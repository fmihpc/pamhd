/*
Particle-assisted magnetohydrodynamics.

Copyright 2015, 2016, 2017 Ilja Honkonen
Copyright 2019 Finnish Meteorological Institute
All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.


Combination of ../particle/test.cpp and ../mhd/two_test.cpp where
particles represent one of the fluids.
*/


#include "array"
#include "cmath"
#include "cstdlib"
#include "fstream"
#include "iostream"
#include "limits"
#include "random"
#include "streambuf"
#include "string"
#include "vector"

#include "boost/filesystem.hpp"
#include "boost/numeric/odeint.hpp"
#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#include "Eigen/Core" // must be included before gensimcell.hpp
#include "Eigen/Geometry"
#include "mpi.h" // must be included before gensimcell.hpp
#include "gensimcell.hpp"
#include "rapidjson/document.h"
#include "rapidjson/error/en.h"

#include "background_magnetic_field.hpp"
#include "boundaries/geometries.hpp"
#include "boundaries/multivariable_boundaries.hpp"
#include "boundaries/multivariable_initial_conditions.hpp"
#include "divergence/remove.hpp"
#include "grid/options.hpp"
#include "mhd/boundaries.hpp"
#include "mhd/common.hpp"
#include "mhd/options.hpp"
#include "mhd/save.hpp"
#include "mhd/N_solve.hpp"
#include "mhd/N_hll_athena.hpp"
#include "mhd/initialize.hpp"
#include "mhd/N_rusanov.hpp"
#include "mhd/variables.hpp"
#include "particle/accumulate_dccrg.hpp"
#include "particle/boundaries.hpp"
#include "particle/common.hpp"
#include "particle/initialize.hpp"
//#include "particle/joiner.hpp"
#include "particle/options.hpp"
#include "particle/save.hpp"
#include "particle/solve_dccrg.hpp"
#include "particle/splitter.hpp"
#include "particle/variables.hpp"
#include "simulation_options.hpp"


using namespace std;
namespace odeint = boost::numeric::odeint;

// counter for assigning unique id to particles
unsigned long long int next_particle_id;


// data stored in every cell of simulation grid
using Cell = gensimcell::Cell<
	gensimcell::Optional_Transfer,
	pamhd::mhd::HD_State_Conservative, // fluid
	pamhd::mhd::HD2_State_Conservative, // particles
	pamhd::Electric_Current_Density,
	pamhd::particle::Solver_Info,
	pamhd::MPI_Rank,
	pamhd::Resistivity,
	pamhd::Magnetic_Field,
	pamhd::Bg_Magnetic_Field_Pos_X,
	pamhd::Bg_Magnetic_Field_Pos_Y,
	pamhd::Bg_Magnetic_Field_Pos_Z,
	pamhd::Magnetic_Field_Resistive,
	pamhd::Magnetic_Field_Temp,
	pamhd::Magnetic_Field_Divergence,
	pamhd::Scalar_Potential_Gradient,
	pamhd::particle::Electric_Field,
	pamhd::particle::Number_Of_Particles,
	pamhd::particle::Bdy_Number_Density,
	pamhd::particle::Bdy_Velocity,
	pamhd::particle::Bdy_Temperature,
	pamhd::particle::Bdy_Species_Mass,
	pamhd::particle::Bdy_Charge_Mass_Ratio,
	pamhd::particle::Bdy_Nr_Particles_In_Cell,
	pamhd::particle::Bulk_Mass,
	pamhd::particle::Bulk_Momentum,
	pamhd::particle::Bulk_Velocity,
	pamhd::particle::Current_Minus_Velocity,
	pamhd::particle::Bulk_Relative_Velocity2,
	pamhd::particle::Nr_Particles_Internal,
	pamhd::particle::Nr_Particles_External,
	pamhd::particle::Nr_Accumulated_To_Cells,
	pamhd::particle::Particles_Internal,
	pamhd::particle::Particles_External,
	pamhd::particle::Accumulated_To_Cells,
	pamhd::mhd::HD_Flux_Conservative,
	pamhd::mhd::HD2_Flux_Conservative,
	pamhd::Magnetic_Field_Flux
>;

// simulation data, see doi:10.1016/j.cpc.2012.12.017 or arxiv.org/abs/1212.3496
using Grid = dccrg::Dccrg<
	Cell,
	dccrg::Cartesian_Geometry,
	std::tuple<>,
	std::tuple<pamhd::particle::Is_Local>
>;

// returns a reference to cell's list of particles not moving to another cell
const auto Part_Int
	= [](Cell& cell_data)->typename pamhd::particle::Particles_Internal::data_type&{
		return cell_data[pamhd::particle::Particles_Internal()];
	};
// particles moving to another cell
const auto Part_Ext
	= [](Cell& cell_data)->typename pamhd::particle::Particles_External::data_type&{
		return cell_data[pamhd::particle::Particles_External()];
	};
// number of particles in above list, for allocating memory for arriving particles
const auto Nr_Ext
	= [](Cell& cell_data)->typename pamhd::particle::Nr_Particles_External::data_type&{
		return cell_data[pamhd::particle::Nr_Particles_External()];
	};

// references to initial condition & boundary data of cell
const auto Bdy_N
	= [](Cell& cell_data)->typename pamhd::particle::Bdy_Number_Density::data_type&{
		return cell_data[pamhd::particle::Bdy_Number_Density()];
	};
const auto Bdy_V
	= [](Cell& cell_data)->typename pamhd::particle::Bdy_Velocity::data_type&{
		return cell_data[pamhd::particle::Bdy_Velocity()];
	};
const auto Bdy_T
	= [](Cell& cell_data)->typename pamhd::particle::Bdy_Temperature::data_type&{
		return cell_data[pamhd::particle::Bdy_Temperature()];
	};
const auto Bdy_Nr_Par
	= [](Cell& cell_data)->typename pamhd::particle::Bdy_Nr_Particles_In_Cell::data_type&{
		return cell_data[pamhd::particle::Bdy_Nr_Particles_In_Cell()];
	};
const auto Bdy_SpM
	= [](Cell& cell_data)->typename pamhd::particle::Bdy_Species_Mass::data_type&{
		return cell_data[pamhd::particle::Bdy_Species_Mass()];
	};
const auto Bdy_C2M
	= [](Cell& cell_data)->typename pamhd::particle::Bdy_Charge_Mass_Ratio::data_type&{
		return cell_data[pamhd::particle::Bdy_Charge_Mass_Ratio()];
	};

// given a particle these return references to particle's parameters
const auto Part_Pos
	= [](pamhd::particle::Particle_Internal& particle)->typename pamhd::particle::Position::data_type&{
		return particle[pamhd::particle::Position()];
	};
const auto Part_Vel
	= [](pamhd::particle::Particle_Internal& particle)->typename pamhd::particle::Velocity::data_type&{
		return particle[pamhd::particle::Velocity()];
	};
// as above but for caller that also provides cell's data
const auto Part_Vel_Cell
	= [](Cell&, pamhd::particle::Particle_Internal& particle)
		-> typename pamhd::particle::Velocity::data_type&
	{
		return particle[pamhd::particle::Velocity()];
	};
const auto Part_C2M
	= [](pamhd::particle::Particle_Internal& particle)->typename pamhd::particle::Charge_Mass_Ratio::data_type&{
		return particle[pamhd::particle::Charge_Mass_Ratio()];
	};
// copy of number of real particles represented by simulation particle
const auto Part_Nr
	= [](pamhd::particle::Particle_Internal& particle)->typename pamhd::particle::Mass::data_type{
		return particle[pamhd::particle::Mass()] / particle[pamhd::particle::Species_Mass()];
	};
const auto Part_Mas
	= [](pamhd::particle::Particle_Internal& particle)->typename pamhd::particle::Mass::data_type&{
		return particle[pamhd::particle::Mass()];
	};
// as above but for caller that also provides cell's data
const auto Part_Mas_Cell
	= [](Cell&, pamhd::particle::Particle_Internal& particle)->typename pamhd::particle::Mass::data_type&{
		return particle[pamhd::particle::Mass()];
	};
const auto Part_Des
	= [](pamhd::particle::Particle_External& particle)->typename pamhd::particle::Destination_Cell::data_type&{
		return particle[pamhd::particle::Destination_Cell()];
	};
// reference to mass of given particle's species
const auto Part_SpM
	= [](pamhd::particle::Particle_Internal& particle)->typename pamhd::particle::Species_Mass::data_type&{
		return particle[pamhd::particle::Species_Mass()];
	};
// as above but for caller that also provides cell's data
const auto Part_SpM_Cell
	= [](Cell&, pamhd::particle::Particle_Internal& particle)->typename pamhd::particle::Species_Mass::data_type&{
		return particle[pamhd::particle::Species_Mass()];
	};
// copy of particle's kinetic energy relative to pamhd::particle::Bulk_Velocity
const auto Part_Ekin
	= [](Cell& cell_data, pamhd::particle::Particle_Internal& particle)->double{
		return
			0.5 * particle[pamhd::particle::Mass()]
			* (
				particle[pamhd::particle::Velocity()]
				- cell_data[pamhd::particle::Bulk_Velocity()].first
			).squaredNorm();
	};

// reference to accumulated number of particles in given cell
const auto Nr_Particles
	= [](Cell& cell_data)->typename pamhd::particle::Number_Of_Particles::data_type&{
		return cell_data[pamhd::particle::Number_Of_Particles()];
	};

const auto Bulk_Mass_Getter
	= [](Cell& cell_data)->typename pamhd::particle::Bulk_Mass::data_type&{
		return cell_data[pamhd::particle::Bulk_Mass()];
	};

const auto Bulk_Momentum_Getter
	= [](Cell& cell_data)->typename pamhd::particle::Bulk_Momentum::data_type&{
		return cell_data[pamhd::particle::Bulk_Momentum()];
	};

const auto Bulk_Relative_Velocity2_Getter
	= [](Cell& cell_data)->typename pamhd::particle::Bulk_Relative_Velocity2::data_type&{
		return cell_data[pamhd::particle::Bulk_Relative_Velocity2()];
	};

// accumulated momentum / accumulated velocity of particles in given cell
const auto Bulk_Velocity_Getter
	= [](Cell& cell_data)->typename pamhd::particle::Bulk_Velocity::data_type&{
		return cell_data[pamhd::particle::Bulk_Velocity()];
	};

// list of items (variables above) accumulated from particles in given cell
const auto Accu_List_Getter
	= [](Cell& cell_data)->typename pamhd::particle::Accumulated_To_Cells::data_type&{
		return cell_data[pamhd::particle::Accumulated_To_Cells()];
	};

// length of above list (for transferring between processes)
const auto Accu_List_Length_Getter
	= [](Cell& cell_data)->typename pamhd::particle::Nr_Accumulated_To_Cells::data_type&{
		return cell_data[pamhd::particle::Nr_Accumulated_To_Cells()];
	};

// target cell of accumulated particle values in an accumulation list item
const auto Accu_List_Target_Getter
	= [](pamhd::particle::Accumulated_To_Cell& accu_item)
		->typename pamhd::particle::Target::data_type&
	{
		return accu_item[pamhd::particle::Target()];
	};

// accumulated number of particles in an accumulation list item
const auto Accu_List_Number_Of_Particles_Getter
	= [](pamhd::particle::Accumulated_To_Cell& accu_item)
		->typename pamhd::particle::Number_Of_Particles::data_type&
	{
		return accu_item[pamhd::particle::Number_Of_Particles()];
	};

const auto Accu_List_Bulk_Mass_Getter
	= [](pamhd::particle::Accumulated_To_Cell& accu_item)
		->typename pamhd::particle::Bulk_Mass::data_type&
	{
		return accu_item[pamhd::particle::Bulk_Mass()];
	};

const auto Accu_List_Bulk_Velocity_Getter
	= [](pamhd::particle::Accumulated_To_Cell& accu_item)
		->typename pamhd::particle::Bulk_Velocity::data_type&
	{
		return accu_item[pamhd::particle::Bulk_Velocity()];
	};

const auto Accu_List_Bulk_Relative_Velocity2_Getter
	= [](pamhd::particle::Accumulated_To_Cell& accu_item)
		->typename pamhd::particle::Bulk_Relative_Velocity2::data_type&
	{
		return accu_item[pamhd::particle::Bulk_Relative_Velocity2()];
	};

const auto Mag
	= [](Cell& cell_data)->typename pamhd::Magnetic_Field::data_type&{
		return cell_data[pamhd::Magnetic_Field()];
	};

// field before divergence removal in case removal fails
const auto Mag_tmp
	= [](Cell& cell_data)->typename pamhd::Magnetic_Field_Temp::data_type&{
		return cell_data[pamhd::Magnetic_Field_Temp()];
	};
// divergence of magnetic field
const auto Mag_div
	= [](Cell& cell_data)->typename pamhd::Magnetic_Field_Divergence::data_type&{
		return cell_data[pamhd::Magnetic_Field_Divergence()];
	};
// electrical resistivity
const auto Res
	= [](Cell& cell_data)->typename pamhd::Resistivity::data_type&{
		return cell_data[pamhd::Resistivity()];
	};
// adjustment to magnetic field due to resistivity
const auto Mag_res
	= [](Cell& cell_data)->typename pamhd::Magnetic_Field_Resistive::data_type&{
		return cell_data[pamhd::Magnetic_Field_Resistive()];
	};
// curl of magnetic field
const auto Cur
	= [](Cell& cell_data)->typename pamhd::Electric_Current_Density::data_type&{
		return cell_data[pamhd::Electric_Current_Density()];
	};
// electric current minus bulk velocity
const auto J_m_V
	= [](Cell& cell_data)->typename pamhd::particle::Current_Minus_Velocity::data_type&{
		return cell_data[pamhd::particle::Current_Minus_Velocity()];
	};
// electric field for propagating particles
const auto Ele
	= [](Cell& cell_data)->typename pamhd::particle::Electric_Field::data_type&{
		return cell_data[pamhd::particle::Electric_Field()];
	};
// solver info variable for boundary logic
const auto Sol_Info
	= [](Cell& cell_data)->typename pamhd::particle::Solver_Info::data_type&{
		return cell_data[pamhd::particle::Solver_Info()];
	};


// references to background magnetic fields
const auto Bg_B_Pos_X
	= [](Cell& cell_data)->typename pamhd::Bg_Magnetic_Field_Pos_X::data_type&{
		return cell_data[pamhd::Bg_Magnetic_Field_Pos_X()];
	};
const auto Bg_B_Pos_Y
	= [](Cell& cell_data)->typename pamhd::Bg_Magnetic_Field_Pos_Y::data_type&{
		return cell_data[pamhd::Bg_Magnetic_Field_Pos_Y()];
	};
const auto Bg_B_Pos_Z
	= [](Cell& cell_data)->typename pamhd::Bg_Magnetic_Field_Pos_Z::data_type&{
		return cell_data[pamhd::Bg_Magnetic_Field_Pos_Z()];
	};

// flux / total change of magnetic field over one time step
const auto Mag_f
	= [](Cell& cell_data)->typename pamhd::Magnetic_Field_Flux::data_type&{
		return cell_data[pamhd::Magnetic_Field_Flux()];
	};

// reference to mass density of fluid 1 in given cell
const auto Mas1
	= [](Cell& cell_data)->typename pamhd::mhd::Mass_Density::data_type&{
		return cell_data[pamhd::mhd::HD_State_Conservative()][pamhd::mhd::Mass_Density()];
	};
const auto Mom1
	= [](Cell& cell_data)->typename pamhd::mhd::Momentum_Density::data_type&{
		return cell_data[pamhd::mhd::HD_State_Conservative()][pamhd::mhd::Momentum_Density()];
	};
const auto Nrj1
	= [](Cell& cell_data)->typename pamhd::mhd::Total_Energy_Density::data_type&{
		return cell_data[pamhd::mhd::HD_State_Conservative()][pamhd::mhd::Total_Energy_Density()];
	};
// reference to mass density of fluid 2 in given cell
const auto Mas2
	= [](Cell& cell_data)->typename pamhd::mhd::Mass_Density::data_type&{
		return cell_data[pamhd::mhd::HD2_State_Conservative()][pamhd::mhd::Mass_Density()];
	};
const auto Mom2
	= [](Cell& cell_data)->typename pamhd::mhd::Momentum_Density::data_type&{
		return cell_data[pamhd::mhd::HD2_State_Conservative()][pamhd::mhd::Momentum_Density()];
	};
const auto Nrj2
	= [](Cell& cell_data)->typename pamhd::mhd::Total_Energy_Density::data_type&{
		return cell_data[pamhd::mhd::HD2_State_Conservative()][pamhd::mhd::Total_Energy_Density()];
	};

// flux of mass density of fluid 1 over one time step
const auto Mas1_f
	= [](Cell& cell_data)->typename pamhd::mhd::Mass_Density::data_type&{
		return cell_data[pamhd::mhd::HD_Flux_Conservative()][pamhd::mhd::Mass_Density()];
	};
const auto Mom1_f
	= [](Cell& cell_data)->typename pamhd::mhd::Momentum_Density::data_type&{
		return cell_data[pamhd::mhd::HD_Flux_Conservative()][pamhd::mhd::Momentum_Density()];
	};
const auto Nrj1_f
	= [](Cell& cell_data)->typename pamhd::mhd::Total_Energy_Density::data_type&{
		return cell_data[pamhd::mhd::HD_Flux_Conservative()][pamhd::mhd::Total_Energy_Density()];
	};
// flux of mass density of fluid 2 over one time step
const auto Mas2_f
	= [](Cell& cell_data)->typename pamhd::mhd::Mass_Density::data_type&{
		return cell_data[pamhd::mhd::HD2_Flux_Conservative()][pamhd::mhd::Mass_Density()];
	};
const auto Mom2_f
	= [](Cell& cell_data)->typename pamhd::mhd::Momentum_Density::data_type&{
		return cell_data[pamhd::mhd::HD2_Flux_Conservative()][pamhd::mhd::Momentum_Density()];
	};
const auto Nrj2_f
	= [](Cell& cell_data)->typename pamhd::mhd::Total_Energy_Density::data_type&{
		return cell_data[pamhd::mhd::HD2_Flux_Conservative()][pamhd::mhd::Total_Energy_Density()];
	};


int main(int argc, char* argv[])
{
	using std::asin;
	using std::atan2;
	using std::get;
	using std::min;
	using std::sqrt;
	using std::to_string;

	/*
	Initialize MPI
	*/

	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
		std::cerr << "Couldn't initialize MPI." << std::endl;
		abort();
	}

	MPI_Comm comm = MPI_COMM_WORLD;

	int rank = 0, comm_size = 0;
	if (MPI_Comm_rank(comm, &rank) != MPI_SUCCESS) {
		std::cerr << "Couldn't obtain MPI rank." << std::endl;
		abort();
	}
	if (MPI_Comm_size(comm, &comm_size) != MPI_SUCCESS) {
		std::cerr << "Couldn't obtain size of MPI communicator." << std::endl;
		abort();
	}

	next_particle_id = 1 + rank;

	// intialize Zoltan
	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
		std::cerr << "Zoltan_Initialize failed." << std::endl;
		abort();
	}

	/*
	Parse configuration file
	*/

	if (argc != 2) {
		if (argc < 2 and rank == 0) {
			std::cerr
				<< "Name of configuration file required."
				<< std::endl;
		}
		if (argc > 2 and rank == 0) {
			std::cerr
				<< "Too many arguments given to " << argv[0]
				<< ": " << argc - 1 << ", should be 1"
				<< std::endl;
		}
		MPI_Finalize();
		return EXIT_FAILURE;
	}

	std::ifstream json_file(argv[1]);
	if (not json_file.good()) {
		if (rank == 0) {
			std::cerr << "Couldn't open configuration file " << argv[1] << std::endl;
		}
		MPI_Finalize();
		return EXIT_FAILURE;
	}

	std::string json{
		std::istreambuf_iterator<char>(json_file),
		std::istreambuf_iterator<char>()
	};

	rapidjson::Document document;
	document.Parse(json.c_str());
	if (document.HasParseError()) {
		std::cerr << "Couldn't parse json data in file " << argv[1]
			<< " at character position " << document.GetErrorOffset()
			<< ": " << rapidjson::GetParseError_En(document.GetParseError())
			<< std::endl;
		MPI_Finalize();
		return EXIT_FAILURE;
	}

	pamhd::Options options_sim{document};
	pamhd::grid::Options options_grid{document};
	pamhd::mhd::Options options_mhd{document};
	pamhd::particle::Options options_particle{document};

	if (rank == 0 and options_sim.output_directory != "") {
		std::cout << "Saving results into directory " << options_sim.output_directory << std::endl;
		try {
			boost::filesystem::create_directories(options_sim.output_directory);
		} catch (const boost::filesystem::filesystem_error& e) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << ") "
				"Couldn't create output directory "
				<< options_sim.output_directory << ": "
				<< e.what()
				<< std::endl;
			abort();
		}
	}

	const int particle_stepper = [&](){
		if (options_particle.solver == "euler") {
			return 0;
		} else if (options_particle.solver == "midpoint") {
			return 1;
		} else if (options_particle.solver == "rk4") {
			return 2;
		} else if (options_particle.solver == "rkck54") {
			return 3;
		} else if (options_particle.solver == "rkf78") {
			return 4;
		} else {
			std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
				<< "Unsupported solver: " << options_particle.solver
				<< ", should be one of: euler, (modified) midpoint, rk4 (runge_kutta4), rkck54 (runge_kutta_cash_karp54), rkf78 (runge_kutta_fehlberg78), see http://www.boost.org/doc/libs/release/libs/numeric/odeint/doc/html/boost_numeric_odeint/odeint_in_detail/steppers.html#boost_numeric_odeint.odeint_in_detail.steppers.stepper_overview"
				<< std::endl;
			abort();
		}
	}();

	using geometry_id_t = unsigned int;

	pamhd::boundaries::Geometries<
		geometry_id_t,
		std::array<double, 3>,
		double,
		uint64_t
	> geometries;
	geometries.set(document);

	/*
	Initial and boundary conditions for fluid 1 and magnetic field
	*/
	pamhd::boundaries::Multivariable_Initial_Conditions<
		geometry_id_t,
		pamhd::mhd::Number_Density,
		pamhd::mhd::Velocity,
		pamhd::mhd::Pressure,
		pamhd::Magnetic_Field
	> initial_conditions_fluid;
	initial_conditions_fluid.set(document);

	pamhd::boundaries::Multivariable_Boundaries<
		uint64_t,
		geometry_id_t,
		pamhd::mhd::Number_Density,
		pamhd::mhd::Velocity,
		pamhd::mhd::Pressure,
		pamhd::Magnetic_Field
	> boundaries_fluid;
	boundaries_fluid.set(document);

	/*
	Initial and boundary conditions for particles (fluid 2)
	*/
	std::vector<
		pamhd::boundaries::Multivariable_Initial_Conditions<
			geometry_id_t,
			pamhd::particle::Bdy_Number_Density,
			pamhd::particle::Bdy_Temperature,
			pamhd::particle::Bdy_Velocity,
			pamhd::particle::Bdy_Nr_Particles_In_Cell,
			pamhd::particle::Bdy_Charge_Mass_Ratio,
			pamhd::particle::Bdy_Species_Mass
		>
	> initial_conditions_particles;
	for (size_t population_id = 0; population_id < 99; population_id++) {
		const auto& obj_population_i
			= document.FindMember(
				(
					"particle-population-"
					+ to_string(population_id)
				).c_str()
			);
		if (obj_population_i == document.MemberEnd()) {
			if (population_id == 0) {
				continue; // allow population ids to start from 0 and 1
			} else {
				break;
			}
		}
		const auto& obj_population = obj_population_i->value;

		const auto old_size = initial_conditions_particles.size();
		initial_conditions_particles.resize(old_size + 1);
		initial_conditions_particles[old_size].set(obj_population);
	}

	std::vector<
		pamhd::boundaries::Multivariable_Boundaries<
			uint64_t,
			geometry_id_t,
			pamhd::particle::Bdy_Number_Density,
			pamhd::particle::Bdy_Temperature,
			pamhd::particle::Bdy_Velocity,
			pamhd::particle::Bdy_Nr_Particles_In_Cell,
			pamhd::particle::Bdy_Charge_Mass_Ratio,
			pamhd::particle::Bdy_Species_Mass
		>
	> boundaries_particles;
	for (size_t population_id = 0; population_id < 99; population_id++) {
		const auto& obj_population_i
			= document.FindMember(
				(
					"particle-population-"
					+ to_string(population_id)
				).c_str()
			);
		if (obj_population_i == document.MemberEnd()) {
			if (population_id == 0) {
				continue; // allow population ids to start from 0 and 1
			} else {
				break;
			}
		}
		const auto& obj_population = obj_population_i->value;

		const auto old_size = boundaries_particles.size();
		boundaries_particles.resize(old_size + 1);
		boundaries_particles[old_size].set(obj_population);
	}

	pamhd::Background_Magnetic_Field<
		double,
		pamhd::Magnetic_Field::data_type
	> background_B;
	background_B.set(document);

	const auto mhd_solver
		= [&options_mhd, &background_B, &rank](){
			if (options_mhd.solver == "rusanov") {
				return pamhd::mhd::Solver::rusanov;
			} else if (options_mhd.solver == "hll-athena") {
				return pamhd::mhd::Solver::hll_athena;
			} else {
				std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
					<< "Unsupported solver: " << options_mhd.solver
					<< std::endl;
				abort();
			}
		}();

	/*
	Prepare resistivity
	*/

	pamhd::boundaries::Math_Expression<pamhd::Resistivity> resistivity;
	mup::Value J_val;
	mup::Variable J_var(&J_val);
	resistivity.add_expression_variable("J", J_var);

	const auto resistivity_name = pamhd::Resistivity::get_option_name();
	if (not document.HasMember(resistivity_name.c_str())) {
		if (rank == 0) {
			std::cerr << __FILE__ "(" << __LINE__
				<< "): Configuration file doesn't have a "
				<< resistivity_name << " key."
				<< std::endl;
		};
		MPI_Finalize();
		return EXIT_FAILURE;
	}
	const auto& json_resistivity = document[resistivity_name.c_str()];
	if (not json_resistivity.IsString()) {
		if (rank == 0) {
			std::cerr << __FILE__ "(" << __LINE__
				<< "): Resistivity option is not of type string."
				<< std::endl;
		};
		MPI_Finalize();
		return EXIT_FAILURE;
	}

	resistivity.set_expression(json_resistivity.GetString());


	/*
	Initialize simulation grid
	*/
	const unsigned int neighborhood_size = 1;
	const auto& number_of_cells = options_grid.get_number_of_cells();
	const auto& periodic = options_grid.get_periodic();

	Grid grid; grid
		.set_neighborhood_length(neighborhood_size)
		.set_maximum_refinement_level(0)
		.set_load_balancing_method(options_sim.lb_name.c_str())
		.set_periodic(periodic[0], periodic[1], periodic[2])
		.set_initial_length(number_of_cells)
		.initialize(comm)
		.balance_load();

	// set grid geometry
	const std::array<double, 3>
		simulation_volume
			= options_grid.get_volume(),
		cell_volume{
			simulation_volume[0] / number_of_cells[0],
			simulation_volume[1] / number_of_cells[1],
			simulation_volume[2] / number_of_cells[2]
		};

	dccrg::Cartesian_Geometry::Parameters geom_params;
	geom_params.start = options_grid.get_start();
	geom_params.level_0_cell_length = cell_volume;

	grid.set_geometry(geom_params);

	// update owner process of cells for saving into file
	for (auto& cell: grid.local_cells()) {
		(*cell.data)[pamhd::MPI_Rank()] = rank;
	}

	// assign cells into boundary geometries
	for (const auto& cell: grid.local_cells()) {
		const auto
			start = grid.geometry.get_min(cell.id),
			end = grid.geometry.get_max(cell.id);
		geometries.overlaps(start, end, cell.id);
	}


	/*
	Simulate
	*/

	const double time_end = options_sim.time_start + options_sim.time_length;
	double
		max_dt_mhd = 0,
		max_dt_particle_gyro = 0,
		max_dt_particle_flight = 0,
		simulation_time = options_sim.time_start,
		next_particle_save = options_particle.save_n,
		next_mhd_save = options_mhd.save_n;

	if (rank == 0) {
		cout << "Initializing simulation... " << endl;
	}

	// zero B before initializing fluids to get correct total energy
	for (auto& cell: grid.local_cells()) {
		Mag(*cell.data) = {0, 0, 0};
	}

	// fluid 1
	pamhd::mhd::initialize_fluid(
		geometries,
		initial_conditions_fluid,
		grid,
		simulation_time,
		options_sim.adiabatic_index,
		options_sim.vacuum_permeability,
		options_sim.proton_mass,
		true,
		Mas1, Mom1, Nrj1, Mag,
		Mas1_f, Mom1_f, Nrj1_f
	);

	// particles
	std::mt19937_64 random_source;

	unsigned long long int nr_particles_created = 0;
	for (auto& init_cond_part: initial_conditions_particles) {
		nr_particles_created = pamhd::particle::initialize_particles<
			pamhd::particle::Particle_Internal,
			pamhd::particle::Mass,
			pamhd::particle::Charge_Mass_Ratio,
			pamhd::particle::Position,
			pamhd::particle::Velocity,
			pamhd::particle::Particle_ID,
			pamhd::particle::Species_Mass
		>(
			geometries,
			init_cond_part,
			simulation_time,
			grid,
			random_source,
			options_particle.boltzmann,
			next_particle_id,
			grid.get_comm_size(),
			false,
			true,
			Part_Int,
			Bdy_N,
			Bdy_V,
			Bdy_T,
			Bdy_Nr_Par,
			Bdy_SpM,
			Bdy_C2M,
			Sol_Info
		);
		next_particle_id += nr_particles_created * grid.get_comm_size();
	}

	nr_particles_created
		+= pamhd::particle::apply_boundaries<
			pamhd::particle::Particle_Internal,
			pamhd::particle::Mass,
			pamhd::particle::Charge_Mass_Ratio,
			pamhd::particle::Position,
			pamhd::particle::Velocity,
			pamhd::particle::Particle_ID,
			pamhd::particle::Species_Mass
		>(
			geometries,
			boundaries_particles,
			simulation_time,
			0,
			grid,
			random_source,
			options_particle.boltzmann,
			options_sim.vacuum_permeability,
			next_particle_id,
			grid.get_comm_size(),
			true,
			Sol_Info,
			Part_Int,
			Bdy_N,
			Bdy_V,
			Bdy_T,
			Bdy_Nr_Par,
			Bdy_SpM,
			Bdy_C2M
		);
	next_particle_id += nr_particles_created * grid.get_comm_size();

	// fluid 2 from particles
	try {
		pamhd::particle::accumulate_mhd_data(
			grid,
			Part_Int,
			Part_Pos,
			Part_Mas_Cell,
			Part_SpM_Cell,
			Part_Vel_Cell,
			Part_Ekin,
			Nr_Particles,
			Part_Nr,
			Bulk_Mass_Getter,
			Bulk_Momentum_Getter,
			Bulk_Relative_Velocity2_Getter,
			Bulk_Velocity_Getter,
			Accu_List_Number_Of_Particles_Getter,
			Accu_List_Bulk_Mass_Getter,
			Accu_List_Bulk_Velocity_Getter,
			Accu_List_Bulk_Relative_Velocity2_Getter,
			Accu_List_Target_Getter,
			Accu_List_Length_Getter,
			Accu_List_Getter,
			pamhd::particle::Nr_Accumulated_To_Cells(),
			pamhd::particle::Accumulated_To_Cells(),
			pamhd::particle::Bulk_Velocity(),
			Sol_Info
		);
	} catch (const std::exception& e) {
		std::cerr << __FILE__ "(" << __LINE__ << ": "
			<< "Couldn't accumulate MHD data from particles: " << e.what()
			<< std::endl;
		abort();
	}

	try {
		pamhd::particle::fill_mhd_fluid_values(
			grid,
			options_sim.adiabatic_index,
			options_sim.vacuum_permeability,
			options_particle.boltzmann,
			options_mhd.min_pressure,
			Nr_Particles,
			Bulk_Mass_Getter,
			Bulk_Momentum_Getter,
			Bulk_Relative_Velocity2_Getter,
			Part_Int,
			Mas2, Mom2, Nrj2, Mag,
			Sol_Info
		);
	} catch (const std::exception& e) {
		std::cerr << __FILE__ "(" << __LINE__ << ": "
			<< "Couldn't fill MHD fluid values: " << e.what()
			<< std::endl;
		abort();
	}

	// magnetic field
	pamhd::mhd::initialize_magnetic_field<pamhd::Magnetic_Field>(
		geometries,
		initial_conditions_fluid,
		background_B,
		grid,
		simulation_time,
		options_sim.vacuum_permeability,
		Mag, Mag_f,
		Bg_B_Pos_X, Bg_B_Pos_Y, Bg_B_Pos_Z
	);
	// update background field between processes
	Cell::set_transfer_all(
		true,
		pamhd::Bg_Magnetic_Field_Pos_X(),
		pamhd::Bg_Magnetic_Field_Pos_Y(),
		pamhd::Bg_Magnetic_Field_Pos_Z()
	);
	grid.update_copies_of_remote_neighbors();
	Cell::set_transfer_all(
		false,
		pamhd::Bg_Magnetic_Field_Pos_X(),
		pamhd::Bg_Magnetic_Field_Pos_Y(),
		pamhd::Bg_Magnetic_Field_Pos_Z()
	);

	// initialize resistivity
	for (auto& cell: grid.local_cells()) {
		Res(*cell.data) = 0;
	}

	pamhd::mhd::apply_magnetic_field_boundaries(
		grid,
		boundaries_fluid,
		geometries,
		simulation_time,
		Mag
	);

	// add magnetic field contribution to fluids' total energy densities
	for (const auto& cell: grid.local_cells()) {
		const auto
			total_mass = Mas1(*cell.data) + Mas2(*cell.data),
			mass_frac1 = Mas1(*cell.data) / total_mass,
			mass_frac2 = Mas2(*cell.data) / total_mass;
		Nrj1(*cell.data) += mass_frac1 * 0.5 * Mag(*cell.data).squaredNorm() / options_sim.vacuum_permeability;
		Nrj2(*cell.data) += mass_frac2 * 0.5 * Mag(*cell.data).squaredNorm() / options_sim.vacuum_permeability;
	}

	if (rank == 0) {
		cout << "done." << endl;
	}

	/*
	Classify cells into normal, boundary and dont_solve
	*/
	Cell::set_transfer_all(true, pamhd::particle::Solver_Info());
	pamhd::mhd::set_solver_info<pamhd::mhd::Solver_Info>(
		grid, boundaries_fluid, geometries, Sol_Info
	);
	pamhd::particle::set_solver_info<pamhd::particle::Solver_Info>(
		grid, boundaries_particles, geometries, Sol_Info
	);
	Cell::set_transfer_all(false, pamhd::particle::Solver_Info());

	// make lists from above for divergence removal functions
	std::vector<uint64_t> solve_cells, bdy_cells, skip_cells;
	for (const auto& cell: grid.local_cells()) {
		if ((Sol_Info(*cell.data) & pamhd::particle::Solver_Info::dont_solve) > 0) {
			skip_cells.push_back(cell.id);
		} else if (Sol_Info(*cell.data) > 0) {
			bdy_cells.push_back(cell.id);
		} else {
			solve_cells.push_back(cell.id);
		}
	}

	size_t simulated_steps = 0;
	while (simulation_time < time_end) {
		simulated_steps++;

		double
			// don't step over the final simulation time
			until_end = time_end - simulation_time,
			// max allowed step for this rank
			local_time_step = min(min(min(
				options_mhd.time_step_factor * max_dt_mhd,
				options_particle.gyroperiod_time_step_factor * max_dt_particle_gyro),
				options_particle.flight_time_step_factor * max_dt_particle_flight),
				until_end),
			time_step = -1;

		if (
			MPI_Allreduce(
				&local_time_step,
				&time_step,
				1,
				MPI_DOUBLE,
				MPI_MIN,
				comm
			) != MPI_SUCCESS
		) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< ": Couldn't reduce time step."
				<< std::endl;
			abort();
		}

		pamhd::particle::split_particles(
			options_particle.min_particles,
			random_source,
			grid,
			Part_Int,
			Part_Pos,
			Part_Mas,
			Sol_Info,
			pamhd::particle::Solver_Info::normal
		);

		try {
			pamhd::particle::accumulate_mhd_data(
				grid,
				Part_Int,
				Part_Pos,
				Part_Mas_Cell,
				Part_SpM_Cell,
				Part_Vel_Cell,
				Part_Ekin,
				Nr_Particles,
				Part_Nr,
				Bulk_Mass_Getter,
				Bulk_Momentum_Getter,
				Bulk_Relative_Velocity2_Getter,
				Bulk_Velocity_Getter,
				Accu_List_Number_Of_Particles_Getter,
				Accu_List_Bulk_Mass_Getter,
				Accu_List_Bulk_Velocity_Getter,
				Accu_List_Bulk_Relative_Velocity2_Getter,
				Accu_List_Target_Getter,
				Accu_List_Length_Getter,
				Accu_List_Getter,
				pamhd::particle::Nr_Accumulated_To_Cells(),
				pamhd::particle::Accumulated_To_Cells(),
				pamhd::particle::Bulk_Velocity(),
				Sol_Info
			);
		} catch (const std::exception& e) {
			std::cerr << __FILE__ "(" << __LINE__ << ": "
				<< "Couldn't accumulate MHD data from particles: " << e.what()
				<< std::endl;
			abort();
		}

		// B required for E calculation
		Cell::set_transfer_all(true, pamhd::Magnetic_Field());
		grid.start_remote_neighbor_copy_updates();

		try {
			pamhd::particle::fill_mhd_fluid_values(
				grid,
				options_sim.adiabatic_index,
				options_sim.vacuum_permeability,
				options_particle.boltzmann,
				options_mhd.min_pressure,
				Nr_Particles,
				Bulk_Mass_Getter,
				Bulk_Momentum_Getter,
				Bulk_Relative_Velocity2_Getter,
				Part_Int,
				Mas2, Mom2, Nrj2, Mag,
				Sol_Info
			);
		} catch (const std::exception& e) {
			std::cerr << __FILE__ "(" << __LINE__ << ": "
				<< "Couldn't fill MHD fluid values: " << e.what()
				<< std::endl;
			abort();
		}

		// inner: J for E = (J - V) x B
		pamhd::divergence::get_curl(
			grid.inner_cells(),
			grid,
			Mag,
			Cur,
			Sol_Info
		);
		// not included in get_curl above
		for (const auto& cell: grid.inner_cells()) {
			Cur(*cell.data) /= options_sim.vacuum_permeability;
		}

		grid.wait_remote_neighbor_copy_update_receives();

		// outer: J for E = (J - V) x B
		pamhd::divergence::get_curl(
			grid.outer_cells(),
			grid,
			Mag,
			Cur,
			Sol_Info
		);
		for (const auto& cell: grid.outer_cells()) {
			Cur(*cell.data) /= options_sim.vacuum_permeability;
		}

		grid.wait_remote_neighbor_copy_update_sends();
		Cell::set_transfer_all(false, pamhd::Magnetic_Field());

		// inner: E = (J - V) x B
		for (const auto& cell: grid.inner_cells()) {
			J_m_V(*cell.data)
				= Cur(*cell.data)
				- (Mom1(*cell.data) + Mom2(*cell.data))
					/ (Mas1(*cell.data) + Mas2(*cell.data));
			Ele(*cell.data) = J_m_V(*cell.data).cross(Mag(*cell.data));
		}

		// outer: E = (J - V) x B
		for (const auto& cell: grid.outer_cells()) {
			J_m_V(*cell.data)
				= Cur(*cell.data)
				- (Mom1(*cell.data) + Mom2(*cell.data))
					/ (Mas1(*cell.data) + Mas2(*cell.data));
			Ele(*cell.data) = J_m_V(*cell.data).cross(Mag(*cell.data));
		}


		Cell::set_transfer_all(true, pamhd::particle::Current_Minus_Velocity());
		grid.update_copies_of_remote_neighbors();
		Cell::set_transfer_all(false, pamhd::particle::Current_Minus_Velocity());


		/*
		Solve
		*/

		max_dt_mhd             =
		max_dt_particle_flight =
		max_dt_particle_gyro   = std::numeric_limits<double>::max();

		if (rank == 0) {
			cout << "Solving at time " << simulation_time
				<< " s with time step " << time_step << " s" << endl;
		}

		// TODO: don't use preprocessor
		#define SOLVE_WITH_STEPPER(given_type, given_cells) \
			pamhd::particle::solve<\
				given_type\
			>(\
				time_step,\
				given_cells,\
				grid,\
				background_B,\
				options_sim.vacuum_permeability,\
				true,\
				J_m_V,\
				Mag,\
				Nr_Ext,\
				Part_Int,\
				Part_Ext,\
				Part_Pos,\
				Part_Vel,\
				Part_C2M,\
				Part_Mas,\
				Part_Des,\
				Sol_Info\
			)

		std::pair<double, double> particle_max_dt{0, 0};

		// outer cells
		switch (particle_stepper) {
		case 0:
			particle_max_dt = SOLVE_WITH_STEPPER(odeint::euler<pamhd::particle::state_t>, grid.outer_cells());
			break;
		case 1:
			particle_max_dt = SOLVE_WITH_STEPPER(odeint::modified_midpoint<pamhd::particle::state_t>, grid.outer_cells());
			break;
		case 2:
			particle_max_dt = SOLVE_WITH_STEPPER(odeint::runge_kutta4<pamhd::particle::state_t>, grid.outer_cells());
			break;
		case 3:
			particle_max_dt = SOLVE_WITH_STEPPER(odeint::runge_kutta_cash_karp54<pamhd::particle::state_t>, grid.outer_cells());
			break;
		case 4:
			particle_max_dt = SOLVE_WITH_STEPPER(odeint::runge_kutta_fehlberg78<pamhd::particle::state_t>, grid.outer_cells());
			break;
		default:
			std::cerr <<  __FILE__ << "(" << __LINE__ << "): " << particle_stepper << std::endl;
			abort();
		}
		max_dt_particle_flight = min(particle_max_dt.first, max_dt_particle_flight);
		max_dt_particle_gyro = min(particle_max_dt.second, max_dt_particle_gyro);

		Cell::set_transfer_all(
			true,
			pamhd::Magnetic_Field(),
			pamhd::mhd::HD_State_Conservative(),
			pamhd::mhd::HD2_State_Conservative(),
			pamhd::particle::Nr_Particles_External()
		);
		grid.start_remote_neighbor_copy_updates();

		// inner MHD
		double solve_max_dt = -1;

		try {
			solve_max_dt = pamhd::mhd::N_solve(
				mhd_solver,
				grid.inner_cells(),
				grid,
				time_step,
				options_sim.adiabatic_index,
				options_sim.vacuum_permeability,
				std::make_pair(Mas1, Mas2),
				std::make_pair(Mom1, Mom2),
				std::make_pair(Nrj1, Nrj2),
				Mag,
				Bg_B_Pos_X, Bg_B_Pos_Y, Bg_B_Pos_Z,
				std::make_pair(Mas1_f, Mas2_f),
				std::make_pair(Mom1_f, Mom2_f),
				std::make_pair(Nrj1_f, Nrj2_f),
				Mag_f,
				Sol_Info
			);
		} catch (const std::exception& e) {
			std::cerr << __FILE__ "(" << __LINE__ << ": "
				<< "MHD solution failed in inner cells: " << e.what()
				<< std::endl;
			abort();
		}
		max_dt_mhd = min(solve_max_dt, max_dt_mhd);

		// inner particles
		switch (particle_stepper) {
		case 0:
			particle_max_dt = SOLVE_WITH_STEPPER(odeint::euler<pamhd::particle::state_t>, grid.inner_cells());
			break;
		case 1:
			particle_max_dt = SOLVE_WITH_STEPPER(odeint::modified_midpoint<pamhd::particle::state_t>, grid.inner_cells());
			break;
		case 2:
			particle_max_dt = SOLVE_WITH_STEPPER(odeint::runge_kutta4<pamhd::particle::state_t>, grid.inner_cells());
			break;
		case 3:
			particle_max_dt = SOLVE_WITH_STEPPER(odeint::runge_kutta_cash_karp54<pamhd::particle::state_t>, grid.inner_cells());
			break;
		case 4:
			particle_max_dt = SOLVE_WITH_STEPPER(odeint::runge_kutta_fehlberg78<pamhd::particle::state_t>, grid.inner_cells());
			break;
		default:
			std::cerr <<  __FILE__ << "(" << __LINE__ << "): " << particle_stepper << std::endl;
			abort();
		}
		#undef SOLVE_WITH_STEPPER
		max_dt_particle_flight = min(particle_max_dt.first, max_dt_particle_flight);
		max_dt_particle_gyro = min(particle_max_dt.second, max_dt_particle_gyro);

		grid.wait_remote_neighbor_copy_update_receives();

		// outer MHD
		try {
			solve_max_dt = pamhd::mhd::N_solve(
				mhd_solver,
				grid.outer_cells(),
				grid,
				time_step,
				options_sim.adiabatic_index,
				options_sim.vacuum_permeability,
				std::make_pair(Mas1, Mas2),
				std::make_pair(Mom1, Mom2),
				std::make_pair(Nrj1, Nrj2),
				Mag,
				Bg_B_Pos_X, Bg_B_Pos_Y, Bg_B_Pos_Z,
				std::make_pair(Mas1_f, Mas2_f),
				std::make_pair(Mom1_f, Mom2_f),
				std::make_pair(Nrj1_f, Nrj2_f),
				Mag_f,
				Sol_Info
			);
		} catch (const std::exception& e) {
			std::cerr << __FILE__ "(" << __LINE__ << ": "
				<< "MHD solution failed in outer cells: " << e.what()
				<< std::endl;
			abort();
		}
		max_dt_mhd = min(solve_max_dt, max_dt_mhd);

		pamhd::divergence::get_curl(
			grid.outer_cells(),
			grid,
			Mag,
			Cur,
			Sol_Info
		);
		for (const auto& cell: grid.outer_cells()) {
			Cur(*cell.data) /= options_sim.vacuum_permeability;
		}

		pamhd::particle::resize_receiving_containers<
			pamhd::particle::Nr_Particles_External,
			pamhd::particle::Particles_External
		>(grid.remote_cells(), grid);

		grid.wait_remote_neighbor_copy_update_sends();
		Cell::set_transfer_all(
			false,
			pamhd::Magnetic_Field(),
			pamhd::mhd::HD_State_Conservative(),
			pamhd::mhd::HD2_State_Conservative(),
			pamhd::particle::Nr_Particles_External()
		);

		// transfer J for calculating additional contributions to B
		Cell::set_transfer_all(true, pamhd::Electric_Current_Density());
		grid.start_remote_neighbor_copy_updates();

		// add contribution to change of B from resistivity
		pamhd::divergence::get_curl(
			grid.inner_cells(),
			grid,
			Cur,
			Mag_res,
			Sol_Info
		);
		for (const auto& cell: grid.inner_cells()) {
			const auto c = grid.geometry.get_center(cell.id);
			const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);

			J_val = Cur(*cell.data).norm();
			Res(*cell.data) = resistivity.evaluate(
				simulation_time,
				c[0], c[1], c[2],
				r, asin(c[2] / r), atan2(c[1], c[0])
			);

			//TODO keep pressure/temperature constant despite electric resistivity
			Mag_res(*cell.data) *= -Res(*cell.data);
			Mag_f(*cell.data) += Mag_res(*cell.data);
		}

		grid.wait_remote_neighbor_copy_update_receives();

		pamhd::divergence::get_curl(
			grid.outer_cells(),
			grid,
			Cur,
			Mag_res,
			Sol_Info
		);
		for (const auto& cell: grid.outer_cells()) {
			const auto c = grid.geometry.get_center(cell.id);
			const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);

			J_val = Cur(*cell.data).norm();
			Res(*cell.data) = resistivity.evaluate(
				simulation_time,
				c[0], c[1], c[2],
				r, asin(c[2] / r), atan2(c[1], c[0])
			);

			Mag_res(*cell.data) *= -Res(*cell.data);
			Mag_f(*cell.data) += Mag_res(*cell.data);
		}

		grid.wait_remote_neighbor_copy_update_sends();
		Cell::set_transfer_all(false, pamhd::Electric_Current_Density());

		Cell::set_transfer_all(true, pamhd::particle::Particles_External());
		grid.start_remote_neighbor_copy_updates();

		try {
			pamhd::mhd::apply_fluxes_N(
				grid,
				options_mhd.min_pressure,
				options_sim.adiabatic_index,
				options_sim.vacuum_permeability,
				std::make_pair(Mas1, Mas2),
				std::make_pair(Mom1, Mom2),
				std::make_pair(Nrj1, Nrj2),
				Mag,
				std::make_pair(Mas1_f, Mas2_f),
				std::make_pair(Mom1_f, Mom2_f),
				std::make_pair(Nrj1_f, Nrj2_f),
				Mag_f,
				Sol_Info
			);
		} catch (const std::exception& e) {
			std::cerr << __FILE__ "(" << __LINE__ << ": "
				<< "Couldn't apply fluxes: " << e.what()
				<< std::endl;
			abort();
		}

		pamhd::particle::incorporate_external_particles<
			pamhd::particle::Nr_Particles_Internal,
			pamhd::particle::Particles_Internal,
			pamhd::particle::Particles_External,
			pamhd::particle::Destination_Cell
		>(grid.inner_cells(), grid);

		grid.wait_remote_neighbor_copy_update_receives();

		pamhd::particle::incorporate_external_particles<
			pamhd::particle::Nr_Particles_Internal,
			pamhd::particle::Particles_Internal,
			pamhd::particle::Particles_External,
			pamhd::particle::Destination_Cell
		>(grid.outer_cells(), grid);

		pamhd::particle::remove_external_particles<
			pamhd::particle::Nr_Particles_External,
			pamhd::particle::Particles_External
		>(grid.inner_cells(), grid);

		grid.wait_remote_neighbor_copy_update_sends();
		Cell::set_transfer_all(false, pamhd::particle::Particles_External());

		pamhd::particle::remove_external_particles<
			pamhd::particle::Nr_Particles_External,
			pamhd::particle::Particles_External
		>(grid.outer_cells(), grid);


		simulation_time += time_step;


		/*
		Update internal particles for setting particle copy boundaries.

		TODO overlap computation and communication in boundary processing
		*/
		for (const auto& cell: grid.local_cells()) {
			// (ab)use external number counter as internal number counter
			(*cell.data)[pamhd::particle::Nr_Particles_External()]
				= (*cell.data)[pamhd::particle::Particles_Internal()].size();
		}
		Cell::set_transfer_all(
			true,
			pamhd::Magnetic_Field(),
			pamhd::mhd::HD_State_Conservative(),
			pamhd::mhd::HD2_State_Conservative(),
			pamhd::particle::Nr_Particles_External()
		);
		grid.update_copies_of_remote_neighbors();
		Cell::set_transfer_all(
			false,
			pamhd::Magnetic_Field(),
			pamhd::mhd::HD_State_Conservative(),
			pamhd::mhd::HD2_State_Conservative(),
			pamhd::particle::Nr_Particles_External()
		);

		pamhd::particle::resize_receiving_containers<
			pamhd::particle::Nr_Particles_External,
			pamhd::particle::Particles_Internal
		>(grid.remote_cells(), grid);
		Cell::set_transfer_all(true, pamhd::particle::Particles_Internal());
		grid.update_copies_of_remote_neighbors();
		Cell::set_transfer_all(false, pamhd::particle::Particles_Internal());

		try {
			pamhd::mhd::apply_magnetic_field_boundaries(
				grid,
				boundaries_fluid,
				geometries,
				simulation_time,
				Mag
			);
		} catch (const std::exception& e) {
			std::cerr << __FILE__ "(" << __LINE__ << ": "
				<< "Couldn't apply magnetic field boundaries: " << e.what()
				<< std::endl;
			abort();
		}

		try {
			pamhd::mhd::apply_fluid_boundaries(
				grid,
				boundaries_fluid,
				geometries,
				simulation_time,
				Mas1, Mom1, Nrj1, Mag,
				options_sim.proton_mass,
				options_sim.adiabatic_index,
				options_sim.vacuum_permeability
			);
		} catch (const std::exception& e) {
			std::cerr << __FILE__ "(" << __LINE__ << ": "
				<< "Couldn't apply fluid boundaries: " << e.what()
				<< std::endl;
			abort();
		}

		nr_particles_created
			+= pamhd::particle::apply_boundaries<
				pamhd::particle::Particle_Internal,
				pamhd::particle::Mass,
				pamhd::particle::Charge_Mass_Ratio,
				pamhd::particle::Position,
				pamhd::particle::Velocity,
				pamhd::particle::Particle_ID,
				pamhd::particle::Species_Mass
			>(
				geometries,
				boundaries_particles,
				simulation_time,
				simulated_steps,
				grid,
				random_source,
				options_particle.boltzmann,
				options_sim.vacuum_permeability,
				next_particle_id,
				grid.get_comm_size(),
				true,
				Sol_Info,
				Part_Int,
				Bdy_N,
				Bdy_V,
				Bdy_T,
				Bdy_Nr_Par,
				Bdy_SpM,
				Bdy_C2M
			);
		next_particle_id += nr_particles_created * grid.get_comm_size();

			// add magnetic nrj to fluid boundary cells
			/*if ((*cell_data)[pamhd::mhd::Cell_Type()] == bdy_classifier_fluid1.value_boundary_cell) {
				const auto& bdy_id = (*cell_data)[pamhd::mhd::Value_Boundary_Id()];

				const auto mass_frac
					= Mas1(*cell_data)
					/ (Mas1(*cell_data) + Mas2(*cell_data));
				Nrj1(*cell_data) += mass_frac * 0.5 * Mag(*cell_data).squaredNorm() / options_sim.vacuum_permeability;
			}*/

		/*
		Save simulation to disk
		*/

		// particles
		if (
			(options_particle.save_n >= 0 and (simulation_time == 0 or simulation_time >= time_end))
			or (options_particle.save_n > 0 and simulation_time >= next_particle_save)
		) {
			if (next_particle_save <= simulation_time) {
				next_particle_save += options_particle.save_n;
			}

			if (rank == 0) {
				cout << "Saving particles at time " << simulation_time << "...";
			}

			if (
				not pamhd::particle::save<
					pamhd::particle::Electric_Field,
					pamhd::Magnetic_Field,
					pamhd::Electric_Current_Density,
					pamhd::particle::Nr_Particles_Internal,
					pamhd::particle::Particles_Internal
				>(
					boost::filesystem::canonical(
						boost::filesystem::path(options_sim.output_directory)
					).append("particle_").generic_string(),
					grid,
					simulation_time,
					options_sim.adiabatic_index,
					options_sim.vacuum_permeability,
					options_particle.boltzmann
				)
			) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
					"Couldn't save particle result."
					<< std::endl;
				MPI_Finalize();
				return EXIT_FAILURE;
			}

			if (rank == 0) {
				cout << "done." << endl;
			}
		}

		// mhd
		if (
			(options_mhd.save_n >= 0 and (simulation_time == 0 or simulation_time >= time_end))
			or (options_mhd.save_n > 0 and simulation_time >= next_mhd_save)
		) {
			if (next_mhd_save <= simulation_time) {
				next_mhd_save += options_mhd.save_n;
			}

			if (rank == 0) {
				cout << "Saving fluid at time " << simulation_time << "... ";
			}

			if (
				not pamhd::mhd::save(
					boost::filesystem::canonical(
						boost::filesystem::path(options_sim.output_directory)
					).append("2mhd_").generic_string(),
					grid,
					2,
					simulation_time,
					options_sim.adiabatic_index,
					options_sim.proton_mass,
					options_sim.vacuum_permeability,
					pamhd::mhd::HD_State_Conservative(),
					pamhd::mhd::HD2_State_Conservative(),
					pamhd::Electric_Current_Density(),
					pamhd::particle::Solver_Info(),
					pamhd::MPI_Rank(),
					pamhd::Resistivity(),
					pamhd::Magnetic_Field(),
					pamhd::Bg_Magnetic_Field_Pos_X(),
					pamhd::Bg_Magnetic_Field_Pos_Y(),
					pamhd::Bg_Magnetic_Field_Pos_Z()
				)
			) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
					"Couldn't save mhd result."
					<< std::endl;
				abort();
			}

			if (rank == 0) {
				cout << "done." << endl;
			}
		}
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
