/*
Hybrid PIC program of PAMHD.

Copyright 2015, 2016, 2017 Ilja Honkonen
Copyright 2019, 2024 Finnish Meteorological Institute
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


Author(s): Ilja Honkonen
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
#include "math/staggered.hpp"
#include "mhd/boundaries.hpp"
#include "mhd/common.hpp"
#include "mhd/hll_athena.hpp"
#include "mhd/hlld_athena.hpp"
#include "mhd/initialize.hpp"
#include "mhd/initialize_staggered.hpp"
#include "mhd/options.hpp"
#include "mhd/roe_athena.hpp"
#include "mhd/rusanov.hpp"
#include "mhd/save.hpp"
#include "mhd/solve_staggered.hpp"
#include "mhd/variables.hpp"
#include "particle/accumulate_dccrg.hpp"
#include "particle/boundaries.hpp"
#include "particle/common.hpp"
#include "particle/initialize.hpp"
#include "particle/options.hpp"
#include "particle/save.hpp"
#include "particle/solve_dccrg.hpp"
#include "particle/splitter.hpp"
#include "particle/variables.hpp"
#include "simulation_options.hpp"
#include "variables.hpp"


namespace odeint = boost::numeric::odeint;

// counter for assigning unique id to particles
unsigned long long int next_particle_id;

using Cell = pamhd::particle::Cell_hyb_particle_staggered;
using Grid = dccrg::Dccrg<
	Cell,
	dccrg::Cartesian_Geometry,
	std::tuple<pamhd::grid::Cell_Is_Local>,
	std::tuple<
		pamhd::grid::Face_Neighbor,
		pamhd::grid::Edge_Neighbor,
		pamhd::grid::Relative_Size,
		pamhd::grid::Neighbor_Is_Local>
>;

// returns reference to background magnetic field on cell faces
const auto Bg_B = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::Bg_Magnetic_Field()];
};

// returns reference to total mass density in given cell
const auto Mas = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::mhd::MHD_State_Conservative()][pamhd::mhd::Mass_Density()];
};
const auto Mom = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::mhd::MHD_State_Conservative()][pamhd::mhd::Momentum_Density()];
};
const auto Nrj = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::mhd::MHD_State_Conservative()][pamhd::mhd::Total_Energy_Density()];
};
const auto Vol_B = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::mhd::MHD_State_Conservative()][pamhd::Magnetic_Field()];
};
const auto Face_B = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::Face_Magnetic_Field()];
};
const auto Face_dB = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::Face_dB()];
};
// divergence of magnetic field
const auto Mag_div = [](Cell& cell_data) ->auto&{
	return cell_data[pamhd::Magnetic_Field_Divergence()];
};
// curl of magnetic field
const auto Vol_J = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::Electric_Current_Density()];
};
// electric current minus bulk velocity
const auto J_m_V = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::particle::Current_Minus_Velocity()];
};
// electric field for propagating particles
const auto Vol_E = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::particle::Electric_Field()];
};

// solver info variable for boundary logic
const auto Sol_Info = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::particle::Solver_Info()];
};
// returns 1 for normal cell, -1 for dont_solve and 0 otherwise
const auto Sol_Info2 = [](Cell& cell_data)->int {
	const auto info = cell_data[pamhd::mhd::Solver_Info()];
	if (info == 0) {
		return 1;
	}
	if ((info & pamhd::mhd::Solver_Info::dont_solve) > 0) {
		return -1;
	}
	return 0;
};

const auto Substep = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::mhd::Substepping_Period()];
};

const auto Substep_Min = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::mhd::Substep_Min()];
};

const auto Substep_Max = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::mhd::Substep_Max()];
};

const auto Timestep = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::mhd::Timestep()];
};

const auto Max_v = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::mhd::Max_Velocity()];
};

// flux of mass density through positive x face of cell
const auto Mas_pfx = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](0, 1)[pamhd::mhd::Mass_Density()];
};
const auto Mas_nfx = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](0, 0)[pamhd::mhd::Mass_Density()];
};
const auto Mas_pfy = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](1, 1)[pamhd::mhd::Mass_Density()];
};
const auto Mas_nfy = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](1, 0)[pamhd::mhd::Mass_Density()];
};
const auto Mas_pfz = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](2, 1)[pamhd::mhd::Mass_Density()];
};
const auto Mas_nfz = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](2, 0)[pamhd::mhd::Mass_Density()];
};
const auto Mom_pfx = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](0, 1)[pamhd::mhd::Momentum_Density()];
};
const auto Mom_nfx = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](0, 0)[pamhd::mhd::Momentum_Density()];
};
const auto Mom_pfy = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](1, 1)[pamhd::mhd::Momentum_Density()];
};
const auto Mom_nfy = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](1, 0)[pamhd::mhd::Momentum_Density()];
};
const auto Mom_pfz = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](2, 1)[pamhd::mhd::Momentum_Density()];
};
const auto Mom_nfz = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](2, 0)[pamhd::mhd::Momentum_Density()];
};
const auto Nrj_pfx = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](0, 1)[pamhd::mhd::Total_Energy_Density()];
};
const auto Nrj_nfx = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](0, 0)[pamhd::mhd::Total_Energy_Density()];
};
const auto Nrj_pfy = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](1, 1)[pamhd::mhd::Total_Energy_Density()];
};
const auto Nrj_nfy = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](1, 0)[pamhd::mhd::Total_Energy_Density()];
};
const auto Nrj_pfz = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](2, 1)[pamhd::mhd::Total_Energy_Density()];
};
const auto Nrj_nfz = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](2, 0)[pamhd::mhd::Total_Energy_Density()];
};
const auto Mag_pfx = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](0, 1)[pamhd::Magnetic_Field()];
};
const auto Mag_nfx = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](0, 0)[pamhd::Magnetic_Field()];
};
const auto Mag_pfy = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](1, 1)[pamhd::Magnetic_Field()];
};
const auto Mag_nfy = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](1, 0)[pamhd::Magnetic_Field()];
};
const auto Mag_pfz = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](2, 1)[pamhd::Magnetic_Field()];
};
const auto Mag_nfz = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](2, 0)[pamhd::Magnetic_Field()];
};

// collections of above to shorten function arguments
const auto Mas_fs = std::make_tuple(
	Mas_nfx, Mas_pfx, Mas_nfy, Mas_pfy, Mas_nfz, Mas_pfz
);
const auto Mom_fs = std::make_tuple(
	Mom_nfx, Mom_pfx, Mom_nfy, Mom_pfy, Mom_nfz, Mom_pfz
);
const auto Nrj_fs = std::make_tuple(
	Nrj_nfx, Nrj_pfx, Nrj_nfy, Nrj_pfy, Nrj_nfz, Nrj_pfz
);
const auto Mag_fs = std::make_tuple(
	Mag_nfx, Mag_pfx, Mag_nfy, Mag_pfy, Mag_nfz, Mag_pfz
);

const auto FInfo = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::mhd::Face_Boundary_Type()];
};
const auto Ref_min = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::grid::Target_Refinement_Level_Min()];
};
const auto Ref_max = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::grid::Target_Refinement_Level_Max()];
};

// list of particles in cell not moving to another cell
const auto Part_Int = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::particle::Particles_Internal()];
};
// particles moving to another cell
const auto Part_Ext = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::particle::Particles_External()];
};
// number of particles in above list, for allocating memory for arriving particles
const auto Nr_Ext = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::particle::Nr_Particles_External()];
};
const auto Nr_Int = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::particle::Nr_Particles_Internal()];
};
// references to initial condition & boundary data of cell
const auto Bdy_N = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::particle::Bdy_Number_Density()];
};
const auto Bdy_V = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::particle::Bdy_Velocity()];
};
const auto Bdy_T = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::particle::Bdy_Temperature()];
};
const auto Bdy_Nr_Par = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::particle::Bdy_Nr_Particles_In_Cell()];
};
const auto Bdy_SpM = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::particle::Bdy_Species_Mass()];
};
const auto Bdy_C2M = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::particle::Bdy_Charge_Mass_Ratio()];
};

// given a particle these return references to particle's parameters
const auto Part_Pos = [](
	pamhd::particle::Particle_Internal& particle
) ->auto& {
	return particle[pamhd::particle::Position()];
};
// as above but for caller that also provides cell's data
const auto Part_Vel_Cell = [](
	Cell&, pamhd::particle::Particle_Internal& particle
) ->auto& {
	return particle[pamhd::particle::Velocity()];
};
const auto Part_Vel = [](
	pamhd::particle::Particle_Internal& particle
) ->auto& {
	return particle[pamhd::particle::Velocity()];
};
const auto Part_C2M = [](
	pamhd::particle::Particle_Internal& particle
) ->auto& {
	return particle[pamhd::particle::Charge_Mass_Ratio()];
};
// copy of number of real particles represented by simulation particle
const auto Part_Nr = [](
	pamhd::particle::Particle_Internal& particle
) ->auto {
	return particle[pamhd::particle::Mass()] / particle[pamhd::particle::Species_Mass()];
};
const auto Part_Mas = [](
	pamhd::particle::Particle_Internal& particle
) ->auto& {
	return particle[pamhd::particle::Mass()];
};
// as above but for caller that also provides cell's data
const auto Part_Mas_Cell = [](
	Cell&, pamhd::particle::Particle_Internal& particle
) ->auto& {
	return particle[pamhd::particle::Mass()];
};
const auto Part_Des = [](
	pamhd::particle::Particle_External& particle
) ->auto& {
	return particle[pamhd::particle::Destination_Cell()];
};
// reference to mass of given particle's species
const auto Part_SpM = [](
	pamhd::particle::Particle_Internal& particle
) ->auto& {
	return particle[pamhd::particle::Species_Mass()];
};
// as above but for caller that also provides cell's data
const auto Part_SpM_Cell = [](
	Cell&, pamhd::particle::Particle_Internal& particle
) ->auto& {
	return particle[pamhd::particle::Species_Mass()];
};
// copy of particle's kinetic energy relative to pamhd::particle::Bulk_Velocity
const auto Part_Ekin = [](
	Cell& cell_data,
	pamhd::particle::Particle_Internal& particle
) ->auto {
	return
		0.5 * particle[pamhd::particle::Mass()]
		* (
			particle[pamhd::particle::Velocity()]
			- cell_data[pamhd::particle::Bulk_Velocity()].first
		).squaredNorm();
};

// reference to accumulated number of particles in given cell
const auto Nr_Particles = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::particle::Number_Of_Particles()];
};

const auto Bulk_Mass_Getter = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::particle::Bulk_Mass()];
};

const auto Bulk_Momentum_Getter = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::particle::Bulk_Momentum()];
};

const auto Bulk_Relative_Velocity2_Getter = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::particle::Bulk_Relative_Velocity2()];
};

const auto Bulk_Velocity_Getter = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::particle::Bulk_Velocity()];
};

// list of items (variables above) accumulated from particles in given cell
const auto Accu_List_Getter = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::particle::Accumulated_To_Cells()];
};

// length of above list (for transferring between processes)
const auto Accu_List_Length_Getter = [](Cell& cell_data) ->auto& {
	return cell_data[pamhd::particle::Nr_Accumulated_To_Cells()];
};

// target cell of accumulated particle values in an accumulation list item
const auto Accu_List_Target_Getter = [](
	pamhd::particle::Accumulated_To_Cell& accu_item
) ->auto& {
	return accu_item[pamhd::particle::Target()];
};

// accumulated number of particles in an accumulation list item
const auto Accu_List_Number_Of_Particles_Getter = [](
	pamhd::particle::Accumulated_To_Cell& accu_item
) ->auto& {
	return accu_item[pamhd::particle::Number_Of_Particles()];
};

const auto Accu_List_Bulk_Mass_Getter = [](
	pamhd::particle::Accumulated_To_Cell& accu_item
) ->auto& {
	return accu_item[pamhd::particle::Bulk_Mass()];
};

const auto Accu_List_Bulk_Velocity_Getter = [](
	pamhd::particle::Accumulated_To_Cell& accu_item
) ->auto& {
	return accu_item[pamhd::particle::Bulk_Velocity()];
};

const auto Accu_List_Bulk_Relative_Velocity2_Getter = [](
	pamhd::particle::Accumulated_To_Cell& accu_item
) ->auto& {
	return accu_item[pamhd::particle::Bulk_Relative_Velocity2()];
};


int main(int argc, char* argv[])
{
	using std::ceil;
	using std::cerr;
	using std::cout;
	using std::endl;
	using std::flush;
	using std::max;
	using std::min;
	using std::to_string;

	/*
	Initialize MPI
	*/

	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
		cerr << "Couldn't initialize MPI." << endl;
		abort();
	}

	MPI_Comm comm = MPI_COMM_WORLD;

	int rank = 0, comm_size = 0;
	if (MPI_Comm_rank(comm, &rank) != MPI_SUCCESS) {
		cerr << "Couldn't obtain MPI rank." << endl;
		abort();
	}
	if (MPI_Comm_size(comm, &comm_size) != MPI_SUCCESS) {
		cerr << "Couldn't obtain size of MPI communicator." << endl;
		abort();
	}

	next_particle_id = 1 + rank;

	// intialize Zoltan
	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
		cerr << "Zoltan_Initialize failed." << endl;
		abort();
	}

	/*
	Parse configuration file
	*/

	if (argc != 2) {
		if (argc < 2 and rank == 0) {
			cerr
				<< "Name of configuration file required."
				<< endl;
		}
		if (argc > 2 and rank == 0) {
			cerr
				<< "Too many arguments given to " << argv[0]
				<< ": " << argc - 1 << ", should be 1"
				<< endl;
		}
		MPI_Finalize();
		return EXIT_FAILURE;
	}

	std::ifstream json_file(argv[1]);
	if (not json_file.good()) {
		if (rank == 0) {
			cerr << "Couldn't open configuration file " << argv[1] << endl;
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
		cerr << "Couldn't parse json data in file " << argv[1]
			<< " at character position " << document.GetErrorOffset()
			<< ": " << rapidjson::GetParseError_En(document.GetParseError())
			<< endl;
		MPI_Finalize();
		return EXIT_FAILURE;
	}

	pamhd::Options options_sim{document};
	pamhd::grid::Options options_grid{document};
	pamhd::mhd::Options options_mhd{document};
	pamhd::particle::Options options_particle{document};

	if (rank == 0 and options_sim.output_directory != "") {
		cout << "Saving results into directory " << options_sim.output_directory << endl;
		try {
			boost::filesystem::create_directories(options_sim.output_directory);
		} catch (const boost::filesystem::filesystem_error& e) {
			cerr <<  __FILE__ << "(" << __LINE__ << ") "
				"Couldn't create output directory "
				<< options_sim.output_directory << ": "
				<< e.what()
				<< endl;
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
			cerr <<  __FILE__ << "(" << __LINE__ << "): "
				<< "Unsupported solver: " << options_particle.solver
				<< ", should be one of: euler, (modified) midpoint, rk4 (runge_kutta4), rkck54 (runge_kutta_cash_karp54), rkf78 (runge_kutta_fehlberg78), see http://www.boost.org/doc/libs/release/libs/numeric/odeint/doc/html/boost_numeric_odeint/odeint_in_detail/steppers.html#boost_numeric_odeint.odeint_in_detail.steppers.stepper_overview"
				<< endl;
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

	pamhd::boundaries::Multivariable_Initial_Conditions<
		geometry_id_t,
		pamhd::Magnetic_Field
	> initial_conditions_fields;
	initial_conditions_fields.set(document);

	pamhd::boundaries::Multivariable_Boundaries<
		uint64_t,
		geometry_id_t,
		pamhd::Magnetic_Field
	> boundaries_fields;
	boundaries_fields.set(document);

	// separate initial and boundary conditions for each particle population
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
			} else if (options_mhd.solver == "hlld-athena") {
				if (background_B.exists() and rank == 0) {
					cout << "NOTE: background magnetic field ignored by hlld-athena solver." << endl;
				}
				return pamhd::mhd::Solver::hlld_athena;
			} else if (options_mhd.solver == "roe-athena") {
				if (background_B.exists() and rank == 0) {
					cout << "NOTE: background magnetic field ignored by roe-athena solver." << endl;
				}
				return pamhd::mhd::Solver::roe_athena;
			} else {
				cerr <<  __FILE__ << "(" << __LINE__ << "): "
					<< "Unsupported solver: " << options_mhd.solver
					<< endl;
				abort();
			}
		}();


	/*
	Initialize simulation grid
	*/
	const unsigned int neighborhood_size = 3;
	const auto& number_of_cells = options_grid.get_number_of_cells();
	const auto& periodic = options_grid.get_periodic();

	Grid grid; grid
		.set_initial_length(number_of_cells)
		.set_neighborhood_length(neighborhood_size)
		.set_periodic(periodic[0], periodic[1], periodic[2])
		.set_load_balancing_method(options_sim.lb_name.c_str())
		.set_maximum_refinement_level(0)
		.initialize(comm);

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

	try {
		grid.set_geometry(geom_params);
	} catch (...) {
		cerr << __FILE__ << ":" << __LINE__
			<< ": Couldn't set grid geometry."
			<< endl;
		abort();
	}

	for (const auto& cell: grid.local_cells()) {
		(*cell.data)[pamhd::MPI_Rank()] = rank;
		Substep(*cell.data) = 1;
		Max_v(*cell.data) = {-1, -1, -1, -1, -1, -1};
	}
	Cell::set_transfer_all(true,
		pamhd::mhd::Max_Velocity(),
		pamhd::mhd::Substepping_Period(),
		pamhd::mhd::Substep_Min(),
		pamhd::mhd::Substep_Max(),
		pamhd::MPI_Rank()
	);
	grid.update_copies_of_remote_neighbors();
	Cell::set_transfer_all(false,
		pamhd::mhd::Max_Velocity(),
		pamhd::mhd::Substepping_Period(),
		pamhd::mhd::Substep_Min(),
		pamhd::mhd::Substep_Max(),
		pamhd::MPI_Rank()
	);

	// assign cells into boundary geometries
	for (const auto& gid: geometries.get_geometry_ids()) {
		geometries.clear_cells(gid);
	}
	for (const auto& cell: grid.local_cells()) {
		geometries.overlaps(
			grid.geometry.get_min(cell.id),
			grid.geometry.get_max(cell.id),
			cell.id);
	}

	/*
	Simulate
	*/

	if (rank == 0) {
		cout << "Initializing particles and magnetic field... " << endl;
	}

	double simulation_time = options_sim.time_start;
	pamhd::mhd::initialize_magnetic_field_staggered<pamhd::Magnetic_Field>(
		geometries,
		initial_conditions_fields,
		background_B,
		grid,
		simulation_time,
		options_sim.vacuum_permeability,
		Face_B, Mag_fs, Bg_B
	);
	Cell::set_transfer_all(true,
		pamhd::Face_Magnetic_Field(),
		pamhd::Bg_Magnetic_Field()
	);
	grid.update_copies_of_remote_neighbors();
	Cell::set_transfer_all(false,
		pamhd::Face_Magnetic_Field(),
		pamhd::Bg_Magnetic_Field()
	);

	pamhd::mhd::update_B_consistency(
		0, grid.local_cells(),
		Mas, Mom, Nrj, Vol_B, Face_B,
		Sol_Info2, Substep,
		options_sim.adiabatic_index,
		options_sim.vacuum_permeability,
		false // fluid not initialized yet
	);
	Cell::set_transfer_all(true,
		pamhd::Face_Magnetic_Field(),
		pamhd::mhd::MHD_State_Conservative()
	);
	grid.update_copies_of_remote_neighbors();
	Cell::set_transfer_all(false, pamhd::Face_Magnetic_Field());

/*	pamhd::mhd::initialize_fluid_staggered(
		geometries,
		initial_conditions_fields,
		grid,
		simulation_time,
		options_sim.adiabatic_index,
		options_sim.vacuum_permeability,
		options_sim.proton_mass,
		true,
		Mas, Mom, Nrj, Vol_B,
		Mas_fs, Mom_fs, Nrj_fs
	);
	grid.update_copies_of_remote_neighbors();
	Cell::set_transfer_all(false, pamhd::mhd::MHD_State_Conservative());*/

/*	pamhd::mhd::apply_magnetic_field_boundaries(
		grid,
		boundaries_fields,
		geometries,
		simulation_time,
		Mag
	);*/

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
			Mas, Mom, Nrj, Vol_B,
			Sol_Info
		);
	} catch (const std::exception& e) {
		std::cerr << __FILE__ "(" << __LINE__ << ": "
			<< "Couldn't fill MHD fluid values: " << e.what()
			<< std::endl;
		abort();
	}

	if (rank == 0) {
		cout << "done initializing particles and fields." << endl;
	}

	/*
	Classify cells into normal, boundary and dont_solve
	*/

	Cell::set_transfer_all(true, pamhd::particle::Solver_Info());
	pamhd::mhd::set_solver_info_magnetic<pamhd::particle::Solver_Info>(
		grid, boundaries_fields, geometries, Sol_Info
	);
	pamhd::particle::set_solver_info<pamhd::particle::Solver_Info>(
		grid, boundaries_particles, geometries, Sol_Info
	);
	Cell::set_transfer_all(false, pamhd::particle::Solver_Info());

	Cell::set_transfer_all(true, pamhd::mhd::Face_Boundary_Type());
	pamhd::mhd::classify_faces(grid, Sol_Info2, FInfo);
	grid.update_copies_of_remote_neighbors();
	Cell::set_transfer_all(false, pamhd::mhd::Face_Boundary_Type());

	/*pamhd::mhd::apply_fluid_boundaries(
		grid,
		boundaries,
		geometries,
		simulation_time,
		Mas, Mom, Nrj, Vol_B, Sol_Info2,
		options_sim.proton_mass,
		options_sim.adiabatic_index,
		options_sim.vacuum_permeability
	);*/
	Cell::set_transfer_all(true, pamhd::mhd::MHD_State_Conservative());
	grid.update_copies_of_remote_neighbors();
	Cell::set_transfer_all(false, pamhd::mhd::MHD_State_Conservative());

	pamhd::mhd::apply_magnetic_field_boundaries_staggered(
		grid,
		boundaries_fields,
		geometries,
		simulation_time,
		Face_B, FInfo
	);
	Cell::set_transfer_all(true, pamhd::Face_Magnetic_Field());
	grid.update_copies_of_remote_neighbors();
	Cell::set_transfer_all(false, pamhd::Face_Magnetic_Field());

	pamhd::mhd::update_B_consistency(
		0, grid.local_cells(),
		Mas, Mom, Nrj, Vol_B, Face_B,
		Sol_Info2, Substep,
		options_sim.adiabatic_index,
		options_sim.vacuum_permeability,
		true
	);
	Cell::set_transfer_all(true,
		pamhd::Face_Magnetic_Field(),
		pamhd::mhd::MHD_State_Conservative()
	);
	grid.update_copies_of_remote_neighbors();
	Cell::set_transfer_all(false,
		pamhd::Face_Magnetic_Field(),
		pamhd::mhd::MHD_State_Conservative()
	);

	for (const auto& cell: grid.local_cells()) {
		Max_v(*cell.data) = {-1, -1, -1, -1, -1, -1};
		Timestep(*cell.data) = 0;
		Substep(*cell.data)     =
		Substep_Max(*cell.data) =
		Substep_Min(*cell.data) = 1;
	}
	Grid::cell_data_type::set_transfer_all(true,
		pamhd::mhd::Timestep(),
		pamhd::mhd::Max_Velocity(),
		pamhd::mhd::Substepping_Period()
	);
	grid.update_copies_of_remote_neighbors();
	Grid::cell_data_type::set_transfer_all(false,
		pamhd::mhd::Timestep(),
		pamhd::mhd::Max_Velocity(),
		pamhd::mhd::Substepping_Period()
	);

	if (rank == 0) {
		cout << "Done initializing" << endl;
	}

	const double time_end = options_sim.time_start + options_sim.time_length;
	double
		max_dt_mhd = 0,
		max_dt_particle_gyro = 0,
		max_dt_particle_flight = 0,
		next_particle_save = 0,
		next_mhd_save = 0,
		next_amr = options_grid.amr_n;
	size_t simulation_step = 0;
	while (simulation_time < time_end) {
		simulation_step++;

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
			cerr << __FILE__ << ":" << __LINE__
				<< ": Couldn't reduce time step."
				<< endl;
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
			cerr << __FILE__ "(" << __LINE__ << ": "
				<< "Couldn't accumulate MHD data from particles: " << e.what()
				<< endl;
			abort();
		}

		// B required for E calculation
		Cell::set_transfer_all(true, pamhd::mhd::MHD_State_Conservative());
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
				Mas, Mom, Nrj, Vol_B,
				Sol_Info
			);
		} catch (const std::exception& e) {
			cerr << __FILE__ "(" << __LINE__ << ": "
				<< "Couldn't fill MHD fluid values: " << e.what()
				<< endl;
			abort();
		}

		// inner: J for E = (J - V) x B
		pamhd::divergence::get_curl(
			grid.inner_cells(),
			grid,
			Vol_B,
			Vol_J,
			Sol_Info
		);
		// not included in get_curl above
		for (const auto& cell: grid.inner_cells()) {
			Vol_J(*cell.data) /= options_sim.vacuum_permeability;
		}

		grid.wait_remote_neighbor_copy_update_receives();

		// outer: J for E = (J - V) x B
		pamhd::divergence::get_curl(
			grid.outer_cells(),
			grid,
			Vol_B,
			Vol_J,
			Sol_Info
		);
		for (const auto& cell: grid.outer_cells()) {
			Vol_J(*cell.data) /= options_sim.vacuum_permeability;
		}

		grid.wait_remote_neighbor_copy_update_sends();
		Cell::set_transfer_all(false, pamhd::mhd::MHD_State_Conservative());

		// inner: E = (J - V) x B
		for (const auto& cell: grid.inner_cells()) {
			J_m_V(*cell.data) = Vol_J(*cell.data) - pamhd::mhd::get_velocity(Mom(*cell.data), Mas(*cell.data));
			// calculate electric field for output file
			Vol_E(*cell.data) = J_m_V(*cell.data).cross(Vol_B(*cell.data));
		}

		// outer: E = (J - V) x B
		for (const auto& cell: grid.outer_cells()) {
			J_m_V(*cell.data) = Vol_J(*cell.data) - pamhd::mhd::get_velocity(Mom(*cell.data), Mas(*cell.data));
			Vol_E(*cell.data) = J_m_V(*cell.data).cross(Vol_B(*cell.data));
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

		if (grid.get_rank() == 0) {
			cout << "Solving at time " << simulation_time
				<< " s with time step " << time_step << " s" << flush;
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
				Vol_B,\
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

		// outer particles
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
			cerr <<  __FILE__ << "(" << __LINE__ << "): " << particle_stepper << endl;
			abort();
		}
		max_dt_particle_flight = min(particle_max_dt.first, max_dt_particle_flight);
		max_dt_particle_gyro = min(particle_max_dt.second, max_dt_particle_gyro);

		Cell::set_transfer_all(
			true,
			pamhd::mhd::MHD_State_Conservative(),
			pamhd::particle::Nr_Particles_External()
		);
		grid.start_remote_neighbor_copy_updates();

		// inner MHD
		pamhd::mhd::get_fluxes(
			mhd_solver,
			grid.inner_cells(),
			grid,
			1,
			options_sim.adiabatic_index,
			options_sim.vacuum_permeability,
			time_step,
			Mas, Mom, Nrj, Vol_B, Face_dB, Bg_B,
			Mas_fs, Mom_fs, Nrj_fs, Mag_fs,
			Sol_Info2, Substep, Max_v
		);

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
			cerr <<  __FILE__ << "(" << __LINE__ << "): " << particle_stepper << endl;
			abort();
		}
		#undef SOLVE_WITH_STEPPER
		max_dt_particle_flight = min(particle_max_dt.first, max_dt_particle_flight);
		max_dt_particle_gyro = min(particle_max_dt.second, max_dt_particle_gyro);

		grid.wait_remote_neighbor_copy_update_receives();

		// outer MHD
		pamhd::mhd::get_fluxes(
			mhd_solver,
			grid.outer_cells(),
			grid,
			1,
			options_sim.adiabatic_index,
			options_sim.vacuum_permeability,
			time_step,
			Mas, Mom, Nrj, Vol_B, Face_dB, Bg_B,
			Mas_fs, Mom_fs, Nrj_fs, Mag_fs,
			Sol_Info2, Substep, Max_v
		);

		pamhd::particle::resize_receiving_containers<
			pamhd::particle::Nr_Particles_External,
			pamhd::particle::Particles_External
		>(grid.remote_cells(), grid);

		grid.wait_remote_neighbor_copy_update_sends();
		Cell::set_transfer_all(
			false,
			pamhd::mhd::MHD_State_Conservative(),
			pamhd::particle::Nr_Particles_External()
		);

		Cell::set_transfer_all(true, pamhd::particle::Particles_External());
		grid.start_remote_neighbor_copy_updates();

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


		Cell::set_transfer_all(true, pamhd::mhd::MHD_Flux());
		grid.update_copies_of_remote_neighbors();
		Cell::set_transfer_all(false, pamhd::mhd::MHD_Flux());

		Cell::set_transfer_all(true,
			// update pressure for B consistency calculation
			pamhd::mhd::MHD_State_Conservative()
		);
		grid.update_copies_of_remote_neighbors();
		Cell::set_transfer_all(false,
			pamhd::mhd::MHD_State_Conservative()
		);

		// constant thermal pressure when updating vol B after solution
		pamhd::mhd::update_B_consistency(
			0, grid.local_cells(),
			Mas, Mom, Nrj, Vol_B, Face_B,
			Sol_Info2, Substep,
			options_sim.adiabatic_index,
			options_sim.vacuum_permeability,
			true
		);
		Cell::set_transfer_all(true,
			pamhd::Face_Magnetic_Field(),
			pamhd::mhd::MHD_State_Conservative()
		);
		grid.update_copies_of_remote_neighbors();
		Cell::set_transfer_all(false,
			pamhd::Face_Magnetic_Field(),
			pamhd::mhd::MHD_State_Conservative()
		);

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
			pamhd::mhd::MHD_State_Conservative(),
			pamhd::particle::Nr_Particles_External()
		);
		grid.update_copies_of_remote_neighbors();
		Cell::set_transfer_all(
			false,
			pamhd::mhd::MHD_State_Conservative(),
			pamhd::particle::Nr_Particles_External()
		);

		pamhd::particle::resize_receiving_containers<
			pamhd::particle::Nr_Particles_External,
			pamhd::particle::Particles_Internal
		>(grid.remote_cells(), grid);
		Cell::set_transfer_all(true, pamhd::particle::Particles_Internal());
		grid.update_copies_of_remote_neighbors();
		Cell::set_transfer_all(false, pamhd::particle::Particles_Internal());

		pamhd::mhd::apply_magnetic_field_boundaries_staggered(
			grid,
			boundaries_fields,
			geometries,
			simulation_time,
			Face_B, FInfo
		);
		Cell::set_transfer_all(true, pamhd::Face_Magnetic_Field());
		grid.update_copies_of_remote_neighbors();
		Cell::set_transfer_all(false, pamhd::Face_Magnetic_Field());

		pamhd::mhd::update_B_consistency(
			0, grid.local_cells(),
			Mas, Mom, Nrj, Vol_B, Face_B,
			Sol_Info2, Substep,
			options_sim.adiabatic_index,
			options_sim.vacuum_permeability,
			true
		);
		Cell::set_transfer_all(true,
			pamhd::Face_Magnetic_Field(),
			pamhd::mhd::MHD_State_Conservative()
		);
		grid.update_copies_of_remote_neighbors();
		Cell::set_transfer_all(false,
			pamhd::Face_Magnetic_Field(),
			pamhd::mhd::MHD_State_Conservative()
		);

		const auto avg_div = pamhd::math::get_divergence_staggered(
			grid.local_cells(), grid,
			Face_B, Mag_div,
			Sol_Info2
		);
		if (rank == 0) {
			cout << " average divergence " << avg_div << endl;
		}
		Cell::set_transfer_all(true, pamhd::Magnetic_Field_Divergence());
		grid.update_copies_of_remote_neighbors();
		Cell::set_transfer_all(false, pamhd::Magnetic_Field_Divergence());

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
				simulation_step,
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

		/*
		Save simulation to disk
		*/

		// particles
		if (
			(options_particle.save_n >= 0 and simulation_time >= time_end)
			or (options_particle.save_n > 0 and simulation_time >= next_particle_save)
		) {
			if (next_particle_save <= simulation_time) {
				next_particle_save
					+= options_particle.save_n
					* ceil(max(options_particle.save_n, simulation_time - next_particle_save) / options_particle.save_n);
			}

			if (rank == 0) {
				cout << "Saving particles at time " << simulation_time << "... " << endl;
			}

			// update number of internal particles
			for (const auto& cell: grid.local_cells()) {
				Nr_Int(*cell.data) = Part_Int(*cell.data).size();
			}
			constexpr uint64_t file_version = 4;
			if (
				not pamhd::particle::save(
					boost::filesystem::canonical(
						boost::filesystem::path(options_sim.output_directory)
					).append("particle_").generic_string(),
					grid,
					file_version,
					simulation_step,
					simulation_time,
					options_sim.adiabatic_index,
					options_sim.proton_mass,
					options_particle.boltzmann
				)
			) {
				cerr <<  __FILE__ << "(" << __LINE__ << "): "
					"Couldn't save particle result."
					<< endl;
				abort();
			}
		}

		// mhd
		if (
			(options_mhd.save_n >= 0 and simulation_time >= time_end)
			or (options_mhd.save_n > 0 and simulation_time >= next_mhd_save)
		) {
			if (next_mhd_save <= simulation_time) {
				next_mhd_save
					+= options_mhd.save_n
					* ceil(max(options_mhd.save_n, simulation_time - next_mhd_save) / options_mhd.save_n);
			}

			if (rank == 0) {
				cout << "Saving MHD at time " << simulation_time << "... " << endl;
			}

			constexpr uint64_t file_version = 4;
			if (
				not pamhd::mhd::save_staggered(
					boost::filesystem::canonical(
						boost::filesystem::path(options_sim.output_directory)
					).append("mhd_staggered_").generic_string(),
					grid,
					file_version,
					simulation_step,
					simulation_time,
					options_sim.adiabatic_index,
					options_sim.proton_mass,
					options_sim.vacuum_permeability
				)
			) {
				cerr <<  __FILE__ << "(" << __LINE__ << "): "
					"Couldn't save mhd result."
					<< endl;
				abort();
			}
		}
	}

	if (rank == 0) {
		cout << "Simulation finished at time " << simulation_time << endl;
	}
	MPI_Finalize();

	return EXIT_SUCCESS;
}
