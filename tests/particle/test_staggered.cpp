/*
Hybrid PIC program of PAMHD.

Copyright 2015, 2016, 2017 Ilja Honkonen
Copyright 2019, 2024, 2025 Finnish Meteorological Institute
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
#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#include "Eigen/Core" // must be included before gensimcell.hpp
#include "mpi.h" // must be included before gensimcell.hpp
#include "gensimcell.hpp"
#include "rapidjson/document.h"
#include "rapidjson/error/en.h"

#include "background_magnetic_field.hpp"
#include "boundaries/geometries.hpp"
#include "boundaries/multivariable_boundaries.hpp"
#include "boundaries/multivariable_initial_conditions.hpp"
#include "grid/options.hpp"
#include "grid/variables.hpp"
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
#include "variable_getter.hpp"
#include "variables.hpp"


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

const auto Bg_B = pamhd::Variable_Getter<pamhd::Bg_Magnetic_Field>();
bool pamhd::Bg_Magnetic_Field::is_stale = true;

const auto Mas = pamhd::Variable_Getter<pamhd::mhd::Mass_Density>();
bool pamhd::mhd::Mass_Density::is_stale = true;

const auto Mom = pamhd::Variable_Getter<pamhd::mhd::Momentum_Density>();
bool pamhd::mhd::Momentum_Density::is_stale = true;

const auto Nrj = pamhd::Variable_Getter<pamhd::mhd::Total_Energy_Density>();
bool pamhd::mhd::Total_Energy_Density::is_stale = true;

const auto Vol_B = pamhd::Variable_Getter<pamhd::Magnetic_Field>();
bool pamhd::Magnetic_Field::is_stale = true;

const auto Face_B = pamhd::Variable_Getter<pamhd::Face_Magnetic_Field>();
bool pamhd::Face_Magnetic_Field::is_stale = true;

const auto Face_dB = pamhd::Variable_Getter<pamhd::Face_dB>();

const auto Div_B = pamhd::Variable_Getter<pamhd::Magnetic_Field_Divergence>();

// curl of magnetic field
const auto Vol_J = pamhd::Variable_Getter<pamhd::Electric_Current_Density>();
bool pamhd::Electric_Current_Density::is_stale = true;

// electric current minus bulk velocity
const auto J_m_V = pamhd::Variable_Getter<pamhd::particle::Current_Minus_Velocity>();
bool pamhd::particle::Current_Minus_Velocity::is_stale = true;

// electric field for propagating particles
const auto Vol_E = pamhd::Variable_Getter<pamhd::particle::Electric_Field>();
bool pamhd::particle::Electric_Field::is_stale = true;

/*! Solver info variable for boundary logic

-1 for cells not to be read nor written
0 for read-only cells
1 for read-write cells
*/
const auto SInfo = pamhd::Variable_Getter<pamhd::Solver_Info>();
bool pamhd::Solver_Info::is_stale = true;

const auto Max_v_part = pamhd::Variable_Getter<pamhd::particle::Max_Spatial_Velocity>();
bool pamhd::particle::Max_Spatial_Velocity::is_stale = true;

const auto Max_ω_part = pamhd::Variable_Getter<pamhd::particle::Max_Angular_Velocity>();
bool pamhd::particle::Max_Angular_Velocity::is_stale = true;

const auto Substep = pamhd::Variable_Getter<pamhd::mhd::Substepping_Period>();
bool pamhd::mhd::Substepping_Period::is_stale = true;

const auto Substep_Min = pamhd::Variable_Getter<pamhd::mhd::Substep_Min>();
bool pamhd::mhd::Substep_Min::is_stale = true;

const auto Substep_Max = pamhd::Variable_Getter<pamhd::mhd::Substep_Max>();
bool pamhd::mhd::Substep_Max::is_stale = true;

const auto Timestep = pamhd::Variable_Getter<pamhd::mhd::Timestep>();
bool pamhd::mhd::Timestep::is_stale = true;

const auto Max_v_wave = pamhd::Variable_Getter<pamhd::mhd::Max_Velocity>();
bool pamhd::mhd::Max_Velocity::is_stale = true;

const auto MHDF = pamhd::Variable_Getter<pamhd::mhd::MHD_Flux>();
const auto Mas_f = [](Cell& cell_data, const int dir)->auto& {
	return MHDF.data(cell_data)(dir)[pamhd::mhd::Mass_Density()];
};
const auto Mom_f = [](Cell& cell_data, const int dir)->auto& {
	return MHDF.data(cell_data)(dir)[pamhd::mhd::Momentum_Density()];
};
const auto Nrj_f = [](Cell& cell_data, const int dir)->auto& {
	return MHDF.data(cell_data)(dir)[pamhd::mhd::Total_Energy_Density()];
};
const auto Mag_f = [](Cell& cell_data, const int dir)->auto& {
	return MHDF.data(cell_data)(dir)[pamhd::Magnetic_Field()];
};

const auto FInfo = pamhd::Variable_Getter<pamhd::mhd::Face_Boundary_Type>();
bool pamhd::mhd::Face_Boundary_Type::is_stale = true;

const auto Ref_min = pamhd::Variable_Getter<pamhd::grid::Target_Refinement_Level_Min>();
bool pamhd::grid::Target_Refinement_Level_Min::is_stale = true;

const auto Ref_max = pamhd::Variable_Getter<pamhd::grid::Target_Refinement_Level_Max>();
bool pamhd::grid::Target_Refinement_Level_Max::is_stale = true;

// list of particles in cell not moving to another cell
const auto Part_Int = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::particle::Particles_Internal()];
};

// particles moving to another cell
const auto Part_Ext = pamhd::Variable_Getter<pamhd::particle::Particles_External>();
bool pamhd::particle::Particles_External::is_stale = true;

// number of particles in above list, for allocating memory for arriving particles
const auto Nr_Ext = pamhd::Variable_Getter<pamhd::particle::Nr_Particles_External>();
bool pamhd::particle::Nr_Particles_External::is_stale = true;

const auto Nr_Int = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::particle::Nr_Particles_Internal()];
};
// references to initial condition & boundary data of cell
const auto Bdy_N = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::particle::Bdy_Number_Density()];
};
const auto Bdy_V = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::particle::Bdy_Velocity()];
};
const auto Bdy_T = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::particle::Bdy_Temperature()];
};
const auto Bdy_Nr_Par = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::particle::Bdy_Nr_Particles_In_Cell()];
};
const auto Bdy_SpM = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::particle::Bdy_Species_Mass()];
};
const auto Bdy_C2M = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::particle::Bdy_Charge_Mass_Ratio()];
};

// given a particle these return references to particle's parameters
const auto Part_Pos = [](
	pamhd::particle::Particle_Internal& particle
)->auto& {
	return particle[pamhd::particle::Position()];
};
// as above but for caller that also provides cell's data
const auto Part_Vel_Cell = [](
	Cell&, pamhd::particle::Particle_Internal& particle
)->auto& {
	return particle[pamhd::particle::Velocity()];
};
const auto Part_Vel = [](
	pamhd::particle::Particle_Internal& particle
)->auto& {
	return particle[pamhd::particle::Velocity()];
};
const auto Part_C2M = [](
	pamhd::particle::Particle_Internal& particle
)->auto& {
	return particle[pamhd::particle::Charge_Mass_Ratio()];
};
// copy of number of real particles represented by simulation particle
const auto Part_Nr = [](
	pamhd::particle::Particle_Internal& particle
)->auto {
	return particle[pamhd::particle::Mass()] / particle[pamhd::particle::Species_Mass()];
};
const auto Part_Mas = [](
	pamhd::particle::Particle_Internal& particle
)->auto& {
	return particle[pamhd::particle::Mass()];
};
// as above but for caller that also provides cell's data
const auto Part_Mas_Cell = [](
	Cell&, pamhd::particle::Particle_Internal& particle
)->auto& {
	return particle[pamhd::particle::Mass()];
};
const auto Part_Des = [](
	pamhd::particle::Particle_External& particle
)->auto& {
	return particle[pamhd::particle::Destination_Cell()];
};
// reference to mass of given particle's species
const auto Part_SpM = [](
	pamhd::particle::Particle_Internal& particle
)->auto& {
	return particle[pamhd::particle::Species_Mass()];
};
// as above but for caller that also provides cell's data
const auto Part_SpM_Cell = [](
	Cell&, pamhd::particle::Particle_Internal& particle
)->auto& {
	return particle[pamhd::particle::Species_Mass()];
};
// copy of particle's kinetic energy relative to pamhd::particle::Bulk_Velocity
const auto Part_Ekin = [](
	Cell& cell_data,
	pamhd::particle::Particle_Internal& particle
)->auto {
	return
		0.5 * particle[pamhd::particle::Mass()]
		* (
			particle[pamhd::particle::Velocity()]
			- cell_data[pamhd::particle::Bulk_Velocity()].first
		).squaredNorm();
};

// reference to accumulated number of particles in given cell
const auto Nr_Particles = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::particle::Number_Of_Particles()];
};

const auto Bulk_Mass_Getter = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::particle::Bulk_Mass()];
};

const auto Bulk_Momentum_Getter = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::particle::Bulk_Momentum()];
};

const auto Bulk_Relative_Velocity2_Getter = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::particle::Bulk_Relative_Velocity2()];
};

const auto Bulk_Velocity_Getter = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::particle::Bulk_Velocity()];
};

// list of items (variables above) accumulated from particles in given cell
const auto Accu_List_Getter = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::particle::Accumulated_To_Cells()];
};

// length of above list (for transferring between processes)
const auto Accu_List_Length_Getter = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::particle::Nr_Accumulated_To_Cells()];
};

// target cell of accumulated particle values in an accumulation list item
const auto Accu_List_Target_Getter = [](
	pamhd::particle::Accumulated_To_Cell& accu_item
)->auto& {
	return accu_item[pamhd::particle::Target()];
};

// accumulated number of particles in an accumulation list item
const auto Accu_List_Number_Of_Particles_Getter = [](
	pamhd::particle::Accumulated_To_Cell& accu_item
)->auto& {
	return accu_item[pamhd::particle::Number_Of_Particles()];
};

const auto Accu_List_Bulk_Mass_Getter = [](
	pamhd::particle::Accumulated_To_Cell& accu_item
)->auto& {
	return accu_item[pamhd::particle::Bulk_Mass()];
};

const auto Accu_List_Bulk_Velocity_Getter = [](
	pamhd::particle::Accumulated_To_Cell& accu_item
)->auto& {
	return accu_item[pamhd::particle::Bulk_Velocity()];
};

const auto Accu_List_Bulk_Relative_Velocity2_Getter = [](
	pamhd::particle::Accumulated_To_Cell& accu_item
)->auto& {
	return accu_item[pamhd::particle::Bulk_Relative_Velocity2()];
};


int main(int argc, char* argv[]) {
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

	// initialize Zoltan
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
		if (rank == 0) {
			cerr << "Couldn't parse json data in file "
				<< argv[1] << " at character position "
				<< document.GetErrorOffset() << ": "
				<< rapidjson::GetParseError_En(document.GetParseError())
				<< endl;
		}
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
		pamhd::mhd::Number_Density,
		pamhd::mhd::Velocity,
		pamhd::mhd::Pressure,
		pamhd::Magnetic_Field
	> initial_conditions_mhd;
	initial_conditions_mhd.set(document);

	pamhd::boundaries::Multivariable_Boundaries<
		uint64_t,
		geometry_id_t,
		pamhd::mhd::Number_Density,
		pamhd::mhd::Velocity,
		pamhd::mhd::Pressure,
		pamhd::Magnetic_Field
	> boundaries_mhd;
	boundaries_mhd.set(document);

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
		Substep.data(*cell.data) = 1;
		Max_v_wave.data(*cell.data) = {-1, -1, -1, -1, -1, -1};
	}
	pamhd::mhd::set_minmax_substepping_period(
		options_sim.time_start, grid,
		options_mhd, Substep_Min, Substep_Max);
	Cell::set_transfer_all(true,
		Max_v_wave.type(), Substep.type(), Substep_Min.type(),
		Substep_Max.type(), pamhd::MPI_Rank());
	grid.update_copies_of_remote_neighbors();
	Cell::set_transfer_all(false,
		Max_v_wave.type(), Substep.type(), Substep_Min.type(),
		Substep_Max.type(), pamhd::MPI_Rank());

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

	const double time_end = options_sim.time_start + options_sim.time_length;
	double
		simulation_time = options_sim.time_start,
		next_particle_save = 0,
		next_mhd_save = 0,
		next_amr = options_grid.amr_n;

	if (rank == 0) {
		cout << "Initializing... " << endl;
	}
	pamhd::mhd::initialize_magnetic_field_staggered<pamhd::Magnetic_Field>(
		geometries, initial_conditions_mhd, background_B,
		grid, simulation_time, options_sim.vacuum_permeability,
		Face_B, Mag_f, Bg_B
	);

	pamhd::mhd::update_B_consistency(
		0, grid.local_cells(), grid,
		Mas, Mom, Nrj, Vol_B, Face_B,
		SInfo, Substep,
		options_sim.adiabatic_index,
		options_sim.vacuum_permeability,
		false // fluid not initialized yet
	);

	pamhd::mhd::initialize_fluid_staggered(
		geometries, initial_conditions_mhd,
		grid, simulation_time,
		options_sim.adiabatic_index,
		options_sim.vacuum_permeability,
		options_sim.proton_mass, true,
		Mas, Mom, Nrj, Vol_B,
		Mas_f, Mom_f, Nrj_f
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
			options_sim.temp2nrj,
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
			SInfo
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
			options_sim.temp2nrj,
			options_sim.vacuum_permeability,
			next_particle_id,
			grid.get_comm_size(),
			true,
			SInfo,
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
			SInfo
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
			options_sim.temp2nrj,
			options_mhd.min_pressure,
			Nr_Particles,
			Bulk_Mass_Getter,
			Bulk_Momentum_Getter,
			Bulk_Relative_Velocity2_Getter,
			Part_Int,
			Mas, Mom, Nrj, Vol_B,
			SInfo
		);
	} catch (const std::exception& e) {
		std::cerr << __FILE__ "(" << __LINE__ << ": "
			<< "Couldn't fill MHD fluid values: " << e.what()
			<< std::endl;
		abort();
	}

	/*
	Classify cells & faces into normal, boundary and dont_solve
	*/

	pamhd::mhd::set_solver_info(grid, boundaries_mhd, geometries, SInfo);
	pamhd::mhd::classify_faces(grid, SInfo, FInfo);

	/*pamhd::mhd::apply_fluid_boundaries(
		grid,
		boundaries,
		geometries,
		simulation_time,
		Mas, Mom, Nrj, Vol_B, SInfo,
		options_sim.proton_mass,
		options_sim.adiabatic_index,
		options_sim.vacuum_permeability
	);*/

	pamhd::mhd::apply_magnetic_field_boundaries_staggered(
		grid,
		boundaries_mhd,
		geometries,
		simulation_time,
		Face_B, FInfo
	);

	pamhd::mhd::update_B_consistency(
		0, grid.local_cells(), grid,
		Mas, Mom, Nrj, Vol_B, Face_B,
		SInfo, Substep,
		options_sim.adiabatic_index,
		options_sim.vacuum_permeability,
		true
	);

	for (const auto& cell: grid.local_cells()) {
		Timestep.data(*cell.data) = -1;
		Substep.data(*cell.data)     =
		Substep_Max.data(*cell.data) =
		Substep_Min.data(*cell.data) = 1;
		Max_v_wave.data(*cell.data) = {-1, -1, -1, -1, -1, -1};
		Max_v_part.data(*cell.data) =
		Max_ω_part.data(*cell.data) = -1;
	}
	Cell::set_transfer_all(true,
		Timestep.type(), Max_v_wave.type(), Substep.type());
	grid.update_copies_of_remote_neighbors();
	Cell::set_transfer_all(false,
		Timestep.type(), Max_v_wave.type(), Substep.type());

	// final init with timestep of 0
	pamhd::particle::timestep(
		options_sim.time_start, grid, options_mhd, Part_Int,
		Part_Pos, Part_Mas, Part_Mas_Cell, Part_SpM,
		Part_SpM_Cell, Part_Vel, Part_Vel_Cell, Part_Ekin,
		Nr_Particles, Part_Nr, Bulk_Mass_Getter,
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
		pamhd::particle::Bulk_Velocity(), SInfo,
		options_sim.adiabatic_index,
		options_sim.vacuum_permeability,
		options_sim.temp2nrj,
		options_mhd.min_pressure,
		Mas, Mom, Nrj, Vol_B, Vol_J, J_m_V, Vol_E, Nr_Ext,
		Max_v_part, Max_ω_part, Part_Ext, Part_C2M, Part_Des,
		Face_dB, Bg_B, Mas_f, Mom_f, Nrj_f, Mag_f, Substep,
		Substep_Min, Substep_Max, Max_v_wave, Face_B, background_B,
		mhd_solver, Timestep, 0, options_mhd.time_step_factor
	);
	if (rank == 0) {
		cout << "done" << endl;
	}

	size_t simulation_step = 0;
	constexpr uint64_t file_version = 4;
	if (options_particle.save_n >= 0) {
		if (rank == 0) {
			cout << "Saving particles at time " << simulation_time << endl;
		}

		// update number of internal particles
		for (const auto& cell: grid.local_cells()) {
			Nr_Int(*cell.data) = Part_Int(*cell.data).size();
		}
		if (
			not pamhd::particle::save(
				boost::filesystem::canonical(
					boost::filesystem::path(options_sim.output_directory)
				).append("particle_").generic_string(),
				grid, file_version,
				simulation_step, simulation_time,
				options_sim.adiabatic_index,
				options_sim.proton_mass,
				options_sim.temp2nrj
			)
		) {
			cerr <<  __FILE__ << "(" << __LINE__ << "): "
				"Couldn't save particle result."
				<< endl;
			abort();
		}
	}

	if (options_mhd.save_n >= 0) {
		if (rank == 0) {
			cout << "Saving MHD at time " << simulation_time << endl;
		}

		if (
			not pamhd::mhd::save_staggered(
				boost::filesystem::canonical(
					boost::filesystem::path(options_sim.output_directory)
				).append("mhd_staggered_").generic_string(),
				grid, file_version,
				simulation_step, simulation_time,
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

	while (simulation_time < time_end) {
		simulation_step++;

		pamhd::particle::split_particles(
			options_particle.min_particles, random_source,
			grid, Part_Int, Part_Pos, Part_Mas, SInfo
		);

		// don't step over the final simulation time
		const double
			until_end = time_end - simulation_time,
			dt = pamhd::particle::timestep(
				simulation_time, grid, options_mhd, Part_Int,
				Part_Pos, Part_Mas, Part_Mas_Cell, Part_SpM,
				Part_SpM_Cell, Part_Vel, Part_Vel_Cell, Part_Ekin,
				Nr_Particles, Part_Nr, Bulk_Mass_Getter,
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
				pamhd::particle::Bulk_Velocity(), SInfo,
				options_sim.adiabatic_index,
				options_sim.vacuum_permeability,
				options_sim.temp2nrj,
				options_mhd.min_pressure,
				Mas, Mom, Nrj, Vol_B, Vol_J, J_m_V, Vol_E, Nr_Ext,
				Max_v_part, Max_ω_part, Part_Ext, Part_C2M, Part_Des,
				Face_dB, Bg_B, Mas_f, Mom_f, Nrj_f, Mag_f, Substep,
				Substep_Min, Substep_Max, Max_v_wave, Face_B, background_B,
				mhd_solver, Timestep, until_end, options_mhd.time_step_factor
			);

		if (rank == 0) {
			cout << "Solution calculated at time " << simulation_time
				<< " s with time step " << dt << " s" << flush;
		}

		simulation_time += dt;

		/*
		Update internal particles for setting particle copy boundaries.

		TODO overlap computation and communication in boundary processing
		*/
		for (const auto& cell: grid.local_cells()) {
			// (ab)use external number counter as internal number counter
			(*cell.data)[pamhd::particle::Nr_Particles_External()]
				= (*cell.data)[pamhd::particle::Particles_Internal()].size();
		}
		Cell::set_transfer_all(true,
			pamhd::particle::Nr_Particles_External()
		);
		grid.update_copies_of_remote_neighbors();
		Cell::set_transfer_all(false,
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
			boundaries_mhd,
			geometries,
			simulation_time,
			Face_B, FInfo
		);

		pamhd::mhd::update_B_consistency(
			0, grid.local_cells(), grid,
			Mas, Mom, Nrj, Vol_B, Face_B,
			SInfo, Substep,
			options_sim.adiabatic_index,
			options_sim.vacuum_permeability,
			true
		);

		const auto avg_div = pamhd::math::get_divergence_staggered(
			grid.local_cells(), grid,
			Face_B, Div_B, SInfo
		);
		if (rank == 0) {
			cout << " average divergence " << avg_div << endl;
		}

		nr_particles_created += pamhd::particle::apply_boundaries<
			pamhd::particle::Particle_Internal,
			pamhd::particle::Mass,
			pamhd::particle::Charge_Mass_Ratio,
			pamhd::particle::Position,
			pamhd::particle::Velocity,
			pamhd::particle::Particle_ID,
			pamhd::particle::Species_Mass
		>(
			geometries, boundaries_particles, simulation_time,
			simulation_step, grid, random_source,
			options_sim.temp2nrj,
			options_sim.vacuum_permeability,
			next_particle_id, grid.get_comm_size(), true,
			SInfo, Part_Int, Bdy_N, Bdy_V, Bdy_T,
			Bdy_Nr_Par, Bdy_SpM, Bdy_C2M
		);
		next_particle_id += nr_particles_created * grid.get_comm_size();

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
					options_sim.temp2nrj
				)
			) {
				cerr <<  __FILE__ << "(" << __LINE__ << "): "
					"Couldn't save particle result."
					<< endl;
				abort();
			}
		}

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
