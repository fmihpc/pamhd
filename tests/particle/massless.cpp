/*
Propagates test particles (0 mass) in prescribed electric and magnetic fields.

Copyright 2015, 2016, 2017 Ilja Honkonen
Copyright 2019, 2024, 2025 Finnish Meteorological Institute
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

* Neither the names of the copyright holders nor the names of their contributors
  may be used to endorse or promote products derived from this software
  without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


Author(s): Ilja Honkonen

Copy boundaries have no effect in this program and shouldn't
be used in config file.
*/

#include "array"
#include "cmath"
#include "cstdlib"
#include "fstream"
#include "functional"
#include "iostream"
#include "random"
#include "string"

#include "boost/filesystem.hpp"
#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#include "Eigen/Core" // must be included before gensimcell
#include "Eigen/Geometry"
#include "mpi.h" // must be included before gensimcell
#include "gensimcell.hpp"
#include "rapidjson/document.h"
#include "rapidjson/error/en.h"

#include "background_magnetic_field.hpp"
#include "boundaries/geometries.hpp"
#include "boundaries/multivariable_boundaries.hpp"
#include "boundaries/multivariable_initial_conditions.hpp"
#include "grid/options.hpp"
#include "mhd/initialize.hpp"
#include "mhd/options.hpp"
#include "mhd/variables.hpp"
#include "particle/boundaries.hpp"
#include "particle/initialize.hpp"
#include "particle/options.hpp"
#include "particle/save.hpp"
#include "particle/solve_dccrg.hpp"
#include "particle/variables.hpp"
#include "simulation_options.hpp"
#include "variable_getter.hpp"


using namespace std;

// counter for assigning unique id to particles
unsigned long long int next_particle_id;

// data stored in every cell of simulation grid
using Cell = pamhd::particle::Cell_test_particle;
// simulation data, see doi:10.1016/j.cpc.2012.12.017 or arxiv.org/abs/1212.3496
using Grid = dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>;

// background magnetic field not stored on cell faces in this program
pamhd::Bg_Magnetic_Field::data_type zero_bg_b;
struct _Bg_B {
	pamhd::Bg_Magnetic_Field::data_type& data(auto&) const {
		return zero_bg_b;
	}
	_Bg_B type() const {
		return _Bg_B();
	}
};
const auto Bg_B = _Bg_B();

const auto Vol_B = pamhd::Variable_Getter<pamhd::Magnetic_Field>();
bool pamhd::Magnetic_Field::is_stale = true;

// electric field for propagating particles
const auto Ele = pamhd::Variable_Getter<pamhd::particle::Electric_Field>();
bool pamhd::particle::Electric_Field::is_stale = true;

const auto Part_Int = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::particle::Particles_Internal()];
};
// particles moving to another cell
const auto Part_Ext = pamhd::Variable_Getter<pamhd::particle::Particles_External>();
bool pamhd::particle::Particles_External::is_stale = true;

// number of particles in above list, for allocating memory for arriving particles
const auto Nr_Ext = pamhd::Variable_Getter<pamhd::particle::Nr_Particles_External>();
bool pamhd::particle::Nr_Particles_External::is_stale = true;

const auto SInfo = pamhd::Variable_Getter<pamhd::Solver_Info>();
bool pamhd::Solver_Info::is_stale = true;

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
const auto Part_Mas = [](
	pamhd::particle::Particle_Internal& particle
)->auto& {
	return particle[pamhd::particle::Mass()];
};
const auto Part_Des = [](
	pamhd::particle::Particle_External& particle
)->auto& {
	return particle[pamhd::particle::Destination_Cell()];
};
// unused
pamhd::Magnetic_Field::data_type zero_mag_flux = {0, 0, 0};
const auto Mag_f = [](Cell& cell_data)->auto& {
	return zero_mag_flux;
};


int main(int argc, char* argv[])
{
	using std::min;

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


	// read and parse json data from configuration file
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
	pamhd::particle::Options options_particle{document};

	if (rank == 0 and options_sim.output_directory != "") {
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
		pamhd::particle::Bdy_Number_Density,
		pamhd::particle::Bdy_Temperature,
		pamhd::particle::Bdy_Velocity,
		pamhd::particle::Bdy_Nr_Particles_In_Cell,
		pamhd::particle::Bdy_Charge_Mass_Ratio,
		pamhd::particle::Bdy_Species_Mass,
		pamhd::particle::Electric_Field,
		pamhd::Magnetic_Field
	> initial_conditions;
	initial_conditions.set(document);

	std::vector<
		pamhd::boundaries::Multivariable_Boundaries<
			uint64_t,
			geometry_id_t,
			pamhd::particle::Bdy_Number_Density,
			pamhd::particle::Bdy_Temperature,
			pamhd::particle::Bdy_Velocity,
			pamhd::particle::Bdy_Nr_Particles_In_Cell,
			pamhd::particle::Bdy_Charge_Mass_Ratio,
			pamhd::particle::Bdy_Species_Mass,
			pamhd::particle::Electric_Field,
			pamhd::Magnetic_Field
		>
	> boundaries(1);
	boundaries[0].set(document);

	pamhd::Background_Magnetic_Field<
		double,
		pamhd::Magnetic_Field::data_type
	> background_B;
	background_B.set(document);


	/*
	Initialize simulation grid
	*/
	Grid grid;

	pamhd::grid::Options grid_options;
	grid_options.set(document);

	const unsigned int neighborhood_size = 1;
	const auto& number_of_cells = grid_options.get_number_of_cells();
	const auto& periodic = grid_options.get_periodic();
	grid
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
			= grid_options.get_volume(),
		cell_volume{
			simulation_volume[0] / number_of_cells[0],
			simulation_volume[1] / number_of_cells[1],
			simulation_volume[2] / number_of_cells[2]
		};

	dccrg::Cartesian_Geometry::Parameters geom_params;
	geom_params.start = grid_options.get_start();
	geom_params.level_0_cell_length = cell_volume;

	grid.set_geometry(geom_params);

	// update owner process of cells for saving into file
	for (const auto& cell: grid.local_cells()) {
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
		max_dt_particle_gyro = 0,
		max_dt_particle_flight = 0,
		simulation_time = options_sim.time_start,
		next_particle_save = options_particle.save_n;

	// set initial condition
	std::mt19937_64 random_source;

	for (auto dir: {-3,-2,-1,+1,+2,+3}) {
		zero_bg_b(dir) = {0, 0, 0};
	}
	pamhd::mhd::initialize_magnetic_field<pamhd::Magnetic_Field>(
		geometries,
		initial_conditions,
		background_B,
		grid,
		simulation_time,
		options_sim.vacuum_permeability,
		Vol_B, Mag_f, Bg_B
	);

	pamhd::particle::initialize_electric_field<pamhd::particle::Electric_Field>(
		geometries,
		initial_conditions,
		simulation_time,
		grid,
		Ele
	);

	auto nr_particles_created
		= pamhd::particle::initialize_particles<
			pamhd::particle::Particle_Internal,
			pamhd::particle::Mass,
			pamhd::particle::Charge_Mass_Ratio,
			pamhd::particle::Position,
			pamhd::particle::Velocity,
			pamhd::particle::Particle_ID,
			pamhd::particle::Species_Mass
		>(
			geometries,
			initial_conditions,
			simulation_time,
			grid,
			random_source,
			options_sim.temp2nrj,
			next_particle_id,
			grid.get_comm_size(),
			true,
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

	nr_particles_created
		+= pamhd::particle::apply_massless_boundaries<
			pamhd::particle::Particle_Internal,
			pamhd::particle::Mass,
			pamhd::particle::Charge_Mass_Ratio,
			pamhd::particle::Position,
			pamhd::particle::Velocity,
			pamhd::particle::Particle_ID,
			pamhd::particle::Species_Mass
		>(
			geometries,
			boundaries,
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
			Ele,
			Vol_B,
			Part_Int,
			Bdy_N,
			Bdy_V,
			Bdy_T,
			Bdy_Nr_Par,
			Bdy_SpM,
			Bdy_C2M
		);
	next_particle_id += nr_particles_created * grid.get_comm_size();

	if (rank == 0) {
		cout << "Done initializing particles" << endl;
	}

	/*
	Classify cells into normal, boundary and dont_solve
	*/

	Cell::set_transfer_all(true, pamhd::Solver_Info());
	pamhd::particle::set_solver_info<pamhd::Solver_Info>(
		grid, boundaries, geometries, SInfo
	);
	Cell::set_transfer_all(false, pamhd::Solver_Info());


	size_t simulated_steps = 0;
	while (simulation_time < time_end) {
		simulated_steps++;

		double
			// don't step over the final simulation time
			until_end = time_end - simulation_time,
			local_time_step = min(min(
				options_particle.gyroperiod_time_step_factor * max_dt_particle_gyro,
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

		/*
		Solve
		*/

		if (rank == 0) {
			cout << "Solving particles at time " << simulation_time
				<< " s with time step " << time_step << " s" << endl;
		}

		max_dt_particle_gyro   =
		max_dt_particle_flight = std::numeric_limits<double>::max();

		Cell::set_transfer_all(
			true,
			pamhd::particle::Electric_Field(),
			pamhd::Magnetic_Field()
		);
		grid.update_copies_of_remote_neighbors();
		Cell::set_transfer_all(
			false,
			pamhd::particle::Electric_Field(),
			pamhd::Magnetic_Field()
		);

		// E is given directly to particle propagator
		// TODO: don't use preprocessor
		std::pair<double, double> particle_max_dt{0, 0};
		particle_max_dt = pamhd::particle::solve(
			time_step, grid.outer_cells(), grid,
			background_B, options_sim.vacuum_permeability,
			false, Ele, Vol_B, Nr_Ext, Part_Int, Part_Ext,
			Part_Pos, Part_Vel, Part_C2M,
			Part_Mas, Part_Des, SInfo
		);

		max_dt_particle_flight = min(particle_max_dt.first, max_dt_particle_flight);
		max_dt_particle_gyro = min(particle_max_dt.second, max_dt_particle_gyro);

		Cell::set_transfer_all(true, pamhd::particle::Nr_Particles_External());
		grid.start_remote_neighbor_copy_updates();

		particle_max_dt = pamhd::particle::solve(
			time_step, grid.inner_cells(), grid,
			background_B, options_sim.vacuum_permeability,
			false, Ele, Vol_B, Nr_Ext, Part_Int, Part_Ext,
			Part_Pos, Part_Vel, Part_C2M,
			Part_Mas, Part_Des, SInfo
		);
		max_dt_particle_flight = min(particle_max_dt.first, max_dt_particle_flight);
		max_dt_particle_gyro = min(particle_max_dt.second, max_dt_particle_gyro);

		simulation_time += time_step;

		grid.wait_remote_neighbor_copy_update_receives();
		pamhd::particle::resize_receiving_containers<
			pamhd::particle::Nr_Particles_External,
			pamhd::particle::Particles_External
		>(grid.remote_cells(), grid);

		grid.wait_remote_neighbor_copy_update_sends();

		Cell::set_transfer_all(false, pamhd::particle::Nr_Particles_External());
		Cell::set_transfer_all(true, pamhd::particle::Particles_External());

		grid.start_remote_neighbor_copy_updates();

		pamhd::particle::incorporate_external_particles<
			pamhd::particle::Nr_Particles_External,
			pamhd::particle::Particles_Internal,
			pamhd::particle::Particles_External,
			pamhd::particle::Destination_Cell
		>(grid.inner_cells(), grid);

		grid.wait_remote_neighbor_copy_update_receives();

		pamhd::particle::incorporate_external_particles<
			pamhd::particle::Nr_Particles_External,
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


		nr_particles_created
			+= pamhd::particle::apply_massless_boundaries<
				pamhd::particle::Particle_Internal,
				pamhd::particle::Mass,
				pamhd::particle::Charge_Mass_Ratio,
				pamhd::particle::Position,
				pamhd::particle::Velocity,
				pamhd::particle::Particle_ID,
				pamhd::particle::Species_Mass
			>(
				geometries,
				boundaries,
				simulation_time,
				simulated_steps,
				grid,
				random_source,
				options_sim.temp2nrj,
				options_sim.vacuum_permeability,
				next_particle_id,
				grid.get_comm_size(),
				false,
				SInfo,
				Ele,
				Vol_B,
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

		if (
			(options_particle.save_n >= 0 and (simulation_time == 0 or simulation_time >= time_end))
			or (options_particle.save_n > 0 and simulation_time >= next_particle_save)
		) {
			if (next_particle_save <= simulation_time) {
				next_particle_save += options_particle.save_n;
			}

			if (rank == 0) {
				cout << "Saving particles at time " << simulation_time << endl;
			}

			constexpr uint64_t file_version = 4;
			if (
				not pamhd::particle::save(
					boost::filesystem::canonical(
						boost::filesystem::path(options_sim.output_directory)
					).append("particle_").generic_string(),
					grid,
					file_version,
					simulated_steps,
					simulation_time,
					0,
					0,
					options_sim.temp2nrj
				)
			) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << "): Couldn't save particle result." << std::endl;
				MPI_Finalize();
				return EXIT_FAILURE;
			}
		}
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
