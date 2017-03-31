/*
Tests parallel particle solver of PAMHD in dipole field.

Copyright 2015, 2016, 2017 Ilja Honkonen
All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "array"
#include "cmath"
#include "cstdlib"
#include "iostream"
#include "utility"
#include "vector"

#include "mhd/background_magnetic_field.hpp"
#include "boost/numeric/odeint.hpp"
#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#include "Eigen/Core" // must be included before gensimcell
#include "Eigen/Geometry"
#include "mpi.h" // must be included before gensimcell
#include "gensimcell.hpp"
#include "rapidjson/document.h"
#include "rapidjson/error/en.h"

#include "mhd/variables.hpp"
#include "particle/save.hpp"
#include "particle/solve_dccrg.hpp"
#include "particle/variables.hpp"

using namespace std;

using Cell = pamhd::particle::Cell_test_particle;
using Grid = dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>;


// returns reference to magnetic field for propagating particles
const auto Mag
	= [](Cell& cell_data)->typename pamhd::particle::Magnetic_Field::data_type&{
		return cell_data[pamhd::particle::Magnetic_Field()];
	};
// electric field for propagating particles
const auto Ele
	= [](Cell& cell_data)->typename pamhd::particle::Electric_Field::data_type&{
		return cell_data[pamhd::particle::Electric_Field()];
	};
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

// given a particle these return references to particle's parameters
const auto Part_Pos
	= [](pamhd::particle::Particle_Internal& particle)->typename pamhd::particle::Position::data_type&{
		return particle[pamhd::particle::Position()];
	};
const auto Part_Vel
	= [](pamhd::particle::Particle_Internal& particle)->typename pamhd::particle::Velocity::data_type&{
		return particle[pamhd::particle::Velocity()];
	};
const auto Part_C2M
	= [](pamhd::particle::Particle_Internal& particle)->typename pamhd::particle::Charge_Mass_Ratio::data_type&{
		return particle[pamhd::particle::Charge_Mass_Ratio()];
	};
const auto Part_Mas
	= [](pamhd::particle::Particle_Internal& particle)->typename pamhd::particle::Mass::data_type&{
		return particle[pamhd::particle::Mass()];
	};
const auto Part_Des
	= [](pamhd::particle::Particle_External& particle)->typename pamhd::particle::Destination_Cell::data_type&{
		return particle[pamhd::particle::Destination_Cell()];
	};

int main(int argc, char* argv[])
{
	constexpr double Re = 6.371e6; // radius of earth

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

	// intialize Zoltan
	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
		std::cerr << "Zoltan_Initialize failed." << std::endl;
		abort();
	}


	Grid grid;

	const unsigned int neighborhood_size = 1;
	const std::array<uint64_t, 3> number_of_cells{{25, 25, 6}};
	if (
		not grid.initialize(
			number_of_cells,
			comm,
			"RANDOM",
			neighborhood_size,
			0,
			false, false, false
		)
	) {
		std::cerr << __FILE__ << "(" << __LINE__ << "): "
			<< "Couldn't initialize one or more grids."
			<< std::endl;
		abort();
	}

	// set grid geometry
	dccrg::Cartesian_Geometry::Parameters geom_params;
	geom_params.start = {{-6.0 * Re, -6.0 * Re, -0.6 * Re}};
	geom_params.level_0_cell_length = {{Re / 5, Re / 5, Re / 5}};

	if (not grid.set_geometry(geom_params)) {
		std::cerr << __FILE__ << "(" << __LINE__ << "): "
			<< "Couldn't set geometry of one or more grids."
			<< std::endl;
		abort();
	}

	grid.balance_load();

	const auto
		inner_cells = grid.get_local_cells_not_on_process_boundary(),
		outer_cells = grid.get_local_cells_on_process_boundary(),
		remote_cells = grid.get_remote_cells_on_process_boundary(),
		cells = grid.get_cells();

	// background magnetic field
	pamhd::mhd::Background_Magnetic_Field<Eigen::Vector3d> bg_B;

	rapidjson::Document document;
	document.Parse(
		"{\"background-magnetic-field\": {"
		"	\"dipoles\": ["
		"		{\"moment\": [0, 0, -7.94e22], \"position\": [0, 0, 0]}"
		"	],"
		"	\"minimum-distance\": 1e5"
		"}}"
	);
	if (document.HasParseError()) {
		std::cerr << "Couldn't parse json data in file " << argv[1]
			<< " at character position " << document.GetErrorOffset()
			<< ": " << rapidjson::GetParseError_En(document.GetParseError())
			<< std::endl;
		MPI_Finalize();
		return EXIT_FAILURE;
	}

	bg_B.set(document);

	// create particles at equator 5 Re from origin
	int initial_particles_local = 0, initial_particles = 0;
	for (const auto& cell_id: cells) {
		auto* const cell_ptr = grid[cell_id];
		if (cell_ptr == nullptr) {
			std::cerr << __FILE__ << "(" << __LINE__ << ")" << std::endl;
			abort();
		}
		auto& cell_data = *cell_ptr;

		const auto cell_center = grid.geometry.get_center(cell_id);
		const Eigen::Vector3d r{
			cell_center[0],
			cell_center[1],
			cell_center[2]
		};
		const auto distance = r.norm();

		Ele(cell_data) = {0, 0, 0};
		Mag(cell_data)
			= bg_B.get_background_field(
				Eigen::Vector3d{cell_center[0], cell_center[1], cell_center[2]},
				1.257e-06
			);

		if (distance > (5 + 1.0/6.0) * Re or distance < (5 - 1.0/6.0) * Re) {
			continue;
		}

		if (std::fabs(cell_center[2]) > Re / 5) {
			continue;
		}

		pamhd::particle::Particle_Internal particle;
		Part_Pos(particle) = r;
		Part_Vel(particle) = 2e6 * r / r.norm(); // m / s
		Part_Mas(particle) = 0;
		Part_C2M(particle) = 95788335.8;
		Part_Int(cell_data).push_back(particle);

		Part_Pos(particle) += Eigen::Vector3d{0.05 * Re, 0.05 * Re, 0.05 * Re};
		Part_C2M(particle) += 1e6;
		Part_Int(cell_data).push_back(particle);

		Nr_Ext(cell_data) = Part_Ext(cell_data).size();

		initial_particles_local += Part_Int(cell_data).size();
	}
	// allocate copies of remote neighbor cells
	grid.update_copies_of_remote_neighbors();

	MPI_Allreduce(
		&initial_particles_local,
		&initial_particles,
		1,
		MPI_INT,
		MPI_SUM,
		comm
	);


	double
		max_dt = 0,
		start_time = 0,
		end_time = 12,
		save_particle_n = 2,
		next_particle_save = save_particle_n,
		simulation_time = start_time;

	size_t simulated_steps = 0;
	while (simulation_time < end_time) {
		simulated_steps++;

		double
			// don't step over the final simulation time
			until_end = end_time - simulation_time,
			local_time_step = std::min(0.5 * max_dt, until_end),
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
				<< ": Couldn't set reduce time step."
				<< std::endl;
			abort();
		}

		if (grid.get_rank() == 0) {
			/*cout << "Solving particles at time " << simulation_time
				<< " s with time step " << time_step << " s" << endl;*/
		}

		max_dt = std::numeric_limits<double>::max();

		max_dt = std::min(
			max_dt,
			pamhd::particle::solve<
				boost::numeric::odeint::runge_kutta_fehlberg78<pamhd::particle::state_t>
			>(
				time_step,
				outer_cells,
				grid,
				false,
				Ele,
				Mag,
				Nr_Ext,
				Part_Int,
				Part_Ext,
				Part_Pos,
				Part_Vel,
				Part_C2M,
				Part_Mas,
				Part_Des
			)
		);

		Cell::set_transfer_all(
			true,
			pamhd::particle::Electric_Field(),
			pamhd::particle::Magnetic_Field(),
			pamhd::particle::Nr_Particles_External()
		);
		grid.start_remote_neighbor_copy_updates();

		max_dt = std::min(
			max_dt,
			pamhd::particle::solve<
				boost::numeric::odeint::runge_kutta_fehlberg78<pamhd::particle::state_t>
			>(
				time_step,
				inner_cells,
				grid,
				false,
				Ele,
				Mag,
				Nr_Ext,
				Part_Int,
				Part_Ext,
				Part_Pos,
				Part_Vel,
				Part_C2M,
				Part_Mas,
				Part_Des
			)
		);

		simulation_time += time_step;

		grid.wait_remote_neighbor_copy_update_receives();
		pamhd::particle::resize_receiving_containers<
			pamhd::particle::Nr_Particles_External,
			pamhd::particle::Particles_External
		>(remote_cells, grid);

		grid.wait_remote_neighbor_copy_update_sends();

		Cell::set_transfer_all(
			false,
			pamhd::particle::Electric_Field(),
			pamhd::particle::Magnetic_Field(),
			pamhd::particle::Nr_Particles_External()
		);
		Cell::set_transfer_all(
			true,
			pamhd::particle::Particles_External()
		);

		grid.start_remote_neighbor_copy_updates();

		pamhd::particle::incorporate_external_particles<
			pamhd::particle::Nr_Particles_External,
			pamhd::particle::Particles_Internal,
			pamhd::particle::Particles_External,
			pamhd::particle::Destination_Cell
		>(inner_cells, grid);

		grid.wait_remote_neighbor_copy_update_receives();

		pamhd::particle::incorporate_external_particles<
			pamhd::particle::Nr_Particles_External,
			pamhd::particle::Particles_Internal,
			pamhd::particle::Particles_External,
			pamhd::particle::Destination_Cell
		>(outer_cells, grid);

		pamhd::particle::remove_external_particles<
			pamhd::particle::Nr_Particles_External,
			pamhd::particle::Particles_External
		>(inner_cells, grid);

		grid.wait_remote_neighbor_copy_update_sends();
		Cell::set_transfer_all(
			false,
			pamhd::particle::Particles_External()
		);

		pamhd::particle::remove_external_particles<
			pamhd::particle::Nr_Particles_External,
			pamhd::particle::Particles_External
		>(outer_cells, grid);


		if (
			(save_particle_n >= 0 and (simulation_time == 0 or simulation_time >= end_time))
			or (save_particle_n > 0 and simulation_time >= next_particle_save)
		) {
			if (next_particle_save <= simulation_time) {
				next_particle_save += save_particle_n;
			}

			if (rank == 0) {
				//cout << "Saving particles at time " << simulation_time << endl;
			}

			if (
				not pamhd::particle::save<
					pamhd::particle::Electric_Field,
					pamhd::particle::Magnetic_Field,
					pamhd::mhd::Electric_Current_Density,
					pamhd::particle::Nr_Particles_Internal,
					pamhd::particle::Particles_Internal
				>(
					"tests/particle/",
					grid,
					simulation_time,
					0,
					0,
					0
				)
			) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
					"Couldn't save particle result, did you run this "
					"from the root pamhd directory?"
					<< std::endl;
				abort();
			}
		}

		// check for errors
		int total_particles_local = 0, total_particles = 0;

		for (const auto& cell_id: cells) {
			auto* const cell_ptr = grid[cell_id];
			if (cell_ptr == nullptr) {
				std::cerr << __FILE__ << "(" << __LINE__ << ")" << std::endl;
				abort();
			}
			auto& cell_data = *cell_ptr;

			total_particles_local
				+= cell_data[pamhd::particle::Particles_Internal()].size()
				+ cell_data[pamhd::particle::Particles_External()].size();
		}

		MPI_Allreduce(
			&total_particles_local,
			&total_particles,
			1,
			MPI_INT,
			MPI_SUM,
			comm
		);
		if (total_particles != initial_particles) {
			if (grid.get_rank() == 0) {
				std::cerr << __FILE__ << "(" << __LINE__ << ") "
					<< "Incorrect total number of particles at end of step "
					<< simulated_steps << ": " << total_particles
					<< std::endl;
			}
			abort();
		}
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
