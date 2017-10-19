/*
Tests parallel particle solver of PAMHD in 3 dimensions with periodic grid.

Copyright 2015, 2016, 2017 Ilja Honkonen
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
*/

#include "array"
#include "cmath"
#include "cstdlib"
#include "iostream"

#include "boost/numeric/odeint.hpp"
#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#include "Eigen/Core" // must be included before gensimcell
#include "Eigen/Geometry"
#include "mpi.h" // must be included before gensimcell
#include "gensimcell.hpp"

#include "background_magnetic_field.hpp"
#include "particle/solve_dccrg.hpp"
#include "particle/variables.hpp"

using namespace std;

using Cell = pamhd::particle::Cell_test_particle;
using Grid = dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>;


// returns reference to magnetic field for propagating particles
const auto Mag
	= [](Cell& cell_data)->typename pamhd::Magnetic_Field::data_type&{
		return cell_data[pamhd::Magnetic_Field()];
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
	const std::array<uint64_t, 3> number_of_cells{{ 10, 10, 10}};
	if (
		not grid.initialize(
			number_of_cells,
			comm,
			"RANDOM",
			neighborhood_size,
			0,
			true, true, true
		)
	) {
		std::cerr << __FILE__ << "(" << __LINE__ << "): "
			<< "Couldn't initialize one or more grids."
			<< std::endl;
		abort();
	}

	// set grid geometry
	dccrg::Cartesian_Geometry::Parameters geom_params;
	geom_params.start = {{0.0, 0.0, 0.0}};
	geom_params.level_0_cell_length = {{1.0, 1.0, 1.0}};

	if (not grid.set_geometry(geom_params)) {
		std::cerr << __FILE__ << "(" << __LINE__ << "): "
			<< "Couldn't set geometry of one or more grids."
			<< std::endl;
		abort();
	}

	grid.balance_load();

	const auto
		inner_cell_ids = grid.get_local_cells_not_on_process_boundary(),
		outer_cell_ids = grid.get_local_cells_on_process_boundary(),
		remote_cell_ids = grid.get_remote_cells_on_process_boundary(),
		cell_ids = grid.get_cells();

	// initial condition
	for (const auto& cell_id: cell_ids) {
		auto* const cell_ptr = grid[cell_id];
		if (cell_ptr == nullptr) {
			std::cerr << __FILE__ << "(" << __LINE__ << ")" << std::endl;
			abort();
		}
		auto& cell_data = *cell_ptr;

		const auto cell_center = grid.geometry.get_center(cell_id);

		pamhd::particle::Particle_Internal particle;
		Part_Pos(particle) = {
			cell_center[0],
			cell_center[1],
			cell_center[2]
		};
		Part_Vel(particle) = {-1.0, 1.0, 1.0};

		Part_Mas(particle) =
		Part_C2M(particle) = 0;

		Part_Int(cell_data).push_back(particle);

		Ele(cell_data) =
		Mag(cell_data) = {0, 0, 0};

		Nr_Ext(cell_data) = Part_Ext(cell_data).size();
	}
	// allocate copies of remote neighbor cells
	grid.update_copies_of_remote_neighbors();

	pamhd::Background_Magnetic_Field<double, Eigen::Vector3d> bg_B;

	// short hand notation for calling solvers
	auto solve = [&bg_B](
		const std::vector<uint64_t>& cell_ids,
		Grid& grid
	) {
		pamhd::particle::solve<
			boost::numeric::odeint::runge_kutta_fehlberg78<pamhd::particle::state_t>
		>(
			1.0,
			cell_ids,
			grid,
			bg_B,
			1,
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
		);
	};

	auto resize_receiving = [](
		const std::vector<uint64_t>& cell_ids,
		Grid& grid
	) {
		pamhd::particle::resize_receiving_containers<
			pamhd::particle::Nr_Particles_External,
			pamhd::particle::Particles_External
		>(cell_ids, grid);
	};

	auto incorporate_external = [](
		const std::vector<uint64_t>& cell_ids,
		Grid& grid
	) {
		pamhd::particle::incorporate_external_particles<
			pamhd::particle::Nr_Particles_External,
			pamhd::particle::Particles_Internal,
			pamhd::particle::Particles_External,
			pamhd::particle::Destination_Cell
		>(cell_ids, grid);
	};

	auto remove_external = [](
		const std::vector<uint64_t>& cell_ids,
		Grid& grid
	) {
		pamhd::particle::remove_external_particles<
			pamhd::particle::Nr_Particles_External,
			pamhd::particle::Particles_External
		>(cell_ids, grid);
	};


	for (size_t step = 0; step < 10; step++) {
		solve(outer_cell_ids, grid);

		Cell::set_transfer_all(
			true,
			pamhd::particle::Electric_Field(),
			pamhd::Magnetic_Field(),
			pamhd::particle::Nr_Particles_External()
		);
		grid.start_remote_neighbor_copy_updates();

		solve(inner_cell_ids, grid);

		grid.wait_remote_neighbor_copy_update_receives();
		resize_receiving(remote_cell_ids, grid);

		grid.wait_remote_neighbor_copy_update_sends();

		Cell::set_transfer_all(
			false,
			pamhd::particle::Electric_Field(),
			pamhd::Magnetic_Field(),
			pamhd::particle::Nr_Particles_External()
		);
		Cell::set_transfer_all(
			true,
			pamhd::particle::Particles_External()
		);

		grid.start_remote_neighbor_copy_updates();

		incorporate_external(inner_cell_ids, grid);

		grid.wait_remote_neighbor_copy_update_receives();

		incorporate_external(outer_cell_ids, grid);

		remove_external(inner_cell_ids, grid);

		grid.wait_remote_neighbor_copy_update_sends();
		Cell::set_transfer_all(
			false,
			pamhd::particle::Particles_External()
		);

		remove_external(outer_cell_ids, grid);


		// check that solution is correct
		int total_particles_local = 0, total_particles = 0;

		for (const auto& cell_id: cell_ids) {
			auto* const cell_ptr = grid[cell_id];
			if (cell_ptr == nullptr) {
				std::cerr << __FILE__ << "(" << __LINE__ << ")" << std::endl;
				abort();
			}
			auto& cell_data = *cell_ptr;

			total_particles_local
				+= cell_data[pamhd::particle::Particles_Internal()].size()
				+ cell_data[pamhd::particle::Particles_External()].size();

			if (cell_data[pamhd::particle::Particles_Internal()].size() > 1) {
				std::cerr << __FILE__ << "(" << __LINE__ << ") "
					<< "Incorrect number of internal particles in cell " << cell_id
					<< ": " << cell_data[pamhd::particle::Particles_Internal()].size()
					<< std::endl;
				abort();
			}
		}

		MPI_Allreduce(
			&total_particles_local,
			&total_particles,
			1,
			MPI_INT,
			MPI_SUM,
			comm
		);
		if (total_particles != 1000) {
			if (grid.get_rank() == 0) {
				std::cerr << __FILE__ << "(" << __LINE__ << ") "
					<< "Incorrect total number of particles at end of step "
					<< step << ": " << total_particles
					<< std::endl;
			}
			abort();
		}
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
