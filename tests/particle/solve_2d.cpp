/*
Tests parallel particle solver of PAMHD in 2 dimensions.

Copyright 2015, 2016, 2017 Ilja Honkonen
Copyright 2019, 2025 Finnish Meteorological Institute
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
*/

#include "array"
#include "cmath"
#include "cstdlib"
#include "iostream"

#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#include "Eigen/Core" // must be included before gensimcell
#include "Eigen/Geometry"
#include "mpi.h" // must be included before gensimcell
#include "gensimcell.hpp"

#include "background_magnetic_field.hpp"
#include "particle/solve_dccrg.hpp"
#include "particle/variables.hpp"
#include "variable_getter.hpp"

using namespace std;

using Cell = pamhd::particle::Cell_test_particle;
using Grid = dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>;


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


	Grid grid_x, grid_y, grid_z;

	const unsigned int neighborhood_size = 1;
	const std::array<uint64_t, 3>
		number_of_cells_x{{ 1, 10, 10}},
		number_of_cells_y{{10,  1, 10}},
		number_of_cells_z{{10, 10,  1}};
	grid_x
		.set_neighborhood_length(neighborhood_size)
		.set_maximum_refinement_level(0)
		.set_load_balancing_method("RCB")
		.set_initial_length(number_of_cells_x)
		.initialize(comm);
	grid_y
		.set_neighborhood_length(neighborhood_size)
		.set_maximum_refinement_level(0)
		.set_load_balancing_method("RCB")
		.set_initial_length(number_of_cells_y)
		.initialize(comm);
	grid_z
		.set_neighborhood_length(neighborhood_size)
		.set_maximum_refinement_level(0)
		.set_load_balancing_method("RCB")
		.set_initial_length(number_of_cells_z)
		.initialize(comm);

	// set grid geometry
	dccrg::Cartesian_Geometry::Parameters geom_params;
	geom_params.start = {{0.0, 0.0, 0.0}};
	geom_params.level_0_cell_length = {{1.0, 1.0, 1.0}};

	grid_x.set_geometry(geom_params);
	grid_y.set_geometry(geom_params);
	grid_z.set_geometry(geom_params);

	// use same domain decomposition in all grids
	for (size_t cell_id = 1; cell_id <= 100; cell_id++) {
		grid_x.pin(cell_id, cell_id % grid_x.get_comm_size());
		grid_y.pin(cell_id, cell_id % grid_y.get_comm_size());
		grid_z.pin(cell_id, cell_id % grid_z.get_comm_size());
	}
	grid_x.balance_load(false);
	grid_y.balance_load(false);
	grid_z.balance_load(false);

	// initial condition
	for (const auto& cell: grid_x.local_cells()) {
		auto
			*const cell_ptr_x = grid_x[cell.id],
			*const cell_ptr_y = grid_y[cell.id],
			*const cell_ptr_z = grid_z[cell.id];
		if (
			cell_ptr_x == nullptr
			or cell_ptr_y == nullptr
			or cell_ptr_z == nullptr
		) {
			std::cerr << __FILE__ << "(" << __LINE__ << ")" << std::endl;
			abort();
		}
		auto
			&cell_data_x = *cell_ptr_x,
			&cell_data_y = *cell_ptr_y,
			&cell_data_z = *cell_ptr_z;

		const auto
			cell_center_x = grid_x.geometry.get_center(cell.id),
			cell_center_y = grid_y.geometry.get_center(cell.id),
			cell_center_z = grid_z.geometry.get_center(cell.id);

		pamhd::particle::Particle_Internal particle_x, particle_y, particle_z;
		Part_Pos(particle_x) = {
			cell_center_x[0],
			cell_center_x[1],
			cell_center_x[2]
		};
		Part_Pos(particle_y) = {
			cell_center_y[0],
			cell_center_y[1],
			cell_center_y[2]
		};
		Part_Pos(particle_z) = {
			cell_center_z[0],
			cell_center_z[1],
			cell_center_z[2]
		};
		Part_Vel(particle_x) = {  0, 1.0, 1.0};
		Part_Vel(particle_y) = {1.0,   0, 1.0};
		Part_Vel(particle_z) = {1.0, 1.0,   0};

		Part_Mas(particle_x) =
		Part_Mas(particle_y) =
		Part_Mas(particle_z) =
		Part_C2M(particle_x) =
		Part_C2M(particle_y) =
		Part_C2M(particle_z) = 0;

		Part_Int(cell_data_x).push_back(particle_x);
		Part_Int(cell_data_y).push_back(particle_y);
		Part_Int(cell_data_z).push_back(particle_z);

		Ele.data(cell_data_x) =
		Vol_B.data(cell_data_x) =
		Ele.data(cell_data_y) =
		Vol_B.data(cell_data_y) =
		Ele.data(cell_data_z) =
		Vol_B.data(cell_data_z) = {0, 0, 0};

		Nr_Ext.data(cell_data_x) = Part_Ext.data(cell_data_x).size();
		Nr_Ext.data(cell_data_y) = Part_Ext.data(cell_data_y).size();
		Nr_Ext.data(cell_data_z) = Part_Ext.data(cell_data_z).size();
	}
	// allocate copies of remote neighbor cells
	grid_x.update_copies_of_remote_neighbors();
	grid_y.update_copies_of_remote_neighbors();
	grid_z.update_copies_of_remote_neighbors();

	pamhd::Background_Magnetic_Field<double, Eigen::Vector3d> bg_B;

	// short hand notation for calling solvers
	auto solve = [&bg_B](
		const auto& cells,
		Grid& grid
	) {
		pamhd::particle::solve(
			1.0, cells, grid, bg_B, 1, false,
			Ele, Vol_B, Nr_Ext, Part_Int, Part_Ext,
			Part_Pos, Part_Vel, Part_C2M,
			Part_Mas, Part_Des, SInfo
		);
	};

	using namespace pamhd::particle;
	using NPE = Nr_Particles_External;
	using PE = Particles_External;
	using PI = Particles_Internal;
	using DC = Destination_Cell;

	for (size_t step = 0; step < 10; step++) {
		solve(grid_x.outer_cells(), grid_x);
		solve(grid_y.outer_cells(), grid_y);
		solve(grid_z.outer_cells(), grid_z);

		Cell::set_transfer_all(true, Electric_Field(), pamhd::Magnetic_Field(), NPE());
		grid_x.start_remote_neighbor_copy_updates();
		grid_y.start_remote_neighbor_copy_updates();
		grid_z.start_remote_neighbor_copy_updates();

		solve(grid_x.inner_cells(), grid_x);
		solve(grid_y.inner_cells(), grid_y);
		solve(grid_z.inner_cells(), grid_z);

		grid_x.wait_remote_neighbor_copy_update_receives();
		grid_y.wait_remote_neighbor_copy_update_receives();
		grid_z.wait_remote_neighbor_copy_update_receives();
		resize_receiving_containers<NPE, PE>(grid_x.remote_cells(), grid_x);
		resize_receiving_containers<NPE, PE>(grid_y.remote_cells(), grid_y);
		resize_receiving_containers<NPE, PE>(grid_z.remote_cells(), grid_z);

		grid_x.wait_remote_neighbor_copy_update_sends();
		grid_y.wait_remote_neighbor_copy_update_sends();
		grid_z.wait_remote_neighbor_copy_update_sends();

		Cell::set_transfer_all(false, Electric_Field(), pamhd::Magnetic_Field(), NPE());
		Cell::set_transfer_all(true, PE());

		grid_x.start_remote_neighbor_copy_updates();
		grid_y.start_remote_neighbor_copy_updates();
		grid_z.start_remote_neighbor_copy_updates();

		incorporate_external_particles<NPE, PI, PE, DC>(grid_x.inner_cells(), grid_x);
		incorporate_external_particles<NPE, PI, PE, DC>(grid_y.inner_cells(), grid_y);
		incorporate_external_particles<NPE, PI, PE, DC>(grid_z.inner_cells(), grid_z);

		grid_x.wait_remote_neighbor_copy_update_receives();
		grid_y.wait_remote_neighbor_copy_update_receives();
		grid_z.wait_remote_neighbor_copy_update_receives();

		incorporate_external_particles<NPE, PI, PE, DC>(grid_x.outer_cells(), grid_x);
		incorporate_external_particles<NPE, PI, PE, DC>(grid_y.outer_cells(), grid_y);
		incorporate_external_particles<NPE, PI, PE, DC>(grid_z.outer_cells(), grid_z);

		remove_external_particles<NPE, PE>(grid_x.inner_cells(), grid_x);
		remove_external_particles<NPE, PE>(grid_y.inner_cells(), grid_y);
		remove_external_particles<NPE, PE>(grid_z.inner_cells(), grid_z);

		grid_x.wait_remote_neighbor_copy_update_sends();
		grid_y.wait_remote_neighbor_copy_update_sends();
		grid_z.wait_remote_neighbor_copy_update_sends();
		Cell::set_transfer_all(false, PE());

		remove_external_particles<NPE, PE>(grid_x.outer_cells(), grid_x);
		remove_external_particles<NPE, PE>(grid_y.outer_cells(), grid_y);
		remove_external_particles<NPE, PE>(grid_z.outer_cells(), grid_z);


		// check that solution is correct
		std::array<int, 3>
			total_particles_local{{0, 0, 0}},
			total_particles{{0, 0, 0}};

		for (const auto& cell: grid_x.local_cells()) {
			auto
				*const cell_ptr_x = grid_x[cell.id],
				*const cell_ptr_y = grid_y[cell.id],
				*const cell_ptr_z = grid_z[cell.id];
			if (
				cell_ptr_x == nullptr
				or cell_ptr_y == nullptr
				or cell_ptr_z == nullptr
			) {
				std::cerr << __FILE__ << "(" << __LINE__ << ")" << std::endl;
				abort();
			}
			auto
				&cell_data_x = *cell_ptr_x,
				&cell_data_y = *cell_ptr_y,
				&cell_data_z = *cell_ptr_z;

			total_particles_local[0]
				+= cell_data_x[pamhd::particle::Particles_Internal()].size()
				+ cell_data_x[pamhd::particle::Particles_External()].size();
			total_particles_local[1]
				+= cell_data_y[pamhd::particle::Particles_Internal()].size()
				+ cell_data_y[pamhd::particle::Particles_External()].size();
			total_particles_local[2]
				+= cell_data_z[pamhd::particle::Particles_Internal()].size()
				+ cell_data_z[pamhd::particle::Particles_External()].size();

			if (cell_data_x[pamhd::particle::Particles_Internal()].size() > 1) {
				std::cerr << __FILE__ << "(" << __LINE__ << ") "
					<< "Incorrect number of internal particles in cell " << cell.id
					<< ": " << cell_data_x[pamhd::particle::Particles_Internal()].size()
					<< std::endl;
				abort();
			}
			if (cell_data_y[pamhd::particle::Particles_Internal()].size() > 1) {
				std::cerr << __FILE__ << "(" << __LINE__ << ") "
					<< "Incorrect number of internal particles in cell " << cell.id
					<< ": " << cell_data_y[pamhd::particle::Particles_Internal()].size()
					<< std::endl;
				abort();
			}
			if (cell_data_z[pamhd::particle::Particles_Internal()].size() > 1) {
				std::cerr << __FILE__ << "(" << __LINE__ << ") "
					<< "Incorrect number of internal particles in cell " << cell.id
					<< ": " << cell_data_z[pamhd::particle::Particles_Internal()].size()
					<< std::endl;
				abort();
			}
		}

		MPI_Allreduce(
			total_particles_local.data(),
			total_particles.data(),
			3,
			MPI_INT,
			MPI_SUM,
			comm
		);
		if (total_particles[0] != int(std::pow(9 - step, 2))) {
			if (grid_x.get_rank() == 0) {
				std::cerr << __FILE__ << "(" << __LINE__ << ") "
					<< "Incorrect total number of particles at end of step "
					<< step << ": " << total_particles[0]
					<< std::endl;
			}
			abort();
		}
		if (total_particles[1] != int(std::pow(9 - step, 2))) {
			if (grid_y.get_rank() == 0) {
				std::cerr << __FILE__ << "(" << __LINE__ << ") "
					<< "Incorrect total number of particles at end of step "
					<< step << ": " << total_particles[1]
					<< std::endl;
			}
			abort();
		}
		if (total_particles[2] != int(std::pow(9 - step, 2))) {
			if (grid_z.get_rank() == 0) {
				std::cerr << __FILE__ << "(" << __LINE__ << ") "
					<< "Incorrect total number of particles at end of step "
					<< step << ": " << total_particles[2]
					<< std::endl;
			}
			abort();
		}
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
