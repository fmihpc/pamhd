/*
Tests parallel particle solver of PAMHD in 3 dimensions with periodic grid.

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

const auto Max_v_part = pamhd::Variable_Getter<pamhd::particle::Max_Spatial_Velocity>();
bool pamhd::particle::Max_Spatial_Velocity::is_stale = true;

const auto Max_ω_part = pamhd::Variable_Getter<pamhd::particle::Max_Angular_Velocity>();
bool pamhd::particle::Max_Angular_Velocity::is_stale = true;

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


	Grid grid;

	const unsigned int neighborhood_size = 1;
	const std::array<uint64_t, 3> number_of_cells{{ 10, 10, 10}};
	grid
		.set_neighborhood_length(neighborhood_size)
		.set_maximum_refinement_level(0)
		.set_load_balancing_method("RANDOM")
		.set_periodic(true, true, true)
		.set_initial_length(number_of_cells)
		.initialize(comm)
		.balance_load();

	// set grid geometry
	dccrg::Cartesian_Geometry::Parameters geom_params;
	geom_params.start = {{0.0, 0.0, 0.0}};
	geom_params.level_0_cell_length = {{1.0, 1.0, 1.0}};

	grid.set_geometry(geom_params);

	// initial condition
	for (const auto& cell: grid.local_cells()) {
		const auto cell_center = grid.geometry.get_center(cell.id);

		pamhd::particle::Particle_Internal particle;
		Part_Pos(particle) = {
			cell_center[0],
			cell_center[1],
			cell_center[2]
		};
		Part_Vel(particle) = {-1.0, 1.0, 1.0};

		Part_Mas(particle) =
		Part_C2M(particle) = 0;

		Part_Int(*cell.data).push_back(particle);

		Ele.data(*cell.data) =
		Vol_B.data(*cell.data) = {0, 0, 0};

		Nr_Ext.data(*cell.data) = Part_Ext.data(*cell.data).size();
	}
	// allocate copies of remote neighbor cells
	grid.update_copies_of_remote_neighbors();

	pamhd::Background_Magnetic_Field<double, Eigen::Vector3d> bg_B;

	// short hand notation for calling solvers
	auto solve = [&bg_B](
		const auto& cells,
		Grid& grid
	) {
		pamhd::particle::solve(
			1.0, cells, grid, bg_B, 1, false, Ele, Vol_B,
			Nr_Ext, Part_Int, Part_Ext, Max_v_part,
			Max_ω_part, Part_Pos, Part_Vel, Part_C2M,
			Part_Mas, Part_Des, SInfo
		);
	};

	using namespace pamhd::particle;
	using NPE = Nr_Particles_External;
	using PE = Particles_External;
	using PI = Particles_Internal;
	using DC = Destination_Cell;

	for (size_t step = 0; step < 10; step++) {
		solve(grid.outer_cells(), grid);

		Cell::set_transfer_all(true, Electric_Field(), pamhd::Magnetic_Field(), NPE());
		grid.start_remote_neighbor_copy_updates();

		solve(grid.inner_cells(), grid);

		grid.wait_remote_neighbor_copy_update_receives();
		resize_receiving_containers<NPE, PE>(grid.remote_cells(), grid);

		grid.wait_remote_neighbor_copy_update_sends();

		Cell::set_transfer_all(false, Electric_Field(), pamhd::Magnetic_Field(), NPE());
		Cell::set_transfer_all(true, PE());

		grid.start_remote_neighbor_copy_updates();

		incorporate_external_particles<NPE, PI, PE, DC>(grid.inner_cells(), grid);

		grid.wait_remote_neighbor_copy_update_receives();

		incorporate_external_particles<NPE, PI, PE, DC>(grid.outer_cells(), grid);

		remove_external_particles<NPE, PE>(grid.inner_cells(), grid);

		grid.wait_remote_neighbor_copy_update_sends();
		Cell::set_transfer_all(false, PE());

		remove_external_particles<NPE, PE>(grid.outer_cells(), grid);


		// check that solution is correct
		int total_particles_local = 0, total_particles = 0;

		for (const auto& cell: grid.local_cells()) {
			total_particles_local
				+= (*cell.data)[pamhd::particle::Particles_Internal()].size()
				+ (*cell.data)[pamhd::particle::Particles_External()].size();

			if ((*cell.data)[pamhd::particle::Particles_Internal()].size() > 1) {
				std::cerr << __FILE__ << "(" << __LINE__ << ") "
					<< "Incorrect number of internal particles in cell " << cell.id
					<< ": " << (*cell.data)[pamhd::particle::Particles_Internal()].size()
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
