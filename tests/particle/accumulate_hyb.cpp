/*
Test for plain hybrid PIC accumulator.

Copyright 2025 Finnish Meteorological Institute
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
#include "algorithm"
#include "cmath"
#include "cstdlib"
#include "exception"
#include "iostream"
#include "utility"
#include "vector"

#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#include "mpi.h" // must be included before gensimcell.hpp
#include "gensimcell.hpp"
#include "prettyprint.hpp"

#include "common_functions.hpp"
#include "particle/accumulate_dccrg.hpp"
#include "particle/variables.hpp"
#include "variable_getter.hpp"


double charge_density(const std::array<double, 3>& r) {
	return std::sin(2*M_PI*r[0]/7)*std::cos(M_PI*r[1]/5)+r[2]*r[2];
}

std::array<double, 3> current_density(const std::array<double, 3>& r) {
	return {
		(charge_density(r) - r[2]*r[2]) / std::cos(M_PI*r[1]/5),
		(charge_density(r) - r[2]*r[2]) / std::sin(2*M_PI*r[0]/7),
		charge_density(r) - std::sin(2*M_PI*r[0]/7)*std::cos(M_PI*r[1]/5)
	};
}

bool pamhd::particle::Volume_Ion_Charge_Density::is_stale = true;
const auto Vol_Qi = pamhd::Variable_Getter<
	pamhd::particle::Volume_Ion_Charge_Density>();

bool pamhd::particle::Volume_Ion_Current_Density::is_stale = true;
const auto Vol_Ji = pamhd::Variable_Getter<
	pamhd::particle::Volume_Ion_Current_Density>();

bool pamhd::Cell_Type::is_stale = true;
const auto CType = pamhd::Variable_Getter<pamhd::Cell_Type>();

const auto Nr_Int = [](auto& cell_data)->auto& {
	return cell_data[pamhd::particle::Nr_Particles_Internal()];
};

const auto Part_Int = [](auto& cell_data)->auto& {
	return cell_data[pamhd::particle::Particles_Internal()];
};

const auto PPos = [](
	pamhd::particle::Particle_Internal& particle
)->auto& {
	return particle[pamhd::particle::Position()];
};

const auto PVel = [](
	pamhd::particle::Particle_Internal& particle
)->auto& {
	return particle[pamhd::particle::Velocity()];
};

const auto PC2M = [](
	pamhd::particle::Particle_Internal& particle
)->auto& {
	return particle[pamhd::particle::Charge_Mass_Ratio()];
};

const auto PMas = [](
	pamhd::particle::Particle_Internal& particle
)->auto& {
	return particle[pamhd::particle::Mass()];
};

const auto PSMas = [](
	pamhd::particle::Particle_Internal& particle
)->auto& {
	return particle[pamhd::particle::Species_Mass()];
};

using Cell = gensimcell::Cell<
	gensimcell::Optional_Transfer,
	pamhd::particle::Particles_Internal,
	pamhd::particle::Nr_Particles_Internal,
	pamhd::particle::Volume_Ion_Charge_Density,
	pamhd::particle::Volume_Ion_Current_Density,
	pamhd::Cell_Type
>;

using Grid = dccrg::Dccrg<
	Cell,
	dccrg::Cartesian_Geometry,
	std::tuple<pamhd::grid::Cell_Is_Local>,
	std::tuple<
		pamhd::grid::Face_Neighbor,
		pamhd::grid::Edge_Neighbor,
		pamhd::grid::Vertex_Neighbor,
		pamhd::grid::Relative_Size,
		pamhd::grid::Neighbor_Is_Local>
>;


/*
Creates evenly spaced particles in given cells.
*/
void create_particles(
	const size_t values_per_cell,
	const size_t dimension,
	Grid& grid
) {
	for (const auto& cell: grid.local_cells()) {
		const auto
			cell_min = grid.geometry.get_min(cell.id),
			cell_length = grid.geometry.get_length(cell.id),
			cell_center = grid.geometry.get_center(cell.id);

		const auto indices = grid.mapping.get_indices(cell.id);
		if (
			indices[dimension] > 0
			and indices[dimension] < grid.mapping.length.get()[dimension] - 1
		) {
			CType.data(*cell.data) = 1;
		} else {
			CType.data(*cell.data) = 0;
		}

		for (size_t i = 0; i < values_per_cell; i++) {
			pamhd::particle::Particle_Internal new_particle;

			auto& ppos = new_particle[pamhd::particle::Position()];
			ppos = cell_center;
			ppos[dimension]
				= cell_min[dimension]
				+ (double(i) + 0.5) * cell_length[dimension] / values_per_cell;

			const auto
				charge = charge_density(ppos),
				mass = charge,
				c2m = 1.0;
			// accumulate mass variable
			new_particle[pamhd::particle::Mass()] = mass / values_per_cell;
			new_particle[pamhd::particle::Charge_Mass_Ratio()] = c2m;

			const auto current = current_density(ppos);
			new_particle[pamhd::particle::Velocity()] = pamhd::mul(current, 1 / charge);

			(*cell.data)[pamhd::particle::Particles_Internal()].push_back(new_particle);
		}

		Nr_Int(*cell.data) = Part_Int(*cell.data).size();
	}
}


// returns infinite norm between analytic and grid results
double get_norm(
	const size_t dimension,
	const Grid& grid,
	MPI_Comm& comm
) {
	using std::abs;
	using std::max;

	double
		norm_local = 0,
		norm_global = 0;
	for (const auto& cell: grid.local_cells()) {
		if (CType.data(*cell.data) < 1) continue;

		const auto cell_center = grid.geometry.get_center(cell.id);
		const auto charge = charge_density(cell_center);
		const auto current = current_density(cell_center);
		const auto abs_qi = abs(Vol_Qi.data(*cell.data) - charge);
		const auto abs_ji = abs(Vol_Ji.data(*cell.data)[dimension] - current[dimension]);
		norm_local = max(max(norm_local, abs_qi), abs_ji);
	}
	if (
		MPI_Allreduce(
			&norm_local,
			&norm_global,
			1,
			MPI_DOUBLE,
			MPI_MAX,
			comm
		) != MPI_SUCCESS
	) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< ": Couldn't set reduce norm."
			<< std::endl;
		abort();
	}

	return norm_global;
}


int main(int argc, char* argv[])
{
	using std::array;

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

	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
		std::cerr << "Zoltan_Initialize failed." << std::endl;
		abort();
	}


	constexpr size_t
		nr_of_values = 10000,
		max_nr_of_cells = 640;

	// same as accumulate.cpp but for all dimensions
	double
		old_norm_x_np = std::numeric_limits<double>::max(),
		old_norm_y_np = std::numeric_limits<double>::max(),
		old_norm_z_np = std::numeric_limits<double>::max();
	size_t old_nr_of_cells = 0;
	for (size_t nr_of_cells = 10; nr_of_cells <= max_nr_of_cells; nr_of_cells *= 2) {

		const array<uint64_t, 3>
			nr_of_cells_x{nr_of_cells,           1,           1},
			nr_of_cells_y{          1, nr_of_cells,           1},
			nr_of_cells_z{          1,           1, nr_of_cells};

		const unsigned int neighborhood_size = 1;
		const array<bool, 3> non_periodic{false, false, false};

		Grid grid_x, grid_x_np, grid_y, grid_y_np, grid_z, grid_z_np;

		grid_x_np
			.set_neighborhood_length(neighborhood_size)
			.set_maximum_refinement_level(0)
			.set_load_balancing_method("RANDOM")
			.set_periodic(non_periodic[0], non_periodic[1], non_periodic[2])
			.set_initial_length(nr_of_cells_x)
			.initialize(comm)
			.balance_load();

		grid_y_np
			.set_neighborhood_length(neighborhood_size)
			.set_maximum_refinement_level(0)
			.set_load_balancing_method("RANDOM")
			.set_periodic(non_periodic[0], non_periodic[1], non_periodic[2])
			.set_initial_length(nr_of_cells_y)
			.initialize(comm)
			.balance_load();

		grid_z_np
			.set_neighborhood_length(neighborhood_size)
			.set_maximum_refinement_level(0)
			.set_load_balancing_method("RANDOM")
			.set_periodic(non_periodic[0], non_periodic[1], non_periodic[2])
			.set_initial_length(nr_of_cells_z)
			.initialize(comm)
			.balance_load();

		dccrg::Cartesian_Geometry::Parameters
			geom_params_x, geom_params_y, geom_params_z;
		geom_params_x.start = {-3,  0,  0};
		geom_params_y.start = { 0, -3,  0};
		geom_params_z.start = { 0,  0, -3};
		geom_params_x.level_0_cell_length = {9.0 / nr_of_cells_x[0],   1,   1};
		geom_params_y.level_0_cell_length = {  1, 9.0 / nr_of_cells_y[1],   1};
		geom_params_z.level_0_cell_length = {  1,   1, 9.0 / nr_of_cells_z[2]};

		grid_x_np.set_geometry(geom_params_x);
		grid_y_np.set_geometry(geom_params_y);
		grid_z_np.set_geometry(geom_params_z);

		const size_t values_per_cell = nr_of_values / nr_of_cells;
		create_particles(values_per_cell, 0, grid_x_np);
		create_particles(values_per_cell, 1, grid_y_np);
		create_particles(values_per_cell, 2, grid_z_np);

		Cell::set_transfer_all(true, pamhd::particle::Nr_Particles_Internal());
		grid_x_np.update_copies_of_remote_neighbors();
		grid_y_np.update_copies_of_remote_neighbors();
		grid_z_np.update_copies_of_remote_neighbors();
		Cell::set_transfer_all(false, pamhd::particle::Nr_Particles_Internal());

		for (const auto& cell: grid_x_np.all_cells()) {
			if (cell.is_local) continue;
			Part_Int(*cell.data).resize(Nr_Int(*cell.data));
		}
		for (const auto& cell: grid_y_np.all_cells()) {
			if (cell.is_local) continue;
			Part_Int(*cell.data).resize(Nr_Int(*cell.data));
		}
		for (const auto& cell: grid_z_np.all_cells()) {
			if (cell.is_local) continue;
			Part_Int(*cell.data).resize(Nr_Int(*cell.data));
		}

		Cell::set_transfer_all(true, pamhd::particle::Particles_Internal());
		grid_x_np.update_copies_of_remote_neighbors();
		grid_y_np.update_copies_of_remote_neighbors();
		grid_z_np.update_copies_of_remote_neighbors();
		Cell::set_transfer_all(false, pamhd::particle::Particles_Internal());

		accumulate_hyb(
			grid_x_np.local_cells(),
			grid_x_np,
			Part_Int,
			PPos,
			PVel,
			PMas,
			PSMas,
			PC2M,
			Vol_Qi,
			Vol_Ji,
			CType
		);
		accumulate_hyb(
			grid_y_np.local_cells(),
			grid_y_np,
			Part_Int,
			PPos,
			PVel,
			PMas,
			PSMas,
			PC2M,
			Vol_Qi,
			Vol_Ji,
			CType
		);
		accumulate_hyb(
			grid_z_np.local_cells(),
			grid_z_np,
			Part_Int,
			PPos,
			PVel,
			PMas,
			PSMas,
			PC2M,
			Vol_Qi,
			Vol_Ji,
			CType
		);

		Cell::set_transfer_all(true, Vol_Qi.type(), Vol_Ji.type());
		grid_x_np.update_copies_of_remote_neighbors();
		grid_y_np.update_copies_of_remote_neighbors();
		grid_z_np.update_copies_of_remote_neighbors();
		Cell::set_transfer_all(false, Vol_Qi.type(), Vol_Ji.type());

		const double
			norm_x_np = get_norm(0, grid_x_np, comm) / nr_of_cells,
			norm_y_np = get_norm(1, grid_y_np, comm) / nr_of_cells,
			norm_z_np = get_norm(2, grid_z_np, comm) / nr_of_cells;

		if (old_nr_of_cells > 0) {
			const double
				order_of_accuracy_x_np
					= -log(norm_x_np / old_norm_x_np)
					/ log(double(nr_of_cells) / old_nr_of_cells),
				order_of_accuracy_y_np
					= -log(norm_y_np / old_norm_y_np)
					/ log(double(nr_of_cells) / old_nr_of_cells),
				order_of_accuracy_z_np
					= -log(norm_z_np / old_norm_z_np)
					/ log(double(nr_of_cells) / old_nr_of_cells);

			if (order_of_accuracy_x_np < 2.9) {
				if (rank == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< ": Order of accuracy in x dimension from "
						<< old_nr_of_cells << " to " << nr_of_cells
						<< " is too low: " << order_of_accuracy_x_np
						<< std::endl;
				}
				MPI_Finalize();
				return EXIT_FAILURE;
			}

			if (order_of_accuracy_y_np < 2.9) {
				if (rank == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< ": Order of accuracy in y dimension from "
						<< old_nr_of_cells << " to " << nr_of_cells
						<< " is too low: " << order_of_accuracy_y_np
						<< std::endl;
				}
				MPI_Finalize();
				return EXIT_FAILURE;
			}

			if (order_of_accuracy_z_np < 2.9) {
				if (rank == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< ": Order of accuracy in z dimension from "
						<< old_nr_of_cells << " to " << nr_of_cells
						<< " is too low: " << order_of_accuracy_z_np
						<< std::endl;
				}
				MPI_Finalize();
				return EXIT_FAILURE;
			}
		}

		old_norm_x_np = norm_x_np;
		old_norm_y_np = norm_y_np;
		old_norm_z_np = norm_z_np;
		old_nr_of_cells = nr_of_cells;
	}

	MPI_Finalize();
	return EXIT_SUCCESS;
}
