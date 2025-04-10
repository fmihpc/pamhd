/*
Test for particle accumulator of PAMHD built on top of DCCRG.

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

#include "algorithm"
#include "cmath"
#include "cstdlib"
#include "exception"
#include "iostream"
#include "utility"
#include "vector"

#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#include "Eigen/Core" // must be included before gensimcell.hpp
#include "mpi.h" // must be included before gensimcell.hpp
#include "gensimcell.hpp"

#include "particle/accumulate_dccrg.hpp"
#include "particle/accumulation_variables.hpp"
#include "particle/variables.hpp"
#include "variable_getter.hpp"


Eigen::Vector3d f(const Eigen::Vector3d& r)
{
	return {
		1.5 * sin(2 * M_PI * (r[0] + 3) / 9),
		1.5 * sin(2 * M_PI * (r[1] + 3) / 9),
		1.5 * sin(2 * M_PI * (r[2] + 3) / 9)
	};
}

//! average number of particles in a cell
struct Count {
	using data_type = double;
};

using Accumulated_To_Cell = pamhd::particle::Accumulated_To_Cell_T<Count>;
using Accumulated_To_Cells = pamhd::particle::Accumulated_To_Cells_T<Accumulated_To_Cell>;

using Cell = gensimcell::Cell<
	gensimcell::Optional_Transfer,
	pamhd::particle::Particles_Internal,
	Count,
	pamhd::Solver_Info,
	// only following two are transferred between processes
	pamhd::particle::Nr_Accumulated_To_Cells,
	Accumulated_To_Cells
>;
using Grid = dccrg::Dccrg<
	Cell,
	dccrg::Cartesian_Geometry,
	std::tuple<>,
	std::tuple<pamhd::particle::Is_Local>
>;


// returns a reference to data accumulated from particles in given cell
const auto bulk_value_getter = pamhd::Variable_Getter<Count>();

// returns a reference to data accumulated from particles in accumulation list
const auto list_bulk_value_getter
	= [](Accumulated_To_Cell& accu_item)->Count::data_type&{
		return accu_item[Count()];
	};

// returns a reference to target of accumulated values in accumulation list
const auto list_target_getter
	= [](Accumulated_To_Cell& accu_item)->pamhd::particle::Target::data_type&{
		return accu_item[pamhd::particle::Target()];
	};

/*
Returns a reference to given cell's list of data accumulated
from particles in given cell into other cells
*/
const auto accumulation_list_getter
	= [](Cell& cell_data)->Accumulated_To_Cells::data_type&{
		return cell_data[Accumulated_To_Cells()];
	};

// returns a reference to length of list above incoming from other processes
const auto accumulation_list_length_getter
	= [](Cell& cell_data)->pamhd::particle::Nr_Accumulated_To_Cells::data_type&{
		return cell_data[pamhd::particle::Nr_Accumulated_To_Cells()];
	};

const auto solver_info_getter = pamhd::Variable_Getter<pamhd::Solver_Info>();

const auto accumulate_from_remote_neighbors
	= [](Grid& grid){
		pamhd::particle::accumulate_from_remote_neighbors(
			grid,
			bulk_value_getter,
			list_bulk_value_getter,
			list_target_getter,
			accumulation_list_getter,
			solver_info_getter
		);
	};

const auto allocate_accumulation_lists
	= [](Grid& grid){
		pamhd::particle::allocate_accumulation_lists(
			grid,
			accumulation_list_getter,
			accumulation_list_length_getter
		);
	};


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

		solver_info_getter.data(*cell.data) = 0;

		for (size_t i = 0; i < values_per_cell; i++) {
			pamhd::particle::Particle_Internal new_particle;
			new_particle[pamhd::particle::Velocity()][0] =
			new_particle[pamhd::particle::Velocity()][1] =
			new_particle[pamhd::particle::Velocity()][2] =
			new_particle[pamhd::particle::Charge_Mass_Ratio()] = 0;

			new_particle[pamhd::particle::Position()][0] = cell_center[0];
			new_particle[pamhd::particle::Position()][1] = cell_center[1];
			new_particle[pamhd::particle::Position()][2] = cell_center[2];

			new_particle[pamhd::particle::Position()][dimension]
				= cell_min[dimension]
				+ (double(i) + 0.5) * cell_length[dimension] / values_per_cell;

			// accumulate mass variable
			new_particle[pamhd::particle::Mass()]
				= f(new_particle[pamhd::particle::Position()])[dimension]
				/ values_per_cell;

			(*cell.data)[pamhd::particle::Particles_Internal()].push_back(new_particle);
		}
	}
}


// returns infinite norm between analytic results and that in given cells
double get_norm(
	const size_t dimension,
	const Grid& grid,
	MPI_Comm& comm
) {
	double
		norm_local = 0,
		norm_global = 0;
	for (const auto& cell: grid.local_cells()) {
		const Eigen::Vector3d cell_center{
			grid.geometry.get_center(cell.id)[0],
			grid.geometry.get_center(cell.id)[1],
			grid.geometry.get_center(cell.id)[2]
		};

		norm_local
			= std::max(
				norm_local,
				std::fabs(
					(*cell.data)[Count()]
					- f(cell_center)[dimension]
				)
			);
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
		old_norm_x = std::numeric_limits<double>::max(),
		old_norm_x_np = std::numeric_limits<double>::max(),
		old_norm_y = std::numeric_limits<double>::max(),
		old_norm_y_np = std::numeric_limits<double>::max(),
		old_norm_z = std::numeric_limits<double>::max(),
		old_norm_z_np = std::numeric_limits<double>::max();
	size_t old_nr_of_cells = 0;
	for (size_t nr_of_cells = 10; nr_of_cells <= max_nr_of_cells; nr_of_cells *= 2) {

		const std::array<uint64_t, 3>
			nr_of_cells_x{{nr_of_cells,           1,           1}},
			nr_of_cells_y{{          1, nr_of_cells,           1}},
			nr_of_cells_z{{          1,           1, nr_of_cells}};

		const unsigned int neighborhood_size = 1;
		const std::array<bool, 3>
			periodic_x = {{true, false, false}},
			periodic_y = {{false, true, false}},
			periodic_z = {{false, false, true}},
			non_periodic = {{false, false, false}};

		Grid grid_x, grid_x_np, grid_y, grid_y_np, grid_z, grid_z_np;

		grid_x
			.set_neighborhood_length(neighborhood_size)
			.set_maximum_refinement_level(0)
			.set_load_balancing_method("RANDOM")
			.set_periodic(periodic_x[0], periodic_x[1], periodic_x[2])
			.set_initial_length(nr_of_cells_x)
			.initialize(comm)
			.balance_load();

		grid_x_np
			.set_neighborhood_length(neighborhood_size)
			.set_maximum_refinement_level(0)
			.set_load_balancing_method("RANDOM")
			.set_periodic(non_periodic[0], non_periodic[1], non_periodic[2])
			.set_initial_length(nr_of_cells_x)
			.initialize(comm)
			.balance_load();

		grid_y
			.set_neighborhood_length(neighborhood_size)
			.set_maximum_refinement_level(0)
			.set_load_balancing_method("RANDOM")
			.set_periodic(periodic_y[0], periodic_y[1], periodic_y[2])
			.set_initial_length(nr_of_cells_y)
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

		grid_z
			.set_neighborhood_length(neighborhood_size)
			.set_maximum_refinement_level(0)
			.set_load_balancing_method("RANDOM")
			.set_periodic(periodic_z[0], periodic_z[1], periodic_z[2])
			.set_initial_length(nr_of_cells_z)
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
		geom_params_x.start = {{-3,  0,  0}};
		geom_params_y.start = {{ 0, -3,  0}};
		geom_params_z.start = {{ 0,  0, -3}};
		geom_params_x.level_0_cell_length = {{9.0 / nr_of_cells_x[0],   1,   1}};
		geom_params_y.level_0_cell_length = {{  1, 9.0 / nr_of_cells_y[1],   1}};
		geom_params_z.level_0_cell_length = {{  1,   1, 9.0 / nr_of_cells_z[2]}};

		grid_x.set_geometry(geom_params_x);
		grid_x_np.set_geometry(geom_params_x);
		grid_y.set_geometry(geom_params_y);
		grid_y_np.set_geometry(geom_params_y);
		grid_z.set_geometry(geom_params_z);
		grid_z_np.set_geometry(geom_params_z);

		// create copies of remote neighbors
		grid_x.update_copies_of_remote_neighbors();
		grid_x_np.update_copies_of_remote_neighbors();
		grid_y.update_copies_of_remote_neighbors();
		grid_y_np.update_copies_of_remote_neighbors();
		grid_z.update_copies_of_remote_neighbors();
		grid_z_np.update_copies_of_remote_neighbors();

		const size_t values_per_cell = nr_of_values / nr_of_cells;
		create_particles(values_per_cell, 0, grid_x);
		create_particles(values_per_cell, 0, grid_x_np);
		create_particles(values_per_cell, 1, grid_y);
		create_particles(values_per_cell, 1, grid_y_np);
		create_particles(values_per_cell, 2, grid_z);
		create_particles(values_per_cell, 2, grid_z_np);

		const auto accumulate_particles
			= [](Grid& grid) {
				pamhd::particle::accumulate(
					grid.local_cells(),
					grid,
					[](Cell& cell)->pamhd::particle::Particles_Internal::data_type&{
						return cell[pamhd::particle::Particles_Internal()];
					},
					[](pamhd::particle::Particle_Internal& particle)
						->pamhd::particle::Position::data_type&
					{
						return particle[pamhd::particle::Position()];
					},
					[](Cell&, pamhd::particle::Particle_Internal& particle)
						->pamhd::particle::Mass::data_type&
					{
						return particle[pamhd::particle::Mass()];
					},
					bulk_value_getter,
					list_bulk_value_getter,
					list_target_getter,
					accumulation_list_length_getter,
					accumulation_list_getter,
					solver_info_getter,
					true
				);
			};

		accumulate_particles(grid_x);
		accumulate_particles(grid_x_np);
		accumulate_particles(grid_y);
		accumulate_particles(grid_y_np);
		accumulate_particles(grid_z);
		accumulate_particles(grid_z_np);

		// transfer number accumulated values between processes
		Cell::set_transfer_all(true, pamhd::particle::Nr_Accumulated_To_Cells());
		grid_x.update_copies_of_remote_neighbors();
		grid_x_np.update_copies_of_remote_neighbors();
		grid_y.update_copies_of_remote_neighbors();
		grid_y_np.update_copies_of_remote_neighbors();
		grid_z.update_copies_of_remote_neighbors();
		grid_z_np.update_copies_of_remote_neighbors();
		Cell::set_transfer_all(false, pamhd::particle::Nr_Accumulated_To_Cells());

		allocate_accumulation_lists(grid_x);
		allocate_accumulation_lists(grid_x_np);
		allocate_accumulation_lists(grid_y);
		allocate_accumulation_lists(grid_y_np);
		allocate_accumulation_lists(grid_z);
		allocate_accumulation_lists(grid_z_np);

		// transfer accumulated values between processes
		Cell::set_transfer_all(true, Accumulated_To_Cells());
		grid_x.update_copies_of_remote_neighbors();
		grid_x_np.update_copies_of_remote_neighbors();
		grid_y.update_copies_of_remote_neighbors();
		grid_y_np.update_copies_of_remote_neighbors();
		grid_z.update_copies_of_remote_neighbors();
		grid_z_np.update_copies_of_remote_neighbors();
		Cell::set_transfer_all(false, Accumulated_To_Cells());

		accumulate_from_remote_neighbors(grid_x);
		accumulate_from_remote_neighbors(grid_x_np);
		accumulate_from_remote_neighbors(grid_y);
		accumulate_from_remote_neighbors(grid_y_np);
		accumulate_from_remote_neighbors(grid_z);
		accumulate_from_remote_neighbors(grid_z_np);

		const double
			norm_x = get_norm(0, grid_x, comm) / nr_of_cells,
			norm_x_np = get_norm(0, grid_x_np, comm) / nr_of_cells,
			norm_y = get_norm(1, grid_y, comm) / nr_of_cells,
			norm_y_np = get_norm(1, grid_y_np, comm) / nr_of_cells,
			norm_z = get_norm(2, grid_z, comm) / nr_of_cells,
			norm_z_np = get_norm(2, grid_z_np, comm) / nr_of_cells;

		if (old_nr_of_cells > 0) {
			const double
				order_of_accuracy_x
					= -log(norm_x / old_norm_x)
					/ log(double(nr_of_cells) / old_nr_of_cells),
				order_of_accuracy_x_np
					= -log(norm_x_np / old_norm_x_np)
					/ log(double(nr_of_cells) / old_nr_of_cells),
				order_of_accuracy_y
					= -log(norm_y / old_norm_y)
					/ log(double(nr_of_cells) / old_nr_of_cells),
				order_of_accuracy_y_np
					= -log(norm_y_np / old_norm_y_np)
					/ log(double(nr_of_cells) / old_nr_of_cells),
				order_of_accuracy_z
					= -log(norm_z / old_norm_z)
					/ log(double(nr_of_cells) / old_nr_of_cells),
				order_of_accuracy_z_np
					= -log(norm_z_np / old_norm_z_np)
					/ log(double(nr_of_cells) / old_nr_of_cells);

			if (order_of_accuracy_x < 1.9) {
				if (rank == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< ": Order of accuracy in x dimension from "
						<< old_nr_of_cells << " to " << nr_of_cells
						<< " is too low: " << order_of_accuracy_x
						<< std::endl;
				}
				MPI_Finalize();
				return EXIT_FAILURE;
			}
			if (order_of_accuracy_x_np < 1.9) {
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

			if (order_of_accuracy_y < 1.9) {
				if (rank == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< ": Order of accuracy in y dimension from "
						<< old_nr_of_cells << " to " << nr_of_cells
						<< " is too low: " << order_of_accuracy_y
						<< std::endl;
				}
				MPI_Finalize();
				return EXIT_FAILURE;
			}
			if (order_of_accuracy_y_np < 1.9) {
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

			if (order_of_accuracy_z < 1.9) {
				if (rank == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< ": Order of accuracy in z dimension from "
						<< old_nr_of_cells << " to " << nr_of_cells
						<< " is too low: " << order_of_accuracy_z
						<< std::endl;
				}
				MPI_Finalize();
				return EXIT_FAILURE;
			}
			if (order_of_accuracy_z_np < 1.9) {
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

		old_norm_x    = norm_x;
		old_norm_x_np = norm_x_np;
		old_norm_y    = norm_y;
		old_norm_y_np = norm_y_np;
		old_norm_z    = norm_z;
		old_norm_z_np = norm_z_np;
		old_nr_of_cells = nr_of_cells;
	}

	MPI_Finalize();
	return EXIT_SUCCESS;
}
