/*
Tests vector field divergence removal of PAMHD in 3d.

Copyright 2014, 2015, 2016, 2017 Ilja Honkonen
Copyright 2018, 2019 Finnish Meteorological Institute
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
#include "cstdlib"
#include "iostream"
#include "iomanip"
#include "limits"
#include "vector"

#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#include "gensimcell.hpp"
#include "prettyprint.hpp"

#include "divergence/remove.hpp"


int Poisson_Cell::transfer_switch = Poisson_Cell::INIT;


std::array<double, 3> function(const std::array<double, 3>& r)
{
	return {0, 0, 1 + 0.1 * std::sin(r[2])};
}

double div_removed_function()
{
	return 1;
}


struct Vector_Field {
	using data_type = std::array<double, 3>;
};

struct Divergence {
	using data_type = double;
};

struct Gradient {
	using data_type = std::array<double, 3>;
};

struct Type {
	using data_type = int;
};

using Cell = gensimcell::Cell<
	gensimcell::Always_Transfer,
	Vector_Field,
	Divergence,
	Gradient,
	Type
>;


template<class Grid> double get_max_norm(const Grid& grid)
{
	double local_norm = 0, global_norm = 0;
	for (const auto& cell: grid.local_cells()) {
		const auto c = grid.geometry.get_center(cell.id);
		const auto l = grid.geometry.get_length(cell.id);
		if (
			std::abs(c[0] - M_PI) > l[0]
			or std::abs(c[1] - M_PI) > l[1]
			or std::abs(c[2] - M_PI) > l[2]
		) {
			continue;
		}

		local_norm = std::max(
			local_norm,
			std::fabs((*cell.data)[Vector_Field()][2] - div_removed_function())
		);
	}

	MPI_Comm comm = grid.get_communicator();
	MPI_Allreduce(&local_norm, &global_norm, 1, MPI_DOUBLE, MPI_MAX, comm);
	MPI_Comm_free(&comm);
	return global_norm;
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


	// intialize Zoltan
	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
		std::cerr << "Zoltan_Initialize failed." << std::endl;
		abort();
	}

	const unsigned int neighborhood_size = 0;
	const int max_refinement_level = 1;

	double old_norm = std::numeric_limits<double>::max();
	size_t old_nr_of_cells = 0;
	for (size_t nr_of_cells = 8; nr_of_cells <= 32; nr_of_cells *= 2) {

		dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry> grid;

		const std::array<uint64_t, 3> grid_size{{
			nr_of_cells + 2,
			nr_of_cells + 2,
			nr_of_cells + 2
		}};

		grid
			.set_load_balancing_method("RANDOM")
			.set_initial_length(grid_size)
			.set_maximum_refinement_level(max_refinement_level)
			.set_neighborhood_length(neighborhood_size)
			.initialize(comm)
			.balance_load();

		const std::array<double, 3>
			cell_length{{
				2 * M_PI / (grid_size[0] - 2),
				2 * M_PI / (grid_size[1] - 2),
				2 * M_PI / (grid_size[2] - 2)
			}},
			grid_start{{
				-cell_length[0], -cell_length[1], -cell_length[2]
			}};

		dccrg::Cartesian_Geometry::Parameters geom_params;
		geom_params.start = grid_start;
		geom_params.level_0_cell_length = cell_length;
		grid.set_geometry(geom_params);

		for (int i = 0; i < max_refinement_level; i++) {
			for (const auto& cell: grid.local_cells()) {
				const auto center = grid.geometry.get_center(cell.id);
				if (
					center[0] > 2.0*M_PI*3/8 and center[0] < 2.0*M_PI*6/8
					and center[1] > 2.0*M_PI*3/8 and center[1] < 2.0*M_PI*6/8
					and center[2] > 2.0*M_PI*3/8 and center[2] < 2.0*M_PI*6/8
				) {
					grid.refine_completely(cell.id);
				}
			}
			grid.stop_refining();
		}

		for (const auto& cell: grid.local_cells()) {
			const auto center = grid.geometry.get_center(cell.id);
			(*cell.data)[Vector_Field()] = function(center);
		}
		grid.update_copies_of_remote_neighbors();

		// classify cells
		uint64_t solve_cells_local = 0, solve_cells_global = 0;
		for (const auto& cell: grid.local_cells()) {
			const auto center = grid.geometry.get_center(cell.id);
			if (
				center[0] > 0 and center[0] < 2 * M_PI
				and center[1] > 0 and center[1] < 2 * M_PI
				and center[2] > 0 and center[2] < 2 * M_PI
			) {
				(*cell.data)[Type()] = 1;
				solve_cells_local++;
			} else {
				(*cell.data)[Type()] = 0;
			}
		}
		if (
			MPI_Allreduce(
				&solve_cells_local,
				&solve_cells_global,
				1,
				MPI_UINT64_T,
				MPI_SUM,
				comm
			) != MPI_SUCCESS
		) {
			std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}

		// apply copy boundaries
		for (const auto& cell: grid.local_cells()) {
			if ((*cell.data)[Type()] != 0) {
				continue;
			}

			const auto index = grid.mapping.get_indices(cell.id);
			auto neighbor_index = index;

			// length of ref lvl 0 cells in indices
			const auto init_cell_size = (1 << max_refinement_level);
			// assume cells close to boundaries haven't been refined
			if (index[0] == 0) {
				neighbor_index[0] = index[0] + init_cell_size;
			} else if (index[0] == init_cell_size * (grid_size[0] - 1)) {
				neighbor_index[0] = index[0] - init_cell_size;
			} else if (index[1] == 0) {
				neighbor_index[1] = index[1] + init_cell_size;
			} else if (index[1] == init_cell_size * (grid_size[1] - 1)) {
				neighbor_index[1] = index[1] - init_cell_size;
			} else if (index[2] == 0) {
				neighbor_index[2] = index[2] + init_cell_size;
			} else if (index[2] == init_cell_size * (grid_size[2] - 1)) {
				neighbor_index[2] = index[2] - init_cell_size;
			}
			const auto neighbor = grid.mapping.get_cell_from_indices(neighbor_index, 0);

			const auto* const neighbor_data = grid[neighbor];
			if (neighbor_data == nullptr) {
				std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}

			(*cell.data)[Vector_Field()] = (*neighbor_data)[Vector_Field()];
		}
		grid.update_copies_of_remote_neighbors();

		auto Vector_Getter = [](Cell& cell_data) -> Vector_Field::data_type& {
			return cell_data[Vector_Field()];
		};
		auto Div_Getter = [](Cell& cell_data) -> Divergence::data_type& {
			return cell_data[Divergence()];
		};
		auto Gradient_Getter = [](Cell& cell_data) -> Gradient::data_type& {
			return cell_data[Gradient()];
		};
		auto Type_Getter = [](Cell& cell_data) -> Type::data_type& {
			return cell_data[Type()];
		};
		const double div_before = pamhd::divergence::get_divergence(
			grid.local_cells(),
			grid,
			Vector_Getter,
			Div_Getter,
			Type_Getter
		);

		pamhd::divergence::remove(
			grid.local_cells(),
			grid,
			Vector_Getter,
			Div_Getter,
			Gradient_Getter,
			Type_Getter,
			2000, 0, 1e-10, 2, 100, 10, false, false
		);
		grid.update_copies_of_remote_neighbors();

		const double norm = get_max_norm(grid);
		// update copy boundaries to correspond to removed divergence
		for (const auto& cell: grid.local_cells()) {
			if ((*cell.data)[Type()] != 0) {
				continue;
			}

			const auto index = grid.mapping.get_indices(cell.id);
			auto neighbor_index = index;

			const auto init_cell_size = (1 << max_refinement_level);
			if (index[0] == 0) {
				neighbor_index[0] = index[0] + init_cell_size;
			} else if (index[0] == init_cell_size * (grid_size[0] - 1)) {
				neighbor_index[0] = index[0] - init_cell_size;
			} else if (index[1] == 0) {
				neighbor_index[1] = index[1] + init_cell_size;
			} else if (index[1] == init_cell_size * (grid_size[1] - 1)) {
				neighbor_index[1] = index[1] - init_cell_size;
			} else if (index[2] == 0) {
				neighbor_index[2] = index[2] + init_cell_size;
			} else if (index[2] == init_cell_size * (grid_size[2] - 1)) {
				neighbor_index[2] = index[2] - init_cell_size;
			}
			const auto neighbor = grid.mapping.get_cell_from_indices(neighbor_index, 0);

			const auto* const neighbor_data = grid[neighbor];
			if (neighbor_data == nullptr) {
				std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}

			(*cell.data)[Vector_Field()] = (*neighbor_data)[Vector_Field()];
		}
		grid.update_copies_of_remote_neighbors();

		const double div_after = pamhd::divergence::get_divergence(
			grid.local_cells(),
			grid,
			Vector_Getter,
			Div_Getter,
			Type_Getter
		);

		if (div_after > div_before) {
			if (grid.get_rank() == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< ": Divergence after removal " << div_after
					<< " is larger than before " << div_before
					<< " with " << solve_cells_global << " cells."
					<< std::endl;
			}
			abort();
		}

		if (old_nr_of_cells > 0) {
			const double
				order_of_accuracy
					= -log(norm / old_norm)
					/ log(double(solve_cells_global) / old_nr_of_cells);

			if (order_of_accuracy < 0.3) {
				if (grid.get_rank() == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< ": Order of accuracy from "
						<< old_nr_of_cells << " to " << solve_cells_global
						<< " is too low: " << order_of_accuracy
						<< "  " << old_norm << "->" << norm
						<< std::endl;
				}
				abort();
			}
		}

		old_nr_of_cells = solve_cells_global;
		old_norm = norm;
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
