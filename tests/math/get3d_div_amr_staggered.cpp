/*
Tests vector field divergence calculation of PAMHD in 3d.

Copyright 2014, 2015, 2016, 2017 Ilja Honkonen
Copyright 2018, 2022, 2023 Finnish Meteorological Institute
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
#include "limits"
#include "vector"

#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#include "gensimcell.hpp"

#include "grid/variables.hpp"
#include "math/staggered.hpp"
#include "tests/math/common.hpp"


std::array<double, 3> function(const std::array<double, 3>& r)
{
	const double x = r[0], y = r[1], z = r[2];

	return {{
		std::exp(2*x) + 5*x + 6,
		std::sin(y) * std::cos(y) + 4*y + 5,
		3*z*z*z - 2*z*z + z - 2
	}};
}

double div_of_function(const std::array<double, 3>& r)
{
	const double x = r[0], y = r[1], z = r[2];

	return 9*z*z - 4*z + 2 * std::pow(std::cos(y), 2) + 2 * std::exp(2*x) + 9;
}


/*!
Returns maximum norm.
*/
template<class Grid> double get_max_norm(
	const Grid& grid
) {
	double local_norm = 0, global_norm = 0;
	for (const auto& cell: grid.local_cells()) {
		if ((*cell.data)[Type()] != 1) {
			continue;
		}

		const auto div_of = div_of_function(grid.geometry.get_center(cell.id));

		local_norm = std::max(
			local_norm,
			std::fabs((*cell.data)[Divergence()] - div_of)
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

		dccrg::Dccrg<
			Cell,
			dccrg::Cartesian_Geometry,
			std::tuple<>,
			std::tuple<
				pamhd::grid::Is_Face_Neighbor,
				pamhd::grid::Relative_Size>
		> grid;

		const std::array<uint64_t, 3> grid_size{{nr_of_cells + 2, nr_of_cells + 2, nr_of_cells + 2}};

		grid
			.set_load_balancing_method("RANDOM")
			.set_initial_length(grid_size)
			.set_neighborhood_length(neighborhood_size)
			.set_maximum_refinement_level(max_refinement_level)
			.initialize(comm)
			.balance_load();

		const std::array<double, 3>
			cell_length{{
				double(3) / (grid_size[0] - 2),
				1.5 / (grid_size[1] - 2),
				double(4) / (grid_size[2] - 2)
			}},
			grid_start{{
				-1 - cell_length[0],
				-M_PI / 4 - cell_length[1],
				-2 - cell_length[2]
			}};

		dccrg::Cartesian_Geometry::Parameters geom_params;
		geom_params.start = grid_start;
		geom_params.level_0_cell_length = cell_length;
		grid.set_geometry(geom_params);

		for (int i = 0; i < max_refinement_level; i++) {
			for (const auto& cell: grid.local_cells()) {
				const auto center = grid.geometry.get_center(cell.id);
				if (
					center[0] > grid_start[0] + 3.0 * 1 / 4
					and center[0] < grid_start[0] + 3.0 * 3 / 4
					and center[1] > grid_start[1] + 1.5 * 1 / 4
					and center[1] < grid_start[1] + 1.5 * 3 / 4
					and center[2] > grid_start[2] + 4.0 * 1 / 4
					and center[2] < grid_start[2] + 4.0 * 3 / 4
				) {
					grid.refine_completely(cell.id);
				}
			}
			grid.stop_refining();
		}

		for (const auto& cell: grid.local_cells()) {
			const auto
				center = grid.geometry.get_center(cell.id),
				length = grid.geometry.get_length(cell.id);
			(*cell.data)[Vector_Pos()] = function({
				center[0] + length[0]/2,
				center[1] + length[1]/2,
				center[2] + length[2]/2
			});
			(*cell.data)[Vector_Neg()] = function({
				center[0] - length[0]/2,
				center[1] - length[1]/2,
				center[2] - length[2]/2
			});
		}
		grid.update_copies_of_remote_neighbors();

		uint64_t solve_cells_local = 0, solve_cells_global = 0;
		for (const auto& cell: grid.local_cells()) {
			const auto center = grid.geometry.get_center(cell.id);
			if (
				center[0] > -1 and center[0] < -1 + 3
				and center[1] > -M_PI / 4 and center[1] < -M_PI / 4 + 1.5
				and center[2] > -2 and center[2] < -2 + 4
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

		pamhd::math::get_divergence_staggered(
			grid.local_cells(),
			grid,
			[](Cell& cell_data)->auto& {
				return cell_data[Vector_Pos()];
			},
			[](Cell& cell_data)->auto& {
				return cell_data[Vector_Neg()];
			},
			[](Cell& cell_data)->auto& {
				return cell_data[Divergence()];
			},
			[](Cell& cell_data)->auto& {
				return cell_data[Type()];
			}
		);

		const double norm = get_max_norm(grid);

		if (norm > old_norm) {
			if (grid.get_rank() == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< ": Norm with " << solve_cells_global
					<< " cells " << norm
					<< " is larger than with " << old_nr_of_cells
					<< " cells " << old_norm
					<< std::endl;
			}
			abort();
		}

		if (old_nr_of_cells > 0) {
			const double order_of_accuracy
				= -log(norm / old_norm)
				/ log(double(solve_cells_global) / old_nr_of_cells);

			if (order_of_accuracy < 0.5) {
				if (grid.get_rank() == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< ": Order of accuracy from "
						<< old_nr_of_cells << " to " << solve_cells_global
						<< " is too low: " << order_of_accuracy
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
