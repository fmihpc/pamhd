/*
vertex2volume test in non-periodic grid.

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
#include "cstdlib"
#include "iostream"
#include "limits"
#include "utility"

#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#include "gensimcell.hpp"

#include "grid/variables.hpp"
#include "math/interpolation.hpp"
#include "tests/math/common.hpp"


double function(const std::array<double, 3>& r)
{
	return -std::sin(r[0] / 2) * std::cos(2 * r[1]) / (r[2] + 10);
}

// target cannot be vertex type
const auto Tgt = pamhd::Variable_Getter<Source>();
bool Source::is_stale = true;

// source must be vertex type
const auto Src = pamhd::Variable_Getter<Target>();
bool Target::is_stale = true;

const auto Typ = pamhd::Variable_Getter<Type>();
bool Type::is_stale = true;

/*!
Returns maximum norm if p == 0
*/
template<class Grid> double get_diff_lp_norm(
	const Grid& grid, const double cell_volume
) {
	double local_norm = 0, global_norm = 0;
	for (const auto& cell: grid.local_cells()) {
		if (Typ.data(*cell.data) < 1) {
			continue;
		}

		const auto center = grid.geometry.get_center(cell.id);
		local_norm = std::max(
			local_norm,
			std::fabs(Tgt.data(*cell.data) - function(center))
		);
	}
	local_norm *= cell_volume;

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

	const unsigned int neighborhood_size = 1;
	const int max_refinement_level = 0;

	double old_norm = std::numeric_limits<double>::max();
	size_t old_nr_of_cells = 0;
	for (size_t nr_of_cells = 8; nr_of_cells <= 64; nr_of_cells *= 2) {

		dccrg::Dccrg<
			Cell_interp,
			dccrg::Cartesian_Geometry,
			std::tuple<>,
			std::tuple<
				pamhd::grid::Face_Neighbor,
				pamhd::grid::Edge_Neighbor,
				pamhd::grid::Vertex_Neighbor>
		> grid;

		const std::array<uint64_t, 3> grid_size{
			nr_of_cells + 2, nr_of_cells + 2, nr_of_cells + 2};

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

		const double cell_volume
			= cell_length[0] * cell_length[1] * cell_length[2];

		dccrg::Cartesian_Geometry::Parameters geom_params;
		geom_params.start = grid_start;
		geom_params.level_0_cell_length = cell_length;
		grid.set_geometry(geom_params);

		for (const auto& cell: grid.local_cells()) {
			const auto
				center = grid.geometry.get_center(cell.id),
				length = grid.geometry.get_length(cell.id);
			for (int dx: {-1, +1})
			for (int dy: {-1, +1})
			for (int dz: {-1, +1}) {
				const std::array<double, 3> r{
					center[0] + 0.5*dx*length[0],
					center[1] + 0.5*dy*length[1],
					center[2] + 0.5*dz*length[2]
				};
				Src.data(*cell.data)(dx,dy,dz) = function(r);
			}
		}
		Src.type().is_stale = true;

		uint64_t solve_cells_local = 0, solve_cells_global = 0;
		for (const auto& cell: grid.local_cells()) {
			const auto index = grid.mapping.get_indices(cell.id);
			if (
				index[0] > 1
				and index[0] < grid_size[0] - 2
				and index[1] > 1
				and index[1] < grid_size[1] - 2
				and index[2] > 1
				and index[2] < grid_size[2] - 2
			) {
				Typ.data(*cell.data) = 1;
				solve_cells_local++;
			} else if (
				index[0] > 0
				and index[0] < grid_size[0] - 1
				and index[1] > 0
				and index[1] < grid_size[1] - 1
				and index[2] > 0
				and index[2] < grid_size[2] - 1
			) {
				Typ.data(*cell.data) = 0;
			} else {
				Typ.data(*cell.data) = -1;
			}
		}
		if (
			MPI_Allreduce(
				&solve_cells_local, &solve_cells_global,
				1, MPI_UINT64_T, MPI_SUM, comm
			) != MPI_SUCCESS
		) {
			std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}

		pamhd::math::vertex2volume(grid, Src, Tgt, Typ);

		const double norm = get_diff_lp_norm(grid, cell_volume);

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

			if (order_of_accuracy < 1.2) {
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
