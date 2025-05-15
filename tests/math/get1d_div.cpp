/*
Staggered version of get1d_div.cpp.

Copyright 2014, 2015, 2016, 2017 Ilja Honkonen
Copyright 2018, 2022, 2023,
          2024, 2025 Finnish Meteorological Institute
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
#include "vector"

#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#include "gensimcell.hpp"

#include "common_variables.hpp"
#include "grid/variables.hpp"
#include "math/staggered.hpp"
#include "tests/math/common.hpp"
#include "variable_getter.hpp"


double function(const double x)
{
	return std::pow(std::sin(x / 2), 2);
}

double div_of_function(const double x)
{
	return std::sin(x / 2) * std::cos(x / 2);
}

/*!
Returns maximum norm if p == 0
*/
template<class Grid> double get_diff_lp_norm(
	const Grid& grid,
	const double p,
	const double cell_volume,
	const size_t dimension
) {
	double local_norm = 0, global_norm = 0;
	for (const auto& cell: grid.local_cells()) {
		if ((*cell.data)[Type()] != 1) {
			continue;
		}

		const auto center = grid.geometry.get_center(cell.id);

		if (p == 0) {
			local_norm = std::max(
				local_norm,
				std::fabs((*cell.data)[Divergence()] - div_of_function(center[dimension]))
			);
		} else {
			local_norm += std::pow(
				std::fabs((*cell.data)[Divergence()] - div_of_function(center[dimension])),
				p
			);
		}
	}
	local_norm *= cell_volume;

	if (p == 0) {
		MPI_Comm comm = grid.get_communicator();
		MPI_Allreduce(&local_norm, &global_norm, 1, MPI_DOUBLE, MPI_MAX, comm);
		MPI_Comm_free(&comm);
		return global_norm;
	} else {
		MPI_Comm comm = grid.get_communicator();
		MPI_Allreduce(&local_norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM, comm);
		MPI_Comm_free(&comm);
		return std::pow(global_norm, 1.0 / p);
	}
}


template<
	class Vector,
	class Type,
	class Grid
> void initialize(
	Grid& grid,
	MPI_Comm& comm,
	const uint64_t nr_of_cells,
	const unsigned int dimension
) {
	std::array<uint64_t, 3> grid_size;
	std::array<double, 3> cell_length, grid_start;
	switch (dimension) {
	case 0:
		grid_size = {nr_of_cells + 2, 1, 1};
		cell_length = {2 * M_PI / (grid_size[dimension] - 2), 1, 1};
		grid_start = {-cell_length[dimension], 0, 0};
		break;
	case 1:
		grid_size = {1, nr_of_cells + 2, 1};
		cell_length = {1, 2 * M_PI / (grid_size[dimension] - 2), 1};
		grid_start = {0, -cell_length[dimension], 0};
		break;
	case 2:
		grid_size = {1, 1, nr_of_cells + 2};
		cell_length = {1, 1, 2 * M_PI / (grid_size[dimension] - 2)};
		grid_start = {0, 0, -cell_length[dimension]};
		break;
	default:
		abort();
		break;
	}
	grid
		.set_load_balancing_method("RANDOM")
		.set_initial_length(grid_size)
		.set_neighborhood_length(0)
		.set_maximum_refinement_level(0)
		.initialize(comm);

	dccrg::Cartesian_Geometry::Parameters geom_params;
	geom_params.start = grid_start;
	geom_params.level_0_cell_length = cell_length;
	grid.set_geometry(geom_params);

	for (const auto& cell: grid.local_cells()) {
		const auto
			center = grid.geometry.get_center(cell.id),
			length = grid.geometry.get_length(cell.id);

		auto& vec = (*cell.data)[Vector()];
		vec = {0, 0, 0, 0, 0, 0};
		vec(dimension, +1) = function(center[dimension] + length[dimension]/2);
		vec(dimension, -1) = function(center[dimension] - length[dimension]/2);

		// exclude one layer of boundary cells
		const auto index = grid.mapping.get_indices(cell.id);
		if (index[dimension] > 0 and index[dimension] < grid_size[dimension] - 1) {
			(*cell.data)[Type()] = 1;
		} else {
			(*cell.data)[Type()] = 0;
		}
	}
	grid.update_copies_of_remote_neighbors();
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

	double
		old_norm_x = std::numeric_limits<double>::max(),
		old_norm_y = std::numeric_limits<double>::max(),
		old_norm_z = std::numeric_limits<double>::max();

	size_t old_nr_of_cells = 0;
	for (size_t nr_of_cells = 8; nr_of_cells <= 4096; nr_of_cells *= 2) {

		dccrg::Dccrg<
			Cell,
			dccrg::Cartesian_Geometry,
			std::tuple<>,
			std::tuple<
				pamhd::grid::Face_Neighbor,
				pamhd::grid::Relative_Size>
		> grid_x, grid_y, grid_z;

		initialize<Vector, Type>(grid_x, comm, nr_of_cells, 0);
		initialize<Vector, Type>(grid_y, comm, nr_of_cells, 1);
		initialize<Vector, Type>(grid_z, comm, nr_of_cells, 2);

		const auto Vector_Getter = pamhd::Variable_Getter<Vector>();
		const auto Divergence_Getter = pamhd::Variable_Getter<Divergence>();
		const auto Type_Getter = pamhd::Variable_Getter<Type>();
		pamhd::math::get_divergence_face2volume(
			grid_x.local_cells(),
			grid_x,
			Vector_Getter,
			Divergence_Getter,
			Type_Getter
		);
		pamhd::math::get_divergence_face2volume(
			grid_y.local_cells(),
			grid_y,
			Vector_Getter,
			Divergence_Getter,
			Type_Getter
		);
		pamhd::math::get_divergence_face2volume(
			grid_z.local_cells(),
			grid_z,
			Vector_Getter,
			Divergence_Getter,
			Type_Getter
		);

		const auto cell_length = grid_x.geometry.get_level_0_cell_length();
		const double
			cell_volume = cell_length[0] * cell_length[1] * cell_length[2],
			p_of_norm = 2,
			norm_x = get_diff_lp_norm(grid_x, p_of_norm, cell_volume, 0),
			norm_y = get_diff_lp_norm(grid_y, p_of_norm, cell_volume, 1),
			norm_z = get_diff_lp_norm(grid_z, p_of_norm, cell_volume, 2);

		if (norm_x > old_norm_x) {
			if (grid_x.get_rank() == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< ": X norm with " << nr_of_cells
					<< " cells " << norm_x
					<< " is larger than with " << nr_of_cells / 2
					<< " cells " << old_norm_x
					<< std::endl;
			}
			abort();
		}
		if (norm_y > old_norm_y) {
			if (grid_y.get_rank() == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< ": Y norm with " << nr_of_cells
					<< " cells " << norm_x
					<< " is larger than with " << nr_of_cells / 2
					<< " cells " << old_norm_x
					<< std::endl;
			}
			abort();
		}
		if (norm_z > old_norm_z) {
			if (grid_z.get_rank() == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< ": Z norm with " << nr_of_cells
					<< " cells " << norm_x
					<< " is larger than with " << nr_of_cells / 2
					<< " cells " << old_norm_x
					<< std::endl;
			}
			abort();
		}

		if (old_nr_of_cells > 0) {
			const double
				order_of_accuracy_x
					= -log(norm_x / old_norm_x)
					/ log(double(nr_of_cells) / old_nr_of_cells),
				order_of_accuracy_y
					= -log(norm_y / old_norm_y)
					/ log(double(nr_of_cells) / old_nr_of_cells),
				order_of_accuracy_z
					= -log(norm_z / old_norm_z)
					/ log(double(nr_of_cells) / old_nr_of_cells);

			if (order_of_accuracy_x < 1.99) {
				if (grid_x.get_rank() == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< ": X order of accuracy from "
						<< old_nr_of_cells << " to " << nr_of_cells
						<< " is too low: " << order_of_accuracy_x
						<< std::endl;
				}
				abort();
			}
			if (order_of_accuracy_y < 1.99) {
				if (grid_y.get_rank() == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< ": Y order of accuracy from "
						<< old_nr_of_cells << " to " << nr_of_cells
						<< " is too low: " << order_of_accuracy_y
						<< " (" << norm_y << ", " << old_norm_y
						<< ", " << nr_of_cells << ", " << old_nr_of_cells
						<< ")" << std::endl;
				}
				abort();
			}
			if (order_of_accuracy_z < 1.99) {
				if (grid_z.get_rank() == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< ": Z order of accuracy from "
						<< old_nr_of_cells << " to " << nr_of_cells
						<< " is too low: " << order_of_accuracy_z
						<< std::endl;
				}
				abort();
			}
		}

		old_nr_of_cells = nr_of_cells;
		old_norm_x = norm_x;
		old_norm_y = norm_y;
		old_norm_z = norm_z;
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
