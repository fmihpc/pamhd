/*
Tests vector field double curl calculation of PAMHD in 3d.

Copyright 2015, 2016, 2017 Ilja Honkonen
Copyright 2018 Finnish Meteorological Institute
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


#include "algorithm"
#include "array"
#include "cstdlib"
#include "iostream"
#include "limits"
#include "vector"

#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#include "gensimcell.hpp"
#include "prettyprint.hpp"

#include "divergence/remove.hpp"

using namespace std;

std::array<double, 3> f(const std::array<double, 3>& r)
{
	using std::sin;

	const double x = r[0], y = r[1], z = r[2];

	return {
		sin(x) + sin(y) + sin(z),
		sin(x) + sin(y) + sin(z),
		sin(x) + sin(y) + sin(z)
	};
}

std::array<double, 3> result(const std::array<double, 3>& r)
{
	const double x = r[0], y = r[1], z = r[2];

	return {
		sin(z) + cos(y),
		sin(z) + sin(x),
		cos(y) + sin(x)
	};
}


struct Vector {
	using data_type = std::array<double, 3>;
};

struct Result {
	using data_type = std::array<double, 3>;
};

struct Type {
	using data_type = int;
};

using Cell = gensimcell::Cell<
	gensimcell::Always_Transfer,
	Vector,
	Result,
	Type
>;


/*!
Returns maximum norm.
*/
template<class Grid> double get_diff_max_norm(const Grid& grid) {
	double local_norm = 0, global_norm = 0;
	uint64_t nr_cells_local = 0, nr_cells_global = 0;
	for (const auto& cell: grid.local_cells) {
		if ((*cell.data)[Type()] != 1) {
			continue;
		}
		nr_cells_local++;

		const auto res = result(grid.geometry.get_center(cell.id));

		local_norm =
			std::max(
			std::max(
			std::max(
				local_norm,
				std::fabs((*cell.data)[Result()][0] - res[0])),
				std::fabs((*cell.data)[Result()][1] - res[1])),
				std::fabs((*cell.data)[Result()][2] - res[2])
			);
	}

	MPI_Comm comm = grid.get_communicator();
	MPI_Allreduce(&local_norm, &global_norm, 1, MPI_DOUBLE, MPI_MAX, comm);
	MPI_Allreduce(&nr_cells_local, &nr_cells_global, 1, MPI_UINT64_T, MPI_SUM, comm);
	MPI_Comm_free(&comm);

	return global_norm / nr_cells_global;
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
	const int max_refinement_level = 0;

	double old_norm = std::numeric_limits<double>::max();
	size_t old_nr_of_cells = 0;
	for (size_t nr_of_cells = 4; nr_of_cells <= 64; nr_of_cells *= 2) {

		dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry> grid;

		const std::array<uint64_t, 3> grid_size{{nr_of_cells, nr_of_cells, nr_of_cells}};

		grid
			.set_load_balancing_method("RANDOM")
			.set_initial_length(grid_size)
			.set_neighborhood_length(neighborhood_size)
			.set_maximum_refinement_level(max_refinement_level)
			.set_periodic(true, true, true)
			.initialize(comm)
			.balance_load();

		const std::array<double, 3>
			cell_length{2*M_PI/grid_size[0], 2*M_PI/grid_size[1], 2*M_PI/grid_size[2]},
			grid_start{-M_PI, -M_PI, -M_PI};

		dccrg::Cartesian_Geometry::Parameters geom_params;
		geom_params.start = grid_start;
		geom_params.level_0_cell_length = cell_length;
		grid.set_geometry(geom_params);

		for (const auto& cell: grid.local_cells) {
			(*cell.data)[Vector()] = f(grid.geometry.get_center(cell.id));
			(*cell.data)[Result()] = {0, 0, 0};
		}
		grid.update_copies_of_remote_neighbors();

		pamhd::divergence::get_curl_curl(
			grid.local_cells,
			grid,
			[](Cell& cell_data) -> Vector::data_type& {
				return cell_data[Vector()];
			},
			[](Cell& cell_data) -> Result::data_type& {
				return cell_data[Result()];
			},
			[](Cell& cell_data) -> Type::data_type& {
				return cell_data[Type()];
			}
		);

		const double norm = get_diff_max_norm(grid);

		if (old_nr_of_cells > 0) {
			if (norm > old_norm) {
				if (grid.get_rank() == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< ": X norm with " << nr_of_cells
						<< " cells " << norm
						<< " is larger than with " << old_nr_of_cells
						<< " cells " << old_norm
						<< std::endl;
				}
				abort();
			}
		}

		if (old_nr_of_cells > 0) {
			const double
				order_of_accuracy
					= -log(norm / old_norm)
					/ log(double(nr_of_cells) / old_nr_of_cells);
			if (order_of_accuracy < 2.3) {
				if (grid.get_rank() == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< ": Order of accuracy from " << old_nr_of_cells
						<< " cells to " << nr_of_cells
						<< " cells is too low: " << order_of_accuracy
						<< std::endl;
				}
				abort();
			}
		}

		old_nr_of_cells = nr_of_cells;
		old_norm = norm;
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
