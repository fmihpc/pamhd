/*
Staggered version of get1d_div_amr.cpp.

Copyright 2014, 2015, 2016, 2017 Ilja Honkonen
Copyright 2018, 2022 Finnish Meteorological Institute
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
#include "set"
#include "vector"

#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#include "gensimcell.hpp"

#include "grid/variables.hpp"
#include "math/staggered.hpp"


double function(const double x)
{
	return std::pow(std::sin(x / 2), 2);
}

double div_of_function(const double x)
{
	return std::sin(x / 2) * std::cos(x / 2);
}


//! each component is normal to and at corresponding + dir face
struct Vector {
	using data_type = std::array<double, 3>;
};

struct Divergence {
	using data_type = double;
};

struct Type {
	using data_type = int;
};

using Cell = gensimcell::Cell<
	gensimcell::Always_Transfer,
	Vector,
	Divergence,
	Type
>;


/*!
Returns maximum norm if p == 0
*/
template<class Grid> double get_max_norm(
	const Grid& grid,
	const size_t dimension
) {
	double local_norm = 0, global_norm = 0;
	for (const auto& cell: grid.local_cells()) {
		if ((*cell.data)[Type()] != 1) {
			continue;
		}

		const auto center = grid.geometry.get_center(cell.id);

		local_norm = std::max(
			local_norm,
			std::fabs((*cell.data)[Divergence()] - div_of_function(center[dimension]))
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

	std::array<double, 3> old_norms{{
		std::numeric_limits<double>::max(),
		std::numeric_limits<double>::max(),
		std::numeric_limits<double>::max()
	}};

	size_t old_nr_of_cells = 0;
	for (size_t nr_of_cells = 8; nr_of_cells <= 512; nr_of_cells *= 2) {

		std::array<dccrg::Dccrg<
			Cell,
			dccrg::Cartesian_Geometry,
			std::tuple<>,
			std::tuple<pamhd::Is_Face_Neighbor>
		>, 3> grids;

		const std::array<std::array<uint64_t, 3>, 3>
			grid_sizes{{
				{{nr_of_cells + 2, 1, 1}},
				{{1, nr_of_cells + 2, 1}},
				{{1, 1, nr_of_cells + 2}}
			}};

		const int max_ref_lvl = 2;
		for (size_t dim = 0; dim < 3; dim++) {
			grids[dim]
				.set_load_balancing_method("RANDOM")
				.set_initial_length(grid_sizes[dim])
				.set_neighborhood_length(0)
				.set_maximum_refinement_level(max_ref_lvl)
				.initialize(comm);
		}

		const std::array<std::array<double, 3>, 3>
			cell_lengths{{
				{{2 * M_PI / (grid_sizes[0][0] - 2), 1, 1}},
				{{1, 2 * M_PI / (grid_sizes[1][1] - 2), 1}},
				{{1, 1, 2 * M_PI / (grid_sizes[2][2] - 2)}}
			}},
			grid_starts{{
				{{-cell_lengths[0][0], 0, 0}},
				{{0, -cell_lengths[1][1], 0}},
				{{0, 0, -cell_lengths[2][2]}}
			}};

		std::array<dccrg::Cartesian_Geometry::Parameters, 3> geom_params;
		for (size_t dim = 0; dim < 3; dim++) {
			geom_params[dim].start = grid_starts[dim];
			geom_params[dim].level_0_cell_length = cell_lengths[dim];
			grids[dim].set_geometry(geom_params[dim]);
		}

		for (int i = 0; i < max_ref_lvl; i++) {
			for (size_t dim = 0; dim < 3; dim++) {
				for (const auto& cell: grids[dim].local_cells()) {
					if (
						(grids[dim].geometry.get_center(cell.id)[dim] > 2 * M_PI / 8
						and grids[dim].geometry.get_center(cell.id)[dim] < 6 * M_PI / 8)
						or (grids[dim].geometry.get_center(cell.id)[dim] > 10 * M_PI / 8
						and grids[dim].geometry.get_center(cell.id)[dim] < 14 * M_PI / 8)
					) {
						grids[dim].refine_completely(cell.id);
					}
				}
				grids[dim].stop_refining();
			}
		}

		// real number of cells
		uint64_t local_real_nr_cells = 0, real_nr_cells = 0;
		std::set<uint64_t> uniq_indices;
		for (const auto& cell: grids[0].local_cells()) {
			const auto index = grids[0].mapping.get_indices(cell.id)[0];
			uniq_indices.insert(index);
		}
		local_real_nr_cells = uniq_indices.size();
		MPI_Comm comm = grids[0].get_communicator();
		MPI_Allreduce(&local_real_nr_cells, &real_nr_cells, 1, MPI_UINT64_T, MPI_SUM, comm);
		MPI_Comm_free(&comm);

		for (size_t dim = 0; dim < 3; dim++) {
			for (const auto& cell: grids[dim].local_cells()) {
				auto& vec = (*cell.data)[Vector()];
				vec[0] = vec[1] = vec[2] = 0;

				const auto
					center = grids[dim].geometry.get_center(cell.id)[dim],
					length = grids[dim].geometry.get_length(cell.id)[dim];
				vec[dim] = function(center + length/2);
			}
			grids[dim].update_copies_of_remote_neighbors();
		}

		// exclude one layer of boundary cells
		for (size_t dim = 0; dim < 3; dim++) {
			for (const auto& cell: grids[dim].local_cells()) {
				const auto center = grids[dim].geometry.get_center(cell.id)[dim];
				if (center > 0 and center < 2 * M_PI) {
					(*cell.data)[Type()] = 1;
				} else {
					(*cell.data)[Type()] = 0;
				}
			}
		}

		auto Vector_Getter = [](Cell& cell_data) -> Vector::data_type& {
			return cell_data[Vector()];
		};
		auto Divergence_Getter = [](Cell& cell_data) -> Divergence::data_type& {
			return cell_data[Divergence()];
		};
		auto Type_Getter = [](Cell& cell_data) -> Type::data_type& {
			return cell_data[Type()];
		};

		for (size_t dim = 0; dim < 3; dim++) {
			pamhd::math::get_divergence_staggered(
				grids[dim].local_cells(),
				grids[dim],
				Vector_Getter,
				Divergence_Getter,
				Type_Getter
			);
		}

		const std::array<double, 3> norms{{
			get_max_norm(grids[0], 0),
			get_max_norm(grids[1], 1),
			get_max_norm(grids[2], 2)
		}};

		for (size_t dim = 0; dim < 3; dim++) {
			if (norms[dim] > old_norms[dim]) {
				if (grids[dim].get_rank() == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " in dim " << dim
						<< ": max norm with " << real_nr_cells
						<< " cells " << norms[dim]
						<< " is larger than with " << old_nr_of_cells
						<< " cells " << old_norms[dim]
						<< std::endl;
				}
				abort();
			}
		}

		for (size_t dim = 0; dim < 3; dim++) {
			if (old_nr_of_cells > 0) {
				const double order_of_accuracy
					= -log(norms[dim] / old_norms[dim])
					/ log(double(real_nr_cells) / old_nr_of_cells);

				if (real_nr_cells > 300 and order_of_accuracy < 1.9) {
					if (grids[dim].get_rank() == 0) {
						std::cerr << __FILE__ << ":" << __LINE__
							<< " in dim " << dim
							<< ": order of accuracy from "
							<< old_nr_of_cells << " to " << real_nr_cells
							<< " is too low: " << order_of_accuracy
							<< std::endl;
					}
					abort();
				}
			}
		}

		old_nr_of_cells = real_nr_cells;
		old_norms = norms;
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
