/*
Functions for working with divergence of vector field.

Copyright 2014, 2015, 2016, 2017 Ilja Honkonen
Copyright 2018, 2019, 2022 Finnish Meteorological Institute
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

#ifndef PAMHD_MATH_STAGGERED_HPP
#define PAMHD_MATH_STAGGERED_HPP


#include "string"

#include "mpi.h"


namespace pamhd {
namespace math {


/*!
Same as get_divergence() but for staggered vector variable...

TODO
*/
template <
	class Cell_Iterator,
	class Grid,
	class Vector_Getter,
	class Divergence_Getter,
	class Cell_Type_Getter
> double get_divergence_staggered(
	const Cell_Iterator& cells,
	Grid& grid,
	Vector_Getter Vector,
	Divergence_Getter Divergence,
	Cell_Type_Getter Cell_Type
) {
	double local_divergence = 0, global_divergence = 0;
	uint64_t local_calculated_cells = 0, global_calculated_cells = 0;
	for (const auto& cell: cells) {
		if (Cell_Type(*cell.data) < 0) {
			continue;
		}
		local_calculated_cells++;

		const int cell_length_i = grid.mapping.get_cell_length_in_indices(cell.id);
		const auto cell_length_g = grid.geometry.get_length(cell.id);
		auto& div = Divergence(*cell.data);
		div = 0.0;

		for (const auto& neighbor: cell.neighbors_of) {
			if (Cell_Type(*neighbor.data) < 0) {
				continue;
			}

			const int neigh_length_i = grid.mapping.get_cell_length_in_indices(neighbor.id);
			if (neigh_length_i > 2*cell_length_i or cell_length_i > 2*neigh_length_i) {
				throw std::runtime_error(
					__FILE__"(" + std::to_string(__LINE__)
					+ ") Logical size difference too large between cells "
					+ std::to_string(cell.id) + " and " + std::to_string(neighbor.id)
					+ ": " + std::to_string(cell_length_i) + " vs "
					+ std::to_string(neigh_length_i)
				);
			}

			if (
				neighbor.x == -neigh_length_i
				and (neighbor.y < cell_length_i)
				and neighbor.y > -neigh_length_i
				and (neighbor.z < cell_length_i)
				and neighbor.z > -neigh_length_i
			) {
				const auto diff = (Vector(*cell.data)[0] - Vector(*neighbor.data)[0]) / cell_length_g[0];
				if (neigh_length_i < cell_length_i) {
					div += diff / 4;
				} else {
					div += diff;
				}
			}
			if (
				neighbor.y == -neigh_length_i
				and (neighbor.x < cell_length_i)
				and neighbor.x > -neigh_length_i
				and (neighbor.z < cell_length_i)
				and neighbor.z > -neigh_length_i
			) {
				const auto diff = (Vector(*cell.data)[1] - Vector(*neighbor.data)[1]) / cell_length_g[1];
				if (neigh_length_i < cell_length_i) {
					div += diff / 4;
				} else {
					div += diff;
				}
			}
			if (
				neighbor.z == -neigh_length_i
				and (neighbor.x < cell_length_i)
				and neighbor.x > -neigh_length_i
				and (neighbor.y < cell_length_i)
				and neighbor.y > -neigh_length_i
			) {
				const auto diff = (Vector(*cell.data)[2] - Vector(*neighbor.data)[2]) / cell_length_g[2];
				if (neigh_length_i < cell_length_i) {
					div += diff / 4;
				} else {
					div += diff;
				}
			}
		}
		local_divergence += std::fabs(div);
	}

	MPI_Comm comm = grid.get_communicator();
	MPI_Allreduce(
		&local_divergence,
		&global_divergence,
		1,
		MPI_DOUBLE,
		MPI_SUM,
		comm
	);
	MPI_Allreduce(
		&local_calculated_cells,
		&global_calculated_cells,
		1,
		MPI_UINT64_T,
		MPI_SUM,
		comm
	);
	MPI_Comm_free(&comm);

	return global_divergence / global_calculated_cells;
}

}} // namespaces

#endif // ifndef PAMHD_MATH_STAGGERED_HPP
