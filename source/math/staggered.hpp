/*
Functions for working with divergence of vector field.

Copyright 2014, 2015, 2016, 2017 Ilja Honkonen
Copyright 2018, 2019, 2022, 2023 Finnish Meteorological Institute
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
Same as get_divergence() but for vector variable that's stored on cell faces.

Vector_Pos returns vector whose components are at cell faces on positive
side from cell center and are normal to said faces, Vector_Neg returns
components on negative side faces pointing in positive direction.
*/
template <
	class Cell_Iterator,
	class Grid,
	class Vector_Pos_Getter,
	class Vector_Neg_Getter,
	class Divergence_Getter,
	class Is_Primary_Face_Getter,
	class Cell_Type_Getter
> double get_divergence_staggered(
	const Cell_Iterator& cells,
	Grid& grid,
	const Vector_Pos_Getter Vector_Pos,
	const Vector_Neg_Getter Vector_Neg,
	const Divergence_Getter Divergence,
	const Is_Primary_Face_Getter PFace,
	const Cell_Type_Getter Cell_Type
) {
	double local_divergence = 0, global_divergence = 0;
	uint64_t local_calculated_cells = 0, global_calculated_cells = 0;
	for (const auto& cell: cells) {
		if (Cell_Type(*cell.data) <= 0) {
			continue;
		}
		local_calculated_cells++;

		const auto cell_length = grid.geometry.get_length(cell.id);
		auto& div = Divergence(*cell.data);
		div = 0.0;

		const auto primary = PFace(*cell.data);
		for (size_t i = 0; i < 3; i++) {
			if (primary[2*i + 0]) div -= Vector_Neg(*cell.data)[i] / cell_length[i];
			if (primary[2*i + 1]) div += Vector_Pos(*cell.data)[i] / cell_length[i];
		}

		for (const auto& neighbor: cell.neighbors_of) {
			const auto n = neighbor.face_neighbor;
			if (n == 0) continue;
			if (n == -1 and primary[0]) continue;
			if (n == +1 and primary[1]) continue;
			if (n == -2 and primary[2]) continue;
			if (n == +2 and primary[3]) continue;
			if (n == -3 and primary[4]) continue;
			if (n == +3 and primary[5]) continue;

			if (Cell_Type(*neighbor.data) < 0) {
				continue;
			}

			const size_t dim = std::abs(n) - 1;
			const double factor = [&]{
				if (neighbor.relative_size > 0) {
					return 0.25 / cell_length[dim];
				} else {
					return 1.0 / cell_length[dim];
				}
			}();
			if (neighbor.face_neighbor < 0) {
				div -= Vector_Pos(*neighbor.data)[dim] * factor;
			} else {
				div += Vector_Neg(*neighbor.data)[dim] * factor;
			}
		}
		local_divergence += std::abs(div);
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
