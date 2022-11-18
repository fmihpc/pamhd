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


#include "mpi.h"


namespace pamhd {
namespace math {


/*!
Calculates divergence of vector variable.

Divergence is calculated in all @cells for
which @Cell_Type returns 0 or 1.

Returns the total divergence i.e. sum of absolute
divergence in given cells of all processes divided
by number of cells on all processes in which divergence
was calculated.

Variable returned by @Vector must have a defined value
in face neighbors of cells where divergence is calculated.
Dimension(s) where at least one neighbor is missing
doesn't contribute to divergence.

@Cell_Type must be an object that, when given a reference
to data of one grid cell, returns 0 or 1 if divergence
should be calculated for that cell and -1 otherwise.
@Vector must be an object that, when given a reference
to the data of one grid cell, returns a *reference* to
the data from which divergence should be calculated.
Similarly @Divergence should return a *reference* to
data in which calculated divergence should be stored.

Example:

struct Cell_Data {
	std::array<double, 3> vector_data;
	double scalar_data;
	int calculate_div;
	std::tuple<...> get_mpi_datatype() {...}
};

dccrg::Dccrg<Cell_Data, ...> grid;

pamhd::math::get_divergence_staggered(
	grid.local_cells(),
	grid,
	[](Cell_Data& cell_data)->auto& {
		return cell_data.vector_data;
	},
	[](Cell_Data& cell_data)->auto& {
		return cell_data.scalar_data;
	},
	[](Cell_Data& cell_data)->auto& {
		return cell_data.calculate_div;
	}
);
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

		const auto cell_length = grid.geometry.get_length(cell.id);
		auto& div = Divergence(*cell.data);
		div = 0.0;

		for (const auto& neighbor: cell.neighbors_of) {
			if (Cell_Type(*neighbor.data) < 0) {
				continue;
			}
			if (neighbor.x == -1 and neighbor.y == 0 and neighbor.z == 0) {
				div += (Vector(*cell.data)[0] - Vector(*neighbor.data)[0]) / cell_length[0];
			}
			if (neighbor.x == 0 and neighbor.y == -1 and neighbor.z == 0) {
				div += (Vector(*cell.data)[1] - Vector(*neighbor.data)[1]) / cell_length[1];
			}
			if (neighbor.x == 0 and neighbor.y == 0 and neighbor.z == -1) {
				div += (Vector(*cell.data)[2] - Vector(*neighbor.data)[2]) / cell_length[2];
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
