/*
Grid stuff related to adaptive mesh refinement of PAMHD.

Copyright 2023 Finnish Meteorological Institute
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

* Neither the name of copyright holders nor the names of their contributors
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

#ifndef PAMHD_GRID_AMR_HPP
#define PAMHD_GRID_AMR_HPP

namespace pamhd {
namespace grid {


/*! Variable indicating refinement level towards which grid cells should be refined.
*/
struct Targer_Refinement_Level {
	using data_type = int;
	static const std::string get_name() { return {"target refinement level"}; }
	static const std::string get_option_name() { return {"target-ref-lvl"}; }
	static const std::string get_option_help() { return {"Target refinement level of cell"}; }
};


/*! Records which of cell's faces are primary.

By default positive faces are primary, i.e. values in
corresponding negative face of neighbor(s) in positive
direction from cell are ignored.

Faces of smaller neighbors are always primary.

0 == -x, 1 == +x, ..., 5 == +z
*/
struct Is_Primary_Face {
	using data_type = std::array<bool, 6>;
	#ifdef MPI_VERSION
	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype() const {
		return std::make_tuple(nullptr, 0, MPI_BYTE);
	}
	#endif
};


/*! Updates Is_Primary_Face variable of given cells.
*/
template <
	class Cells,
	class Faces_Getter
> void update_primary_faces(
	const Cells& cells,
	const Faces_Getter& Faces
) {
	for (const auto& cell: cells) {
		Faces(*cell.data) = {false, true, false, true, false, true};
		for (const auto& neighbor: cell.neighbors_of) {
			const auto n = neighbor.face_neighbor;
			if (n == 0) {
				continue;
			}
			if (neighbor.relative_size > 0) {
				if (n == +1) Faces(*cell.data)[1] = false;
				if (n == +2) Faces(*cell.data)[3] = false;
				if (n == +3) Faces(*cell.data)[5] = false;
				// default values for n < 0
			}
			if (neighbor.relative_size < 0) {
				if (n == -1) Faces(*cell.data)[0] = true;
				if (n == -2) Faces(*cell.data)[2] = true;
				if (n == -3) Faces(*cell.data)[4] = true;
				// default values for n > 0
			}
			// default values for relative_size == 0
		}
	}
}


}} // namespaces

#endif // ifndef PAMHD_GRID_AMR_HPP
