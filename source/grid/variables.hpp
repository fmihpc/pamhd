/*
TODO

Copyright 2022 Finnish Meteorological Institute
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

#ifndef PAMHD_GRID_VARIABLES_HPP
#define PAMHD_GRID_VARIABLES_HPP


namespace pamhd {

//! Calculate whether cell and neighbor share a face whenever grid changes
struct Is_Face_Neighbor {
	bool is_face_neighbor = false;
	template<
		class Grid, class Cell_Item, class Neighbor_Item
	> void update(
		const Grid& grid,
		const Cell_Item& cell,
		const Neighbor_Item& neighbor,
		const int&,
		const Is_Face_Neighbor&
	) {
		const auto
			temp_cli = grid.mapping.get_cell_length_in_indices(cell.id),
			temp_nli = grid.mapping.get_cell_length_in_indices(neighbor.id);
		if (temp_cli > std::numeric_limits<int>::max()) {
			throw std::runtime_error(__FILE__"(" + std::to_string(__LINE__) + ") Cell length in indices too large for int");
		}
		if (temp_nli > std::numeric_limits<int>::max()) {
			throw std::runtime_error(__FILE__"(" + std::to_string(__LINE__) + ") Neighbor length in indices too large for int");
		}
		const int cell_length_i = temp_cli, neigh_length_i = temp_nli;
		if (
			(neighbor.x == -neigh_length_i
			and (neighbor.y < cell_length_i)
			and neighbor.y > -neigh_length_i
			and (neighbor.z < cell_length_i)
			and neighbor.z > -neigh_length_i)

			or

			(neighbor.y == -neigh_length_i
			and (neighbor.x < cell_length_i)
			and neighbor.x > -neigh_length_i
			and (neighbor.z < cell_length_i)
			and neighbor.z > -neigh_length_i)

			or

			(neighbor.z == -neigh_length_i
			and (neighbor.x < cell_length_i)
			and neighbor.x > -neigh_length_i
			and (neighbor.y < cell_length_i)
			and neighbor.y > -neigh_length_i)
		) {
			is_face_neighbor = true;
		}
	}
};

} // namespace

#endif // ifndef PAMHD_GRID_VARIABLES_HPP
