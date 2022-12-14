/*
Grid-related stuff common to all PAMHD models.

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
namespace grid {

/*! Calculate whether cell and neighbor share a face whenever grid changes.

Non-face neighbors have face_neighbor == 0,
One or more face neighbors on negative side of cell in y dimension
have face_neighbor == -2, positive side of cell in x == +1, etc.
*/
struct Is_Face_Neighbor {
	int face_neighbor = 0;
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
			(neighbor.y < cell_length_i)
			and neighbor.y > -neigh_length_i
			and (neighbor.z < cell_length_i)
			and neighbor.z > -neigh_length_i
		) {
			if (neighbor.x == -neigh_length_i) {
				face_neighbor = -1;
			} else {
				face_neighbor = 1;
			}
		} else if (
			(neighbor.x < cell_length_i)
			and neighbor.x > -neigh_length_i
			and (neighbor.z < cell_length_i)
			and neighbor.z > -neigh_length_i
		) {
			if (neighbor.y == -neigh_length_i) {
				face_neighbor = -2;
			} else {
				face_neighbor = 2;
			}
		} else if (
			(neighbor.x < cell_length_i)
			and neighbor.x > -neigh_length_i
			and (neighbor.y < cell_length_i)
			and neighbor.y > -neigh_length_i
		) {
			if (neighbor.z == -neigh_length_i) {
				face_neighbor = -3;
			} else {
				face_neighbor = 3;
			}
		} else {
			face_neighbor = 0;
		}
	}
};

//! Whether logical size of neighbor is smaller than cell.
struct Is_Smaller {
	bool is_smaller = false;
	template<
		class Grid, class Cell_Item, class Neighbor_Item
	> void update(
		const Grid& grid,
		const Cell_Item& cell,
		const Neighbor_Item& neighbor,
		const int&,
		const Is_Smaller&
	) {
		const auto
			cell_len = grid.mapping.get_cell_length_in_indices(cell.id),
			neigh_len = grid.mapping.get_cell_length_in_indices(neighbor.id);
		if (neigh_len < cell_len) {
			is_smaller = true;
		} else {
			is_smaller = false;
		}
	}
};

}} // namespaces

#endif // ifndef PAMHD_GRID_VARIABLES_HPP
