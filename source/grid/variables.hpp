/*
Grid-related stuff common to all PAMHD models.

Copyright 2022, 2023 Finnish Meteorological Institute
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


#include "limits"
#include "stdexcept"
#include "string"


namespace pamhd {
namespace grid {

/*! Stores whether cell and neighbor share a face.

Non-face neighbors have face_neighbor == 0.
One or more face neighbors on negative side of cell in y dimension
have face_neighbor == -2, positive side of cell in x == +1, etc.
*/
struct Face_Neighbor {
	int face_neighbor = 0;
	template<
		class Grid, class Cell_Item, class Neighbor_Item
	> void update(
		const Grid& grid,
		const Cell_Item& cell,
		const Neighbor_Item& neighbor,
		const int&,
		const Face_Neighbor&
	) {
		using std::numeric_limits;
		using std::range_error;
		using std::to_string;

		const auto
			temp_cli = grid.mapping.get_cell_length_in_indices(cell.id),
			temp_nli = grid.mapping.get_cell_length_in_indices(neighbor.id);
		if (temp_cli > numeric_limits<int>::max()) {
			throw range_error(__FILE__"(" + to_string(__LINE__) + ") Cell length in indices too large for int");
		}
		if (temp_nli > numeric_limits<int>::max()) {
			throw range_error(__FILE__"(" + to_string(__LINE__) + ") Neighbor length in indices too large for int");
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
			} else if (neighbor.x == cell_length_i) {
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
			} else if (neighbor.y == cell_length_i) {
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
			} else if (neighbor.z == cell_length_i) {
				face_neighbor = 3;
			}
		} else {
			face_neighbor = 0;
		}
	}
};

/*! Logical size of neighbor compared to cell.

relative_size == 0 if neighbor's refinement level i.e. size is same,
... == +1 if neighbor's ref lvl is larger by one i.e.
neighbor is 1/2x of cell in each dimension,
... == -1 if neighbor's ref lvl is smaller by one i.e.
neighbor is 2x larger than cell in each dimension, etc.
*/
struct Relative_Size {
	int relative_size = 0;
	template<
		class Grid, class Cell_Item, class Neighbor_Item
	> void update(
		const Grid& grid,
		const Cell_Item& cell,
		const Neighbor_Item& neighbor,
		const int&,
		const Relative_Size&
	) {
		this->relative_size
			= grid.mapping.get_refinement_level(neighbor.id)
			- grid.mapping.get_refinement_level(cell.id);
	}
};

/*! Stores whether cell and neighbor share an edge.

edge_neighbor[0] specifies dimension in which cell and neighbor overlap,
with 0 == x, ..., 2 == z. If neighbor doesn't share an edge or shares
more than one edge with cell, edge_neighbor[0] == -1.

edge_neighbor[1] specifies offset of neighbor relative to cell in
first dimension perpendicular to edge_neighbor[0], for example
if edge_neighbor[0] == +1, edge_neighbor[1] specifies offset in x dim.
Value of -1 means neighbor in on negative side from cell in that dimension,
value of +1 means neighbor is on positive side.

edge_neighbor[2] specifies offset of neighbor in second dim
perpendicular to edge_neighbor[0], for example if edge_neighbor[0] == 2,
edge_neighbor[2] specifies offset in y dimension.
*/
struct Edge_Neighbor {
	std::array<int, 3> edge_neighbor{-1, 0, 0};
	template<
		class Grid, class Cell_Item, class Neighbor_Item
	> void update(
		const Grid& grid,
		const Cell_Item& cell,
		const Neighbor_Item& neighbor,
		const int&,
		const Edge_Neighbor&
	) {
		using std::numeric_limits;
		using std::range_error;
		using std::to_string;

		const auto
			cell_length = grid.mapping.get_cell_length_in_indices(cell.id),
			neigh_length = grid.mapping.get_cell_length_in_indices(neighbor.id);
		if (cell_length > numeric_limits<int>::max()) {
			throw range_error(__FILE__"(" + to_string(__LINE__) + ") Cell length in indices too large for int");
		}
		if (neigh_length > numeric_limits<int>::max()) {
			throw range_error(__FILE__"(" + to_string(__LINE__) + ") Neighbor length in indices too large for int");
		}
		const int cell_len = cell_length, neigh_len = neigh_length;

		// get dimension of return value
		if ((neighbor.x < cell_len) and neighbor.x > -neigh_len) {
			this->edge_neighbor[0] = 0;
		}
		if ((neighbor.y < cell_len) and neighbor.y > -neigh_len) {
			if (this->edge_neighbor[0] > -1) {
				// face neighbor
				this->edge_neighbor[0] = -1;
				return;
			} else {
				this->edge_neighbor[0] = 1;
			}
		}
		if ((neighbor.z < cell_len) and neighbor.z > -neigh_len) {
			if (this->edge_neighbor[0] > -1) {
				// face neighbor
				this->edge_neighbor[0] = -1;
				return;
			} else {
				this->edge_neighbor[0] = 2;
			}
		}
		if (this->edge_neighbor[0] < 0) {
			// not edge neighbor
			return;
		}

		if (neighbor.x == -neigh_len) {
			// 1st perp dim in any case
			this->edge_neighbor[1] = -1;
		}
		if (neighbor.x == cell_len) {
			this->edge_neighbor[1] = +1;
		}

		if (neighbor.y == -neigh_len) {
			switch (this->edge_neighbor[0]) {
			case 0:
				this->edge_neighbor[1] = -1;
				break;
			case 2:
				this->edge_neighbor[2] = -1;
				break;
			default:
				break;
			}
		}
		if (neighbor.y == cell_len) {
			switch (this->edge_neighbor[0]) {
			case 0:
				this->edge_neighbor[1] = +1;
				break;
			case 2:
				this->edge_neighbor[2] = +1;
				break;
			default:
				break;
			}
		}

		if (neighbor.z == -neigh_len) {
			// 2nd perp dim in any case
			this->edge_neighbor[2] = -1;
		}
		if (neighbor.z == cell_len) {
			this->edge_neighbor[2] = +1;
		}

		if (this->edge_neighbor[1] == 0 or this->edge_neighbor[2] == 0) {
			// not edge neighbor
			this->edge_neighbor[0] = -1;
			return;
		}
	}
};

}} // namespaces

#endif // ifndef PAMHD_GRID_VARIABLES_HPP
