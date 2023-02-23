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


#include "string"
#include "type_traits"


namespace pamhd {
namespace grid {


/*! Variable indicating refinement level towards which grid cells should be refined.
*/
struct Target_Refinement_Level {
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
		Faces(*cell.data) = {true, true, true, true, true, true};
		for (const auto& neighbor: cell.neighbors_of) {
			const auto n = neighbor.face_neighbor;
			if (n == 0) {
				continue;
			}
			if (neighbor.relative_size < 0) {
				continue; // larger neighbor never has priority
			}
			if (neighbor.relative_size > 0) {
				// smaller neighbor always has priority
				if (n == -1) Faces(*cell.data)[0] = false;
				if (n == +1) Faces(*cell.data)[1] = false;
				if (n == -2) Faces(*cell.data)[2] = false;
				if (n == +2) Faces(*cell.data)[3] = false;
				if (n == -3) Faces(*cell.data)[4] = false;
				if (n == +3) Faces(*cell.data)[5] = false;
				continue;
			}
			if (n == -1) Faces(*cell.data)[0] = false;
			if (n == -2) Faces(*cell.data)[2] = false;
			if (n == -3) Faces(*cell.data)[4] = false;
		}
	}
}


/*! Common accessors for values stored on cell edges.
*/
template<class Data_Type> struct Edge_Type {
	std::array<Data_Type, 12> edge;

	/*
	par_dim_i = dimension to which electric field is parallel to,
	0 = x, 1 = y, 2 = z.
	first_perp_dim_i = negative or positive side of cell in
	lexically earlier dimension perpendicular to electric field,
	e.g. if par_dim_i = 0, first_perp_dim_i = 1 means positive side of
	cell in y dimension.
	second_perp_dim = neg or pos side in lexically later dimension,
	if par_dim_i = 2, second_perp_dim_i = 0 means negative side of
	cell in y dimension.
	par_dim_i | first_perp_dim_i | second..._i | E on edge of cell
	    0     |         0        |      0      | x directed: -y, -z
	    0     |         0        |      1      | x dir:      -y, +z
	    0     |         1        |      1      | x dir:      +y, +z
	...
	    2     |         1        |      0      | z dir:      +x, -y
	    2     |         1        |      1      | z dir:      +x, +y
	*/
	const Data_Type& operator()(
		const size_t par_dim_i,
		const size_t first_perp_dim_i,
		const size_t second_perp_dim_i
	) const {
		if (par_dim_i > 2) {
			throw std::domain_error("Parallel dimension > 2");
		}
		if (first_perp_dim_i > 1) {
			throw std::domain_error("First perpendicular dimension > 1");
		}
		if (second_perp_dim_i > 1) {
			throw std::domain_error("Second perpendicular dimension > 1");
		}
		return this->edge[par_dim_i*2*2 + first_perp_dim_i*2 + second_perp_dim_i];
	}

	// https://stackoverflow.com/a/123995
	Data_Type& operator()(
		const size_t par_dim_i,
		const size_t first_perp_dim_i,
		const size_t second_perp_dim_i
	) {
		return const_cast<Data_Type&>(static_cast<const Edge_Type<Data_Type>&>(*this).operator()(par_dim_i, first_perp_dim_i, second_perp_dim_i));
	}

	Edge_Type<Data_Type>& operator=(const decltype(edge)& given) {
		this->edge = given;
		return *this;
	}

	#ifdef MPI_VERSION
	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype() const {
		if constexpr (std::is_same_v<Data_Type, double>) {
			return std::make_tuple((void*) this->edge.data(), this->edge.size(), MPI_DOUBLE);
		} else {
			return std::make_tuple(nullptr, 0, MPI_BYTE);
		}
	}
	#endif
};

/*! Records which of cell's edges are primary.

\see Edge_Type for further info.
*/
struct Is_Primary_Edge {
	using data_type = Edge_Type<bool>;
};


/*! Returns info on potentially shared edge with given neighbor.

Returns (d, a, b) where d is dimension (0..2) parallel to edge,
a and b are 1st and 2nd perpendicular dimension and are 0 if
neighbor is on negative side of cell or 1 if on positive. If d < 0
neighbor doesn't share or shares more than one edge with cell.
*/
std::array<int, 3> edge_neighbor(
	const size_t cell_length,
	const size_t neigh_length,
	const int neigh_offset_x,
	const int neigh_offset_y,
	const int neigh_offset_z
) {
	using std::numeric_limits;
	using std::runtime_error;
	using std::to_string;

	if (neigh_offset_x == 0 and neigh_offset_y == 0 and neigh_offset_z == 0) {
		throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
	}

	if (cell_length > numeric_limits<int>::max()) {
		throw runtime_error(__FILE__"(" + to_string(__LINE__) + ") Cell length in indices too large for int");
	}
	if (neigh_length > numeric_limits<int>::max()) {
		throw runtime_error(__FILE__"(" + to_string(__LINE__) + ") Neighbor length in indices too large for int");
	}
	const int cell_len = cell_length, neigh_len = neigh_length;

	std::array<int, 3> ret_val{-1, -1, -1};
	// get dimension of return value
	if (neigh_offset_x < cell_len or neigh_offset_x > -neigh_len) {
		ret_val[0] = 0;
	}
	if (neigh_offset_y < cell_len or neigh_offset_y > -neigh_len) {
		if (ret_val[0] > -1) {
			// face neighbor
			ret_val[0] = -1;
			return ret_val;
		} else {
			ret_val[0] = 1;
		}
	}
	if (neigh_offset_z < cell_len or neigh_offset_z > -neigh_len) {
		if (ret_val[0] > -1) {
			// face neighbor
			ret_val[0] = -1;
			return ret_val;
		} else {
			ret_val[0] = 2;
		}
	}
	if (ret_val[0] < 0) {
		// not edge neighbor
		return ret_val;
	}

	if (neigh_offset_x == -neigh_len) {
		// 1st perp dim in any case
		ret_val[1] = 0;
	}
	if (neigh_offset_x == cell_len) {
		ret_val[1] = 1;
	}

	if (neigh_offset_y == -neigh_len) {
		if      (ret_val[0] == 0) ret_val[1] = 0;
		else if (ret_val[0] == 2) ret_val[2] = 0;
		else throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
	}
	if (neigh_offset_y == cell_len) {
		if      (ret_val[0] == 0) ret_val[1] = 1;
		else if (ret_val[0] == 2) ret_val[2] = 1;
		else throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
	}

	if (neigh_offset_z == -neigh_len) {
		// 2nd perp dim in any case
		ret_val[2] = 0;
	}
	if (neigh_offset_z == cell_len) {
		ret_val[2] = 1;
	}

	return ret_val;
}

/*! Updates Is_Primary_Edge variable of given cells.
*/
template <
	class Cells,
	class Grid,
	class Edges_Getter
> void update_primary_edges(
	const Cells& cells,
	const Grid& grid,
	const Edges_Getter& Edges
) {
	for (const auto& cell: cells) {
		const auto clen = grid.mapping.get_cell_length_in_indices(cell.id);

		Edges(*cell.data) = {
			true, true, true, true, // x-directed edges: -y,-z;-y,+z,...
			true, true, true, true, // y-directed edges: -x,-z;-x,+z,...
			true, true, true, true};// ...

		for (const auto& neighbor: cell.neighbors_of) {
			const auto en = edge_neighbor(
				clen,
				grid.mapping.get_cell_length_in_indices(neighbor.id),
				neighbor.x,
				neighbor.y,
				neighbor.z);
			if (neighbor.face_neighbor != 0 and en[0] >= 0) {
				throw std::runtime_error(__FILE__"(" + std::to_string(__LINE__) + ")");
			}
			if (en[0] < 0) {
				continue;
			}
			if (neighbor.relative_size < 0) {
				continue; // larger neighbor never has priority
			}
			if (neighbor.relative_size > 0) {
				// smaller neighbor always has priority
				Edges(*cell.data)(en[0], en[1], en[2]) = false;
				continue;
			}
			if (en[1] != 1 and en[2] != 1) {
				// identical neighbor usually has priority
				Edges(*cell.data)(en[0], en[1], en[2]) = false;
				continue;
			}
		}
	}
}


}} // namespaces

#endif // ifndef PAMHD_GRID_AMR_HPP
