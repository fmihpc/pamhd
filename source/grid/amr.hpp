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


#include "algorithm"
#include "initializer_list"
#include "stdexcept"
#include "string"
#include "type_traits"

#include "grid/options.hpp"


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


/*! Common accessors for values stored on cell faces.
*/
template<class Data_Type> struct Face_Type {
	std::array<Data_Type, 6> face;

	/*! Returns data of given cell face.
	dir: -1 == -x, +1 == +x, -2 == -y, ..., +3 == +z face
	*/
	const Data_Type& operator()(const int dir) const {
		if (dir == 0 or dir < -3 or dir > +3) {
			throw std::domain_error(__FILE__ "(" + std::to_string(__LINE__) + "): Invalid direction, must be -3,-2,-1,1,2,3 but is " + std::to_string(dir));
		}
		if        (dir == -1) {
			return this->face[0];
		} else if (dir == +1) {
			return this->face[1];
		} else if (dir == -2) {
			return this->face[2];
		} else if (dir == +2) {
			return this->face[3];
		} else if (dir == -3) {
			return this->face[4];
		} else if (dir == +3) {
			return this->face[5];
		} else throw std::runtime_error("Internal error");
	}

	// https://stackoverflow.com/a/123995
	Data_Type& operator()(const int dir) {
		return const_cast<Data_Type&>(static_cast<const Face_Type<Data_Type>&>(*this).operator()(dir));
	}

	/*! Returns data of given cell face.
	dim: 0 == cell face with normal in x direction, 1 == y, 2 == z,
	side: -1 == negative side of cell from center, +1 == positive
	*/
	const Data_Type& operator()(
		const size_t dim,
		const int side
	) const {
		if (dim > 2) {
			throw std::domain_error("Invalid dimension, must be 0..2 but is " + std::to_string(dim));
		}
		if (side != -1 and side != +1) {
			throw std::domain_error("Invalid side, must be -1,1 but is " + std::to_string(side));
		}
		if        (dim == 0) {
			if (side < 0) return this->face[0];
			else          return this->face[1];
		} else if (dim == 1) {
			if (side < 0) return this->face[2];
			else          return this->face[3];
		} else if (dim == 2) {
			if (side < 0) return this->face[4];
			else          return this->face[5];
		} else throw std::runtime_error("Internal error");
	}

	Face_Type<Data_Type>& operator=(const Face_Type<Data_Type>& other) noexcept {
		if (this == &other) {
			return *this;
		}
		this->face = other.face;
		return *this;
	}

	decltype(face)& operator=(const decltype(face)& other) noexcept {
		this->face = other;
		return this->face;
	}

	decltype(face)& operator=(const std::initializer_list<Data_Type>& other) {
		if (other.size() != this->face.size()) {
			throw std::runtime_error(__FILE__"(" + std::to_string(__LINE__) + ")");
		}
		std::copy(other.begin(), other.end(), this->face.begin());
		return this->face;
	}

	#ifdef MPI_VERSION
	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype() const {
		if constexpr (std::is_same_v<Data_Type, double>) {
			return std::make_tuple((void*) this->face.data(), this->face.size(), MPI_DOUBLE);
		} else if constexpr (std::is_same_v<Data_Type, float>) {
			return std::make_tuple((void*) this->face.data(), this->face.size(), MPI_FLOAT);
		} else if constexpr (std::is_same_v<Data_Type, bool>) {
			return std::make_tuple((void*) this->face.data(), this->face.size(), MPI_CXX_BOOL);
		} else {
			static_assert(true, "Unsupported face item type for MPI");
			//return std::make_tuple(nullptr, 0, MPI_BYTE);
		}
	}
	#endif
};

/*! Records which of cell's faces are primary.

By default positive faces are primary, i.e. values in
corresponding negative face of neighbor(s) in positive
direction from cell are ignored.

Faces of smaller neighbors are always primary.

0 == -x, 1 == +x, ..., 5 == +z
*/
struct Is_Primary_Face {
	using data_type = Face_Type<bool>;
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
				Faces(*cell.data)(n) = false;
				continue;
			}
			if (n < 0) Faces(*cell.data)(n) = false;
		}
	}
}


/*! Common accessors for values stored on cell edges.
*/
template<class Data_Type> struct Edge_Type {
	std::array<Data_Type, 12> edge;

	/*
	par_dim_i = dimension which edge is parallel to,
	0 = x, 1 = y, 2 = z.
	first_perp_dim_i = negative or positive side of cell in
	lexically earlier dimension perpendicular to edge,
	e.g. if par_dim_i = 0, first_perp_dim_i = 1 means positive side of
	cell in y dimension.
	second_perp_dim = neg or pos side in lexically later dimension,
	if par_dim_i = 2, second_perp_dim_i = 0 means negative side of
	cell in y dimension.
	par_dim_i | first_perp_dim_i | second..._i | edge of cell
	    0     |         0        |      0      | x directed: -y, -z sides
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

	Edge_Type<Data_Type>& operator=(const Edge_Type<Data_Type>& other) {
		if (this == &other) {
			return *this;
		}
		this->edge = other.edge;
		return *this;
	}

	decltype(edge)& operator=(const decltype(edge)& other) {
		this->edge = other;
		return this->edge;
	}

	decltype(edge)& operator=(const std::initializer_list<Data_Type>& other) {
		if (other.size() != this->edge.size()) {
			throw std::runtime_error(__FILE__"(" + std::to_string(__LINE__) + ")");
		}
		std::copy(other.begin(), other.end(), this->edge.begin());
		return this->edge;
	}

	#ifdef MPI_VERSION
	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype() const {
		if constexpr (std::is_same_v<Data_Type, double>) {
			return std::make_tuple((void*) this->edge.data(), this->edge.size(), MPI_DOUBLE);
		} else if constexpr (std::is_same_v<Data_Type, float>) {
			return std::make_tuple((void*) this->edge.data(), this->edge.size(), MPI_FLOAT);
		} else if constexpr (std::is_same_v<Data_Type, bool>) {
			return std::make_tuple((void*) this->edge.data(), this->edge.size(), MPI_CXX_BOOL);
		} else {
			static_assert(true, "Unsupported edge item type for MPI");
			//return std::make_tuple(nullptr, 0, MPI_BYTE);
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
	if (grid.get_neighborhood_length() == 0) {
		throw std::domain_error(
			__FILE__ "(" + std::to_string(__LINE__)
			+ "): Neighborhood length must be at least 1");
	}
	for (const auto& cell: cells) {
		Edges(*cell.data) = {
			true, true, true, true, // x-directed edges: -y,-z;-y,+z,...
			true, true, true, true, // y-directed edges: -x,-z;-x,+z,...
			true, true, true, true};// ...

		for (const auto& neighbor: cell.neighbors_of) {
			const auto& fn = neighbor.face_neighbor;
			const auto& en = neighbor.edge_neighbor;

			if (fn != 0 and en[0] >= 0) {
				throw std::runtime_error(__FILE__"(" + std::to_string(__LINE__) + ")");
			}

			if (neighbor.relative_size < 0) {
				continue; // larger neighbor never has priority
			}

			if (fn == -1) {
				// this or another neighbor has priority
				Edges(*cell.data)(1, 0, 0) =
				Edges(*cell.data)(1, 0, 1) =
				Edges(*cell.data)(2, 0, 0) =
				Edges(*cell.data)(2, 0, 1) = false;
				continue;
			}
			if (fn == +1 and neighbor.relative_size > 0) {
				// only smaller neighbors have priority
				Edges(*cell.data)(1, 1, 0) =
				Edges(*cell.data)(1, 1, 1) =
				Edges(*cell.data)(2, 1, 0) =
				Edges(*cell.data)(2, 1, 1) = false;
				continue;
			}
			if (fn == -2) {
				Edges(*cell.data)(0, 0, 0) =
				Edges(*cell.data)(0, 0, 1) =
				Edges(*cell.data)(2, 0, 0) =
				Edges(*cell.data)(2, 1, 0) = false;
				continue;
			}
			if (fn == +2 and neighbor.relative_size > 0) {
				Edges(*cell.data)(0, 1, 0) =
				Edges(*cell.data)(0, 1, 1) =
				Edges(*cell.data)(2, 0, 1) =
				Edges(*cell.data)(2, 1, 1) = false;
				continue;
			}
			if (fn == -3) {
				Edges(*cell.data)(0, 0, 0) =
				Edges(*cell.data)(0, 1, 0) =
				Edges(*cell.data)(1, 0, 0) =
				Edges(*cell.data)(1, 1, 0) = false;
				continue;
			}
			if (fn == +3 and neighbor.relative_size > 0) {
				Edges(*cell.data)(0, 0, 1) =
				Edges(*cell.data)(0, 1, 1) =
				Edges(*cell.data)(1, 0, 1) =
				Edges(*cell.data)(1, 1, 1) = false;
				continue;
			}

			if (en[0] < 0) {
				continue;
			}

			if (neighbor.relative_size > 0 or en[1] == 0) {
				Edges(*cell.data)(en[0], en[1], en[2]) = false;
			}
		}
	}
}


//! Returns smallest [0] and largest [1] refinement level for given cell.
std::array<int, 2> get_target_refinement_level(
	Options& options,
	const double& t,
	const std::array<double, 3>& cell_center
) {
	const auto& c = cell_center;
	const auto
		r = std::sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]),
		lat = std::asin(c[2] / r),
		lon = std::atan2(c[1], c[0]);

	const auto min_ref_lvl = options.get_ref_lvl_at_least(
		t, c[0], c[1], c[2], r, lat, lon);
	if (min_ref_lvl < 0) {
		throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): options returned min_ref_lvl < 0");
	}

	const auto max_ref_lvl = options.get_ref_lvl_at_most(
		t, c[0], c[1], c[2], r, lat, lon);
	if (max_ref_lvl < 0) {
		throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): options returned min_ref_lvl < 0");
	}
	if (max_ref_lvl < min_ref_lvl) {
		throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): options returned max_ref_lvl < min_ref_lvl");
	}

	std::array<int, 2> ret_val{min_ref_lvl, max_ref_lvl};
	return ret_val;
}

// for adapt_grid() below
struct Default_Cells_Handler {
	template <class Grid, class Cells> void operator()(const Grid&, const Cells&) const {
		return;
	}
};

/*! Adapts given grid based on given geometrical options

Handlers must have following operator() signature:
nch(grid, new_cells);
rch(grid, removed_cells);
where
const auto new_cells = grid.stop_refining();
const auto removed_cells = grid.get_removed_cells();
*/
template <
	class Grid,
	class New_Cells_Handler = Default_Cells_Handler,
	class Removed_Cells_Handler = Default_Cells_Handler
> void adapt_grid(
	Grid& grid,
	Options& options,
	const double& sim_time,
	const New_Cells_Handler& nch = Default_Cells_Handler(),
	const Removed_Cells_Handler& rch = Default_Cells_Handler()
) {
	for (int i = 0; i < grid.get_maximum_refinement_level(); i++) {
		for (const auto& cell: grid.local_cells()) {
			const auto c = grid.geometry.get_center(cell.id);
			const auto ref_lvl = grid.get_refinement_level(cell.id);
			const auto tgt_ref_lvl = get_target_refinement_level(options, sim_time, c);
			if (ref_lvl < tgt_ref_lvl[0]) {
				grid.refine_completely(cell.id);
			} else if (ref_lvl == tgt_ref_lvl[0]) {
				grid.dont_unrefine(cell.id);
			} else if (ref_lvl > tgt_ref_lvl[1]) {
				grid.unrefine_completely(cell.id);
			}
		}
		const auto new_cells = grid.stop_refining();
		const auto removed_cells = grid.get_removed_cells();
		if (new_cells.size() == 0 and removed_cells.size() == 0) {
			break;
		} else {
			nch(grid, new_cells);
			rch(grid, removed_cells);
			grid.clear_refined_unrefined_data();
		}
	}
}

}} // namespaces

#endif // ifndef PAMHD_GRID_AMR_HPP
