/*
Grid stuff related to adaptive mesh refinement of PAMHD.

Copyright 2023, 2024 Finnish Meteorological Institute
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


Author(s): Ilja Honkonen
*/

#ifndef PAMHD_GRID_AMR_HPP
#define PAMHD_GRID_AMR_HPP


#include "algorithm"
#include "initializer_list"
#include "set"
#include "stdexcept"
#include "string"
#include "tuple"
#include "type_traits"

#include "grid/options.hpp"


namespace pamhd {
namespace grid {


/*! Variable indicating minimum refinement level towards which grid cells should be refined.
*/
struct Target_Refinement_Level_Min {
	using data_type = int;
	static const std::string get_name() { return {"target minimum refinement level"}; }
	static const std::string get_option_name() { return {"target-min-ref-lvl"}; }
	static const std::string get_option_help() { return {"Target minimum refinement level of cell"}; }
};

struct Target_Refinement_Level_Max {
	using data_type = int;
	static const std::string get_name() { return {"target maximum refinement level"}; }
	static const std::string get_option_name() { return {"target-max-ref-lvl"}; }
	static const std::string get_option_help() { return {"Target maximum refinement level of cell"}; }
};


/*! Common accessors for values stored on cell faces.
*/
template<class Data_Type> struct Face_Type {
	static constexpr size_t nr_faces = 6;
	std::array<Data_Type, nr_faces> face;

	/*! Returns data of given cell face.
	dir: -1 == -x, +1 == +x, -2 == -y, ..., +3 == +z face
	*/
	const Data_Type& operator()(const int dir) const {
		using std::to_string;

		if (dir == 0 or dir < -3 or dir > +3) {
			throw std::domain_error(__FILE__ "(" + to_string(__LINE__) + "): Invalid direction, must be -3,-2,-1,1,2,3 but is " + to_string(dir));
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
		using std::domain_error;
		using std::to_string;

		if (dim > 2) {
			throw domain_error("Invalid dimension, must be 0..2 but is " + to_string(dim));
		}
		if (side != -1 and side != +1) {
			throw domain_error("Invalid side, must be -1,1 but is " + to_string(side));
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

	Data_Type& operator()(const size_t dim, const int side) {
		return const_cast<Data_Type&>(static_cast<const Face_Type<Data_Type>&>(*this).operator()(dim, side));
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
		using std::array;
		using std::is_same_v;
		using std::make_tuple;
		using std::runtime_error;

		if        constexpr (is_same_v<Data_Type, double>) {
			return make_tuple((void*) this->face.data(), this->face.size(), MPI_DOUBLE);
		} else if constexpr (is_same_v<Data_Type, float>) {
			return make_tuple((void*) this->face.data(), this->face.size(), MPI_FLOAT);
		} else if constexpr (is_same_v<Data_Type, uint64_t>) {
			return make_tuple((void*) this->face.data(), this->face.size(), MPI_UINT64_T);
		} else if constexpr (is_same_v<Data_Type, uint32_t>) {
			return make_tuple((void*) this->face.data(), this->face.size(), MPI_UINT32_T);
		} else if constexpr (is_same_v<Data_Type, uint16_t>) {
			return make_tuple((void*) this->face.data(), this->face.size(), MPI_UINT16_T);
		} else if constexpr (is_same_v<Data_Type, uint8_t>) {
			return make_tuple((void*) this->face.data(), this->face.size(), MPI_UINT8_T);
		} else if constexpr (is_same_v<Data_Type, int64_t>) {
			return make_tuple((void*) this->face.data(), this->face.size(), MPI_INT64_T);
		} else if constexpr (is_same_v<Data_Type, int32_t>) {
			return make_tuple((void*) this->face.data(), this->face.size(), MPI_INT32_T);
		} else if constexpr (is_same_v<Data_Type, int16_t>) {
			return make_tuple((void*) this->face.data(), this->face.size(), MPI_INT16_T);
		} else if constexpr (is_same_v<Data_Type, int8_t>) {
			return make_tuple((void*) this->face.data(), this->face.size(), MPI_INT8_T);
		} else if constexpr (is_same_v<Data_Type, long long>) {
			return make_tuple((void*) this->face.data(), this->face.size(), MPI_LONG_LONG);
		} else if constexpr (is_same_v<Data_Type, long>) {
			return make_tuple((void*) this->face.data(), this->face.size(), MPI_LONG);
		} else if constexpr (is_same_v<Data_Type, int>) {
			return make_tuple((void*) this->face.data(), this->face.size(), MPI_INT);
		} else if constexpr (is_same_v<Data_Type, short>) {
			return make_tuple((void*) this->face.data(), this->face.size(), MPI_SHORT);
		} else if constexpr (is_same_v<Data_Type, char>) {
			return make_tuple((void*) this->face.data(), this->face.size(), MPI_CHAR);
		} else if constexpr (is_same_v<Data_Type, signed char>) {
			return make_tuple((void*) this->face.data(), this->face.size(), MPI_SIGNED_CHAR);
		} else if constexpr (is_same_v<Data_Type, unsigned char>) {
			return make_tuple((void*) this->face.data(), this->face.size(), MPI_UNSIGNED_CHAR);
		} else if constexpr (is_same_v<Data_Type, bool>) {
			return make_tuple((void*) this->face.data(), this->face.size(), MPI_CXX_BOOL);
		#ifdef EIGEN_WORLD_VERSION
		} else if constexpr (is_same_v<Data_Type, Eigen::Vector3d>) {
			array<void*, nr_faces> addresses;
			array<int, nr_faces> counts;
			array<MPI_Datatype, nr_faces> datatypes;
			array<MPI_Aint, nr_faces> displacements;

			MPI_Datatype final_datatype = MPI_DATATYPE_NULL;
			for (size_t i = 0; i < nr_faces; i++) {
				addresses[i] = (void*)this->face[i].data();
				counts[i] = 3;
				datatypes[i] = MPI_DOUBLE;
				displacements[i]
					= static_cast<char*>(addresses[i])
					- static_cast<char*>(addresses[0]);
			}
			const auto result = MPI_Type_create_struct(
				int(counts.size()),
				counts.data(),
				displacements.data(),
				datatypes.data(),
				&final_datatype
			);
			if (result != MPI_SUCCESS) {
				throw runtime_error("Couldn't create MPI datatype for MHD flux");
			}
			return make_tuple(addresses[0], 1, final_datatype);
		#endif
		} else { // assume Data_Type has get_mpi_datatype()
			array<void*, nr_faces> addresses;
			array<int, nr_faces> counts;
			array<MPI_Datatype, nr_faces> datatypes;
			array<MPI_Aint, nr_faces> displacements;

			MPI_Datatype final_datatype = MPI_DATATYPE_NULL;
			for (size_t i = 0; i < nr_faces; i++) {
				std::tie(
					addresses[i],
					counts[i],
					datatypes[i]
				) = this->face[i].get_mpi_datatype();
				displacements[i]
					= static_cast<char*>(addresses[i])
					- static_cast<char*>(addresses[0]);
			}
			auto result = MPI_Type_create_struct(
				int(counts.size()),
				counts.data(),
				displacements.data(),
				datatypes.data(),
				&final_datatype
			);
			if (result != MPI_SUCCESS) {
				throw runtime_error("Couldn't create MPI datatype for MHD flux");
			}

			// free user-defined component datatypes
			for (size_t i = 0; i < nr_faces; i++) {
				if (datatypes[i] == MPI_DATATYPE_NULL) {
					continue;
				}
				int combiner = -1, tmp1 = -1, tmp2 = -1, tmp3 = -1;
				MPI_Type_get_envelope(datatypes[i], &tmp1, &tmp2, &tmp3, &combiner);
				if (combiner != MPI_COMBINER_NAMED) {
					MPI_Type_free(&datatypes[i]);
				}
			}

			return make_tuple(addresses[0], 1, final_datatype);
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
	e.g. if par_dim_i = -1, first_perp_dim_i = +1 means positive
	side of cell in y dimension.
	second_perp_dim = neg or pos side in lexically later dimension,
	if par_dim_i = 2, second_perp_dim_i = -1 means negative side of
	cell in y dimension.
	par_dim_i | first_perp_dim_i | second..._i | edge of cell
	    0     |        -1        |     -1      | x directed: -y, -z sides
	    0     |        -1        |     +1      | x dir:      -y, +z
	    0     |        +1        |     +1      | x dir:      +y, +z
	...
	    2     |        +1        |     -1      | z dir:      +x, -y
	    2     |        +1        |     +1      | z dir:      +x, +y
	*/
	const Data_Type& operator()(
		const int par_dim_i,
		const int first_perp_dim_i,
		const int second_perp_dim_i
	) const {
		using std::domain_error;
		using std::to_string;

		if (par_dim_i < 0 or par_dim_i > 2) {
			throw domain_error("Parallel dimension != 0,1,2: " + to_string(par_dim_i));
		}
		if (first_perp_dim_i != +1 and first_perp_dim_i != -1) {
			throw domain_error("First perpendicular dimension != +-1: " + to_string(first_perp_dim_i));
		}
		if (second_perp_dim_i != +1 and second_perp_dim_i != -1) {
			throw domain_error("Second perpendicular dimension != +-1: " + to_string(second_perp_dim_i));
		}
		const size_t
			p1 = [&](){if (first_perp_dim_i < 0) return 0; else return 1;}(),
			p2 = [&](){if (second_perp_dim_i < 0) return 0; else return 1;}();
		return this->edge[size_t(par_dim_i)*2*2 + p1*2 + p2];
	}

	// https://stackoverflow.com/a/123995
	Data_Type& operator()(
		const int par_dim_i,
		const int first_perp_dim_i,
		const int second_perp_dim_i
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
		using std::is_same_v;
		using std::make_tuple;

		if constexpr (is_same_v<Data_Type, double>) {
			return make_tuple((void*) this->edge.data(), this->edge.size(), MPI_DOUBLE);
		} else if constexpr (is_same_v<Data_Type, float>) {
			return make_tuple((void*) this->edge.data(), this->edge.size(), MPI_FLOAT);
		} else if constexpr (is_same_v<Data_Type, uint64_t>) {
			return make_tuple((void*) this->edge.data(), this->edge.size(), MPI_UINT64_T);
		} else if constexpr (is_same_v<Data_Type, uint32_t>) {
			return make_tuple((void*) this->edge.data(), this->edge.size(), MPI_UINT32_T);
		} else if constexpr (is_same_v<Data_Type, uint16_t>) {
			return make_tuple((void*) this->edge.data(), this->edge.size(), MPI_UINT16_T);
		} else if constexpr (is_same_v<Data_Type, uint8_t>) {
			return make_tuple((void*) this->edge.data(), this->edge.size(), MPI_UINT8_T);
		} else if constexpr (is_same_v<Data_Type, int64_t>) {
			return make_tuple((void*) this->edge.data(), this->edge.size(), MPI_INT64_T);
		} else if constexpr (is_same_v<Data_Type, int32_t>) {
			return make_tuple((void*) this->edge.data(), this->edge.size(), MPI_INT32_T);
		} else if constexpr (is_same_v<Data_Type, int16_t>) {
			return make_tuple((void*) this->edge.data(), this->edge.size(), MPI_INT16_T);
		} else if constexpr (is_same_v<Data_Type, int8_t>) {
			return make_tuple((void*) this->edge.data(), this->edge.size(), MPI_INT8_T);
		} else if constexpr (is_same_v<Data_Type, long long>) {
			return make_tuple((void*) this->edge.data(), this->edge.size(), MPI_LONG_LONG);
		} else if constexpr (is_same_v<Data_Type, long>) {
			return make_tuple((void*) this->edge.data(), this->edge.size(), MPI_LONG);
		} else if constexpr (is_same_v<Data_Type, int>) {
			return make_tuple((void*) this->edge.data(), this->edge.size(), MPI_INT);
		} else if constexpr (is_same_v<Data_Type, short>) {
			return make_tuple((void*) this->edge.data(), this->edge.size(), MPI_SHORT);
		} else if constexpr (is_same_v<Data_Type, char>) {
			return make_tuple((void*) this->edge.data(), this->edge.size(), MPI_CHAR);
		} else if constexpr (is_same_v<Data_Type, signed char>) {
			return make_tuple((void*) this->edge.data(), this->edge.size(), MPI_SIGNED_CHAR);
		} else if constexpr (is_same_v<Data_Type, unsigned char>) {
			return make_tuple((void*) this->edge.data(), this->edge.size(), MPI_UNSIGNED_CHAR);
		} else if constexpr (is_same_v<Data_Type, bool>) {
			return make_tuple((void*) this->edge.data(), this->edge.size(), MPI_CXX_BOOL);
		} else {
			static_assert(false, "Unsupported edge item type for MPI");
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
				Edges(*cell.data)(1,-1,-1) =
				Edges(*cell.data)(1,-1,+1) =
				Edges(*cell.data)(2,-1,-1) =
				Edges(*cell.data)(2,-1,+1) = false;
				continue;
			}
			if (fn == +1 and neighbor.relative_size > 0) {
				// only smaller neighbors have priority
				Edges(*cell.data)(1,+1,-1) =
				Edges(*cell.data)(1,+1,+1) =
				Edges(*cell.data)(2,+1,-1) =
				Edges(*cell.data)(2,+1,+1) = false;
				continue;
			}
			if (fn == -2) {
				Edges(*cell.data)(0,-1,-1) =
				Edges(*cell.data)(0,-1,+1) =
				Edges(*cell.data)(2,-1,-1) =
				Edges(*cell.data)(2,+1,-1) = false;
				continue;
			}
			if (fn == +2 and neighbor.relative_size > 0) {
				Edges(*cell.data)(0,+1,-1) =
				Edges(*cell.data)(0,+1,+1) =
				Edges(*cell.data)(2,-1,+1) =
				Edges(*cell.data)(2,+1,+1) = false;
				continue;
			}
			if (fn == -3) {
				Edges(*cell.data)(0,-1,-1) =
				Edges(*cell.data)(0,+1,-1) =
				Edges(*cell.data)(1,-1,-1) =
				Edges(*cell.data)(1,+1,-1) = false;
				continue;
			}
			if (fn == +3 and neighbor.relative_size > 0) {
				Edges(*cell.data)(0,-1,+1) =
				Edges(*cell.data)(0,+1,+1) =
				Edges(*cell.data)(1,-1,+1) =
				Edges(*cell.data)(1,+1,+1) = false;
				continue;
			}

			if (en[0] < 0) {
				continue;
			}

			if (neighbor.relative_size > 0 or en[1] == -1) {
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
	using std::runtime_error;
	using std::to_string;

	const auto& c = cell_center;
	const auto
		r = std::sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]),
		lat = std::asin(c[2] / r),
		lon = std::atan2(c[1], c[0]);

	const auto min_ref_lvl = options.get_ref_lvl_at_least(
		t, c[0], c[1], c[2], r, lat, lon);
	if (min_ref_lvl < 0) {
		throw runtime_error(__FILE__ "(" + to_string(__LINE__) + "): options returned min_ref_lvl < 0");
	}

	const auto max_ref_lvl = options.get_ref_lvl_at_most(
		t, c[0], c[1], c[2], r, lat, lon);
	if (max_ref_lvl < 0) {
		throw runtime_error(__FILE__ "(" + to_string(__LINE__) + "): options returned min_ref_lvl < 0");
	}
	if (max_ref_lvl < min_ref_lvl) {
		throw runtime_error(__FILE__ "(" + to_string(__LINE__) + "): options returned max_ref_lvl < min_ref_lvl");
	}

	std::array<int, 2> ret_val{min_ref_lvl, max_ref_lvl};
	return ret_val;
}

// Handlers that only handle target refinement levels.
template<
	class Target_Refinement_Level_Min_Getter,
	class Target_Refinement_Level_Max_Getter
> struct New_Cells_Handler {
	const Target_Refinement_Level_Min_Getter& RLMin;
	const Target_Refinement_Level_Max_Getter& RLMax;

	New_Cells_Handler(
		const Target_Refinement_Level_Min_Getter& RLMin_,
		const Target_Refinement_Level_Max_Getter& RLMax_
	) :
		RLMin(RLMin_), RLMax(RLMax_)
	{};

	template <
		class Grid, class Cells
	> void operator()(
		const Grid& grid, const Cells& new_cells
	) const {
		using std::runtime_error;
		using std::to_string;
		for (const auto& new_cell_id: new_cells) {
			auto* const cell_data = grid[new_cell_id];
			if (cell_data == nullptr) {
				throw runtime_error(__FILE__ ":" + to_string(__LINE__));
			}
			const auto parent_id = grid.mapping.get_parent(new_cell_id);
			auto* const parent_data = grid[parent_id];
			if (parent_data == nullptr) {
				throw runtime_error(__FILE__ ":" + to_string(__LINE__));
			}
			RLMin(*cell_data) = RLMin(*parent_data);
			RLMax(*cell_data) = RLMax(*parent_data);
		}
	}
};

template<
	class Target_Refinement_Level_Min_Getter,
	class Target_Refinement_Level_Max_Getter
> struct Removed_Cells_Handler {
	const Target_Refinement_Level_Min_Getter& RLMin;
	const Target_Refinement_Level_Max_Getter& RLMax;

	Removed_Cells_Handler(
		const Target_Refinement_Level_Min_Getter& RLMin_,
		const Target_Refinement_Level_Max_Getter& RLMax_
	) :
		RLMin(RLMin_), RLMax(RLMax_)
	{};

	template <
		class Grid, class Cells
	> void operator()(
		const Grid& grid, const Cells& removed_cells
	) const {
		using std::runtime_error;
		using std::to_string;

		// process each parent of removed cells only once
		std::set<uint64_t> parents;
		for (const auto& removed_cell: removed_cells) {
			parents.insert(grid.mapping.get_parent(removed_cell));
		}

		// initialize data of parents of removed cells
		for (const auto& parent: parents) {
			auto* const parent_data = grid[parent];
			if (parent_data == nullptr) {
				throw runtime_error(__FILE__ ":" + to_string(__LINE__));
			}
			RLMin(*parent_data) = 999;
			RLMax(*parent_data) = 0;
		}

		// parent's refinement range spans all of childrens'
		for (const auto& removed_cell_id: removed_cells) {
			auto* const removed_cell_data = grid[removed_cell_id];
			if (removed_cell_data == nullptr) {
				throw runtime_error(__FILE__ ":" + to_string(__LINE__));
			}
			const auto parent_id = grid.mapping.get_parent(removed_cell_id);
			auto* const parent_data = grid[parent_id];
			if (parent_data == nullptr) {
				throw runtime_error(__FILE__ ":" + to_string(__LINE__));
			}
			RLMin(*parent_data) = std::min(RLMin(*removed_cell_data), RLMin(*parent_data));
			RLMax(*parent_data) = std::max(RLMax(*removed_cell_data), RLMax(*parent_data));
		}
	}
};

/*! Adapts given grid based on target min and max refinement level.

Grid is adapted at most grid.get_maximum_refinement_level() times.
Cells whose refinement level is smaller than target minimum are
refined.
Cells whose refinement level is larger than target maximum are
unrefined.

Handlers must have following operator() signature:
nch(grid, new_cells);
rch(grid, removed_cells);
where
const auto new_cells = grid.stop_refining();
const auto removed_cells = grid.get_removed_cells();
*/
template <
	class Grid,
	class Target_Refinement_Level_Min_Getter,
	class Target_Refinement_Level_Max_Getter,
	class New_Cells_Handler,
	class Removed_Cells_Handler
> void adapt_grid(
	Grid& grid,
	const Target_Refinement_Level_Min_Getter& RLMin,
	const Target_Refinement_Level_Max_Getter& RLMax,
	const New_Cells_Handler& nch,
	const Removed_Cells_Handler& rch,
	const int max_rounds = 1 << 30
) {
	#ifdef MPI_VERSION
	MPI_Comm comm = grid.get_communicator();
	#endif
	for (
		int i = 0;
		i < std::min(max_rounds, grid.get_maximum_refinement_level());
		i++
	) {
		for (const auto& cell: grid.local_cells()) {
			const auto ref_lvl = grid.get_refinement_level(cell.id);
			if (ref_lvl < RLMin(*cell.data)) {
				grid.refine_completely(cell.id);
			} else if (ref_lvl == RLMin(*cell.data)) {
				grid.dont_unrefine(cell.id);
			} else if (ref_lvl > RLMax(*cell.data)) {
				grid.unrefine_completely(cell.id);
			}
		}
		const auto new_cells = grid.stop_refining();
		const auto removed_cells = grid.get_removed_cells();

		const uint64_t total_size_local = new_cells.size() + removed_cells.size();
		uint64_t total_size_global = total_size_local;
		#ifdef MPI_VERSION
		total_size_global = dccrg::All_Reduce()(total_size_local, comm);
		#endif

		if (total_size_global == 0) {
			break;
		} else {
			if (new_cells.size() > 0) nch(grid, new_cells);
			if (removed_cells.size() > 0) rch(grid, removed_cells);
			grid.clear_refined_unrefined_data();
		}
	}
	#ifdef MPI_VERSION
	MPI_Comm_free(&comm);
	#endif
}


/*! Sets min,max refinement levels based on geometry.

If clamped == true then target minimum refinement level
is not decreased and tgt max ref lvl is not increased.
*/
template <
	class Cells,
	class Grid,
	class Target_Refinement_Level_Min_Getter,
	class Target_Refinement_Level_Max_Getter
> void set_minmax_refinement_level(
	const Cells& cells,
	Grid& grid,
	Options& options,
	const double& sim_time,
	const Target_Refinement_Level_Min_Getter& RLMin,
	const Target_Refinement_Level_Max_Getter& RLMax,
	bool clamped
) {
	using std::clamp;

	for (const auto& cell: cells) {
		const auto c = grid.geometry.get_center(cell.id);
		const auto [rlmin, rlmax] = get_target_refinement_level(options, sim_time, c);
		if (clamped) {
			const auto
				prev_min = RLMin(*cell.data),
				prev_max = RLMax(*cell.data);
			RLMin(*cell.data) = clamp(rlmin, prev_min, prev_max);
			RLMax(*cell.data) = clamp(rlmax, prev_min, prev_max);
		} else {
			RLMin(*cell.data) = rlmin;
			RLMax(*cell.data) = rlmax;
		}
	}
}

}} // namespaces

#endif // ifndef PAMHD_GRID_AMR_HPP
