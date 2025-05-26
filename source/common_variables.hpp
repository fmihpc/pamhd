/*
Variables of PAMHD common to several test programs.

Copyright 2014, 2015, 2016, 2017 Ilja Honkonen
Copyright 2023, 2024, 2025 Finnish Meteorological Institute
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

#ifndef PAMHD_VARIABLES_HPP
#define PAMHD_VARIABLES_HPP


#include "array"


namespace pamhd {


/*! Helper for simulation variables stored on cell faces.

Stores given variable on all cell faces.

When compiled with MPI provides get_mpi_datatype()
covering variable of every face.
*/
template<class Data_Type> struct Face_Type {
	static constexpr size_t nr_faces = 6;
	std::array<Data_Type, nr_faces> face;

	/*! Returns data of given cell face.
	dir: -1 == -x, +1 == +x, -2 == -y, ..., +3 == +z face
	*/
	const Data_Type& operator()(const int& dir) const {
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
	Data_Type& operator()(const int& dir) {
		return const_cast<Data_Type&>(static_cast<const Face_Type<Data_Type>&>(*this).operator()(dir));
	}

	/*! Returns data of given cell face.
	dim: 0 == cell face with normal in x direction, 1 == y, 2 == z,
	side: -1 == negative side of cell from center, +1 == positive
	*/
	const Data_Type& operator()(
		const size_t& dim,
		const int& side
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

	Data_Type& operator()(const size_t& dim, const int& side) {
		return const_cast<Data_Type&>(static_cast<const Face_Type<Data_Type>&>(*this).operator()(dim, side));
	}

	Face_Type<Data_Type>& operator=(
		const Face_Type<Data_Type>& other
	) noexcept {
		if (this == &other) {
			return *this;
		}
		this->face = other.face;
		return *this;
	}

	decltype(Face_Type::face)& operator=(
		const decltype(Face_Type::face)& other
	) noexcept {
		this->face = other;
		return this->face;
	}

	decltype(Face_Type::face)& operator=(
		const std::initializer_list<Data_Type>& other
	) {
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
		} else if constexpr (is_same_v<Data_Type, array<double, 3>>) {
			return make_tuple((void*) this->face.data(), 3*this->face.size(), MPI_DOUBLE);
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


/*! Helper for simulation variables stored on cell edges.

Stores given variable on all cell edges.

When compiled with MPI provides get_mpi_datatype()
covering variable of every edge.
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
		const int& par_dim_i,
		const int& first_perp_dim_i,
		const int& second_perp_dim_i
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
		const int& par_dim_i,
		const int& first_perp_dim_i,
		const int& second_perp_dim_i
	) {
		return const_cast<Data_Type&>(static_cast<const Edge_Type<Data_Type>&>(*this).operator()(par_dim_i, first_perp_dim_i, second_perp_dim_i));
	}

	Edge_Type<Data_Type>& operator=(
		const Edge_Type<Data_Type>& other
	) noexcept {
		if (this == &other) {
			return *this;
		}
		this->edge = other.edge;
		return *this;
	}

	decltype(Edge_Type::edge)& operator=(
		const decltype(Edge_Type::edge)& other
	) noexcept {
		this->edge = other;
		return this->edge;
	}

	decltype(Edge_Type::edge)& operator=(
		const std::initializer_list<Data_Type>& other
	) {
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


/*! Helper for simulation variables stored on cell vertices.

Stores given variable on all cell vertices.

When compiled with MPI provides get_mpi_datatype()
covering variable of every vertex.
*/
template<class Data_Type> struct Vertex_Type {
	std::array<Data_Type, 8> vertex;

	/*
	Vertex defined by direction from cell center in each dimension.

	Negative values mean vertex negative side from cell center,
	positive values mean positive side.
	*/
	const Data_Type& operator()(
		const int& x_offset,
		const int& y_offset,
		const int& z_offset
	) const {
		using std::domain_error;

		if (x_offset == 0) {
			throw domain_error("x offset of vertex cannot be 0");
		}
		if (y_offset == 0) {
			throw domain_error("y offset of vertex cannot be 0");
		}
		if (z_offset == 0) {
			throw domain_error("z offset of vertex cannot be 0");
		}
		const size_t
			i = [&](){if (x_offset < 0) return 0; else return 1;}(),
			j = [&](){if (y_offset < 0) return 0; else return 1;}(),
			k = [&](){if (z_offset < 0) return 0; else return 1;}();
		return this->vertex[k*2*2 + j*2 + i];
	}

	// https://stackoverflow.com/a/123995
	Data_Type& operator()(
		const int& x_offset,
		const int& y_offset,
		const int& z_offset
	) {
		return const_cast<Data_Type&>(
			static_cast<const Vertex_Type<Data_Type>&>(*this)
				.operator()(x_offset, y_offset, z_offset));
	}

	const Data_Type& operator()(
		const std::array<int, 3>& offsets
	) const {
		return this->operator()(offsets[0], offsets[1], offsets[2]);
	}

	Data_Type& operator()(
		const std::array<int, 3>& offsets
	) {
		return const_cast<Data_Type&>(
			static_cast<const Vertex_Type<Data_Type>&>(*this)
				.operator()(offsets));
	}

	Vertex_Type<Data_Type>& operator=(
		const Vertex_Type<Data_Type>& other
	) noexcept {
		if (this == &other) {
			return *this;
		}
		this->vertex = other.vertex;
		return *this;
	}

	decltype(Vertex_Type::vertex)& operator=(
		const decltype(Vertex_Type::vertex)& other
	) noexcept {
		this->vertex = other;
		return this->vertex;
	}

	decltype(Vertex_Type::vertex)& operator=(
		const std::initializer_list<Data_Type>& other
	) {
		if (other.size() != this->vertex.size()) {
			throw std::runtime_error(__FILE__"(" + std::to_string(__LINE__) + ")");
		}
		std::copy(other.begin(), other.end(), this->vertex.begin());
		return this->vertex;
	}

	#ifdef MPI_VERSION
	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype() const {
		using std::is_same_v;
		using std::make_tuple;

		if constexpr (is_same_v<Data_Type, double>) {
			return make_tuple((void*) this->vertex.data(), this->vertex.size(), MPI_DOUBLE);
		} else if constexpr (is_same_v<Data_Type, float>) {
			return make_tuple((void*) this->vertex.data(), this->vertex.size(), MPI_FLOAT);
		} else if constexpr (is_same_v<Data_Type, uint64_t>) {
			return make_tuple((void*) this->vertex.data(), this->vertex.size(), MPI_UINT64_T);
		} else if constexpr (is_same_v<Data_Type, uint32_t>) {
			return make_tuple((void*) this->vertex.data(), this->vertex.size(), MPI_UINT32_T);
		} else if constexpr (is_same_v<Data_Type, uint16_t>) {
			return make_tuple((void*) this->vertex.data(), this->vertex.size(), MPI_UINT16_T);
		} else if constexpr (is_same_v<Data_Type, uint8_t>) {
			return make_tuple((void*) this->vertex.data(), this->vertex.size(), MPI_UINT8_T);
		} else if constexpr (is_same_v<Data_Type, int64_t>) {
			return make_tuple((void*) this->vertex.data(), this->vertex.size(), MPI_INT64_T);
		} else if constexpr (is_same_v<Data_Type, int32_t>) {
			return make_tuple((void*) this->vertex.data(), this->vertex.size(), MPI_INT32_T);
		} else if constexpr (is_same_v<Data_Type, int16_t>) {
			return make_tuple((void*) this->vertex.data(), this->vertex.size(), MPI_INT16_T);
		} else if constexpr (is_same_v<Data_Type, int8_t>) {
			return make_tuple((void*) this->vertex.data(), this->vertex.size(), MPI_INT8_T);
		} else if constexpr (is_same_v<Data_Type, long long>) {
			return make_tuple((void*) this->vertex.data(), this->vertex.size(), MPI_LONG_LONG);
		} else if constexpr (is_same_v<Data_Type, long>) {
			return make_tuple((void*) this->vertex.data(), this->vertex.size(), MPI_LONG);
		} else if constexpr (is_same_v<Data_Type, int>) {
			return make_tuple((void*) this->vertex.data(), this->vertex.size(), MPI_INT);
		} else if constexpr (is_same_v<Data_Type, short>) {
			return make_tuple((void*) this->vertex.data(), this->vertex.size(), MPI_SHORT);
		} else if constexpr (is_same_v<Data_Type, char>) {
			return make_tuple((void*) this->vertex.data(), this->vertex.size(), MPI_CHAR);
		} else if constexpr (is_same_v<Data_Type, signed char>) {
			return make_tuple((void*) this->vertex.data(), this->vertex.size(), MPI_SIGNED_CHAR);
		} else if constexpr (is_same_v<Data_Type, unsigned char>) {
			return make_tuple((void*) this->vertex.data(), this->vertex.size(), MPI_UNSIGNED_CHAR);
		} else if constexpr (is_same_v<Data_Type, bool>) {
			return make_tuple((void*) this->vertex.data(), this->vertex.size(), MPI_CXX_BOOL);
		} else if constexpr (is_same_v<Data_Type, std::array<double, 3>>) {
			return make_tuple((void*) this->vertex.data(), 3*this->vertex.size(), MPI_DOUBLE);
		} else {
			static_assert(false, "Unsupported vertex item type for MPI");
		}
	}
	#endif
};


/*
Variables used by magnetohydrodynamic (MHD) solver
*/

struct Magnetic_Field {
	static bool is_stale;
	using data_type = std::array<double, 3>;
	static const std::string get_name() { return {"magnetic field"}; }
	static const std::string get_option_name() { return {"magnetic-field"}; }
	static const std::string get_option_help() { return {"Plasma magnetic field (T)"}; }
};

struct Magnetic_Field_Flux {
	using data_type = std::array<double, 3>;
	static const std::string get_name() { return {"magnetic field flux"}; }
	static const std::string get_option_name() { return {"magnetic-field-flux"}; }
	static const std::string get_option_help() { return {"Flux of magnetic field (T)"}; }
};

//! stores B before divergence removal so B can be restored after failed removal
struct Magnetic_Field_Temp {
	using data_type = std::array<double, 3>;
	static const std::string get_name() { return {"temporary magnetic field"}; }
	static const std::string get_option_name() { return {"temporary-magnetic-field"}; }
	static const std::string get_option_help() { return {"Temporary value of magnetic field in plasma (T)"}; }
};

//! stores change in B due to resistivity
struct Magnetic_Field_Resistive {
	using data_type = std::array<double, 3>;
	static const std::string get_name() { return {"resistive magnetic field"}; }
	static const std::string get_option_name() { return {"resistive-magnetic-field"}; }
	static const std::string get_option_help() { return {"Change in magnetic field due to resistivity"}; }
};

//! Background magnetic field vector at cell faces
struct Bg_Magnetic_Field {
	static bool is_stale;
	using data_type = pamhd::Face_Type<std::array<double, 3>>;
	static const std::string get_name() { return {"face background magnetic fields"}; }
	static const std::string get_option_name() { return {"bg-b"}; }
	static const std::string get_option_help() { return {"background magnetic field vector on cell faces"}; }
};

/*! Magnetic field component at cell faces

Magnetic field stored on cell's faces.
Component of magnetic field normal to the face is stored for each face.
*/
struct Face_Magnetic_Field {
	static bool is_stale;
	using data_type = pamhd::Face_Type<double>;
	static const std::string get_name() { return {"magnetic field on faces"}; }
	static const std::string get_option_name() { return {"face-b"}; }
	static const std::string get_option_help() { return {"magnetic field on cell faces normal to cell faces"}; }
};

struct Face_dB {
	using data_type = pamhd::Face_Type<double>;
};

//! Electric field along each cell edge
struct Edge_Electric_Field {
	using data_type = pamhd::Edge_Type<double>;
	static const std::string get_name() { return {"electric field on edges"}; }
	static const std::string get_option_name() { return {"edge-e"}; }
	static const std::string get_option_help() { return {"electric field on cell edges parallel to cell edges"}; }
};

struct MPI_Rank {
	using data_type = int;
	static const std::string get_name() { return {"MPI rank"}; }
	static const std::string get_option_name() { return {"mpi-rank"}; }
	static const std::string get_option_help() { return {"Owner (MPI process) of cell"}; }
};

struct Magnetic_Field_Divergence {
	using data_type = double;
	static const std::string get_name() { return {"magnetic field divergence"}; }
	static const std::string get_option_name() { return {"magnetic-field-divergence"}; }
	static const std::string get_option_help() { return {"Divergence of plasma magnetic field (T/m)"}; }
};

//! J in J = ∇×B
struct Electric_Current_Density {
	static bool is_stale;
	using data_type = std::array<double, 3>;
	static const std::string get_name() { return {"current density"}; }
	static const std::string get_option_name() { return {"current-density"}; }
	static const std::string get_option_help() { return {"Density of electric current"}; }
};

//! Electrical resistivity of plasma (n in E = -VxB + nJ + ...)
struct Resistivity {
	using data_type = double;
	static const std::string get_name() { return {"electrical resistivity"}; }
	static const std::string get_option_name() { return {"resistivity"}; }
	static const std::string get_option_help() { return {"Electrical resistivity"}; }
};

/*!
Information for solver(s) on how to handle a simulation cell.

-1 means dont_solve cell that's neither read nor written
0 means boundary read-only cell
1 means normal read-write cell
*/
struct Cell_Type {
	static bool is_stale;
	using data_type = int;
};

struct Timestep {
	static bool is_stale;
	using data_type = double;
};

/*! Determines how often cell is solved during time substepping.

N == solved every Nth substep
*/
struct Substepping_Period {
	static bool is_stale;
	using data_type = int;
};

//! Minimum substep period, solved no more often than every Nth substep
struct Substep_Min {
	static bool is_stale;
	using data_type = int;
	static const std::string get_name() { return {"Minimum substepping period"}; }
};

//! Maximum substep period, solved no less often than every Nth substep
struct Substep_Max {
	static bool is_stale;
	using data_type = int;
	static const std::string get_name() { return {"Maximum substepping period"}; }
};


} // namespace

#endif // ifndef PAMHD_VARIABLES_HPP
