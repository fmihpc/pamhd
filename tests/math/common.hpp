/*
Common definitions for tests in this directory.

Copyright 2023 Finnish Meteorological Institute
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

#ifndef TESTS_MATH_COMMON_HPP
#define TESTS_MATH_COMMON_HPP

#include "array"

#include "gensimcell.hpp"

//! each component is normal to and at corresponding + dir face
struct Vector_Pos {
	using data_type = std::array<double, 3>;
};

//! component at - dir face
struct Vector_Neg {
	using data_type = std::array<double, 3>;
};

struct Divergence {
	using data_type = double;
};

struct Is_Primary_Face {
	using data_type = std::array<bool, 6>;
	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype() const {
		return std::make_tuple(nullptr, 0, MPI_BYTE);
	}
};

struct Type {
	using data_type = int;
};

using Cell = gensimcell::Cell<
	gensimcell::Always_Transfer,
	Vector_Pos,
	Vector_Neg,
	Divergence,
	Is_Primary_Face,
	Type
>;

#endif
