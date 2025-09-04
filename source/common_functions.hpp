/*
Common functions of PAMHD.

Copyright 2025 Finnish Meteorological Institute
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


Author(s): Ilja Honkonen
*/

#ifndef PAMHD_COMMON_FUNCTIONS_HPP
#define PAMHD_COMMON_FUNCTIONS_HPP


#include "cmath"
#include "numeric"
#include "optional"
#include "stdexcept"

#include "common_variables.hpp"


namespace pamhd {

// TODO: en.cppreference.com/w/cpp/numeric/linalg

auto add(const auto& in1, const auto& in2) {
	// vector and vector
	if constexpr (requires {
		in1[0] + in2[0];
		in1.size(); in2.size();
	}) {
		if (in1.size() != in2.size()) {
			throw std::range_error("Operands of different size().");
		}
		auto out = in1;
		for (size_t i = 0; i < size_t(in1.size()); i++) {
			out[i] += in2[i];
		}
		return out;
	// vector and scalar
	} else if constexpr (requires {in1[0] + in2;}) {
		auto out = in1;
		for (auto& o: out) o += in2;
		return out;
	// scalar and vector
	} else if constexpr (requires {in1 + in2[0];}) {
		auto out = in2;
		for (auto& o: out) o += in1;
		return out;
	// scalar and scalar
	} else if constexpr (requires {in1 + in2;}) {
		return in1 + in2;
	} else {
		static_assert(false, "Unsupported operand(s).");
	}
}

auto mul(const auto& in1, const auto& in2) {
	// vector and vector
	if constexpr (requires {
		in1[0] * in2[0];
		in1.size(); in2.size();
	}) {
		if (in1.size() != in2.size()) {
			throw std::range_error("Operands of different size().");
		}
		auto out = in1;
		for (size_t i = 0; i < size_t(in1.size()); i++) {
			out[i] *= in2[i];
		}
		return out;
	// vector and scalar
	} else if constexpr (requires {in1[0] * in2;}) {
		auto out = in1;
		for (auto& o: out) o *= in2;
		return out;
	// scalar and vector
	} else if constexpr (requires {in1 * in2[0];}) {
		auto out = in2;
		for (auto& o: out) o *= in1;
		return out;
	// scalar and scalar
	} else if constexpr (requires {in1 * in2;}) {
		return in1 * in2;
	} else {
		static_assert(false, "Unsupported operand(s).");
	}
}

auto cube(const auto& in) {
	if constexpr (requires {in[0];}) {
		auto out = in[0]; out = 0;
		for (const auto& i: in) out += i*i*i;
		return out;
	} else {
		return in*in*in;
	}
}

auto sum(const auto& in) {
	if constexpr (requires {in[0];}) {
		auto out = in[0]; out = 0;
		for (const auto& i: in) out += i;
		return out;
	} else {
		return in;
	}
}

auto dot(const auto& in1, const auto& in2) {
	// scalar and scalar
	if constexpr (requires {in1 * in2;}) {
		return in1 * in2;
	// rest supported by mul
	} else if constexpr (requires {mul(in1, in2);}) {
		const auto out = mul(in1, in2);
		return sum(out);
	} else {
		static_assert(false, "Unsupported operand(s).");
	}
}

auto norm2(const auto& in) {
	return dot(in, in);
}

auto norm(const auto& in) {
	using std::sqrt;
	return sqrt(norm2(in));
}

auto cross(const auto& in1, const auto& in2) {
	if constexpr (requires {
		in1[0]*in1[2] + in2[0]*in2[2];
	}) {
		auto out = in1; out = {
			in1[1]*in2[2] - in1[2]*in2[1],
			in1[2]*in2[0] - in1[0]*in2[2],
			in1[0]*in2[1] - in1[1]*in2[0]
		};
		return out;
	} else {
		static_assert(false, "Unsupported operand(s).");
	}
}

//! negates (all components of) in
auto neg(const auto& in) {
	if constexpr (requires {in[0];}) {
		auto out = in;
		for (auto& o: out) o *= -1;
		return out;
	} else {
		return -in;
	}
}

//! Makes sure all cells agree on common face Bs
template <
	class Grid,
	class Face_Magnetic_Field_Getter,
	class Volume_Magnetic_Field_Getter,
	class Mag_Error_Getter,
	class Cell_Type_Getter
> void sync_magnetic_field(
	Grid& grid,
	const Face_Magnetic_Field_Getter& Face_B,
	const Volume_Magnetic_Field_Getter& Vol_B,
	const Mag_Error_Getter& B_Error,
	const Cell_Type_Getter& CType
) try {
	using Cell = Grid::cell_data_type;

	bool update_copies = false;
	if (Face_B.type().is_stale) {
		update_copies = true;
		Cell::set_transfer_all(true, Face_B.type());
		Face_B.type().is_stale = false;
	}
	if (CType.type().is_stale) {
		update_copies = true;
		Cell::set_transfer_all(true, CType.type());
		CType.type().is_stale = false;
	}
	if (update_copies) {
		grid.update_copies_of_remote_neighbors();
	}
	Cell::set_transfer_all(false, Face_B.type(), CType.type());

	for (const auto& cell: grid.local_cells()) {
		if (CType.data(*cell.data) < 0) continue;

		// track error between cell and smaller neighbors
		pamhd::Face_Type<std::optional<double>> oldB;
		// prepare face B on side(s) of smaller neighbors
		for (const auto& neighbor: cell.neighbors_of) {
			const auto& fn = neighbor.face_neighbor;
			if (
				fn != 0
				and CType.data(*neighbor.data) >= 0
				and neighbor.relative_size > 0
			) {
				if (oldB(fn)) continue;
				oldB(fn) = Face_B.data(*cell.data)(fn);
				Face_B.data(*cell.data)(fn) = 0;
			}
		}

		for (const auto& neighbor: cell.neighbors_of) {
			const auto& fn = neighbor.face_neighbor;
			if (
				fn == 0
				or CType.data(*neighbor.data) < 0
			) continue;

			if (neighbor.relative_size == 0) {
				const auto diff
					= Face_B.data(*cell.data)(fn)
					- Face_B.data(*neighbor.data)(-fn);
				B_Error.data(*cell.data) += std::abs(diff);
				Face_B.data(*cell.data)(fn) -= diff / 2;
				Face_B.data(*neighbor.data)(-fn) += diff / 2;
			} else if (neighbor.relative_size > 0) {
				// replace with average of smaller neighbors
				Face_B.data(*cell.data)(fn) += Face_B.data(*neighbor.data)(-fn) / 4;
			}
		}

		for (int dir: {-3,-2,-1,+1,+2,+3}) {
			if (not oldB(dir)) continue;
			oldB(dir).value() -= Face_B.data(*cell.data)(dir);
			B_Error.data(*cell.data) += std::abs(oldB(dir).value());
		}
	}

	Face_B.type().is_stale = true;
	Vol_B.type().is_stale = true;

} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


} // namespace

#endif // ifndef PAMHD_COMMON_FUNCTIONS_HPP
