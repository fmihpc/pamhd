/*
Interpolation functions of PAMHD.

Copyright 2015, 2016, 2017 Ilja Honkonen
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

#ifndef PAMHD_INTERPOLATE_HPP
#define PAMHD_INTERPOLATE_HPP


#include "array"

#include "common_functions.hpp"


namespace pamhd {

namespace detail {

template<class T> T ZERO()
{
	return 0;
}

} // namespace pamhd


/*!
Returns value linearly interpolated from data to position coord.

data is assumed to be in the order:
(-x, -y, -z), (0, -y, -z), (+x, -y, -z),
(-x, 0, -z), ... , (0, +y, +z), (+x, +y, +z)

start marks the coordinate where
data at (-x, -y, -z) is located, end
marks location of (+x, +y, +z).

coord must be within start and end.
*/
template<class Coord_T, class Data_T> Data_T interpolate(
	const Coord_T& coord,
	const Coord_T& start,
	const Coord_T& end,
	const std::array<Data_T, 27>& data
) {
	using std::fabs;
	using std::max;

	const Coord_T
		mid{
			0.5 * (end[0] + start[0]),
			0.5 * (end[1] + start[1]),
			0.5 * (end[2] + start[2])
		},
		// length of one cell
		dr{
			0.5 * (end[0] - start[0]),
			0.5 * (end[1] - start[1]),
			0.5 * (end[2] - start[2])
		};

	// weights for corresponding item in data
	const std::array<double, 27> weights{{
			// data point at -x, -y, -z
			max(0.0, 1 - (coord[0] - start[0]) / dr[0])
			* max(0.0, 1 - (coord[1] - start[1]) / dr[1])
			* max(0.0, 1 - (coord[2] - start[2]) / dr[2]),
			// 0, -y, -z
			(1 - fabs(coord[0] - mid[0]) / dr[0])
			* max(0.0, 1 - (coord[1] - start[1]) / dr[1])
			* max(0.0, 1 - (coord[2] - start[2]) / dr[2]),
			// +x, -y, -z
			max(0.0, 1 - (end[0] - coord[0]) / dr[0])
			* max(0.0, 1 - (coord[1] - start[1]) / dr[1])
			* max(0.0, 1 - (coord[2] - start[2]) / dr[2]),
			// -x, 0, -z
			max(0.0, 1 - (coord[0] - start[0]) / dr[0])
			* (1 - fabs(coord[1] - mid[1]) / dr[1])
			* max(0.0, 1 - (coord[2] - start[2]) / dr[2]),
			// 0, 0, -z
			(1 - fabs(coord[0] - mid[0]) / dr[0])
			* (1 - fabs(coord[1] - mid[1]) / dr[1])
			* max(0.0, 1 - (coord[2] - start[2]) / dr[2]),
			// +x, 0, -z
			max(0.0, 1 - (end[0] - coord[0]) / dr[0])
			* (1 - fabs(coord[1] - mid[1]) / dr[1])
			* max(0.0, 1 - (coord[2] - start[2]) / dr[2]),
			// -x, +y, -z
			max(0.0, 1 - (coord[0] - start[0]) / dr[0])
			* max(0.0, 1 - (end[1] - coord[1]) / dr[1])
			* max(0.0, 1 - (coord[2] - start[2]) / dr[2]),
			// 0, +y, -z
			(1 - fabs(coord[0] - mid[0]) / dr[0])
			* max(0.0, 1 - (end[1] - coord[1]) / dr[1])
			* max(0.0, 1 - (coord[2] - start[2]) / dr[2]),
			// +x, +y, -z
			max(0.0, 1 - (end[0] - coord[0]) / dr[0])
			* max(0.0, 1 - (end[1] - coord[1]) / dr[1])
			* max(0.0, 1 - (coord[2] - start[2]) / dr[2]),
			// -x, -y, 0
			max(0.0, 1 - (coord[0] - start[0]) / dr[0])
			* max(0.0, 1 - (coord[1] - start[1]) / dr[1])
			* (1 - fabs(coord[2] - mid[2]) / dr[2]),
			// 0, -y, 0
			(1 - fabs(coord[0] - mid[0]) / dr[0])
			* max(0.0, 1 - (coord[1] - start[1]) / dr[1])
			* (1 - fabs(coord[2] - mid[2]) / dr[2]),
			// +x, -y, 0
			max(0.0, 1 - (end[0] - coord[0]) / dr[0])
			* max(0.0, 1 - (coord[1] - start[1]) / dr[1])
			* (1 - fabs(coord[2] - mid[2]) / dr[2]),
			// -x, 0, 0
			max(0.0, 1 - (coord[0] - start[0]) / dr[0])
			* (1 - fabs(coord[1] - mid[1]) / dr[1])
			* (1 - fabs(coord[2] - mid[2]) / dr[2]),
			// 0, 0, 0
			(1 - fabs(coord[0] - mid[0]) / dr[0])
			* (1 - fabs(coord[1] - mid[1]) / dr[1])
			* (1 - fabs(coord[2] - mid[2]) / dr[2]),
			// +x, 0, 0
			max(0.0, 1 - (end[0] - coord[0]) / dr[0])
			* (1 - fabs(coord[1] - mid[1]) / dr[1])
			* (1 - fabs(coord[2] - mid[2]) / dr[2]),
			// -x, +y, 0
			max(0.0, 1 - (coord[0] - start[0]) / dr[0])
			* max(0.0, 1 - (end[1] - coord[1]) / dr[1])
			* (1 - fabs(coord[2] - mid[2]) / dr[2]),
			// 0, +y, 0
			(1 - fabs(coord[0] - mid[0]) / dr[0])
			* max(0.0, 1 - (end[1] - coord[1]) / dr[1])
			* (1 - fabs(coord[2] - mid[2]) / dr[2]),
			// +x, +y, 0
			max(0.0, 1 - (end[0] - coord[0]) / dr[0])
			* max(0.0, 1 - (end[1] - coord[1]) / dr[1])
			* (1 - fabs(coord[2] - mid[2]) / dr[2]),
			// -x, -y, +z
			max(0.0, 1 - (coord[0] - start[0]) / dr[0])
			* max(0.0, 1 - (coord[1] - start[1]) / dr[1])
			* max(0.0, 1 - (end[2] - coord[2]) / dr[2]),
			// 0, -y, +z
			(1 - fabs(coord[0] - mid[0]) / dr[0])
			* max(0.0, 1 - (coord[1] - start[1]) / dr[1])
			* max(0.0, 1 - (end[2] - coord[2]) / dr[2]),
			// +x, -y, +z
			max(0.0, 1 - (end[0] - coord[0]) / dr[0])
			* max(0.0, 1 - (coord[1] - start[1]) / dr[1])
			* max(0.0, 1 - (end[2] - coord[2]) / dr[2]),
			// -x, 0, +z
			max(0.0, 1 - (coord[0] - start[0]) / dr[0])
			* (1 - fabs(coord[1] - mid[1]) / dr[1])
			* max(0.0, 1 - (end[2] - coord[2]) / dr[2]),
			// 0, 0, +z
			(1 - fabs(coord[0] - mid[0]) / dr[0])
			* (1 - fabs(coord[1] - mid[1]) / dr[1])
			* max(0.0, 1 - (end[2] - coord[2]) / dr[2]),
			// +x, 0, +z
			max(0.0, 1 - (end[0] - coord[0]) / dr[0])
			* (1 - fabs(coord[1] - mid[1]) / dr[1])
			* max(0.0, 1 - (end[2] - coord[2]) / dr[2]),
			// -x, +y, +z
			max(0.0, 1 - (coord[0] - start[0]) / dr[0])
			* max(0.0, 1 - (end[1] - coord[1]) / dr[1])
			* max(0.0, 1 - (end[2] - coord[2]) / dr[2]),
			// 0, +y, +z
			(1 - fabs(coord[0] - mid[0]) / dr[0])
			* max(0.0, 1 - (end[1] - coord[1]) / dr[1])
			* max(0.0, 1 - (end[2] - coord[2]) / dr[2]),
			// +x, +y, +z
			max(0.0, 1 - (end[0] - coord[0]) / dr[0])
			* max(0.0, 1 - (end[1] - coord[1]) / dr[1])
			* max(0.0, 1 - (end[2] - coord[2]) / dr[2])
	}};

	Data_T ret_val = {0, 0, 0};
	for (size_t i = 0; i < data.size(); i++) {
		ret_val = pamhd::add(ret_val, pamhd::mul(weights[i], data[i]));
	}
	return ret_val;
}

} // namespace pamhd

#endif
