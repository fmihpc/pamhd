/*
Tests for interpolate function of PAMHD with Eigen vector type.

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


#include "array"
#include "cmath"
#include "cstdlib"
#include "iostream"

#include "Eigen/Core"
#include "interpolate.hpp"


using namespace std;
using namespace pamhd;

Eigen::Vector3d f(const double x, const double y, const double z)
{
	return {
		-0.5*(x + 3) + 2*(y - 1) + 4*z + 5,
		2*(x + 3) - 3*(y - 1) + 0.25*z + 3,
		1.5*(x + 3) + 2*(y - 1) - z + 2
	};
}

void test(
	const Eigen::Vector3d& start,
	const Eigen::Vector3d& end
) {
	const Eigen::Vector3d mid{
		0.5 * (end[0] + start[0]),
		0.5 * (end[1] + start[1]),
		0.5 * (end[2] + start[2])
	};

	const std::array<Eigen::Vector3d, 27> data{{
		f(start[0], start[1], start[2]),
		f(mid[0], start[1], start[2]),
		f(end[0], start[1], start[2]),
		f(start[0], mid[1], start[2]),
		f(mid[0], mid[1], start[2]),
		f(end[0], mid[1], start[2]),
		f(start[0], end[1], start[2]),
		f(mid[0], end[1], start[2]),
		f(end[0], end[1], start[2]),
		f(start[0], start[1], mid[2]),
		f(mid[0], start[1], mid[2]),
		f(end[0], start[1], mid[2]),
		f(start[0], mid[1], mid[2]),
		f(mid[0], mid[1], mid[2]),
		f(end[0], mid[1], mid[2]),
		f(start[0], end[1], mid[2]),
		f(mid[0], end[1], mid[2]),
		f(end[0], end[1], mid[2]),
		f(start[0], start[1], end[2]),
		f(mid[0], start[1], end[2]),
		f(end[0], start[1], end[2]),
		f(start[0], mid[1], end[2]),
		f(mid[0], mid[1], end[2]),
		f(end[0], mid[1], end[2]),
		f(start[0], end[1], end[2]),
		f(mid[0], end[1], end[2]),
		f(end[0], end[1], end[2]),
	}};

	for (size_t samples = 1; samples <= 5; samples++)
	for (double x = start[0]; x <= end[0]; x += (end[0] - start[0])/samples)
	for (double y = start[1]; y <= end[1]; y += (end[1] - start[1])/samples)
	for (double z = start[2]; z <= end[2]; z += (end[2] - start[2])/samples) {
		const Eigen::Vector3d
			coord{x, y, z},
			interpolated = interpolate(coord, start, end, data),
			exact = f(coord[0], coord[1], coord[2]);
		if ((interpolated - exact).norm() > 1e-10) {
			abort();
		}
	}
}

int main()
{
	test({1., -4., -2.}, {5., 0., 2.});
	test({-1.3, 4.5, 2.1}, {5.5, 6.7, 10.8});

	return EXIT_SUCCESS;
}
