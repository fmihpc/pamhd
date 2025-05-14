/*
face2r test.

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
#include "cstdlib"
#include "iostream"
#include "limits"
#include "random"

#include "math/interpolation.hpp"


double function(const std::array<double, 3>& r)
{
	return -std::sin(r[0] / 2) * std::cos(2 * r[1]) / (r[2] + 10);
}

int main()
{
	using std::abs;
	using std::array;
	using std::max;

	double
		old_norm = std::numeric_limits<double>::max(),
		old_dx = 0;
	for (size_t nr_of_cells = 1; nr_of_cells <= 8; nr_of_cells *= 2) {
		const array<double, 3> length{
			3.0 / nr_of_cells,
			1.5 / nr_of_cells,
			4.0 / nr_of_cells
		};

		const array<double, 3>
			start{-length[0]/2, -length[1]/2, -length[2]/2},
			end{+length[0]/2, +length[1]/2, +length[2]/2};
		pamhd::Face_Type<double> data;
		data(-3) = function({0, 0, start[2]});
		data(-2) = function({0, start[1], 0});
		data(-1) = function({start[0], 0, 0});
		data(+1) = function({end[0], 0, 0});
		data(+2) = function({0, end[1], 0});
		data(+3) = function({0, 0, end[2]});

		std::uniform_real_distribution<>
			xrange(start[0], end[0]),
			yrange(start[1], end[1]),
			zrange(start[2], end[2]);
		constexpr size_t iters = 10;
		std::mt19937 rnd;
		double norm = 0;
		for (size_t i = 0; i < iters; i++) {
			const array<double, 3> r{
				xrange(rnd), yrange(rnd), zrange(rnd)
			};
			const auto result = pamhd::math::face2r(r, end, length, data);
			norm = max(norm, abs(result[0] - function(r)));
			norm = max(norm, abs(result[1] - function(r)));
			norm = max(norm, abs(result[2] - function(r)));
		}

		if (norm > old_norm) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< ": Norm with " << length[1]
				<< " dx " << norm
				<< " is larger than with " << old_dx
				<< " dx " << old_norm
				<< std::endl;
			abort();
		}

		if (old_dx > 0) {
			const double order_of_accuracy
				= -log(norm / old_norm) / log(old_dx / length[1]);

			if (order_of_accuracy < 0.7) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< ": Order of accuracy from "
					<< old_dx << " to " << length[1]
					<< " is too low: " << order_of_accuracy
					<< std::endl;
				abort();
			}
		}

		old_dx = length[1];
		old_norm = norm;
	}

	return EXIT_SUCCESS;
}
