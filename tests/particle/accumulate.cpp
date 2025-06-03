/*
Serial test for particle accumulator of PAMHD.

Copyright 2014, 2015, 2016, 2017 Ilja Honkonen
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
#include "exception"
#include "iostream"
#include "utility"
#include "vector"

#include "particle/accumulate.hpp"


using namespace std;
using namespace pamhd::particle;


constexpr double get_grid_start() { return -3; }

constexpr double get_grid_end() { return 6; }

constexpr double get_grid_length() { return get_grid_end() - get_grid_start(); }

constexpr double get_cell_length(const size_t grid_size)
{
	return (get_grid_end() - get_grid_start()) / grid_size;
}

constexpr double get_cell_min(const size_t index, const size_t grid_size)
{
	return index * get_cell_length(grid_size) + get_grid_start();
}

constexpr double get_cell_max(const size_t index, const size_t grid_size)
{
	return (index + 1) * get_cell_length(grid_size) + get_grid_start();
}

constexpr double get_cell_center(const size_t index, const size_t grid_size)
{
	return (get_cell_max(index, grid_size) + get_cell_min(index, grid_size)) / 2;
}


double f(const double x)
{
	return 1.5 * sin(2 * M_PI * (x - get_grid_start()) / get_grid_length());
}


int main()
{
	using std::array;

	double value = 0;
	array<double, 3>
		value_min{0, 0, 0},
		value_max{0, 0, 0},
		cell_min{0, 0, 0},
		cell_max{0, 0, 0};

	// test error detection
	try {
		value_min = {0, 0, 1};
		value_max = {1, 1, 0};
		cell_min = {0, 0, 0};
		cell_max = {1, 1, 1};
		get_accumulated_value(
			value,
			value_min,
			value_max,
			cell_min,
			cell_max
		);

		std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
			<< "Invalid value volume not caught."
			<< std::endl;
		abort();
	} catch(const std::out_of_range& e) {}

	try {
		value_min = {0, 0, 0};
		value_max = {1, 1, 1};
		cell_min = {0, 1, 0};
		cell_max = {1, 0, 1};
		get_accumulated_value(
			value,
			value_min,
			value_max,
			cell_min,
			cell_max
		);

		std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
			<< "Invalid cell volume not caught."
			<< std::endl;
		abort();
	} catch(const std::out_of_range& e) {}

	// test regular accumulation
	const auto accu1
		= get_accumulated_value(
			1.0,
			array<double, 3>{-1, -1, -1},
			array<double, 3>{ 1,  1,  1},
			array<double, 3>{ 0,  0,  0},
			array<double, 3>{ 2,  2,  2}
		);
	if (fabs(accu1 - 1.0/8.0) > 1e-10) {
		std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
			<< "Difference for accumulated value too large from correct."
			<< std::endl;
		abort();
	}

	// test weighted accumulation
	const auto
		accu21 = get_accumulated_value_weighted(
			-3.0,
			2.0,
			array<double, 3>{-1, -1, -1},
			array<double, 3>{ 1,  1,  1},
			array<double, 3>{ 0,  0,  0},
			array<double, 3>{ 2,  2,  2}
		),
		accu22 = get_accumulated_value_weighted(
			10.0,
			1.0,
			array<double, 3>{ 0,  1,  1},
			array<double, 3>{ 2,  3,  3},
			array<double, 3>{ 0,  0,  0},
			array<double, 3>{ 2,  2,  2}
		);
	if (fabs(accu21.first + accu22.first - 14) > 1e-10) {
		std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
			<< "Difference for accumulated value too large from correct."
			<< std::endl;
		abort();
	}
	if (fabs(accu21.second + accu22.second - 4) > 1e-10) {
		std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
			<< "Difference for total weights too large from correct."
			<< std::endl;
		abort();
	}

	// test vector accumulation
	const auto vec_accu
		= get_accumulated_value(
			array<double, 3>{1, 2, 3},
			array<double, 3>{0, 0, 0},
			array<double, 3>{1, 1, 1},
			array<double, 3>{0, 0, 0},
			array<double, 3>{1, 1, 1}
		);
	if (
		fabs(vec_accu[0] - 1) > 1e-10
		or fabs(vec_accu[1] - 2) > 1e-10
		or fabs(vec_accu[2] - 3) > 1e-10
	) {
		std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
			<< "Difference for accumulated vector too large from correct."
			<< std::endl;
		abort();
	}

	const auto vec_accu3
		= get_accumulated_value(
			array<double, 3>{20, 30},
			array<double, 3>{0, 0, 0},
			array<double, 3>{1, 1, 1},
			array<double, 3>{0, 0, 0},
			array<double, 3>{1, 1, 1}
		);
	if (vec_accu3.size() != 2) {
		std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
			<< "Wrong number of dimensions in accumulated array."
			<< std::endl;
		abort();
	}
	if (
		fabs(vec_accu3[0] - 20) > 1e-10
		or fabs(vec_accu3[1] - 30) > 1e-10
	) {
		std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
			<< "Difference for accumulated vector too large from correct."
			<< std::endl;
		abort();
	}

	// test accumulation of negative values
	const auto vec_accu2
		= get_accumulated_value(
			array<double, 3>{-1, -2, -3},
			array<double, 3>{0, 0, 0},
			array<double, 3>{1, 1, 1},
			array<double, 3>{0, 0, 0},
			array<double, 3>{1, 1, 1}
		);
	if (
		fabs(vec_accu2[0] + 1) > 1e-10
		or fabs(vec_accu2[1] + 2) > 1e-10
		or fabs(vec_accu2[2] + 3) > 1e-10
	) {
		std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
			<< "Difference for accumulated vector too large from correct."
			<< std::endl;
		abort();
	}

	constexpr size_t
		nr_of_values = 10000,
		max_nr_of_cells = 640;

	// test convergence of accumulator
	double old_norm = std::numeric_limits<double>::max();
	size_t old_nr_of_cells = 0;
	for (size_t nr_of_cells = 10; nr_of_cells <= max_nr_of_cells; nr_of_cells *= 2) {

		const double
			values_per_cell = double(nr_of_values) / nr_of_cells,
			cell_length = get_cell_length(nr_of_cells);

		/*
		Create evenly spaced values inside of grid volume.
		values[i].first = value, values[i].second = position
		*/
		std::vector<std::pair<double, double>> values;
		for (size_t i = 0; i < nr_of_values; i++) {
			const auto position
				= get_grid_start()
				+ (i + 0.5) * get_grid_length() / nr_of_values;

			values.push_back(std::make_pair(f(position) / values_per_cell, position));
		}

		std::vector<double> accumulated(nr_of_cells, 0);

		for (const auto& item: values) {
			for (size_t cell_i = 0; cell_i < nr_of_cells; cell_i++) {
				const auto accumulated_value = get_accumulated_value(
					item.first,
					array<double, 3>{item.second - cell_length / 2, 0, 0},
					array<double, 3>{item.second + cell_length / 2, 1, 1},
					array<double, 3>{get_cell_min(cell_i, nr_of_cells), 0, 0},
					array<double, 3>{get_cell_max(cell_i, nr_of_cells), 1, 1}
				);

				accumulated[cell_i] += accumulated_value;
			}

			// accumulate across periodic boundary
			accumulated[0] += get_accumulated_value(
				item.first,
				array<double, 3>{item.second - cell_length / 2, 0, 0},
				array<double, 3>{item.second + cell_length / 2, 1, 1},
				array<double, 3>{get_cell_max(nr_of_cells - 1, nr_of_cells), 0, 0},
				array<double, 3>{get_cell_max(nr_of_cells - 1, nr_of_cells) + cell_length, 1, 1}
			);
			accumulated[nr_of_cells - 1] += get_accumulated_value(
				item.first,
				array<double, 3>{item.second - cell_length / 2, 0, 0},
				array<double, 3>{item.second + cell_length / 2, 1, 1},
				array<double, 3>{get_cell_min(0, nr_of_cells) - cell_length, 0, 0},
				array<double, 3>{get_cell_min(0, nr_of_cells), 1, 1}
			);
		}

		// get inf norm from analytic
		double norm = 0;
		for (size_t cell_i = 0; cell_i < nr_of_cells; cell_i++) {
			const double cell_center = get_cell_center(cell_i, nr_of_cells);
			norm = max(norm, fabs(accumulated[cell_i] - f(cell_center)));
		}
		norm /= nr_of_cells;

		if (old_nr_of_cells > 0 and old_nr_of_cells < 320) {
			const double order_of_accuracy
				= -log(norm / old_norm)
				/ log(double(nr_of_cells) / old_nr_of_cells);

			if (order_of_accuracy < 0.9) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< ": Order of accuracy from "
					<< old_nr_of_cells << " to " << nr_of_cells
					<< " is too low: " << order_of_accuracy
					<< std::endl;
				abort();
			}
		}

		old_norm = norm;
		old_nr_of_cells = nr_of_cells;
	}

	return EXIT_SUCCESS;
}
