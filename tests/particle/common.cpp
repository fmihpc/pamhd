/*
Tests for common particle functions of PAMHD.

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
#include "iostream"
#include "random"

#include "particle/common.hpp"
#include "particle/variables.hpp"


using namespace std;
using namespace pamhd::particle;


int main()
{
	using std::fabs;
	using std::min;
	using std::pow;
	using std::sqrt;

	// short hand notation for variables
	const Position Pos{};
	const Velocity Vel{};
	const Mass Mas{};
	const Charge_Mass_Ratio C2M{};

	Particle_T<> particle;
	std::array<double, 3> electric_field, magnetic_field;

	// not needed
	particle[Mas] = 0;
	particle[Pos] = {0, 0, 0};

	// test get_gyro_info
	particle[Vel] = {6, 2, 3};
	particle[C2M] = -2;
	electric_field = {0, 0, 0};
	magnetic_field = {2, 3, 6};
	const auto B_mag = 7.0, v_perp_mag = sqrt(54145.0) / 49.0, c2m_abs = 2.0;

	const auto gyro_info = get_gyro_info(particle[C2M], particle[Vel], magnetic_field);
	if (fabs(gyro_info.first - v_perp_mag / B_mag / c2m_abs) > 1e-20) {
		std::cerr <<  __FILE__ << ":" << __LINE__
			<< " Incorrect gyro radius: "
			<< gyro_info.first << ", should be: " << v_perp_mag / B_mag / c2m_abs
			<< std::endl;
		abort();
	}
	if (fabs(gyro_info.second - c2m_abs * B_mag / 2 / M_PI) > 1e-20) {
		std::cerr <<  __FILE__ << ":" << __LINE__
			<< " Incorrect gyro frequency: "
			<< gyro_info.second << ", should be: " << c2m_abs * B_mag / 2 / M_PI
			<< std::endl;
		abort();
	}


	// test get_max_step
	double cell_length = 0.1;
	auto step_size
		= get_step_size(
			0.05,
			cell_length,
			particle[C2M],
			particle[Vel],
			electric_field,
			magnetic_field
		);

	if (min(step_size.first, step_size.second) > 1.001 * cell_length / 7.0) {
		std::cerr <<  __FILE__ << ":" << __LINE__
			<< " Incorrect maximum time step for particle: " << min(step_size.first, step_size.second)
			<< std::endl;
		abort();
	}

	cell_length = 10;
	step_size = get_step_size(
		0.05,
		cell_length,
		particle[C2M],
		particle[Vel],
		electric_field,
		magnetic_field
	);

	if (min(step_size.first, step_size.second) > 1.001 * 0.05 / c2m_abs / B_mag * 2 * M_PI) {
		std::cerr <<  __FILE__ << ":" << __LINE__
			<< " Incorrect maximum time step for particle: " << min(step_size.first, step_size.second)
			<< std::endl;
		abort();
	}

	return EXIT_SUCCESS;
}
