/*
Tests for particle propagator of PAMHD.

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

#include "prettyprint.hpp"

#include "particle/variables.hpp"
#include "particle/solve.hpp"


using namespace std;
using namespace pamhd::particle;


int main()
{
	using std::fabs;
	using std::pow;

	// short hand notation for variables
	const Position Pos{};
	const Velocity Vel{};
	const Mass Mas{};
	const Charge_Mass_Ratio C2M{};

	double time, dt, gyro_radius, gyro_period;
	Particle_T<> particle;
	std::array<double, 3> electric_field, magnetic_field;

	// not needed
	particle[Mas] = 0;

	// only E
	dt = 1;
	particle[Pos] = {0, 3, 6};
	particle[Vel] = {-3, 2, -1};
	particle[C2M] = 2;
	electric_field = {1, -2, 10.5};
	magnetic_field = {0, 0, 0};

	time = 0;
	for (size_t step = 0; step < 5; step++) {
		std::tie(
			particle[Pos],
			particle[Vel]
		) = propagate(
			particle[Pos],
			particle[Vel],
			electric_field,
			magnetic_field,
			particle[C2M],
			dt
		);
		time += dt;
	}

	if (
		particle[Pos][0] != 5
		or particle[Pos][1] != -27
		or particle[Pos][2] != 211
	) {
		std::cerr <<  __FILE__ << ":" << __LINE__
			<< " Particle1 at incorrect location: "
			<< particle[Pos]
			<< std::endl;
		abort();
	}

	if (
		particle[Vel][0] != 7
		or particle[Vel][1] != -18
		or particle[Vel][2] != 104
	) {
		std::cerr <<  __FILE__ << ":" << __LINE__
			<< " Particle1 has incorrect velocity: "
			<< particle[Vel]
			<< std::endl;
		abort();
	}


	// only B, rotate 90 degrees around y axis
	dt = 5 * M_PI / 1000;
	particle[Pos] = {0, 0, 0};
	particle[Vel] = {1, 0, 0};
	particle[C2M] = 1;
	electric_field = {0, 0, 0};
	magnetic_field = {0, 0.1, 0};
	gyro_radius = particle[Vel][0] / magnetic_field[1] / particle[C2M];
	gyro_period = 2 * M_PI / particle[C2M] / magnetic_field[1];

	for (size_t step = 0; step < 1000; step++) {
		std::tie(
			particle[Pos],
			particle[Vel]
		) = propagate(
			particle[Pos],
			particle[Vel],
			electric_field,
			magnetic_field,
			particle[C2M],
			dt
		);
	}

	if (
		fabs(particle[Pos][0] - gyro_radius) > 0.01
		or particle[Pos][1] != 0
		or fabs(particle[Pos][2] - gyro_radius) > 0.01
	) {
		std::cerr <<  __FILE__ << ":" << __LINE__
			<< " Particle1 at incorrect location: "
			<< particle[Pos]
			<< ", should be: " << gyro_radius << ", 0, " << gyro_radius
			<< std::endl;
		abort();
	}

	if (
		fabs(particle[Vel][0]) > 1e-10
		or particle[Vel][1] != 0
		or fabs(particle[Vel][2] - 1) > 1e-10
	) {
		std::cerr <<  __FILE__ << ":" << __LINE__
			<< " Particle1 has incorrect velocity: "
			<< particle[Vel]
			<< std::endl;
		abort();
	}

	// rotate another 90 degrees
	for (size_t step = 0; step < 1000; step++) {
		std::tie(
			particle[Pos],
			particle[Vel]
		) = propagate(
			particle[Pos],
			particle[Vel],
			electric_field,
			magnetic_field,
			particle[C2M],
			dt
		);
	}

	if (
		fabs(particle[Pos][0]) > 0.02
		or particle[Pos][1] != 0
		or fabs(particle[Pos][2] - 2 * gyro_radius) > 1e-5
	) {
		std::cerr <<  __FILE__ << ":" << __LINE__
			<< " Particle1 at incorrect location: "
			<< particle[Pos]
			<< std::endl;
		abort();
	}

	if (
		fabs(particle[Vel][0] + 1) > 1e-10
		or particle[Vel][1] != 0
		or fabs(particle[Vel][2]) > 1e-10
	) {
		std::cerr <<  __FILE__ << ":" << __LINE__
			<< " Particle1 has incorrect velocity: "
			<< particle[Vel]
			<< std::endl;
		abort();
	}

	// rotate 180 degrees
	for (size_t step = 0; step < 2000; step++) {
		std::tie(
			particle[Pos],
			particle[Vel]
		) = propagate(
			particle[Pos],
			particle[Vel],
			electric_field,
			magnetic_field,
			particle[C2M],
			dt
		);
	}

	if (
		fabs(particle[Pos][0]) > 1e-10
		or particle[Pos][1] != 0
		or fabs(particle[Pos][2]) > 1e-10
	) {
		std::cerr <<  __FILE__ << ":" << __LINE__
			<< " Particle1 at incorrect location: "
			<< particle[Pos]
			<< std::endl;
		abort();
	}

	if (
		fabs(particle[Vel][0] - 1) > 1e-10
		or particle[Vel][1] != 0
		or fabs(particle[Vel][2]) > 1e-10
	) {
		std::cerr <<  __FILE__ << ":" << __LINE__
			<< " Particle1 has incorrect velocity: "
			<< particle[Vel]
			<< std::endl;
		abort();
	}


	// only B, rotate 360 degrees around x axis
	dt = 4 * M_PI / 1000;
	particle[Pos] = {0, 0, 0};
	particle[Vel] = {0, 1, 0};
	particle[C2M] = 1;
	electric_field = {0, 0, 0};
	magnetic_field = {0.5, 0, 0};

	for (size_t step = 0; step < 1000; step++) {
		std::tie(
			particle[Pos],
			particle[Vel]
		) = propagate(
			particle[Pos],
			particle[Vel],
			electric_field,
			magnetic_field,
			particle[C2M],
			dt
		);
	}

	if (
		particle[Pos][0] != 0
		or fabs(particle[Pos][1]) > 1e-10
		or fabs(particle[Pos][2]) > 1e-10
	) {
		std::cerr <<  __FILE__ << ":" << __LINE__
			<< " Particle1 at incorrect location: "
			<< particle[Pos]
			<< std::endl;
		abort();
	}

	if (
		particle[Vel][0] != 0
		or fabs(particle[Vel][1] - 1) > 1e-10
		or fabs(particle[Vel][2]) > 1e-10
	) {
		std::cerr <<  __FILE__ << ":" << __LINE__
			<< " Particle1 has incorrect velocity: "
			<< particle[Vel]
			<< std::endl;
		abort();
	}


	// only B, rotate 360 degrees around z axis
	dt = 1 * M_PI / 1000;
	particle[Pos] = {0, 0, 0};
	particle[Vel] = {1, 0, 0};
	particle[C2M] = 1;
	electric_field = {0, 0, 0};
	magnetic_field = {0, 0, -2};

	for (size_t step = 0; step < 1000; step++) {
		std::tie(
			particle[Pos],
			particle[Vel]
		) = propagate(
			particle[Pos],
			particle[Vel],
			electric_field,
			magnetic_field,
			particle[C2M],
			dt
		);
	}

	if (
		fabs(particle[Pos][0]) > 1e-10
		or fabs(particle[Pos][1]) > 1e-10
		or particle[Pos][2] != 0
	) {
		std::cerr <<  __FILE__ << ":" << __LINE__
			<< " Particle1 at incorrect location: "
			<< particle[Pos]
			<< std::endl;
		abort();
	}

	if (
		fabs(particle[Vel][0] - 1) > 1e-10
		or fabs(particle[Vel][1]) > 1e-10
		or particle[Vel][2] != 0
	) {
		std::cerr <<  __FILE__ << ":" << __LINE__
			<< " Particle1 has incorrect velocity: "
			<< particle[Vel]
			<< std::endl;
		abort();
	}


	// only B, rotate 180 around z axis, opposite charge
	dt = M_PI / 2000;
	particle[Pos] = {0, 0, 0};
	particle[Vel] = {-1, 0, 0};
	particle[C2M] = -1;
	electric_field = {0, 0, 0};
	magnetic_field = {0, 0, 1};

	for (size_t step = 0; step < 2000; step++) {
		std::tie(
			particle[Pos],
			particle[Vel]
		) = propagate(
			particle[Pos],
			particle[Vel],
			electric_field,
			magnetic_field,
			particle[C2M],
			dt
		);
	}

	if (
		fabs(particle[Pos][0]) > 0.002
		or fabs(particle[Pos][1] + 2) > 1e-5
		or particle[Pos][2] != 0
	) {
		std::cerr <<  __FILE__ << ":" << __LINE__
			<< " Particle1 at incorrect location: "
			<< particle[Pos]
			<< std::endl;
		abort();
	}

	if (
		fabs(particle[Vel][0] - 1) > 1e-10
		or fabs(particle[Vel][1]) > 1e-10
		or particle[Vel][2] != 0
	) {
		std::cerr <<  __FILE__ << ":" << __LINE__
			<< " Particle1 has incorrect velocity: "
			<< particle[Vel]
			<< std::endl;
		abort();
	}


	/*
	E.cross(B) drift, e.g. from https://en.wikipedia.org/wiki/Guiding_center
	*/
	particle[Pos] = {0, 0, 0};
	particle[Vel] = {0, 0, 1};
	particle[C2M] = 1;
	electric_field = {1, 0, 0};
	magnetic_field = {0, 3, 0};
	gyro_radius = particle[Vel][2] / magnetic_field[1] / particle[C2M];
	gyro_period = 2 * M_PI / particle[C2M] / magnetic_field[1];
	dt = gyro_period / 1000;
	const double vel_drift = electric_field[0] / magnetic_field[1];

	for (size_t step = 0; step < 1000; step++) {
		std::tie(
			particle[Pos],
			particle[Vel]
		) = propagate(
			particle[Pos],
			particle[Vel],
			electric_field,
			magnetic_field,
			particle[C2M],
			dt
		);
	}

	if (
		fabs(particle[Pos][0]) > 1e-10
		or particle[Pos][1] != 0
		or fabs(particle[Pos][2] - vel_drift * gyro_period) > 1e-5
	) {
		std::cerr <<  __FILE__ << ":" << __LINE__
			<< " Particle1 at incorrect location: "
			<< particle[Pos]
			<< std::endl;
		abort();
	}

	if (
		fabs(particle[Vel][0]) > 1e-10
		or particle[Vel][1] != 0
		or fabs(particle[Vel][2] - 1) > 1e-10
	) {
		std::cerr <<  __FILE__ << ":" << __LINE__
			<< " Particle1 has incorrect velocity: "
			<< particle[Vel]
			<< std::endl;
		abort();
	}


	// Polarization drift
	particle[Pos] = {0, 0, 0};
	particle[Vel] = {0, 0, 1};
	particle[C2M] = 2;
	electric_field = {0.5, 0, 0};
	magnetic_field = {0, 3, 0};
	gyro_radius = particle[Vel][2] / magnetic_field[1] / particle[C2M];
	gyro_period = fabs(2 * M_PI / particle[C2M] / magnetic_field[1]);
	dt = gyro_period / 1000;
	const auto E_change = pamhd::mul(1/100, electric_field);
	const double
		ExB_drift = electric_field[0] / magnetic_field[1],
		polarization_drift = E_change[0] / particle[C2M] / pow(magnetic_field[1], 2);

	for (size_t step = 0; step < 1000; step++) {
		std::tie(
			particle[Pos],
			particle[Vel]
		) = propagate(
			particle[Pos],
			particle[Vel],
			pamhd::add(electric_field, pamhd::mul(E_change, step * dt / gyro_period)),
			magnetic_field,
			particle[C2M],
			dt
		);
	}

	if (
		fabs(particle[Pos][0] - 0.000277776) > 1e-5
		or particle[Pos][1] != 0
		or fabs(particle[Pos][2] - ExB_drift * gyro_period) > 0.001
	) {
		std::cerr <<  __FILE__ << ":" << __LINE__
			<< " Particle1 at incorrect location: "
			<< particle[Pos]
			<< ", should be: " << polarization_drift * gyro_period
			<< ", 0, " << ExB_drift * gyro_period
			<< std::endl;
		abort();
	}

	if (
		fabs(particle[Vel][0]) > 1e-10
		or particle[Vel][1] != 0
		or fabs(particle[Vel][2] - 1) > 0.002
	) {
		std::cerr <<  __FILE__ << ":" << __LINE__
			<< " Particle1 has incorrect velocity: "
			<< particle[Vel]
			<< std::endl;
		abort();
	}


	return EXIT_SUCCESS;
}
