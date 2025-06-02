/*
Tests bulk pressure accumulation from particles in PAMHD.

Copyright 2017 Ilja Honkonen
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

#include "common_functions.hpp"
#include "particle/accumulate.hpp"
#include "particle/common.hpp"
#include "particle/variables.hpp"


using namespace std;
using namespace pamhd::particle;


bool test_accumulation(
	const size_t nr_particles,
	const std::array<double, 3> particle_min,
	const std::array<double, 3> particle_max,
	const std::array<double, 3> cell_min,
	const std::array<double, 3> cell_max,
	const double particle_temp_nrj_ratio,
	const std::array<double, 3> bulk_velocity,
	const std::array<double, 3> temperature,
	const double total_mass,
	const double species_mass,
	const size_t seed
) {
	using std::fabs;
	using std::min;
	using std::pow;
	using std::sqrt;

	std::mt19937 random_source;
	random_source.seed(seed);

	const auto particles
		= create_particles<
			Particle_Internal,
			Mass,
			Charge_Mass_Ratio,
			Position,
			Velocity,
			Particle_ID,
			Species_Mass
		>(
			bulk_velocity,
			particle_min,
			particle_max,
			temperature,
			nr_particles,
			12.3,
			total_mass,
			species_mass,
			particle_temp_nrj_ratio,
			random_source
		);

	const double
		volume
			= (particle_max[0] - particle_min[0])
			* (particle_max[1] - particle_min[1])
			* (particle_max[2] - particle_min[2]),
		intersection_volume = get_intersection_volume(particle_min, particle_max, cell_min, cell_max),
		volume_fraction = intersection_volume / volume;

	// number of particles
	const double ref_nr_real = volume_fraction * total_mass / species_mass;
	double accu_nr_real = 0;
	for (const auto& part: particles) {
		accu_nr_real
			+= get_accumulated_value(
				part[Mass()] / part[Species_Mass()],
				particle_min, particle_max,
				cell_min, cell_max
			);
	}
	if (abs(accu_nr_real - ref_nr_real) / abs(ref_nr_real) > 1e-2) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
			<< accu_nr_real << " " << ref_nr_real
			<< std::endl;
		return false;
	}

	// bulk mass
	const double ref_mass = volume_fraction * total_mass;
	double accu_mass = 0;
	for (const auto& part: particles) {
		accu_mass
			+= get_accumulated_value(
				part[Mass()],
				particle_min, particle_max,
				cell_min, cell_max
			);
	}
	if (abs(accu_mass - ref_mass) / abs(ref_mass) > 1e-2) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
			<< accu_mass << " " << ref_mass
			<< std::endl;
		return false;
	}

	// bulk velocity
	std::pair<std::array<double, 3>, double> accu_vel{{0, 0, 0}, 0};
	for (const auto& part: particles) {
		const auto temp = get_accumulated_value_weighted(
			part[Velocity()],
			part[Mass()] / part[Species_Mass()],
			particle_min, particle_max,
			cell_min, cell_max
		);
		accu_vel.first = pamhd::add(accu_vel.first, temp.first);
		accu_vel.second += temp.second;
	}
	accu_vel.first = pamhd::mul(accu_vel.first, 1 / accu_vel.second);

	for (size_t dim = 0; dim < 3; dim++) {
		if (abs(accu_vel.first[dim] - bulk_velocity[dim]) / abs(bulk_velocity[dim]) > 1e-1) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << "): " << dim << ", "
				<< accu_vel.first[dim] << " " << bulk_velocity[dim]
				<< std::endl;
			return false;
		}
	}

	// relative kinetic energy
	double accu_ekin = 0;
	for (const auto& part: particles) {
		accu_ekin += get_accumulated_value(
			0.5 * part[Mass()] * pamhd::norm2(pamhd::add(part[Velocity()], pamhd::neg(bulk_velocity))),
			particle_min, particle_max,
			cell_min, cell_max
		);
	}

	// temperature
	const double
		temp = 2 * accu_ekin / 3.0 / accu_nr_real / particle_temp_nrj_ratio,
		ref_temp = pamhd::sum(temperature) / 3;
	if (abs(temp - ref_temp) / abs(ref_temp) > 1e-2) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
			<< temp << " " << ref_temp
			<< std::endl;
		return false;
	}

	// pressure
	const double
		ref_pres = ref_temp * particle_temp_nrj_ratio * ref_nr_real / volume,
		pres = 2 * accu_ekin / 3.0 / volume;
	if (abs(pres - ref_pres) / abs(ref_pres) > 1e-2) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
			<< pres << " " << ref_pres
			<< std::endl;
		return false;
	}

	return true;
}


int main()
{
	if (
		not test_accumulation(
			1e5, // nr_particles
			{-5, 4, 9}, // part vol
			{-1, 8, 33},
			{-5, 4, 9}, // cell vol
			{-1, 8, 33},
			0.25, // boltzmann
			{1, 2, 3}, // bulk velocity
			{1, 1, 1}, // bulk temperature
			123, // total mass
			1, // species mass
			997 // seed
		)
	) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
		return EXIT_FAILURE;
	}

	if (
		not test_accumulation(
			1e5, // nr_particles
			{-5, 4, 9}, // part vol
			{-1, 8, 33},
			{-3, 6, 7}, // cell vol
			{ 1, 9, 11},
			0.25, // boltzmann
			{1, 2, 3}, // bulk velocity
			{1, 1, 1}, // bulk temperature
			123, // total mass
			1, // species mass
			997 // seed
		)
	) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
		return EXIT_FAILURE;
	}

	if (
		not test_accumulation(
			1e5, // nr_particles
			{-5, 4, 9}, // part vol
			{-1, 8, 33},
			{-5, 4, 9}, // cell vol
			{-1, 8, 33},
			3, // boltzmann
			{-3, -2, -1}, // bulk velocity
			{2, 2, 2}, // bulk temperature
			47, // total mass
			0.1, // species mass
			998 // seed
		)
	) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
		return EXIT_FAILURE;
	}

	constexpr auto Re = 6.371e6;
	if (
		not test_accumulation(
			1e5, // nr_particles
			{0, 0, 0}, // part vol
			{Re, Re, Re},
			{0, 0, 0}, // cell vol
			{2*Re, 2*Re, 2*Re},
			1.38064852e-23, // boltzmann
			{-5e5, 1e4, -1e4}, // bulk velocity
			{1e5, 1e5, 1e5}, // bulk temperature
			432534, // total mass
			1.6726217770000001e-27, // species mass
			999 // seed
		)
	) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
