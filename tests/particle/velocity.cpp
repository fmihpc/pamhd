/*
Tests bulk velocity calculation of particles in PAMHD.

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
#include "random"

#include "prettyprint.hpp"

#include "particle/common.hpp"
#include "particle/variables.hpp"


using namespace std;
using namespace pamhd::particle;


int main()
{
	using std::array;
	using std::min;
	using std::pow;
	using std::sqrt;

	const size_t nr_particles = 100000;
	const array<double, 3>
		volume_min{-5, 4, 9},
		volume_max{-1, 8, 33};

	std::mt19937 random_source;
	random_source.seed(997);


	const array<double, 3>
		ref_bulk_vel1{0, 0, 0},
		temperature1{1, 2, 3};
	const double
		species_mass1 = 1,
		total_mass1 = nr_particles;
	const auto particles1
		= create_particles<
			Particle_Internal,
			Mass,
			Charge_Mass_Ratio,
			Position,
			Velocity,
			Particle_ID,
			Species_Mass
		>(
			ref_bulk_vel1,
			volume_min,
			volume_max,
			temperature1,
			nr_particles,
			10,
			total_mass1,
			species_mass1,
			0.25,
			random_source
		);
	const auto bulk_vel1 = get_bulk_velocity<Mass, Velocity, Species_Mass>(particles1);
	if (pamhd::norm(
		pamhd::add(bulk_vel1, pamhd::neg(ref_bulk_vel1))
	) > 1e-2) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
			<< ref_bulk_vel1 << " " << bulk_vel1
			<< std::endl;
		abort();
	}


	const array<double, 3>
		ref_bulk_vel2{1, 2, 3},
		temperature2{3, 4, 5};
	const double
		species_mass2 = 1,
		total_mass2 = nr_particles;
	const auto particles2
		= create_particles<
			Particle_Internal,
			Mass,
			Charge_Mass_Ratio,
			Position,
			Velocity,
			Particle_ID,
			Species_Mass
		>(
			ref_bulk_vel2,
			volume_min,
			volume_max,
			temperature2,
			nr_particles,
			10,
			total_mass2,
			species_mass2,
			0.25,
			random_source
		);
	const auto bulk_vel2 = get_bulk_velocity<Mass, Velocity, Species_Mass>(particles2);
	if (pamhd::norm(
		pamhd::add(bulk_vel2, pamhd::neg(ref_bulk_vel2))
	) > 1e-2) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
			<< ref_bulk_vel2 << " " << bulk_vel2
			<< std::endl;
		abort();
	}

	auto particles3 = particles1;
	particles3.insert(particles3.end(), particles2.cbegin(), particles2.cend());
	const auto ref_bulk_vel3 = pamhd::mul(0.5, pamhd::add(ref_bulk_vel1, ref_bulk_vel2));
	const auto bulk_vel3 = get_bulk_velocity<Mass, Velocity, Species_Mass>(particles3);
	if (pamhd::norm(
		pamhd::add(bulk_vel3, pamhd::neg(ref_bulk_vel3))
	) > 1e-2) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
			<< ref_bulk_vel3 << " " << bulk_vel3
			<< std::endl;
		abort();
	}


	const double
		species_mass4 = 1,
		total_mass4 = nr_particles / 100;
	const auto particles4
		= create_particles<
			Particle_Internal,
			Mass,
			Charge_Mass_Ratio,
			Position,
			Velocity,
			Particle_ID,
			Species_Mass
		>(
			ref_bulk_vel2,
			volume_min,
			volume_max,
			temperature2,
			nr_particles,
			10,
			total_mass4,
			species_mass4,
			0.25,
			random_source
		);
	auto particles5 = particles1;
	particles5.insert(particles5.end(), particles4.cbegin(), particles4.cend());
	const auto ref_bulk_vel5_denominator
		= total_mass1 / species_mass1 + total_mass4 / species_mass4;
	const auto ref_bulk_vel5 = pamhd::add(
		pamhd::mul(ref_bulk_vel1, total_mass1 / species_mass1 / ref_bulk_vel5_denominator),
		pamhd::mul(ref_bulk_vel2, total_mass4 / species_mass4 / ref_bulk_vel5_denominator));
	const auto bulk_vel5 = get_bulk_velocity<Mass, Velocity, Species_Mass>(particles5);
	if (pamhd::norm(
		pamhd::add(bulk_vel5, pamhd::neg(ref_bulk_vel5))
	) > 1e-2) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
			<< ref_bulk_vel5 << " " << bulk_vel5
			<< std::endl;
		abort();
	}


	return EXIT_SUCCESS;
}
