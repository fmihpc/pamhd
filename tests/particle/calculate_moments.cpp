/*
Tests for bulk value functions of PAMHD.

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
	using std::abs;
	using std::array;
	using std::min;
	using std::pow;
	using std::sqrt;

	size_t nr_particles = 10000;
	array<double, 3>
		bulk_velocity_ref{0, -3, 8},
		volume_min_ref{-2, 4, 9},
		volume_max_ref{-1, 200, 9.5},
		temperature_ref{3, 2 * M_PI, 42};
	double
		volume
			= (volume_max_ref[0] - volume_min_ref[0])
			* (volume_max_ref[1] - volume_min_ref[1])
			* (volume_max_ref[2] - volume_min_ref[2]),
		total_mass_ref = 0,
		species_mass = 1.25,
		particle_temp_nrj_ratio = 0.25,
		scalar_temp_ref
			= (
				temperature_ref[0]
				+ temperature_ref[1]
				+ temperature_ref[2]
			) / 3;

	std::mt19937 random_source;
	random_source.seed(997);

	// massless particles
	double pressure_ref
		= nr_particles // of mass species_mass
		/ volume
		* particle_temp_nrj_ratio
		* scalar_temp_ref;
	auto particles
		= create_particles<
			Particle_Internal,
			Mass,
			Charge_Mass_Ratio,
			Position,
			Velocity,
			Particle_ID,
			Species_Mass
		>(
			bulk_velocity_ref,
			volume_min_ref,
			volume_max_ref,
			temperature_ref,
			nr_particles,
			10,
			total_mass_ref,
			species_mass,
			particle_temp_nrj_ratio,
			random_source
		);

	array<double, 3> bulk_velocity = get_bulk_velocity<Mass, Velocity, Species_Mass>(particles);

	if (abs(bulk_velocity[0] - bulk_velocity_ref[0]) > 0.1) {
		std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
			<< "Incorrect x component of bulk velocity for particles: "
			<< bulk_velocity[0] << ", should be " << bulk_velocity_ref[0]
			<< endl;
		abort();
	}
	if (abs(bulk_velocity[1] - bulk_velocity_ref[1]) > 0.1) {
		std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
			<< "Incorrect y component of bulk velocity for particles: "
			<< bulk_velocity[1] << ", should be " << bulk_velocity_ref[1]
			<< endl;
		abort();
	}
	if (abs(bulk_velocity[2] - bulk_velocity_ref[2]) > 0.1) {
		std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
			<< "Incorrect z component of bulk velocity for particles: "
			<< bulk_velocity[2] << ", should be " << bulk_velocity_ref[2]
			<< endl;
		abort();
	}

	auto scalar_temperature
		= get_temperature<
			Mass,
			Velocity,
			Species_Mass
		>(
			particles,
			particle_temp_nrj_ratio
		);

	if (abs(scalar_temperature - scalar_temp_ref) > 0.1) {
		std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
			<< "Incorrect scalar temperature for particles: "
			<< scalar_temperature << ", should be " << scalar_temp_ref
			<< endl;
		abort();
	}

	auto pressure
		= get_pressure<
			Mass,
			Velocity,
			Species_Mass
		>(
			particles,
			particle_temp_nrj_ratio,
			volume
		);

	if (abs(pressure - pressure_ref) > 4) {
		std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
			<< "Incorrect pressure for particles: "
			<< pressure << ", should be " << pressure_ref
			<< endl;
		abort();
	}

	total_mass_ref = 314;
	species_mass = 0.1;
	particle_temp_nrj_ratio = 0.01;

	particles
		= create_particles<
			Particle_Internal,
			Mass,
			Charge_Mass_Ratio,
			Position,
			Velocity,
			Particle_ID,
			Species_Mass
		>(
			bulk_velocity_ref,
			volume_min_ref,
			volume_max_ref,
			temperature_ref,
			nr_particles,
			10,
			total_mass_ref,
			species_mass,
			particle_temp_nrj_ratio,
			random_source
		);

	bulk_velocity = get_bulk_velocity<Mass, Velocity, Species_Mass>(particles);

	if (abs(bulk_velocity[0] - bulk_velocity_ref[0]) > 0.1) {
		std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
			<< "Incorrect x component of bulk velocity for particles: "
			<< bulk_velocity[0] << ", should be " << bulk_velocity_ref[0]
			<< endl;
		abort();
	}
	if (abs(bulk_velocity[1] - bulk_velocity_ref[1]) > 0.1) {
		std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
			<< "Incorrect y component of bulk velocity for particles: "
			<< bulk_velocity[1] << ", should be " << bulk_velocity_ref[1]
			<< endl;
		abort();
	}
	if (abs(bulk_velocity[2] - bulk_velocity_ref[2]) > 0.1) {
		std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
			<< "Incorrect z component of bulk velocity for particles: "
			<< bulk_velocity[2] << ", should be " << bulk_velocity_ref[2]
			<< endl;
		abort();
	}

	scalar_temperature
		= get_temperature<
			Mass,
			Velocity,
			Species_Mass
		>(
			particles,
			particle_temp_nrj_ratio
		);

	if (abs(scalar_temperature - scalar_temp_ref) > 0.5) {
		std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
			<< "Incorrect scalar temperature for particles: "
			<< scalar_temperature << ", should be " << scalar_temp_ref
			<< endl;
		abort();
	}

	pressure_ref
		= total_mass_ref / species_mass
		/ volume
		* particle_temp_nrj_ratio
		* scalar_temp_ref;
	pressure
		= get_pressure<
			Mass,
			Velocity,
			Species_Mass
		>(
			particles,
			particle_temp_nrj_ratio,
			volume
		);
	if (abs(pressure - pressure_ref) > 4) {
		std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
			<< "Incorrect pressure for particles: "
			<< pressure << ", should be " << pressure_ref
			<< endl;
		abort();
	}


	// test that regardless of volume into which particles are created
	volume_min_ref[0] = -1000;
	volume_min_ref[1] = -1000;
	volume_min_ref[2] = -1000;
	volume_max_ref[0] = 1000;
	volume_max_ref[1] = 1000;
	volume_max_ref[2] = 1000;
	volume
		= (volume_max_ref[0] - volume_min_ref[0])
		* (volume_max_ref[1] - volume_min_ref[1])
		* (volume_max_ref[2] - volume_min_ref[2]);
	particles
		= create_particles<
			Particle_Internal,
			Mass,
			Charge_Mass_Ratio,
			Position,
			Velocity,
			Particle_ID,
			Species_Mass
		>(
			bulk_velocity_ref,
			volume_min_ref,
			volume_max_ref,
			temperature_ref,
			nr_particles,
			10,
			total_mass_ref,
			species_mass,
			particle_temp_nrj_ratio,
			random_source
		);

	scalar_temperature
		= get_temperature<
			Mass,
			Velocity,
			Species_Mass
		>(
			particles,
			particle_temp_nrj_ratio
		);

	if (abs(scalar_temperature - scalar_temp_ref) > 0.5) {
		std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
			<< "Incorrect scalar temperature for particles: "
			<< scalar_temperature << ", should be " << scalar_temp_ref
			<< endl;
		abort();
	}

	pressure_ref
		= total_mass_ref / species_mass
		/ volume
		* particle_temp_nrj_ratio
		* scalar_temp_ref;
	pressure
		= get_pressure<
			Mass,
			Velocity,
			Species_Mass
		>(
			particles,
			particle_temp_nrj_ratio,
			volume
		);
	if (abs(pressure - pressure_ref) > 4) {
		std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
			<< "Incorrect pressure for particles: "
			<< pressure << ", should be " << pressure_ref
			<< endl;
		abort();
	}


	return EXIT_SUCCESS;
}
