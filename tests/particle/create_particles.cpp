/*
Tests for create_particles function of PAMHD.

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
	using std::abs;
	using std::array;
	using std::min;
	using std::numeric_limits;
	using std::pow;
	using std::sqrt;

	// short hand notation for variables
	const Position Pos{};
	const Velocity Vel{};
	const Mass Mas{};

	const array<double, 3>
		bulk_velocity_ref{0, -3, 8},
		volume_min_ref{-2, 4, 9},
		volume_max_ref{-1, 200, 9.5},
		temperature_ref{3, 2 * M_PI, 42};
	const double
		total_mass_ref = 250,
		species_mass = 1.25,
		particle_temp_nrj_ratio = 0.25;

	for (size_t nr_of_particles = 4000; nr_of_particles <= 32000; nr_of_particles *= 2) {

		std::mt19937 random_source;
		random_source.seed(1);

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
				bulk_velocity_ref,
				volume_min_ref,
				volume_max_ref,
				temperature_ref,
				nr_of_particles,
				10,
				total_mass_ref,
				species_mass,
				particle_temp_nrj_ratio,
				random_source
			);

		if (particles.size() != nr_of_particles) {
			std::cerr <<  __FILE__ << ":" << __LINE__
				<< " Incorrect number of particles created: " << particles.size()
				<< std::endl;
			abort();
		}

		array<double, 3>
			bulk_velocity{0, 0, 0},
			volume_min{
				numeric_limits<double>::max(),
				numeric_limits<double>::max(),
				numeric_limits<double>::max()
			},
			volume_max{
				numeric_limits<double>::lowest(),
				numeric_limits<double>::lowest(),
				numeric_limits<double>::lowest()
			},
			temperature{0, 0, 0};

		array<array<double, 3>, 3> non_gyrotropic_temperature{0, 0, 0, 0, 0, 0, 0, 0, 0};
		double total_mass = 0;

		for (const auto& particle: particles) {
			bulk_velocity = pamhd::add(bulk_velocity, particle[Vel]);
			total_mass += particle[Mas];

			for (size_t dim = 0; dim < 3; dim++) {
				if (volume_min[dim] > particle[Pos][dim]) {
					volume_min[dim] = particle[Pos][dim];
				}
				if (volume_max[dim] < particle[Pos][dim]) {
					volume_max[dim] = particle[Pos][dim];
				}
			}
		}

		if (abs(total_mass - total_mass_ref) > 1e-20) {
			std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
				<< "Incorrect total mass for created particles: " << total_mass
				<< ", should be " << total_mass_ref
				<< endl;
			abort();
		}

		bulk_velocity = pamhd::mul(bulk_velocity, 1.0 / particles.size());

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

		for (const auto& particle: particles) {
			for (size_t dim1 = 0; dim1 < 3; dim1++) {
				for (size_t dim2 = 0; dim2 < 3; dim2++) {
					if (dim1 == dim2) {
						temperature[dim1]
							+= species_mass
							* pow(particle[Vel][dim1] - bulk_velocity[dim1], 2)
							/ particle_temp_nrj_ratio;
					} else {
						non_gyrotropic_temperature[dim1][dim2]
							+= species_mass
							* (particle[Vel][dim1] - bulk_velocity[dim1])
							* (particle[Vel][dim2] - bulk_velocity[dim2])
							/ particle_temp_nrj_ratio;
					}
				}
			}
		}
		temperature = pamhd::mul(temperature, 1.0 / particles.size());
		for (auto& x: non_gyrotropic_temperature)
		for (auto& y: x) y /= particles.size();

		for (size_t dim1 = 0; dim1 < 3; dim1++) {
			for (size_t dim2 = 0; dim2 < 3; dim2++) {
				if (dim1 == dim2) {
					if (
						abs(temperature[dim1] - temperature_ref[dim1])
						/ temperature_ref[dim1]
						> 0.1
					) {
					std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
							<< "Incorrect component number " << dim1
							<< " of temperature for particles: " << temperature[dim1]
							<< ", should be " << temperature_ref[dim1]
							<< endl;
						abort();
					}
				} else {
					if (abs(non_gyrotropic_temperature[dim1][dim2]) > 0.3) {
					std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
							<< "Too large non-gyrotropic temperature component "
							<< dim1 << ", " << dim2 << " of temperature for particles: "
							<< non_gyrotropic_temperature[dim1][dim2] << ", should be 0"
							<< endl;
						abort();
					}
				}
			}
		}
	}

	return EXIT_SUCCESS;
}
