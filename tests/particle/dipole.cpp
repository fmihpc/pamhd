/*
Internal propagator version of PAMHD dipole_odeint test.

Copyright 2015, 2016, 2017 Ilja Honkonen
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


Propagation test after http://dx.doi.org/10.1029/2005JA011382
*/

#include "array"
#include "cmath"
#include "cstdlib"
#include "iomanip"
#include "iostream"

#include "particle/solve.hpp"


using namespace std;
using namespace pamhd::particle;


constexpr double
	proton_mass = 1.672621777e-27,
	elementary_charge = 1.602176565e-19,
	proton_charge_mass_ratio = elementary_charge / proton_mass,
	earth_radius = 6.371e6;

const std::array<double, 3> dipole_moment{0, 0, -8e22};

// state[0] = particle position, state[1] = velocity
using state_t = std::array<std::array<double, 3>, 2>;


std::array<double, 3> get_earth_dipole(const std::array<double, 3>& position)
{
	const double distance = pamhd::norm(position);
	const std::array<double, 3>
		direction = pamhd::mul(position, 1 / distance),
		projected_dip_mom = pamhd::mul(pamhd::dot(dipole_moment, direction), direction);

	return pamhd::mul(pamhd::add(pamhd::mul(-3, projected_dip_mom), dipole_moment), -1e-7 / std::pow(distance, 3));
}


bool check_result(
	const double initial_energy,
	const double initial_angle,
	const size_t step_divisor,
	const double nrj_error,
	const double mom_error,
	const std::array<double, 3>& position
) {
	if (nrj_error > 1e-12) {
		return false;
	}

	if (pamhd::norm(position) > 7 * earth_radius) {
		if (initial_angle == 90) {
			if (initial_energy == 1e4) {
				if (step_divisor >= 2) {
					return false;
				}
			} else {
				if (step_divisor >= 512) {
					return false;
				}
			}
		} else {
			if (initial_energy == 1e4) {
				if (step_divisor >= 4) {
					return false;
				}
			} else {
				if (step_divisor >= 128) {
					return false;
				}
			}
		}
	}

	if (initial_energy == 1e4) {
		if (initial_angle == 90) {
			if (step_divisor >= 128 and mom_error > 2e-2) {
				return false;
			}

			return true;

		} else {
			if (step_divisor >= 128 and mom_error > 4e-2) {
				return false;
			}

			return true;
		}

	} else {

		if (initial_angle == 90) {
			if (step_divisor >= 1024 and mom_error > 0.6) {
				return false;
			}

			return true;

		} else {
			if (step_divisor >= 512 and mom_error > 2.1) {
				return false;
			}

			return true;
		}
	}

	return true;
}

int main()
{
	using std::cos;
	using std::abs;
	using std::sin;

	const std::array<double, 3>
		guiding_center_start{5 * earth_radius, 0, 0},
		field_at_start{get_earth_dipole(guiding_center_start)};

	const double gyroperiod
		= 2 * M_PI / abs(proton_charge_mass_ratio) / pamhd::norm(field_at_start);

	/*cout <<
		"Energy (eV), pitch angle (deg), steps/gyroperiod: "
		"max relative error in energy, magnetic moment\n";*/
	for (auto kinetic_energy: {1e4, 1e7}) { // in eV
	for (auto pitch_angle: {90.0, 30.0}) { // V from B, degrees

		const double particle_speed
			= sqrt(2 * kinetic_energy * proton_charge_mass_ratio);

		const std::array<double, 3> initial_velocity{
			0,
			sin(pitch_angle / 180 * M_PI) * particle_speed,
			cos(pitch_angle / 180 * M_PI) * particle_speed
		};

		const double gyroradius
			= initial_velocity[1]
			/ proton_charge_mass_ratio
			/ pamhd::norm(field_at_start);

		const std::array<double, 3>
			offset_from_gc{gyroradius, 0, 0},
			initial_position{pamhd::add(guiding_center_start, offset_from_gc)};

		const double
			initial_energy = 0.5 * proton_mass * pamhd::dot(initial_velocity, initial_velocity),
			initial_magnetic_moment
				= proton_mass
				* std::pow(initial_velocity[1], 2)
				/ (2 * pamhd::norm(get_earth_dipole(initial_position)));

		for (size_t step_divisor = 1; step_divisor <= (1 << 11); step_divisor *= 2) {

			state_t state{{initial_position, initial_velocity}};
			const double time_step = gyroperiod / step_divisor;

			constexpr double propagation_time = 424.1; // seconds, 1/100 of doi above

			double
				nrj_error = 0,
				mom_error = 0;
			for (double time = 0; time < propagation_time; time += time_step) {

				std::tie(
					state[0],
					state[1]
				) = propagate(
					state[0],
					state[1],
					std::array<double, 3>{0, 0, 0},
					get_earth_dipole(state[0]),
					proton_charge_mass_ratio,
					time_step
				);

				const std::array<double, 3>
					v{state[1]},
					b{get_earth_dipole(state[0])},
					v_perp_b{pamhd::add(v, pamhd::mul(-pamhd::dot(v, b) / pamhd::dot(b, b), b))};

				const double
					current_energy
						= 0.5 * proton_mass * pamhd::norm2(state[1]),
					current_magnetic_moment
						= proton_mass * pamhd::norm2(v_perp_b) / (2 * pamhd::norm(b));
				nrj_error = max(
					nrj_error,
					abs(current_energy - initial_energy) / initial_energy
				);
				mom_error = max(
					mom_error,
					abs(current_magnetic_moment - initial_magnetic_moment)
						/ initial_magnetic_moment
				);
			}

			/*cout
				<< kinetic_energy << " "
				<< pitch_angle << " "
				<< step_divisor << ": "
				<< nrj_error << " "
				<< mom_error
				<< endl;*/

			if (
				not check_result(
					kinetic_energy,
					pitch_angle,
					step_divisor,
					nrj_error,
					mom_error,
					state[0]
				)
			) {
				std::cerr <<  __FILE__ << " (" << __LINE__ << "): "
					<< "Unexpected failure in propagator."
					<< std::endl;
				return EXIT_FAILURE;
			}
		}
	}}

	return EXIT_SUCCESS;
}
