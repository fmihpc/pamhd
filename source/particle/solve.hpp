/*
Particle propagator of PAMHD.

Copyright 2014, 2015, 2016, 2017 Ilja Honkonen
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
*/

#ifndef PAMHD_PARTICLE_SOLVE_HPP
#define PAMHD_PARTICLE_SOLVE_HPP


#include "utility"

#include "common_functions.hpp"


namespace pamhd {
namespace particle {


/*!
Returns new position and velocity for given particle after given length of time.

Particle accelerator based on translations and a rotation given in
equations 25 and 26 of:
J.P. Boris, Relativistic Plasma Simulation - Optimization of a Hybrid Code,
Proceedings of the conference on the numerical simulation of plasmas...,
available at http://www.dtic.mil/dtic/tr/fulltext/u2/a023511.pdf
Same notation is used here.
*/
auto propagate(
	const auto& position,
	const auto& velocity,
	const auto& electric_field,
	const auto& magnetic_field,
	const double charge_mass_ratio,
	const double time_step
) {
	using std::pow;
	using std::tan;

	const double
		coeff = charge_mass_ratio * time_step / 2.0,
		B_mag = norm(magnetic_field),
		B_mag_non_zero = [B_mag](){if (B_mag != 0) return B_mag; else return 1.0;}(),
		f1 = tan(coeff * B_mag) / B_mag_non_zero,
		f2 = 2 * f1 / (1 + dot(f1, f1) * dot(B_mag, B_mag));

	const auto
		new_pos = add(position, mul(time_step, velocity)),
		v1 = add(velocity, mul(coeff, electric_field)),
		v2 = add(v1, mul(f1, cross(v1, magnetic_field))),
		v3 = add(v1, mul(f2, cross(v2, magnetic_field))),
		new_vel = add(v3, mul(coeff, electric_field));

	return std::make_pair(new_pos, new_vel);
}


}} // namespaces

#endif
