/*
Two population version rusanov.hpp

Copyright 2016, 2017 Ilja Honkonen
Copyright 2025 Finnish Meteorological Institute
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

* Neither the name of copyright holders nor the names of their contributors
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


#ifndef PAMHD_MHD_N_RUSANOV_HPP
#define PAMHD_MHD_N_RUSANOV_HPP


#include "cmath"
#include "limits"
#include "string"
#include "tuple"

#include "mhd/common.hpp"
#include "mhd/variables.hpp"


namespace pamhd {
namespace mhd {


/*!
Multi-population version of get_flux_rusanov() in rusanov.hpp.

Splits returned flux into contributions from state_neg and state_pos
to all populations based on their fraction of mass vs total mass.
*/
template <
	class Mass_Density,
	class Momentum_Density,
	class Total_Energy_Density,
	class Magnetic_Field,
	class Vector,
	class Scalar
> std::tuple<detail::MHD, detail::MHD, Scalar> get_flux_N_rusanov(
	detail::MHD& state_neg,
	detail::MHD& state_pos,
	const Vector& bg_face_magnetic_field,
	const Scalar& area,
	const Scalar& dt,
	const Scalar& adiabatic_index,
	const Scalar& vacuum_permeability
) {
	using std::isnormal;
	using std::isfinite;
	using std::to_string;

	const Mass_Density Mas{};
	const Momentum_Density Mom{};
	const Total_Energy_Density Nrj{};
	const Magnetic_Field Mag{};

	if (not isnormal(state_neg[Mas]) or state_neg[Mas] < 0) {
		throw std::domain_error(
			"Invalid mass density in state_neg: "
			+ to_string(state_neg[Mas])
		);
	}
	if (not isnormal(state_pos[Mas]) or state_pos[Mas] < 0) {
		throw std::domain_error(
			"Invalid mass density in state_pos: "
			+ to_string(state_pos[Mas])
		);
	}

	const auto
		fast_magnetosonic_neg
			= get_fast_magnetosonic_speed(
				state_neg[Mas],
				state_neg[Mom],
				state_neg[Nrj],
				state_neg[Mag],
				bg_face_magnetic_field,
				adiabatic_index,
				vacuum_permeability
			),

		fast_magnetosonic_pos
			= get_fast_magnetosonic_speed(
				state_pos[Mas],
				state_pos[Mom],
				state_pos[Nrj],
				state_pos[Mag],
				bg_face_magnetic_field,
				adiabatic_index,
				vacuum_permeability
			),

		max_fast_ms = std::max(fast_magnetosonic_neg, fast_magnetosonic_pos),
		max_signal = std::max(
			std::abs(get_velocity(state_neg[Mom], state_neg[Mas])[0]) + max_fast_ms,
			std::abs(get_velocity(state_pos[Mom], state_pos[Mas])[0]) + max_fast_ms
		);

	detail::MHD flux_neg, flux_pos;
	std::tie(flux_neg, flux_pos) = N_get_flux<
		Mass_Density,
		Momentum_Density,
		Total_Energy_Density,
		Magnetic_Field
	>(
		state_neg,
		state_pos,
		bg_face_magnetic_field,
		adiabatic_index,
		vacuum_permeability
	);

	flux_neg[Mas] = 0.5 * area * dt
		* (flux_neg[Mas] + state_neg[Mas] * max_signal);
	flux_neg[Mom] = pamhd::mul(0.5 * area * dt,
		pamhd::add(
			flux_neg[Mom],
			pamhd::mul(state_neg[Mom], max_signal)));
	flux_neg[Nrj] = 0.5 * area * dt
		* (flux_neg[Nrj] + state_neg[Nrj] * max_signal);
	flux_neg[Mag] = pamhd::mul(0.5 * area * dt,
		pamhd::add(
			flux_neg[Mag],
			pamhd::mul(state_neg[Mag], max_signal)));

	flux_pos[Mas] = 0.5 * area * dt
		* (flux_pos[Mas] - state_pos[Mas] * max_signal);
	flux_pos[Mom] = pamhd::mul(0.5 * area * dt,
		pamhd::add(
			flux_pos[Mom],
			pamhd::mul(state_pos[Mom], -max_signal)));
	flux_pos[Nrj] = 0.5 * area * dt
		* (flux_pos[Nrj] - state_pos[Nrj] * max_signal);
	flux_pos[Mag] = pamhd::mul(0.5 * area * dt,
		pamhd::add(
			flux_pos[Mag],
			pamhd::mul(state_pos[Mag], -max_signal)));

	return std::make_tuple(flux_neg, flux_pos, max_signal);
}

}} // namespaces

#endif // ifndef PAMHD_MHD_N_RUSANOV_HPP
