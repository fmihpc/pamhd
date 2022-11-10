/*
Rusanov MHD solver of PAMHD.

Copyright 2016, 2017 Ilja Honkonen
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
*/


#ifndef PAMHD_MHD_RUSANOV_HPP
#define PAMHD_MHD_RUSANOV_HPP


#include "cmath"
#include "limits"
#include "string"
#include "tuple"

#include "mhd/common.hpp"


namespace pamhd {
namespace mhd {


/*!
See get_flux_hll() in hll_athena.hpp

Equation 10.55 in ISBN 978-3-540-25202-3
*/
template <
	class Mass_Density,
	class Momentum_Density,
	class Total_Energy_Density,
	class Magnetic_Field,
	class MHD,
	class Vector,
	class Scalar
> std::tuple<MHD, Scalar> get_flux_rusanov(
	MHD& state_neg,
	MHD& state_pos,
	const Vector& bg_face_magnetic_field,
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
			std::abs(state_neg[Mom][0] / state_neg[Mas]) + max_fast_ms,
			std::abs(state_pos[Mom][0] / state_pos[Mas]) + max_fast_ms
		);

	if (not isfinite(max_fast_ms)) {
		throw std::domain_error(
			"Invalid maximum fast magnetosonic speed: "
			+ to_string(max_fast_ms)
		);
	}
	if (not isfinite(max_signal)) {
		throw std::domain_error(
			"Invalid maximum signal speed: "
			+ to_string(max_signal)
		);
	}

	const auto Mas_getter
		= [](MHD& state) -> typename Mass_Density::data_type& {
			return state[Mass_Density()];
		};
	const auto Mom_getter
		= [](MHD& state) -> typename Momentum_Density::data_type& {
			return state[Momentum_Density()];
		};
	const auto Nrj_getter
		= [](MHD& state) -> typename Total_Energy_Density::data_type& {
			return state[Total_Energy_Density()];
		};
	const auto Mag_getter
		= [](MHD& state) -> typename Magnetic_Field::data_type& {
			return state[Magnetic_Field()];
		};

	MHD
		flux_neg
			= get_flux(
				state_neg,
				bg_face_magnetic_field,
				adiabatic_index,
				vacuum_permeability,
				Mas_getter, Mom_getter, Nrj_getter, Mag_getter
			),
		flux_pos
			= get_flux(
				state_pos,
				bg_face_magnetic_field,
				adiabatic_index,
				vacuum_permeability,
				Mas_getter, Mom_getter, Nrj_getter, Mag_getter
			);
	if (not isfinite(flux_neg[Mas])) {
		throw std::domain_error(
			"Invalid mass density flux_neg: "
			+ to_string(flux_neg[Mas])
		);
	}
	if (not isfinite(flux_pos[Mas])) {
		throw std::domain_error(
			"Invalid mass density flux_pos: "
			+ to_string(flux_pos[Mas])
		);
	}

	MHD flux
		= (flux_neg + flux_pos) / 2
		- (state_pos - state_neg) * (max_signal / 2);

	if (not isfinite(flux[Mas])) {
		throw std::domain_error(
			"Invalid mass density flux: "
			+ to_string(flux[Mas])
		);
	}

	return std::make_tuple(flux, max_signal);
}


}} // namespaces

#endif // ifndef PAMHD_MHD_RUSANOV_HPP
