/*
Particle-assisted version of hll_athena.hpp.
Copyright 2003 James M. Stone
Copyright 2003 Thomas A. Gardiner
Copyright 2003 Peter J. Teuben
Copyright 2003 John F. Hawley
Copyright 2014, 2015, 2016, 2017 Ilja Honkonen
Copyright 2025 Finnish Meteorological Institute

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PAMHD_MHD_N_HLL_ATHENA_HPP
#define PAMHD_MHD_N_HLL_ATHENA_HPP


#include "cmath"
#include "limits"
#include "string"
#include "tuple"

#include "common_functions.hpp"
#include "mhd/common.hpp"
#include "mhd/variables.hpp"


namespace pamhd {
namespace mhd {
namespace athena {


/*!
Multi-population version of get_flux_hll() in hll_athena.hpp.

Splits returned flux into contributions from state_neg and state_pos
to all populations based on their fraction of mass vs total mass.

Ignores background magnetic field.
*/
template <
	class Mass_Density,
	class Momentum_Density,
	class Total_Energy_Density,
	class Magnetic_Field,
	class Vector,
	class Scalar
> std::tuple<detail::MHD, detail::MHD, Scalar> get_flux_N_hll(
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

	const auto
		flow_v_neg = get_velocity(state_neg[Mom], state_neg[Mas]),
		flow_v_pos = get_velocity(state_pos[Mom], state_pos[Mas]);

	const auto
		pressure_thermal_neg
			= get_pressure(
				state_neg[Mas],
				state_neg[Mom],
				state_neg[Nrj],
				state_neg[Mag],
				adiabatic_index,
				vacuum_permeability
			),

		pressure_thermal_pos
			= get_pressure(
				state_pos[Mas],
				state_pos[Mom],
				state_pos[Nrj],
				state_pos[Mag],
				adiabatic_index,
				vacuum_permeability
			),

		pressure_magnetic_neg
			= pamhd::norm2(state_neg[Mag]) / (2 * vacuum_permeability),

		pressure_magnetic_pos
			= pamhd::norm2(state_pos[Mag]) / (2 * vacuum_permeability),

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

		max_signal = std::max(fast_magnetosonic_neg, fast_magnetosonic_pos),

		max_signal_neg
			= (flow_v_neg[0] <= flow_v_pos[0])
			? flow_v_neg[0] - max_signal
			: flow_v_pos[0] - max_signal,

		max_signal_pos
			= (flow_v_neg[0] <= flow_v_pos[0])
			? flow_v_pos[0] + max_signal
			: flow_v_neg[0] + max_signal,

		bm = std::min(max_signal_neg, 0.0),
		bp = std::max(max_signal_pos, 0.0);

	if (not isnormal(pressure_thermal_neg) or pressure_thermal_neg < 0) {
		throw std::domain_error(
			"Invalid thermal pressure in state_neg: "
			+ to_string(pressure_thermal_neg)
		);
	}
	if (not isfinite(pressure_magnetic_neg) or pressure_magnetic_neg < 0) {
		throw std::domain_error(
			"Invalid magnetic pressure in state_neg: "
			+ to_string(pressure_magnetic_neg)
		);
	}
	if (not isfinite(max_signal_neg)) {
		throw std::domain_error(
			"Invalid max signal speed in state_neg: " + to_string(max_signal_neg)
			+ ", max signal: " + to_string(max_signal)
			+ ", flow_v_neg[0]: " + to_string(flow_v_neg[0])
		);
	}

	if (not isnormal(pressure_thermal_pos) or pressure_thermal_pos < 0) {
		throw std::domain_error(
			"Invalid thermal pressure in state_pos: "
			+ to_string(pressure_thermal_pos)
		);
	}
	if (not isfinite(pressure_magnetic_pos) or pressure_magnetic_pos < 0) {
		throw std::domain_error(
			"Invalid magnetic pressure in state_pos: "
			+ to_string(pressure_magnetic_pos)
		);
	}
	if (not isfinite(max_signal_pos)) {
		throw std::domain_error(
			"Invalid max signal speed in state_pos: "
			+ to_string(max_signal_pos)
		);
	}

	if (not isnormal(bp - bm) or bp - bm < 0) {
		detail::MHD flux;

		flux[Mas]    =
		flux[Mom][0] =
		flux[Mom][1] =
		flux[Mom][2] =
		flux[Nrj]    =
		flux[Mag][0] =
		flux[Mag][1] =
		flux[Mag][2] = 0;

		return std::make_tuple(flux, flux, 0);
	}


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

	const auto
		neg_factor = bp / (bp - bm) * area * dt,
		pos_factor = bm / (bm - bp) * area * dt;
	flux_neg[Mas] = neg_factor * (flux_neg[Mas] - state_neg[Mas] * bm);
	flux_neg[Mom] = pamhd::mul(neg_factor,
		pamhd::add(flux_neg[Mom], pamhd::mul(state_neg[Mom], -bm)));
	flux_neg[Nrj] = neg_factor * (flux_neg[Nrj] - state_neg[Nrj] * bm);
	flux_neg[Mag] = pamhd::mul(neg_factor,
		pamhd::add(flux_neg[Mag], pamhd::mul(state_neg[Mag], -bm)));

	flux_pos[Mas] = pos_factor * (flux_pos[Mas] - state_pos[Mas] * bp);
	flux_pos[Mom] = pamhd::mul(pos_factor,
		pamhd::add(flux_pos[Mom], pamhd::mul(state_pos[Mom], -bp)));
	flux_pos[Nrj] = pos_factor * (flux_pos[Nrj] - state_pos[Nrj] * bp);
	flux_pos[Mom] = pamhd::mul(pos_factor,
		pamhd::add(flux_pos[Mom], pamhd::mul(state_pos[Mom], -bp)));

	return std::make_tuple(flux_neg, flux_pos, std::max(std::fabs(bp), std::fabs(bm)));
}

}}} // namespaces

#endif // ifndef PAMHD_MHD_N_HLL_ATHENA_HPP
