/*
HLL MHD solver of PAMHD adapted from Athena (https://trac.princeton.edu/Athena/wiki)
Copyright 2003 James M. Stone
Copyright 2003 Thomas A. Gardiner
Copyright 2003 Peter J. Teuben
Copyright 2003 John F. Hawley
Copyright 2014, 2015, 2016, 2017 Ilja Honkonen

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

#ifndef PAMHD_MHD_HLL_ATHENA_HPP
#define PAMHD_MHD_HLL_ATHENA_HPP


#include "cmath"
#include "limits"
#include "string"
#include "tuple"

#include "mhd/common.hpp"


namespace pamhd {
namespace mhd {
namespace athena {


/*!
Returns the flux between states and maximum signal speed from cells' shared face.

\param [Mass_Density] Used to access mass density of MHD states
\param [Momentum_Density] Used to access momentum density of MHD states
\param [Total_Energy_Density] Used to access total energy density of MHD states
\param [Magnetic_Field] Used to access magnetic field of MHD states
\param [state_neg] MHD state in negative x direction from the face
\param [state_pos] MHD state in positive x direction from the face
\param [bg_face_magnetic_field] Background magnetic field at face
\param [adiabatic_index] en.wikipedia.org/wiki/Heat_capacity_ratio
\param [vacuum_permeability] en.wikipedia.org/wiki/Vacuum_permeability
*/
template <
	class Mass_Density,
	class Momentum_Density,
	class Total_Energy_Density,
	class Magnetic_Field,
	class MHD,
	class Vector,
	class Scalar
> std::tuple<MHD, Scalar> get_flux_hll(
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
			= state_neg[Mag].squaredNorm() / (2 * vacuum_permeability),

		pressure_magnetic_pos
			= state_pos[Mag].squaredNorm() / (2 * vacuum_permeability),

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
		MHD flux;

		flux[Mas]    =
		flux[Mom][0] =
		flux[Mom][1] =
		flux[Mom][2] =
		flux[Nrj]    =
		flux[Mag][0] =
		flux[Mag][1] =
		flux[Mag][2] = 0;

		return std::make_tuple(flux, 0);
	}

	constexpr auto Mas_getter
		= [](MHD& state) -> typename Mass_Density::data_type& {
			return state[Mass_Density()];
		};
	constexpr auto Mom_getter
		= [](MHD& state) -> typename Momentum_Density::data_type& {
			return state[Momentum_Density()];
		};
	constexpr auto Nrj_getter
		= [](MHD& state) -> typename Total_Energy_Density::data_type& {
			return state[Total_Energy_Density()];
		};
	constexpr auto Mag_getter
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

	// compute L/R fluxes along the lines bm/bp: F_{L}-S_{L}U_{L}; F_{R}-S_{R}U_{R}
	flux_neg -= state_neg * bm;
	flux_pos -= state_pos * bp;
	flux_pos[Mag][0] =
	flux_neg[Mag][0] = 0;

	MHD flux
		= (flux_neg + flux_pos) / 2
		+ (flux_neg - flux_pos) * (bp + bm) / (bp - bm) / 2.0;

	return std::make_tuple(flux, std::max(std::fabs(bp), std::fabs(bm)));
}


}}} // namespaces

#endif // ifndef PAMHD_MHD_HLL_ATHENA_HPP
