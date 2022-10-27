/*
Common MHD functions of PAMHD.

Copyright 2014, 2015, 2016, 2017 Ilja Honkonen
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

#ifndef PAMHD_MHD_COMMON_HPP
#define PAMHD_MHD_COMMON_HPP


#include "cmath"
#include "stdexcept"
#include "sstream"
#include "string"


namespace pamhd {
namespace mhd {


// available MHD solvers
enum class Solver {
	rusanov,
	rusanov_staggered,
	hll_athena,
	hlld_athena,
	roe_athena
};


/*!
Returns pressure.

Throws std::domain_error if given a state with non-positive mass density.

Returns negative pressure if total energy in given state
is smaller than kinetic + magnetic energies.

Vector must be compatible with std::array<double, 3>.

mag must not include background magnetic field.
*/
template <
	class Vector,
	class Scalar
> Scalar get_pressure(
	const Scalar& mass_density,
	const Vector& mom,
	const Scalar& total_energy_density,
	const Vector& mag,
	const Scalar& adiabatic_index,
	const Scalar& vacuum_permeability
) {
	using std::to_string;

	if (mass_density <= 0) {
		throw std::domain_error(
			std::string("Non-positive mass density given to ")
			+ __func__ + std::string(": ") + to_string(mass_density)
		);
	}

	const auto
		kinetic_energy
			= 0.5 / mass_density
			* (
				mom[0] * mom[0]
				+ mom[1] * mom[1]
				+ mom[2] * mom[2]
			),
		magnetic_energy
			= 0.5 / vacuum_permeability
			* (
				mag[0] * mag[0]
				+ mag[1] * mag[1]
				+ mag[2] * mag[2]
			);

	return
		(total_energy_density - kinetic_energy - magnetic_energy)
		* (adiabatic_index - 1);
}


template<
	class Vector,
	class Scalar
> Vector get_velocity(
	Vector& momentum,
	Scalar& mass
) {
	if (mass > 0) {
		return momentum / mass;
	} else {
		return {0, 0, 0};
	}
}


/*!
Returns flux from given state into positive x direction.

Throws std::domain_error if given a state with non-positive mass density.
*/
template <
	class Container,
	class Vector,
	class Scalar,
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Magnetic_Field_Getter
> Container get_flux(
	Container& data,
	const Vector& bg_face_magnetic_field,
	const Scalar& adiabatic_index,
	const Scalar& vacuum_permeability,
	const Mass_Density_Getter Mas,
	const Momentum_Density_Getter Mom,
	const Total_Energy_Density_Getter Nrj,
	const Magnetic_Field_Getter Mag
) {
	using std::isfinite;
	using std::isnormal;
	using std::to_string;

	if (not isnormal(Mas(data)) or Mas(data) <= 0) {
		throw std::domain_error(
			std::string("Invalid mass density given to ")
			+ __func__ + std::string(": ") + to_string(Mas(data))
		);
	}

	if (not isfinite(Mom(data)[0])) {
		throw std::domain_error(
			std::string("Invalid momentum x density given to ")
			+ __func__ + std::string(": ") + to_string(Mom(data)[0])
		);
	}

	const auto mag_tot = Mag(data) + bg_face_magnetic_field;
	const auto
		inv_permeability = 1.0 / vacuum_permeability,
		pressure_B0
			= 0.5 * inv_permeability
			* (
				bg_face_magnetic_field[0] * bg_face_magnetic_field[0]
				+ bg_face_magnetic_field[1] * bg_face_magnetic_field[1]
				+ bg_face_magnetic_field[2] * bg_face_magnetic_field[2]
			),
		pressure_B1
			= 0.5 * inv_permeability
			* (
				Mag(data)[0] * Mag(data)[0]
				+ Mag(data)[1] * Mag(data)[1]
				+ Mag(data)[2] * Mag(data)[2]
			),
		pressure_B_tot
			= 0.5 * inv_permeability
			* (
				mag_tot[0] * mag_tot[0]
				+ mag_tot[1] * mag_tot[1]
				+ mag_tot[2] * mag_tot[2]
			),
		pressure_thermal
			= get_pressure(
				Mas(data),
				Mom(data),
				Nrj(data),
				Mag(data),
				adiabatic_index,
				vacuum_permeability
			);

	const auto velocity = get_velocity(Mom(data), Mas(data));

	Container flux;

	Mas(flux) = Mom(data)[0];
	if (not isfinite(Mas(flux))) {
		throw std::domain_error(
			"Invalid mass density flux: "
			+ to_string(Mas(flux))
		);
	}

	Mom(flux)
		= Mom(data) * velocity[0]
		- inv_permeability * (
			mag_tot[0] * mag_tot
			- bg_face_magnetic_field[0] * bg_face_magnetic_field
		);
	Mom(flux)[0] += pressure_thermal + pressure_B_tot - pressure_B0;

	Nrj(flux)
		= velocity[0] * (Nrj(data) + pressure_thermal + pressure_B1)
		- Mag(data)[0] * velocity.dot(Mag(data)) * inv_permeability
		+ inv_permeability * (
			Mag(data)[1] * (velocity[0] * bg_face_magnetic_field[1] - velocity[1] * bg_face_magnetic_field[0])
			+ Mag(data)[2] * (velocity[0] * bg_face_magnetic_field[2] - velocity[2] * bg_face_magnetic_field[0])
		);

	Mag(flux) = velocity[0] * mag_tot - mag_tot[0] * velocity;
	Mag(flux)[0] = 0;

	return flux;
}


/*!
Returns total energy density.

Throws std::domain_error if given a state with non-positive mass density.

Vector must be compatible with std::array<double, 3>.
*/
template <
	class Vector,
	class Scalar
> Scalar get_total_energy_density(
	const Scalar& mass_density,
	const Vector& vel,
	const Scalar& pressure,
	const Vector& mag,
	const Scalar& adiabatic_index,
	const Scalar& vacuum_permeability
) {
	using std::to_string;

	if (mass_density <= 0) {
		throw std::domain_error(
			std::string("Non-positive mass density given to ")
			+ __func__ + std::string(": ") + to_string(mass_density)
		);
	}

	if (adiabatic_index <= 1) {
		throw std::domain_error(
			__func__
			+ std::string("Adiabatic index must be > 1: ")
			+ to_string(adiabatic_index)
		);
	}

	if (vacuum_permeability <= 0) {
		throw std::domain_error(
			__func__
			+ std::string("Vacuum permeability must be > 0: ")
			+ to_string(vacuum_permeability)
		);
	}

	const auto
		kinetic_energy
			= 0.5 * mass_density
			* (
				vel[0] * vel[0]
				+ vel[1] * vel[1]
				+ vel[2] * vel[2]
			),
		magnetic_energy
			= 0.5 / vacuum_permeability
			* (
				mag[0] * mag[0]
				+ mag[1] * mag[1]
				+ mag[2] * mag[2]
			);

	return
		pressure / (adiabatic_index - 1)
		+ kinetic_energy
		+ magnetic_energy;
}


/*!
Returns speed of sound wave.

Throws std::domain_error if given a state with
non-positive mass density or pressure.

Mom and mag must be compatible with std::array<double, 3>.

mag must not include background magnetic field.
*/
template <
	class Vector,
	class Scalar
> Scalar get_sound_speed(
	const Scalar& mass_density,
	const Vector& mom,
	const Scalar& total_energy_density,
	const Vector& mag,
	const Scalar adiabatic_index,
	const Scalar vacuum_permeability
) {
	using std::isnormal;
	using std::to_string;

	if (not isnormal(mass_density) or mass_density < 0) {
		throw std::domain_error(
			std::string("Invalid mass density in state_pos: ") + __func__ + ": "
			+ to_string(mass_density)
		);
	}

	const auto pressure
		= get_pressure(
			mass_density,
			mom,
			total_energy_density,
			mag,
			adiabatic_index,
			vacuum_permeability
		);
	if (not isnormal(pressure) or pressure < 0) {
		std::ostringstream error;
		error << std::scientific;
		error << "Non-positive pressure given to " << __func__;
		error << ": " << pressure;
		error << ", from mass, momentum, energy and field: ";
		error << mass_density;
		error << ", [" << mom[0] << ", " << mom[1] << ", " << mom[2] << "], ";
		error << total_energy_density;
		error << ", [" << mag[0] << ", " << mag[1] << ", " << mag[2] << "]";
		throw std::domain_error(error.str());
	}

	return std::sqrt(adiabatic_index * pressure / mass_density);
}


/*!
Returns speed of Alfvén wave.

Throws std::domain_error if given a state with non-positive mass density.

Vector must be compatible with std::array<double, 3>.
*/
template <
	class Vector,
	class Scalar
> Scalar get_alfven_speed(
	const Scalar& mass_density,
	const Vector& mag,
	const Vector& bg_mag,
	const Scalar& vacuum_permeability
) {
	using std::isnormal;
	using std::pow;
	using std::sqrt;
	using std::to_string;

	if (not isnormal(mass_density) or mass_density < 0) {
		throw std::domain_error(
			std::string("Non-positive mass density given to ")
			+ __func__ + std::string(": ") + to_string(mass_density)
		);
	}

	const auto mag_mag = sqrt(
		pow(mag[0] + bg_mag[0], 2)
		+ pow(mag[1] + bg_mag[1], 2)
		+ pow(mag[2] + bg_mag[2], 2)
	);
	return mag_mag / sqrt(vacuum_permeability * mass_density);
}


/*!
Returns speed of fast magnetosonic wave in first dimension.

Returns a non-negative value.

Vector must be compatible with std::array<double, 3>.
*/
template <
	class Vector,
	class Scalar
> Scalar get_fast_magnetosonic_speed(
	const Scalar& mass_density,
	const Vector& momentum_density,
	const Scalar& total_energy_density,
	const Vector& mag,
	const Vector& bg_mag,
	const Scalar& adiabatic_index,
	const Scalar& vacuum_permeability
) {
	using std::isfinite;
	using std::isnormal;
	using std::pow;
	using std::sqrt;
	using std::to_string;

	const auto mag_tot = mag + bg_mag;
	if (not isfinite(mag_tot[0])) {
		throw std::domain_error(
			"Invalid total magnetic field x: "
			+ to_string(mag_tot[0])
		);
	}

	const auto
		mag_mag = sqrt(
			pow(mag_tot[0], 2)
			+ pow(mag_tot[1], 2)
			+ pow(mag_tot[2], 2)
		),
		sound = get_sound_speed(
			mass_density,
			momentum_density,
			total_energy_density,
			mag,
			adiabatic_index,
			vacuum_permeability
		);

	if (mag_mag <= 0) {
		return sound;
	}

	const auto
		sound2 = sound * sound,
		alfven2 = pow(
			get_alfven_speed(mass_density, mag, bg_mag, vacuum_permeability),
			2
		),
		speeds_squared = sound2 + alfven2;

	if (not isnormal(speeds_squared) or speeds_squared < 0) {
		throw std::domain_error(
			std::string("Invalid squared speeds in ") + __func__ + ": "
			+ to_string(speeds_squared)
		);
	}

	const auto to_sqrt
		= sound2*sound2
		+ alfven2*alfven2
		+ 2*sound2*alfven2
			* (1 - 2 * pow(mag_tot[0] / mag_mag, 2));

	if (not isnormal(to_sqrt) or to_sqrt < 0) {
		throw std::domain_error(
			std::string("Invalid to_sqrt in ") + __func__ + ": "
			+ to_string(to_sqrt) + " with s2: " + to_string(sound2)
			+ ", a2: " + to_string(alfven2)
			+ ", B[0]/B: " + to_string(mag_tot[0] / mag_mag)
		);
	}

	return sqrt(0.5 * (speeds_squared + sqrt(to_sqrt)));
}


/*!
Returns a flux that is separated into contributions from state_neg and state_pos.

Throws std::domain_error if given a state with non-positive mass density.

Ignores background magnetic field.
*/
template <
	class Mass_Density,
	class Momentum_Density,
	class Total_Energy_Density,
	class Magnetic_Field,
	class MHD,
	class Vector,
	class Scalar
> std::tuple<MHD, MHD> N_get_flux(
	MHD& state_neg,
	MHD& state_pos,
	const Vector& bg_face_magnetic_field,
	const Scalar& adiabatic_index,
	const Scalar& vacuum_permeability
) {
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

	return std::make_tuple(
		get_flux(
			state_neg,
			bg_face_magnetic_field,
			adiabatic_index,
			vacuum_permeability,
			Mas_getter, Mom_getter, Nrj_getter, Mag_getter
		),
		get_flux(
			state_pos,
			bg_face_magnetic_field,
			adiabatic_index,
			vacuum_permeability,
			Mas_getter, Mom_getter, Nrj_getter, Mag_getter
		)
	);
}


/*!
Throws std::domain_error if given a state with non-positive mass density.
*/
template <
	class Primitive,
	class Conservative,
	class Vector,
	class Scalar,
	class P_Mass_Density_Getter,
	class Velocity_Getter,
	class Pressure_Getter,
	class C_Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter
> Primitive get_primitive(
	Conservative data,
	const Vector& magnetic_field,
	const Scalar& adiabatic_index,
	const Scalar& vacuum_permeability,
	const P_Mass_Density_Getter Mas_p,
	const Velocity_Getter Vel,
	const Pressure_Getter Pre,
	const C_Mass_Density_Getter Mas_c,
	const Momentum_Density_Getter Mom,
	const Total_Energy_Density_Getter Nrj
) {
	using std::to_string;

	if (Mas_c(data) <= 0) {
		throw std::domain_error(
			std::string("Non-positive mass density given to ")
			+ __func__ + std::string(": ") + to_string(Mas_c(data))
		);
	}

	Primitive ret_val;

	Mas_p(ret_val) = Mas_c(data);
	Vel(ret_val) = get_velocity(Mom(data), Mas_c(data));
	Pre(ret_val)
		= get_pressure(
			Mas_c(data),
			Mom(data),
			Nrj(data),
			magnetic_field,
			adiabatic_index,
			vacuum_permeability
		);

	return ret_val;
}


/*!
Rotates components of given vector for get_flux_* functions.

Positive direction rotates vector for get_flux_* (assume
two states have identical y and z coordinates), negative
rotates it back.
*/
template<class Vector> Vector get_rotated_vector(
	const Vector& v,
	const int direction
) {
	Vector ret_val;

	switch(direction) {
	case -3:
		ret_val[0] = v[1];
		ret_val[1] = v[2];
		ret_val[2] = v[0];
		return ret_val;
	case -2:
		ret_val[0] = v[2];
		ret_val[1] = v[0];
		ret_val[2] = v[1];
		return ret_val;
	case -1:
	case 1:
		return v;
	case 2:
		ret_val[0] = v[1];
		ret_val[1] = v[2];
		ret_val[2] = v[0];
		return ret_val;
	case 3:
		ret_val[0] = v[2];
		ret_val[1] = v[0];
		ret_val[2] = v[1];
		return ret_val;
	default:
		std::cerr <<  __FILE__ << "(" << __LINE__ << ") "
			<< "Invalid direction: " << direction
			<< std::endl;
		abort();
	}
}


/*!
Function call signature of all MHD flux functions.

Input:
    - Conservative MHD variables in two cells that share
      a face and are neighbors in the x dimension
    - Area shared between given cells
    - Length of time step for which to calculate flux

state_neg represents the MHD variables in the cell in the
negative x direction from the shared face, state_pos in
the cell in positive x direction from the face.

Output:
    - Flux of conservative MHD variables over time dt
      through area shared_area
    - Absolute value of maximum signal speed from shared face

See for example hll_athena.hpp for an implementation.
*/
template <
	class MHD,
	class Vector
> using solver_t = std::function<
	std::tuple<
		MHD,
		double
	>(
		MHD& /* state_neg */,
		MHD& /* state_pos */,
		const Vector& /* bg_face_magnetic_field */,
		const double& /* shared_area */,
		const double& /* dt */,
		const double& /* adiabatic_index */,
		const double& /* vacuum_permeability */
	)
>;


/*!
Positive flux adds to given cell multiplied with given factor.

Throws std::domain_error if new state has non-positive
mass density or pressure check_new_state == true.
*/
template <
	class Container,
	class Scalar,
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Magnetic_Field_Getter,
	class Mass_Density_Flux_Getter,
	class Momentum_Density_Flux_Getter,
	class Total_Energy_Density_Flux_Getter,
	class Magnetic_Field_Flux_Getter
> void apply_fluxes(
	Container& data,
	const Scalar& factor,
	const Scalar& adiabatic_index,
	const Scalar& vacuum_permeability,
	const Mass_Density_Getter Mas,
	const Momentum_Density_Getter Mom,
	const Total_Energy_Density_Getter Nrj,
	const Magnetic_Field_Getter Mag,
	const Mass_Density_Flux_Getter Mas_f,
	const Momentum_Density_Flux_Getter Mom_f,
	const Total_Energy_Density_Flux_Getter Nrj_f,
	const Magnetic_Field_Flux_Getter Mag_f
) {
	using std::to_string;

	Mas(data) += Mas_f(data) * factor;
	Mom(data) += Mom_f(data) * factor;
	Nrj(data) += Nrj_f(data) * factor;
	Mag(data) += Mag_f(data) * factor;

	if (Mas(data) <= 0) {
		throw std::domain_error(
			"New state has negative mass density: " + to_string(Mas(data))
		);
	}

	const auto pressure
		= get_pressure(
			Mas(data),
			Mom(data),
			Nrj(data),
			Mag(data),
			adiabatic_index,
			vacuum_permeability
		);
	if (pressure <= 0) {
		throw std::domain_error(
			"New state has negative pressure: " + to_string(pressure)
		);
	}
}

template <
	class Container,
	class Scalar,
	class Mass_Density_Getters,
	class Momentum_Density_Getters,
	class Total_Energy_Density_Getters,
	class Magnetic_Field_Getter,
	class Mass_Density_Flux_Getters,
	class Momentum_Density_Flux_Getters,
	class Total_Energy_Density_Flux_Getters,
	class Magnetic_Field_Flux_Getter
> void apply_fluxes_N(
	Container& data,
	const Scalar& factor,
	const Mass_Density_Getters Mas,
	const Momentum_Density_Getters Mom,
	const Total_Energy_Density_Getters Nrj,
	const Magnetic_Field_Getter Mag,
	const Mass_Density_Flux_Getters Mas_f,
	const Momentum_Density_Flux_Getters Mom_f,
	const Total_Energy_Density_Flux_Getters Nrj_f,
	const Magnetic_Field_Flux_Getter Mag_f
) {
	Mag(data) += Mag_f(data) * factor;

	Mas.first(data) += Mas_f.first(data) * factor;
	Mom.first(data) += Mom_f.first(data) * factor;
	Nrj.first(data) += Nrj_f.first(data) * factor;

	Mas.second(data) += Mas_f.second(data) * factor;
	Mom.second(data) += Mom_f.second(data) * factor;
	Nrj.second(data) += Nrj_f.second(data) * factor;
}


}} // namespaces

#endif // ifndef PAMHD_MHD_COMMON_HPP
