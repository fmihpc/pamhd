/*
MHD parts of solar wind box program.

Copyright 2024 Finnish Meteorological Institute
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

#ifndef PAMHD_MHD_SW_BOX_HPP
#define PAMHD_MHD_SW_BOX_HPP


#include "algorithm"
#include "cmath"
#include "iterator"
#include "limits"
#include "map"
#include "set"
#include "string"
#include "utility"
#include "vector"

#include "dccrg.hpp"
#include "rapidjson/document.h"

#include "grid/amr.hpp"
#include "mhd/common.hpp"
#include "mhd/solve_staggered.hpp"
#include "mhd/variables.hpp"
#include "variables.hpp"


namespace pamhd {
namespace mhd {


std::tuple<
	double, double, // mass (kg), pressure (Pa)
	std::array<double, 3>, // velocity (m/s)
	std::array<double, 3> // magnetic field (T)
> get_solar_wind_parameters(
	const rapidjson::Document& json,
	const double& /*sim_time*/,
	const double& proton_mass
) {
	using std::string;
	using std::invalid_argument;
	using std::to_string;

	if (not json.HasMember("number-density")) {
		throw invalid_argument(
			string(__FILE__ "(") + to_string(__LINE__) + "): "
			+ "JSON data doesn't have a number-density key."
		);
	}
	const auto& mass_json = json["number-density"];
	if (not mass_json.IsNumber()) {
		throw invalid_argument(
			string(__FILE__ "(") + to_string(__LINE__) + "): "
			+ "JSON item number-density is not a number."
		);
	}
	const double mass = proton_mass * json["number-density"].GetDouble();

	if (not json.HasMember("pressure")) {
		throw invalid_argument(
			string(__FILE__ "(") + to_string(__LINE__) + "): "
			+ "JSON data doesn't have a pressure key."
		);
	}
	const auto& pressure_json = json["pressure"];
	if (not pressure_json.IsNumber()) {
		throw invalid_argument(
			string(__FILE__ "(") + to_string(__LINE__) + "): "
			+ "JSON item pressure is not a number."
		);
	}
	const double pressure = json["pressure"].GetDouble();

	if (not json.HasMember("velocity")) {
		throw invalid_argument(
			string(__FILE__ "(") + to_string(__LINE__) + "): "
			+ "JSON data doesn't have a velocity key."
		);
	}
	const auto& velocity_json = json["velocity"];
	if (not velocity_json.IsArray()) {
		throw invalid_argument(
			string(__FILE__ "(") + to_string(__LINE__) + "): "
			+ "JSON item velocity is not an array."
		);
	}
	const auto& velocity_array = velocity_json.GetArray();
	if (velocity_array.Size() != 3) {
		throw invalid_argument(
			string(__FILE__ "(") + to_string(__LINE__) + "): "
			+ "Invalid size for velocity, should be 3"
		);
	}
	const std::array<double, 3> velocity{
		velocity_array[0].GetDouble(),
		velocity_array[1].GetDouble(),
		velocity_array[2].GetDouble()
	};

	if (not json.HasMember("magnetic-field")) {
		throw invalid_argument(
			string(__FILE__ "(") + to_string(__LINE__) + "): "
			+ "JSON data doesn't have a magnetic field key."
		);
	}
	const auto& magnetic_field_json = json["magnetic-field"];
	if (not magnetic_field_json.IsArray()) {
		throw invalid_argument(
			string(__FILE__ "(") + to_string(__LINE__) + "): "
			+ "JSON item magnetic field is not an array."
		);
	}
	const auto& magnetic_field_array = magnetic_field_json.GetArray();
	if (magnetic_field_array.Size() != 3) {
		throw invalid_argument(
			string(__FILE__ "(") + to_string(__LINE__) + "): "
			+ "Invalid size for magnetic field, should be 3"
		);
	}
	const std::array<double, 3> magnetic_field{
		magnetic_field_array[0].GetDouble(),
		magnetic_field_array[1].GetDouble(),
		magnetic_field_array[2].GetDouble()
	};

	return std::make_tuple(mass, pressure, velocity, magnetic_field);
}


template<
	class Grid,
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Volume_Magnetic_Field_Getter,
	class Face_Magnetic_Field_Getter,
	class Face_Difference_Magnetic_Field_Getter,
	class Background_Magnetic_Field_Getter
> void initialize_plasma(
	Grid& grid,
	const double& sim_time,
	const rapidjson::Document& json,
	const pamhd::Background_Magnetic_Field<
		double,
		pamhd::Magnetic_Field::data_type
	>& background_B,
	const Mass_Density_Getter& Mas,
	const Momentum_Density_Getter& Mom,
	const Total_Energy_Density_Getter& Nrj,
	const Volume_Magnetic_Field_Getter& Vol_B,
	const Face_Magnetic_Field_Getter& Face_B,
	const Face_Difference_Magnetic_Field_Getter& Face_dB,
	const Background_Magnetic_Field_Getter& Bg_B,
	const double& adiabatic_index,
	const double& vacuum_permeability,
	const double& proton_mass
) try {
	using std::string;
	using std::invalid_argument;
	using std::to_string;

	const auto [
		mass, pressure, velocity, magnetic_field
	] = get_solar_wind_parameters(
		json, sim_time, proton_mass
	);
	if (grid.get_rank() == 0) {
		std::cout << "Initializing run, solar wind: "
		<< mass/proton_mass << " #/m^3, "
		<< velocity << " m/s, "
		<< pressure << " Pa, "
		<< magnetic_field << " T"
		<< std::endl;
	}

	for (const auto& cell: grid.local_cells()) {
		const auto [rx, ry, rz] = grid.geometry.get_center(cell.id);
		const auto [sx, sy, sz] = grid.geometry.get_min(cell.id);
		const auto [ex, ey, ez] = grid.geometry.get_max(cell.id);

		Bg_B(*cell.data)(-1) = background_B.get_background_field(
			{sx, ry, rz},
			vacuum_permeability
		);
		Bg_B(*cell.data)(+1) = background_B.get_background_field(
			{ex, ry, rz},
			vacuum_permeability
		);
		Bg_B(*cell.data)(-2) = background_B.get_background_field(
			{rx, sy, rz},
			vacuum_permeability
		);
		Bg_B(*cell.data)(+2) = background_B.get_background_field(
			{rx, ey, rz},
			vacuum_permeability
		);
		Bg_B(*cell.data)(-3) = background_B.get_background_field(
			{rx, ry, sz},
			vacuum_permeability
		);
		Bg_B(*cell.data)(+3) = background_B.get_background_field(
			{rx, ry, ez},
			vacuum_permeability
		);

		Mas(*cell.data) = mass;
		Mom(*cell.data)   =
		Vol_B(*cell.data) = {0, 0, 0};
		Face_B(*cell.data)  =
		Face_dB(*cell.data) = {0, 0, 0, 0, 0, 0};
		Nrj(*cell.data) = pamhd::mhd::get_total_energy_density(
			Mas(*cell.data),
			Mom(*cell.data),
			pressure,
			Vol_B(*cell.data),
			adiabatic_index,
			vacuum_permeability
		);
	}
} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


}} // namespaces


#endif // ifndef PAMHD_MHD_SW_BOX_HPP
