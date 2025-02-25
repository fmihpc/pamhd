/*
Options for solar wind box test programs of PAMHD.

Copyright 2016, 2017 Ilja Honkonen
Copyright 2024, 2025 Finnish Meteorological Institute
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

#ifndef PAMHD_SOLAR_WIND_BOX_OPTIONS_HPP
#define PAMHD_SOLAR_WIND_BOX_OPTIONS_HPP


#include "array"
#include "stdexcept"
#include "string"

#include "rapidjson/document.h"


namespace pamhd {


struct Solar_Wind_Box_Options
{
	std::string sw_dir_s{""};
	int sw_dir{0};
	size_t sw_dim{3};
	double
		inner_radius{-1},
		inner_nr_density{-1},
		inner_pressure{-1},
		sw_nr_density{-1},
		sw_pressure{-1};
	std::array<double, 3>
		inner_velocity{0, 0, 0},
		sw_velocity{0, 0, 0},
		sw_magnetic_field{0, 0, 0};


	Solar_Wind_Box_Options() = default;
	Solar_Wind_Box_Options(const Solar_Wind_Box_Options& other) = default;
	Solar_Wind_Box_Options(Solar_Wind_Box_Options&& other) = delete;
	Solar_Wind_Box_Options& operator=(Solar_Wind_Box_Options&& other) = delete;

	Solar_Wind_Box_Options(const rapidjson::Value& object)
	{
		this->set(object);
	};

	void set(const rapidjson::Value& object) {
		using std::invalid_argument;
		using std::string;
		using std::to_string;

		if (not object.HasMember("solar-wind")) {
			throw invalid_argument(
				string(__FILE__ "(") + to_string(__LINE__) + "): "
				+ "JSON data doesn't have a 'solar-wind' key."
			);
		}
		const auto& sw = object["solar-wind"];
		if (not sw.IsObject()) {
			throw invalid_argument(
				string(__FILE__ "(") + to_string(__LINE__) + "): "
				+ "JSON item 'solar-wind' isn't an object."
			);
		}

		if (not sw.HasMember("direction")) {
			throw invalid_argument(
				string(__FILE__ "(") + to_string(__LINE__) + "): "
				+ "'solar-wind' data doesn't have a 'direction' key."
			);
		}
		const auto& sw_dir_json = sw["direction"];
		if (not sw_dir_json.IsString()) {
			throw invalid_argument(
				string(__FILE__ "(") + to_string(__LINE__) + "): "
				+ "'solar-wind' item 'direction' is not a string."
			);
		}
		this->sw_dir_s = sw_dir_json.GetString();
		if (this->sw_dir_s == "-x") {
			this->sw_dir = -1;
		} else if (this->sw_dir_s == "+x") {
			this->sw_dir = +1;
		} else if (this->sw_dir_s == "-y") {
			this->sw_dir = -2;
		} else if (this->sw_dir_s == "+y") {
			this->sw_dir = +2;
		} else if (this->sw_dir_s == "-z") {
			this->sw_dir = -3;
		} else if (this->sw_dir_s == "+z") {
			this->sw_dir = +3;
		} else {
			throw invalid_argument(
				string(__FILE__ "(") + to_string(__LINE__) + "): "
				+ "'solar-wind' item 'direction' must be one of: -x,+x,-y,+y,-z,+z."
			);
		}
		this->sw_dim = std::abs(this->sw_dir) - 1;

		if (not sw.HasMember("number-density")) {
			throw invalid_argument(
				string(__FILE__ "(") + to_string(__LINE__) + "): "
				+ "'solar-wind' doesn't have a 'number-density' key."
			);
		}
		const auto& sw_nr_json = sw["number-density"];
		if (not sw_nr_json.IsNumber()) {
			throw invalid_argument(
				string(__FILE__ "(") + to_string(__LINE__) + "): "
				+ "'solar-wind' item 'number-density' is not a number."
			);
		}
		this->sw_nr_density = sw_nr_json.GetDouble();

		if (not sw.HasMember("pressure")) {
			throw invalid_argument(
				string(__FILE__ "(") + to_string(__LINE__) + "): "
				+ "'solar-wind' doesn't have a 'pressure' key."
			);
		}
		const auto& sw_p_json = sw["pressure"];
		if (not sw_p_json.IsNumber()) {
			throw invalid_argument(
				string(__FILE__ "(") + to_string(__LINE__) + "): "
				+ "'solar-wind' item 'pressure' is not a number."
			);
		}
		this->sw_pressure = sw_p_json.GetDouble();

		if (not sw.HasMember("velocity")) {
			throw invalid_argument(
				string(__FILE__ "(") + to_string(__LINE__) + "): "
				+ "'solar-wind' doesn't have a 'velocity' key."
			);
		}
		const auto& sw_v_json = sw["velocity"];
		if (not sw_v_json.IsArray()) {
			throw invalid_argument(
				string(__FILE__ "(") + to_string(__LINE__) + "): "
				+ "'solar-wind' item 'velocity' is not an array."
			);
		}
		const auto& sw_v_array = sw_v_json.GetArray();
		if (sw_v_array.Size() != 3) {
			throw invalid_argument(
				string(__FILE__ "(") + to_string(__LINE__) + "): "
				+ "Invalid number of solar wind velocity components, should be 3."
			);
		}
		this->sw_velocity = {
			sw_v_array[0].GetDouble(),
			sw_v_array[1].GetDouble(),
			sw_v_array[2].GetDouble()
		};

		if (not sw.HasMember("magnetic-field")) {
			throw invalid_argument(
				string(__FILE__ "(") + to_string(__LINE__) + "): "
				+ "'solar-wind' doesn't have a 'velocity' key."
			);
		}
		const auto& sw_B_json = sw["magnetic-field"];
		if (not sw_B_json.IsArray()) {
			throw invalid_argument(
				string(__FILE__ "(") + to_string(__LINE__) + "): "
				+ "'solar-wind' item 'magnetic-field' is not an array."
			);
		}
		const auto& sw_B_array = sw_B_json.GetArray();
		if (sw_B_array.Size() != 3) {
			throw invalid_argument(
				string(__FILE__ "(") + to_string(__LINE__) + "): "
				+ "Invalid number of solar wind magnetic field components, should be 3."
			);
		}
		this->sw_magnetic_field = {
			sw_B_array[0].GetDouble(),
			sw_B_array[1].GetDouble(),
			sw_B_array[2].GetDouble()
		};

		if (not object.HasMember("inner-boundary")) {
			throw invalid_argument(
				string(__FILE__ "(") + to_string(__LINE__) + "): "
				+ "JSON data doesn't have an 'inner-boundary' key."
			);
		}
		const auto& inner = object["inner-boundary"];
		if (not inner.IsObject()) {
			throw invalid_argument(
				string(__FILE__ "(") + to_string(__LINE__) + "): "
				+ "JSON item 'inner-boundary' isn't an object."
			);
		}

		if (not inner.HasMember("radius")) {
			throw invalid_argument(
				string(__FILE__ "(") + to_string(__LINE__) + "): "
				+ "'inner-boundary' doesn't have a 'radius' key."
			);
		}
		const auto& r_json = inner["radius"];
		if (not r_json.IsNumber()) {
			throw invalid_argument(
				string(__FILE__ "(") + to_string(__LINE__) + "): "
				+ "'inner-boundary' item 'radius' is not a number."
			);
		}
		this->inner_radius = r_json.GetDouble();

		if (not inner.HasMember("number-density")) {
			throw invalid_argument(
				string(__FILE__ "(") + to_string(__LINE__) + "): "
				+ "'inner-boundary' doesn't have a 'number-density' key."
			);
		}
		const auto& n_json = inner["number-density"];
		if (not n_json.IsNumber()) {
			throw invalid_argument(
				string(__FILE__ "(") + to_string(__LINE__) + "): "
				+ "'inner-boundary' item 'number-density' is not a number."
			);
		}
		this->inner_nr_density = n_json.GetDouble();

		if (not inner.HasMember("pressure")) {
			throw invalid_argument(
				string(__FILE__ "(") + to_string(__LINE__) + "): "
				+ "'inner-boundary' doesn't have a 'pressure' key."
			);
		}
		const auto& p_json = inner["pressure"];
		if (not p_json.IsNumber()) {
			throw invalid_argument(
				string(__FILE__ "(") + to_string(__LINE__) + "): "
				+ "'inner-boundary' item 'pressure' is not a number."
			);
		}
		this->inner_pressure = p_json.GetDouble();

		if (not inner.HasMember("velocity")) {
			throw invalid_argument(
				string(__FILE__ "(") + to_string(__LINE__) + "): "
				+ "'inner-boundary' doesn't have a 'velocity' key."
			);
		}
		const auto& v_json = inner["velocity"];
		if (not v_json.IsArray()) {
			throw invalid_argument(
				string(__FILE__ "(") + to_string(__LINE__) + "): "
				+ "'inner-boundary' item 'velocity' is not an array."
			);
		}
		this->inner_velocity = {
			v_json[0].GetDouble(),
			v_json[1].GetDouble(),
			v_json[2].GetDouble()
		};
	}
};

} // namespace


#endif // ifndef PAMHD_SOLAR_WIND_BOX_OPTIONS_HPP
