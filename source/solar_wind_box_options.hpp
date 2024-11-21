/*
Options for solar wind box test program of PAMHD.

Copyright 2016, 2017 Ilja Honkonen
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

#ifndef PAMHD_SOLAR_WIND_BOX_OPTIONS_HPP
#define PAMHD_SOLAR_WIND_BOX_OPTIONS_HPP


#include "cmath"
#include "limits"
#include "string"

#include "rapidjson/document.h"


namespace pamhd {


struct Solar_Wind_Box_Options
{
	Solar_Wind_Box_Options() = default;
	Solar_Wind_Box_Options(const Solar_Wind_Box_Options& other) = default;
	Solar_Wind_Box_Options(Solar_Wind_Box_Options&& other) = default;

	Solar_Wind_Box_Options(const rapidjson::Value& object)
	{
		this->set(object);
	};

	std::string sw_dir = "";
	double inner_bdy_radius = -1;

	void set(const rapidjson::Value& object) {
		using std::invalid_argument;
		using std::string;
		using std::to_string;

		if (not object.HasMember("solar-wind-dir")) {
			throw invalid_argument(
				string(__FILE__ "(") + to_string(__LINE__) + "): "
				+ "JSON data doesn't have a solar-wind-dir key."
			);
		}
		const auto& sw_dir_json = object["solar-wind-dir"];
		if (not sw_dir_json.IsString()) {
			throw invalid_argument(
				string(__FILE__ "(") + to_string(__LINE__) + "): "
				+ "JSON item solar-wind-dir is not a string."
			);
		}
		this->sw_dir = sw_dir_json.GetString();
		std::set<string> allowed_dirs{"-x", "+x", "-y", "+y", "-z", "+z"};
		if (allowed_dirs.count(this->sw_dir) == 0) {
			throw invalid_argument(
				string(__FILE__ "(") + to_string(__LINE__) + "): "
				+ "JSON item solar-wind-dir has invalid value."
			);
		}

		if (not object.HasMember("inner-bdy-radius")) {
			throw invalid_argument(
				string(__FILE__ "(") + to_string(__LINE__) + "): "
				+ "JSON data doesn't have a inner-bdy-radius key."
			);
		}
		const auto& bdy_json = object["inner-bdy-radius"];
		if (not bdy_json.IsNumber()) {
			throw invalid_argument(
				string(__FILE__ "(") + to_string(__LINE__) + "): "
				+ "JSON item inner-bdy-radius is not a number."
			);
		}
		this->inner_bdy_radius = bdy_json.GetDouble();
	}
};

} // namespace


#endif // ifndef PAMHD_SOLAR_WIND_BOX_OPTIONS_HPP
