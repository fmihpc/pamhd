/*
Handles options of MHD part of PAMHD.

Copyright 2016, 2017 Ilja Honkonen
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

#ifndef PAMHD_MHD_OPTIONS_HPP
#define PAMHD_MHD_OPTIONS_HPP


#include "cmath"
#include "string"

#include "rapidjson/document.h"


namespace pamhd {
namespace mhd {


struct Options
{
	Options() = default;
	Options(const Options& other) = default;
	Options(Options&& other) = default;

	Options(const rapidjson::Value& object)
	{
		this->set(object);
	};


	std::string solver = "roe_athena";
	double
		save_n = -1,
		min_pressure = 0,
		time_step_factor = 0.5;

	void set(const rapidjson::Value& object) {
		using std::isnormal;

		if (not object.HasMember("save-mhd-n")) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "JSON data doesn't have a save-mhd-n key."
			);
		}
		const auto& save_n_json = object["save-mhd-n"];
		if (not save_n_json.IsNumber()) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "JSON item save-mhd-n is not a number."
			);
		}
		save_n = object["save-mhd-n"].GetDouble();

		if (not object.HasMember("solver-mhd")) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "JSON data doesn't have a solver-mhd key."
			);
		}
		solver = object["solver-mhd"].GetString();
		if (
			solver != "rusanov"
			and solver != "rusanov-staggered"
			and solver != "hll-athena"
			and solver != "hlld-athena"
			and solver != "roe-athena"
		) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Invalid mhd solver: " + solver
				+ ", should be one of rusanov, hll-athena, hlld-athena, roe-athena."
			);
		}

		if (not object.HasMember("minimum-pressure")) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "JSON data doesn't have a minimum-pressure key."
			);
		}
		const auto& min_pressure_json = object["minimum-pressure"];
		if (not min_pressure_json.IsNumber()) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "JSON item minimum-pressure is not a number."
			);
		}
		min_pressure = min_pressure_json.GetDouble();

		if (not object.HasMember("mhd-time-step-factor")) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "JSON data doesn't have a mhd-time-step-factor key."
			);
		}
		const auto& time_step_factor_json = object["mhd-time-step-factor"];
		if (not time_step_factor_json.IsNumber()) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "JSON item mhd-time-step-factor is not a number."
			);
		}
		time_step_factor = object["mhd-time-step-factor"].GetDouble();
	}
};

}} // namespaces


#endif // ifndef PAMHD_MHD_OPTIONS_HPP
