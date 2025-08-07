/*
Handles options of MHD part of PAMHD.

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

#ifndef PAMHD_MHD_OPTIONS_HPP
#define PAMHD_MHD_OPTIONS_HPP


#include "cmath"
#include "limits"
#include "string"

#include "rapidjson/document.h"

#include "math/expression.hpp"
#include "mhd/variables.hpp"


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
		time_step_factor = 0.5,
		number_density_mrl_at = 9e99,
		number_density_min_mrg = std::numeric_limits<double>::epsilon(),
		pressure_mrl_at = 9e99,
		pressure_min_mrg = std::numeric_limits<double>::epsilon(),
		vel_mrl_at = 9e99,
		vel_min_mrg = std::numeric_limits<double>::epsilon(),
		mag_mrl_at = 9e99,
		mag_min_mrg = std::numeric_limits<double>::epsilon();

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
			and solver != "hybrid"
		) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Invalid mhd solver: " + solver
				+ ", should be one of rusanov, hll-athena, hlld-athena, roe-athena, hybrid."
			);
		}

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

		if (object.HasMember("number-density-max-ref-lvl-at")) {
			const auto& number_density_mrl_at_json = object["number-density-max-ref-lvl-at"];
			if (not number_density_mrl_at_json.IsNumber()) {
				throw std::invalid_argument(
					std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
					+ "JSON item number-density-max-ref-lvl-at is not a number."
				);
			}
			number_density_mrl_at = number_density_mrl_at_json.GetDouble();
		}
		if (object.HasMember("number-density-min-mrg")) {
			const auto& number_density_min_mrg_json = object["number-density-min-mrg"];
			if (not number_density_min_mrg_json.IsNumber()) {
				throw std::invalid_argument(
					std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
					+ "JSON item number-density-min-mrg is not a number."
				);
			}
			number_density_min_mrg = number_density_min_mrg_json.GetDouble();
		}

		if (object.HasMember("pressure-max-ref-lvl-at")) {
			const auto& pressure_mrl_at_json = object["pressure-max-ref-lvl-at"];
			if (not pressure_mrl_at_json.IsNumber()) {
				throw std::invalid_argument(
					std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
					+ "JSON item pressure-max-ref-lvl-at is not a number."
				);
			}
			pressure_mrl_at = pressure_mrl_at_json.GetDouble();
		}
		if (object.HasMember("pressure-min-mrg")) {
			const auto& pressure_min_mrg_json = object["pressure-min-mrg"];
			if (not pressure_min_mrg_json.IsNumber()) {
				throw std::invalid_argument(
					std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
					+ "JSON item pressure-min-mrg is not a number."
				);
			}
			pressure_min_mrg = pressure_min_mrg_json.GetDouble();
		}

		if (object.HasMember("magnetic-field-max-ref-lvl-at")) {
			const auto& mag_mrl_at_json = object["magnetic-field-max-ref-lvl-at"];
			if (not mag_mrl_at_json.IsNumber()) {
				throw std::invalid_argument(
					std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
					+ "JSON item magnetic-field-max-ref-lvl-at is not a number."
				);
			}
			mag_mrl_at = mag_mrl_at_json.GetDouble();
		}
		if (object.HasMember("magnetic-field-min-mrg")) {
			const auto& mag_min_mrg_json = object["magnetic-field-min-mrg"];
			if (not mag_min_mrg_json.IsNumber()) {
				throw std::invalid_argument(
					std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
					+ "JSON item magnetic-field-min-mrg is not a number."
				);
			}
			mag_min_mrg = mag_min_mrg_json.GetDouble();
		}

		if (object.HasMember("velocity-max-ref-lvl-at")) {
			const auto& vel_mrl_at_json = object["velocity-max-ref-lvl-at"];
			if (not vel_mrl_at_json.IsNumber()) {
				throw std::invalid_argument(
					std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
					+ "JSON item velocity-max-ref-lvl-at is not a number."
				);
			}
			vel_mrl_at = vel_mrl_at_json.GetDouble();
		}
		if (object.HasMember("velocity-min-mrg")) {
			const auto& vel_min_mrg_json = object["velocity-min-mrg"];
			if (not vel_min_mrg_json.IsNumber()) {
				throw std::invalid_argument(
					std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
					+ "JSON item velocity-min-mrg is not a number."
				);
			}
			vel_min_mrg = vel_min_mrg_json.GetDouble();
		}
	}
};

}} // namespaces


#endif // ifndef PAMHD_MHD_OPTIONS_HPP
