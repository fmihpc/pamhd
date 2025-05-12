/*
Handles options common to MHD, particle and PAMHD test programs.

Copyright 2017 Ilja Honkonen
Copyright 2018, 2025 Finnish Meteorological Institute
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

#ifndef PAMHD_SIMULATION_OPTIONS_HPP
#define PAMHD_SIMULATION_OPTIONS_HPP


#include "stdexcept"
#include "string"

#include "rapidjson/document.h"

#include "common_variables.hpp"


namespace pamhd {


struct Options {
	std::string lb_name{"RCB"}, output_directory{""};
	int
		substep_min_i = 0,
		substep_max_i = 999;
	double
		time_start{0}, time_length{1},
		adiabatic_index{5.0 / 3.0},
		vacuum_permeability{4e-7 * M_PI},
		proton_mass{1.672621777e-27},
		charge2mass{95788332}, // charge to mass ratio (C/kg)
		temp2nrj{1.380649e-23}; // Boltzmann constant (J/K)
	math::Expression<Substep_Min> substep_min_e;
	math::Expression<Substep_Max> substep_max_e;

	Options() = default;
	Options(const Options& other) = default;
	Options& operator=(const Options& other) = default;
	Options(Options&& other) = delete;
	Options& operator=(Options&& other) = delete;

	Options(const rapidjson::Value& object)
	{
		this->set(object);
	};

	void set(const rapidjson::Value& object) {
		using std::invalid_argument;
		using std::isfinite;
		using std::isnormal;
		using std::string;
		using std::to_string;

		if (object.HasMember("time-start")) {
			const auto& time_start_json = object["time-start"];
			if (not time_start_json.IsNumber()) {
				throw invalid_argument(
					string(__FILE__ "(") + to_string(__LINE__) + "): "
					+ "JSON item time-start is not a number."
				);
			}
			this->time_start = object["time-start"].GetDouble();
			if (not isfinite(this->time_start)) {
				throw invalid_argument(
					string(__FILE__ "(") + to_string(__LINE__) + "): "
					+ "Invalid time-start: "
					+ to_string(this->time_start)
				);
			}
		}

		if (object.HasMember("time-length")) {
			const auto& time_length_json = object["time-length"];
			if (not time_length_json.IsNumber()) {
				throw invalid_argument(
					string(__FILE__ "(") + to_string(__LINE__) + "): "
					+ "JSON item time-length is not a number."
				);
			}
			this->time_length = object["time-length"].GetDouble();
			if (
				not isnormal(this->time_length)
				or this->time_length < 0
			) {
				throw invalid_argument(
					string(__FILE__ "(") + to_string(__LINE__) + "): "
					+ "Invalid time-length: "
					+ to_string(this->time_length)
					+ ", should be > 0"
				);
			}
		}

		if (object.HasMember("adiabatic-index")) {
			const auto& adiabatic_index_json = object["adiabatic-index"];
			if (not adiabatic_index_json.IsNumber()) {
				throw invalid_argument(
					string(__FILE__ "(") + to_string(__LINE__) + "): "
					+ "JSON item adiabatic-index is not a number."
				);
			}
			this->adiabatic_index = object["adiabatic-index"].GetDouble();
			if (
				not isnormal(this->adiabatic_index)
				or this->adiabatic_index < 0
			) {
				throw invalid_argument(
					string(__FILE__ "(") + to_string(__LINE__) + "): "
					+ "Invalid adiabatic index: "
					+ to_string(this->adiabatic_index)
					+ ", should be > 0"
				);
			}
		}

		if (object.HasMember("vacuum-permeability")) {
			const auto& vacuum_permeability_json = object["vacuum-permeability"];
			if (not vacuum_permeability_json.IsNumber()) {
				throw invalid_argument(
					string(__FILE__ "(") + to_string(__LINE__) + "): "
					+ "JSON item vacuum-permeability is not a number."
				);
			}
			this->vacuum_permeability = object["vacuum-permeability"].GetDouble();
			if (
				not isnormal(this->vacuum_permeability)
				or this->vacuum_permeability < 0
			) {
				throw invalid_argument(
					string(__FILE__ "(") + to_string(__LINE__) + "): "
					+ "Invalid vacuum permeability: "
					+ to_string(this->vacuum_permeability)
					+ ", should be > 0"
				);
			}
		}

		if (object.HasMember("proton-mass")) {
			const auto& proton_mass_json = object["proton-mass"];
			if (not proton_mass_json.IsNumber()) {
				throw invalid_argument(
					string(__FILE__ "(") + to_string(__LINE__) + "): "
					+ "JSON item proton-mass is not a number."
				);
			}
			this->proton_mass = object["proton-mass"].GetDouble();
			if (
				not isnormal(this->proton_mass)
				or this->proton_mass < 0
			) {
				throw invalid_argument(
					string(__FILE__ "(") + to_string(__LINE__) + "): "
					+ "Invalid proton_mass: "
					+ to_string(this->proton_mass)
					+ ", should be > 0"
				);
			}
		}

		if (object.HasMember("charge-mass-ratio")) {
			const auto& c2m_json = object["charge-mass-ratio"];
			if (not c2m_json.IsNumber()) {
				throw invalid_argument(
					string(__FILE__ "(") + to_string(__LINE__) + "): "
					+ "JSON item charge-mass-ratio is not a number."
				);
			}
			this->charge2mass = object["charge-mass-ratio"].GetDouble();
			if (
				not isnormal(this->charge2mass)
				or this->charge2mass < 0
			) {
				throw invalid_argument(
					string(__FILE__ "(") + to_string(__LINE__) + "): "
					+ "Invalid charge-mass-ratio: "
					+ to_string(this->charge2mass)
					+ ", should be > 0"
				);
			}
		}

		if (object.HasMember("particle-temp-nrj-ratio")) {
			const auto& temp2nrj_json = object["particle-temp-nrj-ratio"];
			if (not temp2nrj_json.IsNumber()) {
				throw invalid_argument(
					string(__FILE__ "(") + to_string(__LINE__) + "): "
					+ "JSON item particle-temp-nrj-ratio is not a number."
				);
			}
			this->temp2nrj = object["particle-temp-nrj-ratio"].GetDouble();
			if (
				not isnormal(this->temp2nrj)
				or this->temp2nrj < 0
			) {
				throw invalid_argument(
					string(__FILE__ "(") + to_string(__LINE__) + "): "
					+ "Invalid particle-temp-nrj-ratio: "
					+ to_string(this->temp2nrj)
					+ ", should be > 0"
				);
			}
		}

		if (object.HasMember("output-directory")) {
			output_directory = object["output-directory"].GetString();
		}

		if (object.HasMember("load-balancer")) {
			this->lb_name = object["load-balancer"].GetString();
		}

		if (object.HasMember("substep-min")) {
			const auto& substep_min_json = object["substep-min"];
			if (substep_min_json.IsInt()) {
				this->substep_min_i = substep_min_json.GetInt();
				if (this->substep_min_i < 0) {
					throw std::invalid_argument(
						__FILE__ "(" + std::to_string(__LINE__) + "): "
						+ "JSON item substep-min cannot be negative."
					);
				}
			} else if (substep_min_json.IsString()) {
				this->substep_min_i = -2;
				this->substep_min_e.set_expression(substep_min_json.GetString());
			} else {
				throw std::invalid_argument(
					__FILE__ "(" + std::to_string(__LINE__) + "): "
					+ "JSON item substep-min must either non-negative integer or string."
				);
			}
		}

		if (object.HasMember("substep-max")) {
			const auto& substep_max_json = object["substep-max"];
			if (substep_max_json.IsInt()) {
				this->substep_max_i = substep_max_json.GetInt();
				if (this->substep_max_i < 0) {
					throw std::invalid_argument(
						__FILE__ "(" + std::to_string(__LINE__) + "): "
						+ "JSON item substep-max cannot be negative."
					);
				}
			} else if (substep_max_json.IsString()) {
				this->substep_max_i = -2;
				this->substep_max_e.set_expression(substep_max_json.GetString());
			} else {
				throw std::invalid_argument(
					__FILE__ "(" + std::to_string(__LINE__) + "): "
					+ "JSON item substep-max must either non-negative integer or string."
				);
			}
		}
	}
};

} // namespace


#endif // ifndef PAMHD_SIMULATION_OPTIONS_HPP
