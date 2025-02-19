/*
Handles options of particle part of PAMHD.

Copyright 2016, 2017 Ilja Honkonen
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

#ifndef PAMHD_PARTICLE_OPTIONS_HPP
#define PAMHD_PARTICLE_OPTIONS_HPP


#include "cstdint"
#include "string"

#include "rapidjson/document.h"


namespace pamhd {
namespace particle {


struct Options
{
	double
		save_n{0},
		gyroperiod_time_step_factor{1},
		flight_time_step_factor{1};
	size_t particles_in_cell{10}, min_particles{0};

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
		using std::isnormal;
		using std::string;
		using std::to_string;

		if (object.HasMember("save-particle-n")) {
			const auto& save_n_json = object["save-particle-n"];
			if (not save_n_json.IsNumber()) {
				throw invalid_argument(
					string(__FILE__ "(") + to_string(__LINE__) + "): "
					+ "JSON item save-particle-n is not a number."
				);
			}
			this->save_n = save_n_json.GetDouble();
		}

		if (object.HasMember("particles-in-cell")) {
			const auto& nr_json = object["particles-in-cell"];
			if (not nr_json.IsUint()) {
				throw invalid_argument(
					string(__FILE__ "(") + to_string(__LINE__) + "): "
					+ "JSON item particles-in-cell is not an unsigned integer."
				);
			}
			this->particles_in_cell = nr_json.GetUint();
			if (this->particles_in_cell <= 0) {
				throw invalid_argument(
					string(__FILE__ "(") + to_string(__LINE__) + "): "
					+ "Invalid particles-in-cell: "
					+ to_string(this->particles_in_cell)
					+ ", should be > 0"
				);
			}
		}

		if (object.HasMember("minimum-particles")) {
			const auto& min_particles_json = object["minimum-particles"];
			if (not min_particles_json.IsNumber()) {
				throw invalid_argument(
					string(__FILE__ "(") + to_string(__LINE__) + "): "
					+ "JSON item minimum-particles is not a number."
				);
			}
			this->min_particles = min_particles_json.GetUint();
		}

		if (object.HasMember("gyroperiod-time-step-factor")) {
			const auto& gyroperiod_json = object["gyroperiod-time-step-factor"];
			if (not gyroperiod_json.IsNumber()) {
				throw invalid_argument(
					string(__FILE__ "(") + to_string(__LINE__) + "): "
					+ "JSON item gyroperiod-time-step-factor is not a number."
				);
			}
			this->gyroperiod_time_step_factor = gyroperiod_json.GetDouble();
		}

		if (object.HasMember("flight-time-step-factor")) {
			const auto& flight_time_json = object["flight-time-step-factor"];
			if (not flight_time_json.IsNumber()) {
				throw invalid_argument(
					string(__FILE__ "(") + to_string(__LINE__) + "): "
					+ "JSON item flight-time-step-factor is not a number."
				);
			}
			this->flight_time_step_factor = flight_time_json.GetDouble();
		}
	}
};

}} // namespaces


#endif // ifndef PAMHD_PARTICLE_OPTIONS_HPP
