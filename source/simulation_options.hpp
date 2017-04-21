/*
Handles options common to MHD, particle and PAMHD test programs.

Copyright 2017 Ilja Honkonen
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

#ifndef PAMHD_SIMULATION_OPTIONS_HPP
#define PAMHD_SIMULATION_OPTIONS_HPP


#include "cmath"
#include "string"

#include "rapidjson/document.h"


namespace pamhd {


struct Options
{
	Options() = default;
	Options(const Options& other) = default;
	Options(Options&& other) = default;

	Options(const rapidjson::Value& object)
	{
		this->set(object);
	};


	std::string
		lb_name = "RCB",
		output_directory = "";

	double
		time_start = 0,
		time_length = 1,
		time_step_factor = 0.5,
		adiabatic_index = 5.0 / 3.0,
		vacuum_permeability = 4e-7 * M_PI,
		proton_mass = 1.672621777e-27,
		boltzmann = 1.38064852e-23;

	void set(const rapidjson::Value& object) {
		using std::isnormal;

		if (not object.HasMember("time-start")) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "JSON data doesn't have a time-start key."
			);
		}
		time_start = object["time-start"].GetDouble();

		if (not object.HasMember("time-length")) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "JSON data doesn't have a time-length key."
			);
		}
		time_length = object["time-length"].GetDouble();

		if (not object.HasMember("time-step-factor")) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "JSON data doesn't have a time-step-factor key."
			);
		}
		time_step_factor = object["time-step-factor"].GetDouble();

		if (not object.HasMember("adiabatic-index")) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "JSON data doesn't have a adiabatic-index key."
			);
		}
		adiabatic_index = object["adiabatic-index"].GetDouble();
		if (not isnormal(adiabatic_index) or adiabatic_index < 0) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Invalid adiabatic index: " + std::to_string(adiabatic_index)
				+ ", should be > 0"
			);
		}

		if (not object.HasMember("vacuum-permeability")) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "JSON data doesn't have a vacuum-permeability key."
			);
		}
		vacuum_permeability = object["vacuum-permeability"].GetDouble();
		if (not isnormal(vacuum_permeability) or vacuum_permeability < 0) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Invalid vacuum permeability: " + std::to_string(vacuum_permeability)
				+ ", should be > 0"
			);
		}

		if (not object.HasMember("proton-mass")) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "JSON data doesn't have a proton-mass key."
			);
		}
		proton_mass = object["proton-mass"].GetDouble();
		if (not isnormal(proton_mass) or proton_mass < 0) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Invalid proton_mass: " + std::to_string(proton_mass)
				+ ", should be > 0"
			);
		}

		if (not object.HasMember("particle-temp-nrj-ratio")) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "JSON data doesn't have a particle-temp-nrj-ratio key."
			);
		}
		boltzmann = object["particle-temp-nrj-ratio"].GetDouble();
		if (not isnormal(boltzmann) or boltzmann < 0) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Invalid particle temperature to energy ratio: " + std::to_string(boltzmann)
				+ ", should be > 0"
			);
		}

		if (object.HasMember("output-directory")) {
			output_directory = object["output-directory"].GetString();
		}

		if (not object.HasMember("load-balancer")) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "JSON data doesn't have a load-balancer key."
			);
		}
		lb_name = object["load-balancer"].GetString();
	}
};

} // namespace


#endif // ifndef PAMHD_SIMULATION_OPTIONS_HPP
