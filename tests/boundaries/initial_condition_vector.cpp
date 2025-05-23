/*
Tests initial condition class of PAMHD with vector simulation variable.

Copyright 2016, 2017 Ilja Honkonen
Copyright 2025 Finnish Meteorological Institute
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

#include "array"
#include "cstdlib"
#include "iostream"
#include "set"
#include "string"
#include "vector"

#include "boundaries/initial_condition.hpp"


using namespace pamhd::boundaries;


struct Momentum_Density {
	using data_type = std::array<double, 2>;
	static const std::string get_name(){ return {"momentum density"}; }
};

int main()
{
	const char json[] = "["
		"{"
			"\"geometry-id\": 1,"
			"\"value\": [1, 2]"
		"},"
		"{"
			"\"geometry-id\": 3,"
			"\"value\": \"{t, t*t}\""
		"},"
		"{"
			"\"geometry-id\": 8,"
			"\"value\": {"
				"\"x\": [1],"
				"\"y\": [-2],"
				"\"z\": [0.1],"
				"\"data\": [[-1, 0]]"
			"}"
		"}"
	"]";

	rapidjson::Document document;
	document.Parse(json);
	if (document.HasParseError()) {
		std::cerr << "Couldn't parse json data: " << json << std::endl;
		return EXIT_FAILURE;
	}

	typename Momentum_Density::data_type momentum;

	Initial_Condition<unsigned int, Momentum_Density> init0, init1, init2;

	/*{
		"geometry-id": 1,
		"value": [1, 2]
	}*/
	init0.set(document[0]);
	momentum = init0.get_data(0, 0, 0, 0, 0, 0, 0);
	if (momentum[0] != 1) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< momentum[0] << ", should be 1"
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (momentum[1] != 2) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< momentum[1] << ", should be 2"
			<< std::endl;
		return EXIT_FAILURE;
	}


	/*{
		"geometry-id": 3,
		"value": "{t, t*t}"
	}*/
	init1.set(document[1]);
	momentum = init1.get_data(0.5, 0, 0, 0, 0, 0, 0);
	if (momentum[0] != 0.5) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< momentum[0] << ", should be 0.5"
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (momentum[1] != 0.25) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< momentum[1] << ", should be 0.25"
			<< std::endl;
		return EXIT_FAILURE;
	}


	/*{
		"geometry-id": 8,
		"value": {
			"x": [1],
			"y": [-2],
			"z": [0.1],
			"data": [[-1, 0]]
		}
	}*/
	init2.set(document[2]);
	momentum = init2.get_data(0, 0, 0, 0, 0, 0, 0);
	if (momentum[0] != -1) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< momentum[0] << ", should be -1"
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (momentum[1] != 0) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong value for variable: "
			<< momentum[1] << ", should be 0"
			<< std::endl;
		return EXIT_FAILURE;
	}


	return EXIT_SUCCESS;
}
