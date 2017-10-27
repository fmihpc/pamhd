/*
Tests initial conditions class of PAMHD with scalar simulation variable.

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

#include "array"
#include "cstdlib"
#include "iostream"
#include "set"
#include "string"
#include "vector"

#include "boundaries/initial_conditions.hpp"


using namespace pamhd::boundaries;


struct Mass_Density {
	using data_type = double;
	static const std::string get_name(){ return {"mass density"}; }
};

int main()
{
	typename Mass_Density::data_type mass;

	const char json1[] =
		"{\"default\": -3,"
		"\"initial-conditions\": ["
			"{"
				"\"geometry-id\": 1,"
				"\"value\": 1"
			"}"
		"]}";
	rapidjson::Document doc1;
	doc1.Parse(json1);
	if (doc1.HasParseError()) {
		std::cerr << "Couldn't parse json data: " << json1 << std::endl;
		return EXIT_FAILURE;
	}

	Initial_Conditions<unsigned int, Mass_Density> init1;
	init1.set(doc1);
	mass = init1.get_default_data(0, 0, 0, 0, 0, 0, 0);
	if (mass != -3) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong default value for variable: "
			<< mass << ", should be -3"
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (init1.get_number_of_regions() != 1) {
		std::cerr << __FILE__ "(" << __LINE__
			<< "): Wrong number of non-default initial conditions: "
			<< init1.get_number_of_regions() << ", should be 1"
			<< std::endl;
		return EXIT_FAILURE;
	}
	mass = init1.get_data(0, 0, 0, 0, 0, 0, 0, 0);
	if (mass != 1) {
		std::cerr << __FILE__ "(" << __LINE__
			<< "): Wrong value for variable in non-default initial condition: "
			<< mass << ", should be 1"
			<< std::endl;
		return EXIT_FAILURE;
	}


	const char json2[] =
		"{\"default\": \"t+2\","
		"\"initial-conditions\": ["
			"{"
				"\"geometry-id\": 1,"
				"\"value\": 1"
			"}"
		"]}";
	rapidjson::Document doc2;
	doc2.Parse(json2);
	if (doc2.HasParseError()) {
		std::cerr << "Couldn't parse json data: " << json2 << std::endl;
		return EXIT_FAILURE;
	}

	Initial_Conditions<unsigned int, Mass_Density> init2;
	init2.set(doc2);
	mass = init2.get_default_data(3, 0, 0, 0, 0, 0, 0);
	if (mass != 5) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong default value for variable: "
			<< mass << ", should be 5"
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (init2.get_number_of_regions() != 1) {
		std::cerr << __FILE__ "(" << __LINE__
			<< "): Wrong number of non-default initial conditions: "
			<< init2.get_number_of_regions() << ", should be 1"
			<< std::endl;
		return EXIT_FAILURE;
	}
	mass = init2.get_data(0, 0, 0, 0, 0, 0, 0, 0);
	if (mass != 1) {
		std::cerr << __FILE__ "(" << __LINE__
			<< "): Wrong value for variable in non-default initial condition: "
			<< mass << ", should be 1"
			<< std::endl;
		return EXIT_FAILURE;
	}


	const char json3[] =
		"{\"default\": {"
			"\"x\": [1],"
			"\"y\": [-2],"
			"\"z\": [0.1],"
			"\"data\": [-1]"
		"},"
		"\"initial-conditions\": ["
			"{"
				"\"geometry-id\": 1,"
				"\"value\": 1"
			"}"
		"]}";
	rapidjson::Document doc3;
	doc3.Parse(json3);
	if (doc3.HasParseError()) {
		std::cerr << "Couldn't parse json data: " << json3 << std::endl;
		return EXIT_FAILURE;
	}

	Initial_Conditions<unsigned int, Mass_Density> init3;
	init3.set(doc3);
	mass = init3.get_default_data(3, 0, 0, 0, 0, 0, 0);
	if (mass != -1) {
		std::cerr << __FILE__ "(" << __LINE__ << "): Wrong default value for variable: "
			<< mass << ", should be -1"
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (init3.get_number_of_regions() != 1) {
		std::cerr << __FILE__ "(" << __LINE__
			<< "): Wrong number of non-default initial conditions: "
			<< init3.get_number_of_regions() << ", should be 1"
			<< std::endl;
		return EXIT_FAILURE;
	}
	mass = init3.get_data(0, 0, 0, 0, 0, 0, 0, 0);
	if (mass != 1) {
		std::cerr << __FILE__ "(" << __LINE__
			<< "): Wrong value for variable in non-default initial condition: "
			<< mass << ", should be 1"
			<< std::endl;
		return EXIT_FAILURE;
	}


	return EXIT_SUCCESS;
}
