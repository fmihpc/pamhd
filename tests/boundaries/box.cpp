/*
Box boundary geometry test of PAMHD.

Copyright 2014, 2015, 2016, 2017 Ilja Honkonen
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

#define MAKE_BOX(a, b, c, d, e, f) {a, b, c}, {d, e, f}

#include "rapidjson/document.h"

#include "boundaries/box.hpp"


using namespace pamhd::boundaries;

int main()
{
	const char json[] = "{"
		"\"asdf\": ["
			"{\"box\": {"
				"\"start\": [-1, -2, -3],"
				"\"end\": [1, 2, 3]"
			"}},"
			"{\"sphere\": {"
				"\"center\": [-1, -2, -3],"
				"\"radius\": 0.01"
			"}},"
			"{\"box\": {"
				"\"start\": [1, 2, 3],"
				"\"end\": [-1, -2, -3]"
			"}}"
		"]"
	"}";

	rapidjson::Document document;
	document.Parse(json);
	if (document.HasParseError()) {
		std::cerr << "Couldn't parse json data: " << json << std::endl;
		return EXIT_FAILURE;
	}

	try {
		const auto box_fail = Box<
			std::array<double, 3>, int
		>(document["asdf"][2]["box"]);

		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Could create box from rapidjson object."
			<< std::endl;
		return EXIT_FAILURE;
	} catch (...) {
	}


	Box<std::array<double, 3>, int> box;

	try {
		box = Box<
			std::array<double, 3>, int
		>(document["asdf"][0]["box"]);
	} catch (...) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Couldn't create box from rapidjson object."
			<< std::endl;
		return EXIT_FAILURE;
	}

	if (box.overlaps(MAKE_BOX(-2, -3, -4, -1.1, -2.1, -3.1), 1)) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Cell overlaps box."
			<< std::endl;
		return EXIT_FAILURE;
	}

	if (box.overlaps(MAKE_BOX(1.1, 2.1, 3.1, 2, 3, 4), 2)) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Cell overlaps box."
			<< std::endl;
		return EXIT_FAILURE;
	}

	if (not box.overlaps(MAKE_BOX(0, 1, 2, 5, 5, 5), 3)) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Cell does not overlap box."
			<< std::endl;
		return EXIT_FAILURE;
	}

	if (box.cells.size() != 1) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Wrong number of cells in box's list of overlapping cells: "
			<< box.cells.size() << ", should be 1."
			<< std::endl;
		return EXIT_FAILURE;
	}
	if (box.cells[0] != 3) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Wrong cell in box's list of overlapping cells: "
			<< box.cells[0] << ", should be 3."
			<< std::endl;
		return EXIT_FAILURE;
	}


	return EXIT_SUCCESS;
}
