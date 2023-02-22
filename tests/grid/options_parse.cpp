/*
Parse test for grid options of PAMHD.

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


#include "prettyprint.hpp"
#include "rapidjson/document.h"

#include "grid/options.hpp"


int main()
{
	const char json1[] = "{"
		"\"grid-options\": {"
			"\"periodic\": \"{false, true, false}\","
			"\"cells\": \"{1+1, 2+2, 3+3}\","
			"\"volume\": \"{4+4 + cells[2], 5+5, 6+6 + cells[0]}\","
			"\"start\": \"{7 - volume[1] + cells[2], 8-8, 9-9}\""
		"}"
	"}";
	rapidjson::Document document;
	document.Parse(json1);
	if (document.HasParseError()) {
		return EXIT_FAILURE;
	}

	pamhd::grid::Options options;
	options.set(document);

	const auto& periodic = options.get_periodic();
	if (periodic[0] != false or periodic[1] != true or periodic[2] != false) {
		std::cerr << "Incorrect number of simulation cells: "
			<< periodic << ", should be [0, 1, 0]"
			<< std::endl;
		return EXIT_FAILURE;
	}

	const auto& nr_cells = options.get_number_of_cells();
	if (nr_cells[0] != 2 or nr_cells[1] != 4 or nr_cells[2] != 6) {
		std::cerr << "Incorrect number of simulation cells: "
			<< nr_cells << ", should be [2, 4, 6]"
			<< std::endl;
		return EXIT_FAILURE;
	}

	const auto& volume = options.get_volume();
	if (volume[0] != 14 or volume[1] != 10 or volume[2] != 14) {
		std::cerr << "Incorrect simulation volume: "
			<< volume << ", should be [14, 10, 14]"
			<< std::endl;
		return EXIT_FAILURE;
	}

	const auto& start = options.get_start();
	if (start[0] != 3 or start[1] != 0 or start[2] != 0) {
		std::cerr << "Incorrect simulation start: "
			<< start << ", should be [3, 0, 0]"
			<< std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
