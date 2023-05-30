/*
AMR test for grid options of PAMHD.

Copyright 2023 Finnish Meteorological Institute
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


#include "stdexcept"
#include "string"

#include "rapidjson/document.h"

#include "grid/amr.hpp"
#include "grid/options.hpp"


int main()
{
	using std::runtime_error;
	using std::to_string;

	{const char json[] = "{"
		"\"grid-options\": {"
			"\"periodic\": \"{false, false, false}\","
			"\"cells\": \"{1, 1, 1}\","
			"\"volume\": \"{1, 1, 1}\","
			"\"start\": \"{0, 0, 0}\","
			"\"max-ref-lvl\": 0,"
			"\"ref-lvl-at-least\": 0,"
			"\"ref-lvl-at-most\": 0"
		"}"
	"}";
	rapidjson::Document document;
	document.Parse(json);
	if (document.HasParseError()) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
	pamhd::grid::Options options;
	options.set(document);

	std::array<double, 3> c{0, 0, 0};
	const auto tgt_ref_lvl = pamhd::grid::get_target_refinement_level(options, 0.0, c);
	if (tgt_ref_lvl[0] != 0) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
	if (tgt_ref_lvl[1] != 0) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
	}

	{const char json[] = "{"
		"\"grid-options\": {"
			"\"periodic\": \"{false, false, false}\","
			"\"cells\": \"{1, 1, 1}\","
			"\"volume\": \"{1, 1, 1}\","
			"\"start\": \"{0, 0, 0}\","
			"\"max-ref-lvl\": 2,"
			"\"ref-lvl-at-least\": 1,"
			"\"ref-lvl-at-most\": 2"
		"}"
	"}";
	rapidjson::Document document;
	document.Parse(json);
	if (document.HasParseError()) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
	pamhd::grid::Options options;
	options.set(document);

	std::array<double, 3> c{0, 0, 0};
	const auto tgt_ref_lvl = pamhd::grid::get_target_refinement_level(options, 0.0, c);
	if (tgt_ref_lvl[0] != 1) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
	if (tgt_ref_lvl[1] != 2) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
	}

	{const char json[] = "{"
		"\"grid-options\": {"
			"\"periodic\": \"{false, false, false}\","
			"\"cells\": \"{1, 1, 1}\","
			"\"volume\": \"{1, 1, 1}\","
			"\"start\": \"{0, 0, 0}\","
			"\"max-ref-lvl\": 2,"
			"\"ref-lvl-at-least\": \"1\","
			"\"ref-lvl-at-most\": 2"
		"}"
	"}";
	rapidjson::Document document;
	document.Parse(json);
	if (document.HasParseError()) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
	pamhd::grid::Options options;
	options.set(document);

	std::array<double, 3> c{0, 0, 0};
	const auto tgt_ref_lvl = pamhd::grid::get_target_refinement_level(options, 0.0, c);
	if (tgt_ref_lvl[0] != 1) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
	if (tgt_ref_lvl[1] != 2) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
	}

	{const char json[] = "{"
		"\"grid-options\": {"
			"\"periodic\": \"{false, false, false}\","
			"\"cells\": \"{1, 1, 1}\","
			"\"volume\": \"{1, 1, 1}\","
			"\"start\": \"{0, 0, 0}\","
			"\"max-ref-lvl\": 2,"
			"\"ref-lvl-at-least\": 1,"
			"\"ref-lvl-at-most\": \"2\""
		"}"
	"}";
	rapidjson::Document document;
	document.Parse(json);
	if (document.HasParseError()) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
	pamhd::grid::Options options;
	options.set(document);

	std::array<double, 3> c{0, 0, 0};
	const auto tgt_ref_lvl = pamhd::grid::get_target_refinement_level(options, 0.0, c);
	if (tgt_ref_lvl[0] != 1) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
	if (tgt_ref_lvl[1] != 2) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
	}

	{const char json[] = "{"
		"\"grid-options\": {"
			"\"periodic\": \"{false, false, false}\","
			"\"cells\": \"{1, 1, 1}\","
			"\"volume\": \"{1, 1, 1}\","
			"\"start\": \"{0, 0, 0}\","
			"\"max-ref-lvl\": 2,"
			"\"ref-lvl-at-least\": \"1\","
			"\"ref-lvl-at-most\": \"2\""
		"}"
	"}";
	rapidjson::Document document;
	document.Parse(json);
	if (document.HasParseError()) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
	pamhd::grid::Options options;
	options.set(document);

	std::array<double, 3> c{0, 0, 0};
	const auto tgt_ref_lvl = pamhd::grid::get_target_refinement_level(options, 0.0, c);
	if (tgt_ref_lvl[0] != 1) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
	if (tgt_ref_lvl[1] != 2) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
	}

	{const char json[] = "{"
		"\"grid-options\": {"
			"\"periodic\": \"{false, false, false}\","
			"\"cells\": \"{1, 1, 1}\","
			"\"volume\": \"{1, 1, 1}\","
			"\"start\": \"{0, 0, 0}\","
			"\"max-ref-lvl\": 3,"
			"\"ref-lvl-at-least\": \"(radius < 2) ? 1 : 3\","
			"\"ref-lvl-at-most\": \"2\""
		"}"
	"}";
	rapidjson::Document document;
	document.Parse(json);
	if (document.HasParseError()) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
	pamhd::grid::Options options;
	options.set(document);

	std::array<double, 3> c{0, 0, 0};
	const auto tgt_ref_lvl = pamhd::grid::get_target_refinement_level(options, 0.0, c);
	if (tgt_ref_lvl[0] != 1) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
	if (tgt_ref_lvl[1] != 2) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
	}

	{const char json[] = "{"
		"\"grid-options\": {"
			"\"periodic\": \"{false, false, false}\","
			"\"cells\": \"{1, 1, 1}\","
			"\"volume\": \"{1, 1, 1}\","
			"\"start\": \"{0, 0, 0}\","
			"\"max-ref-lvl\": 3,"
			"\"ref-lvl-at-least\": \"(radius < 2) ? 1 : 3\","
			"\"ref-lvl-at-most\": \"(x > -1) ? 2 : 0\""
		"}"
	"}";
	rapidjson::Document document;
	document.Parse(json);
	if (document.HasParseError()) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
	pamhd::grid::Options options;
	options.set(document);

	std::array<double, 3> c{0, 0, 0};
	const auto tgt_ref_lvl = pamhd::grid::get_target_refinement_level(options, 0.0, c);
	if (tgt_ref_lvl[0] != 1) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
	if (tgt_ref_lvl[1] != 2) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
	}

	return EXIT_SUCCESS;
}
