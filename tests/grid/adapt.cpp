/*
AMR test for grid options of PAMHD.

Copyright 2023, 2024 Finnish Meteorological Institute
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
#include "stdexcept"
#include "string"

#include "dccrg.hpp"
#include "mpi.h"
#include "rapidjson/document.h"
#include "prettyprint.hpp"

#include "grid/amr.hpp"
#include "grid/options.hpp"

using Cell = std::array<int, 2>;

const auto Ref_Lvl_Min = [](Cell& cell_data)->auto& {
	return cell_data[0];
};
const auto Ref_Lvl_Max = [](Cell& cell_data)->auto& {
	return cell_data[1];
};



int main(int argc, char* argv[])
{
	using std::runtime_error;
	using std::to_string;

	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
		std::cerr << "Couldn't initialize MPI." << std::endl;
		abort();
	}
	MPI_Comm comm = MPI_COMM_WORLD;

	{const char json[] = "{"
		"\"grid-options\": {"
			"\"periodic\": \"{false, false, false}\","
			"\"cells\": \"{2, 2, 2}\","
			"\"volume\": \"{cells[0], cells[1], cells[2]}\","
			"\"start\": \"{0, 0, 0}\","
			"\"max-ref-lvl\": 1,"
			"\"ref-lvl-at-least\": \"(radius < 1) ? 1 : 0\","
			"\"ref-lvl-at-most\": \"(radius < 1) ? 1 : 0\""
		"}"
	"}";
	rapidjson::Document document;
	document.Parse(json);
	if (document.HasParseError()) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
	pamhd::grid::Options options;
	options.set(document);

	dccrg::Dccrg<Cell> grid; grid
		.set_initial_length(options.get_number_of_cells())
		.set_neighborhood_length(0)
		.set_periodic(options.get_periodic()[0], options.get_periodic()[1], options.get_periodic()[2])
		.set_load_balancing_method("RANDOM")
		.set_maximum_refinement_level(options.get_max_ref_lvl())
		.initialize(comm);
	pamhd::grid::set_minmax_refinement_level(
		grid.local_cells(), grid, options,
		0, Ref_Lvl_Min, Ref_Lvl_Max);
	pamhd::grid::adapt_grid(
		grid, Ref_Lvl_Min, Ref_Lvl_Max,
		pamhd::grid::New_Cells_Handler(Ref_Lvl_Min, Ref_Lvl_Max),
		pamhd::grid::Removed_Cells_Handler(Ref_Lvl_Min, Ref_Lvl_Max)
	);
	for (const auto& cell: grid.local_cells()) {
		if (cell.id == 1) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
	}}

	{const char json[] = "{"
		"\"grid-options\": {"
			"\"periodic\": \"{false, false, false}\","
			"\"cells\": \"{1, 1, 1}\","
			"\"volume\": \"{cells[0], cells[1], cells[2]}\","
			"\"start\": \"{0, 0, 0}\","
			"\"max-ref-lvl\": 1,"
			"\"ref-lvl-at-least\": 0,"
			"\"ref-lvl-at-most\": 0"
		"}"
	"}";
	rapidjson::Document document;
	document.Parse(json);
	if (document.HasParseError()) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
	pamhd::grid::Options options;
	options.set(document);

	dccrg::Dccrg<Cell> grid; grid
		.set_initial_length(options.get_number_of_cells())
		.set_neighborhood_length(0)
		.set_periodic(options.get_periodic()[0], options.get_periodic()[1], options.get_periodic()[2])
		.set_load_balancing_method("RANDOM")
		.set_maximum_refinement_level(options.get_max_ref_lvl())
		.initialize(comm);
	grid.refine_completely(1);
	grid.stop_refining();
	pamhd::grid::set_minmax_refinement_level(
		grid.local_cells(), grid, options,
		0, Ref_Lvl_Min, Ref_Lvl_Max);
	pamhd::grid::adapt_grid(
		grid, Ref_Lvl_Min, Ref_Lvl_Max,
		pamhd::grid::New_Cells_Handler(Ref_Lvl_Min, Ref_Lvl_Max),
		pamhd::grid::Removed_Cells_Handler(Ref_Lvl_Min, Ref_Lvl_Max)
	);
	for (const auto& cell: grid.local_cells()) {
		if (cell.id != 1) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
	}}

	{const char json[] = "{"
		"\"grid-options\": {"
			"\"periodic\": \"{false, false, false}\","
			"\"cells\": \"{1, 1, 1}\","
			"\"volume\": \"{cells[0], cells[1], cells[2]}\","
			"\"start\": \"{0, 0, 0}\","
			"\"max-ref-lvl\": 2,"
			"\"ref-lvl-at-least\": \"(radius < 2) ? 2 : 0\","
			"\"ref-lvl-at-most\": \"(radius < 2) ? 2 : 0\""
		"}"
	"}";
	rapidjson::Document document;
	document.Parse(json);
	if (document.HasParseError()) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
	pamhd::grid::Options options;
	options.set(document);

	dccrg::Dccrg<Cell> grid; grid
		.set_initial_length(options.get_number_of_cells())
		.set_neighborhood_length(0)
		.set_periodic(options.get_periodic()[0], options.get_periodic()[1], options.get_periodic()[2])
		.set_load_balancing_method("RANDOM")
		.set_maximum_refinement_level(options.get_max_ref_lvl())
		.initialize(comm);
	pamhd::grid::set_minmax_refinement_level(
		grid.local_cells(), grid, options,
		0, Ref_Lvl_Min, Ref_Lvl_Max);
	pamhd::grid::adapt_grid(
		grid, Ref_Lvl_Min, Ref_Lvl_Max,
		pamhd::grid::New_Cells_Handler(Ref_Lvl_Min, Ref_Lvl_Max),
		pamhd::grid::Removed_Cells_Handler(Ref_Lvl_Min, Ref_Lvl_Max)
	);
	for (const auto& cell: grid.local_cells()) {
		const auto ref_lvl = grid.get_refinement_level(cell.id);
		if (ref_lvl < 2) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
	}}

	{const char json[] = "{"
		"\"grid-options\": {"
			"\"periodic\": \"{false, false, false}\","
			"\"cells\": \"{1, 1, 1}\","
			"\"volume\": \"{cells[0], cells[1], cells[2]}\","
			"\"start\": \"{0, 0, 0}\","
			"\"max-ref-lvl\": 3,"
			"\"ref-lvl-at-least\": \"(radius < 2) ? 3 : 0\","
			"\"ref-lvl-at-most\": \"(radius < 2) ? 3 : 0\""
		"}"
	"}";
	rapidjson::Document document;
	document.Parse(json);
	if (document.HasParseError()) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
	pamhd::grid::Options options;
	options.set(document);

	dccrg::Dccrg<Cell> grid; grid
		.set_initial_length(options.get_number_of_cells())
		.set_neighborhood_length(0)
		.set_periodic(options.get_periodic()[0], options.get_periodic()[1], options.get_periodic()[2])
		.set_load_balancing_method("RANDOM")
		.set_maximum_refinement_level(options.get_max_ref_lvl())
		.initialize(comm);
	pamhd::grid::set_minmax_refinement_level(
		grid.local_cells(), grid, options,
		0, Ref_Lvl_Min, Ref_Lvl_Max);
	pamhd::grid::adapt_grid(
		grid, Ref_Lvl_Min, Ref_Lvl_Max,
		pamhd::grid::New_Cells_Handler(Ref_Lvl_Min, Ref_Lvl_Max),
		pamhd::grid::Removed_Cells_Handler(Ref_Lvl_Min, Ref_Lvl_Max)
	);
	for (const auto& cell: grid.local_cells()) {
		const auto ref_lvl = grid.get_refinement_level(cell.id);
		if (ref_lvl < 3) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
	}}

	MPI_Finalize();
	return EXIT_SUCCESS;
}
