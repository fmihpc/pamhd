/*
Tests primary cell edge calculation of PAMHD.

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


#include "iostream"
#include "stdexcept"
#include "string"

#include "dccrg.hpp"
#include "dccrg_no_geometry.hpp"
#include "mpi.h" // must be included before gensimcell.hpp

#include "common.hpp"


using namespace std;

const auto PEdge = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::grid::Is_Primary_Edge()];
};

int main(int argc, char* argv[])
{
	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
		cerr << "Couldn't initialize MPI." << endl;
		abort();
	}
	MPI_Comm comm = MPI_COMM_WORLD;
	float zoltan_version = -1.0;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
		cerr << "Zoltan_Initialize failed." << endl;
		abort();
	}

	// 2x2 cells of which one is refined
	{Grid grid; grid
		.set_periodic(true, true, true)
		.set_initial_length({2, 2, 1})
		.set_neighborhood_length(1)
		.set_maximum_refinement_level(1)
		.initialize(comm);
	for (const auto& cell: grid.local_cells()) if (cell.id == 1) grid.refine_completely(cell.id);
	grid.stop_refining();
	pamhd::grid::update_primary_edges(grid.local_cells(), grid, PEdge);
	for (const auto& cell: grid.local_cells()) {
		if (cell.id < 2 or cell.id > 26) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
		const auto pe = PEdge(*cell.data);
		// test  only diagonal neighbor logic
		if (cell.id != 4) continue;
		for (auto d: {0, 1, 2})
		for (auto a: {-1, +1})
		for (auto b: {-1, +1})
			if (d < 2 and a == 1 and b == 1) {
				if (pe(d,a,b) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			} else {
				if (pe(d,a,b) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			}
	}}

	{Grid grid; grid
		.set_periodic(true, true, true)
		.set_initial_length({2, 1, 2})
		.set_neighborhood_length(1)
		.set_maximum_refinement_level(1)
		.initialize(comm);
	for (const auto& cell: grid.local_cells()) if (cell.id == 1) grid.refine_completely(cell.id);
	grid.stop_refining();
	pamhd::grid::update_primary_edges(grid.local_cells(), grid, PEdge);
	for (const auto& cell: grid.local_cells()) {
		if (cell.id < 2 or cell.id > 26) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
		const auto pe = PEdge(*cell.data);
		// test  only diagonal neighbor logic
		if (cell.id != 4) continue;
		for (auto d: {0, 1, 2})
		for (auto a: {-1, +1})
		for (auto b: {-1, +1})
			if (d != 1 and a == 1 and b == 1) {
				if (pe(d,a,b) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			} else {
				if (pe(d,a,b) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			}
	}}

	{Grid grid; grid
		.set_periodic(true, true, true)
		.set_initial_length({1, 2, 2})
		.set_neighborhood_length(1)
		.set_maximum_refinement_level(1)
		.initialize(comm);
	for (const auto& cell: grid.local_cells()) if (cell.id == 1) grid.refine_completely(cell.id);
	grid.stop_refining();
	pamhd::grid::update_primary_edges(grid.local_cells(), grid, PEdge);
	for (const auto& cell: grid.local_cells()) {
		if (cell.id < 2 or cell.id > 26) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
		const auto pe = PEdge(*cell.data);
		// test  only diagonal neighbor logic
		if (cell.id != 4) continue;
		for (auto d: {0, 1, 2})
		for (auto a: {-1, +1})
		for (auto b: {-1, +1})
			if (d > 0 and a == 1 and b == 1) {
				if (pe(d,a,b) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			} else {
				if (pe(d,a,b) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			}
	}}

	// 2x2 grid of which two diagonal neighbors are refined
	{Grid grid; grid
		.set_periodic(true, true, true)
		.set_initial_length({2, 2, 1})
		.set_neighborhood_length(1)
		.set_maximum_refinement_level(1)
		.initialize(comm);
	for (const auto& cell: grid.local_cells()) if (cell.id == 1 or cell.id == 4) grid.refine_completely(cell.id);
	grid.stop_refining();
	pamhd::grid::update_primary_edges(grid.local_cells(), grid, PEdge);
	for (const auto& cell: grid.local_cells()) {
		if (cell.id < 2 or cell.id > 36) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
		const auto pe = PEdge(*cell.data);
		// test  only diagonal neighbor logic
		switch (cell.id) {
		case 5:
		case 15:
		case 21:
		case 31:
			if (pe(2,-1,-1) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			break;
		case 6:
		case 16:
		case 22:
		case 32:
			if (pe(2,+1,-1) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			break;
		case 9:
		case 19:
		case 25:
		case 35:
			if (pe(2,-1,+1) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			break;
		case 10:
		case 20:
		case 26:
		case 36:
			if (pe(2,+1,+1) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			break;
		default:
			break;
		}
	}}

	{Grid grid; grid
		.set_periodic(true, true, true)
		.set_initial_length({2, 1, 2})
		.set_neighborhood_length(1)
		.set_maximum_refinement_level(1)
		.initialize(comm);
	for (const auto& cell: grid.local_cells()) if (cell.id == 1 or cell.id == 4) grid.refine_completely(cell.id);
	grid.stop_refining();
	pamhd::grid::update_primary_edges(grid.local_cells(), grid, PEdge);
	for (const auto& cell: grid.local_cells()) {
		if (cell.id < 2 or cell.id > 36) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
		const auto pe = PEdge(*cell.data);
		// test  only diagonal neighbor logic
		switch (cell.id) {
		case 5:
		case 9:
		case 23:
		case 27:
			if (pe(1,-1,-1) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			break;
		case 6:
		case 10:
		case 24:
		case 28:
			if (pe(1,+1,-1) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			break;
		case 13:
		case 17:
		case 31:
		case 35:
			if (pe(1,-1,+1) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			break;
		case 14:
		case 18:
		case 32:
		case 36:
			if (pe(1,+1,+1) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			break;
		default:
			break;
		}
	}}

	{Grid grid; grid
		.set_periodic(true, true, true)
		.set_initial_length({1, 2, 2})
		.set_neighborhood_length(1)
		.set_maximum_refinement_level(1)
		.initialize(comm);
	for (const auto& cell: grid.local_cells()) if (cell.id == 1 or cell.id == 4) grid.refine_completely(cell.id);
	grid.stop_refining();
	pamhd::grid::update_primary_edges(grid.local_cells(), grid, PEdge);
	for (const auto& cell: grid.local_cells()) {
		if (cell.id < 2 or cell.id > 36) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
		const auto pe = PEdge(*cell.data);
		// test  only diagonal neighbor logic
		switch (cell.id) {
		case 5:
		case 6:
		case 25:
		case 26:
			if (pe(0,-1,-1) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			break;
		case 7:
		case 8:
		case 27:
		case 28:
			if (pe(0,+1,-1) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			break;
		case 13:
		case 14:
		case 33:
		case 34:
			if (pe(0,-1,+1) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			break;
		case 15:
		case 16:
		case 35:
		case 36:
			if (pe(0,+1,+1) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			break;
		default:
			break;
		}
	}}

	MPI_Finalize();
	return EXIT_SUCCESS;
}
