/*
Tests primary cell edge calculation of PAMHD.

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

	/*
	Two-cell grids with different direction, periodicities with AMR
	*/
	{Grid grid; grid
		.set_periodic(true, true, true)
		.set_initial_length({2, 1, 1})
		.set_neighborhood_length(1)
		.set_maximum_refinement_level(1)
		.initialize(comm);
	for (const auto& cell: grid.local_cells()) if (cell.id == 1) grid.refine_completely(cell.id);
	grid.stop_refining();
	pamhd::grid::update_primary_edges(grid.local_cells(), grid, PEdge);
	for (const auto& cell: grid.local_cells()) {
		if (cell.id < 2 or cell.id > 18) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
		const auto pe = PEdge(*cell.data);
		if (pe(0,0,0) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
		if (pe(0,0,1) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
		if (pe(0,1,0) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
		if (pe(0,1,1) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
		if (cell.id == 2) {
			for (size_t d: {1, 2})
			for (size_t a: {0, 1})
			for (size_t b: {0, 1})
				if (pe(d,a,b) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			continue;
		}
		if (cell.id % 2 == 1) {
			for (size_t d: {1, 2})
			for (size_t a: {0, 1})
			for (size_t b: {0, 1})
				if (b == 1) {
					if (pe(d,a,b) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
				} else {
					if (pe(d,a,b) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
				}
		} else {
			for (size_t d: {1, 2})
			for (size_t a: {0, 1})
			for (size_t b: {0, 1})
				if (a == 1 and b == 1) {
					if (pe(d,a,b) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
				} else {
					if (pe(d,a,b) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
				}
		}
	}}

	{Grid grid; grid
		.set_periodic(false, false, false)
		.set_initial_length({2, 1, 1})
		.set_neighborhood_length(1)
		.set_maximum_refinement_level(1)
		.initialize(comm);
	for (const auto& cell: grid.local_cells()) if (cell.id == 1) grid.refine_completely(cell.id);
	grid.stop_refining();
	pamhd::grid::update_primary_edges(grid.local_cells(), grid, PEdge);
	for (const auto& cell: grid.local_cells()) {
		if (cell.id < 2 or cell.id > 18) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
		const auto pe = PEdge(*cell.data);
		if (cell.id == 2) {
			for (size_t d: {0, 1, 2})
			for (size_t a: {0, 1})
			for (size_t b: {0, 1})
				if (d == 0) {
					if (pe(d,a,b) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
				} else if (a == 0) {
					if (pe(d,a,b) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
				} else {
					if (pe(d,a,b) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
				}
		} else {
			if (cell.id < 5) {
				if (pe(0,0,0) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			} else {
				if (pe(0,0,0) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			}
			if (cell.id % 8 == 3 or cell.id % 8 == 4) {
				if (pe(0,0,1) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			} else {
				if (pe(0,0,1) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			}
			if (cell.id < 9) {
				if (pe(0,1,0) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			} else {
				if (pe(0,1,0) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			}
			if (pe(0,1,1) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			if (cell.id == 3 or cell.id == 7) {
				if (pe(1,0,0) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			} else {
				if (pe(1,0,0) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			}
			if (cell.id % 2 == 1) {
				if (pe(1,0,1) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			} else {
				if (pe(1,0,1) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			}
			if (cell.id < 9) {
				if (pe(1,1,0) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			} else {
				if (pe(1,1,0) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			}
			if (pe(1,1,1) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			if (cell.id == 3 or cell.id == 11) {
				if (pe(2,0,0) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			} else {
				if (pe(2,0,0) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			}
			if (cell.id % 2 == 1) {
				if (pe(2,0,1) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			} else {
				if (pe(2,0,1) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			}
			if (cell.id % 8 == 3 or cell.id % 8 == 4) {
				if (pe(2,1,0) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			} else {
				if (pe(2,1,0) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			}
			if (pe(2,1,1) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
		}
	}}

	{Grid grid; grid
		.set_periodic(true, true, false)
		.set_initial_length({2, 1, 1})
		.set_neighborhood_length(1)
		.set_maximum_refinement_level(1)
		.initialize(comm);
	for (const auto& cell: grid.local_cells()) if (cell.id == 1) grid.refine_completely(cell.id);
	grid.stop_refining();
	pamhd::grid::update_primary_edges(grid.local_cells(), grid, PEdge);
	for (const auto& cell: grid.local_cells()) {
		if (cell.id < 2 or cell.id > 18) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
		const auto pe = PEdge(*cell.data);
		if (cell.id == 2) {
			for (size_t d: {0, 1, 2})
			for (size_t a: {0, 1})
			for (size_t b: {0, 1})
				if (d == 0 and a != 0) {
					if (pe(d,a,b) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
				} else {
					if (pe(d,a,b) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
				}
		} else {
			if (pe(0,0,0) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			if (pe(0,0,1) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			if (cell.id < 9) {
				if (pe(0,1,0) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			} else {
				if (pe(0,1,0) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			}
			if (pe(0,1,1) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			if (cell.id == 3 or cell.id == 7) {
				if (pe(1,0,0) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			} else {
				if (pe(1,0,0) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			}
			if (cell.id % 2 == 1) {
				if (pe(1,0,1) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			} else {
				if (pe(1,0,1) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			}
			if (cell.id < 9) {
				if (pe(1,1,0) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			} else {
				if (pe(1,1,0) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			}
			if (pe(1,1,1) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			if (pe(2,0,0) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			if (cell.id % 2 == 1) {
				if (pe(2,0,1) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			} else {
				if (pe(2,0,1) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			}
			if (pe(2,1,0) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			if (pe(2,1,1) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
		}
	}}

	{Grid grid; grid
		.set_periodic(true, false, true)
		.set_initial_length({2, 1, 1})
		.set_neighborhood_length(1)
		.set_maximum_refinement_level(1)
		.initialize(comm);
	for (const auto& cell: grid.local_cells()) if (cell.id == 1) grid.refine_completely(cell.id);
	grid.stop_refining();
	pamhd::grid::update_primary_edges(grid.local_cells(), grid, PEdge);
	for (const auto& cell: grid.local_cells()) {
		if (cell.id < 2 or cell.id > 18) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
		const auto pe = PEdge(*cell.data);
		if (cell.id == 2) {
			for (size_t d: {0, 1, 2})
			for (size_t a: {0, 1})
			for (size_t b: {0, 1})
				if (d == 0 and b != 0) {
					if (pe(d,a,b) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
				} else {
					if (pe(d,a,b) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
				}
		} else {
			if (pe(0,0,0) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			if (cell.id % 8 == 3 or cell.id % 8 == 4) {
				if (pe(0,0,1) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			} else {
				if (pe(0,0,1) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			}
			if (pe(0,1,0) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			if (pe(0,1,1) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			if (pe(1,0,0) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			if (cell.id % 2 == 1) {
				if (pe(1,0,1) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			} else {
				if (pe(1,0,1) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			}
			if (pe(1,1,0) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			if (pe(1,1,1) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			if (cell.id == 3 or cell.id == 11) {
				if (pe(2,0,0) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			} else {
				if (pe(2,0,0) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			}
			if (cell.id % 2 == 1) {
				if (pe(2,0,1) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			} else {
				if (pe(2,0,1) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			}
			if (cell.id % 8 == 3 or cell.id % 8 == 4) {
				if (pe(2,1,0) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			} else {
				if (pe(2,1,0) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			}
			if (pe(2,1,1) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
		}
	}}

	{Grid grid; grid
		.set_periodic(false, true, true)
		.set_initial_length({2, 1, 1})
		.set_neighborhood_length(1)
		.set_maximum_refinement_level(1)
		.initialize(comm);
	for (const auto& cell: grid.local_cells()) if (cell.id == 1) grid.refine_completely(cell.id);
	grid.stop_refining();
	pamhd::grid::update_primary_edges(grid.local_cells(), grid, PEdge);
	for (const auto& cell: grid.local_cells()) {
		if (cell.id < 2 or cell.id > 18) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
		const auto pe = PEdge(*cell.data);
		for (size_t d: {0, 1, 2})
		for (size_t a: {0, 1})
		for (size_t b: {0, 1})
			if (
				(a == 1 and b == 1)
				or (d > 0 and a == 0 and b == 1 and cell.id % 2 == 1)
			) {
				if (pe(d,a,b) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			} else {
				if (pe(d,a,b) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			}
	}}

	{Grid grid; grid
		.set_periodic(true, false, true)
		.set_initial_length({1, 2, 1})
		.set_neighborhood_length(1)
		.set_maximum_refinement_level(1)
		.initialize(comm);
	for (const auto& cell: grid.local_cells()) if (cell.id == 1) grid.refine_completely(cell.id);
	grid.stop_refining();
	pamhd::grid::update_primary_edges(grid.local_cells(), grid, PEdge);
	for (const auto& cell: grid.local_cells()) {
		if (cell.id < 2 or cell.id > 18) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
		const auto pe = PEdge(*cell.data);
		for (size_t d: {0, 1, 2})
		for (size_t a: {0, 1})
		for (size_t b: {0, 1})
			if (
				(a == 1 and b == 1)
				or
					((cell.id % 4 == 0 or cell.id % 4 == 3)
					and
						((d == 0 and a == 0 and b == 1)
						or (d == 2 and a == 1 and b == 0)))
			) {
				if (pe(d,a,b) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			} else {
				if (pe(d,a,b) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			}
	}}

	{Grid grid; grid
		.set_periodic(true, true, false)
		.set_initial_length({1, 1, 2})
		.set_neighborhood_length(1)
		.set_maximum_refinement_level(1)
		.initialize(comm);
	for (const auto& cell: grid.local_cells()) if (cell.id == 1) grid.refine_completely(cell.id);
	grid.stop_refining();
	pamhd::grid::update_primary_edges(grid.local_cells(), grid, PEdge);
	for (const auto& cell: grid.local_cells()) {
		if (cell.id < 2 or cell.id > 18) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
		const auto pe = PEdge(*cell.data);
		for (size_t d: {0, 1, 2})
		for (size_t a: {0, 1})
		for (size_t b: {0, 1})
			if (
				(a == 1 and b == 1)
				or (
					d < 2 and a == 1 and b == 0
					and cell.id > 2 and cell.id < 7)
			) {
				if (pe(d,a,b) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			} else {
				if (pe(d,a,b) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			}
	}}

	// shouldn't depend on neighborhood length as long as > 0
	{Grid grid; grid
		.set_periodic(true, true, false)
		.set_initial_length({1, 1, 2})
		.set_neighborhood_length(2)
		.set_maximum_refinement_level(1)
		.initialize(comm);
	for (const auto& cell: grid.local_cells()) if (cell.id == 1) grid.refine_completely(cell.id);
	grid.stop_refining();
	pamhd::grid::update_primary_edges(grid.local_cells(), grid, PEdge);
	for (const auto& cell: grid.local_cells()) {
		if (cell.id < 2 or cell.id > 18) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
		const auto pe = PEdge(*cell.data);
		for (size_t d: {0, 1, 2})
		for (size_t a: {0, 1})
		for (size_t b: {0, 1})
			if (
				(a == 1 and b == 1)
				or (
					d < 2 and a == 1 and b == 0
					and cell.id > 2 and cell.id < 7)
			) {
				if (pe(d,a,b) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			} else {
				if (pe(d,a,b) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			}
	}}

	// shouldn't depend on max ref level as long as > 0
	{Grid grid; grid
		.set_periodic(true, true, false)
		.set_initial_length({1, 1, 2})
		.set_neighborhood_length(1)
		.set_maximum_refinement_level(2)
		.initialize(comm);
	for (const auto& cell: grid.local_cells()) if (cell.id == 1) grid.refine_completely(cell.id);
	grid.stop_refining();
	pamhd::grid::update_primary_edges(grid.local_cells(), grid, PEdge);
	for (const auto& cell: grid.local_cells()) {
		if (cell.id < 2 or cell.id > 18) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
		const auto pe = PEdge(*cell.data);
		for (size_t d: {0, 1, 2})
		for (size_t a: {0, 1})
		for (size_t b: {0, 1})
			if (
				(a == 1 and b == 1)
				or (
					d < 2 and a == 1 and b == 0
					and cell.id > 2 and cell.id < 7)
			) {
				if (pe(d,a,b) == false) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			} else {
				if (pe(d,a,b) == true) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			}
	}}

	MPI_Finalize();
	return EXIT_SUCCESS;
}
