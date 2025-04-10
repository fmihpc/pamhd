/*
Tests boundaries class for all simulation variables of PAMHD.

Copyright 2016, 2017 Ilja Honkonen
Copyright 2019, 2025 Finnish Meteorological Institute
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


#include "cstdlib"
#include "iostream"
#include "set"
#include "utility"

#include "dccrg.hpp"
#include "gensimcell.hpp"
#include "mpi.h"

#include "boundaries/geometries.hpp"
#include "boundaries/multivariable_boundaries.hpp"
#include "variable_getter.hpp"


using namespace pamhd::boundaries;


// storage for classifier algorithm
struct Cell_Type {
	using data_type = int;
	static const std::string get_name(){ return {"cell type"}; }
};
const auto Cell_Type_Getter = pamhd::Variable_Getter<Cell_Type>();

struct Mass_Density {
	using data_type = int;
	static const std::string get_name(){ return {"mass density"}; }
	static const std::string get_option_name(){ return {"mass-density"}; }
};
struct Momentum_Density {
	using data_type = std::array<int, 2>;
	static const std::string get_name(){ return {"momentum density"}; }
	static const std::string get_option_name(){ return {"momentum-density"}; }
};

using Cell = gensimcell::Cell<
	gensimcell::Optional_Transfer,
	Cell_Type, Mass_Density, Momentum_Density
>;


int main(int argc, char* argv[])
{
	const char json[] = "{"
		"\"geometries\": ["
			// 1 layer of cells on negative x side of simulation volume
			"{\"box\": {"
				"\"start\": [-9, -9, -9],"
				"\"end\": [0.8, 9, 9]"
			"}},"
			// close to center of sim volume
			"{\"sphere\": {"
				"\"center\": [2.5, 2.5, 2.5],"
				"\"radius\": 0.4"
			"}},"
			// positive z side of sim vol
			"{\"box\": {"
				"\"start\": [-9, -9, 3.2],"
				"\"end\": [9, 9, 9]"
			"}},"
			// negative y side of simulation volume
			"{\"box\": {"
				"\"start\": [-9, -9, -9],"
				"\"end\": [9, 0.8, 9]"
			"}},"
			// positive x side of sim vol
			"{\"box\": {"
				"\"start\": [3.2, -9, -9],"
				"\"end\": [9, 9, 9]"
			"}}"
		"],"
		"\"mass-density\": {"
			"\"value-boundaries\": ["
				"{"
					"\"geometry-id\": 0,"
					"\"time-stamps\": [-3],"
					"\"values\": [2]"
				"},"
				"{"
					"\"geometry-id\": 1,"
					"\"time-stamps\": [1],"
					"\"values\": [\"t\"]"
				"}"
			"],"
			"\"copy-boundaries\": ["
				"{\"geometry-id\": 2}"
			"]"
		"},"
		"\"momentum-density\": {"
			"\"value-boundaries\": ["
				"{"
					"\"geometry-id\": 3,"
					"\"time-stamps\": [-3],"
					"\"values\": [[3, 4]]"
				"},"
				"{"
					"\"geometry-id\": 1,"
					"\"time-stamps\": [1],"
					"\"values\": [\"{t+1, t+2}\"]"
				"}"
			"],"
			"\"copy-boundaries\": ["
				"{\"geometry-id\": 4}"
			"]"
		"}"
	"}";

	rapidjson::Document document;
	document.Parse(json);
	if (document.HasParseError()) {
		std::cerr << "Couldn't parse json data: " << json << std::endl;
		return EXIT_FAILURE;
	}

	Geometries<int, std::array<double, 3>, double, uint64_t> geometries;
	geometries.set(document);

	Multivariable_Boundaries<
		uint64_t,
		int,
		Mass_Density,
		Momentum_Density
	> boundaries;
	boundaries.set(document);

	MPI_Init(&argc, &argv);

	dccrg::Dccrg<Cell> grid; grid
		.set_neighborhood_length(0)
		.set_maximum_refinement_level(0)
		.set_load_balancing_method("RANDOM")
		.set_periodic(false, false, false)
		.set_initial_length({4, 4, 4})
		.initialize(MPI_COMM_WORLD);

	for (const auto& cell: grid.local_cells()) {
		const auto
			start = grid.geometry.get_min(cell.id),
			end = grid.geometry.get_max(cell.id);
		geometries.overlaps(start, end, cell.id);
	}

	// check that cells have been assigned to correct geometries
	const std::set<uint64_t>
		ref_geom0_cells{1, 5, 9, 13, 17, 21, 25, 29, 33, 37, 41, 45, 49, 53, 57, 61},
		ref_geom1_cells{43},
		ref_geom2_cells{49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64},
		ref_geom3_cells{1, 2, 3, 4, 17, 18, 19, 20, 33, 34, 35, 36, 49, 50, 51, 52},
		ref_geom4_cells{4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64};

	const auto geom_ids = geometries.get_geometry_ids();
	if (geom_ids.size() != 5) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< ": Wrong number of geometries: " << geom_ids.size()
			<< ", should be 5."
			<< std::endl;
		return EXIT_FAILURE;
	}


	const auto& geom0_cells = geometries.get_cells(0);
	unsigned int local_geom0_cells_size = geom0_cells.size(), geom0_cells_size = 0;
	MPI_Allreduce(
		&local_geom0_cells_size,
		&geom0_cells_size,
		1,
		MPI_UNSIGNED,
		MPI_SUM,
		MPI_COMM_WORLD
	);

	if (geom0_cells_size != ref_geom0_cells.size()) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< ": Too few cells in geometry 0: " << geom0_cells_size
			<< ", should be " << ref_geom0_cells.size()
			<< std::endl;
		return EXIT_FAILURE;
	}

	for (const auto& cell: geom0_cells) {
		if (ref_geom0_cells.count(cell) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< ": Cell " << cell << " not in geometry 0"
				<< std::endl;
			return EXIT_FAILURE;
		}
	}


	const auto& geom1_cells = geometries.get_cells(1);
	unsigned int local_geom1_cells_size = geom1_cells.size(), geom1_cells_size = 0;
	MPI_Allreduce(
		&local_geom1_cells_size,
		&geom1_cells_size,
		1,
		MPI_UNSIGNED,
		MPI_SUM,
		MPI_COMM_WORLD
	);

	if (geom1_cells_size != ref_geom1_cells.size()) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< ": Wrong number of cells in geometry 1: " << geom1_cells_size
			<< ", should be " << ref_geom1_cells.size()
			<< std::endl;
		return EXIT_FAILURE;
	}

	for (const auto& cell: geom1_cells) {
		if (ref_geom1_cells.count(cell) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< ": Cell " << cell << " not in geometry 1"
				<< std::endl;
			return EXIT_FAILURE;
		}
	}


	const auto& geom2_cells = geometries.get_cells(2);
	unsigned int local_geom2_cells_size = geom2_cells.size(), geom2_cells_size = 0;
	MPI_Allreduce(
		&local_geom2_cells_size,
		&geom2_cells_size,
		1,
		MPI_UNSIGNED,
		MPI_SUM,
		MPI_COMM_WORLD
	);

	if (geom2_cells_size != ref_geom2_cells.size()) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< ": Too few cells in geometry 2: " << geom2_cells_size
			<< ", should be " << ref_geom2_cells.size()
			<< std::endl;
		return EXIT_FAILURE;
	}

	for (const auto& cell: geom2_cells) {
		if (ref_geom2_cells.count(cell) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< ": Cell " << cell << " not in geometry 2"
				<< std::endl;
			return EXIT_FAILURE;
		}
	}


	const auto& geom3_cells = geometries.get_cells(3);
	unsigned int local_geom3_cells_size = geom3_cells.size(), geom3_cells_size = 0;
	MPI_Allreduce(
		&local_geom3_cells_size,
		&geom3_cells_size,
		1,
		MPI_UNSIGNED,
		MPI_SUM,
		MPI_COMM_WORLD
	);

	if (geom3_cells_size != ref_geom3_cells.size()) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< ": Too few cells in geometry 3: " << geom3_cells_size
			<< ", should be " << ref_geom3_cells.size()
			<< std::endl;
		return EXIT_FAILURE;
	}

	for (const auto& cell: geom3_cells) {
		if (ref_geom3_cells.count(cell) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< ": Cell " << cell << " not in geometry 3"
				<< std::endl;
			return EXIT_FAILURE;
		}
	}


	const auto& geom4_cells = geometries.get_cells(4);
	unsigned int local_geom4_cells_size = geom4_cells.size(), geom4_cells_size = 0;
	MPI_Allreduce(
		&local_geom4_cells_size,
		&geom4_cells_size,
		1,
		MPI_UNSIGNED,
		MPI_SUM,
		MPI_COMM_WORLD
	);

	if (geom4_cells_size != ref_geom4_cells.size()) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< ": Too few cells in geometry 4: " << geom4_cells_size
			<< ", should be " << ref_geom4_cells.size()
			<< std::endl;
		return EXIT_FAILURE;
	}

	for (const auto& cell: geom4_cells) {
		if (ref_geom4_cells.count(cell) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< ": Cell " << cell << " not in geometry 4"
				<< std::endl;
			return EXIT_FAILURE;
		}
	}



	boundaries.classify(grid, geometries, Cell_Type_Getter);

	constexpr Mass_Density Mass{};
	constexpr Momentum_Density Momentum{};

	// check that mass density cells have been classified correctly
	unsigned int
		local_value_cells_size = boundaries.get_value_boundary_cells(Mass).size(),
		value_cells_size = 0;
	MPI_Allreduce(
		&local_value_cells_size,
		&value_cells_size,
		1,
		MPI_UNSIGNED,
		MPI_SUM,
		MPI_COMM_WORLD
	);

	if (value_cells_size != 13) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< ": Wrong number of value boundary cells: "
			<< value_cells_size << ", should be 13"
			<< std::endl;
		return EXIT_FAILURE;
	}

	for (const auto& cell: boundaries.get_value_boundary_cells(Mass)) {
		if (
			ref_geom0_cells.count(cell) == 0
			and ref_geom1_cells.count(cell) == 0
		) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< ": Value boundary cell " << cell << " not in geometries 0 or 1."
				<< std::endl;
			return EXIT_FAILURE;
		}
	}


	unsigned int
		local_copy_cells_size = boundaries.get_copy_boundary_cells(Mass).size(),
		copy_cells_size = 0;
	MPI_Allreduce(
		&local_copy_cells_size,
		&copy_cells_size,
		1,
		MPI_UNSIGNED,
		MPI_SUM,
		MPI_COMM_WORLD
	);

	if (copy_cells_size != 11) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< ": Wrong number of copy boundary cells: "
			<< copy_cells_size << ", should be 11"
			<< std::endl;
		return EXIT_FAILURE;
	}

	const std::set<std::vector<uint64_t>> ref_mass_copy_cells{
		{49, 33}, // target, source
		{50, 34}, {51, 35}, {52, 36}, {53, 37}, {54, 38},
		{55, 39}, {56, 40}, {57, 41}, {58, 42},
		//{59, 43}, 43 is value boundary so 59 should be dont_solve
		{60, 44}, {61, 45}, {62, 46}, {63, 47}, {64, 48}
	};

	for (const auto& item: boundaries.get_copy_boundary_cells(Mass)) {
		if (ref_mass_copy_cells.count(item) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< ": Copy boundary item " << item[0] << "<-" << item[1]
				<< " not in copy boundary cells."
				<< std::endl;
			return EXIT_FAILURE;
		}
	}


	unsigned int
		local_dont_cells_size = boundaries.get_dont_solve_cells(Mass).size(),
		dont_cells_size = 0;
	MPI_Allreduce(
		&local_dont_cells_size,
		&dont_cells_size,
		1,
		MPI_UNSIGNED,
		MPI_SUM,
		MPI_COMM_WORLD
	);

	if (dont_cells_size != 5) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< ": Wrong number of dont solve cells: "
			<< dont_cells_size << ", should be 5"
			<< std::endl;
		return EXIT_FAILURE;
	}

	const std::set<uint64_t> ref_dont_solve_cells{49, 53, 57, 59, 61};
	for (const auto& cell: boundaries.get_dont_solve_cells(Mass))
		if (ref_dont_solve_cells.count(cell) == 0) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< ": Invalid dont solve cell: " << cell
			<< std::endl;
		return EXIT_FAILURE;
	}


	// check that momentum density cells have been classified correctly
	local_value_cells_size = boundaries.get_value_boundary_cells(Momentum).size(),
	value_cells_size = 0;
	MPI_Allreduce(
		&local_value_cells_size,
		&value_cells_size,
		1,
		MPI_UNSIGNED,
		MPI_SUM,
		MPI_COMM_WORLD
	);

	if (value_cells_size != 13) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< ": Wrong number of value boundary cells: "
			<< value_cells_size << ", should be 13"
			<< std::endl;
		return EXIT_FAILURE;
	}

	for (const auto& cell: boundaries.get_value_boundary_cells(Momentum)) {
		if (
			ref_geom1_cells.count(cell) == 0
			and ref_geom3_cells.count(cell) == 0
		) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< ": Value boundary cell " << cell << " not in geometries 1 or 3."
				<< std::endl;
			return EXIT_FAILURE;
		}
	}


	local_copy_cells_size = boundaries.get_copy_boundary_cells(Momentum).size(),
	copy_cells_size = 0;
	MPI_Allreduce(
		&local_copy_cells_size,
		&copy_cells_size,
		1,
		MPI_UNSIGNED,
		MPI_SUM,
		MPI_COMM_WORLD
	);

	if (copy_cells_size != 11) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< ": Wrong number of copy boundary cells: "
			<< copy_cells_size << ", should be 11"
			<< std::endl;
		return EXIT_FAILURE;
	}

	const std::set<std::vector<uint64_t>> ref_momentum_copy_cells{
		{8, 7}, {12, 11}, // {target, source}
		{16, 15}, {24, 23}, {28, 27},
		{32, 31}, {40, 39}, {48, 47},
		{56, 55}, {60, 59}, {64, 63}
	};

	for (const auto& item: boundaries.get_copy_boundary_cells(Momentum)) {
		if (ref_momentum_copy_cells.count(item) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< ": Copy boundary item " << item[0] << "<-" << item[1]
				<< " not in copy boundary cells."
				<< std::endl;
			return EXIT_FAILURE;
		}
	}


	MPI_Finalize();
	return EXIT_SUCCESS;
}
