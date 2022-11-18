/*
Extracts subvolumes of MHD output of PAMHD.

Copyright 2016, 2017 Ilja Honkonen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

* Neither the name of copyright holders nor the names of their contributors
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
#include "cstdint"
#include "cstdio"
#include "iostream"
#include "string"
#include "vector"

#include "boost/program_options.hpp"
#include "mpi.h"


using namespace std;


/*
Writes given volume from file named infile_name to outfile_name.
*/
void extract(
	const string& infile_name,
	const string& outfile_name,
	const array<double, 3>& vol_min,
	const array<double, 3>& vol_max
) {
	FILE
		*infile = fopen(infile_name.c_str(), "rb"),
		*outfile = fopen(outfile_name.c_str(), "wb");

	if (infile == nullptr) {
		throw runtime_error("Couldn't open file " + infile_name);
	}
	if (outfile == nullptr) {
		throw runtime_error("Couldn't open file " + outfile_name);
	}

	size_t items_read;

	uint64_t file_version;
	items_read = fread(&file_version, sizeof(file_version), 1, infile);
	if (items_read != 1) { throw runtime_error("Couldn't read file version"); }

	array<double, 4> sim_params;
	items_read = fread(sim_params.data(), sizeof(double), sim_params.size(), infile);
	if (items_read != sim_params.size()) { throw runtime_error("Couldn't read simulation parameters"); }

	uint64_t endianness;
	items_read = fread(&endianness, sizeof(endianness), 1, infile);
	if (items_read != 1) { throw runtime_error("Couldn't read endianness"); }

	array<uint64_t, 3> ref_lvl_0_cells;
	items_read = fread(ref_lvl_0_cells.data(), sizeof(uint64_t), ref_lvl_0_cells.size(), infile);
	if (items_read != ref_lvl_0_cells.size()) { throw runtime_error("Couldn't read grid length in cells"); }

	int max_ref_lvl;
	items_read = fread(&max_ref_lvl, sizeof(max_ref_lvl), 1, infile);
	if (items_read != 1) { throw runtime_error("Couldn't read maximum refinement level"); }
	if (max_ref_lvl > 0) { throw runtime_error("Maximum refinement level > 0 not supported"); }

	unsigned int neigh_len;
	items_read = fread(&neigh_len, sizeof(neigh_len), 1, infile);
	if (items_read != 1) { throw runtime_error("Couldn't read neighborhood length"); }

	array<uint8_t, 3> periodicity;
	items_read = fread(periodicity.data(), sizeof(uint8_t), periodicity.size(), infile);
	if (items_read != periodicity.size()) { throw runtime_error("Couldn't read grid periodicity"); }

	int geom_id;
	items_read = fread(&geom_id, sizeof(geom_id), 1, infile);
	if (items_read != 1) { throw runtime_error("Couldn't read geometry id"); }

	array<double, 3> grid_start;
	items_read = fread(grid_start.data(), sizeof(double), grid_start.size(), infile);
	if (items_read != grid_start.size()) { throw runtime_error("Couldn't read grid starting coordinate"); }

	array<double, 3> lvl_0_cell_len;
	items_read = fread(lvl_0_cell_len.data(), sizeof(double), lvl_0_cell_len.size(), infile);
	if (items_read != lvl_0_cell_len.size()) { throw runtime_error("Couldn't read initial cell length"); }

	uint64_t total_cells;
	items_read = fread(&total_cells, sizeof(total_cells), 1, infile);
	if (items_read != 1) { throw runtime_error("Couldn't read total number of cells"); }

	// TODO: don't read in all at once if not needed
	vector<uint64_t> cell_ids_offsets(2 * total_cells);
	items_read = fread(cell_ids_offsets.data(), sizeof(uint64_t), cell_ids_offsets.size(), infile);
	if (items_read != cell_ids_offsets.size()) { throw runtime_error("Couldn't read cell ids and data offsets"); }

	array<uint64_t, 3>
		min_index{
			numeric_limits<uint64_t>::max(),
			numeric_limits<uint64_t>::max(),
			numeric_limits<uint64_t>::max()
		},
		max_index{
			numeric_limits<uint64_t>::min(),
			numeric_limits<uint64_t>::min(),
			numeric_limits<uint64_t>::min()
		};

	vector<uint64_t> out_cells;
	using cell_data_t = array<uint8_t, sizeof(double) * 24 + sizeof(int) * 2>;
	vector<cell_data_t> out_data;
	for (size_t i = 0; i < cell_ids_offsets.size(); i += 2) {
		auto cell_id = cell_ids_offsets[i];

		cell_id--;
		const array<uint64_t, 3> cell_index{
			cell_id % ref_lvl_0_cells[0],
			cell_id / ref_lvl_0_cells[0] % ref_lvl_0_cells[1],
			cell_id / (ref_lvl_0_cells[0] * ref_lvl_0_cells[1])
		};
		cell_id++;

		const array<double, 3> cell_center{
			grid_start[0] + lvl_0_cell_len[0] * (0.5 + cell_index[0]),
			grid_start[1] + lvl_0_cell_len[1] * (0.5 + cell_index[1]),
			grid_start[2] + lvl_0_cell_len[2] * (0.5 + cell_index[2])
		};

		if (
			cell_center[0] < vol_min[0]
			or cell_center[0] > vol_max[0]
			or cell_center[1] < vol_min[1]
			or cell_center[1] > vol_max[1]
			or cell_center[2] < vol_min[2]
			or cell_center[2] > vol_max[2]
		) {
			continue;
		}

		for (size_t dim = 0; dim < 3; dim++) {
			min_index[dim] = min(min_index[dim], cell_index[dim]);
			max_index[dim] = max(max_index[dim], cell_index[dim]);
		}

		out_cells.push_back(cell_id);

		const auto offset = cell_ids_offsets[i + 1];
		fseek(infile, offset, SEEK_SET);

		cell_data_t data;
		items_read = fread(data.data(), sizeof(cell_data_t::value_type), data.size(), infile);
		if (items_read != data.size()) {
			throw runtime_error(
				"Couldn't read data of cell " + to_string(cell_id)
				+ " at offset " + to_string(offset)
			);
		}
		out_data.push_back(data);
	}

	size_t items_written;

	items_written = fwrite(&file_version, sizeof(file_version), 1, outfile);
	if (items_written != 1) { throw runtime_error("Couldn't write file version"); }

	items_written = fwrite(sim_params.data(), sizeof(double), sim_params.size(), outfile);
	if (items_written != sim_params.size()) { throw runtime_error("Couldn't write simulation parameters"); }

	items_written = fwrite(&endianness, sizeof(endianness), 1, outfile);
	if (items_written != 1) { throw runtime_error("Couldn't write endianness"); }

	const decltype(ref_lvl_0_cells) new_lvl_0_cells{
		max_index[0] - min_index[0] + 1,
		max_index[1] - min_index[1] + 1,
		max_index[2] - min_index[2] + 1
	};
	items_written = fwrite(new_lvl_0_cells.data(), sizeof(uint64_t), new_lvl_0_cells.size(), outfile);
	if (items_written != new_lvl_0_cells.size()) { throw runtime_error("Couldn't write grid length in cells"); }

	items_written = fwrite(&max_ref_lvl, sizeof(max_ref_lvl), 1, outfile);
	if (items_written != 1) { throw runtime_error("Couldn't write maximum refinement level"); }
	if (max_ref_lvl > 0) { throw runtime_error("Maximum refinement level > 0 not supported"); }

	items_written = fwrite(&neigh_len, sizeof(neigh_len), 1, outfile);
	if (items_written != 1) { throw runtime_error("Couldn't write neighborhood length"); }

	items_written = fwrite(periodicity.data(), sizeof(uint8_t), periodicity.size(), outfile);
	if (items_written != periodicity.size()) { throw runtime_error("Couldn't write grid periodicity"); }

	items_written = fwrite(&geom_id, sizeof(geom_id), 1, outfile);
	if (items_written != 1) { throw runtime_error("Couldn't write geometry id"); }

	const decltype(grid_start) new_grid_start{
		grid_start[0] + lvl_0_cell_len[0] * min_index[0],
		grid_start[1] + lvl_0_cell_len[1] * min_index[1],
		grid_start[2] + lvl_0_cell_len[2] * min_index[2]
	};
	items_written = fwrite(new_grid_start.data(), sizeof(double), new_grid_start.size(), outfile);
	if (items_written != new_grid_start.size()) { throw runtime_error("Couldn't write grid starting coordinate"); }

	items_written = fwrite(lvl_0_cell_len.data(), sizeof(double), lvl_0_cell_len.size(), outfile);
	if (items_written != lvl_0_cell_len.size()) { throw runtime_error("Couldn't write initial cell length"); }

	uint64_t new_total_cells = out_cells.size();
	items_written = fwrite(&new_total_cells, sizeof(new_total_cells), 1, outfile);
	if (items_written != 1) { throw runtime_error("Couldn't write total number of cells"); }

	const auto cell_data_start = ftell(outfile);
	vector<uint64_t> new_cell_ids_offsets;
	for (size_t i = 0; i < out_cells.size(); i++) {
		const auto old = out_cells[i] - 1;

		const array<uint64_t, 3>
			old_index{
				old % ref_lvl_0_cells[0],
				old / ref_lvl_0_cells[0] % ref_lvl_0_cells[1],
				old / (ref_lvl_0_cells[0] * ref_lvl_0_cells[1])
			},
			new_index{
				old_index[0] - min_index[0],
				old_index[1] - min_index[1],
				old_index[2] - min_index[2]
			};
		const auto new_id
			= new_index[2] * new_lvl_0_cells[0] * new_lvl_0_cells[1]
			+ new_index[1] * new_lvl_0_cells[0]
			+ new_index[0] + 1;

		new_cell_ids_offsets.push_back(new_id);
		new_cell_ids_offsets.push_back(
			cell_data_start
			+ 2 * out_cells.size() * sizeof(uint64_t)
			+ i * sizeof(cell_data_t::value_type) * tuple_size<cell_data_t>::value
		);
	}

	items_written = fwrite(
		new_cell_ids_offsets.data(),
		sizeof(uint64_t),
		new_cell_ids_offsets.size(),
		outfile
	);
	if (items_written != new_cell_ids_offsets.size()) { throw runtime_error("Couldn't write cell ids and data offsets"); }

	items_written = fwrite(
		out_data.data(),
		tuple_size<cell_data_t>::value * sizeof(cell_data_t::value_type),
		out_data.size(),
		outfile
	);
	if (items_written != out_data.size()) { throw runtime_error("Couldn't write cell data"); }
}


int main(int argc, char* argv[])
{
	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
		cerr << "Coudln't initialize MPI." << endl;
		abort();
	}

	MPI_Comm comm = MPI_COMM_WORLD;

	int rank = 0, comm_size = 0;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &comm_size);

	// program options
	vector<string> input_files;
	double
		min_x = -numeric_limits<double>::infinity(),
		max_x = +numeric_limits<double>::infinity(),
		min_y = -numeric_limits<double>::infinity(),
		max_y = +numeric_limits<double>::infinity(),
		min_z = -numeric_limits<double>::infinity(),
		max_z = +numeric_limits<double>::infinity();

	boost::program_options::options_description
		options("Usage: program_name [options], where options are");
	options.add_options()
		("help", "Print usage instructions")
		("input-file",
			boost::program_options::value<vector<string>>(&input_files),
			"Files to process (in1 out1 in2 out2 ..., the --input-file part can be omitted")
		("min-x",
			boost::program_options::value<double>(&min_x)
				->default_value(min_x),
			"Minimum x coordinate of volume to extract")
		("max-x",
			boost::program_options::value<double>(&max_x)
				->default_value(max_x),
			"Maximum x coordinate of volume to extract")
		("min-y",
			boost::program_options::value<double>(&min_y)
				->default_value(min_y),
			"Minimum y coordinate of volume to extract")
		("max-y",
			boost::program_options::value<double>(&max_y)
				->default_value(max_y),
			"Maximum y coordinate of volume to extract")
		("min-z",
			boost::program_options::value<double>(&min_z)
				->default_value(min_z),
			"Minimum z coordinate of volume to extract")
		("max-z",
			boost::program_options::value<double>(&max_z)
				->default_value(max_z),
			"Maximum z coordinate of volume to extract");

	boost::program_options::positional_options_description positional_options;
	positional_options.add("input-file", -1);

	boost::program_options::variables_map var_map;
	try {
		boost::program_options::store(
			boost::program_options::command_line_parser(argc, argv)
				.options(options)
				.positional(positional_options)
				.run(),
			var_map
		);
	} catch (exception& e) {
		cerr <<  __FILE__ << "(" << __LINE__ << "): "
			<< "Couldn't parse command line options: " << e.what()
			<< endl;
		abort();
	}
	boost::program_options::notify(var_map);

	if (var_map.count("help") > 0) {
		cout << options << endl;
		return EXIT_SUCCESS;
	}

	for (size_t i = 1; i < input_files.size(); i += 2) {
		if (int((i-1)/2) % comm_size != rank) {
			continue;
		}

		extract(
			input_files[i - 1],
			input_files[i],
			{min_x, min_y, min_z},
			{max_x, max_y, max_z}
		);
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
