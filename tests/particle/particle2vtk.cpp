/*
Program for converting particle output of PAMHD to vtk format.

Copyright 2015, 2016, 2017 Ilja Honkonen
Copyright 2025 Finnish Meteorological Institute
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


Author(s): Ilja Honkonen
*/

#include "array"
#include "cstdint"
#include "cstdlib"
#include "functional"
#include "fstream"
#include "string"
#include "tuple"
#include "unordered_map"
#include "vector"

#include "boost/filesystem.hpp"
#include "boost/optional.hpp"
#include "boost/program_options.hpp"
#include "dccrg_cartesian_geometry.hpp"
#include "dccrg_mapping.hpp"
#include "dccrg_topology.hpp"
#include "mpi.h" // must be included before gensimcell
#include "Eigen/Core" // must be included before gensimcell
#include "gensimcell.hpp"
//#include "prettyprint.hpp"

#include "common_functions.hpp"
#include "particle/save.hpp"
#include "particle/variables.hpp"

using namespace std;
using namespace pamhd::particle;

using Cell = pamhd::particle::Cell_test_particle;

/*
Reads simulation data from given file.

Fills out grid info and simulation data.

On success returns vacuum permeability.
*/
boost::optional<std::array<double, 4>> read_data(
	dccrg::Mapping& cell_id_mapping,
	dccrg::Grid_Topology& topology,
	dccrg::Cartesian_Geometry& geometry,
	unordered_map<uint64_t, Cell>& simulation_data,
	const std::string& file_name,
	const int mpi_rank
) {
	MPI_File file;
	if (
		MPI_File_open(
			MPI_COMM_SELF,
			const_cast<char*>(file_name.c_str()),
			MPI_MODE_RDONLY,
			MPI_INFO_NULL,
			&file
		) != MPI_SUCCESS
	) {
		cerr << "Process " << mpi_rank
			<< " couldn't open file " << file_name
			<< endl;
		return boost::optional<std::array<double, 4>>();
	}

	MPI_Offset offset = 0;

	// read metadata
	std::array<double, 4> metadata;
	MPI_File_read_at(
		file,
		offset,
		(void*) metadata.data(),
		metadata.size(),
		MPI_DOUBLE,
		MPI_STATUS_IGNORE
	);
	offset += sizeof(double);

	// skip endianness check data
	offset += sizeof(uint64_t);

	if (not cell_id_mapping.read(file, offset)) {
		cerr << "Process " << mpi_rank
			<< " couldn't set cell id mapping from file " << file_name
			<< endl;
		return boost::optional<std::array<double, 4>>();
	}

	offset
		+= cell_id_mapping.data_size()
		+ sizeof(unsigned int)
		+ topology.data_size();

	if (not geometry.read(file, offset)) {
		cerr << "Process " << mpi_rank
			<< " couldn't read geometry from file " << file_name
			<< endl;
		return boost::optional<std::array<double, 4>>();
	}
	offset += geometry.data_size();

	// read number of cells
	uint64_t total_cells = 0;
	MPI_File_read_at(
		file,
		offset,
		&total_cells,
		1,
		MPI_UINT64_T,
		MPI_STATUS_IGNORE
	);
	offset += sizeof(uint64_t);

	if (total_cells == 0) {
		MPI_File_close(&file);
		return boost::optional<std::array<double, 4>>(metadata);
	}

	// read cell ids and data offsets
	vector<pair<uint64_t, uint64_t>> cells_offsets(total_cells);
	MPI_File_read_at(
		file,
		offset,
		cells_offsets.data(),
		2 * total_cells,
		MPI_UINT64_T,
		MPI_STATUS_IGNORE
	);

	// read cell data
	for (const auto& item: cells_offsets) {
		uint64_t
			cell_id = item.first,
			file_address = item.second;

		simulation_data[cell_id];
		auto& cell_data = simulation_data.at(cell_id);

		// read fixed size data first
		void* memory_address = NULL;
		int memory_count = -1;
		MPI_Datatype
			memory_datatype = MPI_DATATYPE_NULL,
			file_datatype = MPI_DATATYPE_NULL;

		Cell::set_transfer_all(
			true,
			Electric_Field(),
			pamhd::Magnetic_Field(),
			Nr_Particles_Internal()
		);
		tie(
			memory_address,
			memory_count,
			memory_datatype
		) = cell_data.get_mpi_datatype();

		int sizeof_memory_datatype;
		MPI_Type_size(memory_datatype, &sizeof_memory_datatype);
		MPI_Type_contiguous(sizeof_memory_datatype, MPI_BYTE, &file_datatype);

		// interpret data from the file using the non-padded type
		MPI_Type_commit(&file_datatype);
		MPI_File_set_view(
			file,
			file_address,
			MPI_BYTE,
			file_datatype,
			const_cast<char*>("native"),
			MPI_INFO_NULL
		);

		MPI_Type_commit(&memory_datatype);
		MPI_File_read_at(
			file,
			0,
			memory_address,
			memory_count,
			memory_datatype,
			MPI_STATUS_IGNORE
		);

		MPI_Type_free(&memory_datatype);
		MPI_Type_free(&file_datatype);

		if (cell_data[Nr_Particles_Internal()] == 0) {
			continue;
		}

		// read particle data
		file_address += size_t(sizeof_memory_datatype);

		cell_data[Particles_Internal()].resize(cell_data[Nr_Particles_Internal()]);

		Cell::set_transfer_all(
			false,
			Electric_Field(),
			pamhd::Magnetic_Field(),
			Nr_Particles_Internal()
		);
		Cell::set_transfer_all(true, Particles_Internal());
		tie(
			memory_address,
			memory_count,
			memory_datatype
		) = simulation_data.at(cell_id).get_mpi_datatype();
		Cell::set_transfer_all(false, Particles_Internal());

		MPI_Type_size(memory_datatype, &sizeof_memory_datatype);
		MPI_Type_contiguous(sizeof_memory_datatype, MPI_BYTE, &file_datatype);

		MPI_Type_commit(&file_datatype);
		MPI_File_set_view(
			file,
			file_address,
			MPI_BYTE,
			file_datatype,
			const_cast<char*>("native"),
			MPI_INFO_NULL
		);

		MPI_Type_commit(&memory_datatype);
		MPI_File_read_at(
			file,
			0,
			memory_address,
			memory_count,
			memory_datatype,
			MPI_STATUS_IGNORE
		);

		MPI_Type_free(&memory_datatype);
		MPI_Type_free(&file_datatype);
	}

	MPI_File_close(&file);

	return boost::optional<std::array<double, 4>>(metadata);
}


/*
Writes given data in vtk format to given file appended with .vtk.
*/
void convert(
	const dccrg::Cartesian_Geometry& geometry,
	const unordered_map<uint64_t, Cell>& simulation_data,
	const std::string& output_file_name_prefix
) {
	std::vector<uint64_t> cells;
	for (const auto& i: simulation_data) {
		const auto cell_id = i.first;
		if (geometry.mapping.get_refinement_level(cell_id) > 0) {
			std::cerr << "Refined mesh not supported." << std::endl;
			abort();
		}
		cells.push_back(cell_id);
	}
	std::sort(cells.begin(), cells.end());

	std::ofstream
		grid_file(output_file_name_prefix + "_grid.vtk"),
		particle_file(output_file_name_prefix + "_particles.vtk");

	grid_file <<
		"# vtk DataFile Version 2.0\n"
		"Grid data from particle results of PAMHD\n"
		"ASCII\n"
		"DATASET STRUCTURED_POINTS\n";
	particle_file <<
		"# vtk DataFile Version 2.0\n"
		"Particle data from PAMHD\n"
		"ASCII\n"
		"DATASET POLYDATA\n";

	const auto& grid_size = geometry.length.get();
	grid_file
		<< "DIMENSIONS "
		<< grid_size[0] + 1 << " "
		<< grid_size[1] + 1 << " "
		<< grid_size[2] + 1 << "\n";

	const auto grid_start = geometry.get_start();
	grid_file
		<< "ORIGIN "
		<< grid_start[0] << " "
		<< grid_start[1] << " "
		<< grid_start[2] << "\n";

	const auto cell_length = geometry.get_level_0_cell_length();
	grid_file
		<< "SPACING "
		<< cell_length[0] << " "
		<< cell_length[1] << " "
		<< cell_length[2] << "\n";

	grid_file << "CELL_DATA " << cells.size() << "\n";

	const pamhd::Magnetic_Field Mag{};
	grid_file << "VECTORS magnetic_field float\n";
	for (const auto& cell: cells) {
		const auto magnetic_field = simulation_data.at(cell)[Mag];
		grid_file
			<< magnetic_field[0] << " "
			<< magnetic_field[1] << " "
			<< magnetic_field[2] << "\n";
	}

	const Electric_Field Ele{};
	grid_file << "VECTORS electric_field float\n";
	for (const auto& cell: cells) {
		const auto electric_field = simulation_data.at(cell)[Ele];
		grid_file
			<< electric_field[0] << " "
			<< electric_field[1] << " "
			<< electric_field[2] << "\n";
	}


	size_t total_particles = 0;
	for (const auto& cell: cells) {
		total_particles
			+= simulation_data.at(cell)[Particles_Internal()].size();
	}

	particle_file << "POINTS " << total_particles << " float\n";
	for (const auto& cell: cells) {
		for (const auto& particle:
			simulation_data.at(cell)[Particles_Internal()]
		) {
			particle_file
				<< particle[Position()][0] << " "
				<< particle[Position()][1] << " "
				<< particle[Position()][2] << "\n";
		}
	}

	particle_file
		<< "VERTICES "
		<< total_particles << " "
		<< 2 * total_particles << "\n";
	for (size_t i = 0; i < total_particles; i++) {
		particle_file << "1 " << i << "\n";
	}

	particle_file << "POINT_DATA " << total_particles << "\n";

	particle_file << "SCALARS id float\nLOOKUP_TABLE default\n";
	for (const auto& cell: cells) {
		for (const auto& particle:
			simulation_data.at(cell)[Particles_Internal()]
		) {
			particle_file << particle[Particle_ID()] << "\n";
		}
	}

	particle_file << "SCALARS velocity_magnitude float\nLOOKUP_TABLE default\n";
	for (const auto& cell: cells) {
		for (const auto& particle:
			simulation_data.at(cell)[Particles_Internal()]
		) {
			particle_file << pamhd::norm(particle[Velocity()]) << "\n";
		}
	}
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
	bool verbose = false;
	std::vector<std::string> input_files;

	boost::program_options::options_description
		options("Usage: program_name [options], where options are");
	options.add_options()
		("help", "Print this help message")
		("verbose", "Print run time information")
		("input-file",
			boost::program_options::value<std::vector<std::string>>(&input_files),
			"Results to plot, the --input-file part can be omitted");

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
	} catch (std::exception& e) {
		if (rank == 0) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
				<< "Couldn't parse command line options: " << e.what()
				<< std::endl;
		}
		abort();
	}
	boost::program_options::notify(var_map);

	if (var_map.count("help") > 0) {
		if (rank == 0) {
			cout
				<<
					"Converts given .dc files from MHD test program to "
					"VTK format and writes the names of all created files "
					"to a mhd.visit file in the directory of the first "
					"given .dc file.\n"
				<< options
				<< endl;
		}
		MPI_Finalize();
		return EXIT_SUCCESS;
	}

	if (input_files.size() == 0) {
		MPI_Finalize();
		return EXIT_SUCCESS;
	}

	if (var_map.count("verbose") > 0) {
		verbose = true;
	}

	if (rank == 0) {
		boost::filesystem::path input_dir(input_files[0]);
		input_dir = input_dir.parent_path();

		std::ofstream
			visit_grid_file(
				(input_dir / boost::filesystem::path("particle_grid.visit")).string()
			),
			visit_particle_file(
				(input_dir / boost::filesystem::path("particle_particles.visit")).string()
			);

		visit_grid_file << "!NBLOCKS 1\n";
		for (size_t i = 0; i < input_files.size(); i++) {
			visit_grid_file << boost::filesystem::path(
				input_files[i].substr(0, input_files[i].size() - 3) + "_grid.vtk\n"
			).filename().string();
		}
		visit_grid_file.close();

		visit_particle_file << "!NBLOCKS 1\n";
		for (size_t i = 0; i < input_files.size(); i++) {
			visit_particle_file << boost::filesystem::path(
				input_files[i].substr(0, input_files[i].size() - 3) + "_particles.vtk\n"
			).filename().string();
		}
		visit_particle_file.close();
	}

	for (size_t i = 0; i < input_files.size(); i++) {
		if (int(i) % comm_size != rank) {
			continue;
		}

		if (verbose) {
			std::cout << "MPI rank " << rank
				<< " processing file " << input_files[i]
				<< std::endl;
		}

		dccrg::Mapping cell_id_mapping;
		dccrg::Grid_Topology topology;
		dccrg::Cartesian_Geometry geometry(cell_id_mapping.length, cell_id_mapping, topology);
		unordered_map<uint64_t, Cell> simulation_data;

		boost::optional<std::array<double, 4>> metadata
			= read_data(
				cell_id_mapping,
				topology,
				geometry,
				simulation_data,
				input_files[i],
				rank
			);
		if (not metadata) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< "Couldn't read simulation data from file " << input_files[i]
				<< std::endl;
			continue;
		}

		convert(
			geometry,
			simulation_data,
			input_files[i].substr(0, input_files[i].size() - 3)
		);
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
