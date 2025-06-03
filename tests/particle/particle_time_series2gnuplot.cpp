/*
Program for plotting time series of particle output of PAMHD with gnuplot.

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
#include "limits"
#include "string"
#include "tuple"
#include "unordered_map"
#include "vector"

#include "boost/optional.hpp"
#include "boost/program_options.hpp"
#include "dccrg_cartesian_geometry.hpp"
#include "dccrg_mapping.hpp"
#include "dccrg_topology.hpp"
#include "mpi.h" // must be included before gensimcell
#include "Eigen/Core" // must be included before gensimcell
#include "gensimcell.hpp"
#include "prettyprint.hpp"

#include "common_functions.hpp"
#include "particle/save.hpp"
#include "particle/variables.hpp"

using namespace std;
using namespace pamhd::particle;

/*
Reads simulation data from given file.

Fills out grid info and simulation data withing given volume.

On success returns vacuum permeability.
*/
boost::optional<std::array<double, 4>> read_data(
	const Eigen::Vector3d& volume_start,
	const Eigen::Vector3d& volume_end,
	dccrg::Mapping& cell_id_mapping,
	dccrg::Grid_Topology& topology,
	dccrg::Cartesian_Geometry& geometry,
	unordered_map<uint64_t, Cell_test_particle>& simulation_data,
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

	// read simulation parameters
	std::array<double, 4> metadata;
	MPI_File_read_at(
		file,
		offset,
		(void*) metadata.data(),
		metadata.size(),
		MPI_DOUBLE,
		MPI_STATUS_IGNORE
	);
	offset += sizeof(double) * metadata.size();

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

		const auto
			cell_start = geometry.get_min(cell_id),
			cell_end = geometry.get_max(cell_id);

		// skip data outside of volume to be plotted
		if (
			cell_start[0] > volume_end[0]
			or cell_start[1] > volume_end[1]
			or cell_start[2] > volume_end[2]
			or cell_end[0] < volume_start[0]
			or cell_end[1] < volume_start[1]
			or cell_end[2] < volume_start[2]
		) {
			continue;
		}

		simulation_data[cell_id];
		auto& cell_data = simulation_data.at(cell_id);

		// read fixed size data first
		void* memory_address = NULL;
		int memory_count = -1;
		MPI_Datatype
			memory_datatype = MPI_DATATYPE_NULL,
			file_datatype = MPI_DATATYPE_NULL;

		Cell_test_particle::set_transfer_all(
			true,
			Electric_Field(),
			pamhd::Magnetic_Field(),
			pamhd::Electric_Current_Density(),
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

		Cell_test_particle::set_transfer_all(
			false,
			Electric_Field(),
			pamhd::Magnetic_Field(),
			pamhd::Electric_Current_Density(),
			Nr_Particles_Internal()
		);
		Cell_test_particle::set_transfer_all(true, Particles_Internal());
		tie(
			memory_address,
			memory_count,
			memory_datatype
		) = simulation_data.at(cell_id).get_mpi_datatype();
		Cell_test_particle::set_transfer_all(false, Particles_Internal());

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


int plot(
	const std::string& output_file_name,
	const std::unordered_map<
		Particle_ID::data_type,
		std::vector<std::pair<double, double>>
	>& plot_data,
	const std::string& plot_command_prefix,
	const std::string& xlabel,
	const std::string& ylabel
) {
	const std::string gnuplot_file_name(output_file_name + ".gnuplot");
	std::ofstream gnuplot_file(gnuplot_file_name);

	gnuplot_file
		<< plot_command_prefix
		<< "\nset output '" << output_file_name
		<< "'\nset title '..." << output_file_name.substr(output_file_name.size() - 47, output_file_name.size())
		<< "' noenhanced\nset xlabel '" << xlabel
		<< "'\nset ylabel '" << ylabel
		<< "\nplot '-' using 1:2 with line title ''";
	for (size_t i = 1; i < plot_data.size(); i++) {
		gnuplot_file << ", '-' using 1:2 with line title ''";
	}
	gnuplot_file << "\n";

	gnuplot_file.precision(20);
	for (const auto& item: plot_data) {
		for (const auto& value: item.second) {
			gnuplot_file << value.first << " " << value.second << "\n";
		}
		gnuplot_file << "end\n";
	}

	if (not gnuplot_file.good()) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Writing gnuplot file probably failed..."
			<< std::endl;
	}

	gnuplot_file.close();

	return system(("gnuplot '" + gnuplot_file_name + "'").c_str());
}


int main(int argc, char* argv[])
{
	using std::to_string;

	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
		cerr << "Coudln't initialize MPI." << endl;
		abort();
	}

	MPI_Comm comm = MPI_COMM_WORLD;

	int rank = 0, comm_size = 0;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &comm_size);

	// program options
	bool relative = false;
	std::vector<std::string> input_files;
	std::string
		vertical_variable("r"),
		horizontal_variable("fileno"),
		reference_file(""),
		plot_command_prefix(
			"set term png enhanced size 800, 600"
		);
	constexpr double inf = std::numeric_limits<double>::infinity();
	double
		r_mag_start = 0,
		r_mag_end = inf,
		v_mag_start = 0,
		v_mag_end = inf;
	Eigen::Vector3d
		r_start{-inf, -inf, -inf},
		r_end{inf, inf, inf},
		v_start{-inf, -inf, -inf},
		v_end{inf, inf, inf};

	boost::program_options::options_description
		options("Usage: program_name [options], where options are");
	options.add_options()
		("help", "Print this help message")
		("input-file",
			boost::program_options::value<>(&input_files)
				->default_value(input_files),
			"Results to plot, the --input-file part can be omitted")
		("reference-file",
			boost::program_options::value<>(&reference_file)
				->default_value(reference_file),
			"Use arg as reference result (1st given input file used if empty)")
		("vertical-variable",
			boost::program_options::value<>(&vertical_variable)
				->default_value(vertical_variable),
			"Plot variable arg as a function horizontal variable for "
			"selected particles (one of r, rx, ry, rz, v, vx, vy, vz, E, B)")
		("horizontal-variable",
			boost::program_options::value<>(&horizontal_variable)
				->default_value(horizontal_variable),
			"Plot vertical variable as a function arg for selected particles "
			"(one of fileno, r, rx, ry, rz, v, vx, vy, vz)")
		("plot-command-prefix",
			boost::program_options::value<>(&plot_command_prefix)
				->default_value(plot_command_prefix),
			"Prepend arg to plot command given to gnuplot, "
			"last newline is added automatically")
		("relative", "Plot variable using initial value as the unit of magnitude")
		("r-start",
			boost::program_options::value<>(&r_mag_start)
				->default_value(r_mag_start),
			"Select particles with distance from origin > arg")
		("rx-start",
			boost::program_options::value<>(&r_start[0])
				->default_value(r_start[0]),
			"Select particles with x coordinate of position > arg")
		("ry-start",
			boost::program_options::value<>(&r_start[1])
				->default_value(r_start[1]),
			"Select particles with y coordinate of position > arg")
		("rz-start",
			boost::program_options::value<>(&r_start[2])
				->default_value(r_start[2]),
			"Select particles with z coordinate of position > arg")
		("r-end",
			boost::program_options::value<>(&r_mag_end)
				->default_value(r_mag_end),
			"Select particles with distance from origin < arg")
		("rx-end",
			boost::program_options::value<>(&r_end[0])
				->default_value(r_end[0]),
			"Select particles with x coordinate of position < arg")
		("ry-end",
			boost::program_options::value<>(&r_end[1])
				->default_value(r_end[1]),
			"Select particles with y coordinate of position < arg")
		("rz-end",
			boost::program_options::value<>(&r_end[2])
				->default_value(r_end[2]),
			"Select particles with z coordinate of position < arg")
		("v-start",
			boost::program_options::value<>(&v_mag_start)
				->default_value(v_mag_start),
			"Select particles with velocity > arg")
		("vx-start",
			boost::program_options::value<>(&v_start[0])
				->default_value(v_start[0]),
			"Select particles with x coordinate of velocity > arg")
		("vy-start",
			boost::program_options::value<>(&v_start[1])
				->default_value(v_start[1]),
			"Select particles with y coordinate of velocity > arg")
		("vz-start",
			boost::program_options::value<>(&v_start[2])
				->default_value(v_start[2]),
			"Select particles with z coordinate of velocity > arg")
		("v-end",
			boost::program_options::value<>(&v_mag_end)
				->default_value(v_mag_end),
			"Select particles with velocity < arg")
		("vx-end",
			boost::program_options::value<>(&v_end[0])
				->default_value(v_end[0]),
			"Select particles with x coordinate of velocity < arg")
		("vy-end",
			boost::program_options::value<>(&v_end[1])
				->default_value(v_end[1]),
			"Select particles with y coordinate of velocity < arg")
		("vz-end",
			boost::program_options::value<>(&v_end[2])
				->default_value(v_end[2]),
			"Select particles with z coordinate of velocity < arg");

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
			cout << options <<
				"Particle selections are based only on the first given input file."
				<< endl;
		}
		MPI_Finalize();
		return EXIT_SUCCESS;
	}

	if (
		vertical_variable != "r"
		and vertical_variable != "rx"
		and vertical_variable != "ry"
		and vertical_variable != "rz"
		and vertical_variable != "v"
		and vertical_variable != "vx"
		and vertical_variable != "vy"
		and vertical_variable != "vz"
		and vertical_variable != "E"
		and vertical_variable != "B"
	) {
		if (rank == 0) {
			std::cerr << "Unsupported plot vertical variable: " << vertical_variable << std::endl;
		}
		MPI_Finalize();
		return EXIT_FAILURE;
	}


	if (var_map.count("relative") > 0) {
		relative = true;
	}

	if (rank != 0) {
		MPI_Finalize();
		return EXIT_SUCCESS;
	}

	if (
		input_files.size() == 0
		or (reference_file == "" and input_files.size() == 1)
	) {
		std::cout << "At least two files required."
			<< std::endl;
		MPI_Finalize();
		return EXIT_FAILURE;
	}

	if (reference_file == "") {
		reference_file = input_files[0];
	}

	std::unordered_map<
		Particle_ID::data_type,
		std::vector<std::pair<double, double>> // horizontal, vertical variable
	> plot_data;

	/*
	Fill plot data with selected particles from reference file
	*/

	dccrg::Mapping cell_id_mapping;
	dccrg::Grid_Topology topology;
	dccrg::Cartesian_Geometry geometry(cell_id_mapping.length, cell_id_mapping, topology);
	unordered_map<uint64_t, Cell_test_particle> simulation_data;

	boost::optional<std::array<double, 4>> metadata = read_data(
		r_start,
		r_end,
		cell_id_mapping,
		topology,
		geometry,
		simulation_data,
		reference_file,
		rank
	);
	if (not metadata) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Couldn't read simulation data from file " << input_files[0]
			<< std::endl;
		MPI_Finalize();
		return EXIT_FAILURE;
	}

	for (const auto& item: simulation_data) {
		const auto& cell_data = item.second;

		for (const auto& particle: cell_data[Particles_Internal()]) {

			const auto& pos = particle[Position()];
			if (
				pamhd::norm(pos) < r_mag_start
				or pamhd::norm(pos) > r_mag_end
				or pos[0] < r_start[0]
				or pos[1] < r_start[1]
				or pos[2] < r_start[2]
				or pos[0] > r_end[0]
				or pos[1] > r_end[1]
				or pos[2] > r_end[2]
			) {
				continue;
			}

			const auto& vel = particle[Velocity()];
			if (
				pamhd::norm(vel) < v_mag_start
				or pamhd::norm(vel) > v_mag_end
				or vel[0] < v_start[0]
				or vel[1] < v_start[1]
				or vel[2] < v_start[2]
				or vel[0] > v_end[0]
				or vel[1] > v_end[1]
				or vel[2] > v_end[2]
			) {
				continue;
			}

			double vertical_value = 0, horizontal_value = 0;
			if (vertical_variable == "r") {
				vertical_value = pamhd::norm(pos);
			} else if (vertical_variable == "rx") {
				vertical_value = pos[0];
			} else if (vertical_variable == "ry") {
				vertical_value = pos[1];
			} else if (vertical_variable == "rz") {
				vertical_value = pos[2];
			} else if (vertical_variable == "v") {
				vertical_value = pamhd::norm(vel);
			} else if (vertical_variable == "vx") {
				vertical_value = vel[0];
			} else if (vertical_variable == "vy") {
				vertical_value = vel[1];
			} else if (vertical_variable == "vz") {
				vertical_value = vel[2];
			} else if (vertical_variable == "E") {
				vertical_value = pamhd::norm(cell_data[Electric_Field()]);
			} else if (vertical_variable == "B") {
				vertical_value = pamhd::norm(cell_data[pamhd::Magnetic_Field()]);
			}
			if (horizontal_variable == "fileno") {
				horizontal_value = 0;
			} else if (horizontal_variable == "r") {
				horizontal_value = pamhd::norm(pos);
			} else if (horizontal_variable == "rx") {
				horizontal_value = pos[0];
			} else if (horizontal_variable == "ry") {
				horizontal_value = pos[1];
			} else if (horizontal_variable == "rz") {
				horizontal_value = pos[2];
			} else if (horizontal_variable == "v") {
				horizontal_value = pamhd::norm(vel);
			} else if (horizontal_variable == "vx") {
				horizontal_value = vel[0];
			} else if (horizontal_variable == "vy") {
				horizontal_value = vel[1];
			} else if (horizontal_variable == "vz") {
				horizontal_value = vel[2];
			}

			const auto id = particle[Particle_ID()];
			plot_data[id].emplace_back(horizontal_value, vertical_value);
		}
	}

	/*
	Append to data of selected particles from rest of given input files
	*/
	for (size_t i = 1; i < input_files.size(); i++) {
		if (input_files[i] == reference_file) {
			continue;
		}

		simulation_data.clear();
		boost::optional<std::array<double, 4>> metadata = read_data(
			r_start,
			r_end,
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
			MPI_Finalize();
			return EXIT_FAILURE;
		}

		for (const auto& item: simulation_data) {
			const auto& cell_data = item.second;

			for (const auto& particle: cell_data[Particles_Internal()]) {
				const auto id = particle[Particle_ID()];

				if (plot_data.count(id) == 0) {
					continue;
				}

				const auto& pos = particle[Position()];
				const auto& vel = particle[Velocity()];

				double vertical_value = 0, horizontal_value = 0;
				if (vertical_variable == "r") {
					vertical_value = pamhd::norm(pos);
				} else if (vertical_variable == "rx") {
					vertical_value = pos[0];
				} else if (vertical_variable == "ry") {
					vertical_value = pos[1];
				} else if (vertical_variable == "rz") {
					vertical_value = pos[2];
				} else if (vertical_variable == "v") {
					vertical_value = pamhd::norm(vel);
				} else if (vertical_variable == "vx") {
					vertical_value = vel[0];
				} else if (vertical_variable == "vy") {
					vertical_value = vel[1];
				} else if (vertical_variable == "vz") {
					vertical_value = vel[2];
				} else if (vertical_variable == "E") {
					vertical_value = pamhd::norm(cell_data[Electric_Field()]);
				} else if (vertical_variable == "B") {
					vertical_value = pamhd::norm(cell_data[pamhd::Magnetic_Field()]);
				}
				if (horizontal_variable == "fileno") {
					horizontal_value = double(i);
				} else if (horizontal_variable == "r") {
					horizontal_value = pamhd::norm(pos);
				} else if (horizontal_variable == "rx") {
					horizontal_value = pos[0];
				} else if (horizontal_variable == "ry") {
					horizontal_value = pos[1];
				} else if (horizontal_variable == "rz") {
					horizontal_value = pos[2];
				} else if (horizontal_variable == "v") {
					horizontal_value = pamhd::norm(vel);
				} else if (horizontal_variable == "vx") {
					horizontal_value = vel[0];
				} else if (horizontal_variable == "vy") {
					horizontal_value = vel[1];
				} else if (horizontal_variable == "vz") {
					horizontal_value = vel[2];
				}
				if (relative) {
					vertical_value /= plot_data.at(id)[0].second;
				}
				plot_data.at(id).emplace_back(
					horizontal_value,
					vertical_value
				);
			}
		}
	}

	if (relative) {
		for (auto& item: plot_data) {
			item.second[0].second = 1;
		}
	}

	std::string vertical_variable_min, vertical_variable_max;
	if (vertical_variable == "r") {
		vertical_variable_min += "_" + to_string(r_mag_start);
		vertical_variable_max += "_" + to_string(r_mag_end);
	} else if (vertical_variable == "rx") {
		vertical_variable_min += "_" + to_string(r_start[0]);
		vertical_variable_max += "_" + to_string(r_end[0]);
	} else if (vertical_variable == "ry") {
		vertical_variable_min += "_" + to_string(r_start[1]);
		vertical_variable_max += "_" + to_string(r_end[1]);
	} else if (vertical_variable == "rz") {
		vertical_variable_min += "_" + to_string(r_start[2]);
		vertical_variable_max += "_" + to_string(r_end[2]);
	} else if (vertical_variable == "v") {
		vertical_variable_min += "_" + to_string(v_mag_start);
		vertical_variable_max += "_" + to_string(v_mag_end);
	} else if (vertical_variable == "vx") {
		vertical_variable_min += "_" + to_string(v_start[0]);
		vertical_variable_max += "_" + to_string(v_end[0]);
	} else if (vertical_variable == "vy") {
		vertical_variable_min += "_" + to_string(v_start[1]);
		vertical_variable_max += "_" + to_string(v_end[1]);
	} else if (vertical_variable == "vz") {
		vertical_variable_min += "_" + to_string(v_start[2]);
		vertical_variable_max += "_" + to_string(v_end[2]);
	}

	plot(
		reference_file + "_"
			+ horizontal_variable + "_"
			+ vertical_variable
			+ vertical_variable_min
			+ vertical_variable_max
			+ (relative ? "_rel" : "")
			+ "_trajectories.png",
		plot_data,
		plot_command_prefix,
		horizontal_variable,
		vertical_variable
	);

	MPI_Finalize();

	return EXIT_SUCCESS;
}
