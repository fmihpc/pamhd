/*
Program for plotting particle statistics from particle program of PAMHD with gnuplot.

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
#include "cmath"
#include "cstdint"
#include "cstdlib"
#include "functional"
#include "fstream"
#include "limits"
#include "string"
#include "tuple"
#include "unordered_map"
#include "vector"

#include "boost/lexical_cast.hpp"
#include "boost/optional.hpp"
#include "boost/program_options.hpp"
#include "dccrg_cartesian_geometry.hpp"
#include "dccrg_mapping.hpp"
#include "dccrg_topology.hpp"
#include "mpi.h" // must be included before gensimcell
#include "Eigen/Core"
#include "gensimcell.hpp"
#include "prettyprint.hpp"

#include "common_functions.hpp"
#include "particle/common.hpp"
#include "particle/save.hpp"
#include "particle/variables.hpp"

using namespace std;
using namespace pamhd::particle;

/*
Reads simulation data from given file.

Fills out grid info and simulation data withing given volume.

On success returns vacuum permeability and Boltzmann constant.
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
	const std::vector<
		std::tuple<
			double, // bin center
			double, // value to plot
			size_t // number of particles in bin
		>
	>& plot_data,
	const std::string& plot_command_prefix,
	const std::string& xlabel,
	const std::string& ylabel
) {
	const std::string
		gnuplot_file_name(output_file_name + ".gnuplot"),
		plot_title(
			output_file_name.substr(
				output_file_name.size() - min(size_t(47), output_file_name.size()),
				output_file_name.size()
			)
		);
	std::ofstream gnuplot_file(gnuplot_file_name);

	gnuplot_file
		<< plot_command_prefix
		<< "\nset output '" << output_file_name
		<< "'\nset title '..." << plot_title
		<< "' noenhanced\nset xlabel '" << xlabel
		<< "'\nset format x '%1.2e'\nset format y '%1.2e'\nset format y2 '%1.2e'"
		<< "\nset ylabel '" << ylabel
		<< "\nset y2label 'Number of particles'\n"
			"set y2tics\nset ytics nomirror\nset y2tics nomirror\n"
			"plot '-' using 1:2 with lines title '', "
			"'' using 1:2 axis x1y2 w l t '#'\n";
	for (const auto& item: plot_data) {
		gnuplot_file << get<0>(item) << " " << get<1>(item) << "\n";
	}
	gnuplot_file << "end\n";
	for (const auto& item: plot_data) {
		gnuplot_file << get<0>(item) << " " << get<2>(item) << "\n";
	}
	gnuplot_file << "end\n";

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
	using std::isinf;

	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
		cerr << "Coudln't initialize MPI." << endl;
		abort();
	}

	MPI_Comm comm = MPI_COMM_WORLD;

	int rank = 0, comm_size = 0;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &comm_size);

	// program options
	std::vector<std::string> input_files;
	std::string
		vertical_variable("N"),
		horizontal_variable("r"),
		plot_command_prefix(
			"set term png enhanced size 800, 600"
		);
	constexpr double inf = std::numeric_limits<double>::infinity();
	double
		r_mag_start = 0,
		r_mag_end = inf,
		boltzmann = 1.3806488e-23,
		charge = 1.602176565e-19;
	Eigen::Vector3d
		r_start{-inf, -inf, -inf},
		r_end{inf, inf, inf};
	size_t horizontal_resolution = 50;

	boost::program_options::options_description
		options("Usage: program_name [options], where options are");
	options.add_options()
		("help", "Print this help message")
		("input-file",
			boost::program_options::value<>(&input_files)
				->default_value(input_files),
			"Results to plot, the --input-file part can be omitted")
		("vertical-variable",
			boost::program_options::value<>(&vertical_variable)
				->default_value(vertical_variable),
			"Plot variable arg as a function horizontal variable for "
			"selected particles [one of N (nr of particles), "
			"V (velocity), P (pressure) or T (temperature)]")
		("horizontal-variable",
			boost::program_options::value<>(&horizontal_variable)
				->default_value(horizontal_variable),
			"Plot vertical variable as a function arg for selected particles "
			"(one of r, rx, ry, rz)")
		("horizontal-resolution",
			boost::program_options::value<>(&horizontal_resolution)
				->default_value(horizontal_resolution),
			"Gather particles in arg bins of selected horizontal variable for plotting")
		("plot-command-prefix",
			boost::program_options::value<>(&plot_command_prefix)
				->default_value(plot_command_prefix),
			"Prepend arg to plot command given to gnuplot, "
			"last newline is added automatically")
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
		("species-charge",
			boost::program_options::value<>(&charge)
				->default_value(charge),
			"Assume every particle has a charge arg")
		("boltzmann",
			boost::program_options::value<>(&boltzmann)
				->default_value(boltzmann),
			"Assume Boltzmann constant has value arg");

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
		vertical_variable != "N"
		and vertical_variable != "V"
		and vertical_variable != "T"
		and vertical_variable != "P"
	) {
		if (rank == 0) {
			std::cerr << "Unsupported vertical variable for plot: " << vertical_variable << std::endl;
		}
		MPI_Finalize();
		return EXIT_FAILURE;
	}
	if (
		vertical_variable == "P"
		and horizontal_variable == "r"
	) {
		if (rank == 0) {
			std::cerr << "Pressure is not supported as a vertical variable when plotting as a function of r"
				<< std::endl;
		}
		MPI_Finalize();
		return EXIT_FAILURE;
	}


	if (rank != 0) {
		MPI_Finalize();
		return EXIT_SUCCESS;
	}

	if (input_files.size() == 0) {
		MPI_Finalize();
		return EXIT_SUCCESS;
	}


	dccrg::Mapping cell_id_mapping;
	dccrg::Grid_Topology topology;
	dccrg::Cartesian_Geometry geometry(
		cell_id_mapping.length,
		cell_id_mapping,
		topology
	);
	unordered_map<uint64_t, Cell_test_particle> simulation_data;
	boost::optional<std::array<double, 4>> metadata = read_data(
		r_start,
		r_end,
		cell_id_mapping,
		topology,
		geometry,
		simulation_data,
		input_files[0],
		rank
	);
	if (not metadata) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Couldn't read simulation data from file " << input_files[0]
			<< std::endl;
		MPI_Finalize();
		return EXIT_FAILURE;
	}

	const auto
		grid_start = geometry.get_start(),
		grid_end = geometry.get_end();
	decltype(grid_start)
		grid_length{{
			grid_end[0] - grid_start[0],
			grid_end[1] - grid_start[1],
			grid_end[2] - grid_start[2]
		}};

	if (isinf(r_start[0])) {
		r_start[0] = grid_start[0];
	}
	if (isinf(r_start[1])) {
		r_start[1] = grid_start[1];
	}
	if (isinf(r_start[2])) {
		r_start[2] = grid_start[2];
	}
	if (isinf(r_end[0])) {
		r_end[0] = grid_end[0];
	}
	if (isinf(r_end[1])) {
		r_end[1] = grid_end[1];
	}
	if (isinf(r_end[2])) {
		r_end[2] = grid_end[2];
	}
	if (isinf(r_mag_end)) {
		r_mag_end =
			max(pow(grid_start[0], 2)+pow(grid_start[1], 2)+pow(grid_start[2], 2),
			max(pow(grid_start[0], 2)+pow(grid_start[1], 2)+pow(  grid_end[2], 2),
			max(pow(grid_start[0], 2)+pow(  grid_end[1], 2)+pow(grid_start[2], 2),
			max(pow(grid_start[0], 2)+pow(  grid_end[1], 2)+pow(  grid_end[2], 2),
			max(pow(  grid_end[0], 2)+pow(grid_start[1], 2)+pow(grid_start[2], 2),
			max(pow(  grid_end[0], 2)+pow(grid_start[1], 2)+pow(  grid_end[2], 2),
			max(pow(  grid_end[0], 2)+pow(  grid_end[1], 2)+pow(grid_start[2], 2),
			    pow(  grid_end[0], 2)+pow(  grid_end[1], 2)+pow(  grid_end[2], 2)
			)))))));
		r_mag_end = sqrt(r_mag_end);
	}

	/*
	Gather data for plotting
	*/

	double bin_start = 0, bin_length = 0;
	if (horizontal_variable == "r") {
		bin_start = r_mag_start;
		bin_length = (r_mag_end - r_mag_start) / horizontal_resolution;
	} else if (horizontal_variable == "rx") {
		bin_start = r_start[0];
		bin_length = (r_end[0] - r_start[0]) / horizontal_resolution;
	} else if (horizontal_variable == "ry") {
		bin_start = r_start[1];
		bin_length = (r_end[1] - r_start[1]) / horizontal_resolution;
	} else if (horizontal_variable == "rz") {
		bin_start = r_start[2];
		bin_length = (r_end[2] - r_start[2]) / horizontal_resolution;
	}

	// data of each bin, from neg to pos direction
	std::vector<
		std::unordered_map<
			Particle_ID::data_type, // id of particle contributing to this bin
			Particle_Internal
		>
	> binned_data(horizontal_resolution);

	// assign particle data to bins
	for (size_t i = 0; i < input_files.size(); i++) {
		if (i > 0) {
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
		}

		for (const auto& item: simulation_data) {
			const auto& cell_data = item.second;

			for (const auto& particle: cell_data[Particles_Internal()]) {

				const auto& pos = particle[Position()];
				double horizontal_value = 0;
				if (horizontal_variable == "r") {
					horizontal_value = pamhd::norm(pos);
				} else if (horizontal_variable == "rx") {
					horizontal_value = pos[0];
				} else if (horizontal_variable == "ry") {
					horizontal_value = pos[1];
				} else if (horizontal_variable == "rz") {
					horizontal_value = pos[2];
				}
				const size_t bin_i = min(
					horizontal_resolution - 1,
					size_t(floor((horizontal_value - bin_start) / bin_length))
				);

				const auto& id = particle[Particle_ID()];
				if (binned_data[bin_i].count(id) == 0) {
					binned_data[bin_i][id] = particle;
				}
				// TODO: replace existing particle if one closer to
				// center of bin, or calculate an average, or...
			}
		}
	}

	std::vector<
		std::tuple<
			double, // bin center
			double, // value to plot
			size_t // number of particles in bin
		>
	> plot_data(horizontal_resolution);
	for (size_t bin_i = 0; bin_i < binned_data.size(); bin_i++) {

		// copy binned particles into a vector for processing
		std::vector<Particle_Internal> particles;
		typename Species_Mass::data_type species_mass = 0;
		for (const auto& bin_data: binned_data[bin_i]) {
			particles.push_back(bin_data.second);
			species_mass += bin_data.second[Species_Mass()];
		}
		species_mass /= particles.size();

		get<0>(plot_data[bin_i]) = bin_start + bin_length * (bin_i + 0.5);
		get<2>(plot_data[bin_i]) = particles.size();

		if (vertical_variable == "N") {
			get<1>(plot_data[bin_i]) = particles.size();
		} else if (vertical_variable == "V") {
			const auto vel = get_bulk_velocity<Mass, Velocity, Species_Mass>(particles);
			get<1>(plot_data[bin_i]) = pamhd::norm(vel);
		} else if (vertical_variable == "T") {
			get<1>(plot_data[bin_i])
				= get_temperature<
					Mass,
					Velocity,
					Species_Mass
				>(
					particles,
					(*metadata)[3]
				);
		} else if (vertical_variable == "P") {
			const auto volume
				= grid_length[0]
				* grid_length[1]
				* grid_length[2]
				/ horizontal_resolution;

			get<1>(plot_data[bin_i])
				= get_pressure<
					Mass,
					Velocity,
					Species_Mass
				>(
					particles,
					(*metadata)[3],
					volume
				);
		}
	}

	plot(
		input_files[0] + "_"
			+ horizontal_variable + "_"
			+ vertical_variable + "_"
			+ "averages.png",
		plot_data,
		plot_command_prefix,
		horizontal_variable,
		vertical_variable
	);

	MPI_Finalize();
	return EXIT_SUCCESS;
}
