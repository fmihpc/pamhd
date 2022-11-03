/*
Program for plotting MHD output of PAMHD with gnuplot.

Copyright 2014, 2015, 2016, 2017 Ilja Honkonen
Copyright 2022 Finnish Meteorological Institute
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
#include "cstdlib"
#include "functional"
#include "fstream"
#include "string"
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

#include "mhd/common.hpp"
#include "mhd/save.hpp"
#include "mhd/variables.hpp"
#include "variables.hpp"


using namespace std;


/*
Reads simulation data from given file.

Fills out grid info and simulation data.

On success returns simulation time, adiabatic index,
proton mass, vacuum permeability.
Returns uninitialized value on error.
*/
boost::optional<std::array<double, 4>> read_data(
	dccrg::Mapping& cell_id_mapping,
	dccrg::Grid_Topology& topology,
	dccrg::Cartesian_Geometry& geometry,
	unordered_map<uint64_t, pamhd::mhd::Cell>& simulation_data,
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

	// read file version
	uint64_t file_version = 0;
	MPI_File_read_at( // TODO: add error checking
		file,
		offset,
		&file_version,
		1,
		MPI_UINT64_T,
		MPI_STATUS_IGNORE
	);
	offset += sizeof(uint64_t);
	if (file_version != 2) {
		cerr << "Process " << mpi_rank
			<< " Unsupported file version: " << file_version
			<< endl;
		return boost::optional<std::array<double, 4>>();
	}

	// read simulation parameters
	std::array<double, 4> simulation_parameters;
	MPI_File_read_at(
		file,
		offset,
		simulation_parameters.data(),
		4,
		MPI_DOUBLE,
		MPI_STATUS_IGNORE
	);
	offset += 4 * sizeof(double);

	// check endianness
	uint64_t endianness = 0;
	MPI_File_read_at(
		file,
		offset,
		&endianness,
		1,
		MPI_UINT64_T,
		MPI_STATUS_IGNORE
	);
	offset += sizeof(uint64_t);
	if (endianness != 0x1234567890abcdef) {
		cerr << "Process " << mpi_rank
			<< " Unsupported endianness: " << endianness
			<< ", should be " << 0x1234567890abcdef
			<< endl;
		return boost::optional<std::array<double, 4>>();
	}

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
		return boost::optional<std::array<double, 4>>(simulation_parameters);
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
	pamhd::mhd::Cell::set_transfer_all(
		true,
		pamhd::mhd::HD_State_Conservative(),
		pamhd::Electric_Current_Density(),
		pamhd::mhd::Solver_Info(),
		pamhd::MPI_Rank(),
		pamhd::Resistivity(),
		pamhd::Magnetic_Field(),
		pamhd::Face_Magnetic_Field(),
		pamhd::Edge_Electric_Field(),
		pamhd::Bg_Magnetic_Field_Pos_X(),
		pamhd::Bg_Magnetic_Field_Pos_Y(),
		pamhd::Bg_Magnetic_Field_Pos_Z()
	);
	for (const auto& item: cells_offsets) {
		const uint64_t
			cell_id = item.first,
			file_address = item.second;

		simulation_data[cell_id];

		/*
		dccrg writes cell data without padding so store the
		non-padded version of the cell's datatype into file_datatype
		*/
		void* memory_address = NULL;
		int memory_count = -1;
		MPI_Datatype
			memory_datatype = MPI_DATATYPE_NULL,
			file_datatype = MPI_DATATYPE_NULL;

		tie(
			memory_address,
			memory_count,
			memory_datatype
		) = simulation_data.at(cell_id).get_mpi_datatype();

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
	}

	MPI_File_close(&file);

	return boost::optional<std::array<double, 4>>(simulation_parameters);
}


/*
Plots 1d data from given list in that order.
*/
int plot_1d(
	const dccrg::Cartesian_Geometry& geometry,
	const unordered_map<uint64_t, pamhd::mhd::Cell>& simulation_data,
	const std::vector<uint64_t>& cells,
	const std::string& output_file_name_prefix,
	const double simulation_time,
	const double adiabatic_index,
	const double vacuum_permeability,
	const std::string& common_cmd,
	const std::string& density_pressure_cmd,
	const std::string& velocity_cmd,
	const std::string& magnetic_field_cmd,
	const std::string& electric_field_cmd,
	const std::string& current_density_cmd
) {
	using std::to_string;

	const string gnuplot_file_name(output_file_name_prefix + ".dat");
	ofstream gnuplot_file(gnuplot_file_name);

	const size_t tube_dim
		= [&](){
			for (size_t i = 0; i < geometry.length.get().size(); i++) {
				if (geometry.length.get()[i] > 1) {
					return i;
				}
			}
			return size_t(999);
		}();

	const double
		tube_start = geometry.get_start()[tube_dim],
		tube_end = geometry.get_end()[tube_dim];

	// mass density & pressure
	gnuplot_file
		<< common_cmd
		<< "\nset title 'Time " << simulation_time << "'\n"
		<< "set output '"
		<< output_file_name_prefix + ".png"
		<< "'\nset xlabel 'Dimension "
		<< to_string(tube_dim + 1)
		<< "'\nset xrange [" << to_string(tube_start)
		<< " : " << to_string(tube_end)
		<< "]\n" << density_pressure_cmd
		<< "\nplot "
		     "'-' using 1:2 axes x1y1 with line linewidth 2, "
		     "'-' u 1:2 axes x1y2 w l lt 3 lw 2\n";

	for (const auto& cell_id: cells) {
		const auto& hd_data
			= simulation_data.at(cell_id)[pamhd::mhd::HD_State_Conservative()];
		const double x = geometry.get_center(cell_id)[tube_dim];
		gnuplot_file << x << " " << hd_data[pamhd::mhd::Mass_Density()] << "\n";
	}
	gnuplot_file << "end\n";

	for (const auto& cell_id: cells) {
		auto& hd_data = simulation_data.at(cell_id)[pamhd::mhd::HD_State_Conservative()];
		const double
			x = geometry.get_center(cell_id)[tube_dim],
			pressure = [&](){
				try {
					return pamhd::mhd::get_pressure(
						hd_data[pamhd::mhd::Mass_Density()],
						hd_data[pamhd::mhd::Momentum_Density()],
						hd_data[pamhd::mhd::Total_Energy_Density()],
						simulation_data.at(cell_id)[pamhd::Magnetic_Field()],
						adiabatic_index,
						vacuum_permeability
					);
				} catch (std::domain_error& e) {
					return 0.0;
				}
			}();
		gnuplot_file << x << " " << pressure << "\n";
	}
	gnuplot_file << "end\nreset\n";

	// velocity
	gnuplot_file
		<< common_cmd
		<< "\nset output '"
		<< output_file_name_prefix + "_V.png"
		<< "'\n" << velocity_cmd
		<< "\nset key horizontal bottom outside\nplot "
		     "'-' u 1:2 w l lw 2 t 'component 1', "
		     "'-' u 1:2 w l lw 2 t 'component 2', "
		     "'-' u 1:2 w l lw 2 t 'component 3'\n";

	for (const auto& cell_id: cells) {
		const auto& hd_data
			= simulation_data.at(cell_id)[pamhd::mhd::HD_State_Conservative()];
		const auto& m = hd_data[pamhd::mhd::Momentum_Density()];
		const auto& rho = hd_data[pamhd::mhd::Mass_Density()];
		const double x = geometry.get_center(cell_id)[tube_dim];
		if (rho > 0) {
			gnuplot_file << x << " " << pamhd::mhd::get_velocity(m, rho)[0] << "\n";
		} else {
			gnuplot_file << x << " 0\n";
		}
	}
	gnuplot_file << "end\n";

	for (const auto& cell_id: cells) {
		const auto& hd_data
			= simulation_data.at(cell_id)[pamhd::mhd::HD_State_Conservative()];
		const auto& m = hd_data[pamhd::mhd::Momentum_Density()];
		const auto& rho = hd_data[pamhd::mhd::Mass_Density()];
		const double x = geometry.get_center(cell_id)[tube_dim];
		if (rho > 0) {
			gnuplot_file << x << " " << pamhd::mhd::get_velocity(m, rho)[1] << "\n";
		} else {
			gnuplot_file << x << " 0\n";
		}
	}
	gnuplot_file << "end\n";

	for (const auto& cell_id: cells) {
		const auto& hd_data
			= simulation_data.at(cell_id)[pamhd::mhd::HD_State_Conservative()];
		const auto& m = hd_data[pamhd::mhd::Momentum_Density()];
		const auto& rho = hd_data[pamhd::mhd::Mass_Density()];
		const double x = geometry.get_center(cell_id)[tube_dim];
		if (rho > 0) {
			gnuplot_file << x << " " << pamhd::mhd::get_velocity(m, rho)[2] << "\n";
		} else {
			gnuplot_file << x << " 0\n";
		}
	}
	gnuplot_file << "end\nreset\n";

	// magnetic fields
	gnuplot_file
		<< common_cmd
		<< "\nset output '"
		<< output_file_name_prefix + "_B1.png"
		<< "'\n" << magnetic_field_cmd
		<< "\nset key horizontal bottom outside\nplot "
		     "'-' u 1:2 w l lw 2 t 'component 1', "
		     "'-' u 1:2 w l lw 2 t 'component 2', "
		     "'-' u 1:2 w l lw 2 t 'component 3'\n";

	for (const auto& cell_id: cells) {
		const auto& B = simulation_data.at(cell_id)[pamhd::Magnetic_Field()];
		const double x = geometry.get_center(cell_id)[tube_dim];
		gnuplot_file << x << " " << B[0] << "\n";
	}
	gnuplot_file << "end\n";

	for (const auto& cell_id: cells) {
		const auto& B = simulation_data.at(cell_id)[pamhd::Magnetic_Field()];
		const double x = geometry.get_center(cell_id)[tube_dim];
		gnuplot_file << x << " " << B[1] << "\n";
	}
	gnuplot_file << "end\n";

	for (const auto& cell_id: cells) {
		const auto& B = simulation_data.at(cell_id)[pamhd::Magnetic_Field()];
		const double x = geometry.get_center(cell_id)[tube_dim];
		gnuplot_file << x << " " << B[2] << "\n";
	}
	gnuplot_file << "end\nreset\n";

	gnuplot_file
		<< common_cmd
		<< "\nset output '"
		<< output_file_name_prefix + "_B0px.png"
		<< "'\n" << magnetic_field_cmd
		<< "\nset key horizontal bottom outside\nplot "
		     "'-' u 1:2 w l lw 2 t 'component 1', "
		     "'-' u 1:2 w l lw 2 t 'component 2', "
		     "'-' u 1:2 w l lw 2 t 'component 3'\n";

	for (const auto& cell_id: cells) {
		const auto& B = simulation_data.at(cell_id)[pamhd::Bg_Magnetic_Field_Pos_X()];
		const double x = geometry.get_center(cell_id)[tube_dim];
		gnuplot_file << x << " " << B[0] << "\n";
	}
	gnuplot_file << "end\n";

	for (const auto& cell_id: cells) {
		const auto& B = simulation_data.at(cell_id)[pamhd::Bg_Magnetic_Field_Pos_X()];
		const double x = geometry.get_center(cell_id)[tube_dim];
		gnuplot_file << x << " " << B[1] << "\n";
	}
	gnuplot_file << "end\n";

	for (const auto& cell_id: cells) {
		const auto& B = simulation_data.at(cell_id)[pamhd::Bg_Magnetic_Field_Pos_X()];
		const double x = geometry.get_center(cell_id)[tube_dim];
		gnuplot_file << x << " " << B[2] << "\n";
	}
	gnuplot_file << "end\nreset\n";

	// current density
	gnuplot_file
		<< common_cmd
		<< "\nset output '"
		<< output_file_name_prefix + "_J.png"
		<< "'\n" << current_density_cmd
		<< "\nset key horizontal bottom outside\nplot "
		     "'-' u 1:2 w l lw 2 t 'component 1', "
		     "'-' u 1:2 w l lw 2 t 'component 2', "
		     "'-' u 1:2 w l lw 2 t 'component 3'\n";

	for (const auto& cell_id: cells) {
		const auto& J = simulation_data.at(cell_id)[pamhd::Electric_Current_Density()];
		const double x = geometry.get_center(cell_id)[tube_dim];
		gnuplot_file << x << " " << J[0] << "\n";
	}
	gnuplot_file << "end\n";

	for (const auto& cell_id: cells) {
		const auto& J = simulation_data.at(cell_id)[pamhd::Electric_Current_Density()];
		const double x = geometry.get_center(cell_id)[tube_dim];
		gnuplot_file << x << " " << J[1] << "\n";
	}
	gnuplot_file << "end\n";

	for (const auto& cell_id: cells) {
		const auto& J = simulation_data.at(cell_id)[pamhd::Electric_Current_Density()];
		const double x = geometry.get_center(cell_id)[tube_dim];
		gnuplot_file << x << " " << J[2] << "\n";
	}
	gnuplot_file << "end\nreset\n";

	// resistivity & mpi rank
	gnuplot_file
		<< common_cmd
		<< "\nset output '"
		<< output_file_name_prefix + "_R.png"
		<< "'\nset xlabel 'Dimension "
		<< to_string(tube_dim + 1)
		<< "'\nset xrange [" << to_string(tube_start)
		<< " : " << to_string(tube_end)
		<< "]\nset ylabel \"Resistivity\" textcolor lt 1\n"
		   "set y2label \"MPI rank\" textcolor lt 3"
		   "\nunset key\n"
		   "set ytics nomirror\n"
		   "set format x '%.2e'\n"
		   "set format y '%.2e'\n"
		   "set format y2 '%.2e'\n"
		   "set y2tics auto\n"
		   "plot "
		     "'-' using 1:2 axes x1y1 with line linewidth 2, "
		     "'-' u 1:2 axes x1y2 w l lt 3 lw 2\n";

	for (const auto& cell_id: cells) {
		const double x = geometry.get_center(cell_id)[tube_dim];
		gnuplot_file << x << " " << simulation_data.at(cell_id)[pamhd::Resistivity()] << "\n";
	}
	gnuplot_file << "end\n";

	for (const auto& cell_id: cells) {
		const double x = geometry.get_center(cell_id)[tube_dim];
		gnuplot_file
			<< x << " " << simulation_data.at(cell_id)[pamhd::MPI_Rank()] << "\n";
	}
	gnuplot_file << "end\nreset\n";

	// face magnetic field
	gnuplot_file
		<< common_cmd
		<< "\nset output '"
		<< output_file_name_prefix + "_Bf.png"
		<< "'\n" << magnetic_field_cmd
		<< "\nset key horizontal bottom outside\nplot "
		     "'-' u 1:2 w l lw 2 t 'component 1', "
		     "'-' u 1:2 w l lw 2 t 'component 2', "
		     "'-' u 1:2 w l lw 2 t 'component 3'\n";

	for (const auto& cell_id: cells) {
		const auto& B = simulation_data.at(cell_id)[pamhd::Face_Magnetic_Field()];
		const double x = geometry.get_center(cell_id)[tube_dim];
		gnuplot_file << x << " " << B[0] << "\n";
	}
	gnuplot_file << "end\n";

	for (const auto& cell_id: cells) {
		const auto& B = simulation_data.at(cell_id)[pamhd::Face_Magnetic_Field()];
		const double x = geometry.get_center(cell_id)[tube_dim];
		gnuplot_file << x << " " << B[1] << "\n";
	}
	gnuplot_file << "end\n";

	for (const auto& cell_id: cells) {
		const auto& B = simulation_data.at(cell_id)[pamhd::Face_Magnetic_Field()];
		const double x = geometry.get_center(cell_id)[tube_dim];
		gnuplot_file << x << " " << B[2] << "\n";
	}
	gnuplot_file << "end\nreset\n";

	// edge electric field
	gnuplot_file
		<< common_cmd
		<< "\nset output '"
		<< output_file_name_prefix + "_Ee.png"
		<< "'\n" << electric_field_cmd
		<< "\nset key horizontal bottom outside\nplot "
		     "'-' u 1:2 w l lw 2 t 'component 1', "
		     "'-' u 1:2 w l lw 2 t 'component 2', "
		     "'-' u 1:2 w l lw 2 t 'component 3'\n";

	for (const auto& cell_id: cells) {
		const auto& B = simulation_data.at(cell_id)[pamhd::Edge_Electric_Field()];
		const double x = geometry.get_center(cell_id)[tube_dim];
		gnuplot_file << x << " " << B[0] << "\n";
	}
	gnuplot_file << "end\n";

	for (const auto& cell_id: cells) {
		const auto& B = simulation_data.at(cell_id)[pamhd::Edge_Electric_Field()];
		const double x = geometry.get_center(cell_id)[tube_dim];
		gnuplot_file << x << " " << B[1] << "\n";
	}
	gnuplot_file << "end\n";

	for (const auto& cell_id: cells) {
		const auto& B = simulation_data.at(cell_id)[pamhd::Edge_Electric_Field()];
		const double x = geometry.get_center(cell_id)[tube_dim];
		gnuplot_file << x << " " << B[2] << "\n";
	}
	gnuplot_file << "end\nreset\n";


	gnuplot_file.close();

	return system(("gnuplot " + gnuplot_file_name).c_str());
}


/*!
Prints gnuplot commands and data to plot given variable in 2d.
*/
void write_gnuplot_cmd_2d(
	ofstream& gnuplot_file,
	const std::array<double, 3>& grid_geom_start,
	const std::array<double, 3>& grid_geom_length,
	const std::array<uint64_t, 3>& grid_size,
	const std::array<size_t, 2>& dimensions,
	const unordered_map<uint64_t, pamhd::mhd::Cell>& simulation_data,
	const double simulation_time,
	const std::vector<uint64_t>& cells,
	const std::string& common_cmd,
	const std::string& output_file_name_prefix,
	const std::string& output_file_name_suffix,
	const std::string& var_name,
	const std::string& var_cmd,
	const std::function<double(const pamhd::mhd::Cell& cell_data)>& value_for_plot
) {
	const auto
		grid_size1 = grid_size[dimensions[0]],
		grid_size2 = grid_size[dimensions[1]];
	const auto
		grid_start1 = grid_geom_start[dimensions[0]],
		grid_start2 = grid_geom_start[dimensions[1]],
		grid_length1 = grid_geom_length[dimensions[0]],
		grid_length2 = grid_geom_length[dimensions[1]],
		cell_size1 = grid_length1 / grid_size1,
		cell_size2 = grid_length2 / grid_size2;

	gnuplot_file
		<< common_cmd
		<< "set output '"
		<< output_file_name_prefix + "_"
		<< var_name << "." << output_file_name_suffix
		<< "'\nset ylabel 'Dimension "
		<< to_string(dimensions[1] + 1)
		<< "'\nset xlabel 'Dimension "
		<< to_string(dimensions[0] + 1)
		<< ", simulation time " << to_string(simulation_time)
		<< "'\nset xrange[" << to_string(grid_start1) << " : "
		<< to_string(grid_start1 + grid_length1) << "]"
		<< "\nset yrange[" << to_string(grid_start2) << " : "
		<< to_string(grid_start2 + grid_length2) << "]\n"
		<< var_cmd
		<< "\nplot '-' using ($1 / "
		<< grid_size1 << " * " << grid_length1 << " + "
		<< grid_start1 << " + " << cell_size1 / 2 << "):($2 / "
		<< grid_size2 << " * " << grid_length2 << " + "
		<< grid_start2 << " + " << cell_size2 / 2
		<< "):3 matrix with image title ''\n";

	for (size_t i = 0; i < cells.size(); i++) {
		gnuplot_file << value_for_plot(simulation_data.at(cells[i])) << " ";

		if (i % grid_size1 == grid_size1 - 1) {
			gnuplot_file << "\n";
		}
	}
	gnuplot_file << "\nend\nreset\n";
}


/*
Plots data of given cells as a 2d color map.
*/
int plot_2d(
	const dccrg::Cartesian_Geometry& geometry,
	const unordered_map<uint64_t, pamhd::mhd::Cell>& simulation_data,
	const std::vector<uint64_t>& cells,
	const std::string& output_file_name_prefix,
	const std::string& output_file_name_suffix,
	const double simulation_time,
	const double adiabatic_index,
	const double vacuum_permeability,
	const std::string& common_cmd,
	const std::string& mass_density_cmd,
	const std::string& pressure_cmd,
	const std::string& velocity_cmd,
	const std::string& magnetic_field_cmd,
	const std::string& electric_field_cmd,
	const std::string& current_density_cmd,
	const std::string& rank_cmd,
	const std::string& resistivity_cmd,
	const std::string& bdy_cmd
) {
	using std::to_string;

	const auto& grid_size = geometry.length.get();

	// indices of dimensions with more than one cell
	const std::array<size_t, 2> dimensions
		= [&grid_size](){

			std::array<size_t, 2> ret_val{{0, 0}};
			size_t dims = 0;

			for (size_t dim = 0; dim < grid_size.size(); dim++) {
				if (grid_size[dim] > 1) {
					ret_val[dims] = dim;
					dims++;
				}
			}
			if (dims != 2) {
				std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
					<< "Grid with != 2 dimensions given to plot_2d: " << grid_size
					<< std::endl;
				abort();
			}

			return ret_val;
		}();

	const auto
		grid_start = geometry.get_start(),
		grid_end = geometry.get_end();
	const decltype(grid_start) grid_geom_length{{
		grid_end[0] - grid_start[0],
		grid_end[1] - grid_start[1],
		grid_end[2] - grid_start[2]
	}};
	const string gnuplot_file_name(output_file_name_prefix + ".dat");

	std::ofstream gnuplot_file(gnuplot_file_name);

	auto write_gnuplot_cmd_current = std::bind(
		write_gnuplot_cmd_2d,
		std::ref(gnuplot_file),
		grid_start,
		grid_geom_length,
		grid_size,
		dimensions,
		std::cref(simulation_data),
		simulation_time,
		cells,
		common_cmd + "\n",
		output_file_name_prefix,
		output_file_name_suffix,
		std::placeholders::_1,
		std::placeholders::_2,
		std::placeholders::_3
	);

	// mass density
	if (mass_density_cmd != "") {
		write_gnuplot_cmd_current(
			"rho",
			"\n" + mass_density_cmd + "\n",
			[](const pamhd::mhd::Cell& cell_data){
				return cell_data[pamhd::mhd::HD_State_Conservative()][pamhd::mhd::Mass_Density()];
			}
		);
	}

	// pressure
	if (pressure_cmd != "") {
		write_gnuplot_cmd_current(
			"P",
			"\n" + pressure_cmd + "\n",
			[&](const pamhd::mhd::Cell& cell_data){
				try {
					return pamhd::mhd::get_pressure(
						cell_data[pamhd::mhd::HD_State_Conservative()][pamhd::mhd::Mass_Density()],
						cell_data[pamhd::mhd::HD_State_Conservative()][pamhd::mhd::Momentum_Density()],
						cell_data[pamhd::mhd::HD_State_Conservative()][pamhd::mhd::Total_Energy_Density()],
						cell_data[pamhd::Magnetic_Field()],
						adiabatic_index,
						vacuum_permeability
					);
				} catch (std::domain_error& e) {
					return 0.0;
				}
			}
		);
	}

	// velocity
	if (velocity_cmd != "") {
		write_gnuplot_cmd_current(
			"Vx",
			"\n" + velocity_cmd + " 1\"\n",
			[](const pamhd::mhd::Cell& cell_data)->double {
				return pamhd::mhd::get_velocity(
					cell_data[pamhd::mhd::HD_State_Conservative()][pamhd::mhd::Momentum_Density()],
					cell_data[pamhd::mhd::HD_State_Conservative()][pamhd::mhd::Mass_Density()]
				)[0];
			}
		);

		write_gnuplot_cmd_current(
			"Vy",
			"\n" + velocity_cmd + " 2\"\n",
			[](const pamhd::mhd::Cell& cell_data)->double {
				return pamhd::mhd::get_velocity(
					cell_data[pamhd::mhd::HD_State_Conservative()][pamhd::mhd::Momentum_Density()],
					cell_data[pamhd::mhd::HD_State_Conservative()][pamhd::mhd::Mass_Density()]
				)[1];
			}
		);

		write_gnuplot_cmd_current(
			"Vz",
			"\n" + velocity_cmd + " 3\"\n",
			[](const pamhd::mhd::Cell& cell_data)->double {
				return pamhd::mhd::get_velocity(
					cell_data[pamhd::mhd::HD_State_Conservative()][pamhd::mhd::Momentum_Density()],
					cell_data[pamhd::mhd::HD_State_Conservative()][pamhd::mhd::Mass_Density()]
				)[2];
			}
		);
	}

	// magnetic field
	if (magnetic_field_cmd != "") {
		write_gnuplot_cmd_current(
			"Bx",
			"\n" + magnetic_field_cmd + " 1\"\n",
			[](const pamhd::mhd::Cell& cell_data){
				return cell_data[pamhd::Magnetic_Field()][0];
			}
		);

		write_gnuplot_cmd_current(
			"By",
			"\n" + magnetic_field_cmd + " 2\"\n",
			[](const pamhd::mhd::Cell& cell_data){
				return cell_data[pamhd::Magnetic_Field()][1];
			}
		);

		write_gnuplot_cmd_current(
			"Bz",
			"\n" + magnetic_field_cmd + " 3\"\n",
			[](const pamhd::mhd::Cell& cell_data){
				return cell_data[pamhd::Magnetic_Field()][2];
			}
		);
	}

	// total magnetic field
	if (magnetic_field_cmd != "") {
		write_gnuplot_cmd_current(
			"Bx_tot",
			"\n" + magnetic_field_cmd + " 1\"\n",
			[](const pamhd::mhd::Cell& cell_data){
				return cell_data[pamhd::Magnetic_Field()][0] + cell_data[pamhd::Bg_Magnetic_Field_Pos_X()][0];
			}
		);

		write_gnuplot_cmd_current(
			"By_tot",
			"\n" + magnetic_field_cmd + " 2\"\n",
			[](const pamhd::mhd::Cell& cell_data){
				return cell_data[pamhd::Magnetic_Field()][1] + cell_data[pamhd::Bg_Magnetic_Field_Pos_X()][1];
			}
		);

		write_gnuplot_cmd_current(
			"Bz_tot",
			"\n" + magnetic_field_cmd + " 3\"\n",
			[](const pamhd::mhd::Cell& cell_data){
				return cell_data[pamhd::Magnetic_Field()][2] + cell_data[pamhd::Bg_Magnetic_Field_Pos_X()][2];
			}
		);
	}

	// background magnetic field
	if (magnetic_field_cmd != "") {
		write_gnuplot_cmd_current(
			"Bx_bg",
			"\n" + magnetic_field_cmd + " 1\"\n",
			[](const pamhd::mhd::Cell& cell_data){
				return cell_data[pamhd::Bg_Magnetic_Field_Pos_X()][0];
			}
		);

		write_gnuplot_cmd_current(
			"By_bg",
			"\n" + magnetic_field_cmd + " 2\"\n",
			[](const pamhd::mhd::Cell& cell_data){
				return cell_data[pamhd::Bg_Magnetic_Field_Pos_X()][1];
			}
		);

		write_gnuplot_cmd_current(
			"Bz_bg",
			"\n" + magnetic_field_cmd + " 3\"\n",
			[](const pamhd::mhd::Cell& cell_data){
				return cell_data[pamhd::Bg_Magnetic_Field_Pos_X()][2];
			}
		);
	}

	// electric field
	if (electric_field_cmd != "") {
		write_gnuplot_cmd_current(
			"Ex",
			"\n" + electric_field_cmd + " 1\"\n",
			[](const pamhd::mhd::Cell& cell_data){
				return cell_data[pamhd::Edge_Electric_Field()][0];
			}
		);

		write_gnuplot_cmd_current(
			"Ey",
			"\n" + electric_field_cmd + " 2\"\n",
			[](const pamhd::mhd::Cell& cell_data){
				return cell_data[pamhd::Edge_Electric_Field()][1];
			}
		);

		write_gnuplot_cmd_current(
			"Ez",
			"\n" + electric_field_cmd + " 3\"\n",
			[](const pamhd::mhd::Cell& cell_data){
				return cell_data[pamhd::Edge_Electric_Field()][2];
			}
		);
	}

	// current density
	if (current_density_cmd != "") {
		write_gnuplot_cmd_current(
			"Jx",
			"\n" + current_density_cmd + " 1\"\n",
			[](const pamhd::mhd::Cell& cell_data){
				return cell_data[pamhd::Electric_Current_Density()][0];
			}
		);

		write_gnuplot_cmd_current(
			"Jy",
			"\n" + current_density_cmd + " 2\"\n",
			[](const pamhd::mhd::Cell& cell_data){
				return cell_data[pamhd::Electric_Current_Density()][1];
			}
		);

		write_gnuplot_cmd_current(
			"Jz",
			"\n" + current_density_cmd + " 3\"\n",
			[](const pamhd::mhd::Cell& cell_data){
				return cell_data[pamhd::Electric_Current_Density()][2];
			}
		);
	}

	// MPI rank
	if (rank_cmd != "") {
		write_gnuplot_cmd_current(
			"rank",
			"\n" + rank_cmd + "\n",
			[](const pamhd::mhd::Cell& cell_data){
				return cell_data[pamhd::MPI_Rank()];
			}
		);
	}

	// resistivity
	if (resistivity_cmd != "") {
		write_gnuplot_cmd_current(
			"res",
			"\n" + resistivity_cmd + "\n",
			[](const pamhd::mhd::Cell& cell_data){
				return cell_data[pamhd::Resistivity()];
			}
		);
	}

	// boundary info
	if (bdy_cmd != "") {
		write_gnuplot_cmd_current(
			"bdy",
			"\n" + bdy_cmd + "\n",
			[](const pamhd::mhd::Cell& cell_data){
				return cell_data[pamhd::mhd::Solver_Info()];
			}
		);
	}


	/*
	Scale vector lengths
	*/

	double max_v = 0, max_B = 0, max_B0 = 0;
	for (size_t i = 0; i < cells.size(); i++) {
		const auto& hd_data
			= simulation_data.at(cells[i])[pamhd::mhd::HD_State_Conservative()];

		const auto rho = hd_data[pamhd::mhd::Mass_Density()];
		if (rho > 0) {
			const double v
				= pamhd::mhd::get_velocity(
					hd_data[pamhd::mhd::Momentum_Density()],
					hd_data[pamhd::mhd::Mass_Density()]
				).norm();

			if (max_v < v) {
				max_v = v;
			}
		}

		const double B = simulation_data.at(cells[i])[pamhd::Magnetic_Field()].norm();
		if (max_B < B) {
			max_B = B;
		}

		const double B0 = simulation_data.at(cells[i])[pamhd::Bg_Magnetic_Field_Pos_X()].norm();
		if (max_B0 < B0) {
			max_B0 = B0;
		}
	}
	if (max_v <= 0) {
		max_v = 1;
	}
	if (max_B <= 0) {
		max_B = 1;
	}
	if (max_B0 <= 0) {
		max_B0 = 1;
	}

	constexpr size_t vectors_per_screen = 1000;
	const double
		screen_length = geometry.get_end()[0] - geometry.get_start()[0],
		vector_length = screen_length / sqrt(vectors_per_screen),
		velocity_scaling = vector_length / max_v,
		B_scaling = vector_length / max_B,
		B0_scaling = 10 * vector_length / max_B0;

	// velocity
	gnuplot_file
		<< common_cmd
		<< "\nset output '"
		<< output_file_name_prefix + "_V." << output_file_name_suffix
		<< "'\nset ylabel 'Dimension "
		<< to_string(dimensions[1] + 1)
		<< "'\nset xlabel 'Dimension "
		<< to_string(dimensions[0] + 1)
		<< "'\nset title 'Velocity'\n"
		   "plot '-' u 1:2:3:4 w vectors t ''\n";

	for (size_t i = 0; i < cells.size(); i++) {
		if (
			cells.size() > vectors_per_screen
			and i % (cells.size() / vectors_per_screen) != 0
		) {
			continue;
		}

		const auto& hd_data
			= simulation_data.at(cells[i])[pamhd::mhd::HD_State_Conservative()];

		const auto& m = hd_data[pamhd::mhd::Momentum_Density()];
		const auto& rho = hd_data[pamhd::mhd::Mass_Density()];
		const auto cell_center = geometry.get_center(cells[i]);
		gnuplot_file
			<< cell_center[dimensions[0]] << " "
			<< cell_center[dimensions[1]] << " ";
		if (rho > 0) {
			gnuplot_file
				<< velocity_scaling * pamhd::mhd::get_velocity(m, rho)[dimensions[0]] << " "
				<< velocity_scaling * pamhd::mhd::get_velocity(m, rho)[dimensions[1]] << "\n";
		} else {
			gnuplot_file << "0 0\n";
		}
	}
	gnuplot_file << "\nend\n";

	// magnetic field
	gnuplot_file
		<< "\nset output '"
		<< output_file_name_prefix + "_B." << output_file_name_suffix
		<< "'\nset title 'Magnetic field'\n"
		   "plot '-' u 1:2:3:4 w vectors t ''\n";

	for (size_t i = 0; i < cells.size(); i++) {
		if (
			cells.size() > vectors_per_screen
			and i % (cells.size() / vectors_per_screen) != 0
		) {
			continue;
		}

		const auto& B = simulation_data.at(cells[i])[pamhd::Magnetic_Field()];
		const auto cell_center = geometry.get_center(cells[i]);
		gnuplot_file
			<< cell_center[dimensions[0]] << " "
			<< cell_center[dimensions[1]] << " "
			<< B_scaling * B[dimensions[0]] << " "
			<< B_scaling * B[dimensions[1]] << "\n";
	}
	gnuplot_file << "\nend\n";

	// background magnetic field on one face
	gnuplot_file
		<< "\nset output '"
		<< output_file_name_prefix + "_B0." << output_file_name_suffix
		<< "'\nset title 'Background magnetic field'\n"
		   "plot '-' u 1:2:3:4 w vectors t ''\n";

	for (size_t i = 0; i < cells.size(); i++) {
		if (
			cells.size() > vectors_per_screen
			and i % (cells.size() / vectors_per_screen) != 0
		) {
			continue;
		}

		const auto& B0 = simulation_data.at(cells[i])[pamhd::Bg_Magnetic_Field_Pos_X()];
		const auto cell_center = geometry.get_center(cells[i]);
		gnuplot_file
			<< cell_center[dimensions[0]] << " "
			<< cell_center[dimensions[1]] << " "
			<< B0_scaling * B0[dimensions[0]] << " "
			<< B0_scaling * B0[dimensions[1]] << "\n";
	}
	gnuplot_file << "\nend\n";


	gnuplot_file.close();

	return system(("gnuplot " + gnuplot_file_name).c_str());
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
	std::string
		common_plot_1d("set term png enhanced size 800, 600"),
		density_pressure_plot_1d(
			"set ylabel \"Mass density\" textcolor lt 1\n"
			"set y2label \"Pressure\" textcolor lt 3\n"
			"unset key\n"
			"set ytics nomirror\n"
			"set format x '%.2e'\n"
			"set format y '%.2e'\n"
			"set format y2 '%.2e'\n"
			"set y2tics auto"
		),
		velocity_plot_1d("set ylabel \"Velocity\""),
		magnetic_field_plot_1d("set ylabel \"Magnetic field\""),
		electric_field_plot_1d("set ylabel \"Electric field\""),
		current_density_plot_1d("set ylabel \"Current density\""),
		common_plot_2d(
			"set term png enhanced size 800, 600\n"
			"set pal gray\n"
			"set format cb \"%.2e\""
		),
		mass_density_plot_2d("set title \"Mass density\""),
		pressure_plot_2d("set title \"Pressure\""),
		velocity_plot_2d("set title \"Velocity"),
		magnetic_field_plot_2d("set title \"Magnetic field"),
		electric_field_plot_2d("set title \"Electric field"),
		current_density_plot_2d("set title \"Current density"),
		rank_plot_2d("set title \"MPI rank\""),
		resistivity_plot_2d("set title \"Resistivity\""),
		bdy_plot_2d("set title \"Type\"");

	boost::program_options::options_description
		options("Usage: program_name [options], where options are");
	options.add_options()
		("help", "Print this help message")
		("verbose", "Print run time information")
		("input-file",
			boost::program_options::value<std::vector<std::string>>(&input_files)
				->default_value(input_files),
			"Results to plot, the --input-file part can be omitted")
		("common-1d",
			boost::program_options::value<std::string>(&common_plot_1d)
				->default_value(common_plot_1d),
			"Gnuplot command(s) common to all variables in 1d")
		("density-pressure-1d",
			boost::program_options::value<std::string>(&density_pressure_plot_1d)
				->default_value(density_pressure_plot_1d),
			"Gnuplot command(s) for plotting mass density and pressure in 1d")
		("velocity-1d",
			boost::program_options::value<std::string>(&velocity_plot_1d)
				->default_value(velocity_plot_1d),
			"Gnuplot command(s) for plotting each component of velocity in 1d")
		("magnetic-field-1d",
			boost::program_options::value<std::string>(&magnetic_field_plot_1d)
				->default_value(magnetic_field_plot_1d),
			"Gnuplot command(s) for plotting each component of magnetic field in 1d")
		("current-density-1d",
			boost::program_options::value<std::string>(&current_density_plot_1d)
				->default_value(current_density_plot_1d),
			"Gnuplot command(s) for plotting each component of current density in 1d")
		("common-2d",
			boost::program_options::value<std::string>(&common_plot_2d)
				->default_value(common_plot_2d),
			"Gnuplot command(s) common to all variables in 2d")
		("mass-density-2d",
			boost::program_options::value<std::string>(&mass_density_plot_2d)
				->default_value(mass_density_plot_2d),
			"Gnuplot command(s) for plotting mass density in 2d")
		("pressure-2d",
			boost::program_options::value<std::string>(&pressure_plot_2d)
				->default_value(pressure_plot_2d),
			"Gnuplot command(s) for plotting pressure in 2d")
		("velocity-2d",
			boost::program_options::value<std::string>(&velocity_plot_2d)
				->default_value(velocity_plot_2d),
			"Gnuplot command(s) for plotting each component of velocity in 2d "
			"(component number and closing \" added automatically)")
		("magnetic-field-2d",
			boost::program_options::value<std::string>(&magnetic_field_plot_2d)
				->default_value(magnetic_field_plot_2d),
			"Gnuplot command(s) for plotting each component of magnetic field in 2d "
			"(component number and closing \" added automatically)")
		("current-density-2d",
			boost::program_options::value<std::string>(&current_density_plot_2d)
				->default_value(current_density_plot_2d),
			"Gnuplot command(s) for plotting each component of current density in 2d "
			"(component number and closing \" added automatically)")
		("rank-2d",
			boost::program_options::value<std::string>(&rank_plot_2d)
				->default_value(rank_plot_2d),
			"Gnuplot command(s) for plotting MPI rank in 2d")
		("resistivity-2d",
			boost::program_options::value<std::string>(&resistivity_plot_2d)
				->default_value(resistivity_plot_2d),
			"Gnuplot command(s) for plotting electric resistivity in 2d")
		("bdy-2d",
			boost::program_options::value<std::string>(&bdy_plot_2d)
				->default_value(bdy_plot_2d),
			"Gnuplot command(s) for plotting cell boundary type in 2d");

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
				"Each variable's plot command is given to gnuplot after "
				"the common plot command block. To include newlines in "
				"arguments use e.g. in bash ./mhd2gnuplot --opt "
				"$'set title \"T\"\\nset xrange...' which expands \\n etc. "
				"Variables with empty plot command are not plotted."
				<< endl;
		}
		MPI_Finalize();
		return EXIT_SUCCESS;
	}

	if (var_map.count("verbose") > 0) {
		verbose = true;
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
		unordered_map<uint64_t, pamhd::mhd::Cell> simulation_data;

		boost::optional<std::array<double, 4>> header = read_data(
			cell_id_mapping,
			topology,
			geometry,
			simulation_data,
			input_files[i],
			rank
		);
		if (not header) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< "Couldn't read simulation data from file " << input_files[i]
				<< std::endl;
			continue;
		}


		// cell coordinates increase with id (without AMR)
		std::vector<uint64_t> cells;
		cells.reserve(simulation_data.size());

		for (const auto& item: simulation_data) {
			cells.push_back(item.first);
		}
		if (cells.size() == 0) {
			continue;
		}

		std::sort(cells.begin(), cells.end());

		// get number of dimensions
		const auto grid_length = cell_id_mapping.length.get();
		const size_t dimensions
			= [&](){
				size_t ret_val = 0;
				for (size_t i = 0; i < grid_length.size(); i++) {
					if (grid_length[i] > 1) {
						ret_val++;
					}
				}
				return ret_val;
			}();

		switch(dimensions) {
		case 1:
			plot_1d(
				geometry,
				simulation_data,
				cells,
				input_files[i].substr(0, input_files[i].size() - 3),
				(*header)[0],
				(*header)[1],
				(*header)[3],
				common_plot_1d,
				density_pressure_plot_1d,
				velocity_plot_1d,
				magnetic_field_plot_1d,
				electric_field_plot_1d,
				current_density_plot_1d
			);
			break;
		case 2:
			plot_2d(
				geometry,
				simulation_data,
				cells,
				input_files[i].substr(0, input_files[i].size() - 3),
				"png",
				(*header)[0],
				(*header)[1],
				(*header)[3],
				common_plot_2d,
				mass_density_plot_2d,
				pressure_plot_2d,
				velocity_plot_2d,
				magnetic_field_plot_2d,
				electric_field_plot_2d,
				current_density_plot_2d,
				rank_plot_2d,
				resistivity_plot_2d,
				bdy_plot_2d
			);
			break;
		default:
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< "Unsupported number of dimensions in file " << input_files[i]
				<< std::endl;
			continue;
			break;
		}
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
