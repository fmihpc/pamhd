/*
Program for plotting particle output of PAMHD with gnuplot.

Copyright 2015, 2016, 2017 Ilja Honkonen
Copyright 2024 Finnish Meteorological Institute
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
#include "Eigen/Core" // must be included before gensimcell
#include "gensimcell.hpp"
#include "prettyprint.hpp"

#include "mhd/variables.hpp"
#include "particle/common.hpp"
#include "particle/save.hpp"
#include "particle/variables.hpp"

using namespace std;
using namespace pamhd::particle;

/*
Reads simulation data from given file.

Fills out grid info and simulation data withing given volume.

On success returns simulation_time, adiabatic index,
vacuum permeability and particle temperature to energy ratio
(Boltzmann constant).
Returns uninitialized value on error.
*/
boost::optional<std::array<double, 4>> read_data(
	const Eigen::Vector3d& volume_start,
	const Eigen::Vector3d& volume_end,
	dccrg::Mapping& cell_id_mapping,
	dccrg::Grid_Topology& topology,
	dccrg::Cartesian_Geometry& geometry,
	unordered_map<uint64_t, Cell_test_particle>& simulation_data,
	uint64_t& simulation_step,
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
	if (file_version != 1) {
		cerr << "Process " << mpi_rank
			<< " Unsupported file version: " << file_version
			<< endl;
		return boost::optional<std::array<double, 4>>();
	}

	uint64_t read_simulation_step = ~0ull;
	MPI_File_read_at(
		file,
		offset,
		&read_simulation_step,
		1,
		MPI_UINT64_T,
		MPI_STATUS_IGNORE
	);
	offset += sizeof(uint64_t);
	if (read_simulation_step > 999999999) {
		cerr << "Process " << mpi_rank
			<< " Suspicious simulation step: " << read_simulation_step
			<< endl;
		return boost::optional<std::array<double, 4>>();
	}
	simulation_step = read_simulation_step;

	// read physical constants
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


/*!
Returns data of given variable required for plot_2d.

Data is collected from particles within given ranges
in real and velocity spaces.

Data is returned for plot_2d in the same order as
its arguments where applicable.

Infinite components of r_start and r_end are modified
to coincide with boundaries of given geometry.

If any component of v_start or v_end is inf then goes
through particle data twice to get the min/max value.
*/
std::tuple<Eigen::MatrixXd, double, double, double, double> prepare_plot_data(
	const std::string& horizontal_variable,
	const std::string& vertical_variable,
	const std::string& plot_variable,
	const Eigen::Vector3d& r_start,
	const Eigen::Vector3d& r_end,
	const Eigen::Vector3d& v_start,
	const Eigen::Vector3d& v_end,
	const size_t& horizontal_resolution,
	const size_t& vertical_resolution,
	const unordered_map<uint64_t, Cell_test_particle>& simulation_data,
	const dccrg::Mapping& cell_id_mapping,
	const dccrg::Grid_Topology& topology,
	const dccrg::Cartesian_Geometry& geometry,
	const double particle_temp_nrj_ratio
) {
	static_assert(
		std::numeric_limits<double>::is_iec559,
		"IEC 559 / IEEE 754 required"
	);

	using std::isinf;

	Eigen::MatrixXd plot_data
		= Eigen::MatrixXd::Zero(vertical_resolution, horizontal_resolution);
	double horiz_max = -1, horiz_min = 1, vert_max = -1, vert_min = 1;

	if (
		horizontal_variable == vertical_variable
		or horizontal_resolution == 0
		or vertical_resolution == 0
		or r_start[0] >= r_end[0]
		or r_start[1] >= r_end[1]
		or r_start[2] >= r_end[2]
		or v_start[0] >= v_end[0]
		or v_start[1] >= v_end[1]
		or v_start[2] >= v_end[2]
	) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Invalid resolution or data range."
			<< std::endl;
		return std::make_tuple(plot_data, horiz_min, horiz_max, vert_min, vert_max);
	}

	if (
		not (
			plot_variable == "rx"
			or plot_variable == "ry"
			or plot_variable == "rz"
			or plot_variable == "vx"
			or plot_variable == "vy"
			or plot_variable == "vz"
			or plot_variable == "r"
			or plot_variable == "v"
			or plot_variable == "mass"
			or plot_variable == "count"
			or plot_variable == "temperature"
		)
	) {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Unsupported plot variable."
			<< std::endl;
		return std::make_tuple(plot_data, horiz_min, horiz_max, vert_min, vert_max);
	}

	// cap r_* at grid boundaries unless capped by user
	const auto
		grid_start = geometry.get_start(),
		grid_end = geometry.get_end();

	auto
		real_r_start = r_start,
		real_r_end = r_end;
	if (isinf(real_r_start[0])) {
		real_r_start[0] = grid_start[0];
	}
	if (isinf(real_r_start[1])) {
		real_r_start[1] = grid_start[1];
	}
	if (isinf(real_r_start[2])) {
		real_r_start[2] = grid_start[2];
	}
	if (isinf(real_r_end[0])) {
		real_r_end[0] = grid_end[0];
	}
	if (isinf(real_r_end[1])) {
		real_r_end[1] = grid_end[1];
	}
	if (isinf(real_r_end[2])) {
		real_r_end[2] = grid_end[2];
	}

	// cap v_* at data boundaries unless capped by user
	bool need_v_range = false;
	if (
		(horizontal_variable == "v" or vertical_variable == "v")
		or (
			(horizontal_variable == "vx" or vertical_variable == "vx")
			and (isinf(v_start[0]) or isinf(v_end[0]))
		) or (
			(horizontal_variable == "vy" or vertical_variable == "vy")
			and (isinf(v_start[1]) or isinf(v_end[1]))
		) or (
			(horizontal_variable == "vz" or vertical_variable == "vz")
			and (isinf(v_start[2]) or isinf(v_end[2]))
		)
	) {
		need_v_range = true;
	}

	auto
		real_v_start = v_start,
		real_v_end = v_end;
	if (need_v_range) {
		if (isinf(real_v_start[0])) {
			real_v_start[0] = std::numeric_limits<double>::max();
		}
		if (isinf(real_v_start[1])) {
			real_v_start[1] = std::numeric_limits<double>::max();
		}
		if (isinf(real_v_start[2])) {
			real_v_start[2] = std::numeric_limits<double>::max();
		}
		if (isinf(real_v_end[0])) {
			real_v_end[0] = std::numeric_limits<double>::lowest();
		}
		if (isinf(real_v_end[1])) {
			real_v_end[1] = std::numeric_limits<double>::lowest();
		}
		if (isinf(real_v_end[2])) {
			real_v_end[2] = std::numeric_limits<double>::lowest();
		}

		for (const auto& item: simulation_data) {
			const auto cell_id = item.first;
			const auto
				cell_start = geometry.get_min(cell_id),
				cell_end = geometry.get_max(cell_id);

			if (
				cell_end[0] < r_start[0]
				or cell_end[1] < r_start[1]
				or cell_end[2] < r_start[2]
				or cell_start[0] > r_end[0]
				or cell_start[1] > r_end[1]
				or cell_start[2] > r_end[2]
			) {
				continue;
			}

			const auto& cell_data = item.second;
			for (const auto& particle: cell_data[Particles_Internal()]) {
				const auto& vel = particle[Velocity()];
				// get range if not set by user
				if (isinf(v_start[0])) {
					real_v_start[0] = min(real_v_start[0], vel[0]);
				}
				if (isinf(v_start[1])) {
					real_v_start[1] = min(real_v_start[1], vel[1]);
				}
				if (isinf(v_start[2])) {
					real_v_start[2] = min(real_v_start[2], vel[2]);
				}
				if (isinf(v_end[0])) {
					real_v_end[0] = max(real_v_end[0], vel[0]);
				}
				if (isinf(v_end[1])) {
					real_v_end[1] = max(real_v_end[1], vel[1]);
				}
				if (isinf(v_end[2])) {
					real_v_end[2] = max(real_v_end[2], vel[2]);
				}
			}
		}
	}

	// get actual plot range
	if (horizontal_variable == "rx") {
		horiz_min = real_r_start[0];
		horiz_max = real_r_end[0];
	} else if (horizontal_variable == "ry") {
		horiz_min = real_r_start[1];
		horiz_max = real_r_end[1];
	} else if (horizontal_variable == "rz") {
		horiz_min = real_r_start[2];
		horiz_max = real_r_end[2];
	} else if (horizontal_variable == "vx") {
		horiz_min = real_v_start[0];
		horiz_max = real_v_end[0];
	} else if (horizontal_variable == "vy") {
		horiz_min = real_v_start[1];
		horiz_max = real_v_end[1];
	} else if (horizontal_variable == "vz") {
		horiz_min = real_v_start[2];
		horiz_max = real_v_end[2];
	} else if (horizontal_variable == "v") {
		horiz_min = 0;
		horiz_max = max(real_v_start.norm(), real_v_end.norm());
	} else if (horizontal_variable == "r") {
		horiz_min = 0;
		horiz_max = max(real_r_start.norm(), real_r_end.norm());
	} else {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Unsupported horizontal variable for plotting."
			<< std::endl;
		abort();
	}

	if (vertical_variable == "rx") {
		vert_min = real_r_start[0];
		vert_max = real_r_end[0];
	} else if (vertical_variable == "ry") {
		vert_min = real_r_start[1];
		vert_max = real_r_end[1];
	} else if (vertical_variable == "rz") {
		vert_min = real_r_start[2];
		vert_max = real_r_end[2];
	} else if (vertical_variable == "vx") {
		vert_min = real_v_start[0];
		vert_max = real_v_end[0];
	} else if (vertical_variable == "vy") {
		vert_min = real_v_start[1];
		vert_max = real_v_end[1];
	} else if (vertical_variable == "vz") {
		vert_min = real_v_start[2];
		vert_max = real_v_end[2];
	} else if (vertical_variable == "v") {
		vert_min = 0;
		vert_max = max(real_v_start.norm(), real_v_end.norm());
	} else if (vertical_variable == "r") {
		vert_min = 0;
		vert_max = max(real_r_start.norm(), real_r_end.norm());
	} else {
		std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
			<< "Unsupported vertical variable for plotting."
			<< std::endl;
		abort();
	}

	size_t nr_of_particles = 0;

	const double
		horiz_cell_length = (horiz_max - horiz_min) / horizontal_resolution,
		vert_cell_length = (vert_max - vert_min) / vertical_resolution;

	std::vector<std::vector<std::vector<Particle_Internal>>> plot_data_particles(
		vertical_resolution,
		std::vector<std::vector<Particle_Internal>>(
			horizontal_resolution,
			std::vector<Particle_Internal>()
		)
	);

	for (const auto& item: simulation_data) {
		const auto cell_id = item.first;
		const auto
			cell_start = geometry.get_min(cell_id),
			cell_end = geometry.get_max(cell_id);

		if (
			cell_end[0] < real_r_start[0]
			or cell_end[1] < real_r_start[1]
			or cell_end[2] < real_r_start[2]
			or cell_start[0] > real_r_end[0]
			or cell_start[1] > real_r_end[1]
			or cell_start[2] > real_r_end[2]
		) {
			continue;
		}

		const auto& cell_data = item.second;

		for (const auto& particle: cell_data[Particles_Internal()]) {
			const auto
				&r = particle[Position()],
				&v = particle[Velocity()];

			if (
				   r[0] < real_r_start[0]
				or r[1] < real_r_start[1]
				or r[2] < real_r_start[2]
				or r[0] > real_r_end[0]
				or r[1] > real_r_end[1]
				or r[2] > real_r_end[2]
				or v[0] < real_v_start[0]
				or v[1] < real_v_start[1]
				or v[2] < real_v_start[2]
				or v[0] > real_v_end[0]
				or v[1] > real_v_end[1]
				or v[2] > real_v_end[2]
			) {
				continue;
			}

			double value = 0;

			if (plot_variable == "rx") {
				value = r[0];
			} else if (plot_variable == "ry") {
				value = r[1];
			} else if (plot_variable == "rz") {
				value = r[2];
			} else if (plot_variable == "vx") {
				value = v[0];
			} else if (plot_variable == "vy") {
				value = v[1];
			} else if (plot_variable == "vz") {
				value = v[2];
			} else if (plot_variable == "r") {
				value = r.norm();
			} else if (plot_variable == "v") {
				value = v.norm();
			} else if (plot_variable == "mass") {
				value = particle[Mass()];
			} else if (plot_variable == "count") {
				value = 1;
			}

			// find where to assign value in plot_data
			double horiz = horiz_min, vert = vert_min;
			if (horizontal_variable == "rx") {
				if (r[0] < horiz_min or r[0] > horiz_max) { continue; }
				horiz = r[0];
			} else if (horizontal_variable == "ry") {
				if (r[1] < horiz_min or r[1] > horiz_max) { continue; }
				horiz = r[1];
			} else if (horizontal_variable == "rz") {
				if (r[2] < horiz_min or r[2] > horiz_max) { continue; }
				horiz = r[2];
			} else if (horizontal_variable == "vx") {
				if (v[0] < horiz_min or v[0] > horiz_max) { continue; }
					horiz = v[0];
			} else if (horizontal_variable == "vy") {
				if (v[1] < horiz_min or v[1] > horiz_max) { continue; }
				horiz = v[1];
			} else if (horizontal_variable == "vz") {
				if (v[2] < horiz_min or v[2] > horiz_max) { continue; }
				horiz = v[2];
			} else if (horizontal_variable == "v") {
			if (v.norm() < horiz_min or v.norm() > horiz_max) { continue; }
				horiz = v.norm();
			} else if (horizontal_variable == "r") {
				if (r.norm() < horiz_min or r.norm() > horiz_max) { continue; }
				horiz = r.norm();
			}
			if (vertical_variable == "rx") {
				if (r[0] < vert_min or r[0] > vert_max) { continue; }
				vert = r[0];
			} else if (vertical_variable == "ry") {
				if (r[1] < vert_min or r[1] > vert_max) { continue; }
					vert = r[1];
			} else if (vertical_variable == "rz") {
				if (r[2] < vert_min or r[2] > vert_max) { continue; }
				vert = r[2];
			} else if (vertical_variable == "vx") {
				if (v[0] < vert_min or v[0] > vert_max) { continue; }
			vert = v[0];
				} else if (vertical_variable == "vy") {
			if (v[1] < vert_min or v[1] > vert_max) { continue; }
				vert = v[1];
			} else if (vertical_variable == "vz") {
				if (v[2] < vert_min or v[2] > vert_max) { continue; }
				vert = v[2];
			} else if (vertical_variable == "v") {
				if (v.norm() < vert_min or v.norm() > vert_max) { continue; }
				vert = v.norm();
			} else if (vertical_variable == "r") {
				if (r.norm() < vert_min or r.norm() > vert_max) { continue; }
				vert = r.norm();
			}

			size_t
				horiz_index = min(
					horizontal_resolution - 1,
					size_t(floor((horiz - horiz_min) / horiz_cell_length))
				),
				vert_index = min(
					vertical_resolution - 1,
					size_t(floor((vert - vert_min) / vert_cell_length))
				);

			if (plot_variable == "temperature") {
				plot_data_particles[vert_index][horiz_index].push_back(particle);
			}

			plot_data(vert_index, horiz_index) += value;
			nr_of_particles++;
		}
	}

	if (fabs(horiz_min - horiz_max) / fabs(horiz_min + horiz_max) < 1e-4) {
		horiz_max += 1e-4 * fabs(horiz_max);
	}
	if (fabs(vert_min - vert_max) / fabs(vert_min + vert_max) < 1e-4) {
		vert_max += 1e-4 * fabs(vert_max);
	}

	if (plot_variable == "temperature") {
		for (size_t v_i = 0; v_i < vertical_resolution; v_i++) {
		for (size_t h_i = 0; h_i < horizontal_resolution; h_i++) {
			plot_data(v_i, h_i)
				= get_temperature<Mass, Velocity, Species_Mass>(
					plot_data_particles[v_i][h_i],
					particle_temp_nrj_ratio
				);
		}}
	}

	return std::make_tuple(plot_data, horiz_min, horiz_max, vert_min, vert_max);
}


/*
Plots data of given cells as a 2d color map.
*/
int plot_2d(
	const std::string& output_file_name,
	const Eigen::MatrixXd& plot_data,
	const std::string& plot_command_prefix,
	const std::string& xlabel,
	const std::string& ylabel,
	const std::string& cblabel,
	const double xrange_start,
	const double xrange_end,
	const double yrange_start,
	const double yrange_end
) {
	const double
		cell_length_x = (xrange_end - xrange_start) / plot_data.cols(),
		cell_length_y = (yrange_end - yrange_start) / plot_data.rows();

	const std::string gnuplot_file_name(output_file_name + ".gnuplot");
	std::ofstream gnuplot_file(gnuplot_file_name);

	gnuplot_file
		<< plot_command_prefix
		<< "\nset output '" << output_file_name
		<< "'\nset xlabel '" << xlabel
		<< "'\nset ylabel '" << ylabel
		<< "'\nset cblabel '" << cblabel
		<< "'\nset xrange[" << xrange_start
		<< " : " << xrange_end
		<< "]\nset yrange[" << yrange_start
		<< " : " << yrange_end
		<< "]\nplot '-' using ($1 / "
		<< plot_data.cols() << " * " << xrange_end - xrange_start << " + "
		// marks center of first data point in this dimension
		<< xrange_start + cell_length_x / 2 << "):($2 / "
		<< plot_data.rows() << " * " << yrange_end - yrange_start << " + "
		<< yrange_start + cell_length_y / 2
		<< "):3 matrix with image title ''\n";

	for (size_t row = 0; row < size_t(plot_data.rows()); row++) {
		for (size_t col = 0; col < size_t(plot_data.cols()); col++) {
			gnuplot_file << plot_data(row, col) << " ";
		}
		gnuplot_file << "\n";
	}
	gnuplot_file << "\nend";

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
		plot_variable("count"),
		horizontal_variable("rx"),
		vertical_variable("ry"),
		plot_command_prefix(
			"set term png enhanced size 800, 600\n"
			"set pal gray\n"
			"set format cb \"%.2e\""
		);
	constexpr double inf = std::numeric_limits<double>::infinity();
	Eigen::Vector3d
		r_start{-inf, -inf, -inf},
		r_end{inf, inf, inf},
		v_start{-inf, -inf, -inf},
		v_end{inf, inf, inf};
	size_t horizontal_resolution = 25, vertical_resolution = 25;

	boost::program_options::options_description
		options("Usage: program_name [options], where options are");
	options.add_options()
		("help", "Print this help message")
		("verbose", "Print run time information")
		("input-file",
			boost::program_options::value<>(&input_files)
				->default_value(input_files),
			"Results to plot, the --input-file part can be omitted")
		("plot-variable",
			boost::program_options::value<>(&plot_variable)
				->default_value(plot_variable),
			"Plot variable arg as a function of horizontal and vertical variables "
			"(one of r, rx, ry, rz, v, vx, vy, vz, mass, count, temperature)")
		("horizontal-variable",
			boost::program_options::value<>(&horizontal_variable)
				->default_value(horizontal_variable),
			"Use arg as horizontal variable (one of r, rx, ry, rz, v, vx, vy, vz)")
		("vertical-variable",
			boost::program_options::value<>(&vertical_variable)
				->default_value(vertical_variable),
			"Use arg as vertical variable (same choices as in horizontal-variable)")
		("horizontal-resolution",
			boost::program_options::value<>(&horizontal_resolution)
				->default_value(horizontal_resolution),
			"Use arg number of cells in horizontal direction of plot")
		("vertical-resolution",
			boost::program_options::value<>(&vertical_resolution)
				->default_value(vertical_resolution),
			"Use arg number of cells in vertical direction of plot")
		("rx-start",
			boost::program_options::value<>(&r_start[0])
				->default_value(r_start[0]),
			"Process particles with x coordinate of position > arg")
		("ry-start",
			boost::program_options::value<>(&r_start[1])
				->default_value(r_start[1]),
			"Process particles with y coordinate of position > arg")
		("rz-start",
			boost::program_options::value<>(&r_start[2])
				->default_value(r_start[2]),
			"Process particles with z coordinate of position > arg")
		("rx-end",
			boost::program_options::value<>(&r_end[0])
				->default_value(r_end[0]),
			"Process particles with x coordinate of position < arg")
		("ry-end",
			boost::program_options::value<>(&r_end[1])
				->default_value(r_end[1]),
			"Process particles with y coordinate of position < arg")
		("rz-end",
			boost::program_options::value<>(&r_end[2])
				->default_value(r_end[2]),
			"Process particles with z coordinate of position < arg")
		("vx-start",
			boost::program_options::value<>(&v_start[0])
				->default_value(v_start[0]),
			"Process particles with x coordinate of velocity > arg")
		("vy-start",
			boost::program_options::value<>(&v_start[1])
				->default_value(v_start[1]),
			"Process particles with y coordinate of velocity > arg")
		("vz-start",
			boost::program_options::value<>(&v_start[2])
				->default_value(v_start[2]),
			"Process particles with z coordinate of velocity > arg")
		("vx-end",
			boost::program_options::value<>(&v_end[0])
				->default_value(v_end[0]),
			"Process particles with x coordinate of velocity < arg")
		("vy-end",
			boost::program_options::value<>(&v_end[1])
				->default_value(v_end[1]),
			"Process particles with y coordinate of velocity < arg")
		("vz-end",
			boost::program_options::value<>(&v_end[2])
				->default_value(v_end[2]),
			"Process particles with z coordinate of velocity < arg")
		("common-2d",
			boost::program_options::value<std::string>(&plot_command_prefix)
				->default_value(plot_command_prefix),
			"String to prepend to gnuplot file");

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
				"To include newlines in arguments use e.g. in bash "
				"./mhd2gnuplot --opt $'set title \"T\"\\nset xrange...' "
				"which expands \\n etc."
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
		unordered_map<uint64_t, Cell_test_particle> simulation_data;
		uint64_t simulation_step;

		boost::optional<std::array<double, 4>> metadata
			= read_data(
				r_start,
				r_end,
				cell_id_mapping,
				topology,
				geometry,
				simulation_data,
				simulation_step,
				input_files[i],
				rank
			);
		if (not metadata) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< "Couldn't read simulation data from file " << input_files[i]
				<< std::endl;
			continue;
		}

		Eigen::MatrixXd plot_data;
		double x_start, x_end, y_start, y_end;
		std::tie(
			plot_data, x_start, x_end, y_start, y_end
		) = prepare_plot_data(
			horizontal_variable,
			vertical_variable,
			plot_variable,
			r_start,
			r_end,
			v_start,
			v_end,
			horizontal_resolution,
			vertical_resolution,
			simulation_data,
			cell_id_mapping,
			topology,
			geometry,
			(*metadata)[3]
		);

		// error or nothing to plot
		if (x_start > x_end or y_start > y_end) {
			std::cerr <<  __FILE__ << "(" << __LINE__<< "): "
				<< "Couldn't prepare plot data from file " << input_files[i]
				<< ", either there was nothing to plot or an error occurred."
				<< std::endl;
			continue;
		}

		const std::string output_name
			= input_files[i] + "_"
			+ horizontal_variable + "_"
			+ vertical_variable + "_"
			+ plot_variable + ".png";
		plot_2d(
			output_name,
			plot_data,
			plot_command_prefix,
			horizontal_variable,
			vertical_variable,
			plot_variable,
			x_start,
			x_end,
			y_start,
			y_end
		);
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
