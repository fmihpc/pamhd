/*
Saves particle solution of PAMHD.

Copyright 2015, 2016, 2017 Ilja Honkonen
Copyright 2024, 2025 Finnish Meteorological Institute
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

#ifndef PAMHD_PARTICLE_SAVE_HPP
#define PAMHD_PARTICLE_SAVE_HPP


#include "iomanip"
#include "set"
#include "vector"

#include "dccrg.hpp"
#include "mpi.h"

#include "particle/variables.hpp"
#include "variables.hpp"


namespace pamhd {
namespace particle {


/*!
Saves particle and related data to path derived from prefix.

Return true on success, false otherwise.
*/
template <class Grid> bool save(
	const std::string& path_name_prefix,
	Grid& grid,
	const uint64_t file_version,
	const uint64_t simulation_step,
	const double simulation_time,
	const double adiabatic_index,
	const double vacuum_permeability,
	const double particle_temp_nrj_ratio,
	std::set<std::string> given_variables = std::set<std::string>()
) {
	using std::get;
	using std::string;
	using std::vector;

	using Cell = Grid::cell_data_type;

	std::set<string>
		variables{},
		allowed_variables{"volE", "volJ", "nr ipart", "ipart"};
	if (given_variables.size() == 0) {
		given_variables = allowed_variables;
	}
	std::set_intersection(
		given_variables.cbegin(), given_variables.cend(),
		allowed_variables.cbegin(), allowed_variables.cend(),
		std::inserter(variables, variables.begin())
	);
	const uint8_t nr_var_offsets = variables.size();
	vector<uint64_t> variable_offsets(nr_var_offsets, 0);

	const vector<double> simulation_parameters{
		simulation_time,
		adiabatic_index,
		-1,
		vacuum_permeability,
		particle_temp_nrj_ratio
	};
	const int nr_sim_params = simulation_parameters.size();

	const vector<int> counts{1, 1, nr_sim_params, 1, nr_var_offsets};
	const vector<MPI_Aint> displacements{
		0,
		reinterpret_cast<char*>(const_cast<uint64_t*>(&simulation_step))
			- reinterpret_cast<char*>(const_cast<uint64_t*>(&file_version)),
		reinterpret_cast<char*>(const_cast<double*>(simulation_parameters.data()))
			- reinterpret_cast<char*>(const_cast<uint64_t*>(&file_version)),
		reinterpret_cast<char*>(const_cast<uint8_t*>(&nr_var_offsets))
			- reinterpret_cast<char*>(const_cast<uint64_t*>(&file_version)),
		reinterpret_cast<char*>(variable_offsets.data())
			- reinterpret_cast<char*>(const_cast<uint64_t*>(&file_version))
	};
	const vector<MPI_Datatype> datatypes{MPI_UINT64_T, MPI_UINT64_T, MPI_DOUBLE, MPI_UINT8_T, MPI_UINT64_T};

	MPI_Datatype header_datatype;
	if (
		MPI_Type_create_struct(
			counts.size(),
			counts.data(),
			displacements.data(),
			datatypes.data(),
			&header_datatype
		) != MPI_SUCCESS
	) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< " Couldn't create header datatype"
			<< std::endl;
		abort();
	}

	std::tuple<void*, int, MPI_Datatype> header{
		(void*) &file_version,
		1,
		header_datatype
	};

	std::ostringstream step_string;
	step_string << std::setw(9) << std::setfill('0') << simulation_step;

	// make sure data is written in same order to all files
	vector<uint64_t> cells = grid.get_cells();

	// assume transfer of all variables has been switched off
	bool ret_val = grid.save_grid_data(
		path_name_prefix + step_string.str() + ".dc",
		0, header, cells, true, true, false
	);
	variable_offsets.resize(0);

	// append variables to file one by one
	MPI_File outfile;
	MPI_File_open(
		MPI_COMM_WORLD,
		(path_name_prefix + step_string.str() + ".dc").data(),
		MPI_MODE_RDWR, MPI_INFO_NULL, &outfile
	);

	MPI_Offset outsize = 0;
	get<1>(header) = 8;
	get<2>(header) = MPI_BYTE;

	if (variables.count("volE") > 0) {
		MPI_File_get_size(outfile, &outsize);
		variable_offsets.push_back(outsize);
		Cell::set_transfer_all(true, pamhd::particle::Electric_Field());
		const string varname = "volE    ";
		get<0>(header) = (void*)varname.data();
		ret_val = ret_val and grid.save_grid_data(
			path_name_prefix + step_string.str() + ".dc",
			outsize, header, cells, false, false, false);
		Cell::set_transfer_all(false, pamhd::particle::Electric_Field());
	}

	if (variables.count("volJ") > 0) {
		MPI_File_get_size(outfile, &outsize);
		variable_offsets.push_back(outsize);
		Cell::set_transfer_all(true, pamhd::Electric_Current_Density());
		const string varname = "volJ    ";
		get<0>(header) = (void*)varname.data();
		ret_val = ret_val and grid.save_grid_data(
			path_name_prefix + step_string.str() + ".dc",
			outsize, header, cells, false, false, false);
		Cell::set_transfer_all(false, pamhd::Electric_Current_Density());
	}

	if (variables.count("nr ipart") > 0) {
		MPI_File_get_size(outfile, &outsize);
		variable_offsets.push_back(outsize);
		Cell::set_transfer_all(true, pamhd::particle::Nr_Particles_Internal());
		const string varname = "nr ipart";
		get<0>(header) = (void*)varname.data();
		ret_val = ret_val and grid.save_grid_data(
			path_name_prefix + step_string.str() + ".dc",
			outsize, header, cells, false, false, false);
		Cell::set_transfer_all(false, pamhd::particle::Nr_Particles_Internal());
	}

	if (variables.count("ipart") > 0) {
		MPI_File_get_size(outfile, &outsize);
		variable_offsets.push_back(outsize);
		Cell::set_transfer_all(true, pamhd::particle::Particles_Internal());
		const string varname = "ipart   ";
		get<0>(header) = (void*)varname.data();
		ret_val = ret_val and grid.save_grid_data(
			path_name_prefix + step_string.str() + ".dc",
			outsize, header, cells, false, false, false);
		Cell::set_transfer_all(false, pamhd::particle::Particles_Internal());
	}

	if (grid.get_rank() == 0) {
		if (nr_var_offsets != variable_offsets.size()) {
			throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
		}
		MPI_File_write_at(
			outfile,
			8 * (1+1+nr_sim_params) + 1,
			(void*)variable_offsets.data(),
			variable_offsets.size(),
			MPI_UINT64_T,
			MPI_STATUS_IGNORE
		);
	}
	MPI_File_close(&outfile);

	return ret_val;
}


/*!
Saves particle and related data to path derived from prefix.

Return true on success, false otherwise.
*/
template <class Grid> bool save_hyb(
	const std::string& path_name_prefix,
	Grid& grid,
	const uint64_t file_version,
	const uint64_t simulation_step,
	const double simulation_time,
	const double adiabatic_index,
	const double vacuum_permeability,
	const double particle_temp_nrj_ratio,
	std::set<std::string> given_variables = std::set<std::string>()
) {
	using std::get;
	using std::string;
	using std::vector;

	using Cell = Grid::cell_data_type;

	std::set<string>
		variables{},
		allowed_variables{"volJ", "nr ipart", "ipart"};
	if (given_variables.size() == 0) {
		given_variables = allowed_variables;
	}
	std::set_intersection(
		given_variables.cbegin(), given_variables.cend(),
		allowed_variables.cbegin(), allowed_variables.cend(),
		std::inserter(variables, variables.begin())
	);
	const uint8_t nr_var_offsets = variables.size();
	vector<uint64_t> variable_offsets(nr_var_offsets, 0);

	const vector<double> simulation_parameters{
		simulation_time,
		adiabatic_index,
		-1,
		vacuum_permeability,
		particle_temp_nrj_ratio
	};
	const int nr_sim_params = simulation_parameters.size();

	const vector<int> counts{1, 1, nr_sim_params, 1, nr_var_offsets};
	const vector<MPI_Aint> displacements{
		0,
		reinterpret_cast<char*>(const_cast<uint64_t*>(&simulation_step))
			- reinterpret_cast<char*>(const_cast<uint64_t*>(&file_version)),
		reinterpret_cast<char*>(const_cast<double*>(simulation_parameters.data()))
			- reinterpret_cast<char*>(const_cast<uint64_t*>(&file_version)),
		reinterpret_cast<char*>(const_cast<uint8_t*>(&nr_var_offsets))
			- reinterpret_cast<char*>(const_cast<uint64_t*>(&file_version)),
		reinterpret_cast<char*>(variable_offsets.data())
			- reinterpret_cast<char*>(const_cast<uint64_t*>(&file_version))
	};
	const vector<MPI_Datatype> datatypes{
		MPI_UINT64_T, MPI_UINT64_T,
		MPI_DOUBLE, MPI_UINT8_T, MPI_UINT64_T};

	MPI_Datatype header_datatype;
	if (
		MPI_Type_create_struct(
			counts.size(),
			counts.data(),
			displacements.data(),
			datatypes.data(),
			&header_datatype
		) != MPI_SUCCESS
	) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< " Couldn't create header datatype"
			<< std::endl;
		abort();
	}

	std::tuple<void*, int, MPI_Datatype> header{
		(void*) &file_version,
		1,
		header_datatype
	};

	std::ostringstream step_string;
	step_string << std::setw(9) << std::setfill('0') << simulation_step;

	// make sure data is written in same order to all files
	vector<uint64_t> cells = grid.get_cells();

	// assume transfer of all variables has been switched off
	bool ret_val = grid.save_grid_data(
		path_name_prefix + step_string.str() + ".dc",
		0, header, cells, true, true, false
	);
	variable_offsets.resize(0);

	// append variables to file one by one
	MPI_File outfile;
	MPI_File_open(
		MPI_COMM_WORLD,
		(path_name_prefix + step_string.str() + ".dc").data(),
		MPI_MODE_RDWR, MPI_INFO_NULL, &outfile
	);

	MPI_Offset outsize = 0;
	get<1>(header) = 8;
	get<2>(header) = MPI_BYTE;

	if (variables.count("volJ") > 0) {
		MPI_File_get_size(outfile, &outsize);
		variable_offsets.push_back(outsize);
		Cell::set_transfer_all(true, pamhd::Electric_Current_Density());
		const string varname = "volJ    ";
		get<0>(header) = (void*)varname.data();
		ret_val = ret_val and grid.save_grid_data(
			path_name_prefix + step_string.str() + ".dc",
			outsize, header, cells, false, false, false);
		Cell::set_transfer_all(false, pamhd::Electric_Current_Density());
	}

	if (variables.count("nr ipart") > 0) {
		MPI_File_get_size(outfile, &outsize);
		variable_offsets.push_back(outsize);
		Cell::set_transfer_all(true, pamhd::particle::Nr_Particles_Internal());
		const string varname = "nr ipart";
		get<0>(header) = (void*)varname.data();
		ret_val = ret_val and grid.save_grid_data(
			path_name_prefix + step_string.str() + ".dc",
			outsize, header, cells, false, false, false);
		Cell::set_transfer_all(false, pamhd::particle::Nr_Particles_Internal());
	}

	if (variables.count("ipart") > 0) {
		MPI_File_get_size(outfile, &outsize);
		variable_offsets.push_back(outsize);
		Cell::set_transfer_all(true, pamhd::particle::Particles_Internal());
		const string varname = "ipart   ";
		get<0>(header) = (void*)varname.data();
		ret_val = ret_val and grid.save_grid_data(
			path_name_prefix + step_string.str() + ".dc",
			outsize, header, cells, false, false, false);
		Cell::set_transfer_all(false, pamhd::particle::Particles_Internal());
	}

	if (grid.get_rank() == 0) {
		if (nr_var_offsets != variable_offsets.size()) {
			throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
		}
		MPI_File_write_at(
			outfile,
			8 * (1+1+nr_sim_params) + 1,
			(void*)variable_offsets.data(),
			variable_offsets.size(),
			MPI_UINT64_T,
			MPI_STATUS_IGNORE
		);
	}
	MPI_File_close(&outfile);

	return ret_val;
}


}} // namespaces

#endif // ifndef PAMHD_PARTICLE_SAVE_HPP
