/*
Saves the MHD solution of PAMHD.

Copyright 2014, 2015, 2016, 2017 Ilja Honkonen
Copyright 2022, 2024 Finnish Meteorological Institute
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

#ifndef PAMHD_MHD_SAVE_HPP
#define PAMHD_MHD_SAVE_HPP


#include "cstdint"
#include "cstdio"
#include "iomanip"
#include "iterator"
#include "set"
#include "stdexcept"
#include "string"
#include "tuple"
#include "vector"

#include "mpi.h"

#include "grid/variables.hpp"
#include "mhd/variables.hpp"
#include "variables.hpp"


namespace pamhd {
namespace mhd {


/*!
Saves the MHD solution into files with names derived from simulation time.

path_name_prefix is added to the beginning of file names.

The transfer of all first level variables must be switched
off before this function is called. After save returns the
transfer of all first level variables is off.
Transfer of variables of 2nd etc levels must be switched on
in order to get written to file(s).

Return true on success, false otherwise.
*/
template <class Grid> bool save_staggered(
	const std::string& path_name_prefix,
	Grid& grid,
	const uint64_t file_version,
	const uint64_t simulation_step,
	const double simulation_time,
	const double adiabatic_index,
	const double proton_mass,
	const double vacuum_permeability,
	std::set<std::string> given_variables = std::set<std::string>()
) {
	using std::get;
	using std::string;
	using std::vector;

	std::set<string>
		variables{},
		allowed_variables{
			"mhd", "divfaceB", "bgB", "rank", "mhd info",
			"timestep", "substep", "substmin", "substmax",
			"ref lvls", "faceB", "fluxes"
		};
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
		proton_mass,
		vacuum_permeability,
		-1
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

	if (variables.count("mhd") > 0) {
		MPI_File_get_size(outfile, &outsize);
		variable_offsets.push_back(outsize);
		Grid::cell_data_type::set_transfer_all(true, pamhd::mhd::MHD_State_Conservative());
		const string varname = "mhd     ";
		get<0>(header) = (void*)varname.data();
		ret_val = ret_val and grid.save_grid_data(
			path_name_prefix + step_string.str() + ".dc",
			outsize, header, cells, false, false, false);
		Grid::cell_data_type::set_transfer_all(false, pamhd::mhd::MHD_State_Conservative());
	}

	if (variables.count("divfaceB") > 0) {
		MPI_File_get_size(outfile, &outsize);
		variable_offsets.push_back(outsize);
		Grid::cell_data_type::set_transfer_all(true, pamhd::Magnetic_Field_Divergence());
		const string varname = "divfaceB";
		get<0>(header) = (void*)varname.data();
		ret_val = ret_val and grid.save_grid_data(
			path_name_prefix + step_string.str() + ".dc",
			outsize, header, cells, false, false, false);
		Grid::cell_data_type::set_transfer_all(false, pamhd::Magnetic_Field_Divergence());
	}

	if (variables.count("bgB") > 0) {
		MPI_File_get_size(outfile, &outsize);
		variable_offsets.push_back(outsize);
		Grid::cell_data_type::set_transfer_all(true, pamhd::Bg_Magnetic_Field());
		const string varname = "bgB     ";
		get<0>(header) = (void*)varname.data();
		ret_val = ret_val and grid.save_grid_data(
			path_name_prefix + step_string.str() + ".dc",
			outsize, header, cells, false, false, false);
		Grid::cell_data_type::set_transfer_all(false, pamhd::Bg_Magnetic_Field());
	}

	if (variables.count("rank") > 0) {
		MPI_File_get_size(outfile, &outsize);
		variable_offsets.push_back(outsize);
		Grid::cell_data_type::set_transfer_all(true, pamhd::MPI_Rank());
		const string varname = "rank    ";
		get<0>(header) = (void*)varname.data();
		ret_val = ret_val and grid.save_grid_data(
			path_name_prefix + step_string.str() + ".dc",
			outsize, header, cells, false, false, false);
		Grid::cell_data_type::set_transfer_all(false, pamhd::MPI_Rank());
	}

	if (variables.count("mhd info") > 0) {
		MPI_File_get_size(outfile, &outsize);
		variable_offsets.push_back(outsize);
		Grid::cell_data_type::set_transfer_all(true, pamhd::Solver_Info());
		const string varname = "mhd info";
		get<0>(header) = (void*)varname.data();
		ret_val = ret_val and grid.save_grid_data(
			path_name_prefix + step_string.str() + ".dc",
			outsize, header, cells, false, false, false);
		Grid::cell_data_type::set_transfer_all(false, pamhd::Solver_Info());
	}

	if (variables.count("ref lvls") > 0) {
		MPI_File_get_size(outfile, &outsize);
		variable_offsets.push_back(outsize);
		Grid::cell_data_type::set_transfer_all(true,
			pamhd::grid::Target_Refinement_Level_Min(),
			pamhd::grid::Target_Refinement_Level_Max());
		const string varname = "ref lvls";
		get<0>(header) = (void*)varname.data();
		ret_val = ret_val and grid.save_grid_data(
			path_name_prefix + step_string.str() + ".dc",
			outsize, header, cells, false, false, false);
		Grid::cell_data_type::set_transfer_all(false,
			pamhd::grid::Target_Refinement_Level_Min(),
			pamhd::grid::Target_Refinement_Level_Max());
	}

	if (variables.count("faceB") > 0) {
		MPI_File_get_size(outfile, &outsize);
		variable_offsets.push_back(outsize);
		Grid::cell_data_type::set_transfer_all(true, pamhd::Face_Magnetic_Field());
		const string varname = "faceB   ";
		get<0>(header) = (void*)varname.data();
		ret_val = ret_val and grid.save_grid_data(
			path_name_prefix + step_string.str() + ".dc",
			outsize, header, cells, false, false, false);
		Grid::cell_data_type::set_transfer_all(false, pamhd::Face_Magnetic_Field());
	}

	if (variables.count("fluxes") > 0) {
		MPI_File_get_size(outfile, &outsize);
		variable_offsets.push_back(outsize);
		Grid::cell_data_type::set_transfer_all(true, pamhd::mhd::MHD_Flux());
		const string varname = "fluxes  ";
		get<0>(header) = (void*)varname.data();
		ret_val = ret_val and grid.save_grid_data(
			path_name_prefix + step_string.str() + ".dc",
			outsize, header, cells, false, false, false);
		Grid::cell_data_type::set_transfer_all(false, pamhd::mhd::MHD_Flux());
	}

	if (variables.count("substep") > 0) {
		MPI_File_get_size(outfile, &outsize);
		variable_offsets.push_back(outsize);
		Grid::cell_data_type::set_transfer_all(true, pamhd::mhd::Substepping_Period());
		const string varname = "substep ";
		get<0>(header) = (void*)varname.data();
		ret_val = ret_val and grid.save_grid_data(
			path_name_prefix + step_string.str() + ".dc",
			outsize, header, cells, false, false, false);
		Grid::cell_data_type::set_transfer_all(false, pamhd::mhd::Substepping_Period());
	}

	if (variables.count("substmin") > 0) {
		MPI_File_get_size(outfile, &outsize);
		variable_offsets.push_back(outsize);
		Grid::cell_data_type::set_transfer_all(true, pamhd::mhd::Substep_Min());
		const string varname = "substmin";
		get<0>(header) = (void*)varname.data();
		ret_val = ret_val and grid.save_grid_data(
			path_name_prefix + step_string.str() + ".dc",
			outsize, header, cells, false, false, false);
		Grid::cell_data_type::set_transfer_all(false, pamhd::mhd::Substep_Min());
	}

	if (variables.count("substmax") > 0) {
		MPI_File_get_size(outfile, &outsize);
		variable_offsets.push_back(outsize);
		Grid::cell_data_type::set_transfer_all(true, pamhd::mhd::Substep_Max());
		const string varname = "substmax";
		get<0>(header) = (void*)varname.data();
		ret_val = ret_val and grid.save_grid_data(
			path_name_prefix + step_string.str() + ".dc",
			outsize, header, cells, false, false, false);
		Grid::cell_data_type::set_transfer_all(false, pamhd::mhd::Substep_Max());
	}

	if (variables.count("timestep") > 0) {
		MPI_File_get_size(outfile, &outsize);
		variable_offsets.push_back(outsize);
		Grid::cell_data_type::set_transfer_all(true, pamhd::mhd::Timestep());
		const string varname = "timestep";
		get<0>(header) = (void*)varname.data();
		ret_val = ret_val and grid.save_grid_data(
			path_name_prefix + step_string.str() + ".dc",
			outsize, header, cells, false, false, false);
		Grid::cell_data_type::set_transfer_all(false, pamhd::mhd::Timestep());
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


template <class Grid> bool save(
	const std::string& path_name_prefix,
	Grid& grid,
	const uint64_t file_version,
	const uint64_t simulation_step,
	const double simulation_time,
	const double adiabatic_index,
	const double proton_mass,
	const double vacuum_permeability,
	std::set<std::string> variables_ = std::set<std::string>()
) {
	using std::get;
	using std::string;
	using std::vector;

	std::set<string>
		variables{},
		allowed_variables{"hd", "volume B", "Ecurrent", "res", "bgB", "rank", "mhd info"};
	if (variables_.size() == 0) {
		variables_ = allowed_variables;
	}
	std::set_intersection(
		variables_.cbegin(), variables_.cend(),
		allowed_variables.cbegin(), allowed_variables.cend(),
		std::inserter(variables, variables.begin())
	);

	const vector<double> simulation_parameters{
		simulation_time,
		adiabatic_index,
		proton_mass,
		vacuum_permeability
	};
	const int nr_sim_params = simulation_parameters.size();

	const uint8_t nr_var_offsets = variables.size();
	const vector<int> counts{1, 1, nr_sim_params, 1, nr_var_offsets};
	vector<uint64_t> variable_offsets(nr_var_offsets, 0);
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

	// make sure data if written in same order to all files
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

	if (variables.count("hd") > 0) {
		MPI_File_get_size(outfile, &outsize);
		variable_offsets.push_back(outsize);
		Grid::cell_data_type::set_transfer_all(true, pamhd::mhd::HD_State_Conservative());
		const string varname = "hd      ";
		get<0>(header) = (void*)varname.data();
		ret_val = ret_val and grid.save_grid_data(
			path_name_prefix + step_string.str() + ".dc",
			outsize, header, cells, false, false, false
		);
		Grid::cell_data_type::set_transfer_all(false, pamhd::mhd::HD_State_Conservative());
	}

	if (variables.count("bgB") > 0) {
		MPI_File_get_size(outfile, &outsize);
		variable_offsets.push_back(outsize);
		Grid::cell_data_type::set_transfer_all(true, pamhd::Bg_Magnetic_Field());
		const string varname = "bgB     ";
		get<0>(header) = (void*)varname.data();
		ret_val = ret_val and grid.save_grid_data(
			path_name_prefix + step_string.str() + ".dc",
			outsize, header, cells, false, false, false);
		Grid::cell_data_type::set_transfer_all(false, pamhd::Bg_Magnetic_Field());
	}

	if (variables.count("rank") > 0) {
		MPI_File_get_size(outfile, &outsize);
		variable_offsets.push_back(outsize);
		Grid::cell_data_type::set_transfer_all(true, pamhd::MPI_Rank());
		const string varname = "rank    ";
		get<0>(header) = (void*)varname.data();
		ret_val = ret_val and grid.save_grid_data(
			path_name_prefix + step_string.str() + ".dc",
			outsize, header, cells, false, false, false);
		Grid::cell_data_type::set_transfer_all(false, pamhd::MPI_Rank());
	}

	if (variables.count("mhd info") > 0) {
		MPI_File_get_size(outfile, &outsize);
		variable_offsets.push_back(outsize);
		Grid::cell_data_type::set_transfer_all(true, pamhd::Solver_Info());
		const string varname = "mhd info";
		get<0>(header) = (void*)varname.data();
		ret_val = ret_val and grid.save_grid_data(
			path_name_prefix + step_string.str() + ".dc",
			outsize, header, cells, false, false, false);
		Grid::cell_data_type::set_transfer_all(false, pamhd::Solver_Info());
	}

	if (variables.count("volume B") > 0) {
		MPI_File_get_size(outfile, &outsize);
		variable_offsets.push_back(outsize);
		Grid::cell_data_type::set_transfer_all(true, pamhd::Magnetic_Field());
		const string varname = "volume B";
		get<0>(header) = (void*)varname.data();
		ret_val = ret_val and grid.save_grid_data(
			path_name_prefix + step_string.str() + ".dc",
			outsize, header, cells, false, false, false);
		Grid::cell_data_type::set_transfer_all(false, pamhd::Magnetic_Field());
	}

	if (variables.count("res") > 0) {
		MPI_File_get_size(outfile, &outsize);
		variable_offsets.push_back(outsize);
		Grid::cell_data_type::set_transfer_all(true, pamhd::Resistivity());
		const string varname = "res     ";
		get<0>(header) = (void*)varname.data();
		ret_val = ret_val and grid.save_grid_data(
			path_name_prefix + step_string.str() + ".dc",
			outsize, header, cells, false, false, false);
		Grid::cell_data_type::set_transfer_all(false, pamhd::Resistivity());
	}

	if (variables.count("Ecurrent") > 0) {
		MPI_File_get_size(outfile, &outsize);
		variable_offsets.push_back(outsize);
		Grid::cell_data_type::set_transfer_all(true, pamhd::Electric_Current_Density());
		const string varname = "Ecurrent";
		get<0>(header) = (void*)varname.data();
		ret_val = ret_val and grid.save_grid_data(
			path_name_prefix + step_string.str() + ".dc",
			outsize, header, cells, false, false, false);
		Grid::cell_data_type::set_transfer_all(false, pamhd::Electric_Current_Density());
	}

	if (grid.get_rank() == 0) {
		if (nr_var_offsets != variable_offsets.size()) {
			throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
		}
		MPI_File_write_at(
			outfile,
			8 * (1+1+4) + 1,
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

#endif // ifndef PAMHD_MHD_SAVE_HPP
