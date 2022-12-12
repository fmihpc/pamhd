/*
Saves the MHD solution of PAMHD.

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


#include "iomanip"


namespace pamhd {
namespace mhd {


/*!
Saves the MHD solution into a file with name derived from simulation time.

path_name_prefix is added to the beginning of the file name.

The transfer of all first level variables must be switched
off before this function is called. After save returns the
transfer of all first level variables is switched off.
Transfer of variables in MHD_T must be switched on.

Return true on success, false otherwise.
*/
template <
	class Grid,
	class... Variables
> bool save(
	const std::string& path_name_prefix,
	Grid& grid,
	const uint64_t file_version,
	const uint64_t simulation_step,
	const double simulation_time,
	const double adiabatic_index,
	const double proton_mass,
	const double vacuum_permeability,
	const Variables&... variables
) {
	const std::array<double, 4> simulation_parameters{{
		simulation_time,
		adiabatic_index,
		proton_mass,
		vacuum_permeability
	}};

	const std::array<int, 3> counts{{1, 1, 4}};
	const std::array<MPI_Aint, 3> displacements{{
		0,
		reinterpret_cast<char*>(const_cast<uint64_t*>(&simulation_step))
			- reinterpret_cast<char*>(const_cast<uint64_t*>(&file_version)),
		reinterpret_cast<char*>(const_cast<double*>(simulation_parameters.data()))
			- reinterpret_cast<char*>(const_cast<uint64_t*>(&file_version))
	}};
	std::array<MPI_Datatype, 3> datatypes{{MPI_UINT64_T, MPI_UINT64_T, MPI_DOUBLE}};

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

	Grid::cell_data_type::set_transfer_all(true, variables...);
	const bool ret_val = grid.save_grid_data(
		path_name_prefix + step_string.str() + ".dc",
		0, // start offset of grid data
		header
	);
	Grid::cell_data_type::set_transfer_all(false, variables...);

	return ret_val;
}


}} // namespaces

#endif // ifndef PAMHD_MHD_SAVE_HPP
