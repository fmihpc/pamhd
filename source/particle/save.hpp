/*
Saves particle solution of PAMHD.

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

#ifndef PAMHD_PARTICLE_SAVE_HPP
#define PAMHD_PARTICLE_SAVE_HPP


#include "iomanip"

#include "dccrg.hpp"
#include "mpi.h"

#include "particle/variables.hpp"


namespace pamhd {
namespace particle {


/*!
Saves the particle solution into a file with name derived from simulation time.

file_name_prefix is added to the beginning of the file name.

The transfer of all first level variables must be switched
off before this function is called. After save returns the
transfer of all first level variables is switched off.

Cell must be compatible with gensimcell (github.com/nasailja/gensimcell)

Return true on success, false otherwise.
*/
template <
	class Electric_Field_T,
	class Magnetic_Field_T,
	class Electric_Current_Density_T,
	class Nr_Particles_T,
	class Particles_T,
	class Grid
> bool save(
	const std::string& file_name_prefix,
	Grid& grid,
	const uint64_t file_version,
	const uint64_t simulation_step,
	const double simulation_time,
	const double adiabatic_index,
	const double vacuum_permeability,
	const double particle_temp_nrj_ratio
) {
	const std::array<double, 4> simulation_parameters{{
		simulation_time,
		adiabatic_index,
		vacuum_permeability,
		particle_temp_nrj_ratio
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

	// update number of internal particles
	for (const auto& cell_id: grid.get_cells()) {
		auto* const cell_data = grid[cell_id];
		if (cell_data == nullptr) {
			std::cerr << __FILE__ << "(" << __LINE__ << ")" << std::endl;
			abort();
		}
		(*cell_data)[Nr_Particles_T()]
			= (*cell_data)[Particles_T()].size();
	}

	Grid::cell_data_type::set_transfer_all(
		true,
		Electric_Field_T(),
		Magnetic_Field_T(),
		Electric_Current_Density_T(),
		Nr_Particles_T(),
		Particles_T()
	);
	const bool ret_val = grid.save_grid_data(
		file_name_prefix + step_string.str() + ".dc",
		0,
		header
	);
	Grid::cell_data_type::set_transfer_all(
		false,
		Electric_Field_T(),
		Magnetic_Field_T(),
		Electric_Current_Density_T(),
		Nr_Particles_T(),
		Particles_T()
	);

	return ret_val;
}


}} // namespaces

#endif // ifndef PAMHD_PARTICLE_SAVE_HPP
