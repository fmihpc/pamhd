/*
Initializes the MHD solution of PAMHD.

Copyright 2014, 2015, 2016, 2017 Ilja Honkonen
Copyright 2019, 2024, 2025 Finnish Meteorological Institute
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


Author(s): Ilja Honkonen
*/

#ifndef PAMHD_MHD_INITIALIZE_HPP
#define PAMHD_MHD_INITIALIZE_HPP


#include "cmath"
#include "iostream"
#include "limits"
#include "tuple"

#include "dccrg.hpp"

#include "grid/amr.hpp"
#include "common_functions.hpp"
#include "mhd/common.hpp"
#include "mhd/variables.hpp"


namespace pamhd {
namespace mhd {


/*! Sets initial state of MHD simulation and zeroes fluxes

Emulates 1d/2d system if grid length is initially 1 cell
in some dimension(s), even if it has been refined, in which
case some cells will overlap in affected dimension(s).
*/
template <
	class Geometries,
	class Init_Cond,
	class Grid,
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Volume_Magnetic_Field_Getter,
	class Mass_Density_Flux_Getter,
	class Momentum_Density_Flux_Getter,
	class Total_Energy_Density_Flux_Getter
> void initialize_fluid_staggered(
	const Geometries& geometries,
	Init_Cond& initial_conditions,
	Grid& grid,
	const double& time,
	const double& adiabatic_index,
	const double& vacuum_permeability,
	const double& proton_mass,
	const bool& verbose,
	const Mass_Density_Getter& Mas,
	const Momentum_Density_Getter& Mom,
	const Total_Energy_Density_Getter& Nrj,
	const Volume_Magnetic_Field_Getter& Vol_B,
	const Mass_Density_Flux_Getter& Mas_f,
	const Momentum_Density_Flux_Getter& Mom_f,
	const Total_Energy_Density_Flux_Getter& Nrj_f
) {
	using std::get;

	if (verbose and grid.get_rank() == 0) {
		std::cout << "Setting default MHD state... ";
		std::cout.flush();
	}

	// set default state
	for (const auto& cell: grid.local_cells()) {
		// zero fluxes
		Mas_f(*cell.data, -1) =
		Mas_f(*cell.data, +1) =
		Mas_f(*cell.data, -2) =
		Mas_f(*cell.data, +2) =
		Mas_f(*cell.data, -3) =
		Mas_f(*cell.data, +3) =
		Nrj_f(*cell.data, -1) =
		Nrj_f(*cell.data, +1) =
		Nrj_f(*cell.data, -2) =
		Nrj_f(*cell.data, +2) =
		Nrj_f(*cell.data, -3) =
		Nrj_f(*cell.data, +3) = 0;
		Mom_f(*cell.data, -1) =
		Mom_f(*cell.data, +1) =
		Mom_f(*cell.data, -2) =
		Mom_f(*cell.data, +2) =
		Mom_f(*cell.data, -3) =
		Mom_f(*cell.data, +3) = {0, 0, 0};

		const auto [c, s, e] = pamhd::grid::get_cell_geom_emulated(grid, cell.id);
		const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
		const auto
			lat = [&](){
				if (r == 0) return 0.0;
				else return asin(c[2] / r);
			}(),
			lon = atan2(c[1], c[0]);

		const auto mass_density
			= proton_mass
			* initial_conditions.get_default_data(
				Number_Density(),
				time,
				c[0], c[1], c[2],
				r, lat, lon
			);
		const auto velocity
			= initial_conditions.get_default_data(
				Velocity(),
				time,
				c[0], c[1], c[2],
				r, lat, lon
			);
		const auto pressure
			= initial_conditions.get_default_data(
				Pressure(),
				time,
				c[0], c[1], c[2],
				r, lat, lon
			);

		Mas.data(*cell.data) = mass_density;
		Mom.data(*cell.data) = pamhd::mul(mass_density, velocity);
		if (mass_density > 0 and pressure > 0) {
			Nrj.data(*cell.data) = get_total_energy_density(
				mass_density,
				velocity,
				pressure,
				Vol_B.data(*cell.data),
				adiabatic_index,
				vacuum_permeability
			);
		} else {
			Nrj.data(*cell.data) = 0;
		}
	}

	if (verbose and grid.get_rank() == 0) {
		std::cout << "done\nSetting non-default initial MHD state... ";
		std::cout.flush();
	}

	// mass density
	for (
		size_t i = 0;
		i < initial_conditions.get_number_of_regions(Number_Density());
		i++
	) {
		const auto& init_cond = initial_conditions.get_initial_condition(Number_Density(), i);
		const auto& geometry_id = init_cond.get_geometry_id();
		const auto& cells = geometries.get_cells(geometry_id);
		for (const auto& cell: cells) {
			const auto [c, s, e] = pamhd::grid::get_cell_geom_emulated(grid, cell);
			const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
			const auto
				lat = asin(c[2] / r),
				lon = atan2(c[1], c[0]);

			const auto mass_density
				= proton_mass
				* initial_conditions.get_data(
					Number_Density(),
					geometry_id,
					time,
					c[0], c[1], c[2],
					r, lat, lon
				);

			auto* const cell_data = grid[cell];
			if (cell_data == nullptr) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << std::endl;
				abort();
			}

			Mas.data(*cell_data) = mass_density;
		}
	}

	// velocity
	for (
		size_t i = 0;
		i < initial_conditions.get_number_of_regions(Velocity());
		i++
	) {
		const auto& init_cond = initial_conditions.get_initial_condition(Velocity(), i);
		const auto& geometry_id = init_cond.get_geometry_id();
		const auto& cells = geometries.get_cells(geometry_id);
		for (const auto& cell: cells) {
			const auto [c, s, e] = pamhd::grid::get_cell_geom_emulated(grid, cell);
			const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
			const auto
				lat = asin(c[2] / r),
				lon = atan2(c[1], c[0]);

			const auto velocity = initial_conditions.get_data(
				Velocity(),
				geometry_id,
				time,
				c[0], c[1], c[2],
				r, lat, lon
			);

			auto* const cell_data = grid[cell];
			if (cell_data == nullptr) {
				std::cerr <<  __FILE__ << "(" << __LINE__
					<< ") No data for cell: " << cell
					<< std::endl;
				abort();
			}

			Mom.data(*cell_data) = pamhd::mul(Mas.data(*cell_data), velocity);
		}
	}

	// pressure
	for (
		size_t i = 0;
		i < initial_conditions.get_number_of_regions(Pressure());
		i++
	) {
		const auto& init_cond = initial_conditions.get_initial_condition(Pressure(), i);
		const auto& geometry_id = init_cond.get_geometry_id();
		const auto& cells = geometries.get_cells(geometry_id);
		for (const auto& cell: cells) {
			const auto [c, s, e] = pamhd::grid::get_cell_geom_emulated(grid, cell);
			const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
			const auto
				lat = asin(c[2] / r),
				lon = atan2(c[1], c[0]);

			const auto pressure = initial_conditions.get_data(
				Pressure(),
				geometry_id,
				time,
				c[0], c[1], c[2],
				r, lat, lon
			);

			auto* const cell_data = grid[cell];
			if (cell_data == nullptr) {
				std::cerr <<  __FILE__ << "(" << __LINE__
					<< ") No data for cell: " << cell
					<< std::endl;
				abort();
			}

			if (Mas.data(*cell_data) > 0 and pressure > 0) {
				Nrj.data(*cell_data) = get_total_energy_density(
					Mas.data(*cell_data),
					get_velocity(Mom.data(*cell_data), Mas.data(*cell_data)),
					pressure,
					Vol_B.data(*cell_data),
					adiabatic_index,
					vacuum_permeability
				);
			} else {
				Nrj.data(*cell_data) = 0;
			}
		}
	}

	Mas.type().is_stale = true;
	Mom.type().is_stale = true;
	Nrj.type().is_stale = true;

	if (verbose and grid.get_rank() == 0) {
		std::cout << "done" << std::endl;
	}
}

}} // namespaces

#endif // ifndef PAMHD_MHD_INITIALIZE_HPP
