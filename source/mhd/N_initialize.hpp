/*
Two-population version of PAMHD MHD initialization.

Copyright 2015, 2016, 2017 Ilja Honkonen
Copyright 2019, 2023 Finnish Meteorological Institute
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

#ifndef PAMHD_MHD_N_INITIALIZE_HPP
#define PAMHD_MHD_N_INITIALIZE_HPP


#include "cmath"
#include "iostream"
#include "limits"

#include "dccrg.hpp"

#include "mhd/common.hpp"
#include "mhd/variables.hpp"


namespace pamhd {
namespace mhd {


/*!
Sets the initial state of MHD simulation and zeroes fluxes.

Getters should return a reference to data of corresponding variable
when given a simulation cell's data.
\param [Init_Cond] Initial condition class defined in MHD test program
\param [grid] Grid storing the cells to initialize
\param [cells] List of cells to initialize
\param [adiabatic_index] https://en.wikipedia.org/wiki/Heat_capacity_ratio
\param [vacuum_permeability] https://en.wikipedia.org/wiki/Vacuum_permeability
*/
template <
	class Geometries,
	class Init_Cond,
	class Background_Magnetic_Field,
	class Cell_Iterator,
	class Cell,
	class Geometry,
	class Mass_Density_Getter1,
	class Mass_Density_Getter2,
	class Momentum_Density_Getter1,
	class Momentum_Density_Getter2,
	class Total_Energy_Density_Getter1,
	class Total_Energy_Density_Getter2,
	class Magnetic_Field_Getter,
	class Background_Magnetic_Field_Getter,
	class Mass_Density_Flux_Getter1,
	class Mass_Density_Flux_Getter2,
	class Momentum_Density_Flux_Getter1,
	class Momentum_Density_Flux_Getter2,
	class Total_Energy_Density_Flux_Getter1,
	class Total_Energy_Density_Flux_Getter2,
	class Magnetic_Field_Flux_Getter
> void N_initialize(
	const Geometries& geometries,
	Init_Cond& initial_conditions,
	const Background_Magnetic_Field& bg_B,
	const Cell_Iterator& cells,
	dccrg::Dccrg<Cell, Geometry>& grid,
	const double time,
	const double adiabatic_index,
	const double vacuum_permeability,
	const double proton_mass,
	const bool verbose,
	const Mass_Density_Getter1 Mas1,
	const Mass_Density_Getter2 Mas2,
	const Momentum_Density_Getter1 Mom1,
	const Momentum_Density_Getter2 Mom2,
	const Total_Energy_Density_Getter1 Nrj1,
	const Total_Energy_Density_Getter2 Nrj2,
	const Magnetic_Field_Getter Mag,
	const Background_Magnetic_Field_Getter Bg_B,
	const Mass_Density_Flux_Getter1 Mas_f1,
	const Mass_Density_Flux_Getter2 Mas_f2,
	const Momentum_Density_Flux_Getter1 Mom_f1,
	const Momentum_Density_Flux_Getter2 Mom_f2,
	const Total_Energy_Density_Flux_Getter1 Nrj_f1,
	const Total_Energy_Density_Flux_Getter2 Nrj_f2,
	const Magnetic_Field_Flux_Getter Mag_f
) {
	if (verbose and grid.get_rank() == 0) {
		std::cout << "Setting default MHD state... ";
		std::cout.flush();
	}
	// set default state
	for (const auto& cell: cells) {
		// zero fluxes and background fields
		Mas_f1(*cell.data)         =
		Mas_f2(*cell.data)         =
		Nrj_f1(*cell.data)         =
		Nrj_f2(*cell.data)         =
		Mom_f1(*cell.data)[0]      =
		Mom_f1(*cell.data)[1]      =
		Mom_f1(*cell.data)[2]      =
		Mom_f2(*cell.data)[0]      =
		Mom_f2(*cell.data)[1]      =
		Mom_f2(*cell.data)[2]      =
		Mag_f(*cell.data)[0]       =
		Mag_f(*cell.data)[1]       =
		Mag_f(*cell.data)[2]       = 0;

		const auto c = grid.geometry.get_center(cell.id);
		const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
		const auto
			lat = asin(c[2] / r),
			lon = atan2(c[1], c[0]);

		const auto
			mass_density1
				= proton_mass
				* initial_conditions.get_default_data(
					pamhd::mhd::Number_Density(),
					time,
					c[0], c[1], c[2],
					r, lat, lon
				),
			mass_density2
				= proton_mass
				* initial_conditions.get_default_data(
					pamhd::mhd::Number_Density2(),
					time,
					c[0], c[1], c[2],
					r, lat, lon
				);

		const auto
			velocity1
				= initial_conditions.get_default_data(
					pamhd::mhd::Velocity(),
					time,
					c[0], c[1], c[2],
					r, lat, lon
				),
			velocity2
				= initial_conditions.get_default_data(
					pamhd::mhd::Velocity2(),
					time,
					c[0], c[1], c[2],
					r, lat, lon
				);
		const auto
			pressure1
				= initial_conditions.get_default_data(
					pamhd::mhd::Pressure(),
					time,
					c[0], c[1], c[2],
					r, lat, lon
				),
			pressure2
				= initial_conditions.get_default_data(
					pamhd::mhd::Pressure2(),
					time,
					c[0], c[1], c[2],
					r, lat, lon
				);
		const auto magnetic_field
			= initial_conditions.get_default_data(
				pamhd::Magnetic_Field(),
				time,
				c[0], c[1], c[2],
				r, lat, lon
			);

		Mas1(*cell.data) = mass_density1;
		Mas2(*cell.data) = mass_density2;
		Mom1(*cell.data) = mass_density1 * velocity1;
		Mom2(*cell.data) = mass_density2 * velocity2;
		Mag(*cell.data) = magnetic_field;
		if (mass_density1 > 0 and pressure1 > 0) {
			Nrj1(*cell.data) = get_total_energy_density(
				mass_density1,
				velocity1,
				pressure1,
				decltype(velocity1){0, 0, 0},
				adiabatic_index,
				vacuum_permeability
			);
		} else {
			Nrj1(*cell.data) = 0;
		}
		if (mass_density2 > 0 and pressure2 > 0) {
			Nrj2(*cell.data) = get_total_energy_density(
				mass_density2,
				velocity2,
				pressure2,
				decltype(velocity2){0, 0, 0},
				adiabatic_index,
				vacuum_permeability
			);
		} else {
			Nrj2(*cell.data) = 0;
		}

		const auto cell_end = grid.geometry.get_max(cell.id);
		Bg_B(*cell.data)(0,0) = bg_B.get_background_field(
			{cell_end[0], c[1], c[2]},
			vacuum_permeability
		);
		Bg_B(*cell.data)(1,0) = bg_B.get_background_field(
			{c[0], cell_end[1], c[2]},
			vacuum_permeability
		);
		Bg_B(*cell.data)(2,0) = bg_B.get_background_field(
			{c[0], c[1], cell_end[2]},
			vacuum_permeability
		);
	}


	if (verbose and grid.get_rank() == 0) {
		std::cout << "done\nSetting non-default initial state... ";
		std::cout.flush();
	}


	// mass densities
	for (
		size_t i = 0;
		i < initial_conditions.get_number_of_regions(pamhd::mhd::Number_Density());
		i++
	) {
		const auto& init_cond = initial_conditions.get_initial_condition(pamhd::mhd::Number_Density(), i);
		const auto& geometry_id = init_cond.get_geometry_id();
		const auto& cells = geometries.get_cells(geometry_id);
		for (const auto& cell: cells) {
			const auto c = grid.geometry.get_center(cell);
			const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
			const auto
				lat = asin(c[2] / r),
				lon = atan2(c[1], c[0]);

			const auto mass_density
				= proton_mass
				* initial_conditions.get_data(
					pamhd::mhd::Number_Density(),
					geometry_id,
					time,
					c[0], c[1], c[2],
					r, lat, lon
				);

			auto* const cell_data = grid[cell];
			if (cell_data == nullptr) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
				abort();
			}

			Mas1(*cell_data) = mass_density;
		}
	}
	for (
		size_t i = 0;
		i < initial_conditions.get_number_of_regions(pamhd::mhd::Number_Density2());
		i++
	) {
		const auto& init_cond = initial_conditions.get_initial_condition(pamhd::mhd::Number_Density2(), i);
		const auto& geometry_id = init_cond.get_geometry_id();
		const auto& cells = geometries.get_cells(geometry_id);
		for (const auto& cell: cells) {
			const auto c = grid.geometry.get_center(cell);
			const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
			const auto
				lat = asin(c[2] / r),
				lon = atan2(c[1], c[0]);

			const auto mass_density
				= proton_mass
				* initial_conditions.get_data(
					pamhd::mhd::Number_Density2(),
					geometry_id,
					time,
					c[0], c[1], c[2],
					r, lat, lon
				);

			auto* const cell_data = grid[cell];
			if (cell_data == nullptr) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
				abort();
			}

			Mas2(*cell_data) = mass_density;
		}
	}

	// velocities
	for (
		size_t i = 0;
		i < initial_conditions.get_number_of_regions(pamhd::mhd::Velocity());
		i++
	) {
		const auto& init_cond = initial_conditions.get_initial_condition(pamhd::mhd::Velocity(), i);
		const auto& geometry_id = init_cond.get_geometry_id();
		const auto& cells = geometries.get_cells(geometry_id);
		for (const auto& cell: cells) {
			const auto c = grid.geometry.get_center(cell);
			const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
			const auto
				lat = asin(c[2] / r),
				lon = atan2(c[1], c[0]);

			const auto velocity = initial_conditions.get_data(
				pamhd::mhd::Velocity(),
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

			Mom1(*cell_data) = Mas1(*cell_data) * velocity;
		}
	}
	for (
		size_t i = 0;
		i < initial_conditions.get_number_of_regions(pamhd::mhd::Velocity2());
		i++
	) {
		const auto& init_cond = initial_conditions.get_initial_condition(pamhd::mhd::Velocity2(), i);
		const auto& geometry_id = init_cond.get_geometry_id();
		const auto& cells = geometries.get_cells(geometry_id);
		for (const auto& cell: cells) {
			const auto c = grid.geometry.get_center(cell);
			const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
			const auto
				lat = asin(c[2] / r),
				lon = atan2(c[1], c[0]);

			const auto velocity = initial_conditions.get_data(
				pamhd::mhd::Velocity2(),
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

			Mom2(*cell_data) = Mas2(*cell_data) * velocity;
		}
	}

	// pressures, assuming no magnetic field
	for (
		size_t i = 0;
		i < initial_conditions.get_number_of_regions(pamhd::mhd::Pressure());
		i++
	) {
		std::cout << std::endl;
		const auto& init_cond = initial_conditions.get_initial_condition(pamhd::mhd::Pressure(), i);
		const auto& geometry_id = init_cond.get_geometry_id();
		std::cout << geometry_id << std::endl;
		const auto& cells = geometries.get_cells(geometry_id);
		std::cout << cells.size() << std::endl;
		for (const auto& cell: cells) {
			const auto c = grid.geometry.get_center(cell);
			const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
			const auto
				lat = asin(c[2] / r),
				lon = atan2(c[1], c[0]);

			const auto pressure = initial_conditions.get_data(
				pamhd::mhd::Pressure(),
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

			if (Mas1(*cell_data) > 0 and pressure > 0) {
				Nrj1(*cell_data) = get_total_energy_density(
					Mas1(*cell_data),
					get_velocity(Mom1(*cell_data), Mas1(*cell_data)),
					pressure,
					Magnetic_Field::data_type{0, 0, 0},
					adiabatic_index,
					vacuum_permeability
				);
			} else {
				Nrj1(*cell_data) = 0;
			}
		}
	}
	for (
		size_t i = 0;
		i < initial_conditions.get_number_of_regions(pamhd::mhd::Pressure2());
		i++
	) {
		std::cout << std::endl;
		const auto& init_cond = initial_conditions.get_initial_condition(pamhd::mhd::Pressure2(), i);
		const auto& geometry_id = init_cond.get_geometry_id();
		std::cout << geometry_id << std::endl;
		const auto& cells = geometries.get_cells(geometry_id);
		std::cout << cells.size() << std::endl;
		for (const auto& cell: cells) {
			const auto c = grid.geometry.get_center(cell);
			const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
			const auto
				lat = asin(c[2] / r),
				lon = atan2(c[1], c[0]);

			const auto pressure = initial_conditions.get_data(
				pamhd::mhd::Pressure2(),
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

			if (Mas2(*cell_data) > 0 and pressure > 0) {
				Nrj2(*cell_data) = get_total_energy_density(
					Mas2(*cell_data),
					get_velocity(Mom2(*cell_data), Mas2(*cell_data)),
					pressure,
					Magnetic_Field::data_type{0, 0, 0},
					adiabatic_index,
					vacuum_permeability
				);
			} else {
				Nrj2(*cell_data) = 0;
			}
		}
	}

	// magnetic field
	for (
		size_t i = 0;
		i < initial_conditions.get_number_of_regions(pamhd::Magnetic_Field());
		i++
	) {
		const auto& init_cond = initial_conditions.get_initial_condition(pamhd::Magnetic_Field(), i);
		const auto& geometry_id = init_cond.get_geometry_id();
		const auto& cells = geometries.get_cells(geometry_id);
		for (const auto& cell: cells) {
			const auto c = grid.geometry.get_center(cell);
			const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
			const auto
				lat = asin(c[2] / r),
				lon = atan2(c[1], c[0]);

			const auto magnetic_field = initial_conditions.get_data(
				pamhd::Magnetic_Field(),
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

			Mag(*cell_data) = magnetic_field;
		}
	}

	// add magnetic field contribution to total energy densities
	for (const auto& cell: cells) {
		const auto
			total_mass = Mas1(*cell.data) + Mas2(*cell.data),
			mass_frac1 = Mas1(*cell.data) / total_mass,
			mass_frac2 = Mas2(*cell.data) / total_mass;
		Nrj1(*cell.data) += mass_frac1 * 0.5 * Mag(*cell.data).squaredNorm() / vacuum_permeability;
		Nrj2(*cell.data) += mass_frac2 * 0.5 * Mag(*cell.data).squaredNorm() / vacuum_permeability;
	}
}

}} // namespaces

#endif // ifndef PAMHD_MHD_N_INITIALIZE_HPP
