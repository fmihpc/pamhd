/*
Initializes the MHD solution of PAMHD.

Copyright 2014, 2015, 2016, 2017 Ilja Honkonen
Copyright 2019, 2024 Finnish Meteorological Institute
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

#ifndef PAMHD_MHD_INITIALIZE_HPP
#define PAMHD_MHD_INITIALIZE_HPP


#include "cmath"
#include "iostream"
#include "limits"
#include "tuple"

#include "dccrg.hpp"

#include "mhd/common.hpp"
#include "mhd/variables.hpp"


namespace pamhd {
namespace mhd {


// as initialize() below but for magnetic field only
template <
	class Boundary_Magnetic_Field,
	class Geometries,
	class Init_Cond,
	class Background_Magnetic_Field,
	class Grid,
	class Magnetic_Field_Getter,
	class Magnetic_Field_Flux_Getter,
	class Background_Magnetic_Field_Getter
> void initialize_magnetic_field(
	const Geometries& geometries,
	Init_Cond& initial_conditions,
	const Background_Magnetic_Field& bg_B,
	Grid& grid,
	const double time,
	const double vacuum_permeability,
	const Magnetic_Field_Getter Mag,
	const Magnetic_Field_Flux_Getter Mag_f,
	const Background_Magnetic_Field_Getter Bg_B
) {
	// set default magnetic field
	for (const auto& cell: grid.local_cells()) {
		Mag_f(*cell.data)[0] =
		Mag_f(*cell.data)[1] =
		Mag_f(*cell.data)[2] = 0;

		const auto c = grid.geometry.get_center(cell.id);
		const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
		const auto
			lat = asin(c[2] / r),
			lon = atan2(c[1], c[0]);

		const auto magnetic_field
			= initial_conditions.get_default_data(
				Boundary_Magnetic_Field(),
				time,
				c[0], c[1], c[2],
				r, lat, lon
			);

		Mag(*cell.data) = magnetic_field;

		const auto cell_end = grid.geometry.get_max(cell.id);
		Bg_B(*cell.data)(0, 1) = bg_B.get_background_field(
			{cell_end[0], c[1], c[2]},
			vacuum_permeability
		);
		Bg_B(*cell.data)(1, 1) = bg_B.get_background_field(
			{c[0], cell_end[1], c[2]},
			vacuum_permeability
		);
		Bg_B(*cell.data)(2, 1) = bg_B.get_background_field(
			{c[0], c[1], cell_end[2]},
			vacuum_permeability
		);
	}

	// set non-default magnetic field
	for (
		size_t i = 0;
		i < initial_conditions.get_number_of_regions(Boundary_Magnetic_Field());
		i++
	) {
		const auto& init_cond = initial_conditions.get_initial_condition(Boundary_Magnetic_Field(), i);
		const auto& geometry_id = init_cond.get_geometry_id();
		const auto& cells = geometries.get_cells(geometry_id);
		for (const auto& cell: cells) {
			const auto c = grid.geometry.get_center(cell);
			const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
			const auto
				lat = asin(c[2] / r),
				lon = atan2(c[1], c[0]);

			const auto magnetic_field = initial_conditions.get_data(
				Boundary_Magnetic_Field(),
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
}

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
	class Grid,
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Magnetic_Field_Getter,
	class Mass_Density_Flux_Getter,
	class Momentum_Density_Flux_Getter,
	class Total_Energy_Density_Flux_Getter
> void initialize_fluid(
	const Geometries& geometries,
	Init_Cond& initial_conditions,
	Grid& grid,
	const double time,
	const double adiabatic_index,
	const double vacuum_permeability,
	const double proton_mass,
	const bool verbose,
	const Mass_Density_Getter Mas,
	const Momentum_Density_Getter Mom,
	const Total_Energy_Density_Getter Nrj,
	const Magnetic_Field_Getter Mag,
	const Mass_Density_Flux_Getter Mas_f,
	const Momentum_Density_Flux_Getter Mom_f,
	const Total_Energy_Density_Flux_Getter Nrj_f
) {
	if (verbose and grid.get_rank() == 0) {
		std::cout << "Setting default MHD state... ";
		std::cout.flush();
	}

	// set default state
	for (const auto& cell: grid.local_cells()) {
		// zero fluxes
		Mas_f(*cell.data)         =
		Nrj_f(*cell.data)         =
		Mom_f(*cell.data)[0]      =
		Mom_f(*cell.data)[1]      =
		Mom_f(*cell.data)[2]      = 0;

		const auto c = grid.geometry.get_center(cell.id);
		const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
		const auto
			lat = asin(c[2] / r),
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

		Mas(*cell.data) = mass_density;
		Mom(*cell.data) = mass_density * velocity;
		if (mass_density > 0 and pressure > 0) {
			Nrj(*cell.data) = get_total_energy_density(
				mass_density,
				velocity,
				pressure,
				Mag(*cell.data),
				adiabatic_index,
				vacuum_permeability
			);
		} else {
			Nrj(*cell.data) = 0;
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
			const auto c = grid.geometry.get_center(cell);
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

			Mas(*cell_data) = mass_density;
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
			const auto c = grid.geometry.get_center(cell);
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

			Mom(*cell_data) = Mas(*cell_data) * velocity;
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
			const auto c = grid.geometry.get_center(cell);
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

			if (Mas(*cell_data) > 0 and pressure > 0) {
				Nrj(*cell_data) = get_total_energy_density(
					Mas(*cell_data),
					get_velocity(Mom(*cell_data), Mas(*cell_data)),
					pressure,
					Mag(*cell_data),
					adiabatic_index,
					vacuum_permeability
				);
			} else {
				Nrj(*cell_data) = 0;
			}
		}
	}

	if (verbose and grid.get_rank() == 0) {
		std::cout << "done" << std::endl;
	}
}


template <
	class Geometries,
	class Init_Cond,
	class Grid,
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Magnetic_Field_Getter,
	class Mass_Density_Flux_Getter,
	class Momentum_Density_Flux_Getter,
	class Total_Energy_Density_Flux_Getter
> void initialize_fluid_staggered(
	const Geometries& geometries,
	Init_Cond& initial_conditions,
	Grid& grid,
	const double time,
	const double adiabatic_index,
	const double vacuum_permeability,
	const double proton_mass,
	const bool verbose,
	const Mass_Density_Getter Mas,
	const Momentum_Density_Getter Mom,
	const Total_Energy_Density_Getter Nrj,
	const Magnetic_Field_Getter Mag,
	const Mass_Density_Flux_Getter Mas_f,
	const Momentum_Density_Flux_Getter Mom_f,
	const Total_Energy_Density_Flux_Getter Nrj_f
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

		const auto c = grid.geometry.get_center(cell.id);
		const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
		const auto
			lat = asin(c[2] / r),
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

		Mas(*cell.data) = mass_density;
		Mom(*cell.data) = mass_density * velocity;
		if (mass_density > 0 and pressure > 0) {
			Nrj(*cell.data) = get_total_energy_density(
				mass_density,
				velocity,
				pressure,
				Mag(*cell.data),
				adiabatic_index,
				vacuum_permeability
			);
		} else {
			Nrj(*cell.data) = 0;
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
			const auto c = grid.geometry.get_center(cell);
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

			Mas(*cell_data) = mass_density;
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
			const auto c = grid.geometry.get_center(cell);
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

			Mom(*cell_data) = Mas(*cell_data) * velocity;
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
			const auto c = grid.geometry.get_center(cell);
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

			if (Mas(*cell_data) > 0 and pressure > 0) {
				Nrj(*cell_data) = get_total_energy_density(
					Mas(*cell_data),
					get_velocity(Mom(*cell_data), Mas(*cell_data)),
					pressure,
					Mag(*cell_data),
					adiabatic_index,
					vacuum_permeability
				);
			} else {
				Nrj(*cell_data) = 0;
			}
		}
	}

	if (verbose and grid.get_rank() == 0) {
		std::cout << "done" << std::endl;
	}
}

}} // namespaces

#endif // ifndef PAMHD_MHD_INITIALIZE_HPP
