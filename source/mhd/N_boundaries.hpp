/*
Handles boundary cell classification logic of MHD part of PAMHD.

Copyright 2015, 2016, 2017 Ilja Honkonen
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

#ifndef PAMHD_MHD_N_BOUNDARIES_HPP
#define PAMHD_MHD_N_BOUNDARIES_HPP


#include "algorithm"
#include "cmath"
#include "limits"
#include "map"
#include "string"
#include "utility"
#include "vector"

#include "dccrg.hpp"

#include "mhd/common.hpp"
#include "mhd/variables.hpp"


namespace pamhd {
namespace mhd {


/*!
Prepares boundary information needed for MHD solver about each simulation cell.
*/
template<
	class Solver_Info,
	class Cell_Data,
	class Geometry,
	class Boundaries,
	class Boundary_Geometries,
	class Solver_Info_Getter
> void N_set_solver_info(
	dccrg::Dccrg<Cell_Data, Geometry>& grid,
	Boundaries& boundaries,
	const Boundary_Geometries& geometries,
	const Solver_Info_Getter& SInfo
) {
	boundaries.classify(grid, geometries, SInfo);

	for (const auto& cell: grid.cells) {
		SInfo(*cell.data) = 1;
	}

	// number densities
	constexpr pamhd::mhd::Number_Density N{};
	for (const auto& cell: boundaries.get_value_boundary_cells(N)) {
		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		SInfo(*cell_data) = 0;
	}
	for (const auto& item: boundaries.get_copy_boundary_cells(N)) {
		auto* const cell_data = grid[item[0]];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		SInfo(*cell_data) = 0;
	}
	const std::set<uint64_t> dont_solve_mass(
		boundaries.get_dont_solve_cells(N).cbegin(),
		boundaries.get_dont_solve_cells(N).cend()
	);

	constexpr pamhd::mhd::Number_Density2 N2{};
	for (const auto& cell: boundaries.get_value_boundary_cells(N2)) {
		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		SInfo(*cell_data) = 0;
	}
	for (const auto& item: boundaries.get_copy_boundary_cells(N2)) {
		auto* const cell_data = grid[item[0]];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		SInfo(*cell_data) = 0;
	}
	const std::set<uint64_t> dont_solve_mass2(
		boundaries.get_dont_solve_cells(N2).cbegin(),
		boundaries.get_dont_solve_cells(N2).cend()
	);

	// velocities
	constexpr pamhd::mhd::Velocity V{};
	for (const auto& cell: boundaries.get_value_boundary_cells(V)) {
		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		SInfo(*cell_data) = 0;
	}
	for (const auto& item: boundaries.get_copy_boundary_cells(V)) {
		auto* const cell_data = grid[item[0]];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		SInfo(*cell_data) = 0;
	}
	const std::set<uint64_t> dont_solve_velocity(
		boundaries.get_dont_solve_cells(V).cbegin(),
		boundaries.get_dont_solve_cells(V).cend()
	);

	constexpr pamhd::mhd::Velocity2 V2{};
	for (const auto& cell: boundaries.get_value_boundary_cells(V2)) {
		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		SInfo(*cell_data) = 0;
	}
	for (const auto& item: boundaries.get_copy_boundary_cells(V2)) {
		auto* const cell_data = grid[item[0]];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		SInfo(*cell_data) = 0;
	}
	const std::set<uint64_t> dont_solve_velocity2(
		boundaries.get_dont_solve_cells(V2).cbegin(),
		boundaries.get_dont_solve_cells(V2).cend()
	);

	// pressures
	constexpr pamhd::mhd::Pressure P{};
	for (const auto& cell: boundaries.get_value_boundary_cells(P)) {
		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		SInfo(*cell_data) = 0;
	}
	for (const auto& item: boundaries.get_copy_boundary_cells(P)) {
		auto* const cell_data = grid[item[0]];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		SInfo(*cell_data) = 0;
	}
	const std::set<uint64_t> dont_solve_pressure(
		boundaries.get_dont_solve_cells(P).cbegin(),
		boundaries.get_dont_solve_cells(P).cend()
	);

	constexpr pamhd::mhd::Pressure2 P2{};
	for (const auto& cell: boundaries.get_value_boundary_cells(P2)) {
		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		SInfo(*cell_data) = 0;
	}
	for (const auto& item: boundaries.get_copy_boundary_cells(P2)) {
		auto* const cell_data = grid[item[0]];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		SInfo(*cell_data) = 0;
	}
	const std::set<uint64_t> dont_solve_pressure2(
		boundaries.get_dont_solve_cells(P2).cbegin(),
		boundaries.get_dont_solve_cells(P2).cend()
	);

	// magnetic field
	constexpr pamhd::Magnetic_Field B{};
	for (const auto& cell: boundaries.get_value_boundary_cells(B)) {
		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		SInfo(*cell_data) = 0;
	}
	for (const auto& item: boundaries.get_copy_boundary_cells(B)) {
		auto* const cell_data = grid[item[0]];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		SInfo(*cell_data) = 0;
	}
	const std::set<uint64_t> dont_solve_mag(
		boundaries.get_dont_solve_cells(B).cbegin(),
		boundaries.get_dont_solve_cells(B).cend()
	);

	if (
		not std::equal(
			dont_solve_mass.cbegin(),
			dont_solve_mass.cend(),
			dont_solve_mass2.cbegin()
		)
	) {
		throw std::invalid_argument(
			std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
			+ "Mass density and mass density2 dont_solves aren't equal."
		);
	}
	if (
		not std::equal(
			dont_solve_mass.cbegin(),
			dont_solve_mass.cend(),
			dont_solve_velocity.cbegin()
		)
	) {
		throw std::invalid_argument(
			std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
			+ "Mass density and velocity dont_solves aren't equal."
		);
	}
	if (
		not std::equal(
			dont_solve_mass.cbegin(),
			dont_solve_mass.cend(),
			dont_solve_velocity2.cbegin()
		)
	) {
		throw std::invalid_argument(
			std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
			+ "Mass density and velocity2 dont_solves aren't equal."
		);
	}
	if (
		not std::equal(
			dont_solve_mass.cbegin(),
			dont_solve_mass.cend(),
			dont_solve_pressure.cbegin()
		)
	) {
		throw std::invalid_argument(
			std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
			+ "Mass density and pressure dont_solves aren't equal."
		);
	}
	if (
		not std::equal(
			dont_solve_mass.cbegin(),
			dont_solve_mass.cend(),
			dont_solve_pressure2.cbegin()
		)
	) {
		throw std::invalid_argument(
			std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
			+ "Mass density and pressure2 dont_solves aren't equal."
		);
	}
	if (
		not std::equal(
			dont_solve_mass.cbegin(),
			dont_solve_mass.cend(),
			dont_solve_mag.cbegin()
		)
	) {
		throw std::invalid_argument(
			std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
			+ "Mass density and magnetic field dont_solves aren't equal."
		);
	}

	// don't solve cells in which no variable is solved
	for (auto& cell: dont_solve_mass) {
		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		SInfo(*cell_data) = -1;
	}

	grid.update_copies_of_remote_neighbors();
}


/*!
Applies boundaries of all simulation variables.

Value boundaries are applied to all cells
within that bundary's geometry. Value
boundaries are applied in order given in json
data, in case several overlap in geometry so
last one remains in effect.

Copy boundaries are applied after all value
boundaries. In case of more than one normal
neighbor their average is copied, vector
variables are processed by component.
*/
template<
	class Cell_Data,
	class Grid_Geometry,
	class Boundaries,
	class Boundary_Geometries,
	class Mass_Getter,
	class Mass_Getter2,
	class Momentum_Getter,
	class Momentum_Getter2,
	class Energy_Getter,
	class Energy_Getter2,
	class Magnetic_Field_Getter
> void N_apply_boundaries(
	dccrg::Dccrg<Cell_Data, Grid_Geometry>& grid,
	Boundaries& boundaries,
	const Boundary_Geometries& bdy_geoms,
	const double simulation_time,
	const Mass_Getter& Mas,
	const Mass_Getter2& Mas2,
	const Momentum_Getter& Mom,
	const Momentum_Getter2& Mom2,
	const Energy_Getter& Nrj,
	const Energy_Getter2& Nrj2,
	const Magnetic_Field_Getter& Mag,
	const double proton_mass,
	const double adiabatic_index,
	const double vacuum_permeability
) {
	// number densities
	constexpr pamhd::mhd::Number_Density N{};
	for (
		size_t i = 0;
		i < boundaries.get_number_of_value_boundaries(N);
		i++
	) {
		auto& value_bdy = boundaries.get_value_boundary(N, i);
		const auto& geometry_id = value_bdy.get_geometry_id();
		const auto& cells = bdy_geoms.get_cells(geometry_id);
		for (const auto& cell: cells) {
			const auto c = grid.geometry.get_center(cell);
			const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
			const auto
				lat = asin(c[2] / r),
				lon = atan2(c[1], c[0]);

			const auto mass_density
				= proton_mass
				* value_bdy.get_data(
					simulation_time,
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
	for (const auto& item: boundaries.get_copy_boundary_cells(N)) {
		if (item.size() < 2) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}

		pamhd::mhd::Number_Density::data_type source_value = 0.0;
		for (size_t i = 1; i < item.size(); i++) {
			auto* source_data = grid[item[i]];
			if (source_data == nullptr) {
				std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}

			source_value += Mas(*source_data);
		}
		source_value /= item.size() - 1;

		auto *target_data = grid[item[0]];
		if (target_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}

		Mas(*target_data) = source_value;
	}

	constexpr pamhd::mhd::Number_Density2 N2{};
	for (
		size_t i = 0;
		i < boundaries.get_number_of_value_boundaries(N2);
		i++
	) {
		auto& value_bdy = boundaries.get_value_boundary(N2, i);
		const auto& geometry_id = value_bdy.get_geometry_id();
		const auto& cells = bdy_geoms.get_cells(geometry_id);
		for (const auto& cell: cells) {
			const auto c = grid.geometry.get_center(cell);
			const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
			const auto
				lat = asin(c[2] / r),
				lon = atan2(c[1], c[0]);

			const auto mass_density
				= proton_mass
				* value_bdy.get_data(
					simulation_time,
					c[0], c[1], c[2],
					r, lat, lon
				);

			auto* const cell_data = grid[cell];
			if (cell_data == nullptr) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << std::endl;
				abort();
			}

			Mas2(*cell_data) = mass_density;
		}
	}
	for (const auto& item: boundaries.get_copy_boundary_cells(N2)) {
		if (item.size() < 2) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}

		pamhd::mhd::Number_Density2::data_type source_value = 0.0;
		for (size_t i = 1; i < item.size(); i++) {
			auto* source_data = grid[item[i]];
			if (source_data == nullptr) {
				std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}

			source_value += Mas2(*source_data);
		}
		source_value /= item.size() - 1;

		auto *target_data = grid[item[0]];
		if (target_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}

		Mas2(*target_data) = source_value;
	}

	// velocities
	constexpr pamhd::mhd::Velocity V{};
	for (
		size_t i = 0;
		i < boundaries.get_number_of_value_boundaries(V);
		i++
	) {
		auto& value_bdy = boundaries.get_value_boundary(V, i);
		const auto& geometry_id = value_bdy.get_geometry_id();
		const auto& cells = bdy_geoms.get_cells(geometry_id);
		for (const auto& cell: cells) {
			const auto c = grid.geometry.get_center(cell);
			const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
			const auto
				lat = asin(c[2] / r),
				lon = atan2(c[1], c[0]);

			const auto velocity = value_bdy.get_data(
				simulation_time,
				c[0], c[1], c[2],
				r, lat, lon
			);

			auto* const cell_data = grid[cell];
			if (cell_data == nullptr) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << std::endl;
				abort();
			}

			Mom(*cell_data) = Mas(*cell_data) * velocity;
		}
	}
	for (const auto& item: boundaries.get_copy_boundary_cells(V)) {
		if (item.size() < 2) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}

		pamhd::mhd::Velocity::data_type source_value{0, 0, 0};
		for (size_t i = 1; i < item.size(); i++) {
			auto* source_data = grid[item[i]];
			if (source_data == nullptr) {
				std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}

			source_value += pamhd::mhd::get_velocity(
				Mom(*source_data),
				Mas(*source_data)
			);
		}
		source_value /= item.size() - 1;

		auto *target_data = grid[item[0]];
		if (target_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}

		Mom(*target_data) = Mas(*target_data) * source_value;
	}

	constexpr pamhd::mhd::Velocity2 V2{};
	for (
		size_t i = 0;
		i < boundaries.get_number_of_value_boundaries(V2);
		i++
	) {
		auto& value_bdy = boundaries.get_value_boundary(V2, i);
		const auto& geometry_id = value_bdy.get_geometry_id();
		const auto& cells = bdy_geoms.get_cells(geometry_id);
		for (const auto& cell: cells) {
			const auto c = grid.geometry.get_center(cell);
			const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
			const auto
				lat = asin(c[2] / r),
				lon = atan2(c[1], c[0]);

			const auto velocity = value_bdy.get_data(
				simulation_time,
				c[0], c[1], c[2],
				r, lat, lon
			);

			auto* const cell_data = grid[cell];
			if (cell_data == nullptr) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << std::endl;
				abort();
			}

			Mom2(*cell_data) = Mas2(*cell_data) * velocity;
		}
	}
	for (const auto& item: boundaries.get_copy_boundary_cells(V2)) {
		if (item.size() < 2) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}

		pamhd::mhd::Velocity2::data_type source_value{0, 0, 0};
		for (size_t i = 1; i < item.size(); i++) {
			auto* source_data = grid[item[i]];
			if (source_data == nullptr) {
				std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}

			source_value += pamhd::mhd::get_velocity(
				Mom2(*source_data),
				Mas2(*source_data)
			);
		}
		source_value /= item.size() - 1;

		auto *target_data = grid[item[0]];
		if (target_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}

		Mom2(*target_data) = Mas2(*target_data) * source_value;
	}

	// pressures
	constexpr pamhd::mhd::Pressure P{};
	for (
		size_t i = 0;
		i < boundaries.get_number_of_value_boundaries(P);
		i++
	) {
		auto& value_bdy = boundaries.get_value_boundary(P, i);
		const auto& geometry_id = value_bdy.get_geometry_id();
		const auto& cells = bdy_geoms.get_cells(geometry_id);
		for (const auto& cell: cells) {
			const auto c = grid.geometry.get_center(cell);
			const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
			const auto
				lat = asin(c[2] / r),
				lon = atan2(c[1], c[0]);

			const auto pressure = value_bdy.get_data(
				simulation_time,
				c[0], c[1], c[2],
				r, lat, lon
			);

			auto* const cell_data = grid[cell];
			if (cell_data == nullptr) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << std::endl;
				abort();
			}

			if (Mas(*cell_data) > 0 and pressure > 0) {
				Nrj(*cell_data) = pamhd::mhd::get_total_energy_density(
					Mas(*cell_data),
					pamhd::mhd::get_velocity(Mom(*cell_data), Mas(*cell_data)),
					pressure,
					Magnetic_Field::data_type{0, 0, 0},
					adiabatic_index,
					vacuum_permeability
				);
			} else {
				Nrj(*cell_data) = 0;
			}
		}
	}
	for (const auto& item: boundaries.get_copy_boundary_cells(P)) {
		if (item.size() < 2) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}

		pamhd::mhd::Pressure::data_type source_value = 0.0;
		for (size_t i = 1; i < item.size(); i++) {
			auto* source_data = grid[item[i]];
			if (source_data == nullptr) {
				std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}

			source_value += pamhd::mhd::get_pressure(
				Mas(*source_data),
				Mom(*source_data),
				Nrj(*source_data),
				Mag(*source_data),
				adiabatic_index,
				vacuum_permeability
			);
		}
		source_value /= item.size() - 1;

		auto *target_data = grid[item[0]];
		if (target_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}

		if (Mas(*target_data) > 0 and source_value > 0) {
			Nrj(*target_data) = pamhd::mhd::get_total_energy_density(
				Mas(*target_data),
				pamhd::mhd::get_velocity(Mom(*target_data), Mas(*target_data)),
				source_value,
				Magnetic_Field::data_type{0, 0, 0},
				adiabatic_index,
				vacuum_permeability
			);
		} else {
			Nrj(*target_data) = 0;
		}
	}

	constexpr pamhd::mhd::Pressure2 P2{};
	for (
		size_t i = 0;
		i < boundaries.get_number_of_value_boundaries(P2);
		i++
	) {
		auto& value_bdy = boundaries.get_value_boundary(P2, i);
		const auto& geometry_id = value_bdy.get_geometry_id();
		const auto& cells = bdy_geoms.get_cells(geometry_id);
		for (const auto& cell: cells) {
			const auto c = grid.geometry.get_center(cell);
			const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
			const auto
				lat = asin(c[2] / r),
				lon = atan2(c[1], c[0]);

			const auto pressure = value_bdy.get_data(
				simulation_time,
				c[0], c[1], c[2],
				r, lat, lon
			);

			auto* const cell_data = grid[cell];
			if (cell_data == nullptr) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << std::endl;
				abort();
			}

			if (Mas2(*cell_data) > 0 and pressure > 0) {
				Nrj2(*cell_data) = pamhd::mhd::get_total_energy_density(
					Mas2(*cell_data),
					pamhd::mhd::get_velocity(Mom2(*cell_data), Mas2(*cell_data)),
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
	for (const auto& item: boundaries.get_copy_boundary_cells(P2)) {
		if (item.size() < 2) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}

		pamhd::mhd::Pressure::data_type source_value = 0.0;
		for (size_t i = 1; i < item.size(); i++) {
			auto* source_data = grid[item[i]];
			if (source_data == nullptr) {
				std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}

			source_value += pamhd::mhd::get_pressure(
				Mas2(*source_data),
				Mom2(*source_data),
				Nrj2(*source_data),
				Mag(*source_data),
				adiabatic_index,
				vacuum_permeability
			);
		}
		source_value /= item.size() - 1;

		auto *target_data = grid[item[0]];
		if (target_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}

		if (Mas2(*target_data) > 0 and source_value > 0) {
			Nrj2(*target_data) = pamhd::mhd::get_total_energy_density(
				Mas2(*target_data),
				pamhd::mhd::get_velocity(Mom2(*target_data), Mas2(*target_data)),
				source_value,
				Magnetic_Field::data_type{0, 0, 0},
				adiabatic_index,
				vacuum_permeability
			);
		} else {
			Nrj2(*target_data) = 0;
		}
	}

	// magnetic field
	constexpr pamhd::Magnetic_Field B{};
	for (
		size_t i = 0;
		i < boundaries.get_number_of_value_boundaries(B);
		i++
	) {
		auto& value_bdy = boundaries.get_value_boundary(B, i);
		const auto& geometry_id = value_bdy.get_geometry_id();
		const auto& cells = bdy_geoms.get_cells(geometry_id);
		for (const auto& cell: cells) {
			const auto c = grid.geometry.get_center(cell);
			const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
			const auto
				lat = asin(c[2] / r),
				lon = atan2(c[1], c[0]);

			const auto magnetic_field = value_bdy.get_data(
				simulation_time,
				c[0], c[1], c[2],
				r, lat, lon
			);

			auto* const cell_data = grid[cell];
			if (cell_data == nullptr) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << std::endl;
				abort();
			}

			Mag(*cell_data) = magnetic_field;

			// add magnetic field contribution to total energy densities
			const auto
				total_mass = Mas(*cell_data) + Mas2(*cell_data),
				mass_frac1 = Mas(*cell_data) / total_mass,
				mass_frac2 = Mas2(*cell_data) / total_mass;
			Nrj(*cell_data) += mass_frac1 * 0.5 * Mag(*cell_data).squaredNorm() / vacuum_permeability;
			Nrj2(*cell_data) += mass_frac2 * 0.5 * Mag(*cell_data).squaredNorm() / vacuum_permeability;
		}
	}
	for (const auto& item: boundaries.get_copy_boundary_cells(B)) {
		if (item.size() < 2) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}

		pamhd::Magnetic_Field::data_type source_value{0, 0, 0};
		for (size_t i = 1; i < item.size(); i++) {
			auto* source_data = grid[item[i]];
			if (source_data == nullptr) {
				std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}

			source_value += Mag(*source_data);
		}
		source_value /= item.size() - 1;

		auto *target_data = grid[item[0]];
		if (target_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}

		Mag(*target_data) = source_value;

		const auto
			total_mass = Mas(*target_data ) + Mas2(*target_data ),
			mass_frac1 = Mas(*target_data ) / total_mass,
			mass_frac2 = Mas2(*target_data ) / total_mass;
		Nrj(*target_data ) += mass_frac1 * 0.5 * Mag(*target_data ).squaredNorm() / vacuum_permeability;
		Nrj2(*target_data ) += mass_frac2 * 0.5 * Mag(*target_data ).squaredNorm() / vacuum_permeability;
	}
}


}} // namespaces


#endif // ifndef PAMHD_MHD_N_BOUNDARIES_HPP
