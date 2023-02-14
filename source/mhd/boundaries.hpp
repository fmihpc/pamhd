/*
Handles boundary logic of MHD part of PAMHD.

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

#ifndef PAMHD_MHD_BOUNDARIES_HPP
#define PAMHD_MHD_BOUNDARIES_HPP


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


// as set_solver_info() below but for magnetic field only
// TODO: remove magnetic field from set_solver_info()
template<
	class Solver_Info,
	class Grid,
	class Boundaries,
	class Boundary_Geometries,
	class Solver_Info_Getter
> void set_solver_info_magnetic(
	Grid& grid,
	Boundaries& boundaries,
	const Boundary_Geometries& geometries,
	const Solver_Info_Getter& Sol_Info
) {
	boundaries.classify(grid, geometries, Sol_Info);

	for (const auto& cell: grid.cells) {
		Sol_Info(*cell.data) = 0;
	}

	// magnetic field
	constexpr pamhd::Magnetic_Field B{};
	for (const auto& cell: boundaries.get_value_boundary_cells(B)) {
		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		Sol_Info(*cell_data) |= Solver_Info::magnetic_field_bdy;
	}
	for (const auto& item: boundaries.get_copy_boundary_cells(B)) {
		auto* const cell_data = grid[item[0]];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		Sol_Info(*cell_data) |= Solver_Info::magnetic_field_bdy;
	}
	const std::set<uint64_t> dont_solve_mag(
		boundaries.get_dont_solve_cells(B).cbegin(),
		boundaries.get_dont_solve_cells(B).cend()
	);

	// don't solve cells in which no variable is solved
	for (auto& cell: dont_solve_mag) {
		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		Sol_Info(*cell_data) |= Solver_Info::dont_solve;
	}

	grid.update_copies_of_remote_neighbors();
}


/*!
Prepares boundary information needed for MHD solver about each simulation cell.

MPI transfer of Sol_Info variable must be switched on before calling this.
*/
template<
	class Solver_Info,
	class Grid,
	class Boundaries,
	class Boundary_Geometries,
	class Solver_Info_Getter
> void set_solver_info(
	Grid& grid,
	Boundaries& boundaries,
	const Boundary_Geometries& geometries,
	const Solver_Info_Getter& Sol_Info
) {
	boundaries.classify(grid, geometries, Sol_Info);

	for (const auto& cell: grid.cells) {
		Sol_Info(*cell.data) = 0;
	}

	// number density
	constexpr pamhd::mhd::Number_Density N{};
	for (const auto& cell: boundaries.get_value_boundary_cells(N)) {
		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		Sol_Info(*cell_data) |= Solver_Info::mass_density_bdy;
	}
	for (const auto& item: boundaries.get_copy_boundary_cells(N)) {
		auto* const cell_data = grid[item[0]];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		Sol_Info(*cell_data) |= Solver_Info::mass_density_bdy;
	}
	const std::set<uint64_t> dont_solve_mass(
		boundaries.get_dont_solve_cells(N).cbegin(),
		boundaries.get_dont_solve_cells(N).cend()
	);

	// velocity
	constexpr pamhd::mhd::Velocity V{};
	for (const auto& cell: boundaries.get_value_boundary_cells(V)) {
		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		Sol_Info(*cell_data) |= Solver_Info::velocity_bdy;
	}
	for (const auto& item: boundaries.get_copy_boundary_cells(V)) {
		auto* const cell_data = grid[item[0]];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		Sol_Info(*cell_data) |= Solver_Info::velocity_bdy;
	}
	const std::set<uint64_t> dont_solve_velocity(
		boundaries.get_dont_solve_cells(V).cbegin(),
		boundaries.get_dont_solve_cells(V).cend()
	);

	// pressure
	constexpr pamhd::mhd::Pressure P{};
	for (const auto& cell: boundaries.get_value_boundary_cells(P)) {
		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		Sol_Info(*cell_data) |= Solver_Info::pressure_bdy;
	}
	for (const auto& item: boundaries.get_copy_boundary_cells(P)) {
		auto* const cell_data = grid[item[0]];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		Sol_Info(*cell_data) |= Solver_Info::pressure_bdy;
	}
	const std::set<uint64_t> dont_solve_pressure(
		boundaries.get_dont_solve_cells(P).cbegin(),
		boundaries.get_dont_solve_cells(P).cend()
	);

	// magnetic field
	constexpr pamhd::Magnetic_Field B{};
	for (const auto& cell: boundaries.get_value_boundary_cells(B)) {
		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		Sol_Info(*cell_data) |= Solver_Info::magnetic_field_bdy;
	}
	for (const auto& item: boundaries.get_copy_boundary_cells(B)) {
		auto* const cell_data = grid[item[0]];
		if (cell_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		Sol_Info(*cell_data) |= Solver_Info::magnetic_field_bdy;
	}
	const std::set<uint64_t> dont_solve_mag(
		boundaries.get_dont_solve_cells(B).cbegin(),
		boundaries.get_dont_solve_cells(B).cend()
	);

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
		Sol_Info(*cell_data) |= Solver_Info::dont_solve;
	}

	grid.update_copies_of_remote_neighbors();
}


// as apply_boundaries() below but for magnetic field only.
template<
	class Grid,
	class Boundaries,
	class Boundary_Geometries,
	class Magnetic_Field_Pos_Getter,
	class Magnetic_Field_Neg_Getter
> void apply_magnetic_field_boundaries(
	Grid& grid,
	Boundaries& boundaries,
	const Boundary_Geometries& bdy_geoms,
	const double simulation_time,
	const Magnetic_Field_Pos_Getter& Magp,
	const Magnetic_Field_Neg_Getter& /*Magn*/,
	const bool staggered
) {
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
			const auto c = [&](){
				auto c_ = grid.geometry.get_center(cell);
				if (staggered) {
					const auto len = grid.geometry.get_length(cell);
					c_[0] += len[0]/2;
					c_[1] += len[1]/2;
					c_[2] += len[2]/2;
				}
				return c_;
			}();
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

			Magp(*cell_data) = magnetic_field;
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

			source_value += Magp(*source_data);
		}
		source_value /= item.size() - 1;

		auto *target_data = grid[item[0]];
		if (target_data == nullptr) {
			std::cerr <<  __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}

		Magp(*target_data) = source_value;
	}
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
	class Grid,
	class Boundaries,
	class Boundary_Geometries,
	class Mass_Getter,
	class Momentum_Getter,
	class Energy_Getter,
	class Magnetic_Field_Getter
> void apply_fluid_boundaries(
	Grid& grid,
	Boundaries& boundaries,
	const Boundary_Geometries& bdy_geoms,
	const double simulation_time,
	const Mass_Getter& Mas,
	const Momentum_Getter& Mom,
	const Energy_Getter& Nrj,
	const Magnetic_Field_Getter& Mag,
	const double proton_mass,
	const double adiabatic_index,
	const double vacuum_permeability
) {
	// number density
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

	// velocity
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

	// pressure
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

			Nrj(*cell_data) = pamhd::mhd::get_total_energy_density(
				Mas(*cell_data),
				pamhd::mhd::get_velocity(
					Mom(*cell_data),
					Mas(*cell_data)
				),
				pressure,
				Mag(*cell_data),
				adiabatic_index,
				vacuum_permeability
			);
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

		Nrj(*target_data) = pamhd::mhd::get_total_energy_density(
			Mas(*target_data),
			pamhd::mhd::get_velocity(
				Mom(*target_data),
				Mas(*target_data)
			),
			source_value,
			Mag(*target_data),
			adiabatic_index,
			vacuum_permeability
		);
	}
}

}} // namespaces


#endif // ifndef PAMHD_MHD_BOUNDARIES_HPP
