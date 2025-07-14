/*
Initializes staggered MHD solution of PAMHD.

Copyright 2014, 2015, 2016, 2017 Ilja Honkonen
Copyright 2019, 2023, 2024, 2025 Finnish Meteorological Institute
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

#ifndef PAMHD_MHD_INITIALIZE_STAGGERED_HPP
#define PAMHD_MHD_INITIALIZE_STAGGERED_HPP


#include "cmath"
#include "iostream"
#include "limits"
#include "string"
#include "tuple"

#include "dccrg.hpp"

#include "grid/amr.hpp"
#include "mhd/common.hpp"
#include "mhd/variables.hpp"


namespace pamhd {
namespace mhd {


/*! Sets initial state of magnetic field and zeroes fluxes

Emulates 1d/2d system if grid length is initially 1 cell
in some dimension(s), even if it has been refined, in which
case some cells will overlap in affected dimension(s).
*/
template <
	class Boundary_Magnetic_Field,
	class Geometries,
	class Init_Cond,
	class Background_Magnetic_Field,
	class Grid,
	class Face_Magnetic_Field_Getter,
	class Magnetic_Field_Flux_Getters,
	class Background_Magnetic_Field_Getter
> void initialize_magnetic_field_staggered(
	const Geometries& geometries,
	Init_Cond& initial_conditions,
	const Background_Magnetic_Field& bg_B,
	Grid& grid,
	const double& time,
	const double& vacuum_permeability,
	const Face_Magnetic_Field_Getter& Face_B,
	const Magnetic_Field_Flux_Getters& Mag_f,
	const Background_Magnetic_Field_Getter& Bg_B
) {
	using std::get;
	using std::make_tuple;
	using std::runtime_error;
	using std::to_string;

	// set default magnetic field
	for (const auto& cell: grid.local_cells()) {
		Mag_f(*cell.data, -1) =
		Mag_f(*cell.data, +1) =
		Mag_f(*cell.data, -2) =
		Mag_f(*cell.data, +2) =
		Mag_f(*cell.data, -3) =
		Mag_f(*cell.data, +3) = {0, 0, 0};

		const auto [
			center, start, end
		] = pamhd::grid::get_cell_geom_emulated(grid, cell.id);
		const auto& [rx, ry, rz] = center;
		const auto& [sx, sy, sz] = start;
		const auto& [ex, ey, ez] = end;

		Bg_B.data(*cell.data)(-1) = bg_B.get_background_field(
			{sx, ry, rz},
			vacuum_permeability
		);
		Bg_B.data(*cell.data)(+1) = bg_B.get_background_field(
			{ex, ry, rz},
			vacuum_permeability
		);
		Bg_B.data(*cell.data)(-2) = bg_B.get_background_field(
			{rx, sy, rz},
			vacuum_permeability
		);
		Bg_B.data(*cell.data)(+2) = bg_B.get_background_field(
			{rx, ey, rz},
			vacuum_permeability
		);
		Bg_B.data(*cell.data)(-3) = bg_B.get_background_field(
			{rx, ry, sz},
			vacuum_permeability
		);
		Bg_B.data(*cell.data)(+3) = bg_B.get_background_field(
			{rx, ry, ez},
			vacuum_permeability
		);

		for (const int dir: {-3,-2,-1,+1,+2,+3}) {
			const auto [x, y, z] = [&](){
				// center of face
				switch (dir) {
				case -1:
					return make_tuple(sx, ry, rz);
					break;
				case +1:
					return make_tuple(ex, ry, rz);
					break;
				case -2:
					return make_tuple(rx, sy, rz);
					break;
				case +2:
					return make_tuple(rx, ey, rz);
					break;
				case -3:
					return make_tuple(rx, ry, sz);
					break;
				case +3:
					return make_tuple(rx, ry, ez);
					break;
				default:
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
				}
			}();
			const auto
				r = sqrt(x*x + y*y + z*z),
				lat = asin(z / r),
				lon = atan2(y, x);

			const auto magnetic_field
				= initial_conditions.get_default_data(
					Boundary_Magnetic_Field(),
					time,
					x, y, z,
					r, lat, lon
				);

			const size_t dim = std::abs(dir) - 1;
			Face_B.data(*cell.data)(dir) = magnetic_field[dim];
		}
	}
	Bg_B.type().is_stale = true;
	Face_B.type().is_stale = true;

	// set non-default magnetic field
	for (
		size_t i = 0;
		i < initial_conditions.get_number_of_regions(Boundary_Magnetic_Field());
		i++
	) {
		const auto& init_cond = initial_conditions.get_initial_condition(Boundary_Magnetic_Field(), i);
		const auto& geometry_id = init_cond.get_geometry_id();
		const auto& cells = geometries.get_cells(geometry_id);
		for (const auto& cell_id: cells) {
			auto* const cell_data = grid[cell_id];
			if (cell_data == nullptr) {
				throw runtime_error(__FILE__ "(" + to_string(__LINE__)
					+ ") No data for cell: " + to_string(cell_id));
			}

			const auto [
				center, start, end
			] = pamhd::grid::get_cell_geom_emulated(grid, cell_id);
			const auto& [rx, ry, rz] = center;
			const auto& [sx, sy, sz] = start;
			const auto& [ex, ey, ez] = end;
			for (const int dir: {-3,-2,-1,+1,+2,+3}) {
				// center of face
				const auto [x, y, z] = [&](){
					switch (dir) {
					case -1:
						return make_tuple(sx, ry, rz);
						break;
					case +1:
						return make_tuple(ex, ry, rz);
						break;
					case -2:
						return make_tuple(rx, sy, rz);
						break;
					case +2:
						return make_tuple(rx, ey, rz);
						break;
					case -3:
						return make_tuple(rx, ry, sz);
						break;
					case +3:
						return make_tuple(rx, ry, ez);
						break;
					default:
						throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					}
				}();
				const auto
					r = sqrt(x*x + y*y + z*z),
					lat = asin(z / r),
					lon = atan2(y, x);

				const auto magnetic_field
					= initial_conditions.get_default_data(
						Boundary_Magnetic_Field(),
						time,
						x, y, z,
						r, lat, lon
					);

				const size_t dim = std::abs(dir) - 1;
				Face_B.data(*cell_data)(dir) = magnetic_field[dim];
			}
		}
	}
}

}} // namespaces

#endif // ifndef PAMHD_MHD_INITIALIZE_STAGGERED_HPP
