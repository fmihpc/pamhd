/*
Initializes staggered MHD solution of PAMHD.

Copyright 2014, 2015, 2016, 2017 Ilja Honkonen
Copyright 2019, 2023, 2024 Finnish Meteorological Institute
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
	class Face_Magnetic_Field_Getter,
	class Magnetic_Field_Flux_Getters,
	class Background_Magnetic_Field_Getter,
	class Primary_Face_Getter
> void initialize_magnetic_field_staggered(
	const Geometries& geometries,
	Init_Cond& initial_conditions,
	const Background_Magnetic_Field& bg_B,
	Grid& grid,
	const double time,
	const double vacuum_permeability,
	const Face_Magnetic_Field_Getter Face_B,
	const Magnetic_Field_Flux_Getters Mag_f,
	const Background_Magnetic_Field_Getter Bg_B,
	const Primary_Face_Getter PFace
) {
	using std::get;
	using std::make_tuple;

	const auto
		Mag_fnx = get<0>(Mag_f), Mag_fpx = get<1>(Mag_f),
		Mag_fny = get<2>(Mag_f), Mag_fpy = get<3>(Mag_f),
		Mag_fnz = get<4>(Mag_f), Mag_fpz = get<5>(Mag_f);

	// set default magnetic field
	for (const auto& cell: grid.local_cells()) {
		Mag_fnx(*cell.data) =
		Mag_fny(*cell.data) =
		Mag_fnz(*cell.data) =
		Mag_fpx(*cell.data) =
		Mag_fpy(*cell.data) =
		Mag_fpz(*cell.data) = {0, 0, 0};

		const auto [rx, ry, rz] = grid.geometry.get_center(cell.id);
		const auto [sx, sy, sz] = grid.geometry.get_min(cell.id);
		const auto [ex, ey, ez] = grid.geometry.get_max(cell.id);

		Bg_B(*cell.data)(0, -1) = bg_B.get_background_field(
			{sx, ry, rz},
			vacuum_permeability
		);
		Bg_B(*cell.data)(0, +1) = bg_B.get_background_field(
			{ex, ry, rz},
			vacuum_permeability
		);
		Bg_B(*cell.data)(1, -1) = bg_B.get_background_field(
			{rx, sy, rz},
			vacuum_permeability
		);
		Bg_B(*cell.data)(1, +1) = bg_B.get_background_field(
			{rx, ey, rz},
			vacuum_permeability
		);
		Bg_B(*cell.data)(2, -1) = bg_B.get_background_field(
			{rx, ry, sz},
			vacuum_permeability
		);
		Bg_B(*cell.data)(2, +1) = bg_B.get_background_field(
			{rx, ry, ez},
			vacuum_permeability
		);

		const auto& pface = PFace(*cell.data);
		for (const int dir: {-1,+1,-2,+2,-3,+3}) {
			if (not pface(dir)) continue;

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
					std::cout << __FILE__ "(" << __LINE__ << ")" << std::endl;
					abort();
					break;
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
			Face_B(*cell.data)(dir) = magnetic_field[dim];
		}
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
		for (const auto& cell_id: cells) {
			auto* const cell_data = grid[cell_id];
			if (cell_data == nullptr) {
				std::cerr <<  __FILE__ << "(" << __LINE__
					<< ") No data for cell: " << cell_id
					<< std::endl;
				abort();
			}
			const auto& pface = PFace(*cell_data);

			const auto [rx, ry, rz] = grid.geometry.get_center(cell_id);
			const auto [sx, sy, sz] = grid.geometry.get_min(cell_id);
			const auto [ex, ey, ez] = grid.geometry.get_max(cell_id);
			for (const int dir: {-1,+1,-2,+2,-3,+3}) {
				if (not pface(dir)) continue;

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
						std::cout << __FILE__ "(" << __LINE__ << ")" << std::endl;
						abort();
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
				Face_B(*cell_data)(dir) = magnetic_field[dim];
			}
		}
	}
}

}} // namespaces

#endif // ifndef PAMHD_MHD_INITIALIZE_STAGGERED_HPP
