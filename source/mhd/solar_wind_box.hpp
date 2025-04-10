/*
MHD parts of solar wind box program.

Copyright 2024, 2025 Finnish Meteorological Institute
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

#ifndef PAMHD_MHD_SW_BOX_HPP
#define PAMHD_MHD_SW_BOX_HPP


#include "cmath"
#include "stdexcept"
#include "string"
#include "vector"

#include "dccrg.hpp"
#include "rapidjson/document.h"

#include "grid/amr.hpp"
#include "mhd/common.hpp"
#include "mhd/solve_staggered.hpp"
#include "mhd/variables.hpp"
#include "solar_wind_box_options.hpp"
#include "variables.hpp"


namespace pamhd {
namespace mhd {


template<
	class Grid,
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Volume_Magnetic_Field_Getter,
	class Face_Magnetic_Field_Getter,
	class Face_Difference_Magnetic_Field_Getter,
	class Background_Magnetic_Field_Getter
> void initialize_plasma(
	Grid& grid,
	const double& /*sim_time*/,
	const Solar_Wind_Box_Options& options,
	const pamhd::Background_Magnetic_Field<
		double,
		pamhd::Magnetic_Field::data_type
	>& background_B,
	const Mass_Density_Getter& Mas,
	const Momentum_Density_Getter& Mom,
	const Total_Energy_Density_Getter& Nrj,
	const Volume_Magnetic_Field_Getter& Vol_B,
	const Face_Magnetic_Field_Getter& Face_B,
	const Face_Difference_Magnetic_Field_Getter& Face_dB,
	const Background_Magnetic_Field_Getter& Bg_B,
	const double& adiabatic_index,
	const double& vacuum_permeability,
	const double& proton_mass
) try {
	using std::string;
	using std::invalid_argument;
	using std::runtime_error;
	using std::to_string;

	for (const size_t dim: {0, 1, 2}) {
		if (
			options.sw_dim == dim
			and options.sw_magnetic_field[dim] != 0
		) {
			throw runtime_error(
				__FILE__ "(" + to_string(__LINE__) + "): solar "
				"wind boundary in " + to_string(dim) + " dimension"
				" but corresponding B component non-zero: "
				+ to_string(options.sw_magnetic_field[dim]));
		}
	}

	if (grid.get_rank() == 0) {
		std::cout << "Initializing run, solar wind: "
		<< options.sw_nr_density << " #/m^3, "
		<< options.sw_velocity << " m/s, "
		<< options.sw_pressure << " Pa, "
		<< options.sw_magnetic_field << " T"
		<< std::endl;
	}

	const auto lvl0 = grid.mapping.length.get();
	const auto grid_end = grid.geometry.get_end();
	for (const auto& cell: grid.local_cells()) {
		const auto [rx, ry, rz] = grid.geometry.get_center(cell.id);
		const auto [sx, sy, sz] = grid.geometry.get_min(cell.id);
		const auto [ex, ey, ez] = grid.geometry.get_max(cell.id);

		if (lvl0[0] > 1) {
			Bg_B.data(*cell.data)(-1) = background_B.get_background_field(
				{sx, ry, rz},
				vacuum_permeability
			);
			Bg_B.data(*cell.data)(+1) = background_B.get_background_field(
				{ex, ry, rz},
				vacuum_permeability
			);
		} else {
			Bg_B.data(*cell.data)(-1) =
			Bg_B.data(*cell.data)(+1) =
				background_B.get_background_field(
					{rx, ry, rz},
					vacuum_permeability
				);
		}

		if (lvl0[1] > 1) {
			Bg_B.data(*cell.data)(-2) = background_B.get_background_field(
				{rx, sy, rz},
				vacuum_permeability
			);
			Bg_B.data(*cell.data)(+2) = background_B.get_background_field(
				{rx, ey, rz},
				vacuum_permeability
			);
		} else {
			Bg_B.data(*cell.data)(-2) =
			Bg_B.data(*cell.data)(+2) =
				background_B.get_background_field(
					{rx, ry, rz},
					vacuum_permeability
				);
		}

		if (lvl0[2] > 1) {
			Bg_B.data(*cell.data)(-3) = background_B.get_background_field(
				{rx, ry, sz},
				vacuum_permeability
			);
			Bg_B.data(*cell.data)(+3) = background_B.get_background_field(
				{rx, ry, ez},
				vacuum_permeability
			);
		} else {
			Bg_B.data(*cell.data)(-3) =
			Bg_B.data(*cell.data)(+3) =
				background_B.get_background_field(
					{rx, ry, rz},
					vacuum_permeability
				);
		}

		const auto mass = options.sw_nr_density * proton_mass;
		Mas.data(*cell.data) = mass;

		const auto cell_center = grid.geometry.get_center(cell.id);
		const auto v_factor = std::max(0.0, cell_center[0]/grid_end[0]);
		const Eigen::Vector3d velocity{
			v_factor * options.sw_velocity[0],
			v_factor * options.sw_velocity[1],
			v_factor * options.sw_velocity[2]};
		Mom.data(*cell.data) = mass * velocity;

		Vol_B.data(*cell.data)[0]   =
		Face_B.data(*cell.data)(-1) =
		Face_B.data(*cell.data)(+1) = options.sw_magnetic_field[0];
		Vol_B.data(*cell.data)[1]   =
		Face_B.data(*cell.data)(-2) =
		Face_B.data(*cell.data)(+2) = options.sw_magnetic_field[1];
		Vol_B.data(*cell.data)[2]   =
		Face_B.data(*cell.data)(-3) =
		Face_B.data(*cell.data)(+3) = options.sw_magnetic_field[2];
		Face_dB.data(*cell.data) = {0, 0, 0, 0, 0, 0};

		Nrj.data(*cell.data) = pamhd::mhd::get_total_energy_density(
			Mas.data(*cell.data),
			velocity,
			options.sw_pressure,
			Vol_B.data(*cell.data),
			adiabatic_index,
			vacuum_permeability
		);
	}
	Mas.type().is_stale = true;
	Mom.type().is_stale = true;
	Nrj.type().is_stale = true;
	Vol_B.type().is_stale = true;
	Face_B.type().is_stale = true;
	Bg_B.type().is_stale = true;
} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


template<
	class Cells,
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Energy_Density_Getter,
	class Volume_Magnetic_Field_Getter,
	class Face_Magnetic_Field_Getter,
	class Face_Magnetic_Field_Change_Getter
> void apply_solar_wind_boundaries(
	const Solar_Wind_Box_Options& options,
	const Cells& solar_wind_cells,
	const double& /*sim_time*/,
	const double& adiabatic_index,
	const double& vacuum_permeability,
	const double& proton_mass,
	const Mass_Density_Getter& Mas,
	const Momentum_Density_Getter& Mom,
	const Energy_Density_Getter& Nrj,
	const Volume_Magnetic_Field_Getter& Vol_B,
	const Face_Magnetic_Field_Getter& Face_B,
	const Face_Magnetic_Field_Change_Getter& Face_dB
) try {
	for (const auto& cell: solar_wind_cells) {
		if (cell.is_local) {
			const auto mass = options.sw_nr_density * proton_mass;
			Mas.data(*cell.data) = mass;
			Mom.data(*cell.data) = {
				mass*options.sw_velocity[0],
				mass*options.sw_velocity[1],
				mass*options.sw_velocity[2]
			};
			Face_dB.data(*cell.data) = {0, 0, 0, 0, 0, 0};

			Vol_B.data(*cell.data)[0]   =
			Face_B.data(*cell.data)(-1) =
			Face_B.data(*cell.data)(+1) = options.sw_magnetic_field[0];
			Vol_B.data(*cell.data)[1]   =
			Face_B.data(*cell.data)(-2) =
			Face_B.data(*cell.data)(+2) = options.sw_magnetic_field[1];
			Vol_B.data(*cell.data)[2]   =
			Face_B.data(*cell.data)(-3) =
			Face_B.data(*cell.data)(+3) = options.sw_magnetic_field[2];
			Nrj.data(*cell.data) = pamhd::mhd::get_total_energy_density(
				mass, options.sw_velocity,
				options.sw_pressure, Vol_B.data(*cell.data),
				adiabatic_index, vacuum_permeability
			);
		}

		for (const auto& neighbor: cell.neighbors_of) {
			if (not neighbor.is_local) continue;
			if (options.sw_dir != -neighbor.face_neighbor) continue;

			const auto& mas = Mas.data(*neighbor.data);
			const auto& mom = Mom.data(*neighbor.data);
			const auto velocity = mom / mas;
			const auto pressure = pamhd::mhd::get_pressure(
				mas, mom,
				Nrj.data(*neighbor.data),
				Vol_B.data(*neighbor.data),
				adiabatic_index, vacuum_permeability
			);

			Face_B.data(*neighbor.data)(options.sw_dir)
				= options.sw_magnetic_field[options.sw_dim];
			Vol_B.data(*neighbor.data)[options.sw_dim]
				= Face_B.data(*neighbor.data)(+options.sw_dir) / 2
				+ Face_B.data(*neighbor.data)(-options.sw_dir) / 2;
			Nrj.data(*neighbor.data) = pamhd::mhd::get_total_energy_density(
				mas, velocity,
				pressure, Vol_B.data(*neighbor.data),
				adiabatic_index, vacuum_permeability
			);
		}
	}
	Mas.type().is_stale = true;
	Mom.type().is_stale = true;
	Nrj.type().is_stale = true;
	Vol_B.type().is_stale = true;
	Face_B.type().is_stale = true;

} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


template<
	class Grid,
	class SW_Cells,
	class Face_Cells,
	class Edge_Cells,
	class Vert_Cells,
	class Planet_Cells,
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Energy_Density_Getter,
	class Volume_Magnetic_Field_Getter,
	class Face_Magnetic_Field_Getter,
	class Face_Magnetic_Field_Change_Getter,
	class Solver_Info_Getter
> void apply_boundaries_sw_box(
	Grid& grid,
	const double& sim_time,
	const Solar_Wind_Box_Options& options,
	const SW_Cells& solar_wind_cells,
	const Face_Cells& face_cells,
	const Edge_Cells& edge_cells,
	const Vert_Cells& vert_cells,
	const Planet_Cells& planet_cells,
	const double& adiabatic_index,
	const double& vacuum_permeability,
	const double& proton_mass,
	const Mass_Density_Getter& Mas,
	const Momentum_Density_Getter& Mom,
	const Energy_Density_Getter& Nrj,
	const Volume_Magnetic_Field_Getter& Vol_B,
	const Face_Magnetic_Field_Getter& Face_B,
	const Face_Magnetic_Field_Change_Getter& Face_dB,
	const Solver_Info_Getter& SInfo
) try {
	using std::abs;
	using std::runtime_error;
	using std::to_string;

	using Cell = Grid::cell_data_type;

	bool update_copies = false;
	if (
		Mas.type().is_stale
		or Mom.type().is_stale
		or Nrj.type().is_stale
		or Vol_B.type().is_stale
	) {
		update_copies = true;
		Cell::set_transfer_all(true,
			Mas.type(), Mom.type(), Nrj.type(), Vol_B.type());
	}
	if (Face_B.type().is_stale) {
		update_copies = true;
		Cell::set_transfer_all(true, Face_B.type());
	}
	if (SInfo.type().is_stale) {
		update_copies = true;
		Cell::set_transfer_all(true, SInfo.type());
	}
	if (update_copies) {
		grid.update_copies_of_remote_neighbors();
	}
	Cell::set_transfer_all(false,
		Mas.type(), Mom.type(), Nrj.type(),
		Vol_B.type(), Face_B.type(), SInfo.type());

	Mas.type().is_stale = true;
	Mom.type().is_stale = true;
	Nrj.type().is_stale = true;
	Vol_B.type().is_stale = true;
	Face_B.type().is_stale = true;

	// boundary and normal cell share face
	for (int dir: {-3,-2,-1,+1,+2,+3}) {
		for (const auto& cell: face_cells(dir)) {
			Face_dB.data(*cell.data) = {0, 0, 0, 0, 0, 0};
			double pressure = -1;
			for (const auto& neighbor: cell.neighbors_of) {
				const auto& fn = neighbor.face_neighbor;
				if (fn != -dir) continue;
				if (SInfo.data(*neighbor.data) != 1) {
					throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
				}
				if (pressure > 0) {
					throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
				}
				pressure = pamhd::mhd::get_pressure(
					Mas.data(*neighbor.data),
					Mom.data(*neighbor.data),
					Nrj.data(*neighbor.data),
					Vol_B.data(*neighbor.data),
					adiabatic_index, vacuum_permeability
				);
				Mas.data(*cell.data) = Mas.data(*neighbor.data);
				Mom.data(*cell.data) = Mom.data(*neighbor.data);

				for (int dir2: {-3,-2,-1,+1,+2,+3}) {
					if (dir2 == fn or dir2 == -fn) {
						Face_B.data(*cell.data)(dir2) = Face_B.data(*neighbor.data)(-fn);
					} else {
						Face_B.data(*cell.data)(dir2) = Face_B.data(*neighbor.data)(dir2);
					}
				}
			}
			Vol_B.data(*cell.data) = {
				0.5*(Face_B.data(*cell.data)(-1) + Face_B.data(*cell.data)(+1)),
				0.5*(Face_B.data(*cell.data)(-2) + Face_B.data(*cell.data)(+2)),
				0.5*(Face_B.data(*cell.data)(-3) + Face_B.data(*cell.data)(+3))
			};
			Nrj.data(*cell.data) = pamhd::mhd::get_total_energy_density(
				Mas.data(*cell.data),
				Mom.data(*cell.data) / Mas.data(*cell.data),
				pressure,
				Vol_B.data(*cell.data),
				adiabatic_index,
				vacuum_permeability
			);
		}
	}

	// boundary and normal cell share edge
	for (int dim: {0, 1, 2})
	for (int dir1: {-1, +1})
	for (int dir2: {-1, +1}) {
		for (const auto& cell: edge_cells(dim, dir1, dir2)) {
			Face_dB.data(*cell.data) = {0, 0, 0, 0, 0, 0};
			double pressure = -1;
			for (const auto& neighbor: cell.neighbors_of) {
				const auto& en = neighbor.edge_neighbor;
				if (en[0] != dim or en[1] != -dir1 or en[2] != -dir2) continue;
				if (SInfo.data(*neighbor.data) != 1) {
					throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
				}
				if (pressure > 0) {
					throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
				}
				pressure = pamhd::mhd::get_pressure(
					Mas.data(*neighbor.data),
					Mom.data(*neighbor.data),
					Nrj.data(*neighbor.data),
					Vol_B.data(*neighbor.data),
					adiabatic_index, vacuum_permeability
				);
				Mas.data(*cell.data) = Mas.data(*neighbor.data);
				Mom.data(*cell.data) = Mom.data(*neighbor.data);

				if (dim == 0) {
					Face_B.data(*cell.data)(-1) = Face_B.data(*neighbor.data)(-1);
					Face_B.data(*cell.data)(+1) = Face_B.data(*neighbor.data)(+1);
					if (dir1 < 0) {
						Face_B.data(*cell.data)(-2) =
						Face_B.data(*cell.data)(+2) = Face_B.data(*neighbor.data)(-2);
					} else {
						Face_B.data(*cell.data)(-2) =
						Face_B.data(*cell.data)(+2) = Face_B.data(*neighbor.data)(+2);
					}
					if (dir2 < 0) {
						Face_B.data(*cell.data)(-3) =
						Face_B.data(*cell.data)(+3) = Face_B.data(*neighbor.data)(-3);
					} else {
						Face_B.data(*cell.data)(-3) =
						Face_B.data(*cell.data)(+3) = Face_B.data(*neighbor.data)(+3);
					}
				}
				if (dim == 1) {
					Face_B.data(*cell.data)(-2) = Face_B.data(*neighbor.data)(-2);
					Face_B.data(*cell.data)(+2) = Face_B.data(*neighbor.data)(+2);
					if (dir1 < 0) {
						Face_B.data(*cell.data)(-1) =
						Face_B.data(*cell.data)(+1) = Face_B.data(*neighbor.data)(-1);
					} else {
						Face_B.data(*cell.data)(-1) =
						Face_B.data(*cell.data)(+1) = Face_B.data(*neighbor.data)(+1);
					}
					if (dir2 < 0) {
						Face_B.data(*cell.data)(-3) =
						Face_B.data(*cell.data)(+3) = Face_B.data(*neighbor.data)(-3);
					} else {
						Face_B.data(*cell.data)(-3) =
						Face_B.data(*cell.data)(+3) = Face_B.data(*neighbor.data)(+3);
					}
				}
				if (dim == 2) {
					Face_B.data(*cell.data)(-3) = Face_B.data(*neighbor.data)(-3);
					Face_B.data(*cell.data)(+3) = Face_B.data(*neighbor.data)(+3);
					if (dir1 < 0) {
						Face_B.data(*cell.data)(-1) =
						Face_B.data(*cell.data)(+1) = Face_B.data(*neighbor.data)(-1);
					} else {
						Face_B.data(*cell.data)(-1) =
						Face_B.data(*cell.data)(+1) = Face_B.data(*neighbor.data)(+1);
					}
					if (dir2 < 0) {
						Face_B.data(*cell.data)(-2) =
						Face_B.data(*cell.data)(+2) = Face_B.data(*neighbor.data)(-2);
					} else {
						Face_B.data(*cell.data)(-2) =
						Face_B.data(*cell.data)(+2) = Face_B.data(*neighbor.data)(+2);
					}
				}
			}
			Vol_B.data(*cell.data) = {
				0.5*(Face_B.data(*cell.data)(-1) + Face_B.data(*cell.data)(+1)),
				0.5*(Face_B.data(*cell.data)(-2) + Face_B.data(*cell.data)(+2)),
				0.5*(Face_B.data(*cell.data)(-3) + Face_B.data(*cell.data)(+3))
			};
			Nrj.data(*cell.data) = pamhd::mhd::get_total_energy_density(
				Mas.data(*cell.data),
				Mom.data(*cell.data) / Mas.data(*cell.data),
				pressure,
				Vol_B.data(*cell.data),
				adiabatic_index,
				vacuum_permeability
			);
		}
	}

	// boundary and normal cell share vertex
	for (const auto& cell: vert_cells) {
		Face_dB.data(*cell.data) = {0, 0, 0, 0, 0, 0};
		const auto cilen = grid.mapping.get_cell_length_in_indices(cell.id);
		double pressure = -1;
		for (const auto& neighbor: cell.neighbors_of) {
			const auto& fn = neighbor.face_neighbor;
			if (fn != 0) continue;
			const auto& en = neighbor.edge_neighbor;
			if (en[0] >= 0) continue;
			if (
				abs(neighbor.x) > cilen
				or abs(neighbor.y) > cilen
				or abs(neighbor.z) > cilen
			) continue;
			if (neighbor.relative_size != 0) {
				throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			}
			if (SInfo.data(*neighbor.data) != 1) {
				throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			}

			if (pressure > 0) {
				throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			}
			pressure = pamhd::mhd::get_pressure(
				Mas.data(*neighbor.data),
				Mom.data(*neighbor.data),
				Nrj.data(*neighbor.data),
				Vol_B.data(*neighbor.data),
				adiabatic_index, vacuum_permeability
			);
			Mas.data(*cell.data) = Mas.data(*neighbor.data);
			Mom.data(*cell.data) = Mom.data(*neighbor.data);

			if (neighbor.x < 0) {
				Face_B.data(*cell.data)(-1) =
				Face_B.data(*cell.data)(+1) = Face_B.data(*neighbor.data)(+1);
			} else {
				Face_B.data(*cell.data)(-1) =
				Face_B.data(*cell.data)(+1) = Face_B.data(*neighbor.data)(-1);
			}
			if (neighbor.y < 0) {
				Face_B.data(*cell.data)(-2) =
				Face_B.data(*cell.data)(+2) = Face_B.data(*neighbor.data)(+2);
			} else {
				Face_B.data(*cell.data)(-2) =
				Face_B.data(*cell.data)(+2) = Face_B.data(*neighbor.data)(-2);
			}
			if (neighbor.z < 0) {
				Face_B.data(*cell.data)(-3) =
				Face_B.data(*cell.data)(+3) = Face_B.data(*neighbor.data)(+3);
			} else {
				Face_B.data(*cell.data)(-3) =
				Face_B.data(*cell.data)(+3) = Face_B.data(*neighbor.data)(-3);
			}
		}
		Vol_B.data(*cell.data) = {
			0.5*(Face_B.data(*cell.data)(-1) + Face_B.data(*cell.data)(+1)),
			0.5*(Face_B.data(*cell.data)(-2) + Face_B.data(*cell.data)(+2)),
			0.5*(Face_B.data(*cell.data)(-3) + Face_B.data(*cell.data)(+3))
		};
		Nrj.data(*cell.data) = pamhd::mhd::get_total_energy_density(
			Mas.data(*cell.data),
			Mom.data(*cell.data) / Mas.data(*cell.data),
			pressure,
			Vol_B.data(*cell.data),
			adiabatic_index,
			vacuum_permeability
		);
	}

	// planetary boundary cells
	for (const auto& cell: planet_cells) {
		Face_dB.data(*cell.data) = {0, 0, 0, 0, 0, 0};

		for (const int dir: {-3,-2,-1,+1,+2,+3}) {
			Face_B.data(*cell.data)(dir) = 0;
		}
		Vol_B.data(*cell.data) = {0, 0, 0};
		Mas.data(*cell.data) = options.sw_nr_density * proton_mass;
		Mom.data(*cell.data) = {0,0,0};
		Nrj.data(*cell.data) = pamhd::mhd::get_total_energy_density(
			Mas.data(*cell.data),
			Mom.data(*cell.data),
			options.sw_pressure,
			Vol_B.data(*cell.data),
			adiabatic_index,
			vacuum_permeability
		);
	}

	apply_solar_wind_boundaries(
		options, solar_wind_cells, sim_time,
		adiabatic_index, vacuum_permeability, proton_mass,
		Mas, Mom, Nrj, Vol_B, Face_B, Face_dB
	);

} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


}} // namespaces


#endif // ifndef PAMHD_MHD_SW_BOX_HPP
