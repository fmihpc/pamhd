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


#include "algorithm"
#include "cmath"
#include "iterator"
#include "limits"
#include "map"
#include "set"
#include "string"
#include "utility"
#include "vector"

#include "dccrg.hpp"
#include "rapidjson/document.h"

#include "grid/amr.hpp"
#include "mhd/common.hpp"
#include "mhd/solve_staggered.hpp"
#include "mhd/variables.hpp"
#include "variables.hpp"


namespace pamhd {
namespace mhd {


std::tuple<
	double, double, // mass (kg), pressure (Pa)
	std::array<double, 3>, // velocity (m/s)
	std::array<double, 3> // magnetic field (T)
> get_solar_wind_parameters(
	const rapidjson::Document& json,
	const double& /*sim_time*/,
	const double& proton_mass
) {
	using std::string;
	using std::invalid_argument;
	using std::to_string;

	if (not json.HasMember("number-density")) {
		throw invalid_argument(
			string(__FILE__ "(") + to_string(__LINE__) + "): "
			+ "JSON data doesn't have a number-density key."
		);
	}
	const auto& mass_json = json["number-density"];
	if (not mass_json.IsNumber()) {
		throw invalid_argument(
			string(__FILE__ "(") + to_string(__LINE__) + "): "
			+ "JSON item number-density is not a number."
		);
	}
	const double mass = proton_mass * json["number-density"].GetDouble();

	if (not json.HasMember("pressure")) {
		throw invalid_argument(
			string(__FILE__ "(") + to_string(__LINE__) + "): "
			+ "JSON data doesn't have a pressure key."
		);
	}
	const auto& pressure_json = json["pressure"];
	if (not pressure_json.IsNumber()) {
		throw invalid_argument(
			string(__FILE__ "(") + to_string(__LINE__) + "): "
			+ "JSON item pressure is not a number."
		);
	}
	const double pressure = json["pressure"].GetDouble();

	if (not json.HasMember("velocity")) {
		throw invalid_argument(
			string(__FILE__ "(") + to_string(__LINE__) + "): "
			+ "JSON data doesn't have a velocity key."
		);
	}
	const auto& velocity_json = json["velocity"];
	if (not velocity_json.IsArray()) {
		throw invalid_argument(
			string(__FILE__ "(") + to_string(__LINE__) + "): "
			+ "JSON item velocity is not an array."
		);
	}
	const auto& velocity_array = velocity_json.GetArray();
	if (velocity_array.Size() != 3) {
		throw invalid_argument(
			string(__FILE__ "(") + to_string(__LINE__) + "): "
			+ "Invalid size for velocity, should be 3"
		);
	}
	const std::array<double, 3> velocity{
		velocity_array[0].GetDouble(),
		velocity_array[1].GetDouble(),
		velocity_array[2].GetDouble()
	};

	if (not json.HasMember("magnetic-field")) {
		throw invalid_argument(
			string(__FILE__ "(") + to_string(__LINE__) + "): "
			+ "JSON data doesn't have a magnetic field key."
		);
	}
	const auto& magnetic_field_json = json["magnetic-field"];
	if (not magnetic_field_json.IsArray()) {
		throw invalid_argument(
			string(__FILE__ "(") + to_string(__LINE__) + "): "
			+ "JSON item magnetic field is not an array."
		);
	}
	const auto& magnetic_field_array = magnetic_field_json.GetArray();
	if (magnetic_field_array.Size() != 3) {
		throw invalid_argument(
			string(__FILE__ "(") + to_string(__LINE__) + "): "
			+ "Invalid size for magnetic field, should be 3"
		);
	}
	const std::array<double, 3> magnetic_field{
		magnetic_field_array[0].GetDouble(),
		magnetic_field_array[1].GetDouble(),
		magnetic_field_array[2].GetDouble()
	};

	return std::make_tuple(mass, pressure, velocity, magnetic_field);
}


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
	const double& sim_time,
	const rapidjson::Document& json,
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
	using std::to_string;

	const auto [
		mass, pressure, velocity, magnetic_field
	] = get_solar_wind_parameters(
		json, sim_time, proton_mass
	);
	if (grid.get_rank() == 0) {
		std::cout << "Initializing run, solar wind: "
		<< mass/proton_mass << " #/m^3, "
		<< velocity << " m/s, "
		<< pressure << " Pa, "
		<< magnetic_field << " T"
		<< std::endl;
	}

	for (const auto& cell: grid.local_cells()) {
		const auto [rx, ry, rz] = grid.geometry.get_center(cell.id);
		const auto [sx, sy, sz] = grid.geometry.get_min(cell.id);
		const auto [ex, ey, ez] = grid.geometry.get_max(cell.id);

		Bg_B.data(*cell.data)(-1) = background_B.get_background_field(
			{sx, ry, rz},
			vacuum_permeability
		);
		Bg_B.data(*cell.data)(+1) = background_B.get_background_field(
			{ex, ry, rz},
			vacuum_permeability
		);
		Bg_B.data(*cell.data)(-2) = background_B.get_background_field(
			{rx, sy, rz},
			vacuum_permeability
		);
		Bg_B.data(*cell.data)(+2) = background_B.get_background_field(
			{rx, ey, rz},
			vacuum_permeability
		);
		Bg_B.data(*cell.data)(-3) = background_B.get_background_field(
			{rx, ry, sz},
			vacuum_permeability
		);
		Bg_B.data(*cell.data)(+3) = background_B.get_background_field(
			{rx, ry, ez},
			vacuum_permeability
		);

		Mas.data(*cell.data) = mass;
		Mom.data(*cell.data)   =
		Vol_B.data(*cell.data) = {0, 0, 0};
		Face_B.data(*cell.data)  =
		Face_dB.data(*cell.data) = {0, 0, 0, 0, 0, 0};
		Nrj.data(*cell.data) = pamhd::mhd::get_total_energy_density(
			Mas.data(*cell.data),
			Mom.data(*cell.data),
			pressure,
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
	class JSON,
	class Cells,
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Energy_Density_Getter,
	class Volume_Magnetic_Field_Getter,
	class Face_Magnetic_Field_Getter
> void apply_solar_wind_boundaries(
	const JSON& json,
	const Cells& solar_wind_cells,
	const int& solar_wind_dir,
	const double& sim_time,
	const double& adiabatic_index,
	const double& vacuum_permeability,
	const double& proton_mass,
	const Mass_Density_Getter& Mas,
	const Momentum_Density_Getter& Mom,
	const Energy_Density_Getter& Nrj,
	const Volume_Magnetic_Field_Getter& Vol_B,
	const Face_Magnetic_Field_Getter& Face_B
) try {
	const auto [
		mass, pressure, velocity, magnetic_field
	] = pamhd::mhd::get_solar_wind_parameters(
		json, sim_time, proton_mass
	);
	for (const auto& cell: solar_wind_cells) {
		Mas.data(*cell.data) = mass;
		Mom.data(*cell.data) = {
			mass*velocity[0],
			mass*velocity[1],
			mass*velocity[2]
		};
		// prevent div(B) at solar wind boundary
		if (solar_wind_dir != abs(1)) {
			Vol_B.data(*cell.data)[0]   =
			Face_B.data(*cell.data)(-1) =
			Face_B.data(*cell.data)(+1) = magnetic_field[0];
		}
		if (solar_wind_dir != abs(2)) {
			Vol_B.data(*cell.data)[1]   =
			Face_B.data(*cell.data)(-2) =
			Face_B.data(*cell.data)(+2) = magnetic_field[1];
		}
		if (solar_wind_dir != abs(3)) {
			Vol_B.data(*cell.data)[2]   =
			Face_B.data(*cell.data)(-3) =
			Face_B.data(*cell.data)(+3) = magnetic_field[2];
		}
		Nrj.data(*cell.data) = pamhd::mhd::get_total_energy_density(
			mass, velocity, pressure, magnetic_field,
			adiabatic_index, vacuum_permeability
		);
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
	class JSON,
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
	class Solver_Info_Getter
> void apply_boundaries(
	Grid& grid,
	const double& sim_time,
	const JSON& json,
	const int& solar_wind_dir,
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
	const Solver_Info_Getter& SInfo
) try {
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
		// everything will become stale again below
	}
	if (Face_B.type().is_stale) {
		update_copies = true;
		Cell::set_transfer_all(true, Face_B.type());
	}
	if (update_copies) {
		grid.update_copies_of_remote_neighbors();
	}
	Cell::set_transfer_all(false,
		Mas.type(), Mom.type(), Nrj.type(),
		Vol_B.type(), Face_B.type());
	// boundary and normal cell share face
	for (int dir: {-3,-2,-1,+1,+2,+3}) {
		for (const auto& cell: face_cells(dir)) {
			for (const auto& neighbor: cell.neighbors_of) {
				const auto& fn = neighbor.face_neighbor;
				if (fn != -dir) continue;
				Mas.data(*cell.data) = Mas.data(*neighbor.data);
				Mom.data(*cell.data) = Mom.data(*neighbor.data);
				Nrj.data(*cell.data) = Nrj.data(*neighbor.data);
				Vol_B.data(*cell.data) = Vol_B.data(*neighbor.data);
				for (int dir2: {-3,-2,-1,+1,+2,+3}) {
					if (dir2 == fn or dir2 == -fn) {
						Face_B.data(*cell.data)(dir2) = Face_B.data(*neighbor.data)(-fn);
					} else {
						Face_B.data(*cell.data)(dir2) = Face_B.data(*neighbor.data)(dir2);
					}
				}
			}
		}
	}

	// boundary and normal cell share edge
	for (int dim: {0, 1, 2})
	for (int dir1: {-1, +1})
	for (int dir2: {-1, +1}) {
		for (const auto& cell: edge_cells(dim, dir1, dir2)) {
			for (const auto& neighbor: cell.neighbors_of) {
				const auto& en = neighbor.edge_neighbor;
				if (en[0] != dim or en[1] != -dir1 or en[2] != -dir2) continue;
				Mas.data(*cell.data) = Mas.data(*neighbor.data);
				Mom.data(*cell.data) = Mom.data(*neighbor.data);
				Nrj.data(*cell.data) = Nrj.data(*neighbor.data);
				Vol_B.data(*cell.data) = Vol_B.data(*neighbor.data);
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
		}
	}

	// boundary and normal cell share vertex
	for (const auto& cell: vert_cells) {
		for (const auto& neighbor: cell.neighbors_of) {
			const auto& fn = neighbor.face_neighbor;
			if (fn != 0) continue;
			const auto& en = neighbor.edge_neighbor;
			if (en[0] >= 0) continue;

			Mas.data(*cell.data) = Mas.data(*neighbor.data);
			Mom.data(*cell.data) = Mom.data(*neighbor.data);
			Nrj.data(*cell.data) = Nrj.data(*neighbor.data);
			Vol_B.data(*cell.data) = Vol_B.data(*neighbor.data);
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
	}

	// planetary boundary cells
	for (const auto& cell: planet_cells) {
		Mas.data(*cell.data) = proton_mass * 1e9;
		Mom.data(*cell.data) = {0, 0, 0};
		Vol_B.data(*cell.data) = {0, 0, 0};
		Nrj.data(*cell.data) = pamhd::mhd::get_total_energy_density(
			Mas.data(*cell.data),
			Mom.data(*cell.data),
			1e-11,
			Vol_B.data(*cell.data),
			adiabatic_index,
			vacuum_permeability
		);
		Face_B.data(*cell.data)(-1) =
		Face_B.data(*cell.data)(+1) =
		Face_B.data(*cell.data)(-2) =
		Face_B.data(*cell.data)(+2) =
		Face_B.data(*cell.data)(-3) =
		Face_B.data(*cell.data)(+3) = 0;

		pamhd::grid::Face_Type<bool> have_value{false, false, false, false, false, false};
		// corrections to Face_B from normal neighbors
		for (const auto& neighbor: cell.neighbors_of) {
			if (SInfo.data(*neighbor.data) < 1) continue;

			const auto& fn = neighbor.face_neighbor;
			const auto& en = neighbor.edge_neighbor;

			if (fn != 0) {
				Face_B.data(*cell.data)(fn) = Face_B.data(*neighbor.data)(-fn);
				have_value(fn) = true;
			} else if (en[0] >= 0) {
				// TODO
			} else {
				// TODO
			}
		}
		for (int dir: {-3,-2,-1,+1,+2,+3}) {
			if (not have_value(dir)) {
				if (have_value(-dir)) {
					Face_B.data(*cell.data)(dir) = Face_B.data(*cell.data)(-dir);
				} else {
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
			Mom.data(*cell.data),
			1e-11,
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

	apply_solar_wind_boundaries(
		json, solar_wind_cells, solar_wind_dir, sim_time,
		adiabatic_index, vacuum_permeability, proton_mass,
		Mas, Mom, Nrj, Vol_B, Face_B
	);

} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


}} // namespaces


#endif // ifndef PAMHD_MHD_SW_BOX_HPP
