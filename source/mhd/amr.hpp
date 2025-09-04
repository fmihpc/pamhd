/*
Adaptive mesh refinement logic of MHD part of PAMHD

Copyright 2023, 2024, 2025 Finnish Meteorological Institute
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

#ifndef PAMHD_MHD_AMR_HPP
#define PAMHD_MHD_AMR_HPP


#include "array"
#include "cmath"
#include "limits"
#include "set"
#include "string"
#include "tuple"
#include "vector"

#include "dccrg.hpp"
#include "gensimcell.hpp"
#include "prettyprint.hpp"

#include "grid/amr.hpp"
#include "grid/options.hpp"
#include "grid/variables.hpp"
#include "mhd/boundaries.hpp"
#include "mhd/common.hpp"
#include "mhd/options.hpp"
#include "mhd/solve.hpp"
#include "mhd/variables.hpp"
#include "variables.hpp"


namespace pamhd {
namespace mhd {


/*! Constrains min,max refinement levels based on physics.

Dont_solve, boundary cells and their normal face neighbors
are kept at current refinement level.
*/
template <
	class Cells,
	class Grid,
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Volume_Magnetic_Field_Getter,
	class Solver_Info_Getter,
	class Target_Refinement_Level_Min_Getter,
	class Target_Refinement_Level_Max_Getter
> void set_minmax_refinement_level(
	const Cells& cells,
	Grid& grid,
	const Options& options_mhd,
	const Mass_Density_Getter& Mas,
	const Momentum_Density_Getter& Mom,
	const Total_Energy_Density_Getter& Nrj,
	const Volume_Magnetic_Field_Getter& Vol_B,
	const Solver_Info_Getter& SInfo,
	const Target_Refinement_Level_Min_Getter& RLMin,
	const Target_Refinement_Level_Max_Getter& RLMax,
	const double& adiabatic_index,
	const double& vacuum_permeability,
	const double& proton_mass
) try {
	using std::abs;
	using std::clamp;
	using std::max;
	using std::min;
	using std::round;
	using std::runtime_error;
	using std::to_string;

	const auto max_ref_lvl = grid.get_maximum_refinement_level();
	if (max_ref_lvl == 0) return;

	for (const auto& cell: cells) {

		bool skip = false;
		if (SInfo.data(*cell.data) <= 0) {
			skip = true;
		} else {
			for (const auto& neighbor: cell.neighbors_of) {
				if (SInfo.data(*neighbor.data) <= 0) {
					skip = true;
					break;
				}
			}
		}
		if (skip) {
			RLMin.data(*cell.data) =
			RLMax.data(*cell.data) = grid.get_refinement_level(cell.id);
			continue;
		}

		if (Mas.data(*cell.data) <= 0) {
			throw runtime_error(__FILE__ ":" + to_string(__LINE__));
		}
		if (Nrj.data(*cell.data) <= 0) {
			throw runtime_error(__FILE__ ":" + to_string(__LINE__));
		}

		const auto [cdx, cdy, cdz] = grid.geometry.get_length(cell.id);
		const auto
			cpre = max(options_mhd.pressure_min_mrg,
				get_pressure(
					Mas.data(*cell.data), Mom.data(*cell.data),
					Nrj.data(*cell.data), Vol_B.data(*cell.data),
					adiabatic_index, vacuum_permeability)),
			cmas = max(options_mhd.number_density_min_mrg,
				Mas.data(*cell.data) / proton_mass),
			cvel = max(options_mhd.vel_min_mrg,
				pamhd::norm(Mom.data(*cell.data)) / Mas.data(*cell.data)),
			cmag = max(options_mhd.mag_min_mrg,
				pamhd::norm(Vol_B.data(*cell.data)));

		// maximum relative gradients w.r.t. neighbors
		double
			mrg_mas = 0,
			mrg_vel = 0,
			mrg_pre = 0,
			mrg_mag = 0;
		for (const auto& neighbor: cell.neighbors_of) {
			if (SInfo.data(*neighbor.data) <= 0) {
				continue;
			}

			if (neighbor.face_neighbor == 0) {
				continue;
			}

			const auto [ndx, ndy, ndz] = grid.geometry.get_length(neighbor.id);
			const auto
				npre = max(options_mhd.pressure_min_mrg,
					get_pressure(
						Mas.data(*neighbor.data), Mom.data(*neighbor.data),
						Nrj.data(*neighbor.data), Vol_B.data(*neighbor.data),
						adiabatic_index, vacuum_permeability)),
				nmas = max(options_mhd.number_density_min_mrg,
					Mas.data(*neighbor.data) / proton_mass),
				nvel = max(options_mhd.vel_min_mrg,
					pamhd::norm(Mom.data(*neighbor.data)) / Mas.data(*neighbor.data)),
				nmag = max(options_mhd.mag_min_mrg,
					pamhd::norm(Vol_B.data(*neighbor.data)));

			const auto len = [&](){
				switch (neighbor.face_neighbor) {
				case -1:
				case +1:
					return 0.5 * (cdx + ndx);
					break;
				case -2:
				case +2:
					return 0.5 * (cdy + ndy);
					break;
				case -3:
				case +3:
					return 0.5 * (cdz + ndz);
					break;
				default:
					throw runtime_error(__FILE__ ":" + to_string(__LINE__));
					break;
				}
				return -1.0;
			}();

			mrg_mas = max(mrg_mas, abs(cmas-nmas) / max(cmas, nmas) / len);
			mrg_vel = max(mrg_vel, abs(cvel-nvel) / max(cvel, nvel) / len);
			mrg_pre = max(mrg_pre, abs(cpre-npre) / max(cpre, npre) / len);
			mrg_mag = max(mrg_mag, abs(cmag-nmag) / max(cmag, nmag) / len);
		}

		const double tgt_ref_lvl = max_ref_lvl * max(max(max(
			mrg_mas / options_mhd.number_density_mrl_at,
			mrg_vel / options_mhd.vel_mrl_at),
			mrg_pre / options_mhd.pressure_mrl_at),
			mrg_mag / options_mhd.mag_mrl_at);

		const auto
			prev_min = RLMin.data(*cell.data),
			prev_max = RLMax.data(*cell.data);
		RLMin.data(*cell.data) = clamp(
			max(int(round(tgt_ref_lvl-0.25)), RLMin.data(*cell.data)),
			prev_min, prev_max);
		RLMax.data(*cell.data) = clamp(
			min(int(round(tgt_ref_lvl+0.25)), RLMax.data(*cell.data)),
			prev_min, prev_max);
	}
} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


// Handlers that handle all parameters required by MHD test program.
template <
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Face_Magnetic_Field_Getter,
	class Volume_Magnetic_Field_Getter,
	class Background_Magnetic_Field_Getter,
	class Background_Magnetic_Field,
	class Target_Refinement_Level_Min_Getter,
	class Target_Refinement_Level_Max_Getter,
	class Maximum_Signal_Velocity_Getter,
	class Face_B_Error_Getter
> struct New_Cells_Handler {
	const Mass_Density_Getter& Mas;
	const Momentum_Density_Getter& Mom;
	const Total_Energy_Density_Getter& Nrj;
	const Face_Magnetic_Field_Getter& Face_B;
	const Volume_Magnetic_Field_Getter& Vol_B;
	const Background_Magnetic_Field_Getter& Bg_B;
	const Background_Magnetic_Field& bg_B;
	const Target_Refinement_Level_Min_Getter& RLMin;
	const Target_Refinement_Level_Max_Getter& RLMax;
	const Maximum_Signal_Velocity_Getter& Max_v;
	const Face_B_Error_Getter& Berror;
	const double& adiabatic_index;
	const double& vacuum_permeability;

	New_Cells_Handler(
		const Mass_Density_Getter& Mas_,
		const Momentum_Density_Getter& Mom_,
		const Total_Energy_Density_Getter& Nrj_,
		const Face_Magnetic_Field_Getter& Face_B_,
		const Volume_Magnetic_Field_Getter& Vol_B_,
		const Background_Magnetic_Field_Getter& Bg_B_,
		const Background_Magnetic_Field& bg_B_,
		const Target_Refinement_Level_Min_Getter& RLMin_,
		const Target_Refinement_Level_Max_Getter& RLMax_,
		const Maximum_Signal_Velocity_Getter& Max_v_,
		const Face_B_Error_Getter& Berror_,
		const double& adiabatic_index_,
		const double& vacuum_permeability_
	) :
		Mas(Mas_), Mom(Mom_), Nrj(Nrj_),
		Face_B(Face_B_), Vol_B(Vol_B_), Bg_B(Bg_B_), bg_B(bg_B_),
		RLMin(RLMin_), RLMax(RLMax_), Max_v(Max_v_),
		Berror(Berror_), adiabatic_index(adiabatic_index_),
		vacuum_permeability(vacuum_permeability_)
	{};

	template<
		class Grid,
		class Cells
	> void operator()(
		const Grid& grid,
		const Cells& new_cells
	) const {
		using std::array;
		using std::runtime_error;
		using std::to_string;

		for (const auto& new_cell_id: new_cells) {
			auto* const cell_data = grid[new_cell_id];
			if (cell_data == nullptr) {
				throw runtime_error(__FILE__ ":" + to_string(__LINE__));
			}

			const auto parent_id = grid.mapping.get_parent(new_cell_id);
			auto* const parent_data = grid[parent_id];
			if (parent_data == nullptr) {
				throw runtime_error(__FILE__ ":" + to_string(__LINE__));
			}
			Mas.data(*cell_data) = Mas.data(*parent_data);
			Mom.data(*cell_data) = Mom.data(*parent_data);
			RLMin.data(*cell_data) = RLMin.data(*parent_data);
			RLMax.data(*cell_data) = RLMax.data(*parent_data);
			Max_v.data(*cell_data) = Max_v.data(*parent_data);
			Berror.data(*cell_data) = Berror.data(*parent_data);

			// inherit thermal pressure
			const double parent_pressure = [&](){
				try {
					return get_pressure(
						Mas.data(*parent_data),
						Mom.data(*parent_data),
						Nrj.data(*parent_data),
						Vol_B.data(*parent_data),
						adiabatic_index,
						vacuum_permeability
					);
				} catch (...) {
					throw runtime_error(__FILE__ ":" + to_string(__LINE__));
				}
			}();

			const auto
				cindex = grid.mapping.get_indices(new_cell_id),
				pindex = grid.mapping.get_indices(parent_id);
			for (auto dim: {0, 1, 2}) {
				const int side = [&](){
					if (cindex[dim] == pindex[dim]) {
						// neg faces coincide
						return -1;
					} else {
						// pos faces coincide
						return +1;
					}
				}();
				Face_B.data(*cell_data)(dim, side) = Face_B.data(*parent_data)(dim, side);
				Face_B.data(*cell_data)(dim, -side) = Vol_B.data(*parent_data)[dim];
			}

			for (auto dim: {0, 1, 2}) {
				Vol_B.data(*cell_data)[dim] = 0.5 * (Face_B.data(*cell_data)(dim, -1) + Face_B.data(*cell_data)(dim, +1));
			}
			Nrj.data(*cell_data) = get_total_energy_density(
				Mas.data(*cell_data),
				pamhd::mul(Mom.data(*cell_data), 1 / Mas.data(*cell_data)),
				parent_pressure,
				Vol_B.data(*cell_data),
				adiabatic_index,
				vacuum_permeability
			);

			const auto [
				center, start, end
			] = pamhd::grid::get_cell_geom_emulated(grid, new_cell_id);
			const auto [rx, ry, rz] = center;
			const auto [sx, sy, sz] = start;
			const auto [ex, ey, ez] = end;
			Bg_B.data(*cell_data)(0, -1) = bg_B.get_background_field(
				{sx, ry, rz},
				vacuum_permeability
			);
			Bg_B.data(*cell_data)(0, +1) = bg_B.get_background_field(
				{ex, ry, rz},
				vacuum_permeability
			);
			Bg_B.data(*cell_data)(1, -1) = bg_B.get_background_field(
				{rx, sy, rz},
				vacuum_permeability
			);
			Bg_B.data(*cell_data)(1, +1) = bg_B.get_background_field(
				{rx, ey, rz},
				vacuum_permeability
			);
			Bg_B.data(*cell_data)(2, -1) = bg_B.get_background_field(
				{rx, ry, sz},
				vacuum_permeability
			);
			Bg_B.data(*cell_data)(2, +1) = bg_B.get_background_field(
				{rx, ry, ez},
				vacuum_permeability
			);
		}
	}
};


template <
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Face_Magnetic_Field_Getter,
	class Volume_Magnetic_Field_Getter,
	class Background_Magnetic_Field_Getter,
	class Background_Magnetic_Field,
	class Target_Refinement_Level_Min_Getter,
	class Target_Refinement_Level_Max_Getter,
	class Maximum_Signal_Velocity_Getter
> struct Removed_Cells_Handler {
	const Mass_Density_Getter& Mas;
	const Momentum_Density_Getter& Mom;
	const Total_Energy_Density_Getter& Nrj;
	const Face_Magnetic_Field_Getter& Face_B;
	const Volume_Magnetic_Field_Getter& Vol_B;
	const Background_Magnetic_Field_Getter& Bg_B;
	const Background_Magnetic_Field& bg_B;
	const Target_Refinement_Level_Min_Getter& RLMin;
	const Target_Refinement_Level_Max_Getter& RLMax;
	const Maximum_Signal_Velocity_Getter& Max_v;
	const double& adiabatic_index;
	const double& vacuum_permeability;

	Removed_Cells_Handler(
		const Mass_Density_Getter& Mas_,
		const Momentum_Density_Getter& Mom_,
		const Total_Energy_Density_Getter& Nrj_,
		const Face_Magnetic_Field_Getter& Face_B_,
		const Volume_Magnetic_Field_Getter& Vol_B_,
		const Background_Magnetic_Field_Getter& Bg_B_,
		const Background_Magnetic_Field& bg_B_,
		const Target_Refinement_Level_Min_Getter& RLMin_,
		const Target_Refinement_Level_Max_Getter& RLMax_,
		const Maximum_Signal_Velocity_Getter& Max_v_,
		const double& adiabatic_index_,
		const double& vacuum_permeability_
	) :
		Mas(Mas_), Mom(Mom_), Nrj(Nrj_),
		Face_B(Face_B_), Vol_B(Vol_B_), Bg_B(Bg_B_), bg_B(bg_B_),
		RLMin(RLMin_), RLMax(RLMax_), Max_v(Max_v_),
		adiabatic_index(adiabatic_index_),
		vacuum_permeability(vacuum_permeability_)
	{};

	template<
		class Grid,
		class Cells
	> void operator()(
		const Grid& grid,
		const Cells& removed_cells
	) const {
		using std::max;
		using std::min;
		using std::runtime_error;
		using std::to_string;

		// process each parent of removed cells only once
		std::set<uint64_t> parents;
		for (const auto& removed_cell: removed_cells) {
			parents.insert(grid.mapping.get_parent(removed_cell));
		}

		// initialize data of parents of removed cells
		for (const auto& parent: parents) {
			auto* const parent_data = grid[parent];
			if (parent_data == nullptr) {
				throw runtime_error(__FILE__ ":" + to_string(__LINE__));
			}
			RLMin.data(*parent_data) = 999;
			RLMax.data(*parent_data) = 0;
			Max_v.data(*parent_data) = {-1, -1, -1, -1, -1, -1};
			Mas.data(*parent_data) =
			Nrj.data(*parent_data) = 0;
			Mom.data(*parent_data) = {0, 0, 0};
			Face_B.data(*parent_data) = {0, 0, 0, 0, 0, 0};
		}

		// average parents' plasma parameters from their children, etc
		for (const auto& removed_cell_id: removed_cells) {
			auto* const removed_cell_data = grid[removed_cell_id];
			if (removed_cell_data == nullptr) {
				throw runtime_error(__FILE__ ":" + to_string(__LINE__));
			}

			const auto parent_id = grid.mapping.get_parent(removed_cell_id);
			auto* const parent_data = grid[parent_id];
			if (parent_data == nullptr) {
				throw runtime_error(__FILE__ ":" + to_string(__LINE__));
			}

			if (Mas.data(*removed_cell_data) <= 0) {
				throw runtime_error(
					__FILE__ "(" + to_string(__LINE__)
					+ "): " + to_string(removed_cell_id) + ", "
					+ to_string(parent_id)
				);
			}

			// plasma, all children included
			RLMin.data(*parent_data) = min(RLMin.data(*parent_data),
				RLMin.data(*removed_cell_data));
			RLMax.data(*parent_data) = max(RLMax.data(*parent_data),
				RLMax.data(*removed_cell_data));
			for (int dir: {-3,-2,-1,+1,+2,+3}) {
				Max_v.data(*parent_data)(dir) = max(Max_v.data(*parent_data)(dir),
					Max_v.data(*removed_cell_data)(dir));
			}
			Mas.data(*parent_data) += Mas.data(*removed_cell_data) / 8;
			Mom.data(*parent_data) = pamhd::add(Mom.data(*parent_data),
				pamhd::mul(Mom.data(*removed_cell_data), 1.0 / 8.0));
			// temporarily store pressure in energy density variable
			try {
				Nrj.data(*parent_data) += get_pressure(
					Mas.data(*removed_cell_data),
					Mom.data(*removed_cell_data),
					Nrj.data(*removed_cell_data),
					Vol_B.data(*removed_cell_data),
					adiabatic_index,
					vacuum_permeability
				) / 8;
			} catch (...) {
				throw runtime_error(__FILE__ ":" + to_string(__LINE__));
			}

			// magnetic field, only children sharing
			// face(s) with parent included
			const auto
				pindex = grid.mapping.get_indices(parent_id),
				rindex = grid.mapping.get_indices(removed_cell_id);
			for (auto dim: {0, 1, 2}) {
				const int side = [&](){
					if (pindex[dim] == rindex[dim]) {
						return -1;
					} else {
						return +1;
					}
				}();
				Face_B.data(*parent_data)(dim, side) += Face_B.data(*removed_cell_data)(dim, side) / 4;
			}

			const auto [
				center, start, end
			] = pamhd::grid::get_cell_geom_emulated(grid, parent_id);
			const auto [rx, ry, rz] = center;
			const auto [sx, sy, sz] = start;
			const auto [ex, ey, ez] = end;
			Bg_B.data(*parent_data)(0, -1) = bg_B.get_background_field(
				{sx, ry, rz},
				vacuum_permeability
			);
			Bg_B.data(*parent_data)(0, +1) = bg_B.get_background_field(
				{ex, ry, rz},
				vacuum_permeability
			);
			Bg_B.data(*parent_data)(1, -1) = bg_B.get_background_field(
				{rx, sy, rz},
				vacuum_permeability
			);
			Bg_B.data(*parent_data)(1, +1) = bg_B.get_background_field(
				{rx, ey, rz},
				vacuum_permeability
			);
			Bg_B.data(*parent_data)(2, -1) = bg_B.get_background_field(
				{rx, ry, sz},
				vacuum_permeability
			);
			Bg_B.data(*parent_data)(2, +1) = bg_B.get_background_field(
				{rx, ry, ez},
				vacuum_permeability
			);
		}
		// set volume magnetic field and energy density
		for (const auto& parent: parents) {
			auto* const parent_data = grid[parent];
			if (parent_data == nullptr) {
				throw runtime_error(__FILE__ ":" + to_string(__LINE__));
			}

			for (auto dim: {0, 1, 2}) {
				Vol_B.data(*parent_data)[dim] = 0.5 * (Face_B.data(*parent_data)(dim, -1) + Face_B.data(*parent_data)(dim, +1));
			}

			Nrj.data(*parent_data) = get_total_energy_density(
				Mas.data(*parent_data),
				pamhd::mul(Mom.data(*parent_data), 1 / Mas.data(*parent_data)),
				Nrj.data(*parent_data),
				Vol_B.data(*parent_data),
				adiabatic_index,
				vacuum_permeability
			);
		}
	}
};


template<
	class Grid,
	class Geometries,
	class Boundaries,
	class Background_B,
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Volume_Magnetic_Field_Getter,
	class Face_Magnetic_Field_Getter,
	class Background_Magnetic_Field_Getter,
	class Solver_Info_Getter,
	class Face_Info_Getter,
	class Target_Refinement_Level_Min_Getter,
	class Target_Refinement_Level_Max_Getter,
	class Substepping_Period_Getter,
	class Maximum_Signal_Speed_Getter,
	class Face_B_Error_Getter
> void adapt_grid(
	Grid& grid,
	pamhd::grid::Options& options_grid,
	const pamhd::mhd::Options& options_mhd,
	Geometries& geometries,
	Boundaries& boundaries,
	const Background_B& bg_B,
	const double& simulation_time,
	const double& proton_mass,
	const double& adiabatic_index,
	const double& vacuum_permeability,
	const Mass_Density_Getter& Mas,
	const Momentum_Density_Getter& Mom,
	const Total_Energy_Density_Getter& Nrj,
	const Volume_Magnetic_Field_Getter& Vol_B,
	const Face_Magnetic_Field_Getter& Face_B,
	const Background_Magnetic_Field_Getter& Bg_B,
	const Solver_Info_Getter& SInfo,
	const Face_Info_Getter& FInfo,
	const Target_Refinement_Level_Min_Getter& Ref_min,
	const Target_Refinement_Level_Max_Getter& Ref_max,
	const Substepping_Period_Getter& Substep,
	const Maximum_Signal_Speed_Getter& Max_v,
	const Face_B_Error_Getter& Berror
) try {
	using Cell = Grid::cell_data_type;

	pamhd::grid::set_minmax_refinement_level(
		grid.local_cells(), grid, options_grid,
		simulation_time, Ref_min, Ref_max, false);

	pamhd::mhd::set_minmax_refinement_level(
		grid.local_cells(), grid, options_mhd,
		Mas, Mom, Nrj, Vol_B, SInfo, Ref_min, Ref_max,
		adiabatic_index, vacuum_permeability, proton_mass);

	Cell::set_transfer_all(true,
		Mas.type(), Mom.type(), Nrj.type(),
		Vol_B.type(), Face_B.type());
	pamhd::grid::adapt_grid(
		grid, Ref_min, Ref_max,
		pamhd::mhd::New_Cells_Handler(
			Mas, Mom, Nrj, Face_B, Vol_B, Bg_B,
			bg_B, Ref_min, Ref_max, Max_v, Berror,
			adiabatic_index, vacuum_permeability),
		pamhd::mhd::Removed_Cells_Handler(
			Mas, Mom, Nrj, Face_B, Vol_B, Bg_B,
			bg_B, Ref_min, Ref_max, Max_v,
			adiabatic_index, vacuum_permeability)
	);
	Cell::set_transfer_all(true, Bg_B.type());
	grid.update_copies_of_remote_neighbors();
	Cell::set_transfer_all(false,
		Mas.type(), Mom.type(), Nrj.type(),
		Vol_B.type(), Face_B.type(), Bg_B.type());
	Mas.type().is_stale = false;
	Mom.type().is_stale = false;
	Nrj.type().is_stale = false;
	Vol_B.type().is_stale = false;
	Face_B.type().is_stale = false;
	Bg_B.type().is_stale = false;

	for (const auto& cell: grid.local_cells()) {
		(*cell.data)[pamhd::MPI_Rank()] = grid.get_rank();
		Substep.data(*cell.data) = 1;
	}
	Cell::set_transfer_all(true, Substep.type());
	grid.update_copies_of_remote_neighbors();
	Cell::set_transfer_all(false, Substep.type());
	Substep.type().is_stale = false;

	for (const auto& gid: geometries.get_geometry_ids()) {
		geometries.clear_cells(gid);
	}
	for (const auto& cell: grid.local_cells()) {
		geometries.overlaps(
			grid.geometry.get_min(cell.id),
			grid.geometry.get_max(cell.id),
			cell.id);
	}

	pamhd::mhd::set_solver_info(grid, boundaries, geometries, SInfo);
	pamhd::mhd::classify_faces(grid, SInfo, FInfo);

	sync_magnetic_field(grid, Face_B, Vol_B, Berror, SInfo);

	pamhd::mhd::update_vol_B(
		0, grid.local_cells(), Mas, Mom, Nrj, Vol_B, Face_B,
		SInfo, Substep, adiabatic_index, vacuum_permeability, true
	);
} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


}} // namespaces


#endif // ifndef PAMHD_MHD_AMR_HPP
