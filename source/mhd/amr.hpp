/*
Adaptive mesh refinement logic of MHD part of PAMHD

Copyright 2023, 2024 Finnish Meteorological Institute
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
#include "mhd/solve_staggered.hpp"
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
	const Volume_Magnetic_Field_Getter& Mag,
	const Solver_Info_Getter& Solver_Info,
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
		if (Solver_Info(*cell.data) <= 0) {
			skip = true;
		} else {
			for (const auto& neighbor: cell.neighbors_of) {
				if (Solver_Info(*neighbor.data) <= 0) {
					skip = true;
					break;
				}
			}
		}
		if (skip) {
			RLMin(*cell.data) =
			RLMax(*cell.data) = grid.get_refinement_level(cell.id);
			continue;
		}

		if (Mas(*cell.data) <= 0) {
			throw runtime_error(__FILE__ ":" + to_string(__LINE__));
		}
		if (Nrj(*cell.data) <= 0) {
			throw runtime_error(__FILE__ ":" + to_string(__LINE__));
		}

		const auto [cdx, cdy, cdz] = grid.geometry.get_length(cell.id);
		const auto
			cpre = max(options_mhd.pressure_min_mrg, get_pressure(
				Mas(*cell.data), Mom(*cell.data), Nrj(*cell.data),
				Mag(*cell.data), adiabatic_index, vacuum_permeability)),
			cmas = max(Mas(*cell.data) / proton_mass, options_mhd.number_density_min_mrg),
			cvel = max((Mom(*cell.data) / Mas(*cell.data)).norm(), options_mhd.vel_min_mrg),
			cmag = max(Mag(*cell.data).norm(), options_mhd.mag_min_mrg);

		// maximum relative gradients w.r.t. neighbors
		double
			mrg_mas = 0,
			mrg_vel = 0,
			mrg_pre = 0,
			mrg_mag = 0;
		for (const auto& neighbor: cell.neighbors_of) {
			if (Solver_Info(*neighbor.data) <= 0) {
				continue;
			}

			if (neighbor.face_neighbor == 0) {
				continue;
			}

			const auto [ndx, ndy, ndz] = grid.geometry.get_length(neighbor.id);
			const auto
				npre = max(options_mhd.pressure_min_mrg, get_pressure(
					Mas(*neighbor.data), Mom(*neighbor.data), Nrj(*neighbor.data),
					Mag(*neighbor.data), adiabatic_index, vacuum_permeability)),
				nmas = max(Mas(*neighbor.data) / proton_mass, options_mhd.number_density_min_mrg),
				nvel = max((Mom(*neighbor.data) / Mas(*neighbor.data)).norm(), options_mhd.vel_min_mrg),
				nmag = max(Mag(*neighbor.data).norm(), options_mhd.mag_min_mrg);

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
			prev_min = RLMin(*cell.data),
			prev_max = RLMax(*cell.data);
		RLMin(*cell.data) = clamp(
			max(int(round(tgt_ref_lvl-0.25)), RLMin(*cell.data)),
			prev_min, prev_max);
		RLMax(*cell.data) = clamp(
			min(int(round(tgt_ref_lvl+0.25)), RLMax(*cell.data)),
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
	class Maximum_Signal_Velocity_Getter
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
			Mas(*cell_data) = Mas(*parent_data);
			Mom(*cell_data) = Mom(*parent_data);
			RLMin(*cell_data) = RLMin(*parent_data);
			RLMax(*cell_data) = RLMax(*parent_data);
			Max_v(*cell_data) = Max_v(*parent_data);

			// inherit thermal pressure
			const double parent_pressure = [&](){
				try {
					return get_pressure(
						Mas(*parent_data),
						Mom(*parent_data),
						Nrj(*parent_data),
						Vol_B(*parent_data),
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
				Face_B(*cell_data)(dim, side) = Face_B(*parent_data)(dim, side);
				Face_B(*cell_data)(dim, -side) = Vol_B(*parent_data)[dim];
			}

			for (auto dim: {0, 1, 2}) {
				Vol_B(*cell_data)[dim] = 0.5 * (Face_B(*cell_data)(dim, -1) + Face_B(*cell_data)(dim, +1));
			}
			Nrj(*cell_data) = get_total_energy_density(
				Mas(*cell_data),
				(Mom(*cell_data) / Mas(*cell_data)).eval(),
				parent_pressure,
				Vol_B(*cell_data),
				adiabatic_index,
				vacuum_permeability
			);

			const auto [rx, ry, rz] = grid.geometry.get_center(new_cell_id);
			const auto [sx, sy, sz] = grid.geometry.get_min(new_cell_id);
			const auto [ex, ey, ez] = grid.geometry.get_max(new_cell_id);
			Bg_B(*cell_data)(0, -1) = bg_B.get_background_field(
				{sx, ry, rz},
				vacuum_permeability
			);
			Bg_B(*cell_data)(0, +1) = bg_B.get_background_field(
				{ex, ry, rz},
				vacuum_permeability
			);
			Bg_B(*cell_data)(1, -1) = bg_B.get_background_field(
				{rx, sy, rz},
				vacuum_permeability
			);
			Bg_B(*cell_data)(1, +1) = bg_B.get_background_field(
				{rx, ey, rz},
				vacuum_permeability
			);
			Bg_B(*cell_data)(2, -1) = bg_B.get_background_field(
				{rx, ry, sz},
				vacuum_permeability
			);
			Bg_B(*cell_data)(2, +1) = bg_B.get_background_field(
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
			RLMin(*parent_data) = 999;
			RLMax(*parent_data) = 0;
			Max_v(*parent_data) = {-1, -1, -1, -1, -1, -1};
			Mas(*parent_data) =
			Nrj(*parent_data) = 0;
			Mom(*parent_data) = {0, 0, 0};
			Face_B(*parent_data) = {0, 0, 0, 0, 0, 0};
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

			if (Mas(*removed_cell_data) <= 0) {
				throw runtime_error(
					__FILE__ "(" + to_string(__LINE__)
					+ "): " + to_string(removed_cell_id) + ", "
					+ to_string(parent_id)
				);
			}

			// all children included
			RLMin(*parent_data) = min(RLMin(*parent_data),
				RLMin(*removed_cell_data));
			RLMax(*parent_data) = max(RLMax(*parent_data),
				RLMax(*removed_cell_data));
			for (int dir: {-3,-2,-1,+1,+2,+3}) {
				Max_v(*parent_data)(dir) = max(Max_v(*parent_data)(dir),
					Max_v(*removed_cell_data)(dir));
			}
			Mas(*parent_data) += Mas(*removed_cell_data) / 8;
			Mom(*parent_data) += Mom(*removed_cell_data) / 8;
			// temporarily store pressure in energy density variable
			try {
				Nrj(*parent_data) += get_pressure(
					Mas(*removed_cell_data),
					Mom(*removed_cell_data),
					Nrj(*removed_cell_data),
					Vol_B(*removed_cell_data),
					adiabatic_index,
					vacuum_permeability
				) / 8;
			} catch (...) {
				throw runtime_error(__FILE__ ":" + to_string(__LINE__));
			}

			const auto
				pindex = grid.mapping.get_indices(parent_id),
				rindex = grid.mapping.get_indices(removed_cell_id);
			for (auto dim: {0, 1, 2}) {
				// only children sharing a face with parent included
				const int side = [&](){
					if (pindex[dim] == rindex[dim]) {
						return -1;
					} else {
						return +1;
					}
				}();
				Face_B(*parent_data)(dim, side) += Face_B(*removed_cell_data)(dim, side) / 4;
			}

			const auto [rx, ry, rz] = grid.geometry.get_center(parent_id);
			const auto [sx, sy, sz] = grid.geometry.get_min(parent_id);
			const auto [ex, ey, ez] = grid.geometry.get_max(parent_id);
			Bg_B(*parent_data)(0, -1) = bg_B.get_background_field(
				{sx, ry, rz},
				vacuum_permeability
			);
			Bg_B(*parent_data)(0, +1) = bg_B.get_background_field(
				{ex, ry, rz},
				vacuum_permeability
			);
			Bg_B(*parent_data)(1, -1) = bg_B.get_background_field(
				{rx, sy, rz},
				vacuum_permeability
			);
			Bg_B(*parent_data)(1, +1) = bg_B.get_background_field(
				{rx, ey, rz},
				vacuum_permeability
			);
			Bg_B(*parent_data)(2, -1) = bg_B.get_background_field(
				{rx, ry, sz},
				vacuum_permeability
			);
			Bg_B(*parent_data)(2, +1) = bg_B.get_background_field(
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
				Vol_B(*parent_data)[dim] = 0.5 * (Face_B(*parent_data)(dim, -1) + Face_B(*parent_data)(dim, +1));
			}

			Nrj(*parent_data) = get_total_energy_density(
				Mas(*parent_data),
				(Mom(*parent_data) / Mas(*parent_data)).eval(),
				Nrj(*parent_data),
				Vol_B(*parent_data),
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
	class Magnetic_Field_Getter,
	class Face_Magnetic_Field_Getter,
	class Background_Magnetic_Field_Getter,
	class Solver_Info_Getter,
	class Solver_Info_Getter2,
	class Face_Info_Getter,
	class Target_Refinement_Level_Min_Getter,
	class Target_Refinement_Level_Max_Getter,
	class Substepping_Period_Getter,
	class Maximum_Signal_Speed_Getter
> void adapt_grid(
	Grid& grid,
	pamhd::grid::Options& options_grid,
	const pamhd::mhd::Options& options_mhd,
	Geometries& geometries,
	Boundaries& boundaries,
	const Background_B& bg_B,
	const double simulation_time,
	const double proton_mass,
	const double adiabatic_index,
	const double vacuum_permeability,
	const Mass_Density_Getter Mas,
	const Momentum_Density_Getter Mom,
	const Total_Energy_Density_Getter Nrj,
	const Magnetic_Field_Getter Mag,
	const Face_Magnetic_Field_Getter Face_B,
	const Background_Magnetic_Field_Getter Bg_B_Getter,
	const Solver_Info_Getter SInfo,
	const Solver_Info_Getter2 SInfo2,
	const Face_Info_Getter FInfo,
	const Target_Refinement_Level_Min_Getter Ref_min,
	const Target_Refinement_Level_Max_Getter Ref_max,
	const Substepping_Period_Getter Substep,
	const Maximum_Signal_Speed_Getter Max_v
) try {
	pamhd::grid::set_minmax_refinement_level(
		grid.local_cells(), grid, options_grid,
		simulation_time, Ref_min, Ref_max, false);

	pamhd::mhd::set_minmax_refinement_level(
		grid.local_cells(), grid, options_mhd,
		Mas, Mom, Nrj, Mag, SInfo2, Ref_min, Ref_max,
		adiabatic_index, vacuum_permeability, proton_mass);

	Grid::cell_data_type::set_transfer_all(true,
		pamhd::Face_Magnetic_Field(),
		pamhd::mhd::MHD_State_Conservative(),
		pamhd::grid::Target_Refinement_Level_Min(),
		pamhd::grid::Target_Refinement_Level_Max()
	);
	pamhd::grid::adapt_grid(
		grid, Ref_min, Ref_max,
		pamhd::mhd::New_Cells_Handler(
			Mas, Mom, Nrj, Face_B, Mag, Bg_B_Getter,
			bg_B, Ref_min, Ref_max, Max_v,
			adiabatic_index, vacuum_permeability),
		pamhd::mhd::Removed_Cells_Handler(
			Mas, Mom, Nrj, Face_B, Mag, Bg_B_Getter,
			bg_B, Ref_min, Ref_max, Max_v,
			adiabatic_index, vacuum_permeability));
	Grid::cell_data_type::set_transfer_all(false,
		pamhd::grid::Target_Refinement_Level_Min(),
		pamhd::grid::Target_Refinement_Level_Max()
	);
	Grid::cell_data_type::set_transfer_all(true,
		pamhd::Bg_Magnetic_Field()
	);
	grid.update_copies_of_remote_neighbors();
	Grid::cell_data_type::set_transfer_all(false,
		pamhd::Face_Magnetic_Field(),
		pamhd::mhd::MHD_State_Conservative(),
		pamhd::Bg_Magnetic_Field()
	);

	for (const auto& cell: grid.local_cells()) {
		(*cell.data)[pamhd::MPI_Rank()] = grid.get_rank();
		Substep(*cell.data) = 1;
	}
	Grid::cell_data_type::set_transfer_all(true, pamhd::mhd::Substepping_Period());
	grid.update_copies_of_remote_neighbors();
	Grid::cell_data_type::set_transfer_all(false, pamhd::mhd::Substepping_Period());

	for (const auto& gid: geometries.get_geometry_ids()) {
		geometries.clear_cells(gid);
	}
	for (const auto& cell: grid.local_cells()) {
		geometries.overlaps(
			grid.geometry.get_min(cell.id),
			grid.geometry.get_max(cell.id),
			cell.id);
	}

	Grid::cell_data_type::set_transfer_all(true, pamhd::mhd::Solver_Info());
	pamhd::mhd::set_solver_info<pamhd::mhd::Solver_Info>(
		grid, boundaries, geometries, SInfo
	);
	grid.update_copies_of_remote_neighbors();
	Grid::cell_data_type::set_transfer_all(false, pamhd::mhd::Solver_Info());

	Grid::cell_data_type::set_transfer_all(true, pamhd::mhd::Face_Boundary_Type());
	pamhd::mhd::classify_faces(grid, SInfo2, FInfo);
	grid.update_copies_of_remote_neighbors();
	Grid::cell_data_type::set_transfer_all(false, pamhd::mhd::Face_Boundary_Type());

	pamhd::mhd::update_B_consistency(
		0, grid.local_cells(),
		Mas, Mom, Nrj, Mag, Face_B,
		SInfo, Substep,
		adiabatic_index,
		vacuum_permeability,
		true
	);
	Grid::cell_data_type::set_transfer_all(true,
		pamhd::Face_Magnetic_Field(),
		pamhd::mhd::MHD_State_Conservative()
	);
	grid.update_copies_of_remote_neighbors();
	Grid::cell_data_type::set_transfer_all(false,
		pamhd::Face_Magnetic_Field(),
		pamhd::mhd::MHD_State_Conservative()
	);
} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


}} // namespaces


#endif // ifndef PAMHD_MHD_AMR_HPP
