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

#include "grid/variables.hpp"
#include "mhd/common.hpp"
#include "mhd/options.hpp"
#include "mhd/variables.hpp"


namespace pamhd {
namespace mhd {


/*! Calculates target refinement levels for given cells.

Target is based on gradients of plasma parameters similarly to
alpha in section 4.3 of arxiv.org/abs/1212.3496.
*/
template <
	class Solver_Info,
	class Cell_Iter,
	class Grid,
	class Target_Refinement_Level_Getter,
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Magnetic_Field_Getter,
	class Bg_Magnetic_Field_Getter,
	class Is_Primary_Face_Getter,
	class Solver_Info_Getter
> void get_target_refinement_levels(
	const Cell_Iter& cells,
	Grid& grid,
	const Target_Refinement_Level_Getter Ref,
	const Mass_Density_Getter Mas,
	const Momentum_Density_Getter Mom,
	const Total_Energy_Density_Getter Nrj,
	const Magnetic_Field_Getter Mag,
	const Bg_Magnetic_Field_Getter Bg_Face_Mag,
	const Is_Primary_Face_Getter PFace,
	const Solver_Info_Getter Sol_Info,
	const double adiabatic_index,
	const double vacuum_permeability
) {
	using std::abs;
	using std::max;
	using std::min;
	using std::pow;
	using std::sqrt;

	const auto max_ref_lvl = grid.get_maximum_refinement_level();
	if (max_ref_lvl == 0) {
		return;
	}

	for (const auto& cell: cells) {
		if ((Sol_Info(*cell.data) & pamhd::mhd::Solver_Info::dont_solve) > 0) {
			continue;
		}

		const auto cpface = PFace(*cell.data);
		const auto cref_lvl = grid.mapping.get_refinement_level(cell.id);
		const auto
			cmas = Mas(*cell.data),
			cnrj = Nrj(*cell.data);
		const typename std::remove_reference<decltype(Mom(*cell.data))>::type
			cmom = Mom(*cell.data),
			cvel = {cmom[0]/cmas, cmom[1]/cmas, cmom[2]/cmas},
			cvel2 = {cvel[0]*cvel[0], cvel[1]*cvel[1], cvel[2]*cvel[2]},
			cmag = Mag(*cell.data);

		for (const auto& neighbor: cell.neighbors_of) {
			const auto n = neighbor.face_neighbor;
			if (n == 0) {
				continue;
			}

			if ((Sol_Info(*neighbor.data) & pamhd::mhd::Solver_Info::dont_solve) > 0) {
				continue;
			}

			const size_t dim = size_t(abs(n) - 1);

			const auto bg_face_mag = [&]{
				if (n < 0) {
					if (cpface[2*dim]) {
						return Bg_Face_Mag(*cell.data)(dim, 0);
					} else {
						return Bg_Face_Mag(*neighbor.data)(dim, 1);
					}
				} else {
					if (cpface[1 + 2*dim]) {
						return Bg_Face_Mag(*cell.data)(dim, 1);
					} else {
						return Bg_Face_Mag(*neighbor.data)(dim, 0);
					}
				}
			}();

			const auto
				nmas = Mas(*neighbor.data),
				dmas = cmas - nmas,
				mas_min = min(cmas, nmas),
				nnrj = Nrj(*neighbor.data),
				dnrj = cnrj - nnrj,
				nrj_min = min(cnrj, nnrj);
			const decltype(cmom)
				nmom = Mom(*neighbor.data),
				dmom = {cmom[0]-nmom[0], cmom[1]-nmom[1], cmom[2]-nmom[2]},
				nvel = {nmom[0]/nmas, nmom[1]/nmas, nmom[2]/nmas},
				nvel2 = {nvel[0]*nvel[0], nvel[1]*nvel[1], nvel[2]*nvel[2]},
				dvel = {cvel[0]-nvel[0], cvel[1]-nvel[1], cvel[2]-nvel[2]},
				nmag = Mag(*neighbor.data),
				dmag = {cmag[0]-nmag[0], cmag[1]-nmag[1], cmag[2]-nmag[2]};
			const auto
				dmom2 = pow(dmom[0],2) + pow(dmom[1],2) + pow(dmom[2],2),
				mom_min = 2*min(cmas*cnrj, nmas*nnrj),
				mag_min1 = 2*vacuum_permeability*min(cnrj, nnrj),
				dmag2 = sqrt(pow(dmag[0],2) + pow(dmag[1],2) + pow(dmag[2],2)),
				mag_min2 = sqrt(min(
					cmag[0]*cmag[0] + cmag[1]*cmag[1] + cmag[2]*cmag[2],
					nmag[0]*nmag[0] + nmag[1]*nmag[1] + nmag[2]*nmag[2])),
				dmag2_mag_min2 = [&]{
					if (mag_min2 <= 0) return 0.0;
					else return dmag2/mag_min2;}(),
				dvel2 = pow(dvel[0],2) + pow(dvel[1],2) + pow(dvel[2],2),
				cfast = get_fast_magnetosonic_speed(
					cmas, cmom, cnrj, cmag, bg_face_mag,
					adiabatic_index, vacuum_permeability),
				nfast = get_fast_magnetosonic_speed(
					nmas, nmom, nnrj, nmag, bg_face_mag,
					adiabatic_index, vacuum_permeability),
				vel_min
					= min(
						cvel2[0]+cvel2[1]+cvel2[2],
						nvel2[0]+nvel2[1]+nvel2[2])
					+ pow(10*max(cfast, nfast),2);

			const auto ref_index = max(max(max(max(max(
				dmas/mas_min,
				dnrj/nrj_min),
				dmom2/mom_min),
				dmag2*dmag2/mag_min1),
				dmag2_mag_min2),
				dvel2/vel_min);
			const auto tgt_ref_lvl = [&]{
				if (max_ref_lvl <= 0) {
					return max_ref_lvl;
				}
				const auto rel_lvls = 0.2 * (cref_lvl + 1) / max_ref_lvl;
				if (ref_index > rel_lvls) {
					return min(cref_lvl + 1, max_ref_lvl);
				}
				if (ref_index < 0.5 * rel_lvls) {
					return max(0, cref_lvl - 1);
				}
				return cref_lvl;
			}();
			Ref(*cell.data) = max(tgt_ref_lvl, Ref(*cell.data));
		}
	}
}


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
> void get_minmax_refinement_level(
	const Cells& cells,
	Grid& grid,
	Options& options_mhd,
	const Mass_Density_Getter& Mas,
	const Momentum_Density_Getter& Mom,
	const Total_Energy_Density_Getter& Nrj,
	const Volume_Magnetic_Field_Getter& Mag,
	const Solver_Info_Getter& Solver_Info,
	const Target_Refinement_Level_Min_Getter& RLMin,
	const Target_Refinement_Level_Max_Getter& RLMax,
	const double& adiabatic_index,
	const double& vacuum_permeability
) {
	using std::abs;
	using std::clamp;
	using std::max;
	using std::min;
	using std::round;
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

		const auto [cdx, cdy, cdz] = grid.geometry.get_length(cell.id);
		const auto
			cpre = max(options_mhd.pressure_min_mrg, get_pressure(
				Mas(*cell.data), Mom(*cell.data), Nrj(*cell.data),
				Mag(*cell.data), adiabatic_index, vacuum_permeability)),
			cmas = max(Mas(*cell.data), options_mhd.number_density_min_mrg),
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
				nmas = max(Mas(*neighbor.data), options_mhd.number_density_min_mrg),
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
					throw std::runtime_error(__FILE__ ":" + to_string(__LINE__));
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
		RLMin(*cell.data) = clamp(
			max(int(round(tgt_ref_lvl-0.25)), RLMin(*cell.data)),
			0, max_ref_lvl);
		RLMax(*cell.data) = clamp(
			min(int(round(tgt_ref_lvl+0.25)), RLMax(*cell.data)),
			RLMin(*cell.data), max_ref_lvl);
	}
}


// Handlers that handle all parameters required by MHD test program.
template <
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Face_Magnetic_Field_Getter_Pos,
	class Face_Magnetic_Field_Getter_Neg,
	class Volume_Magnetic_Field_Getter,
	class Background_Magnetic_Field_Getter,
	class Background_Magnetic_Field,
	class Target_Refinement_Level_Min_Getter,
	class Target_Refinement_Level_Max_Getter
> struct New_Cells_Handler {
	const Mass_Density_Getter& Mas;
	const Momentum_Density_Getter& Mom;
	const Total_Energy_Density_Getter& Nrj;
	const Face_Magnetic_Field_Getter_Pos& Face_Bp;
	const Face_Magnetic_Field_Getter_Neg& Face_Bn;
	const Volume_Magnetic_Field_Getter& Vol_B;
	const Background_Magnetic_Field_Getter& Bg_B;
	const Background_Magnetic_Field& bg_B;
	const Target_Refinement_Level_Min_Getter& RLMin;
	const Target_Refinement_Level_Max_Getter& RLMax;
	const double& adiabatic_index;
	const double& vacuum_permeability;

	New_Cells_Handler(
		const Mass_Density_Getter& Mas_,
		const Momentum_Density_Getter& Mom_,
		const Total_Energy_Density_Getter& Nrj_,
		const Face_Magnetic_Field_Getter_Pos& Face_Bp_,
		const Face_Magnetic_Field_Getter_Neg& Face_Bn_,
		const Volume_Magnetic_Field_Getter& Vol_B_,
		const Background_Magnetic_Field_Getter& Bg_B_,
		const Background_Magnetic_Field& bg_B_,
		const Target_Refinement_Level_Min_Getter& RLMin_,
		const Target_Refinement_Level_Max_Getter& RLMax_,
		const double& adiabatic_index_,
		const double& vacuum_permeability_
	) :
		Mas(Mas_), Mom(Mom_), Nrj(Nrj_),
		Face_Bp(Face_Bp_), Face_Bn(Face_Bn_),
		Vol_B(Vol_B_), Bg_B(Bg_B_), bg_B(bg_B_),
		RLMin(RLMin_), RLMax(RLMax_),
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

			// inherit thermal pressure
			const double parent_pressure = get_pressure(
				Mas(*parent_data),
				Mom(*parent_data),
				Nrj(*parent_data),
				Vol_B(*parent_data),
				adiabatic_index,
				vacuum_permeability
			);

			const auto
				cindex = grid.mapping.get_indices(new_cell_id),
				pindex = grid.mapping.get_indices(parent_id);
			for (size_t dim = 0; dim < 3; dim++) {
				if (cindex[dim] == pindex[dim]) { // neg faces coincide
					Face_Bn(*cell_data)[dim] = Face_Bn(*parent_data)[dim];
					Face_Bp(*cell_data)[dim] = Vol_B(*parent_data)[dim];
				} else { // pos faces coincide
					Face_Bp(*cell_data)[dim] = Face_Bp(*parent_data)[dim];
					Face_Bn(*cell_data)[dim] = Vol_B(*parent_data)[dim];
				}
			}

			Vol_B(*cell_data) = 0.5*(Face_Bn(*cell_data) + Face_Bp(*cell_data));
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
			Bg_B(*cell_data)(0, 0) = bg_B.get_background_field(
				{sx, ry, rz},
				vacuum_permeability
			);
			Bg_B(*cell_data)(0, 1) = bg_B.get_background_field(
				{ex, ry, rz},
				vacuum_permeability
			);
			Bg_B(*cell_data)(1, 0) = bg_B.get_background_field(
				{rx, sy, rz},
				vacuum_permeability
			);
			Bg_B(*cell_data)(1, 1) = bg_B.get_background_field(
				{rx, ey, rz},
				vacuum_permeability
			);
			Bg_B(*cell_data)(2, 0) = bg_B.get_background_field(
				{rx, ry, sz},
				vacuum_permeability
			);
			Bg_B(*cell_data)(2, 1) = bg_B.get_background_field(
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
	class Face_Magnetic_Field_Getter_Pos,
	class Face_Magnetic_Field_Getter_Neg,
	class Volume_Magnetic_Field_Getter,
	class Background_Magnetic_Field_Getter,
	class Background_Magnetic_Field,
	class Target_Refinement_Level_Min_Getter,
	class Target_Refinement_Level_Max_Getter
> struct Removed_Cells_Handler {
	const Mass_Density_Getter& Mas;
	const Momentum_Density_Getter& Mom;
	const Total_Energy_Density_Getter& Nrj;
	const Face_Magnetic_Field_Getter_Pos& Face_Bp;
	const Face_Magnetic_Field_Getter_Neg& Face_Bn;
	const Volume_Magnetic_Field_Getter& Vol_B;
	const Background_Magnetic_Field_Getter& Bg_B;
	const Background_Magnetic_Field& bg_B;
	const Target_Refinement_Level_Min_Getter& RLMin;
	const Target_Refinement_Level_Max_Getter& RLMax;
	const double& adiabatic_index;
	const double& vacuum_permeability;

	Removed_Cells_Handler(
		const Mass_Density_Getter& Mas_,
		const Momentum_Density_Getter& Mom_,
		const Total_Energy_Density_Getter& Nrj_,
		const Face_Magnetic_Field_Getter_Pos& Face_Bp_,
		const Face_Magnetic_Field_Getter_Neg& Face_Bn_,
		const Volume_Magnetic_Field_Getter& Vol_B_,
		const Background_Magnetic_Field_Getter& Bg_B_,
		const Background_Magnetic_Field& bg_B_,
		const Target_Refinement_Level_Min_Getter& RLMin_,
		const Target_Refinement_Level_Max_Getter& RLMax_,
		const double& adiabatic_index_,
		const double& vacuum_permeability_
	) :
		Mas(Mas_), Mom(Mom_), Nrj(Nrj_),
		Face_Bp(Face_Bp_), Face_Bn(Face_Bn_),
		Vol_B(Vol_B_), Bg_B(Bg_B_), bg_B(bg_B_),
		RLMin(RLMin_), RLMax(RLMax_),
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
			Mas(*parent_data) =
			Nrj(*parent_data) = 0;
			Mom(*parent_data)     =
			Face_Bn(*parent_data) =
			Face_Bp(*parent_data) = {0, 0, 0};
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

			// all children included
			RLMin(*parent_data) = std::min(RLMin(*removed_cell_data), RLMin(*parent_data));
			RLMax(*parent_data) = std::max(RLMax(*removed_cell_data), RLMax(*parent_data));
			Mas(*parent_data) += Mas(*removed_cell_data) / 8;
			Mom(*parent_data) += Mom(*removed_cell_data) / 8;
			// temporarily store pressure in energy density variable
			Nrj(*parent_data) += get_pressure(
				Mas(*removed_cell_data),
				Mom(*removed_cell_data),
				Nrj(*removed_cell_data),
				Vol_B(*removed_cell_data),
				adiabatic_index,
				vacuum_permeability
			) / 8;;

			const auto
				pindex = grid.mapping.get_indices(parent_id),
				rindex = grid.mapping.get_indices(removed_cell_id);
			for (size_t dim = 0; dim < 3; dim++) {
				// only children sharing a face with parent included
				if (pindex[dim] == rindex[dim]) {
					Face_Bn(*parent_data)[dim] += Face_Bn(*removed_cell_data)[dim] / 4;
				} else {
					Face_Bp(*parent_data)[dim] += Face_Bp(*removed_cell_data)[dim] / 4;
				}
			}

			const auto [rx, ry, rz] = grid.geometry.get_center(parent_id);
			const auto [sx, sy, sz] = grid.geometry.get_min(parent_id);
			const auto [ex, ey, ez] = grid.geometry.get_max(parent_id);
			Bg_B(*parent_data)(0, 0) = bg_B.get_background_field(
				{sx, ry, rz},
				vacuum_permeability
			);
			Bg_B(*parent_data)(0, 1) = bg_B.get_background_field(
				{ex, ry, rz},
				vacuum_permeability
			);
			Bg_B(*parent_data)(1, 0) = bg_B.get_background_field(
				{rx, sy, rz},
				vacuum_permeability
			);
			Bg_B(*parent_data)(1, 1) = bg_B.get_background_field(
				{rx, ey, rz},
				vacuum_permeability
			);
			Bg_B(*parent_data)(2, 0) = bg_B.get_background_field(
				{rx, ry, sz},
				vacuum_permeability
			);
			Bg_B(*parent_data)(2, 1) = bg_B.get_background_field(
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

			Vol_B(*parent_data) = 0.5*(Face_Bn(*parent_data) + Face_Bp(*parent_data));

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


}} // namespaces


#endif // ifndef PAMHD_MHD_AMR_HPP
