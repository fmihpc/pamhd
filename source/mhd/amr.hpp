/*
Adaptive mesh refinement logic of MHD part of PAMHD

Copyright 2023 Finnish Meteorological Institute
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

#ifndef PAMHD_MHD_AMR_HPP
#define PAMHD_MHD_AMR_HPP


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

		const auto primary = PFace(*cell.data);
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
					if (primary[2*dim]) {
						return Bg_Face_Mag(*cell.data)(dim, 0);
					} else {
						return Bg_Face_Mag(*neighbor.data)(dim, 1);
					}
				} else {
					if (primary[1 + 2*dim]) {
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


/*! Refines/unrefines given cells.
*/
template <
	class Solver_Info,
	class Cell_Iter,
	class Grid,
	class Target_Refinement_Level_Getter,
	class Solver_Info_Getter
> std::vector<uint64_t> adapt_grid(
	const Cell_Iter& cells,
	Grid& grid,
	const Target_Refinement_Level_Getter Ref,
	const Solver_Info_Getter Sol_Info
) {
	using std::to_string;

	for (const auto& cell: cells) {
		if ((Sol_Info(*cell.data) & pamhd::mhd::Solver_Info::dont_solve) > 0) {
			continue;
		}

		const auto cref_lvl = grid.mapping.get_refinement_level(cell.id);
		const auto tgt_lvl = Ref(*cell.data);
		if (tgt_lvl < cref_lvl) {
			grid.unrefine_completely(cell.id);
		} else if (tgt_lvl > cref_lvl) {
			grid.refine_completely(cell.id);
		} else {
			grid.dont_unrefine(cell.id);
		}
	}

	return grid.stop_refining();
}


/*! Initializes new cells from adaptation.

Sets data of new cells based on their parent if parent
was refined or their children if some were unrefined.

-Magnetic fields should be consistent before calling adapt_grid().
-Should be called soon after adapt_grid().
-PFace variable should be initialized before calling this.
*/
template <
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Face_Magnetic_Field_Getter_Pos,
	class Face_Magnetic_Field_Getter_Neg,
	class Target_Refinement_Level_Getter,
	class Solver_Info_Getter
> struct New_Cells_Handler {
	const Mass_Density_Getter& Mas;
	const Momentum_Density_Getter& Mom;
	const Total_Energy_Density_Getter& Nrj;
	const Face_Magnetic_Field_Getter_Pos& Face_Bp;
	const Face_Magnetic_Field_Getter_Neg& Face_Bn;
	const Target_Refinement_Level_Getter& Ref;
	const Solver_Info_Getter& Sol_Info;

	New_Cells_Handler(
		const Mass_Density_Getter& Mas_,
		const Momentum_Density_Getter& Mom_,
		const Total_Energy_Density_Getter& Nrj_,
		const Face_Magnetic_Field_Getter_Pos& Face_Bp_,
		const Face_Magnetic_Field_Getter_Neg& Face_Bn_,
		const Target_Refinement_Level_Getter& Ref_,
		const Solver_Info_Getter& Sol_Info_
	) :
		Mas(Mas_), Mom(Mom_), Nrj(Nrj_),
		Face_Bp(Face_Bp_), Face_Bn(Face_Bn_),
		Ref(Ref_), Sol_Info(Sol_Info_)
	{};

	template<
		class Grid,
		class Cells
	> void operator()(
		const Grid& grid,
		const Cells& new_cells
	) const {
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

			Sol_Info(*cell_data) = Sol_Info(*parent_data);
			Mas(*cell_data) = Mas(*parent_data);
			Mom(*cell_data) = Mom(*parent_data);
			Nrj(*cell_data) = Nrj(*parent_data);
			Ref(*cell_data) = Ref(*parent_data);

			const auto
				cindex = grid.mapping.get_indices(new_cell_id),
				pindex = grid.mapping.get_indices(parent_id);
			for (size_t dim = 0; dim < 3; dim++) {
				const auto avg_b = 0.5*(Face_Bn(*parent_data)[dim] + Face_Bp(*parent_data)[dim]);
				if (cindex[dim] == pindex[dim]) { // neg faces coincide
					Face_Bn(*cell_data)[dim] = Face_Bn(*parent_data)[dim];
					Face_Bp(*cell_data)[dim] = avg_b;
				} else { // pos faces coincide
					Face_Bp(*cell_data)[dim] = Face_Bp(*parent_data)[dim];
					Face_Bn(*cell_data)[dim] = avg_b;
				}
			}
		}
	}
};


template <
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Face_Magnetic_Field_Getter_Pos,
	class Face_Magnetic_Field_Getter_Neg,
	class Target_Refinement_Level_Getter
> struct Removed_Cells_Handler {
	const Mass_Density_Getter& Mas;
	const Momentum_Density_Getter& Mom;
	const Total_Energy_Density_Getter& Nrj;
	const Face_Magnetic_Field_Getter_Pos& Face_Bp;
	const Face_Magnetic_Field_Getter_Neg& Face_Bn;
	const Target_Refinement_Level_Getter& Ref;

	Removed_Cells_Handler(
		const Mass_Density_Getter& Mas_,
		const Momentum_Density_Getter& Mom_,
		const Total_Energy_Density_Getter& Nrj_,
		const Face_Magnetic_Field_Getter_Pos& Face_Bp_,
		const Face_Magnetic_Field_Getter_Neg& Face_Bn_,
		const Target_Refinement_Level_Getter& Ref_
	) :
		Mas(Mas_), Mom(Mom_), Nrj(Nrj_),
		Face_Bp(Face_Bp_), Face_Bn(Face_Bn_),
		Ref(Ref_)
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
			Mas(*parent_data) += Mas(*removed_cell_data) / 8;
			Mom(*parent_data) += Mom(*removed_cell_data) / 8;
			Nrj(*parent_data) += Nrj(*removed_cell_data) / 8;

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
		}
	}
};


}} // namespaces


#endif // ifndef PAMHD_MHD_AMR_HPP
