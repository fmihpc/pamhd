/*
Grid stuff related to adaptive mesh refinement of PAMHD.

Copyright 2023, 2024, 2025 Finnish Meteorological Institute
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

* Neither the name of copyright holders nor the names of their contributors
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

#ifndef PAMHD_GRID_AMR_HPP
#define PAMHD_GRID_AMR_HPP


#include "algorithm"
#include "initializer_list"
#include "set"
#include "stdexcept"
#include "string"
#include "tuple"
#include "type_traits"

#include "grid/options.hpp"


namespace pamhd {
namespace grid {


/*! Variable indicating minimum refinement level towards which grid cells should be refined.
*/
struct Target_Refinement_Level_Min {
	static bool is_stale;
	using data_type = int;
	static const std::string get_name() { return {"target minimum refinement level"}; }
	static const std::string get_option_name() { return {"target-min-ref-lvl"}; }
	static const std::string get_option_help() { return {"Target minimum refinement level of cell"}; }
};

struct Target_Refinement_Level_Max {
	static bool is_stale;
	using data_type = int;
	static const std::string get_name() { return {"target maximum refinement level"}; }
	static const std::string get_option_name() { return {"target-max-ref-lvl"}; }
	static const std::string get_option_help() { return {"Target maximum refinement level of cell"}; }
};


//! Returns smallest [0] and largest [1] refinement level for given cell.
std::array<int, 2> get_target_refinement_level(
	Options& options,
	const double& t,
	const std::array<double, 3>& cell_center
) {
	using std::runtime_error;
	using std::to_string;

	const auto& c = cell_center;
	const auto
		r = std::sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]),
		lat = std::asin(c[2] / r),
		lon = std::atan2(c[1], c[0]);

	const auto min_ref_lvl = options.get_ref_lvl_at_least(
		t, c[0], c[1], c[2], r, lat, lon);
	if (min_ref_lvl < 0) {
		throw runtime_error(__FILE__ "(" + to_string(__LINE__) + "): options returned min_ref_lvl < 0");
	}

	const auto max_ref_lvl = options.get_ref_lvl_at_most(
		t, c[0], c[1], c[2], r, lat, lon);
	if (max_ref_lvl < 0) {
		throw runtime_error(__FILE__ "(" + to_string(__LINE__) + "): options returned min_ref_lvl < 0");
	}
	if (max_ref_lvl < min_ref_lvl) {
		throw runtime_error(__FILE__ "(" + to_string(__LINE__) + "): options returned max_ref_lvl < min_ref_lvl");
	}

	std::array<int, 2> ret_val{min_ref_lvl, max_ref_lvl};
	return ret_val;
}

// Handlers that only handle target refinement levels.
template<
	class Target_Refinement_Level_Min_Getter,
	class Target_Refinement_Level_Max_Getter
> struct New_Cells_Handler {
	const Target_Refinement_Level_Min_Getter& RLMin;
	const Target_Refinement_Level_Max_Getter& RLMax;

	New_Cells_Handler(
		const Target_Refinement_Level_Min_Getter& RLMin_,
		const Target_Refinement_Level_Max_Getter& RLMax_
	) :
		RLMin(RLMin_), RLMax(RLMax_)
	{};

	template <
		class Grid, class Cells
	> void operator()(
		const Grid& grid, const Cells& new_cells
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
			RLMin.data(*cell_data) = RLMin.data(*parent_data);
			RLMax.data(*cell_data) = RLMax.data(*parent_data);
		}
	}
};

template<
	class Target_Refinement_Level_Min_Getter,
	class Target_Refinement_Level_Max_Getter
> struct Removed_Cells_Handler {
	const Target_Refinement_Level_Min_Getter& RLMin;
	const Target_Refinement_Level_Max_Getter& RLMax;

	Removed_Cells_Handler(
		const Target_Refinement_Level_Min_Getter& RLMin_,
		const Target_Refinement_Level_Max_Getter& RLMax_
	) :
		RLMin(RLMin_), RLMax(RLMax_)
	{};

	template <
		class Grid, class Cells
	> void operator()(
		const Grid& grid, const Cells& removed_cells
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
			RLMin.data(*parent_data) = 999;
			RLMax.data(*parent_data) = 0;
		}

		// parent's refinement range spans all of childrens'
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
			RLMin.data(*parent_data) = std::min(RLMin.data(*removed_cell_data), RLMin.data(*parent_data));
			RLMax.data(*parent_data) = std::max(RLMax.data(*removed_cell_data), RLMax.data(*parent_data));
		}
	}
};

/*! Adapts given grid based on target min and max refinement level.

Grid is adapted at most grid.get_maximum_refinement_level() times.
Cells whose refinement level is smaller than target minimum are
refined.
Cells whose refinement level is larger than target maximum are
unrefined.

Handlers must have following operator() signature:
nch(grid, new_cells);
rch(grid, removed_cells);
where
const auto new_cells = grid.stop_refining();
const auto removed_cells = grid.get_removed_cells();
*/
template <
	class Grid,
	class Target_Refinement_Level_Min_Getter,
	class Target_Refinement_Level_Max_Getter,
	class New_Cells_Handler,
	class Removed_Cells_Handler
> void adapt_grid(
	Grid& grid,
	const Target_Refinement_Level_Min_Getter& RLMin,
	const Target_Refinement_Level_Max_Getter& RLMax,
	const New_Cells_Handler& nch,
	const Removed_Cells_Handler& rch,
	const int max_rounds = 1 << 30
) {
	#ifdef MPI_VERSION
	using Cell = Grid::cell_data_type;
	Cell::set_transfer_all(true, RLMin.type(), RLMax.type());
	MPI_Comm comm = grid.get_communicator();
	#endif
	for (
		int i = 0;
		i < std::min(max_rounds, grid.get_maximum_refinement_level());
		i++
	) {
		for (const auto& cell: grid.local_cells()) {
			const auto ref_lvl = grid.get_refinement_level(cell.id);
			if (ref_lvl < RLMin.data(*cell.data)) {
				grid.refine_completely(cell.id);
			} else if (ref_lvl == RLMin.data(*cell.data)) {
				grid.dont_unrefine(cell.id);
			} else if (ref_lvl > RLMax.data(*cell.data)) {
				grid.unrefine_completely(cell.id);
			}
		}
		const auto new_cells = grid.stop_refining();
		const auto removed_cells = grid.get_removed_cells();

		const uint64_t total_size_local = new_cells.size() + removed_cells.size();
		uint64_t total_size_global = total_size_local;
		#ifdef MPI_VERSION
		total_size_global = dccrg::All_Reduce()(total_size_local, comm);
		#endif

		if (total_size_global == 0) {
			break;
		} else {
			if (new_cells.size() > 0) nch(grid, new_cells);
			if (removed_cells.size() > 0) rch(grid, removed_cells);
			grid.clear_refined_unrefined_data();
		}
	}
	#ifdef MPI_VERSION
	MPI_Comm_free(&comm);
	grid.update_copies_of_remote_neighbors();
	Cell::set_transfer_all(false, RLMin.type(), RLMax.type());
	#endif
}

/*! Potentially emulates 1d/2d grid

If length of grid is initially only one cell in a dimension, pretends
that all cells are flat and in middle of grid in that dimension.

Returns
\verbatim
array<
	array<grid.geometry.get_center()>
	array<...get_min()>,
	array<...get_max()>
>.
\endverbatim
*/
std::array<std::array<double, 3>, 3> get_cell_geom_emulated(
	const auto& grid, const auto& cell
) {
	const auto lvl0 = grid.mapping.length.get();
	const auto
		grid_start = grid.geometry.get_start(),
		grid_end = grid.geometry.get_end();
	decltype(grid_start) grid_mid{
		(grid_end[0] + grid_start[0]) / 2,
		(grid_end[1] + grid_start[1]) / 2,
		(grid_end[2] + grid_start[2]) / 2
	};
	auto
		cctr = grid.geometry.get_center(cell),
		cmin = grid.geometry.get_min(cell),
		cmax = grid.geometry.get_max(cell);
	for (size_t dim: {0, 1, 2}) {
		if (lvl0[dim] == 1) {
			cctr[dim] = cmin[dim] = cmax[dim] = grid_mid[dim];
		}
	}
	return {cctr, cmin, cmax};
}


/*! Sets min,max refinement levels based on geometry.

If clamped == true then target minimum refinement level
is not decreased and tgt max ref lvl is not increased.
*/
template <
	class Cells,
	class Grid,
	class Target_Refinement_Level_Min_Getter,
	class Target_Refinement_Level_Max_Getter
> void set_minmax_refinement_level(
	const Cells& cells,
	Grid& grid,
	pamhd::grid::Options& options,
	const double& sim_time,
	const Target_Refinement_Level_Min_Getter& RLMin,
	const Target_Refinement_Level_Max_Getter& RLMax,
	bool clamped
) {
	using std::clamp;

	for (const auto& cell: cells) {
		const auto [center, start, end]
			= get_cell_geom_emulated(grid, cell.id);
		const auto [rlmin, rlmax] = get_target_refinement_level(options, sim_time, center);
		if (clamped) {
			const auto
				prev_min = RLMin.data(*cell.data),
				prev_max = RLMax.data(*cell.data);
			RLMin.data(*cell.data) = clamp(rlmin, prev_min, prev_max);
			RLMax.data(*cell.data) = clamp(rlmax, prev_min, prev_max);
		} else {
			RLMin.data(*cell.data) = rlmin;
			RLMax.data(*cell.data) = rlmax;
		}
	}
}

}} // namespaces

#endif // ifndef PAMHD_GRID_AMR_HPP
