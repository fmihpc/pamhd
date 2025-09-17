/*
Common functions of PAMHD related to temporal substepping.

Copyright 2025 Finnish Meteorological Institute
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

#ifndef PAMHD_SUBSTEPPING_HPP
#define PAMHD_SUBSTEPPING_HPP


#include "cmath"
#include "limits"
#include "string"

#include "mhd/options.hpp"
#include "common_variables.hpp"


namespace pamhd {


/*! Converts substepping periods from 2^N to N format.

Returns new largest N.
*/
template <
	class Grid,
	class Solver_Info_Getter,
	class Substepping_Period_Getter
> int update_substeps(
	Grid& grid,
	const Solver_Info_Getter SInfo,
	const Substepping_Period_Getter Substep
) try {
	Substep.type().is_stale = true;

	int max_local = -1;
	for (const auto& cell: grid.local_cells()) {
		if (SInfo.data(*cell.data) < 0) {
			continue;
		}
		Substep.data(*cell.data) = 1 << Substep.data(*cell.data);
		max_local = std::max(Substep.data(*cell.data), max_local);
	}

	int max_global = -1;
	auto comm = grid.get_communicator();
	if (
		MPI_Allreduce(
			&max_local,
			&max_global,
			1,
			MPI_INT,
			MPI_MAX,
			comm
		) != MPI_SUCCESS
	) {
		std::cerr << __FILE__ "(" << __LINE__
			<< "): Couldn't reduce max substep." << std::endl;
		abort();
	}
	MPI_Comm_free(&comm);

	return max_global;

} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


/*!
Reduces difference in max substep between neighbors to <= 1.

Assumes 2^N format for substep periods.
*/
template <
	class Grid,
	class Substepping_Period_Getter,
	class Substep_Max_Getter,
	class Solver_Info_Getter
> void restrict_substepping_period(
	Grid& grid,
	const Substepping_Period_Getter Substep,
	const Substep_Max_Getter Substep_Max,
	const Solver_Info_Getter SInfo
) try {
	using std::max;
	using std::min;

	using Cell = Grid::cell_data_type;

	for (const auto& cell: grid.local_cells()) {
		if (SInfo.data(*cell.data) < 0) {
			continue;
		}
		Substep.data(*cell.data) = Substep_Max.data(*cell.data);
	}

	if (SInfo.type().is_stale) {
		Cell::set_transfer_all(true, SInfo.type());
		grid.update_copies_of_remote_neighbors();
		Cell::set_transfer_all(false, SInfo.type());
		SInfo.type().is_stale = false;
	}

	auto comm = grid.get_communicator();
	Cell::set_transfer_all(true, Substep.type());
	uint64_t modified_cells = 0;
	do {
		grid.update_copies_of_remote_neighbors();
		uint64_t modified_cells_local = 0;
		for (const auto& cell: grid.local_cells()) {
			if (SInfo.data(*cell.data) < 0) {
				continue;
			}
			for (const auto& neighbor: cell.neighbors_of) {
				if (SInfo.data(*neighbor.data) < 0) {
					continue;
				}
				if (neighbor.face_neighbor == 0 and neighbor.edge_neighbor[0] < 0) {
					continue;
				}

				if (Substep.data(*cell.data) > Substep.data(*neighbor.data) + 1) {
					Substep.data(*cell.data) = Substep.data(*neighbor.data) + 1;
					modified_cells_local++;
				}
			}
		}
		if (
			MPI_Allreduce(
				&modified_cells_local,
				&modified_cells,
				1,
				MPI_UINT64_T,
				MPI_SUM,
				comm
			) != MPI_SUCCESS
		) {
			std::cerr << __FILE__ "(" << __LINE__
				<< "): Couldn't reduce modified cell count."
				<< std::endl;
			abort();
		}
	} while (modified_cells > 0);
	grid.update_copies_of_remote_neighbors();
	Cell::set_transfer_all(false, Substep.type());

	MPI_Comm_free(&comm);

} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


/*!
Sets min and max substepping periods based on geometry.

Offsets minimum substep so that at it will be 0 in at
least one cell.

Assumes that substep is stored in 2^N format.
*/
template <
	class Grid,
	class Options,
	class Substep_Min_Getter,
	class Substep_Max_Getter
> void set_minmax_substepping_period(
	double time,
	Grid& grid,
	Options& options,
	const Substep_Min_Getter Substep_Min,
	const Substep_Max_Getter Substep_Max
) try {
	int min_substep_min_local = std::numeric_limits<int>::max();
	for (const auto& cell: grid.local_cells()) {
		const auto c = grid.geometry.get_center(cell.id);
		const auto
			r = std::sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]),
			lat = std::asin(c[2] / r),
			lon = std::atan2(c[1], c[0]);
		if (options.substep_min_i >= 0) {
			Substep_Min.data(*cell.data) = options.substep_min_i;
		} else {
			Substep_Min.data(*cell.data) = options.substep_min_e.evaluate(
				time, c[0], c[1], c[2], r, lat, lon);
		}
		if (options.substep_max_i >= 0) {
			Substep_Max.data(*cell.data) = options.substep_max_i;
		} else {
			Substep_Max.data(*cell.data) = options.substep_max_e.evaluate(
				time, c[0], c[1], c[2], r, lat, lon);
		}
		min_substep_min_local = std::min(Substep_Min.data(*cell.data), min_substep_min_local);
	}
	Substep_Max.type().is_stale = true;
	Substep_Min.type().is_stale = true;

	int min_substep_min_global = std::numeric_limits<int>::max();
	auto comm = grid.get_communicator();
	if (
		MPI_Allreduce(
			&min_substep_min_local,
			&min_substep_min_global,
			1,
			MPI_INT,
			MPI_MIN,
			comm
		) != MPI_SUCCESS
	) {
		std::cerr << __FILE__ "(" << __LINE__
			<< "): Couldn't reduce min_substep_min."
			<< std::endl;
		abort();
	}
	MPI_Comm_free(&comm);

	if (min_substep_min_global > 0) {
		for (const auto& cell: grid.local_cells()) {
			Substep_Min.data(*cell.data) -= min_substep_min_global;
		}
	}
} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


/*! TODO

Assumes substep is in 2^N format indicating period of solving
a cell, or in other words, indicating length of cell's timestep
in substeps.

If max_dt is smaller than smallest allowed timestep (including
dt_factor) in any cell, returns max_dt and sets substep range
to [0, 0] in all cells.
Otherwise returns largest substep length satisfying each cell's
combination of timestep and minimum substep period.
For example if one cell's timestep is dt and minimum substep is
0 and another cell's timestep is 1.5*dt and minimum substep is
1, returns 1.5*dt/2.
*/
template <
	class Grid,
	class Solver_Info_Getter,
	class Timestep_Getter,
	class Substep_Min_Getter,
	class Substep_Max_Getter
> double set_minmax_substepping_period(
	Grid& grid,
	const double& max_dt,
	const Solver_Info_Getter& SInfo,
	const Timestep_Getter& Timestep,
	const Substep_Min_Getter& Substep_Min,
	const Substep_Max_Getter& Substep_Max
) try {
	using std::cerr;
	using std::endl;
	using std::max;
	using std::min;
	using std::numeric_limits;

	Substep_Max.type().is_stale = true;
	Substep_Min.type().is_stale = true;

	// calculate range of dts allowed in cells
	double
		smallest_dt_local = numeric_limits<double>::max(),
		largest_dt_local = 0;

	for (const auto& cell: grid.local_cells()) {
		if (SInfo.data(*cell.data) < 0) {
			continue;
		}
		smallest_dt_local = min(Timestep.data(*cell.data), smallest_dt_local);
		largest_dt_local = max(Timestep.data(*cell.data), largest_dt_local);
	}
	auto comm = grid.get_communicator();
	double smallest_dt_global = -1;
	if (
		MPI_Allreduce(
			&smallest_dt_local,
			&smallest_dt_global,
			1,
			MPI_DOUBLE,
			MPI_MIN,
			comm
		) != MPI_SUCCESS
	) {
		cerr << __FILE__ "(" << __LINE__
			<< "): Couldn't reduce smallest_dt." << endl;
		abort();
	}

	if (max_dt <= smallest_dt_global) {
		for (const auto& cell: grid.local_cells()) {
			Substep_Max.data(*cell.data) =
			Substep_Min.data(*cell.data) = 0;
		}
		MPI_Comm_free(&comm);
		return max_dt;
	}

	// minimum substep length
	double ret_val_local = numeric_limits<double>::max();
	for (const auto& cell: grid.local_cells()) {
		if (SInfo.data(*cell.data) < 0) {
			continue;
		}
		ret_val_local = min(ret_val_local,
			Timestep.data(*cell.data) / (1 << Substep_Min.data(*cell.data)));
	}

	double ret_val_global = numeric_limits<double>::max();
	if (
		MPI_Allreduce(
			&ret_val_local,
			&ret_val_global,
			1,
			MPI_DOUBLE,
			MPI_MIN,
			comm
		) != MPI_SUCCESS
	) {
		cerr << __FILE__ "(" << __LINE__
			<< "): Couldn't reduce ret_val." << endl;
		abort();
	}
	MPI_Comm_free(&comm);

	// decrease too large max substep periods
	for (const auto& cell: grid.local_cells()) {
		if (SInfo.data(*cell.data) < 0) {
			continue;
		}

		const auto& dt = Timestep.data(*cell.data);
		if (dt > 0) {
			int step = max(0, (int)std::floor(std::log2(min(dt, max_dt) / ret_val_global)));
			if (Substep_Max.data(*cell.data) > step) {
				Substep_Max.data(*cell.data) = step;
			}
			if (Substep_Min.data(*cell.data) > Substep_Max.data(*cell.data)) {
				cerr << "Unexpected substeps in cell "
					<< cell.id << ": "
					<< Substep_Min.data(*cell.data) << ", "
					<< Substep_Max.data(*cell.data) << ", "
					<< dt << ", " << ret_val_global << endl;
				abort();
			}
		} else {
			Substep_Min.data(*cell.data) =
			Substep_Max.data(*cell.data) = 0;
		}
	}

	return ret_val_global;

} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


} // namespace


#endif // ifndef PAMHD_SUBSTEPPING_HPP
