/*
Particle data accumulator of PAMHD built on top of DCCRG.

Copyright 2015, 2016, 2017 Ilja Honkonen
Copyright 2019 Finnish Meteorological Institute
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

#ifndef PAMHD_PARTICLE_ACCUMULATE_DCCRG_HPP
#define PAMHD_PARTICLE_ACCUMULATE_DCCRG_HPP


#include "algorithm"
#include "iostream"
#include "utility"
#include "vector"

#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#include "Eigen/Core"

#include "mhd/common.hpp"
#include "particle/accumulate.hpp"
#include "particle/common.hpp"
#include "particle/variables.hpp"


namespace pamhd {
namespace particle {


/*! Assumes geometry API is compatible with cartesian dccrg geometry */
template<class Geometry> std::tuple<
	Eigen::Vector3d,
	Eigen::Vector3d,
	Eigen::Vector3d,
	Eigen::Vector3d
> get_cell_geometry(
	const uint64_t cell_id,
	const Geometry& geometry
) {
	const auto
		cell_min_tmp = geometry.get_min(cell_id),
		cell_max_tmp = geometry.get_max(cell_id),
		cell_length_tmp = geometry.get_length(cell_id),
		cell_center_tmp = geometry.get_center(cell_id);
	const Eigen::Vector3d
		cell_min{cell_min_tmp[0], cell_min_tmp[1], cell_min_tmp[2]},
		cell_max{cell_max_tmp[0], cell_max_tmp[1], cell_max_tmp[2]},
		cell_length{cell_length_tmp[0], cell_length_tmp[1], cell_length_tmp[2]},
		cell_center{cell_center_tmp[0], cell_center_tmp[1], cell_center_tmp[2]};
	return std::make_tuple(cell_min, cell_max, cell_length, cell_center);
}


/*!
Accumulates particle data in given cells to those cells and their neighbors.

Accumulated data from particles in the same cell is not zeroed before accumulating
unless clear_at_start == true;

Accumulated data is written to local cells directly, accumulated data to
remote neighbors is stored in local cells' accumulation list.

Updates the number of remote accumulated values.

Particle data to/from dont_solve cells isn't accumulated.

Grid is assumed to be dccrg with Cartesian geometry.
*/
template<
	class Cell_Iterator,
	class Grid,
	class Particles_Getter,
	class Particle_Position_Getter,
	class Particle_Value_Getter,
	class Bulk_Value_Getter,
	class Bulk_Value_In_List_Getter,
	class Target_In_List_Getter,
	class Accumulation_List_Length_Getter,
	class Accumulation_List_Getter,
	class Solver_Info_Getter
> void accumulate(
	const Cell_Iterator& cells,
	Grid& grid,
	Particles_Getter Part,
	Particle_Position_Getter Part_Pos,
	Particle_Value_Getter Part_Val,
	Bulk_Value_Getter Bulk_Val,
	Bulk_Value_In_List_Getter List_Bulk_Val,
	Target_In_List_Getter List_Target,
	Accumulation_List_Length_Getter List_Len,
	Accumulation_List_Getter Accu_List,
	Solver_Info_Getter SInfo,
	const bool clear_at_start = true
) {
	const auto
		grid_start = grid.geometry.get_start(),
		grid_end = grid.geometry.get_end();
	const Eigen::Vector3d grid_len{
		grid_end[0] - grid_start[0],
		grid_end[1] - grid_start[1],
		grid_end[2] - grid_start[2]
	};

	for (const auto& cell: cells) {
		if (SInfo(*cell.data) < 0) {
			Bulk_Val(*cell.data) = {};
			continue;
		}

		Eigen::Vector3d cell_min, cell_max, cell_length, cell_center;
		std::tie(cell_min, cell_max, cell_length, cell_center) = get_cell_geometry(cell.id, grid.geometry);

		if (clear_at_start) {
			Accu_List(*cell.data).clear();
		}

		// TODO: faster to iterate over neighbors first?
		for (auto& particle: Part(*cell.data)) {
			auto& position = Part_Pos(particle);
			const Eigen::Vector3d
				value_box_min{
					position[0] - cell_length[0] / 2,
					position[1] - cell_length[1] / 2,
					position[2] - cell_length[2] / 2
				},
				value_box_max{
					position[0] + cell_length[0] / 2,
					position[1] + cell_length[1] / 2,
					position[2] + cell_length[2] / 2
				};

			// accumulate to current cell
			Bulk_Val(*cell.data)
				+= get_accumulated_value(
					Part_Val(*cell.data, particle),
					value_box_min,
					value_box_max,
					cell_min,
					cell_max
				);

			// accumulate to neighbors
			for (const auto& neighbor: cell.neighbors_of) {

				Eigen::Vector3d neigh_min, neigh_max, neigh_length, neigh_center;
				std::tie(neigh_min, neigh_max, neigh_length, neigh_center) = get_cell_geometry(neighbor.id, grid.geometry);

				// handle periodic grid
				if (neighbor.x < 0 and neigh_center[0] > cell_min[0]) {
					neigh_min[0] -= grid_len[0];
					neigh_max[0] -= grid_len[0];
				} else if (neighbor.x > 0 and neigh_center[0] < cell_max[0]) {
					neigh_min[0] += grid_len[0];
					neigh_max[0] += grid_len[0];
				}
				if (neighbor.y < 0 and neigh_center[1] > cell_min[1]) {
					neigh_min[1] -= grid_len[1];
					neigh_max[1] -= grid_len[1];
				} else if (neighbor.y > 0 and neigh_center[1] < cell_max[1]) {
					neigh_min[1] += grid_len[1];
					neigh_max[1] += grid_len[1];
				}
				if (neighbor.z < 0 and neigh_center[2] > cell_min[2]) {
					neigh_min[2] -= grid_len[2];
					neigh_max[2] -= grid_len[2];
				} else if (neighbor.z > 0 and neigh_center[2] < cell_max[2]) {
					neigh_min[2] += grid_len[2];
					neigh_max[2] += grid_len[2];
				}

				if (
					value_box_min[0] > neigh_max[0]
					or value_box_min[1] > neigh_max[1]
					or value_box_min[2] > neigh_max[2]
					or value_box_max[0] < neigh_min[0]
					or value_box_max[1] < neigh_min[1]
					or value_box_max[2] < neigh_min[2]
				) {
					continue;
				}

				// don't accumulate into dont_solve cells
				if (SInfo(*neighbor.data) < 0) {
					continue;
				}

				const auto accumulated_value
					= get_accumulated_value(
						Part_Val(*neighbor.data, particle),
						value_box_min,
						value_box_max,
						neigh_min,
						neigh_max
					);

				// same as for current cell above
				if (neighbor.is_local) {
					Bulk_Val(*neighbor.data) += accumulated_value;
				// accumulate values to a list in current cell
				} else {
					// use this index in the list
					size_t accumulation_index = 0;

					// find the index of target neighbor
					auto iter
						= std::find_if(
							Accu_List(*cell.data).begin(),
							Accu_List(*cell.data).end(),
							[&neighbor, &List_Target](
								const decltype(*Accu_List(*cell.data).begin()) candidate_item
							) {
								if (List_Target(candidate_item) == neighbor.id) {
									return true;
								} else {
									return false;
								}
							}
						);

					// found
					if (iter != Accu_List(*cell.data).end()) {
						List_Bulk_Val(*iter) += accumulated_value;
					// create the item
					} else {
						const auto old_size = Accu_List(*cell.data).size();
						Accu_List(*cell.data).resize(old_size + 1);
						auto& new_item = Accu_List(*cell.data)[old_size];
						List_Target(new_item) = neighbor.id;
						List_Bulk_Val(new_item) = accumulated_value;
					}
				}
			}
		}

		List_Len(*cell.data) = Accu_List(*cell.data).size();
	}
}


/*!
Same as accumulate() but uses get_accumulated_value_weighted() instead of
get_accumulated_value() and also records the total weight from all particles.
*/
template<
	class Cell_Iterator,
	class Grid,
	class Particles_Getter,
	class Particle_Position_Getter,
	class Particle_Value_Getter,
	class Particle_Weight_Getter,
	class Bulk_Value_Getter,
	class Bulk_Value_In_List_Getter,
	class Target_In_List_Getter,
	class Accumulation_List_Length_Getter,
	class Accumulation_List_Getter,
	class Solver_Info_Getter
> void accumulate_weighted(
	const Cell_Iterator& cells,
	Grid& grid,
	Particles_Getter Part,
	Particle_Position_Getter Part_Pos,
	Particle_Value_Getter Part_Val,
	Particle_Weight_Getter Part_Wei,
	Bulk_Value_Getter Bulk_Val,
	Bulk_Value_In_List_Getter List_Bulk_Val,
	Target_In_List_Getter List_Target,
	Accumulation_List_Length_Getter List_Len,
	Accumulation_List_Getter Accu_List,
	Solver_Info_Getter SInfo,
	const bool clear_at_start = true
) {
	using std::is_same_v;
	using Eigen::Vector3d;

	const auto
		grid_start = grid.geometry.get_start(),
		grid_end = grid.geometry.get_end();
	const Vector3d grid_len{
		grid_end[0] - grid_start[0],
		grid_end[1] - grid_start[1],
		grid_end[2] - grid_start[2]
	};

	for (const auto& cell: cells) {
		if (SInfo(*cell.data) < 0) {
			using bulk_val_t = std::remove_reference_t<
				decltype(Bulk_Val(*cell.data))>;
			if constexpr (is_same_v<bulk_val_t, Vector3d>) {
				Bulk_Val(*cell.data) = Vector3d::Zero();
			} else if constexpr (
				is_same_v<bulk_val_t, std::pair<Vector3d, double>>
			) {
				Bulk_Val(*cell.data) = {Vector3d::Zero(), 0};
			} else {
				Bulk_Val(*cell.data) = {};
			}
			continue;
		}

		Vector3d cell_min, cell_max, cell_length, cell_center;
		std::tie(
			cell_min, cell_max, cell_length, cell_center
		) = get_cell_geometry(cell.id, grid.geometry);

		if (clear_at_start) {
			Accu_List(*cell.data).clear();
		}

		for (auto& particle: Part(*cell.data)) {
			auto& position = Part_Pos(particle);
			const Vector3d
				value_box_min{
					position[0] - cell_length[0] / 2,
					position[1] - cell_length[1] / 2,
					position[2] - cell_length[2] / 2
				},
				value_box_max{
					position[0] + cell_length[0] / 2,
					position[1] + cell_length[1] / 2,
					position[2] + cell_length[2] / 2
				};

			// accumulate to current cell
			const auto accu_to_cell
				= get_accumulated_value_weighted(
					Part_Val(*cell.data, particle),
					Part_Wei(particle),
					value_box_min,
					value_box_max,
					cell_min,
					cell_max
				);
			Bulk_Val(*cell.data).first += accu_to_cell.first;
			// final weight of this particle's data in this cell
			Bulk_Val(*cell.data).second += accu_to_cell.second;

			// accumulate to neighbors
			for (const auto& neighbor: cell.neighbors_of) {

				Vector3d neigh_min, neigh_max, neigh_length, neigh_center;
				std::tie(neigh_min, neigh_max, neigh_length, neigh_center) = get_cell_geometry(neighbor.id, grid.geometry);

				// handle periodic grid
				if (neighbor.x < 0 and neigh_center[0] > cell_min[0]) {
					neigh_min[0] -= grid_len[0];
					neigh_max[0] -= grid_len[0];
				} else if (neighbor.x > 0 and neigh_center[0] < cell_max[0]) {
					neigh_min[0] += grid_len[0];
					neigh_max[0] += grid_len[0];
				}
				if (neighbor.y < 0 and neigh_center[1] > cell_min[1]) {
					neigh_min[1] -= grid_len[1];
					neigh_max[1] -= grid_len[1];
				} else if (neighbor.y > 0 and neigh_center[1] < cell_max[1]) {
					neigh_min[1] += grid_len[1];
					neigh_max[1] += grid_len[1];
				}
				if (neighbor.z < 0 and neigh_center[2] > cell_min[2]) {
					neigh_min[2] -= grid_len[2];
					neigh_max[2] -= grid_len[2];
				} else if (neighbor.z > 0 and neigh_center[2] < cell_max[2]) {
					neigh_min[2] += grid_len[2];
					neigh_max[2] += grid_len[2];
				}

				if (
					value_box_min[0] > neigh_max[0]
					or value_box_min[1] > neigh_max[1]
					or value_box_min[2] > neigh_max[2]
					or value_box_max[0] < neigh_min[0]
					or value_box_max[1] < neigh_min[1]
					or value_box_max[2] < neigh_min[2]
				) {
					continue;
				}

				// don't accumulate into dont_solve cells
				if (SInfo(*neighbor.data) < 0) {
					continue;
				}

				const auto accumulated_value
					= get_accumulated_value_weighted(
						Part_Val(*neighbor.data, particle),
						Part_Wei(particle),
						value_box_min,
						value_box_max,
						neigh_min,
						neigh_max
					);

				// same as for current cell above
				if (neighbor.is_local) {
					Bulk_Val(*neighbor.data).first += accumulated_value.first;
					Bulk_Val(*neighbor.data).second += accumulated_value.second;
				// accumulate values to a list in current cell
				} else {
					// use this index in the list
					size_t accumulation_index = 0;

					// find the index of target neighbor
					auto iter
						= std::find_if(
							Accu_List(*cell.data).begin(),
							Accu_List(*cell.data).end(),
							[&neighbor, &List_Target](
								const decltype(*Accu_List(*cell.data).begin()) candidate_item
							) {
								if (List_Target(candidate_item) == neighbor.id) {
									return true;
								} else {
									return false;
								}
							}
						);

					// found
					if (iter != Accu_List(*cell.data).end()) {
						List_Bulk_Val(*iter).first += accumulated_value.first;
						List_Bulk_Val(*iter).second += accumulated_value.second;
					// create the item
					} else {
						const auto old_size = Accu_List(*cell.data).size();
						Accu_List(*cell.data).resize(old_size + 1);
						auto& new_item = Accu_List(*cell.data)[old_size];
						List_Target(new_item) = neighbor.id;
						List_Bulk_Val(new_item).first = accumulated_value.first;
						List_Bulk_Val(new_item).second = accumulated_value.second;
					}
				}
			}
		}

		List_Len(*cell.data) = Accu_List(*cell.data).size();
	}
}


/*!
Allocates memory for lists of accumulated particle data in copies of remote neighbors of local cells.
*/
template<
	class Grid,
	class Accumulation_List_Getter,
	class Accumulation_List_Length_Getter
> void allocate_accumulation_lists(
	Grid& grid,
	Accumulation_List_Getter Accu_List,
	Accumulation_List_Length_Getter List_Len
) {
	for (const auto& cell: grid.remote_cells()) {
		Accu_List(*cell.data).resize(List_Len(*cell.data));
	}
}


/*!
Adds accumulated particle data from remote neighbors' to local cells.
*/
template<
	class Grid,
	class Bulk_Value_Getter,
	class Value_In_List_Getter,
	class Target_In_List_Getter,
	class Accumulation_List_Getter,
	class Solver_Info_Getter
> void accumulate_from_remote_neighbors(
	Grid& grid,
	Bulk_Value_Getter Bulk_Val,
	Value_In_List_Getter Value_In_List,
	Target_In_List_Getter Target_In_List,
	Accumulation_List_Getter Accu_List,
	Solver_Info_Getter SInfo
) {
	for (const auto& remote_cell_id: grid.get_remote_cells_on_process_boundary()) {
		auto* const source_data = grid[remote_cell_id];
		if (source_data == nullptr) {
			std::cerr << __FILE__ << "(" << __LINE__ << ")" << std::endl;
			abort();
		}

		if (SInfo(*source_data) < 0) {
			continue;
		}

		for (auto& item: Accu_List(*source_data)) {
			if (not grid.is_local(Target_In_List(item))) {
				continue;
			}

			auto* const target_data = grid[Target_In_List(item)];
			if (target_data == nullptr) {
				std::cerr << __FILE__ << "(" << __LINE__ << ")" << std::endl;
				abort();
			}

			Bulk_Val(*target_data) += Value_In_List(item);
		}
	}
}

template<
	class Grid,
	class Bulk_Value_Getter,
	class Value_In_List_Getter,
	class Target_In_List_Getter,
	class Accumulation_List_Getter,
	class Solver_Info_Getter
> void accumulate_weighted_from_remote_neighbors(
	Grid& grid,
	Bulk_Value_Getter Bulk_Val,
	Value_In_List_Getter Value_In_List,
	Target_In_List_Getter Target_In_List,
	Accumulation_List_Getter Accu_List,
	Solver_Info_Getter SInfo
) {
	for (const auto& remote_cell_id: grid.get_remote_cells_on_process_boundary()) {
		auto* const source_data = grid[remote_cell_id];
		if (source_data == nullptr) {
			std::cerr << __FILE__ << "(" << __LINE__ << ")" << std::endl;
			abort();
		}

		if (SInfo(*source_data) < 0) {
			continue;
		}

		for (auto& item: Accu_List(*source_data)) {
			if (not grid.is_local(Target_In_List(item))) {
				continue;
			}

			auto* const target_data = grid[Target_In_List(item)];
			if (target_data == nullptr) {
				std::cerr << __FILE__ << "(" << __LINE__ << ")" << std::endl;
				abort();
			}

			Bulk_Val(*target_data).first += Value_In_List(item).first;
			Bulk_Val(*target_data).second += Value_In_List(item).second;
		}
	}
}


/*!
Updates Bulk_Momentum.
*/
template<
	class Grid,
	class Particles_Getter,
	class Particle_Position_Getter,
	class Particle_Mass_Getter,
	class Particle_Species_Mass_Getter,
	class Particle_Velocity_Getter,
	class Particle_Relative_Kinetic_Energy_Getter,
	class Number_Of_Particles_Getter,
	class Particle_Weight_Getter,
	class Bulk_Mass_Getter,
	class Bulk_Momentum_Getter,
	class Bulk_Relative_Kinetic_Energy_Getter,
	class Number_Of_Particles_In_List_Getter,
	class Bulk_Mass_In_List_Getter,
	class Bulk_Momentum_In_List_Getter,
	class Bulk_Relative_Kinetic_Energy_In_List_Getter,
	class Bulk_Velocity_Getter,
	class Target_In_List_Getter,
	class Accumulation_List_Length_Getter,
	class Accumulation_List_Getter,
	class Accumulation_List_Length_Variable,
	class Accumulation_List_Variable,
	class Bulk_Velocity_Variable,
	class Solver_Info_Getter
> void accumulate_mhd_data(
	Grid& grid,
	const Particles_Getter& Particles,
	const Particle_Position_Getter& Particle_Position,
	const Particle_Mass_Getter& Particle_Mass,
	const Particle_Species_Mass_Getter& Particle_Species_Mass,
	const Particle_Velocity_Getter& Particle_Velocity,
	const Particle_Relative_Kinetic_Energy_Getter& Particle_Relative_Kinetic_Energy,
	const Number_Of_Particles_Getter& Number_Of_Particles,
	const Particle_Weight_Getter& Particle_Weight,
	const Bulk_Mass_Getter& Bulk_Mass,
	const Bulk_Momentum_Getter& Bulk_Momentum,
	const Bulk_Relative_Kinetic_Energy_Getter& Bulk_Relative_Kinetic_Energy,
	const Bulk_Velocity_Getter& Bulk_Velocity,
	const Number_Of_Particles_In_List_Getter& Accu_List_Number_Of_Particles,
	const Bulk_Mass_In_List_Getter& Accu_List_Bulk_Mass,
	const Bulk_Momentum_In_List_Getter& Accu_List_Bulk_Velocity,
	const Bulk_Relative_Kinetic_Energy_In_List_Getter& Accu_List_Bulk_Relative_Kinetic_Energy,
	const Target_In_List_Getter& Accu_List_Target,
	const Accumulation_List_Length_Getter& Accu_List_Length,
	const Accumulation_List_Getter& Accu_List,
	const Accumulation_List_Length_Variable& accu_list_len_var,
	const Accumulation_List_Variable& accu_list_var,
	const Bulk_Velocity_Variable& bulk_vel_var,
	const Solver_Info_Getter& SInfo
) {
	for (const auto& cell: grid.local_cells()) {
		Bulk_Mass(*cell.data) = 0;
		Bulk_Momentum(*cell.data) = {0, 0, 0};
		Bulk_Velocity(*cell.data).first = {0, 0, 0};
		Bulk_Velocity(*cell.data).second = 0;
	}

	accumulate(
		grid.outer_cells(),
		grid,
		Particles,
		Particle_Position,
		Particle_Mass,
		Bulk_Mass,
		Accu_List_Bulk_Mass,
		Accu_List_Target,
		Accu_List_Length,
		Accu_List,
		SInfo
	);
	accumulate_weighted(
		grid.outer_cells(),
		grid,
		Particles,
		Particle_Position,
		Particle_Velocity,
		Particle_Weight,
		Bulk_Velocity,
		Accu_List_Bulk_Velocity,
		Accu_List_Target,
		Accu_List_Length,
		Accu_List,
		SInfo,
		false // keep previous data in accumulation lists
	);

	Grid::cell_data_type::set_transfer_all(true, accu_list_len_var);
	grid.start_remote_neighbor_copy_updates();

	accumulate(
		grid.inner_cells(),
		grid,
		Particles,
		Particle_Position,
		Particle_Mass,
		Bulk_Mass,
		Accu_List_Bulk_Mass,
		Accu_List_Target,
		Accu_List_Length,
		Accu_List,
		SInfo,
		false
	);
	accumulate_weighted(
		grid.inner_cells(),
		grid,
		Particles,
		Particle_Position,
		Particle_Velocity,
		Particle_Weight,
		Bulk_Velocity,
		Accu_List_Bulk_Velocity,
		Accu_List_Target,
		Accu_List_Length,
		Accu_List,
		SInfo,
		false
	);

	grid.wait_remote_neighbor_copy_update_receives();

	allocate_accumulation_lists(
		grid,
		Accu_List,
		Accu_List_Length
	);

	grid.wait_remote_neighbor_copy_update_sends();
	Grid::cell_data_type::set_transfer_all(false, accu_list_len_var);

	Grid::cell_data_type::set_transfer_all(true, accu_list_var);
	grid.start_remote_neighbor_copy_updates();
	grid.wait_remote_neighbor_copy_update_receives();

	accumulate_from_remote_neighbors(
		grid,
		Bulk_Mass,
		Accu_List_Bulk_Mass,
		Accu_List_Target,
		Accu_List,
		SInfo
	);
	accumulate_weighted_from_remote_neighbors(
		grid,
		Bulk_Velocity,
		Accu_List_Bulk_Velocity,
		Accu_List_Target,
		Accu_List,
		SInfo
	);

	grid.wait_remote_neighbor_copy_update_sends();
	Grid::cell_data_type::set_transfer_all(false, accu_list_var);

	// scale velocities relative to total weights
	for (const auto& cell: grid.local_cells()) {
		if (Bulk_Velocity(*cell.data).second <= 0) {
			Bulk_Velocity(*cell.data).first = {0, 0, 0};
		} else {
			Bulk_Velocity(*cell.data).first /= Bulk_Velocity(*cell.data).second;
		}
		Bulk_Momentum(*cell.data) = Bulk_Velocity(*cell.data).first * Bulk_Mass(*cell.data);
	}


	/*
	Accumulate particle data required for pressure
	*/

	// needs remote neighbors' bulk velocity
	Grid::cell_data_type::set_transfer_all(true, bulk_vel_var);
	grid.start_remote_neighbor_copy_updates();

	for (const auto& cell: grid.local_cells()) {
		Number_Of_Particles(*cell.data)          =
		Bulk_Relative_Kinetic_Energy(*cell.data) = 0;
	}

	grid.wait_remote_neighbor_copy_update_receives();

	// returns number of particles of mass species mass
	const auto particle_counter
		= [
			&Particle_Mass,
			&Particle_Species_Mass
		](
			typename Grid::cell_data_type& cell,
			// TODO:switch to auto in c++17
			decltype(*Particles(*(grid[0])).begin())& particle
		) ->double {
			return Particle_Mass(cell, particle) / Particle_Species_Mass(cell, particle);
		};
	// TODO: merge with bulk mass accumulation
	accumulate(
		grid.outer_cells(),
		grid,
		Particles,
		Particle_Position,
		particle_counter,
		Number_Of_Particles,
		Accu_List_Number_Of_Particles,
		Accu_List_Target,
		Accu_List_Length,
		Accu_List,
		SInfo
	);
	accumulate(
		grid.outer_cells(),
		grid,
		Particles,
		Particle_Position,
		Particle_Relative_Kinetic_Energy,
		Bulk_Relative_Kinetic_Energy,
		Accu_List_Bulk_Relative_Kinetic_Energy,
		Accu_List_Target,
		Accu_List_Length,
		Accu_List,
		SInfo,
		false
	);

	grid.wait_remote_neighbor_copy_update_sends();
	Grid::cell_data_type::set_transfer_all(false, bulk_vel_var);


	Grid::cell_data_type::set_transfer_all(true, accu_list_len_var);
	grid.start_remote_neighbor_copy_updates();

	accumulate(
		grid.inner_cells(),
		grid,
		Particles,
		Particle_Position,
		particle_counter,
		Number_Of_Particles,
		Accu_List_Number_Of_Particles,
		Accu_List_Target,
		Accu_List_Length,
		Accu_List,
		SInfo,
		false
	);
	accumulate(
		grid.inner_cells(),
		grid,
		Particles,
		Particle_Position,
		Particle_Relative_Kinetic_Energy,
		Bulk_Relative_Kinetic_Energy,
		Accu_List_Bulk_Relative_Kinetic_Energy,
		Accu_List_Target,
		Accu_List_Length,
		Accu_List,
		SInfo,
		false
	);

	grid.wait_remote_neighbor_copy_update_receives();

	allocate_accumulation_lists(
		grid,
		Accu_List,
		Accu_List_Length
	);

	grid.wait_remote_neighbor_copy_update_sends();
	Grid::cell_data_type::set_transfer_all(false, accu_list_len_var);

	Grid::cell_data_type::set_transfer_all(true, accu_list_var);
	grid.start_remote_neighbor_copy_updates();
	grid.wait_remote_neighbor_copy_update_receives();

	accumulate_from_remote_neighbors(
		grid,
		Number_Of_Particles,
		Accu_List_Number_Of_Particles,
		Accu_List_Target,
		Accu_List,
		SInfo
	);
	accumulate_from_remote_neighbors(
		grid,
		Bulk_Relative_Kinetic_Energy,
		Accu_List_Bulk_Relative_Kinetic_Energy,
		Accu_List_Target,
		Accu_List,
		SInfo
	);

	grid.wait_remote_neighbor_copy_update_sends();
	Grid::cell_data_type::set_transfer_all(false, accu_list_var);
}


/*!
Converts bulk particle data to conservative MHD form.

Skips a cell if it doesn't have particles and no_particles_allowed == true.
*/
template<
	class Grid,
	class Number_Of_Particles_Getter,
	class Particle_Bulk_Mass_Getter,
	class Particle_Bulk_Momentum_Getter,
	class Particle_Bulk_Relative_Kinetic_Energy_Getter,
	class Particle_List_Getter,
	class MHD_Mass_Getter,
	class MHD_Momentum_Getter,
	class MHD_Energy_Getter,
	class MHD_Magnetic_Field_Getter,
	class Solver_Info_Getter
> void fill_mhd_fluid_values(
	Grid& grid,
	const double adiabatic_index,
	const double vacuum_permeability,
	const double particle_temp_nrj_ratio,
	const double minimum_pressure,
	Number_Of_Particles_Getter Number_Of_Particles,
	Particle_Bulk_Mass_Getter Particle_Bulk_Mass,
	Particle_Bulk_Momentum_Getter Particle_Bulk_Momentum,
	Particle_Bulk_Relative_Kinetic_Energy_Getter Particle_Bulk_Relative_Kinetic_Energy,
	Particle_List_Getter Particle_List,
	MHD_Mass_Getter MHD_Mass,
	MHD_Momentum_Getter MHD_Momentum,
	MHD_Energy_Getter MHD_Energy,
	MHD_Magnetic_Field_Getter MHD_Magnetic_Field,
	Solver_Info_Getter SInfo
) {
	for (const auto& cell: grid.local_cells()) {
		if (SInfo(*cell.data) < 0) {
			MHD_Mass(*cell.data) = 0;
			MHD_Momentum(*cell.data) = {0, 0, 0};
			MHD_Energy(*cell.data) = 0;
			continue;
		}

		if (Particle_Bulk_Mass(*cell.data) <= 0) {
			MHD_Mass(*cell.data) = 0;
			MHD_Momentum(*cell.data) = {0, 0, 0};
			MHD_Energy(*cell.data) = 0;
			continue;
		}

		const auto length = grid.geometry.get_length(cell.id);
		const auto volume = length[0] * length[1] * length[2];

		MHD_Mass(*cell.data) = Particle_Bulk_Mass(*cell.data) / volume;
		MHD_Momentum(*cell.data) = Particle_Bulk_Momentum(*cell.data) / volume;

		const double pressure = [&](){
			if (Particle_Bulk_Mass(*cell.data) <= 0) {
				return 0.0;
			} else {
				return std::max(
					minimum_pressure,
					2 * Particle_Bulk_Relative_Kinetic_Energy(*cell.data) / 3 / volume
				);
			}
		}();

		if (MHD_Mass(*cell.data) <= 0) {
			MHD_Energy(*cell.data)
				= pressure / (adiabatic_index - 1)
				+ MHD_Magnetic_Field(*cell.data).squaredNorm() / vacuum_permeability / 2;
		} else {
			MHD_Energy(*cell.data)
				= pressure / (adiabatic_index - 1)
				+ MHD_Momentum(*cell.data).squaredNorm() / MHD_Mass(*cell.data) / 2
				+ MHD_Magnetic_Field(*cell.data).squaredNorm() / vacuum_permeability / 2;
		}
	}
}


}} // namespaces

#endif // ifndef PAMHD_PARTICLE_ACCUMULATE_DCCRG_HPP
