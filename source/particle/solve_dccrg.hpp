/*
Particle propagator of PAMHD for DCCRG grid.

Copyright 2015, 2016, 2017 Ilja Honkonen
Copyright 2019, 2024, 2025 Finnish Meteorological Institute
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

#ifndef PAMHD_PARTICLE_SOLVE_DCCRG_HPP
#define PAMHD_PARTICLE_SOLVE_DCCRG_HPP


#include "array"
#include "cstdlib"
#include "exception"
#include "type_traits"
#include "utility"

#include "dccrg.hpp"
#include "Eigen/Core"
#include "Eigen/Geometry"
#include "prettyprint.hpp"

#include "accumulate_dccrg.hpp"
#include "common.hpp"
#include "interpolate.hpp"
#include "math/staggered.hpp"
#include "mhd/solve_staggered.hpp"
#include "particle/solve.hpp"
#include "substepping.hpp"
#include "variables.hpp"


namespace pamhd {
namespace particle {


/*!
Propagates particles in given cells for a given amount of time.

Returns the longest allowed time step for given cells
and their neighbors. Particles which end up outside of the
cell in which they are stored are moved to the External_Particles_T
list of their previous cell and added to Particle_Destinations_T
information.

Assumes grid was initialized with neighbhorhood size of 1 and
maximum refinement level of 0.
*/
template<
	class Cell_Iterator,
	class Grid,
	class Background_Magnetic_Field,
	class Current_Minus_Velocity_Getter,
	class Volume_Magnetic_Field_Getter,
	class Nr_Particles_External_Getter,
	class Particles_Internal_Getter,
	class Particles_External_Getter,
	class Particle_Max_Spatial_Velocity_Getter,
	class Particle_Max_Angular_Velocity_Getter,
	class Particle_Position_Getter,
	class Particle_Velocity_Getter,
	class Particle_Charge_Mass_Ratio_Getter,
	class Particle_Mass_Getter,
	class Particle_Destination_Cell_Getter,
	class Solver_Info_Getter
> void solve(
	const double& dt,
	const Cell_Iterator& cells,
	Grid& grid,
	const Background_Magnetic_Field& bg_B,
	const double& vacuum_permeability,
	const bool& E_is_derived_quantity,
	// if E_is_derived_quantity == true: JmV = J - V, else JmV = E
	const Current_Minus_Velocity_Getter& JmV,
	const Volume_Magnetic_Field_Getter& Vol_B,
	const Nr_Particles_External_Getter& Nr_Ext,
	const Particles_Internal_Getter& Part_Int,
	const Particles_External_Getter& Part_Ext,
	const Particle_Max_Spatial_Velocity_Getter& Max_v_part,
	const Particle_Max_Angular_Velocity_Getter& Max_ω_part,
	const Particle_Position_Getter& Part_Pos,
	const Particle_Velocity_Getter& Part_Vel,
	const Particle_Charge_Mass_Ratio_Getter& Part_C2M,
	const Particle_Mass_Getter& Part_Mas,
	const Particle_Destination_Cell_Getter& Part_Des,
	const Solver_Info_Getter& SInfo
) {
	using std::abs;
	using std::isnan;
	using std::is_same;
	using std::max;
	using std::min;

	using Cell = Grid::cell_data_type;

	if (grid.get_maximum_refinement_level() > 0) {
		throw std::runtime_error("Only maximum refinement level 0 supported.");
	}

	bool update_copies = false;
	if (JmV.type().is_stale) {
		update_copies = true;
		Cell::set_transfer_all(true, JmV.type());
		JmV.type().is_stale = false;
	}
	if (Vol_B.type().is_stale) {
		update_copies = true;
		Cell::set_transfer_all(true, Vol_B.type());
		Vol_B.type().is_stale = false;
	}
	if (update_copies) {
		grid.update_copies_of_remote_neighbors();
	}
	Cell::set_transfer_all(false, JmV.type(), Vol_B.type());

	std::pair<double, double> max_time_step{std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
	const std::array<double, 3>
		grid_start{grid.geometry.get_start()},
		grid_end{grid.geometry.get_end()},
		grid_center{
			(grid_end[0]-grid_start[0]) / 2,
			(grid_end[1]-grid_start[1]) / 2,
			(grid_end[2]-grid_start[2]) / 2
		};
	for (const auto& cell: cells) {
		Max_v_part.data(*cell.data) =
		Max_ω_part.data(*cell.data) = 0;
		if (SInfo.data(*cell.data) < 0) {
			Part_Int(*cell.data).clear();
			Part_Ext.data(*cell.data).clear();
			Nr_Ext.data(*cell.data) = 0;
			continue;
		}

		// get field data from neighborhood for interpolation
		std::array<Eigen::Vector3d, 27> current_minus_velocities, magnetic_fields;

		// default to current cell's data
		current_minus_velocities.fill(JmV.data(*cell.data));
		magnetic_fields.fill(Vol_B.data(*cell.data));

		for (const auto& neighbor: cell.neighbors_of) {
			if (
				abs(neighbor.x) > 1
				or abs(neighbor.y) > 1
				or abs(neighbor.z) > 1
			) {
				continue;
			}
			if (SInfo.data(*neighbor.data) < 0) {
				continue;
			}

			const size_t index = (neighbor.z + 1) * 9 + (neighbor.y + 1) * 3 + neighbor.x + 1;
			current_minus_velocities[index] = JmV.data(*neighbor.data);
			magnetic_fields[index] = Vol_B.data(*neighbor.data);
		}

		const auto
			cell_min = grid.geometry.get_min(cell.id),
			cell_max = grid.geometry.get_max(cell.id),
			cell_center = grid.geometry.get_center(cell.id),
			cell_length = grid.geometry.get_length(cell.id);

		const Eigen::Vector3d
			interpolation_start{
				cell_center[0] - cell_length[0],
				cell_center[1] - cell_length[1],
				cell_center[2] - cell_length[2]
			},
			interpolation_end{
				cell_center[0] + cell_length[0],
				cell_center[1] + cell_length[1],
				cell_center[2] + cell_length[2]
			};

		for (size_t part_i = 0; part_i < Part_Int(*cell.data).size(); part_i++) {
			auto
				pos = Part_Pos(Part_Int(*cell.data)[part_i]),
				// emulate 2d,1d,0d sim from B0 perspective
				bg_pos = pos;
			const auto lvl0 = grid.mapping.length.get();
			for (size_t dim = 0; dim < 3; dim++) {
				if (lvl0[dim] == 1) {
					bg_pos[dim] = grid_center[dim];
				}
			}
			const Eigen::Vector3d B_at_pos
				= interpolate(
					pos, interpolation_start,
					interpolation_end, magnetic_fields)
				+ bg_B.get_background_field(
					bg_pos, vacuum_permeability);

			auto vel = Part_Vel(Part_Int(*cell.data)[part_i]);
			Max_v_part.data(*cell.data) = max(
				Max_v_part.data(*cell.data),
				vel.norm());
			const auto& c2m = Part_C2M(Part_Int(*cell.data)[part_i]);
			Max_ω_part.data(*cell.data) = max(
				Max_ω_part.data(*cell.data),
				abs(c2m) * B_at_pos.norm());

			const Eigen::Matrix<double,3,1> E_at_pos = [&](){
				if (E_is_derived_quantity) {
					const auto J_m_V_at_pos = interpolate(
						pos, interpolation_start,
						interpolation_end, current_minus_velocities
					);
					return J_m_V_at_pos.cross(B_at_pos);
				} else {
					return interpolate(
						pos, interpolation_start,
						interpolation_end, current_minus_velocities
					);
				}
			}();
			std::tie(pos, vel) = propagate(
				pos, vel, E_at_pos, B_at_pos, c2m, dt
			);
			Part_Vel(Part_Int(*cell.data)[part_i]) = vel;

			// take into account periodic grid
			const auto real_pos
				= grid.geometry.get_real_coordinate(
					{pos[0], pos[1], pos[2]});

			Part_Pos(Part_Int(*cell.data)[part_i]) = {
				real_pos[0], real_pos[1], real_pos[2]
			};

			// remove from simulation if particle not inside of grid
			if (
				isnan(real_pos[0])
				or isnan(real_pos[1])
				or isnan(real_pos[2])
			) {

				Part_Int(*cell.data).erase(Part_Int(*cell.data).begin() + part_i);
				part_i--;

			// move to ext list if particle outside of current cell
			} else if (
				real_pos[0] < cell_min[0]
				or real_pos[0] > cell_max[0]
				or real_pos[1] < cell_min[1]
				or real_pos[1] > cell_max[1]
				or real_pos[2] < cell_min[2]
				or real_pos[2] > cell_max[2]
			) {
				uint64_t destination = dccrg::error_cell;

				for (const auto& neighbor: cell.neighbors_of) {
					const auto
						neighbor_min = grid.geometry.get_min(neighbor.id),
						neighbor_max = grid.geometry.get_max(neighbor.id);

					if (
						real_pos[0] >= neighbor_min[0]
						and real_pos[0] <= neighbor_max[0]
						and real_pos[1] >= neighbor_min[1]
						and real_pos[1] <= neighbor_max[1]
						and real_pos[2] >= neighbor_min[2]
						and real_pos[2] <= neighbor_max[2]
					) {
						destination = neighbor.id;
						break;
					}
				}

				if (destination != dccrg::error_cell) {

					const auto index = Part_Ext.data(*cell.data).size();

					Part_Ext.data(*cell.data).resize(index + 1);
					assign(
						Part_Ext.data(*cell.data)[index],
						Part_Int(*cell.data)[part_i]
					);
					Part_Des(Part_Ext.data(*cell.data)[index]) = destination;

					Part_Int(*cell.data).erase(Part_Int(*cell.data).begin() + part_i);
					part_i--;

				} else {

					std::cerr << __FILE__ << "(" << __LINE__ << "): "
						<< " No destination found for particle at " << real_pos
						<< " propagated from " << real_pos
						<< " with dt " << dt
						<< " in cell " << cell.id
						<< " of length " << cell_length
						<< " at " << cell_center
						<< " with E " << JmV.data(*cell.data).cross(Vol_B.data(*cell.data))
						<< " and B " << Vol_B.data(*cell.data)
						<< " from neighbors ";
					for (const auto& neighbor: cell.neighbors_of) {
						std::cerr << neighbor.id << " ";
					}
					std::cerr << std::endl;
					abort();
				}
			}
		}

		Nr_Ext.data(*cell.data) = Part_Ext.data(*cell.data).size();
	}
	Part_Ext.type().is_stale = true;
	Nr_Ext.type().is_stale = true;
	Max_v_part.type().is_stale = true;
	Max_ω_part.type().is_stale = true;
}


template<
	class Nr_Particles_T,
	class Particles_T,
	class Cell_Iterator,
	class Grid
> void resize_receiving_containers(
	const Cell_Iterator& cells,
	Grid& grid
) {
	for (const auto& cell: cells) {
		(*cell.data)[Particles_T()].resize((*cell.data)[Nr_Particles_T()]);
	}
}


template<
	class Nr_Particles_Internal_T,
	class Particles_Internal_T,
	class Particles_External_T,
	class Particle_Destination_T,
	class Cell_Iterator,
	class Grid
> void incorporate_external_particles(
	const Cell_Iterator& cells,
	Grid& grid
) {
	constexpr Nr_Particles_Internal_T Nr_Int{};
	constexpr Particles_Internal_T Part_Int{};
	constexpr Particles_External_T Part_Ext{};
	constexpr Particle_Destination_T Dest{};

	for (const auto& cell: cells) {
		for (const auto& neighbor: cell.neighbors_of) {
			for (auto& particle: (*neighbor.data)[Part_Ext]) {
				if (particle[Dest] == dccrg::error_cell) {
					continue;
				}

				if (particle[Dest] == cell.id) {
					particle[Dest] = dccrg::error_cell;

					const auto index = (*cell.data)[Part_Int].size();

					(*cell.data)[Part_Int].resize(index + 1);

					assign((*cell.data)[Part_Int][index], particle);
				}
			}
		}

		(*cell.data)[Nr_Int] = (*cell.data)[Part_Int].size();
	}
}


template<
	class Nr_Particles_External_T,
	class Particles_External_T,
	class Cell_Iterator,
	class Grid
> void remove_external_particles(
	const Cell_Iterator& cells,
	Grid& grid
) {
	for (const auto& cell: cells) {
		(*cell.data)[Nr_Particles_External_T()] = 0;
		(*cell.data)[Particles_External_T()].clear();
	}
}


template <
	class Grid,
	class Solver_Info_Getter,
	class Timestep_Getter,
	class Max_Spatial_Velocity_Getter,
	class Max_Angular_Velocity_Getter
> void minimize_timestep(
	Grid& grid,
	const double& dt_factor,
	const Solver_Info_Getter& SInfo,
	const Timestep_Getter& Timestep,
	const Max_Spatial_Velocity_Getter& Max_v_part,
	const Max_Angular_Velocity_Getter& Max_ω_part
) try {
	using std::cerr;
	using std::endl;
	using std::max;
	using std::min;
	using std::numeric_limits;

	Timestep.type().is_stale = true;

	for (const auto& cell: grid.local_cells()) {
		if (SInfo.data(*cell.data) < 0) {
			continue;
		}

		const auto& max_v = Max_v_part.data(*cell.data);
		const auto& max_ω = Max_ω_part.data(*cell.data);
		const auto len = grid.geometry.get_length(cell.id);
		if (max_v <= 0 or max_ω <= 0) {
			Timestep.data(*cell.data) = 0;
		} else {
			auto min_dt = dt_factor * 2*M_PI/max_ω;
			for (const size_t dim: {0, 1, 2}) {
				min_dt = min(min_dt, len[dim]/max_v/2);
			}
			Timestep.data(*cell.data) = min(min_dt, Timestep.data(*cell.data));
		}
	}

} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


//! returns length of timestep taken
double timestep(
	const auto& simulation_time,
	auto& grid,
	auto& options_mhd,
	const auto& Part_Int,
	const auto& Part_Pos,
	const auto& Part_Mas,
	const auto& Part_Mas_Cell,
	const auto& Part_SpM,
	const auto& Part_SpM_Cell,
	const auto& Part_Vel,
	const auto& Part_Vel_Cell,
	const auto& Part_Ekin,
	const auto& Nr_Particles,
	const auto& Part_Nr,
	const auto& Bulk_Mass_Getter,
	const auto& Bulk_Momentum_Getter,
	const auto& Bulk_Relative_Velocity2_Getter,
	const auto& Bulk_Velocity_Getter,
	const auto& Accu_List_Number_Of_Particles_Getter,
	const auto& Accu_List_Bulk_Mass_Getter,
	const auto& Accu_List_Bulk_Velocity_Getter,
	const auto& Accu_List_Bulk_Relative_Velocity2_Getter,
	const auto& Accu_List_Target_Getter,
	const auto& Accu_List_Length_Getter,
	const auto& Accu_List_Getter,
	const auto& Nr_Accumulated_To_Cells_Getter,
	const auto& Accumulated_To_Cells_Getter,
	const auto& Bulk_Velocity,
	const auto& SInfo,
	const auto& adiabatic_index,
	const auto& vacuum_permeability,
	const auto& boltzmann,
	const auto& min_pressure,
	const auto& Mas,
	const auto& Mom,
	const auto& Nrj,
	const auto& Vol_B,
	const auto& Vol_J,
	const auto& J_m_V,
	const auto& Vol_E,
	const auto& Nr_Ext,
	const auto& Max_v_part,
	const auto& Max_ω_part,
	const auto& Part_Ext,
	const auto& Part_C2M,
	const auto& Part_Des,
	const auto& Face_dB,
	const auto& Bg_B,
	const auto& Mas_f,
	const auto& Mom_f,
	const auto& Nrj_f,
	const auto& Mag_f,
	const auto& Substep,
	const auto& Substep_Min,
	const auto& Substep_Max,
	const auto& Max_v_wave,
	const auto& Face_B,
	const auto& background_B,
	const auto& mhd_solver,
	const auto& Timestep,
	const auto& max_time_step,
	const auto& mhd_time_step_factor,
	const auto& particle_time_step_factor
) try {
	using std::cerr;
	using std::cout;
	using std::endl;
	using std::min;
	using std::numeric_limits;

	using Cell = std::remove_reference_t<decltype(grid)>::cell_data_type;

	set_minmax_substepping_period(
		simulation_time, grid, options_mhd,
		Substep_Min, Substep_Max
	);

	pamhd::mhd::minimize_timestep(
		grid, mhd_time_step_factor, SInfo, Timestep, Max_v_wave
	);

	pamhd::particle::minimize_timestep(
		grid, particle_time_step_factor, SInfo,
		Timestep, Max_v_part, Max_ω_part
	);

	const double sub_dt = set_minmax_substepping_period(
		grid, max_time_step, SInfo,
		Timestep, Substep_Min, Substep_Max
	);

	restrict_substepping_period(
		grid,
		Substep,
		Substep_Max,
		SInfo
	);

	const int max_substep = update_substeps(grid, SInfo, Substep);
	if (max_substep > 1) {
		cerr << "Substep > 1 (level > 0) not supported." << endl;
		MPI_Finalize();
		exit(EXIT_FAILURE);
	}
	if (grid.get_rank() == 0) {
		cout << "Substep: " << sub_dt
			<< ", largest substep period: " << max_substep << endl;
	}

	double total_dt = 0;
	for (int substep = 1; substep <= max_substep; substep += 1) {
		total_dt += sub_dt;

		try {
			pamhd::particle::accumulate_mhd_data(
				grid, Part_Int, Part_Pos, Part_Mas_Cell,
				Part_SpM_Cell, Part_Vel_Cell, Part_Ekin,
				Nr_Particles, Part_Nr, Bulk_Mass_Getter,
				Bulk_Momentum_Getter,
				Bulk_Relative_Velocity2_Getter,
				Bulk_Velocity_Getter,
				Accu_List_Number_Of_Particles_Getter,
				Accu_List_Bulk_Mass_Getter,
				Accu_List_Bulk_Velocity_Getter,
				Accu_List_Bulk_Relative_Velocity2_Getter,
				Accu_List_Target_Getter,
				Accu_List_Length_Getter,
				Accu_List_Getter,
				Nr_Accumulated_To_Cells_Getter,
				Accumulated_To_Cells_Getter,
				Bulk_Velocity, SInfo
			);
		} catch (const std::exception& e) {
			cerr << __FILE__ "(" << __LINE__ << ": "
				<< "Couldn't accumulate MHD data from particles: " << e.what()
				<< endl;
			abort();
		}

		try {
			pamhd::particle::fill_mhd_fluid_values(
				grid, adiabatic_index, vacuum_permeability,
				boltzmann, min_pressure, Nr_Particles,
				Bulk_Mass_Getter, Bulk_Momentum_Getter,
				Bulk_Relative_Velocity2_Getter,
				Part_Int, Mas, Mom, Nrj, Vol_B, SInfo
			);
		} catch (const std::exception& e) {
			cerr << __FILE__ "(" << __LINE__ << ": "
				<< "Couldn't fill MHD fluid values: " << e.what()
				<< endl;
			abort();
		}

		// J for E = (J - V) x B
		pamhd::math::get_curl_B(grid, Vol_B, Vol_J, SInfo);

		// E = (J - V) x B
		for (const auto& cell: grid.local_cells()) {
			J_m_V.data(*cell.data)
				= Vol_J.data(*cell.data) / vacuum_permeability
				- pamhd::mhd::get_velocity(Mom.data(*cell.data), Mas.data(*cell.data));
			// calculate electric field for output file
			Vol_E.data(*cell.data) = J_m_V.data(*cell.data).cross(Vol_B.data(*cell.data));
		}
		J_m_V.type().is_stale = true;


		// outer particles
		pamhd::particle::solve(
			sub_dt, grid.outer_cells(), grid,
			background_B, vacuum_permeability,
			true, J_m_V, Vol_B, Nr_Ext, Part_Int, Part_Ext,
			Max_v_part, Max_ω_part, Part_Pos, Part_Vel,
			Part_C2M, Part_Mas, Part_Des, SInfo
		);

		Cell::set_transfer_all(true,
			Mas.type(), Mom.type(), Nrj.type(), Vol_B.type(),
			Bg_B.type(), Substep.type(), SInfo.type(), Nr_Ext.type());
		grid.start_remote_neighbor_copy_updates();

		// inner MHD
		pamhd::mhd::get_fluxes(
			mhd_solver, grid.inner_cells(), grid, substep,
			adiabatic_index, vacuum_permeability, sub_dt,
			Mas, Mom, Nrj, Vol_B, Face_dB, Bg_B,
			Mas_f, Mom_f, Nrj_f, Mag_f,
			SInfo, Substep, Max_v_wave
		);

		// inner particles
		pamhd::particle::solve(
			sub_dt, grid.inner_cells(), grid,
			background_B, vacuum_permeability,
			true, J_m_V, Vol_B, Nr_Ext, Part_Int, Part_Ext,
			Max_v_part, Max_ω_part, Part_Pos, Part_Vel,
			Part_C2M, Part_Mas, Part_Des, SInfo
		);

		grid.wait_remote_neighbor_copy_update_receives();

		// outer MHD
		pamhd::mhd::get_fluxes(
			mhd_solver, grid.outer_cells(), grid, substep,
			adiabatic_index, vacuum_permeability, sub_dt,
			Mas, Mom, Nrj, Vol_B, Face_dB, Bg_B,
			Mas_f, Mom_f, Nrj_f, Mag_f,
			SInfo, Substep, Max_v_wave
		);

		pamhd::particle::resize_receiving_containers<
			pamhd::particle::Nr_Particles_External,
			pamhd::particle::Particles_External
		>(grid.remote_cells(), grid);

		grid.wait_remote_neighbor_copy_update_sends();
		Cell::set_transfer_all(false,
			Mas.type(), Mom.type(), Nrj.type(), Vol_B.type(),
			Bg_B.type(), Substep.type(), SInfo.type(), Nr_Ext.type());

		pamhd::mhd::get_fluxes(
			mhd_solver, grid.remote_cells(), grid, substep,
			adiabatic_index, vacuum_permeability, sub_dt,
			Mas, Mom, Nrj, Vol_B, Face_dB, Bg_B,
			Mas_f, Mom_f, Nrj_f, Mag_f,
			SInfo, Substep, Max_v_wave
		);
		Max_v_wave.type().is_stale = true;

		pamhd::mhd::update_mhd_state(
			grid.local_cells(), grid, substep,
			Face_B, Face_dB, SInfo, Substep,
			Mas, Mom, Nrj, Vol_B,
			Mas_f, Mom_f, Nrj_f, Mag_f
		);

		pamhd::mhd::update_B_consistency(
			substep, grid.local_cells(), grid,
			Mas, Mom, Nrj, Vol_B, Face_B,
			SInfo, Substep,
			adiabatic_index,
			vacuum_permeability,
			true
		);

		Cell::set_transfer_all(true, pamhd::particle::Particles_External());
		grid.start_remote_neighbor_copy_updates();

		pamhd::particle::incorporate_external_particles<
			pamhd::particle::Nr_Particles_Internal,
			pamhd::particle::Particles_Internal,
			pamhd::particle::Particles_External,
			pamhd::particle::Destination_Cell
		>(grid.inner_cells(), grid);

		grid.wait_remote_neighbor_copy_update_receives();

		pamhd::particle::incorporate_external_particles<
			pamhd::particle::Nr_Particles_Internal,
			pamhd::particle::Particles_Internal,
			pamhd::particle::Particles_External,
			pamhd::particle::Destination_Cell
		>(grid.outer_cells(), grid);

		pamhd::particle::remove_external_particles<
			pamhd::particle::Nr_Particles_External,
			pamhd::particle::Particles_External
		>(grid.inner_cells(), grid);

		grid.wait_remote_neighbor_copy_update_sends();
		Cell::set_transfer_all(false, pamhd::particle::Particles_External());

		pamhd::particle::remove_external_particles<
			pamhd::particle::Nr_Particles_External,
			pamhd::particle::Particles_External
		>(grid.outer_cells(), grid);
	}

	// update internal particles
	for (const auto& cell: grid.local_cells()) {
		// (ab)use external number counter as internal number counter
		(*cell.data)[pamhd::particle::Nr_Particles_External()]
			= (*cell.data)[pamhd::particle::Particles_Internal()].size();
	}
	Cell::set_transfer_all(true,
		pamhd::particle::Nr_Particles_External()
	);
	grid.update_copies_of_remote_neighbors();
	Cell::set_transfer_all(false,
		pamhd::particle::Nr_Particles_External()
	);

	pamhd::particle::resize_receiving_containers<
		pamhd::particle::Nr_Particles_External,
		pamhd::particle::Particles_Internal
	>(grid.remote_cells(), grid);
	Cell::set_transfer_all(true, pamhd::particle::Particles_Internal());
	grid.update_copies_of_remote_neighbors();
	Cell::set_transfer_all(false, pamhd::particle::Particles_Internal());

	for (const auto& cell: grid.local_cells()) {
		(*cell.data)[pamhd::particle::Nr_Particles_External()] = 0;
	}
	Cell::set_transfer_all(true,
		pamhd::particle::Nr_Particles_External()
	);
	grid.update_copies_of_remote_neighbors();
	Cell::set_transfer_all(false,
		pamhd::particle::Nr_Particles_External()
	);

	return total_dt;

} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


}} // namespaces

#endif // ifndef PAMHD_PARTICLE_SOLVE_DCCRG_HPP
