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
#include "divergence/remove.hpp"
#include "interpolate.hpp"
#include "mhd/solve_staggered.hpp"
#include "particle/solve.hpp"
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
	class Magnetic_Field_Getter,
	class Nr_Particles_External_Getter,
	class Particles_Internal_Getter,
	class Particles_External_Getter,
	class Particle_Position_Getter,
	class Particle_Velocity_Getter,
	class Particle_Charge_Mass_Ratio_Getter,
	class Particle_Mass_Getter,
	class Particle_Destination_Cell_Getter,
	class Solver_Info_Getter
> std::pair<double, double> solve(
	const double dt,
	const Cell_Iterator& cells,
	Grid& grid,
	const Background_Magnetic_Field& bg_B,
	const double vacuum_permeability,
	const bool E_is_derived_quantity,
	// if E_is_derived_quantity == true: JmV = J - V, else JmV = E
	const Current_Minus_Velocity_Getter JmV,
	const Magnetic_Field_Getter Mag,
	const Nr_Particles_External_Getter Nr_Ext,
	const Particles_Internal_Getter Part_Int,
	const Particles_External_Getter Part_Ext,
	const Particle_Position_Getter Part_Pos,
	const Particle_Velocity_Getter Part_Vel,
	const Particle_Charge_Mass_Ratio_Getter Part_C2M,
	const Particle_Mass_Getter Part_Mas,
	const Particle_Destination_Cell_Getter Part_Des,
	const Solver_Info_Getter SInfo
) {
	using std::abs;
	using std::isnan;
	using std::is_same;
	using std::max;
	using std::min;

	if (grid.get_maximum_refinement_level() > 0) {
		throw std::runtime_error("Only maximum refinement level 0 supported.");
	}


	std::pair<double, double> max_time_step{std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
	for (const auto& cell: cells) {
		if (SInfo(*cell.data) < 0) {
			Part_Int(*cell.data).clear();
			Part_Ext(*cell.data).clear();
			continue;
		}

		// get field data from neighborhood for interpolation
		std::array<Eigen::Vector3d, 27> current_minus_velocities, magnetic_fields;

		// default to current cell's data
		current_minus_velocities.fill(JmV(*cell.data));
		magnetic_fields.fill(Mag(*cell.data));

		for (const auto& neighbor: cell.neighbors_of) {
			if (
				abs(neighbor.x) > 1
				or abs(neighbor.y) > 1
				or abs(neighbor.z) > 1
			) {
				continue;
			}

			const size_t index = (neighbor.z + 1) * 9 + (neighbor.y + 1) * 3 + neighbor.x + 1;
			current_minus_velocities[index] = JmV(*neighbor.data);
			magnetic_fields[index] = Mag(*neighbor.data);
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

		// calculate max length of time step for next step from cell-centered values
		const Eigen::Vector3d B_centered = Mag(*cell.data) + bg_B.get_background_field(
			{cell_center[0], cell_center[1], cell_center[2]},
			vacuum_permeability
		);
		const auto E_centered = [&](){
			if (E_is_derived_quantity) {
				return JmV(*cell.data).cross(Mag(*cell.data));
			} else {
				return JmV(*cell.data);
			}
		}();

		for (size_t part_i = 0; part_i < Part_Int(*cell.data).size(); part_i++) {
			auto particle = Part_Int(*cell.data)[part_i]; // reference faster?

			// TODO check accurately only for most restrictive particle(s) in each cell?
			const auto step_size
				= get_step_size(
					1.0 / 32.0,
					// only allow particles to propagate through half a
					// cell in order to not break field interpolation
					min(cell_length[0], min(cell_length[1], cell_length[2])) / 2.0,
					Part_C2M(particle),
					Part_Vel(particle),
					E_centered,
					B_centered
				);

			max_time_step.first = min(step_size.first, max_time_step.first);
			max_time_step.second = min(step_size.second, max_time_step.second);

			// propagate for dt with substeps at most 1/32 of gyroperiod
			const int substeps
				= [&dt, &step_size](){
					const auto substeps = std::ceil(max(1.0, dt / step_size.second));
					if (substeps > std::numeric_limits<int>::max()) {
						std::cerr << __FILE__ "(" << __LINE__ << ") too many substeps" << std::endl;
						abort();
					} else {
						return int(substeps);
					}
				}();

			auto pos = Part_Pos(particle);
			auto vel = Part_Vel(particle);
			for (int substep = 0; substep < substeps; substep++) {
				const Eigen::Vector3d B_at_pos
					= interpolate(
						pos,
						interpolation_start,
						interpolation_end,
						magnetic_fields)
					+ bg_B.get_background_field(
						pos, vacuum_permeability);
				const Eigen::Matrix<double,3,1> E_at_pos = [&](){
					if (E_is_derived_quantity) {
						const auto J_m_V_at_pos = interpolate(
							pos,
							interpolation_start,
							interpolation_end,
							current_minus_velocities
						);
						return J_m_V_at_pos.cross(B_at_pos);
					} else {
						return interpolate(
							pos,
							interpolation_start,
							interpolation_end,
							current_minus_velocities
						);
					}
				}();
				std::tie(
					pos,
					vel
				) = propagate(
					pos,
					vel,
					E_at_pos,
					B_at_pos,
					Part_C2M(particle),
					dt / substeps
				);
			}
			Part_Vel(Part_Int(*cell.data)[part_i]) = vel;

			// take into account periodic grid
			const auto real_pos
				= grid.geometry.get_real_coordinate({{pos[0], pos[1], pos[2]}});

			Part_Pos(Part_Int(*cell.data)[part_i]) = {
				real_pos[0],
				real_pos[1],
				real_pos[2]
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

					const auto index = Part_Ext(*cell.data).size();

					Part_Ext(*cell.data).resize(index + 1);
					assign(
						Part_Ext(*cell.data)[index],
						Part_Int(*cell.data)[part_i]
					);
					Part_Des(Part_Ext(*cell.data)[index]) = destination;

					Part_Int(*cell.data).erase(Part_Int(*cell.data).begin() + part_i);
					part_i--;

				} else {

					std::cerr << __FILE__ << "(" << __LINE__ << "): "
						<< " No destination found for particle at " << real_pos
						<< " propagated from " << Part_Pos(particle)
						<< " with dt " << dt
						<< " in cell " << cell.id
						<< " of length " << cell_length
						<< " at " << cell_center
						<< " with E " << JmV(*cell.data).cross(Mag(*cell.data))
						<< " and B " << Mag(*cell.data)
						<< " from neighbors ";
					for (const auto& neighbor: cell.neighbors_of) {
						std::cerr << neighbor.id << " ";
					}
					std::cerr << std::endl;
					abort();
				}
			}
		}

		Nr_Ext(*cell.data) = Part_Ext(*cell.data).size();
	}

	return max_time_step;
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
	const auto& Part_Ext,
	const auto& Part_C2M,
	const auto& Part_Des,
	const auto& Face_dB,
	const auto& Bg_B,
	const auto& Mas_fs,
	const auto& Mom_fs,
	const auto& Nrj_fs,
	const auto& Mag_fs,
	const auto& Substep,
	const auto& Substep_Min,
	const auto& Substep_Max,
	const auto& Max_v,
	const auto& Face_B,
	const auto& background_B,
	const auto& mhd_solver,
	const auto& Timestep,
	const auto& max_time_step,
	const auto& time_step_factor
) try {
	using std::cerr;
	using std::cout;
	using std::endl;
	using std::min;

	using Cell = std::remove_reference_t<decltype(grid)>::cell_data_type;


	pamhd::mhd::set_minmax_substepping_period(
		simulation_time,
		grid,
		options_mhd,
		Substep_Min,
		Substep_Max
	);

	Cell::set_transfer_all(true, pamhd::mhd::Max_Velocity());
	grid.update_copies_of_remote_neighbors();
	Cell::set_transfer_all(false, pamhd::mhd::Max_Velocity());

	const double sub_dt = pamhd::mhd::set_minmax_substepping_period(
		grid,
		max_time_step,
		time_step_factor,
		SInfo,
		Timestep,
		Substep_Min,
		Substep_Max,
		Max_v
	);

	Cell::set_transfer_all(true, pamhd::mhd::Substepping_Period());
	restrict_substepping_period(
		grid,
		Substep,
		Substep_Max,
		SInfo
	);

	const int max_substep = update_substeps(grid, SInfo, Substep);
	grid.update_copies_of_remote_neighbors();
	Cell::set_transfer_all(false, pamhd::mhd::Substepping_Period());
	if (max_substep > 1) {
		std::cerr << "Substep > 1 (level > 0) not supported." << std::endl;
		MPI_Finalize();
		exit(EXIT_FAILURE);
	}
	if (grid.get_rank() == 0) {
		std::cout
			<< "Substep: " << sub_dt << ", largest substep period: "
			<< max_substep << std::endl;
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

		// B required for E calculation
		Cell::set_transfer_all(true, pamhd::mhd::MHD_State_Conservative());
		grid.start_remote_neighbor_copy_updates();

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

		// inner J for E = (J - V) x B
		pamhd::divergence::get_curl(
			grid.inner_cells(), grid, Vol_B, Vol_J, SInfo
		);
		// not included in get_curl above
		for (const auto& cell: grid.inner_cells()) {
			Vol_J(*cell.data) /= vacuum_permeability;
		}

		grid.wait_remote_neighbor_copy_update_receives();

		// outer J for E = (J - V) x B
		pamhd::divergence::get_curl(
			grid.outer_cells(), grid, Vol_B, Vol_J, SInfo
		);
		for (const auto& cell: grid.outer_cells()) {
			Vol_J(*cell.data) /= vacuum_permeability;
		}

		grid.wait_remote_neighbor_copy_update_sends();
		Cell::set_transfer_all(false, pamhd::mhd::MHD_State_Conservative());

		// inner E = (J - V) x B
		for (const auto& cell: grid.inner_cells()) {
			J_m_V(*cell.data) = Vol_J(*cell.data) - pamhd::mhd::get_velocity(Mom(*cell.data), Mas(*cell.data));
			// calculate electric field for output file
			Vol_E(*cell.data) = J_m_V(*cell.data).cross(Vol_B(*cell.data));
		}

		// outer E = (J - V) x B
		for (const auto& cell: grid.outer_cells()) {
			J_m_V(*cell.data) = Vol_J(*cell.data) - pamhd::mhd::get_velocity(Mom(*cell.data), Mas(*cell.data));
			Vol_E(*cell.data) = J_m_V(*cell.data).cross(Vol_B(*cell.data));
		}

		Cell::set_transfer_all(true, pamhd::particle::Current_Minus_Velocity());
		grid.update_copies_of_remote_neighbors();
		Cell::set_transfer_all(false, pamhd::particle::Current_Minus_Velocity());


		double
			max_dt_particle_flight = std::numeric_limits<double>::max(),
			max_dt_particle_gyro = std::numeric_limits<double>::max();

		std::pair<double, double> particle_max_dt{0, 0};

		// outer particles
		particle_max_dt = pamhd::particle::solve(
			sub_dt, grid.outer_cells(), grid,
			background_B, vacuum_permeability,
			true, J_m_V, Vol_B, Nr_Ext,
			Part_Int, Part_Ext, Part_Pos,
			Part_Vel, Part_C2M, Part_Mas,
			Part_Des, SInfo
		);

		max_dt_particle_flight = min(particle_max_dt.first, max_dt_particle_flight);
		max_dt_particle_gyro = min(particle_max_dt.second, max_dt_particle_gyro);

		Cell::set_transfer_all(
			true,
			pamhd::mhd::MHD_State_Conservative(),
			pamhd::particle::Nr_Particles_External()
		);
		grid.start_remote_neighbor_copy_updates();

		// inner MHD
		pamhd::mhd::get_fluxes(
			mhd_solver, grid.inner_cells(), grid, 1,
			adiabatic_index, vacuum_permeability,
			sub_dt, Mas, Mom, Nrj, Vol_B, Face_dB,
			Bg_B, Mas_fs, Mom_fs, Nrj_fs, Mag_fs,
			SInfo, Substep, Max_v
		);

		// inner particles
		particle_max_dt = pamhd::particle::solve(
			sub_dt, grid.inner_cells(), grid,
			background_B, vacuum_permeability,
			true, J_m_V, Vol_B, Nr_Ext,
			Part_Int, Part_Ext, Part_Pos,
			Part_Vel, Part_C2M, Part_Mas,
			Part_Des, SInfo
		);
		max_dt_particle_flight = min(particle_max_dt.first, max_dt_particle_flight);
		max_dt_particle_gyro = min(particle_max_dt.second, max_dt_particle_gyro);

		grid.wait_remote_neighbor_copy_update_receives();

		// outer MHD
		pamhd::mhd::get_fluxes(
			mhd_solver, grid.outer_cells(), grid, 1,
			adiabatic_index, vacuum_permeability,
			sub_dt, Mas, Mom, Nrj, Vol_B, Face_dB,
			Bg_B, Mas_fs, Mom_fs, Nrj_fs, Mag_fs,
			SInfo, Substep, Max_v
		);

		pamhd::particle::resize_receiving_containers<
			pamhd::particle::Nr_Particles_External,
			pamhd::particle::Particles_External
		>(grid.remote_cells(), grid);

		grid.wait_remote_neighbor_copy_update_sends();
		Cell::set_transfer_all(
			false,
			pamhd::mhd::MHD_State_Conservative(),
			pamhd::particle::Nr_Particles_External()
		);

		pamhd::mhd::update_mhd_state(
			grid.local_cells(), grid,
			substep, Face_B, Face_dB, SInfo,
			Substep, Mas, Mom, Nrj, Vol_B,
			Mas_fs, Mom_fs, Nrj_fs, Mag_fs
		);

		Cell::set_transfer_all(true,
			pamhd::Face_Magnetic_Field(),
			// update pressure for B consistency calculation
			pamhd::mhd::MHD_State_Conservative()
		);
		grid.update_copies_of_remote_neighbors();

		// constant thermal pressure when updating vol B after solution
		pamhd::mhd::update_B_consistency(
			substep, grid.local_cells(),
			Mas, Mom, Nrj, Vol_B, Face_B,
			SInfo, Substep,
			adiabatic_index,
			vacuum_permeability,
			true
		);

		grid.update_copies_of_remote_neighbors();
		Cell::set_transfer_all(false,
			pamhd::Face_Magnetic_Field(),
			pamhd::mhd::MHD_State_Conservative()
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


		Cell::set_transfer_all(true, pamhd::mhd::MHD_Flux());
		grid.update_copies_of_remote_neighbors();
		Cell::set_transfer_all(false, pamhd::mhd::MHD_Flux());

		Cell::set_transfer_all(true,
			// update pressure for B consistency calculation
			pamhd::mhd::MHD_State_Conservative()
		);
		grid.update_copies_of_remote_neighbors();
		Cell::set_transfer_all(false,
			pamhd::mhd::MHD_State_Conservative()
		);
	}

	return total_dt;

} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


}} // namespaces

#endif // ifndef PAMHD_PARTICLE_SOLVE_DCCRG_HPP
