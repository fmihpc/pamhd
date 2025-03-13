/*
Particle parts of solar wind box program.

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

#ifndef PAMHD_PARTICLE_SW_BOX_HPP
#define PAMHD_PARTICLE_SW_BOX_HPP


#include "cmath"
#include "random"
#include "string"
#include "vector"

#include "dccrg.hpp"

#include "mhd/solar_wind_box.hpp"
#include "particle/common.hpp"
#include "particle/options.hpp"
#include "particle/solve_dccrg.hpp"
#include "particle/variables.hpp"
#include "simulation_options.hpp"
#include "solar_wind_box_options.hpp"
#include "variables.hpp"


namespace pamhd {
namespace particle {


//! Returns id of next particle of this process
template<
	class Grid,
	class Internal_Particle_Getter
> uint64_t initialize_plasma(
	Grid& grid,
	const double& sim_time,
	const pamhd::Options& options_sim,
	const pamhd::Solar_Wind_Box_Options& options_box,
	const particle::Options& options_part,
	std::mt19937_64& random_source,
	const Internal_Particle_Getter& Part_Int
) try {
	using Cell = std::remove_reference_t<decltype(grid)>::cell_data_type;

	if (grid.get_rank() == 0) {
		std::cout << "Initializing particles: "
			<< options_part.particles_in_cell << " #/cell" << std::endl;
	}

	uint64_t
		id_start = grid.get_rank(),
		id_increase = grid.get_comm_size();
	const auto temperature
		= options_box.sw_pressure
		/ options_box.sw_nr_density
		/ options_sim.temp2nrj;
	if (options_box.sw_dir != +1) {
		throw std::runtime_error(__FILE__":"+std::to_string(__LINE__));
	}
	const auto grid_end = grid.geometry.get_end();
	for (const auto& cell: grid.local_cells()) {
		const auto
			cell_start = grid.geometry.get_min(cell.id),
			cell_end = grid.geometry.get_max(cell.id),
			cell_length = grid.geometry.get_length(cell.id),
			cell_center = grid.geometry.get_center(cell.id);
		const auto v_factor = std::max(0.0, cell_center[0]/grid_end[0]);
		const Eigen::Vector3d velocity{
			v_factor * options_box.sw_velocity[0],
			v_factor * options_box.sw_velocity[1],
			v_factor * options_box.sw_velocity[2]};

		random_source.seed(cell.id);
		Part_Int(*cell.data) = create_particles<
			pamhd::particle::Particle_Internal,
			pamhd::particle::Mass,
			pamhd::particle::Charge_Mass_Ratio,
			pamhd::particle::Position,
			pamhd::particle::Velocity,
			pamhd::particle::Particle_ID,
			pamhd::particle::Species_Mass
		>(
			velocity,
			Eigen::Vector3d{cell_start[0], cell_start[1], cell_start[2]},
			Eigen::Vector3d{cell_end[0], cell_end[1], cell_end[2]},
			Eigen::Vector3d{temperature, temperature, temperature},
			options_part.particles_in_cell,
			options_sim.charge2mass,
			options_sim.proton_mass * options_box.sw_nr_density
				* cell_length[0] * cell_length[1] * cell_length[2],
			options_sim.proton_mass,
			options_sim.temp2nrj,
			random_source,
			id_start,
			id_increase
		);
		id_start += Part_Int(*cell.data).size() * id_increase;
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

	return id_start;

} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


template<
	class Grid,
	class Cells,
	class Internal_Particle_Getter
> uint64_t apply_solar_wind_boundaries(
	uint64_t next_particle_id,
	const uint64_t& simulation_step,
	const pamhd::Options& options_sim,
	const pamhd::Solar_Wind_Box_Options& options_box,
	const particle::Options& options_part,
	const Grid& grid,
	const Cells& solar_wind_cells,
	const double& sim_time,
	std::mt19937_64& random_source,
	const Internal_Particle_Getter& Part_Int
) try {
	const uint64_t id_increase = grid.get_comm_size();
	const auto temperature
		= options_box.sw_pressure
		/ options_box.sw_nr_density
		/ options_sim.temp2nrj;
	for (const auto& cell: solar_wind_cells) {
		if (not cell.is_local) continue;
		random_source.seed(cell.id + simulation_step * 1'000'000);

		const auto
			cell_start = grid.geometry.get_min(cell.id),
			cell_end = grid.geometry.get_max(cell.id),
			cell_length = grid.geometry.get_length(cell.id);

		Part_Int(*cell.data) = create_particles<
			pamhd::particle::Particle_Internal,
			pamhd::particle::Mass,
			pamhd::particle::Charge_Mass_Ratio,
			pamhd::particle::Position,
			pamhd::particle::Velocity,
			pamhd::particle::Particle_ID,
			pamhd::particle::Species_Mass
		>(
			Eigen::Vector3d{
				options_box.sw_velocity[0],
				options_box.sw_velocity[1],
				options_box.sw_velocity[2]},
			Eigen::Vector3d{cell_start[0], cell_start[1], cell_start[2]},
			Eigen::Vector3d{cell_end[0], cell_end[1], cell_end[2]},
			Eigen::Vector3d{temperature, temperature, temperature},
			options_part.particles_in_cell,
			options_sim.charge2mass,
			options_sim.proton_mass * options_box.sw_nr_density
				* cell_length[0] * cell_length[1] * cell_length[2],
			options_sim.proton_mass,
			options_sim.temp2nrj,
			random_source,
			next_particle_id,
			id_increase
		);
		next_particle_id += Part_Int(*cell.data).size() * id_increase;
	}
	return next_particle_id;

} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


template<
	class Grid,
	class SW_Cells,
	class Face_Cells,
	class Edge_Cells,
	class Vert_Cells,
	class Planet_Cells,
	class Internal_Particle_Getter,
	class Solver_Info_Getter
> uint64_t apply_boundaries_sw_box(
	uint64_t next_particle_id,
	const uint64_t& simulation_step,
	Grid& grid,
	const double& sim_time,
	const pamhd::Options& options_sim,
	const pamhd::Solar_Wind_Box_Options& options_box,
	const particle::Options& options_part,
	const SW_Cells& solar_wind_cells,
	const Face_Cells& face_cells,
	const Edge_Cells& edge_cells,
	const Vert_Cells& vert_cells,
	const Planet_Cells& planet_cells,
	std::mt19937_64& random_source,
	const Internal_Particle_Getter& Part_Int,
	const Solver_Info_Getter& SInfo
) try {
	using std::abs;
	using std::runtime_error;
	using std::to_string;
	using std::vector;

	using Cell = std::remove_reference_t<decltype(grid)>::cell_data_type;

	const uint64_t id_increase = grid.get_comm_size();

	// make more if not enough particles to copy
	const auto get_filler_particles = [&](
		const size_t& nr_new_particles,
		const uint64_t& cell_id,
		const size_t& next_particle_id,
		const Eigen::Vector3d& velocity
	){
		const auto
			cell_start = grid.geometry.get_min(cell_id),
			cell_end = grid.geometry.get_max(cell_id),
			cell_length = grid.geometry.get_length(cell_id);
		const auto temperature
			= options_box.sw_pressure
			/ options_box.sw_nr_density
			/ options_sim.temp2nrj;
		return create_particles<
			pamhd::particle::Particle_Internal,
			pamhd::particle::Mass,
			pamhd::particle::Charge_Mass_Ratio,
			pamhd::particle::Position,
			pamhd::particle::Velocity,
			pamhd::particle::Particle_ID,
			pamhd::particle::Species_Mass
		>(
			velocity,
			Eigen::Vector3d{cell_start[0], cell_start[1], cell_start[2]},
			Eigen::Vector3d{cell_end[0], cell_end[1], cell_end[2]},
			Eigen::Vector3d{temperature, temperature, temperature},
			nr_new_particles,
			options_sim.charge2mass,
			options_sim.proton_mass * options_box.sw_nr_density
				* cell_length[0] * cell_length[1] * cell_length[2],
			options_sim.proton_mass,
			options_sim.temp2nrj,
			random_source,
			next_particle_id,
			id_increase
		);
	};

	// boundary and normal cell share face
	for (int dir: {-3,-2,-1,+1,+2,+3}) {
		for (const auto& cell: face_cells(dir)) {
			vector<uint64_t> copy_cells{cell.id};
			for (const auto& neighbor: cell.neighbors_of) {
				const auto& fn = neighbor.face_neighbor;
				if (fn != -dir) continue;
				copy_cells.push_back(neighbor.id);
			}
			if (copy_cells.size() != 2) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");

			random_source.seed(cell.id + simulation_step * 1'000'000);
			next_particle_id += id_increase * copy_particles<
				pamhd::particle::Position,
				pamhd::particle::Particle_ID
			>(
				copy_cells, next_particle_id, id_increase,
				random_source, grid, Part_Int
			);

			// make sure boundary cell has enough particles
			auto& particles = Part_Int(*cell.data);
			if (particles.size() < options_part.particles_in_cell) {
				const auto filler = get_filler_particles(
					options_part.particles_in_cell - Part_Int(*cell.data).size(),
					cell.id, next_particle_id,
					Eigen::Vector3d{
						options_box.sw_velocity[0],
						options_box.sw_velocity[1],
						options_box.sw_velocity[2]}
				);
				next_particle_id += filler.size() * id_increase;
				particles.insert(particles.end(), filler.begin(), filler.end());
			}
		}
	}

	// boundary and normal cell share edge
	for (int dim: {0, 1, 2})
	for (int dir1: {-1, +1})
	for (int dir2: {-1, +1}) {
		for (const auto& cell: edge_cells(dim, dir1, dir2)) {
			vector<uint64_t> copy_cells{cell.id};
			for (const auto& neighbor: cell.neighbors_of) {
				const auto& en = neighbor.edge_neighbor;
				if (en[0] != dim or en[1] != -dir1 or en[2] != -dir2) continue;
				copy_cells.push_back(neighbor.id);
			}
			if (copy_cells.size() != 2) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");

			random_source.seed(cell.id + simulation_step * 1'000'000);
			next_particle_id += id_increase * copy_particles<
				pamhd::particle::Position,
				pamhd::particle::Particle_ID
			>(
				copy_cells, next_particle_id, id_increase,
				random_source, grid, Part_Int
			);

			auto& particles = Part_Int(*cell.data);
			if (particles.size() < options_part.particles_in_cell) {
				const auto filler = get_filler_particles(
					options_part.particles_in_cell - Part_Int(*cell.data).size(),
					cell.id, next_particle_id,
					Eigen::Vector3d{
						options_box.sw_velocity[0],
						options_box.sw_velocity[1],
						options_box.sw_velocity[2]}
				);
				next_particle_id += filler.size() * id_increase;
				particles.insert(particles.end(), filler.begin(), filler.end());
			}
		}
	}

	// boundary and normal cell share vertex
	for (const auto& cell: vert_cells) {
		const auto cilen = grid.mapping.get_cell_length_in_indices(cell.id);
		vector<uint64_t> copy_cells{cell.id};
		for (const auto& neighbor: cell.neighbors_of) {
			const auto& fn = neighbor.face_neighbor;
			if (fn != 0) continue;
			const auto& en = neighbor.edge_neighbor;
			if (en[0] >= 0) continue;
			if (
				abs(neighbor.x) > cilen
				or abs(neighbor.y) > cilen
				or abs(neighbor.z) > cilen
			) continue;
			copy_cells.push_back(neighbor.id);
		}
		if (copy_cells.size() != 2) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");

		random_source.seed(cell.id + simulation_step * 1'000'000);
		next_particle_id += id_increase * copy_particles<
			pamhd::particle::Position,
			pamhd::particle::Particle_ID
		>(
			copy_cells, next_particle_id, id_increase,
			random_source, grid, Part_Int
		);

		auto& particles = Part_Int(*cell.data);
		if (particles.size() < options_part.particles_in_cell) {
			const auto filler = get_filler_particles(
				options_part.particles_in_cell - Part_Int(*cell.data).size(),
				cell.id, next_particle_id,
				Eigen::Vector3d{
					options_box.sw_velocity[0],
					options_box.sw_velocity[1],
					options_box.sw_velocity[2]}
			);
			next_particle_id += filler.size() * id_increase;
			particles.insert(particles.end(), filler.begin(), filler.end());
		}
	}

	// planetary boundary cells
	const auto grid_end = grid.geometry.get_end();
	for (const auto& cell: planet_cells) {
		const auto cilen = grid.mapping.get_cell_length_in_indices(cell.id);
		vector<uint64_t> copy_cells{cell.id};
		for (const auto& neighbor: cell.neighbors_of) {
			if (
				abs(neighbor.x) > cilen
				or abs(neighbor.y) > cilen
				or abs(neighbor.z) > cilen
			) continue; // TODO: AMR
			if (SInfo.data(*neighbor.data) < 1) continue;
			copy_cells.push_back(neighbor.id);
		}
		if (copy_cells.size() < 2) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");

		random_source.seed(cell.id + simulation_step * 1'000'000);
		next_particle_id += id_increase * copy_particles<
			pamhd::particle::Position,
			pamhd::particle::Particle_ID
		>(
			copy_cells, next_particle_id, id_increase,
			random_source, grid, Part_Int
		);

		auto& particles = Part_Int(*cell.data);
		if (particles.size() < options_part.particles_in_cell) {
			const auto cell_center = grid.geometry.get_center(cell.id);
			const auto v_factor = std::max(0.0, cell_center[0]/grid_end[0]);
			const auto filler = get_filler_particles(
				options_part.particles_in_cell - Part_Int(*cell.data).size(),
				cell.id, next_particle_id,
				Eigen::Vector3d{
					v_factor * options_box.sw_velocity[0],
					v_factor * options_box.sw_velocity[1],
					v_factor * options_box.sw_velocity[2]}
			);
			next_particle_id += filler.size() * id_increase;
			particles.insert(particles.end(), filler.begin(), filler.end());
		}
	}

	next_particle_id = apply_solar_wind_boundaries(
		next_particle_id, simulation_step, options_sim,
		options_box, options_part, grid, solar_wind_cells,
		sim_time, random_source, Part_Int
	);

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

	return next_particle_id;

} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


}} // namespaces


#endif // ifndef PAMHD_PARTICLE_SW_BOX_HPP
