/*
Particle splitter of PAMHD.

Copyright 2017 Ilja Honkonen
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

#ifndef PAMHD_PARTICLE_SPLITTER_HPP
#define PAMHD_PARTICLE_SPLITTER_HPP


#include "algorithm"
#include "random"
#include "utility"

#include "dccrg.hpp"
#include "prettyprint.hpp"

#include "common.hpp"


namespace pamhd {
namespace particle {


/*
Splits given particle and returns newly created one.

Given particle is modified in-place:
  -moved to random position within a cube fitting within
   given cell volume centered on given particle
  -particle mass is halved
Returned copy is on opposite side of original position
with respect to given particle's new position.
*/
template<
	class Particle,
	class Vector,
	class Particle_Position_Getter,
	class Particle_Mass_Getter
> Particle split(
	Particle& particle,
	const Vector& cell_min,
	const Vector& cell_max,
	std::mt19937_64& random_source,
	const Particle_Position_Getter Part_Pos,
	const Particle_Mass_Getter Part_Mas
) {
	using std::max;
	using std::min;

	const auto old_pos = Part_Pos(particle);

	// particle must be created within this distance in every
	// coord to stay within cell while preserving center of mass
	const auto max_distance =
		min(min(cell_max[0] - old_pos[0], old_pos[0] - cell_min[0]), // x
		min(min(cell_max[1] - old_pos[1], old_pos[1] - cell_min[1]), // y
		    min(cell_max[2] - old_pos[2], old_pos[2] - cell_min[2]))); // z

	std::uniform_real_distribution<double> offset_gen(-max_distance, max_distance);
	const auto
		offset_x = offset_gen(random_source),
		offset_y = offset_gen(random_source),
		offset_z = offset_gen(random_source);

	Part_Mas(particle) /= 2;
	auto new_particle = particle;

	Part_Pos(particle) = {
		old_pos[0] + offset_x,
		old_pos[1] + offset_y,
		old_pos[2] + offset_z
	};
	Part_Pos(new_particle) = {
		old_pos[0] - offset_x,
		old_pos[1] - offset_y,
		old_pos[2] - offset_z
	};

	return new_particle;
}


template<
	class Particle,
	class Cell_Data,
	class Vector,
	class Particle_Position_Getter,
	class Particle_Mass_Getter
> void split(
	std::vector<Particle>& particles,
	const size_t min_particles,
	const uint64_t cell_id,
	const Cell_Data* cell_data,
	const Vector& cell_min,
	const Vector& cell_max,
	std::mt19937_64& random_source,
	const Particle_Position_Getter Part_Pos,
	const Particle_Mass_Getter Part_Mas
) {
	using std::max;
	using std::min;
	using std::to_string;

	const size_t original_nr_particles = particles.size();
	if (min_particles > 0 and original_nr_particles == 0) {
		throw std::domain_error(
			std::string(__FILE__) + "(" + to_string(__LINE__) + "): "
			+ "No particles in cell " + to_string(cell_id) + " starting at "
			+ to_string(cell_min[0]) + ", " + to_string(cell_min[1]) + ", " + to_string(cell_min[2])
		);
	}

	if (original_nr_particles >= min_particles) {
		return;
	}

	// don't allocate in below loop
	particles.reserve(min_particles);

	std::uniform_int_distribution<size_t> index_generator(0, original_nr_particles - 1);
	for (size_t i = 0; i < min_particles - original_nr_particles; i++) {
		particles.push_back(
			split(
				particles[index_generator(random_source)],
				cell_min,
				cell_max,
				random_source,
				Part_Pos,
				Part_Mas
			)
		);
	}
}


/*
Splits random particles in normal cells until given minimum exist.

Every split preserves center of mass, other parameters also unchanged.
*/
template<
	class Grid,
	class Particles_Getter,
	class Particle_Position_Getter,
	class Particle_Mass_Getter,
	class Solver_Info_Getter
> void split_particles(
	const size_t min_particles,
	std::mt19937_64& random_source,
	Grid& grid,
	const Particles_Getter& Particles,
	const Particle_Position_Getter& Part_Pos,
	const Particle_Mass_Getter& Part_Mas,
	const Solver_Info_Getter& SInfo
) {
	using std::to_string;

	for (const auto& cell: grid.local_cells()) {
		if (SInfo.data(*cell.data) < 1) {
			continue;
		}

		const auto
			cell_min = grid.geometry.get_min(cell.id),
			cell_max = grid.geometry.get_max(cell.id);

		try {
			split(
				Particles(*cell.data),
				min_particles,
				cell.id,
				cell.data,
				cell_min,
				cell_max,
				random_source,
				Part_Pos,
				Part_Mas
			);
		} catch (const std::exception& e) {
			throw std::runtime_error(
				std::string("Couldn't split particles in normal cell ")
				+ to_string(cell.id) + ": "
			);
		}
	}
}


}} // namespaces

#endif // ifndef PAMHD_PARTICLE_SPLITTER_HPP
