/*
Particle joiner of PAMHD.

Copyright 2018, 2025 Finnish Meteorological Institute
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

#ifndef PAMHD_PARTICLE_JOINER_HPP
#define PAMHD_PARTICLE_JOINER_HPP


#include "array"
#include "utility"
#include "vector"

#include "dccrg.hpp"

#include "common.hpp"


namespace pamhd {
namespace particle {


/*
Removes (TODO: Joins) particles in normal cells within given radius from origin.

Ignores cells whose center is outside of given radius.

Returns following information about removed particles in every cell from which
particles were removed: mass, momentum, temperature.

TODO: handle mix of regular and test particles.
*/
template<
	class Mass_T,
	class Velocity_T,
	class Species_Mass_T,
	class Cell,
	class Particles_Getter,
	class Solver_Info_Getter
> std::vector<
	std::tuple<
		uint64_t, // cell id
		double, // total removed mass
		std::array<double, 3>, // total removed momentum
		double // total removed temperature, 0 if < 2 removed particles
	>
> join_particles(
	const double radius,
	dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid,
	const Particles_Getter Particles,
	const Solver_Info_Getter Sol_Info,
	const unsigned int normal_cell,
	const double particle_temp_nrj_ratio
) {
	using std::get;

	std::vector<std::tuple<uint64_t, double, std::array<double, 3>, double>> ret_val;

	for (const auto& cell: grid.cells) {
		if (Sol_Info(*cell.data) != normal_cell) {
			continue;
		}

		const auto c = grid.geometry.get_center(cell.id);
		if (c[0]*c[0] + c[1]*c[1] + c[2]*c[2] > radius*radius) {
			continue;
		}

		if (Particles(*cell.data).size() == 0) {
			continue;
		}

		const auto old_size = ret_val.size();
		ret_val.resize(old_size + 1);
		get<0>(ret_val[old_size]) = cell.id;

		const auto mass
			= [&cell, &Particles](){
				double ret_val = 0;
				for (const auto& particle: Particles(*cell.data)) {
					ret_val += particle[Mass_T()];
				}
				return ret_val;
			}();
		get<1>(ret_val[old_size]) = mass;

		const auto velocity
			= get_bulk_velocity<
				Mass_T, Velocity_T, Species_Mass_T
			>(
				Particles(*cell.data)
			);
		get<2>(ret_val[old_size]) = velocity;

		const auto temperature
			= get_temperature<
				Mass_T, Velocity_T, Species_Mass_T
			>(
				Particles(*cell.data), particle_temp_nrj_ratio
			);
		get<3>(ret_val[old_size]) = temperature;
	}

	return ret_val;
}


}} // namespaces

#endif // ifndef PAMHD_PARTICLE_JOINER_HPP
