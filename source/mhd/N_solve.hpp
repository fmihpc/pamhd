/*
Particle-assisted version of solve.hpp.

Copyright 2014, 2015, 2016, 2017 Ilja Honkonen
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

#ifndef PAMHD_MHD_N_SOLVE_HPP
#define PAMHD_MHD_N_SOLVE_HPP


#include "cmath"
#include "limits"
#include "string"
#include "tuple"
#include "vector"

#include "dccrg.hpp"
#include "prettyprint.hpp"

#include "mhd/variables.hpp"
#include "mhd/N_rusanov.hpp"
#include "mhd/N_hll_athena.hpp"


namespace pamhd {
namespace mhd {


/*!
Advances MHD solution for one time step of length dt with given solver.

*_Getters should be a pair of objects that return a reference to given variable
of population 1 and 2 respectively when given a reference to simulation cell data. 

Returns the maximum allowed length of time step for the next step on this process.
*/
template <
	class Solver,
	class Cell,
	class Geometry,
	class Mass_Density_Getters,
	class Momentum_Density_Getters,
	class Total_Energy_Density_Getters,
	class Magnetic_Field_Getter,
	class Background_Magnetic_Field_Pos_X_Getter,
	class Background_Magnetic_Field_Pos_Y_Getter,
	class Background_Magnetic_Field_Pos_Z_Getter,
	class Mass_Density_Flux_Getters,
	class Momentum_Density_Flux_Getters,
	class Total_Energy_Density_Flux_Getters,
	class Magnetic_Field_Flux_Getter,
	class Solver_Info_Getter
> std::pair<double, size_t> N_solve(
	const Solver solver,
	const size_t solve_start_index,
	dccrg::Dccrg<Cell, Geometry>& grid,
	const double dt,
	const double adiabatic_index,
	const double vacuum_permeability,
	const Mass_Density_Getters Mas,
	const Momentum_Density_Getters Mom,
	const Total_Energy_Density_Getters Nrj,
	const Magnetic_Field_Getter Mag,
	const Background_Magnetic_Field_Pos_X_Getter Bg_B_Pos_X,
	const Background_Magnetic_Field_Pos_Y_Getter Bg_B_Pos_Y,
	const Background_Magnetic_Field_Pos_Z_Getter Bg_B_Pos_Z,
	const Mass_Density_Flux_Getters Mas_f,
	const Momentum_Density_Flux_Getters Mom_f,
	const Total_Energy_Density_Flux_Getters Nrj_f,
	const Magnetic_Field_Flux_Getter Mag_f,
	const Solver_Info_Getter Sol_Info
) {
	using std::get;
	using std::to_string;

	if (not std::isfinite(dt) or dt < 0) {
		throw std::domain_error(
			"Invalid time step: "
			+ to_string(dt)
		);
	}

	// shorthand for referring to variables of internal MHD data type
	const Mass_Density mas_int{};
	const Momentum_Density mom_int{};
	const Total_Energy_Density nrj_int{};
	const Magnetic_Field mag_int{};

	// maximum allowed next time step for cells of this process
	double max_dt = std::numeric_limits<double>::max();

	const auto& cell_data_pointers = grid.get_cell_data_pointers();

	size_t i = solve_start_index;
	for ( ; i < cell_data_pointers.size(); i++) {
		const auto& cell_id = get<0>(cell_data_pointers[i]);

		// process only inner xor outer cells
		if (cell_id == dccrg::error_cell) {
			break;
		}

		const auto& offset = get<2>(cell_data_pointers[i]);
		if (offset[0] != 0 or offset[1] != 0 or offset[2] != 0) {
			throw std::runtime_error("Unexpected neighbor cell.");
		}

		auto* const cell_data = get<1>(cell_data_pointers[i]);

		const std::array<double, 3>
			cell_length = grid.geometry.get_length(cell_id),
			// area of cell perpendicular to each dimension
			cell_area{{
				cell_length[1] * cell_length[2],
				cell_length[0] * cell_length[2],
				cell_length[0] * cell_length[1]
			}};

		i++;
		while (i < cell_data_pointers.size()) {
			const auto& neighbor_id = get<0>(cell_data_pointers[i]);

			if (neighbor_id == dccrg::error_cell) {
				i--;
				break;
			}

			const auto& neigh_offset = get<2>(cell_data_pointers[i]);

			if (neigh_offset[0] == 0 and neigh_offset[1] == 0 and neigh_offset[2] == 0) {
				i--;
				break;
			}

			// don't solve between dont_solve_cell and any other
			if ((Sol_Info(*cell_data) & pamhd::mhd::Solver_Info::dont_solve) > 0) {
				i++;
				continue;
			}

			int neighbor_dir = 0;

			// only solve between face neighbors
			if (neigh_offset[0] == 1 and neigh_offset[1] == 0 and neigh_offset[2] == 0) {
				neighbor_dir = 1;
			}
			if (neigh_offset[0] == -1 and neigh_offset[1] == 0 and neigh_offset[2] == 0) {
				neighbor_dir = -1;
			}
			if (neigh_offset[0] == 0 and neigh_offset[1] == 1 and neigh_offset[2] == 0) {
				neighbor_dir = 2;
			}
			if (neigh_offset[0] == 0 and neigh_offset[1] == -1 and neigh_offset[2] == 0) {
				neighbor_dir = -2;
			}
			if (neigh_offset[0] == 0 and neigh_offset[1] == 0 and neigh_offset[2] == 1) {
				neighbor_dir = 3;
			}
			if (neigh_offset[0] == 0 and neigh_offset[1] == 0 and neigh_offset[2] == -1) {
				neighbor_dir = -3;
			}

			if (neighbor_dir == 0) {
				i++;
				continue;
			}

			if (grid.is_local(neighbor_id) and neighbor_dir < 0) {
				/*
				This case is handled when neighbor is the current cell
				and current cell is neighbor in positive direction
				*/
				i++;
				continue;
			}

			auto* const neighbor_data = get<1>(cell_data_pointers[i]);

			if ((Sol_Info(*neighbor_data) & pamhd::mhd::Solver_Info::dont_solve) > 0) {
				i++;
				continue;
			}

			const size_t neighbor_dim = size_t(abs(neighbor_dir) - 1);

			const std::array<double, 3>
				neighbor_length = grid.geometry.get_length(neighbor_id),
				neighbor_area{{
					neighbor_length[1] * neighbor_length[2],
					neighbor_length[0] * neighbor_length[2],
					neighbor_length[0] * neighbor_length[1]
				}};

			const double shared_area
				= std::min(cell_area[neighbor_dim], neighbor_area[neighbor_dim]);

			if (not std::isnormal(shared_area) or shared_area < 0) {
				throw std::domain_error(
					"Invalid area between cells "
					+ to_string(cell_id)
					+ " and "
					+ to_string(neighbor_id)
					+ ": "
					+ to_string(shared_area)
				);
			}

			// returns total plasma state with rotated vectors for solver
			const auto get_total_state
				= [&](Cell& cell_data) {
					detail::MHD state;
					state[mas_int] = Mas.first(cell_data) + Mas.second(cell_data);
					state[mom_int]
						= get_rotated_vector(Mom.first(cell_data), abs(neighbor_dir))
						+ get_rotated_vector(Mom.second(cell_data), abs(neighbor_dir));
					state[nrj_int] = Nrj.first(cell_data) + Nrj.second(cell_data);
					state[mag_int] = get_rotated_vector(Mag(cell_data), abs(neighbor_dir));
					return state;
				};

			// take into account direction of neighbor from cell
			detail::MHD state_neg, state_pos;
			Magnetic_Field::data_type bg_face_b;
			if (neighbor_dir > 0) {
				state_neg = get_total_state(*cell_data);
				state_pos = get_total_state(*neighbor_data);

				switch (neighbor_dir) {
				case 1:
					bg_face_b = get_rotated_vector(Bg_B_Pos_X(*cell_data), 1);
					break;
				case 2:
					bg_face_b = get_rotated_vector(Bg_B_Pos_Y(*cell_data), 2);
					break;
				case 3:
					bg_face_b = get_rotated_vector(Bg_B_Pos_Z(*cell_data), 3);
					break;
				default:
					abort();
				}
			} else {
				state_pos = get_total_state(*cell_data);
				state_neg = get_total_state(*neighbor_data);

				switch (neighbor_dir) {
				case -1:
					bg_face_b = get_rotated_vector(Bg_B_Pos_X(*neighbor_data), 1);
					break;
				case -2:
					bg_face_b = get_rotated_vector(Bg_B_Pos_Y(*neighbor_data), 2);
					break;
				case -3:
					bg_face_b = get_rotated_vector(Bg_B_Pos_Z(*neighbor_data), 3);
					break;
				default:
					abort();
				}
			}

			detail::MHD flux_neg, flux_pos;
			double max_vel;
			try {
				#define SOLVER(name) \
					name< \
						pamhd::mhd::Mass_Density, \
						pamhd::mhd::Momentum_Density, \
						pamhd::mhd::Total_Energy_Density, \
						pamhd::Magnetic_Field \
					>( \
						state_neg, \
						state_pos, \
						bg_face_b, \
						shared_area, \
						dt, \
						adiabatic_index, \
						vacuum_permeability \
					)
				switch (solver) {
				case pamhd::mhd::Solver::rusanov:
					std::tie(flux_neg, flux_pos, max_vel) = SOLVER(pamhd::mhd::get_flux_N_rusanov);
					break;
				case pamhd::mhd::Solver::hll_athena:
					std::tie(flux_neg, flux_pos, max_vel) = SOLVER(pamhd::mhd::athena::get_flux_N_hll);
					break;
				default:
					abort();
				}
				#undef SOLVER
			} catch (const std::domain_error& error) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << ") "
					<< "Solution failed between cells " << cell_id
					<< " and " << neighbor_id
					<< " at " << grid.geometry.get_center(cell_id)
					<< " and " << grid.geometry.get_center(neighbor_id)
					<< " in direction " << neighbor_dir
					<< " with rotated states (mass, momentum, total energy, magnetic field): "
					<< state_neg[mas_int] << ", "
					<< state_neg[mom_int] << ", "
					<< state_neg[nrj_int] << ", "
					<< state_neg[mag_int] << " and "
					<< state_pos[mas_int] << ", "
					<< state_pos[mom_int] << ", "
					<< state_pos[nrj_int] << ", "
					<< state_pos[mag_int]
					<< " because: " << error.what()
					<< std::endl;
				abort();
			}

			max_dt = std::min(max_dt, cell_length[neighbor_dim] / max_vel);

			// rotate flux back
			flux_neg[mom_int] = get_rotated_vector(flux_neg[mom_int], -abs(neighbor_dir));
			flux_pos[mom_int] = get_rotated_vector(flux_pos[mom_int], -abs(neighbor_dir));
			flux_neg[mag_int] = get_rotated_vector(flux_neg[mag_int], -abs(neighbor_dir));
			flux_pos[mag_int] = get_rotated_vector(flux_pos[mag_int], -abs(neighbor_dir));

			// names assume neighbor is in positive direction
			const auto
				mass_frac_spec1_neg
					= Mas.first(*cell_data)
					/ (Mas.first(*cell_data) + Mas.second(*cell_data)),
				mass_frac_spec2_neg
					= Mas.second(*cell_data)
					/ (Mas.first(*cell_data) + Mas.second(*cell_data)),
				mass_frac_spec1_pos
					= Mas.first(*neighbor_data)
					/ (Mas.first(*neighbor_data) + Mas.second(*neighbor_data)),
				mass_frac_spec2_pos
					= Mas.second(*neighbor_data)
					/ (Mas.first(*neighbor_data) + Mas.second(*neighbor_data));

			if (neighbor_dir > 0) {
				Mag_f(*cell_data) -= flux_neg[mag_int] + flux_pos[mag_int];

				Mas_f.first(*cell_data)
					-= mass_frac_spec1_neg * flux_neg[mas_int]
					+ mass_frac_spec1_pos * flux_pos[mas_int];
				Mom_f.first(*cell_data)
					-= mass_frac_spec1_neg * flux_neg[mom_int]
					+ mass_frac_spec1_pos * flux_pos[mom_int];
				Nrj_f.first(*cell_data)
					-= mass_frac_spec1_neg * flux_neg[nrj_int]
					+ mass_frac_spec1_pos * flux_pos[nrj_int];

				Mas_f.second(*cell_data)
					-= mass_frac_spec2_neg * flux_neg[mas_int]
					+ mass_frac_spec2_pos * flux_pos[mas_int];
				Mom_f.second(*cell_data)
					-= mass_frac_spec2_neg * flux_neg[mom_int]
					+ mass_frac_spec2_pos * flux_pos[mom_int];
				Nrj_f.second(*cell_data)
					-= mass_frac_spec2_neg * flux_neg[nrj_int]
					+ mass_frac_spec2_pos * flux_pos[nrj_int];

				if (grid.is_local(neighbor_id)) {
					Mag_f(*neighbor_data) += flux_neg[mag_int] + flux_pos[mag_int];

					Mas_f.first(*neighbor_data)
						+= mass_frac_spec1_neg * flux_neg[mas_int]
						+ mass_frac_spec1_pos * flux_pos[mas_int];
					Mom_f.first(*neighbor_data)
						+= mass_frac_spec1_neg * flux_neg[mom_int]
						+ mass_frac_spec1_pos * flux_pos[mom_int];
					Nrj_f.first(*neighbor_data)
						+= mass_frac_spec1_neg * flux_neg[nrj_int]
						+ mass_frac_spec1_pos * flux_pos[nrj_int];

					Mas_f.second(*neighbor_data)
						+= mass_frac_spec2_neg * flux_neg[mas_int]
						+ mass_frac_spec2_pos * flux_pos[mas_int];
					Mom_f.second(*neighbor_data)
						+= mass_frac_spec2_neg * flux_neg[mom_int]
						+ mass_frac_spec2_pos * flux_pos[mom_int];
					Nrj_f.second(*neighbor_data)
						+= mass_frac_spec2_neg * flux_neg[nrj_int]
						+ mass_frac_spec2_pos * flux_pos[nrj_int];
				}

			} else {

				Mag_f(*cell_data) += flux_neg[mag_int] + flux_pos[mag_int];

				// swap fractions because neighbor in negative direction
				Mas_f.first(*cell_data)
					+= mass_frac_spec1_pos * flux_neg[mas_int]
					+ mass_frac_spec1_neg * flux_pos[mas_int];
				Mom_f.first(*cell_data)
					+= mass_frac_spec1_pos * flux_neg[mom_int]
					+ mass_frac_spec1_neg * flux_pos[mom_int];
				Nrj_f.first(*cell_data)
					+= mass_frac_spec1_pos * flux_neg[nrj_int]
					+ mass_frac_spec1_neg * flux_pos[nrj_int];

				Mas_f.second(*cell_data)
					+= mass_frac_spec2_pos * flux_neg[mas_int]
					+ mass_frac_spec2_neg * flux_pos[mas_int];
				Mom_f.second(*cell_data)
					+= mass_frac_spec2_pos * flux_neg[mom_int]
					+ mass_frac_spec2_neg * flux_pos[mom_int];
				Nrj_f.second(*cell_data)
					+= mass_frac_spec2_pos * flux_neg[nrj_int]
					+ mass_frac_spec2_neg * flux_pos[nrj_int];

				if (grid.is_local(neighbor_id)) {
					std::cerr <<  __FILE__ << "(" << __LINE__ << ") "
						"Invalid direction for adding flux to local neighbor."
						<< std::endl;
					abort();
				}
			}

			i++;
		}
	}

	return std::make_pair(max_dt, i);
}


/*!
Applies the MHD solution to given cells.

Returns 1 + last index where solution was applied.
*/
template <
	class Cell,
	class Geometry,
	class Mass_Density_Getters,
	class Momentum_Density_Getters,
	class Total_Energy_Density_Getters,
	class Magnetic_Field_Getter,
	class Mass_Density_Flux_Getters,
	class Momentum_Density_Flux_Getters,
	class Total_Energy_Density_Flux_Getters,
	class Magnetic_Field_Flux_Getter,
	class Solver_Info_Getter
> void apply_fluxes_N(
	dccrg::Dccrg<Cell, Geometry>& grid,
	const Mass_Density_Getters Mas,
	const Momentum_Density_Getters Mom,
	const Total_Energy_Density_Getters Nrj,
	const Magnetic_Field_Getter Mag,
	const Mass_Density_Flux_Getters Mas_f,
	const Momentum_Density_Flux_Getters Mom_f,
	const Total_Energy_Density_Flux_Getters Nrj_f,
	const Magnetic_Field_Flux_Getter Mag_f,
	const Solver_Info_Getter Sol_Info
) {
	for (auto& cell: grid.cells) {
		const auto length = grid.geometry.get_length(cell.id);
		const double inverse_volume = 1.0 / (length[0] * length[1] * length[2]);

		if ((Sol_Info(*cell.data) & Solver_Info::mass_density_bdy) == 0) {
			Mas.first(*cell.data) += Mas_f.first(*cell.data) * inverse_volume;
		}
		Mas_f.first(*cell.data) = 0;

		if ((Sol_Info(*cell.data) & Solver_Info::mass_density2_bdy) == 0) {
			Mas.second(*cell.data) += Mas_f.second(*cell.data) * inverse_volume;
		}
		Mas_f.second(*cell.data) = 0;


		if ((Sol_Info(*cell.data) & Solver_Info::velocity_bdy) == 0) {
			Mom.first(*cell.data) += Mom_f.first(*cell.data) * inverse_volume;
		}
		Mom_f.first(*cell.data)[0] =
		Mom_f.first(*cell.data)[1] =
		Mom_f.first(*cell.data)[2] = 0;

		if ((Sol_Info(*cell.data) & Solver_Info::velocity2_bdy) == 0) {
			Mom.second(*cell.data) += Mom_f.second(*cell.data) * inverse_volume;
		}
		Mom_f.second(*cell.data)[0] =
		Mom_f.second(*cell.data)[1] =
		Mom_f.second(*cell.data)[2] = 0;


		if ((Sol_Info(*cell.data) & Solver_Info::pressure_bdy) == 0) {
			Nrj.first(*cell.data) += Nrj_f.first(*cell.data) * inverse_volume;
		}
		Nrj_f.first(*cell.data) = 0;

		if ((Sol_Info(*cell.data) & Solver_Info::pressure2_bdy) == 0) {
			Nrj.second(*cell.data) += Nrj_f.second(*cell.data) * inverse_volume;
		}
		Nrj_f.second(*cell.data) = 0;


		if ((Sol_Info(*cell.data) & Solver_Info::magnetic_field_bdy) == 0) {
			Mag(*cell.data) += Mag_f(*cell.data) * inverse_volume;
		}
		Mag_f(*cell.data)[0] =
		Mag_f(*cell.data)[1] =
		Mag_f(*cell.data)[2] = 0;
	}
}


}} // namespaces


#endif // ifndef PAMHD_MHD_N_SOLVE_HPP
