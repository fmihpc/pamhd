/*
Particle-assisted version of solve.hpp.

Copyright 2014, 2015, 2016, 2017 Ilja Honkonen
Copyright 2019, 2023 Finnish Meteorological Institute
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
#include "type_traits"
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
	class Cell_Iterator,
	class Grid,
	class Mass_Density_Getters,
	class Momentum_Density_Getters,
	class Total_Energy_Density_Getters,
	class Magnetic_Field_Getter,
	class Background_Magnetic_Field_Getter,
	class Mass_Density_Flux_Getters,
	class Momentum_Density_Flux_Getters,
	class Total_Energy_Density_Flux_Getters,
	class Magnetic_Field_Flux_Getter,
	class Solver_Info_Getter
> double N_solve(
	const Solver solver,
	const Cell_Iterator& cells,
	Grid& grid,
	const double dt,
	const double adiabatic_index,
	const double vacuum_permeability,
	const Mass_Density_Getters Mas,
	const Momentum_Density_Getters Mom,
	const Total_Energy_Density_Getters Nrj,
	const Magnetic_Field_Getter Mag,
	const Background_Magnetic_Field_Getter Bg_B,
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

	for (const auto& cell: cells) {
		const std::array<double, 3>
			cell_length = grid.geometry.get_length(cell.id),
			// area of cell perpendicular to each dimension
			cell_area{{
				cell_length[1] * cell_length[2],
				cell_length[0] * cell_length[2],
				cell_length[0] * cell_length[1]
			}};
		const int cell_length_i = grid.mapping.get_cell_length_in_indices(cell.id);

		for (const auto& neighbor: cell.neighbors_of) {
			// don't solve between dont_solve_cell and any other
			if ((Sol_Info(*cell.data) & pamhd::mhd::Solver_Info::dont_solve) > 0) {
				continue;
			}

			// only solve between face neighbors
			const int neighbor_length_i = grid.mapping.get_cell_length_in_indices(neighbor.id);
			int overlaps = 0, direction = 0;

			if ((neighbor.x < cell_length_i) and neighbor.x > -neighbor_length_i) {
				overlaps++;
			} else if (neighbor.x == cell_length_i) {
				direction = 1;
			} else if (neighbor.x == -neighbor_length_i) {
				direction = -1;
			}

			if ((neighbor.y < cell_length_i) and neighbor.y > -neighbor_length_i) {
				overlaps++;
			} else if (neighbor.y == cell_length_i) {
				direction = 2;
			} else if (neighbor.y == -neighbor_length_i) {
				direction = -2;
			}

			if ((neighbor.z < cell_length_i) and neighbor.z > -neighbor_length_i) {
				overlaps++;
			} else if (neighbor.z == cell_length_i) {
				direction = 3;
			} else if (neighbor.z == -neighbor_length_i) {
				direction = -3;
			}

			if (overlaps < 2) {
				continue;
			}

			if (direction == 0) {
				continue;
			}

			if (grid.is_local(neighbor.id) and direction < 0) {
				continue;
			}

			if ((Sol_Info(*neighbor.data) & pamhd::mhd::Solver_Info::dont_solve) > 0) {
				continue;
			}

			const size_t neighbor_dim = size_t(abs(direction) - 1);

			const std::array<double, 3>
				neighbor_length = grid.geometry.get_length(neighbor.id),
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
					+ to_string(cell.id) + " and "
					+ to_string(neighbor.id) + ": "
					+ to_string(shared_area)
				);
			}

			// returns total plasma state with rotated vectors for solver
			const auto get_total_state
				= [&](typename Grid::cell_data_type& cell_data) {
					detail::MHD state;
					state[mas_int] = Mas.first(cell_data) + Mas.second(cell_data);
					const typename std::remove_reference<
						decltype(Mom.first(cell_data))
					>::type total_mom = Mom.first(cell_data) + Mom.second(cell_data);
					state[mom_int] = get_rotated_vector(total_mom, abs(direction));
					state[nrj_int] = Nrj.first(cell_data) + Nrj.second(cell_data);
					state[mag_int] = get_rotated_vector(Mag(cell_data), abs(direction));
					return state;
				};

			detail::MHD state_neg, state_pos;
			Magnetic_Field::data_type bg_face_b;
			// take into account direction of neighbor from cell
			if (direction > 0) {
				state_neg = get_total_state(*cell.data);
				if (state_neg[nrj_int] < 0) {
					std::cerr << __FILE__ "(" << __LINE__ << ") "
						<< "Negative energy in total state: "
						<< Nrj.first(*cell.data) << " + "
						<< Nrj.second(*cell.data) << ", "
						<< cell.id << ", "
						<< direction
						<< std::endl;
					abort();
				}

				state_pos = get_total_state(*neighbor.data);
				if (state_pos[nrj_int] < 0) {
					std::cerr << __FILE__ "(" << __LINE__ << ") "
						<< "Negative energy in total state: "
						<< Nrj.first(*neighbor.data) << " + "
						<< Nrj.second(*neighbor.data) << ", "
						<< neighbor.id << ", "
						<< direction
						<< std::endl;
					abort();
				}

				switch (direction) {
				case 1:
					bg_face_b = get_rotated_vector(Bg_B(*cell.data)(0,0), 1);
					break;
				case 2:
					bg_face_b = get_rotated_vector(Bg_B(*cell.data)(1,0), 2);
					break;
				case 3:
					bg_face_b = get_rotated_vector(Bg_B(*cell.data)(2,0), 3);
					break;
				default:
					abort();
				}
			} else {
				state_pos = get_total_state(*cell.data);
				if (state_pos[nrj_int] < 0) {
					std::cerr << __FILE__ "(" << __LINE__ << ") "
						<< "Negative energy in total state: "
						<< Nrj.first(*cell.data) << " + "
						<< Nrj.second(*cell.data) << ", "
						<< cell.id << ", "
						<< direction
						<< std::endl;
					abort();
				}

				state_neg = get_total_state(*neighbor.data);
				if (state_neg[nrj_int] < 0) {
					std::cerr << __FILE__ "(" << __LINE__ << ") "
						<< "Negative energy in total state: "
						<< Nrj.first(*neighbor.data) << " + "
						<< Nrj.second(*neighbor.data) << ", "
						<< neighbor.id << ", "
						<< direction
						<< std::endl;
					abort();
				}

				switch (direction) {
				case -1:
					bg_face_b = get_rotated_vector(Bg_B(*neighbor.data)(0,0), 1);
					break;
				case -2:
					bg_face_b = get_rotated_vector(Bg_B(*neighbor.data)(1,0), 2);
					break;
				case -3:
					bg_face_b = get_rotated_vector(Bg_B(*neighbor.data)(2,0), 3);
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
					<< "Solution failed between cells " << cell.id
					<< " and " << neighbor.id
					<< " of boundary type " << Sol_Info(*cell.data)
					<< " and " << Sol_Info(*neighbor.data)
					<< " at " << grid.geometry.get_center(cell.id)
					<< " and " << grid.geometry.get_center(neighbor.id)
					<< " in direction " << direction
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
			flux_neg[mom_int] = get_rotated_vector(flux_neg[mom_int], -abs(direction));
			flux_pos[mom_int] = get_rotated_vector(flux_pos[mom_int], -abs(direction));
			flux_neg[mag_int] = get_rotated_vector(flux_neg[mag_int], -abs(direction));
			flux_pos[mag_int] = get_rotated_vector(flux_pos[mag_int], -abs(direction));

			// names assume neighbor is in positive direction
			const auto
				mass_frac_spec1_neg
					= Mas.first(*cell.data)
					/ (Mas.first(*cell.data) + Mas.second(*cell.data)),
				mass_frac_spec2_neg
					= Mas.second(*cell.data)
					/ (Mas.first(*cell.data) + Mas.second(*cell.data)),
				mass_frac_spec1_pos
					= Mas.first(*neighbor.data)
					/ (Mas.first(*neighbor.data) + Mas.second(*neighbor.data)),
				mass_frac_spec2_pos
					= Mas.second(*neighbor.data)
					/ (Mas.first(*neighbor.data) + Mas.second(*neighbor.data));

			if (direction > 0) {
				Mag_f(*cell.data) -= flux_neg[mag_int] + flux_pos[mag_int];

				Mas_f.first(*cell.data)
					-= mass_frac_spec1_neg * flux_neg[mas_int]
					+ mass_frac_spec1_pos * flux_pos[mas_int];
				Mom_f.first(*cell.data)
					-= mass_frac_spec1_neg * flux_neg[mom_int]
					+ mass_frac_spec1_pos * flux_pos[mom_int];
				Nrj_f.first(*cell.data)
					-= mass_frac_spec1_neg * flux_neg[nrj_int]
					+ mass_frac_spec1_pos * flux_pos[nrj_int];

				Mas_f.second(*cell.data)
					-= mass_frac_spec2_neg * flux_neg[mas_int]
					+ mass_frac_spec2_pos * flux_pos[mas_int];
				Mom_f.second(*cell.data)
					-= mass_frac_spec2_neg * flux_neg[mom_int]
					+ mass_frac_spec2_pos * flux_pos[mom_int];
				Nrj_f.second(*cell.data)
					-= mass_frac_spec2_neg * flux_neg[nrj_int]
					+ mass_frac_spec2_pos * flux_pos[nrj_int];

				if (grid.is_local(neighbor.id)) {
					Mag_f(*neighbor.data) += flux_neg[mag_int] + flux_pos[mag_int];

					Mas_f.first(*neighbor.data)
						+= mass_frac_spec1_neg * flux_neg[mas_int]
						+ mass_frac_spec1_pos * flux_pos[mas_int];
					Mom_f.first(*neighbor.data)
						+= mass_frac_spec1_neg * flux_neg[mom_int]
						+ mass_frac_spec1_pos * flux_pos[mom_int];
					Nrj_f.first(*neighbor.data)
						+= mass_frac_spec1_neg * flux_neg[nrj_int]
						+ mass_frac_spec1_pos * flux_pos[nrj_int];

					Mas_f.second(*neighbor.data)
						+= mass_frac_spec2_neg * flux_neg[mas_int]
						+ mass_frac_spec2_pos * flux_pos[mas_int];
					Mom_f.second(*neighbor.data)
						+= mass_frac_spec2_neg * flux_neg[mom_int]
						+ mass_frac_spec2_pos * flux_pos[mom_int];
					Nrj_f.second(*neighbor.data)
						+= mass_frac_spec2_neg * flux_neg[nrj_int]
						+ mass_frac_spec2_pos * flux_pos[nrj_int];
				}

			} else {

				Mag_f(*cell.data) += flux_neg[mag_int] + flux_pos[mag_int];

				// swap fractions because neighbor in negative direction
				Mas_f.first(*cell.data)
					+= mass_frac_spec1_pos * flux_neg[mas_int]
					+ mass_frac_spec1_neg * flux_pos[mas_int];
				Mom_f.first(*cell.data)
					+= mass_frac_spec1_pos * flux_neg[mom_int]
					+ mass_frac_spec1_neg * flux_pos[mom_int];
				Nrj_f.first(*cell.data)
					+= mass_frac_spec1_pos * flux_neg[nrj_int]
					+ mass_frac_spec1_neg * flux_pos[nrj_int];

				Mas_f.second(*cell.data)
					+= mass_frac_spec2_pos * flux_neg[mas_int]
					+ mass_frac_spec2_neg * flux_pos[mas_int];
				Mom_f.second(*cell.data)
					+= mass_frac_spec2_pos * flux_neg[mom_int]
					+ mass_frac_spec2_neg * flux_pos[mom_int];
				Nrj_f.second(*cell.data)
					+= mass_frac_spec2_pos * flux_neg[nrj_int]
					+ mass_frac_spec2_neg * flux_pos[nrj_int];

				if (grid.is_local(neighbor.id)) {
					std::cerr <<  __FILE__ << "(" << __LINE__ << ") "
						"Invalid direction for adding flux to local neighbor."
						<< std::endl;
					abort();
				}
			}
		}
	}

	return max_dt;
}


/*!
Applies the MHD solution to given cells.

Returns 1 + last index where solution was applied.
*/
template <
	class Grid,
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
	Grid& grid,
	const double min_pressure,
	const double adiabatic_index,
	const double vacuum_permeability,
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
	for (auto& cell: grid.local_cells()) {
		if ((Sol_Info(*cell.data) & Solver_Info::dont_solve) > 0) {
			Mas_f.first(*cell.data)     =
			Mas_f.second(*cell.data)    =
			Mom_f.first(*cell.data)[0]  =
			Mom_f.first(*cell.data)[1]  =
			Mom_f.first(*cell.data)[2]  =
			Mom_f.second(*cell.data)[0] =
			Mom_f.second(*cell.data)[1] =
			Mom_f.second(*cell.data)[2] =
			Nrj_f.first(*cell.data)     =
			Nrj_f.second(*cell.data)    =
			Mag_f(*cell.data)[0]        =
			Mag_f(*cell.data)[1]        =
			Mag_f(*cell.data)[2]        = 0;
			continue;
		}

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


		/*
		Enforce minimum pressure
		*/

		const auto
			total_mass = Mas.first(*cell.data) + Mas.second(*cell.data),
			total_energy = Nrj.first(*cell.data) + Nrj.second(*cell.data);
		const typename std::remove_reference<decltype(Mom.first(*cell.data))>::type total_momentum
			= Mom.first(*cell.data) + Mom.second(*cell.data);

		const auto pressure_both = get_pressure(
			total_mass,
			total_momentum,
			total_energy,
			Mag(*cell.data),
			adiabatic_index,
			vacuum_permeability
		);
		if (pressure_both >= min_pressure) {
			continue;
		}

		// divide new energy by mass fraction
		const auto new_total_energy
			= get_total_energy_density(
				total_mass,
				get_velocity(total_momentum, total_mass),
				min_pressure,
				Mag(*cell.data),
				adiabatic_index,
				vacuum_permeability
			);
		Nrj.first(*cell.data) = Mas.first(*cell.data) / total_mass * new_total_energy;
		Nrj.second(*cell.data) = Mas.second(*cell.data) / total_mass * new_total_energy;
	}
}


}} // namespaces


#endif // ifndef PAMHD_MHD_N_SOLVE_HPP
