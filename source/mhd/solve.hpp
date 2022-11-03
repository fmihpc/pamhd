/*
Solves the MHD part of PAMHD using an external flux function.

Copyright 2014, 2015, 2016, 2017 Ilja Honkonen
Copyright 2018, 2019 Finnish Meteorological Institute
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

#ifndef PAMHD_MHD_SOLVE_HPP
#define PAMHD_MHD_SOLVE_HPP


#include "cmath"
#include "limits"
#include "string"
#include "tuple"
#include "vector"

#include "dccrg.hpp"
#include "gensimcell.hpp"
#include "prettyprint.hpp"

#include "mhd/variables.hpp"
#include "variables.hpp"


namespace pamhd {
namespace mhd {


/*!
Advances MHD solution for one time step of length dt with given solver.

Returns the maximum allowed length of time step for the next step on this process.
*/
template <
	class Solver_Info,
	class Cell_Iter,
	class Grid,
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Magnetic_Field_Getter,
	class Background_Magnetic_Field_Pos_X_Getter,
	class Background_Magnetic_Field_Pos_Y_Getter,
	class Background_Magnetic_Field_Pos_Z_Getter,
	class Mass_Density_Flux_Getter,
	class Momentum_Density_Flux_Getter,
	class Total_Energy_Density_Flux_Getter,
	class Magnetic_Field_Flux_Getter,
	class Solver_Info_Getter
> double solve(
	const Solver solver,
	const Cell_Iter& cells,
	Grid& grid,
	const double dt,
	const double adiabatic_index,
	const double vacuum_permeability,
	const Mass_Density_Getter Mas,
	const Momentum_Density_Getter Mom,
	const Total_Energy_Density_Getter Nrj,
	const Magnetic_Field_Getter Mag,
	const Background_Magnetic_Field_Pos_X_Getter Bg_B_Pos_X,
	const Background_Magnetic_Field_Pos_Y_Getter Bg_B_Pos_Y,
	const Background_Magnetic_Field_Pos_Z_Getter Bg_B_Pos_Z,
	const Mass_Density_Flux_Getter Mas_f,
	const Momentum_Density_Flux_Getter Mom_f,
	const Total_Energy_Density_Flux_Getter Nrj_f,
	const Magnetic_Field_Flux_Getter Mag_f,
	const Solver_Info_Getter Sol_Info
) {
	using std::get;
	using std::to_string;

	if (not std::isfinite(dt) or dt < 0) {
		throw std::domain_error("Invalid time step: " + to_string(dt));
	}

	// shorthand for referring to variables of internal MHD data type
	const pamhd::mhd::Mass_Density mas_int{};
	const pamhd::mhd::Momentum_Density mom_int{};
	const pamhd::mhd::Total_Energy_Density nrj_int{};
	const pamhd::Magnetic_Field mag_int{};

	// maximum allowed next time step for cells of this process
	double max_dt = std::numeric_limits<double>::max();

	for (const auto& cell: cells) {
		if ((Sol_Info(*cell.data) & pamhd::mhd::Solver_Info::dont_solve) > 0) {
			continue;
		}

		const std::array<double, 3>
			cell_length = grid.geometry.get_length(cell.id),
			// area of cell perpendicular to each dimension
			cell_area{{
				cell_length[1] * cell_length[2],
				cell_length[0] * cell_length[2],
				cell_length[0] * cell_length[1]
			}};

		for (const auto& neighbor: cell.neighbors_of) {
			// only solve between face neighbors
			int neighbor_dir = 0;
			if (neighbor.x == 1 and neighbor.y == 0 and neighbor.z == 0) {
				neighbor_dir = 1;
			}
			if (neighbor.x == -1 and neighbor.y == 0 and neighbor.z == 0) {
				neighbor_dir = -1;
			}
			if (neighbor.x == 0 and neighbor.y == 1 and neighbor.z == 0) {
				neighbor_dir = 2;
			}
			if (neighbor.x == 0 and neighbor.y == -1 and neighbor.z == 0) {
				neighbor_dir = -2;
			}
			if (neighbor.x == 0 and neighbor.y == 0 and neighbor.z == 1) {
				neighbor_dir = 3;
			}
			if (neighbor.x == 0 and neighbor.y == 0 and neighbor.z == -1) {
				neighbor_dir = -3;
			}
			if (neighbor_dir == 0) {
				continue;
			}

			if (grid.is_local(neighbor.id) and neighbor_dir < 0) {
				/*
				This case is handled when neighbor is the current cell
				and current cell is neighbor in positive direction
				*/
				continue;
			}

			if ((Sol_Info(*neighbor.data) & pamhd::mhd::Solver_Info::dont_solve) > 0) {
				continue;
			}

			const size_t neighbor_dim = size_t(abs(neighbor_dir) - 1);

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

			detail::MHD state_neg, state_pos;
			Magnetic_Field::data_type bg_face_b;
			// take into account direction of neighbor from cell
			if (neighbor_dir > 0) {
				state_neg[mas_int] = Mas(*cell.data);
				state_neg[mom_int] = get_rotated_vector(Mom(*cell.data), abs(neighbor_dir));
				state_neg[nrj_int] = Nrj(*cell.data);
				state_neg[mag_int] = get_rotated_vector(Mag(*cell.data), abs(neighbor_dir));

				state_pos[mas_int] = Mas(*neighbor.data);
				state_pos[mom_int] = get_rotated_vector(Mom(*neighbor.data), abs(neighbor_dir));
				state_pos[nrj_int] = Nrj(*neighbor.data);
				state_pos[mag_int] = get_rotated_vector(Mag(*neighbor.data), abs(neighbor_dir));

				switch (neighbor_dir) {
				case 1:
					bg_face_b = get_rotated_vector(Bg_B_Pos_X(*cell.data), 1);
					break;
				case 2:
					bg_face_b = get_rotated_vector(Bg_B_Pos_Y(*cell.data), 2);
					break;
				case 3:
					bg_face_b = get_rotated_vector(Bg_B_Pos_Z(*cell.data), 3);
					break;
				default:
					abort();
				}
			} else {
				state_pos[mas_int] = Mas(*cell.data);
				state_pos[mom_int] = get_rotated_vector(Mom(*cell.data), abs(neighbor_dir));
				state_pos[nrj_int] = Nrj(*cell.data);
				state_pos[mag_int] = get_rotated_vector(Mag(*cell.data), abs(neighbor_dir));

				state_neg[mas_int] = Mas(*neighbor.data);
				state_neg[mom_int] = get_rotated_vector(Mom(*neighbor.data), abs(neighbor_dir));
				state_neg[nrj_int] = Nrj(*neighbor.data);
				state_neg[mag_int] = get_rotated_vector(Mag(*neighbor.data), abs(neighbor_dir));

				switch (neighbor_dir) {
				case -1:
					bg_face_b = get_rotated_vector(Bg_B_Pos_X(*neighbor.data), 1);
					break;
				case -2:
					bg_face_b = get_rotated_vector(Bg_B_Pos_Y(*neighbor.data), 2);
					break;
				case -3:
					bg_face_b = get_rotated_vector(Bg_B_Pos_Z(*neighbor.data), 3);
					break;
				default:
					abort();
				}
			}

			detail::MHD flux;
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
					std::tie(flux, max_vel) = SOLVER(pamhd::mhd::get_flux_rusanov);
					break;
				case pamhd::mhd::Solver::hll_athena:
					std::tie(flux, max_vel) = SOLVER(pamhd::mhd::athena::get_flux_hll);
					break;
				case pamhd::mhd::Solver::hlld_athena:
					std::tie(flux, max_vel) = SOLVER(pamhd::mhd::athena::get_flux_hlld);
					break;
				case pamhd::mhd::Solver::roe_athena:
					std::tie(flux, max_vel) = SOLVER(pamhd::mhd::athena::get_flux_roe);
					break;
				default:
					std::cerr <<  __FILE__ << "(" << __LINE__ << ") "
						<< "Invalid solver" << std::endl;
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
					<< " in direction " << neighbor_dir
					<< " with states (mass, momentum, total energy, magnetic field): "
					<< Mas(*cell.data) << ", "
					<< Mom(*cell.data) << ", "
					<< Nrj(*cell.data) << ", "
					<< Mag(*cell.data) << " and "
					<< Mas(*neighbor.data) << ", "
					<< Mom(*neighbor.data) << ", "
					<< Nrj(*neighbor.data) << ", "
					<< Mag(*neighbor.data)
					<< " because: " << error.what()
					<< std::endl;
				abort();
			}

			max_dt = std::min(max_dt, cell_length[neighbor_dim] / max_vel);

			// rotate flux back
			flux[mom_int] = get_rotated_vector(flux[mom_int], -abs(neighbor_dir));
			flux[mag_int] = get_rotated_vector(flux[mag_int], -abs(neighbor_dir));

			if (neighbor_dir > 0) {
				Mas_f(*cell.data) -= flux[mas_int];
				Mom_f(*cell.data) -= flux[mom_int];
				Nrj_f(*cell.data) -= flux[nrj_int];
				Mag_f(*cell.data) -= flux[mag_int];

				if (grid.is_local(neighbor.id)) {
					Mas_f(*neighbor.data) += flux[mas_int];
					Mom_f(*neighbor.data) += flux[mom_int];
					Nrj_f(*neighbor.data) += flux[nrj_int];
					Mag_f(*neighbor.data) += flux[mag_int];
				}
			} else {
				Mas_f(*cell.data) += flux[mas_int];
				Mom_f(*cell.data) += flux[mom_int];
				Nrj_f(*cell.data) += flux[nrj_int];
				Mag_f(*cell.data) += flux[mag_int];

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
Applies the MHD solution to normal cells of \p grid.
*/
template <
	class Solver_Info,
	class Grid,
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Magnetic_Field_Getter,
	class Mass_Density_Flux_Getter,
	class Momentum_Density_Flux_Getter,
	class Total_Energy_Density_Flux_Getter,
	class Magnetic_Field_Flux_Getter,
	class Solver_Info_Getter
> void apply_fluxes(
	Grid& grid,
	const double min_pressure,
	const double adiabatic_index,
	const double vacuum_permeability,
	const Mass_Density_Getter Mas,
	const Momentum_Density_Getter Mom,
	const Total_Energy_Density_Getter Nrj,
	const Magnetic_Field_Getter Mag,
	const Mass_Density_Flux_Getter Mas_f,
	const Momentum_Density_Flux_Getter Mom_f,
	const Total_Energy_Density_Flux_Getter Nrj_f,
	const Magnetic_Field_Flux_Getter Mag_f,
	const Solver_Info_Getter Sol_Info,
	const bool check_new_state = true
) {
	using std::to_string;

	for (const auto& cell: grid.local_cells()) {
		if ((Sol_Info(*cell.data) & Solver_Info::dont_solve) > 0) {
			Mas_f(*cell.data)    =
			Mom_f(*cell.data)[0] =
			Mom_f(*cell.data)[1] =
			Mom_f(*cell.data)[2] =
			Nrj_f(*cell.data)    =
			Mag_f(*cell.data)[0] =
			Mag_f(*cell.data)[1] =
			Mag_f(*cell.data)[2] = 0;
			continue;
		}

		const auto length = grid.geometry.get_length(cell.id);
		const double inverse_volume = 1.0 / (length[0] * length[1] * length[2]);

		if (check_new_state and (Sol_Info(*cell.data) & Solver_Info::mass_density_bdy) == 0) {
			Mas(*cell.data) += Mas_f(*cell.data) * inverse_volume;

			if (Mas(*cell.data) <= 0) {
				const auto c = grid.geometry.get_center(cell.id);
				throw std::domain_error(
					"New state in cell " + to_string(cell.id)
					+ " of type " + to_string(Sol_Info(*cell.data))
					+ " at (" + to_string(c[0]) + ", "
					+ to_string(c[1]) + ", " + to_string(c[2])
					+ ") has non-positive mass density: "
					+ std::to_string(Mas(*cell.data)) + " with flux "
					+ std::to_string(Mas_f(*cell.data) * inverse_volume)
				);
			}
		}
		Mas_f(*cell.data) = 0;

		if ((Sol_Info(*cell.data) & Solver_Info::velocity_bdy) == 0) {
			Mom(*cell.data) += Mom_f(*cell.data) * inverse_volume;
		}
		Mom_f(*cell.data)[0] =
		Mom_f(*cell.data)[1] =
		Mom_f(*cell.data)[2] = 0;

		if ((Sol_Info(*cell.data) & Solver_Info::magnetic_field_bdy) == 0) {
			Mag(*cell.data) += Mag_f(*cell.data) * inverse_volume;
		}
		Mag_f(*cell.data)[0] =
		Mag_f(*cell.data)[1] =
		Mag_f(*cell.data)[2] = 0;

		if ((Sol_Info(*cell.data) & Solver_Info::pressure_bdy) == 0) {
			Nrj(*cell.data) += Nrj_f(*cell.data) * inverse_volume;
		}
		Nrj_f(*cell.data) = 0;

		const auto pressure = get_pressure(
			Mas(*cell.data),
			Mom(*cell.data),
			Nrj(*cell.data),
			Mag(*cell.data),
			adiabatic_index,
			vacuum_permeability
		);
		if (pressure < min_pressure) {
			Nrj(*cell.data) = get_total_energy_density(
				Mas(*cell.data),
				get_velocity(Mom(*cell.data), Mas(*cell.data)),
				min_pressure,
				Mag(*cell.data),
				adiabatic_index,
				vacuum_permeability
			);
		}
	}
}

}} // namespaces


#endif // ifndef PAMHD_MHD_SOLVE_HPP
