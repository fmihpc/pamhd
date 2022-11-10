/*
Solves the MHD part of PAMHD using an external flux function.

Copyright 2014, 2015, 2016, 2017 Ilja Honkonen
Copyright 2018, 2019, 2022 Finnish Meteorological Institute
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

#ifndef PAMHD_MHD_SOLVE_STAGGERED_HPP
#define PAMHD_MHD_SOLVE_STAGGERED_HPP


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
	class Mass_Density_Flux_X_Getter,
	class Mass_Density_Flux_Y_Getter,
	class Mass_Density_Flux_Z_Getter,
	class Momentum_Density_Flux_X_Getter,
	class Momentum_Density_Flux_Y_Getter,
	class Momentum_Density_Flux_Z_Getter,
	class Total_Energy_Density_Flux_X_Getter,
	class Total_Energy_Density_Flux_Y_Getter,
	class Total_Energy_Density_Flux_Z_Getter,
	class Magnetic_Field_Flux_X_Getter,
	class Magnetic_Field_Flux_Y_Getter,
	class Magnetic_Field_Flux_Z_Getter,
	class Solver_Info_Getter
> double solve_staggered(
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
	const Mass_Density_Flux_X_Getter Mas_fx,
	const Mass_Density_Flux_Y_Getter Mas_fy,
	const Mass_Density_Flux_Z_Getter Mas_fz,
	const Momentum_Density_Flux_X_Getter Mom_fx,
	const Momentum_Density_Flux_Y_Getter Mom_fy,
	const Momentum_Density_Flux_Z_Getter Mom_fz,
	const Total_Energy_Density_Flux_X_Getter Nrj_fx,
	const Total_Energy_Density_Flux_Y_Getter Nrj_fy,
	const Total_Energy_Density_Flux_Z_Getter Nrj_fz,
	const Magnetic_Field_Flux_X_Getter Mag_fx,
	const Magnetic_Field_Flux_Y_Getter Mag_fy,
	const Magnetic_Field_Flux_Z_Getter Mag_fz,
	const Solver_Info_Getter Sol_Info
) {
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
			if ((Sol_Info(*neighbor.data) & pamhd::mhd::Solver_Info::dont_solve) > 0) {
				continue;
			}

			// only solve flux between face neighbors
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
			if (neighbor_dir <= 0) {
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
			bg_face_b = get_rotated_vector(Bg_B_Pos_X(*cell.data), neighbor_dir);

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

			if (neighbor_dir == 1) {
				Mas_fx(*cell.data) = flux[mas_int];
				Mom_fx(*cell.data) = flux[mom_int];
				Nrj_fx(*cell.data) = flux[nrj_int];
				Mag_fx(*cell.data) = flux[mag_int];
			}
			if (neighbor_dir == 2) {
				Mas_fy(*cell.data) = flux[mas_int];
				Mom_fy(*cell.data) = flux[mom_int];
				Nrj_fy(*cell.data) = flux[nrj_int];
				Mag_fy(*cell.data) = flux[mag_int];
			}
			if (neighbor_dir == 3) {
				Mas_fz(*cell.data) = flux[mas_int];
				Mom_fz(*cell.data) = flux[mom_int];
				Nrj_fz(*cell.data) = flux[nrj_int];
				Mag_fz(*cell.data) = flux[mag_int];
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
	class Face_Magnetic_Field_Getter,
	class Edge_Electric_Field_Getter,
	class Mass_Density_Flux_X_Getter,
	class Mass_Density_Flux_Y_Getter,
	class Mass_Density_Flux_Z_Getter,
	class Momentum_Density_Flux_X_Getter,
	class Momentum_Density_Flux_Y_Getter,
	class Momentum_Density_Flux_Z_Getter,
	class Total_Energy_Density_Flux_X_Getter,
	class Total_Energy_Density_Flux_Y_Getter,
	class Total_Energy_Density_Flux_Z_Getter,
	class Magnetic_Field_Flux_X_Getter,
	class Magnetic_Field_Flux_Y_Getter,
	class Magnetic_Field_Flux_Z_Getter,
	class Solver_Info_Getter
> void apply_fluxes_staggered(
	Grid& grid,
	const double dt,
	const Mass_Density_Getter Mas,
	const Momentum_Density_Getter Mom,
	const Total_Energy_Density_Getter Nrj,
	const Magnetic_Field_Getter Mag,
	const Face_Magnetic_Field_Getter Face_B,
	const Edge_Electric_Field_Getter Edge_E,
	const Mass_Density_Flux_X_Getter Mas_fx,
	const Mass_Density_Flux_Y_Getter Mas_fy,
	const Mass_Density_Flux_Z_Getter Mas_fz,
	const Momentum_Density_Flux_X_Getter Mom_fx,
	const Momentum_Density_Flux_Y_Getter Mom_fy,
	const Momentum_Density_Flux_Z_Getter Mom_fz,
	const Total_Energy_Density_Flux_X_Getter Nrj_fx,
	const Total_Energy_Density_Flux_Y_Getter Nrj_fy,
	const Total_Energy_Density_Flux_Z_Getter Nrj_fz,
	const Magnetic_Field_Flux_X_Getter Mag_fx,
	const Magnetic_Field_Flux_Y_Getter Mag_fy,
	const Magnetic_Field_Flux_Z_Getter Mag_fz,
	const Solver_Info_Getter Sol_Info
) {
	using std::to_string;

	for (const auto& cell: grid.local_cells()) {
		if ((Sol_Info(*cell.data) & Solver_Info::dont_solve) > 0) {
			continue;
		}

		const auto [dx, dy, dz] = grid.geometry.get_length(cell.id);

		if ((Sol_Info(*cell.data) & Solver_Info::mass_density_bdy) == 0) {
			Mas(*cell.data) -= dt*(Mas_fx(*cell.data)/dx + Mas_fy(*cell.data)/dy + Mas_fz(*cell.data)/dz);
		}

		if ((Sol_Info(*cell.data) & Solver_Info::velocity_bdy) == 0) {
			Mom(*cell.data) -= dt*(Mom_fx(*cell.data)/dx + Mom_fy(*cell.data)/dy + Mom_fz(*cell.data)/dz);
		}

		if ((Sol_Info(*cell.data) & Solver_Info::pressure_bdy) == 0) {
			Nrj(*cell.data) -= dt*(Nrj_fx(*cell.data)/dx + Nrj_fy(*cell.data)/dy + Nrj_fz(*cell.data)/dz);
		}

		if ((Sol_Info(*cell.data) & Solver_Info::magnetic_field_bdy) == 0) {
			Mag(*cell.data) -= dt*(Mag_fx(*cell.data)/dx + Mag_fy(*cell.data)/dy + Mag_fz(*cell.data)/dz);
		}

		auto& edge_e = Edge_E(*cell.data);
		edge_e[0] = Mag_fz(*cell.data)[1] - Mag_fy(*cell.data)[2];
		edge_e[1] = Mag_fx(*cell.data)[2] - Mag_fz(*cell.data)[0];
		edge_e[2] = Mag_fy(*cell.data)[0] - Mag_fx(*cell.data)[1];
		int e0_items = 2, e1_items = 2, e2_items = 2;

		for (const auto& neighbor: cell.neighbors_of) {
			if (neighbor.x == -1 and neighbor.y == 0 and neighbor.z == 0) {
				if ((Sol_Info(*cell.data) & Solver_Info::mass_density_bdy) == 0) {
					Mas(*cell.data) += Mas_fx(*neighbor.data)*dt/dx;
				}
				if ((Sol_Info(*cell.data) & Solver_Info::velocity_bdy) == 0) {
					Mom(*cell.data) += Mom_fx(*neighbor.data)*dt/dx;
				}
				if ((Sol_Info(*cell.data) & Solver_Info::pressure_bdy) == 0) {
					Nrj(*cell.data) += Nrj_fx(*neighbor.data)*dt/dx;
				}
				if ((Sol_Info(*cell.data) & Solver_Info::magnetic_field_bdy) == 0) {
					Mag(*cell.data) += Mag_fx(*neighbor.data)*dt/dx;
				}
			}
			if (neighbor.x == 0 and neighbor.y == -1 and neighbor.z == 0) {
				if ((Sol_Info(*cell.data) & Solver_Info::mass_density_bdy) == 0) {
					Mas(*cell.data) += Mas_fy(*neighbor.data)*dt/dy;
				}
				if ((Sol_Info(*cell.data) & Solver_Info::velocity_bdy) == 0) {
					Mom(*cell.data) += Mom_fy(*neighbor.data)*dt/dy;
				}
				if ((Sol_Info(*cell.data) & Solver_Info::pressure_bdy) == 0) {
					Nrj(*cell.data) += Nrj_fy(*neighbor.data)*dt/dy;
				}
				if ((Sol_Info(*cell.data) & Solver_Info::magnetic_field_bdy) == 0) {
					Mag(*cell.data) += Mag_fy(*neighbor.data)*dt/dy;
				}
			}
			if (neighbor.x == 0 and neighbor.y == 0 and neighbor.z == -1) {
				if ((Sol_Info(*cell.data) & Solver_Info::mass_density_bdy) == 0) {
					Mas(*cell.data) += Mas_fz(*neighbor.data)*dt/dz;
				}
				if ((Sol_Info(*cell.data) & Solver_Info::velocity_bdy) == 0) {
					Mom(*cell.data) += Mom_fz(*neighbor.data)*dt/dz;
				}
				if ((Sol_Info(*cell.data) & Solver_Info::pressure_bdy) == 0) {
					Nrj(*cell.data) += Nrj_fz(*neighbor.data)*dt/dz;
				}
				if ((Sol_Info(*cell.data) & Solver_Info::magnetic_field_bdy) == 0) {
					Mag(*cell.data) += Mag_fz(*neighbor.data)*dt/dz;
				}
			}
			if (neighbor.x == 1 and neighbor.y == 0 and neighbor.z == 0) {
				edge_e[1] -= Mag_fz(*neighbor.data)[0];
				e1_items++;
				edge_e[2] += Mag_fy(*neighbor.data)[0];
				e2_items++;
			}
			if (neighbor.x == 0 and neighbor.y == 1 and neighbor.z == 0) {
				edge_e[0] += Mag_fz(*neighbor.data)[1];
				e0_items++;
				edge_e[2] -= Mag_fx(*neighbor.data)[1];
				e2_items++;
			}
			if (neighbor.x == 0 and neighbor.y == 0 and neighbor.z == 1) {
				edge_e[0] -= Mag_fy(*neighbor.data)[2];
				e0_items++;
				edge_e[1] += Mag_fx(*neighbor.data)[2];
				e1_items++;
			}
		}
		edge_e[0] /= e0_items;
		edge_e[1] /= e1_items;
		edge_e[2] /= e2_items;
	}
	for (const auto& cell: grid.local_cells()) {
		Mas_fx(*cell.data) =
		Mas_fy(*cell.data) =
		Mas_fz(*cell.data) =
		Nrj_fx(*cell.data) =
		Nrj_fy(*cell.data) =
		Nrj_fz(*cell.data) = 0;

		Mom_fx(*cell.data) =
		Mom_fy(*cell.data) =
		Mom_fz(*cell.data) =
		Mag_fx(*cell.data) =
		Mag_fy(*cell.data) =
		Mag_fz(*cell.data) = {0, 0, 0};

		/*
		Equations 13-15 of https://doi.org/10.1006/jcph.1998.6153
		*/

		const std::array<double, 3>
			cell_length = grid.geometry.get_length(cell.id),
			cell_area{
				cell_length[1] * cell_length[2],
				cell_length[0] * cell_length[2],
				cell_length[0] * cell_length[1]
			};

		const auto& cedge_e = Edge_E(*cell.data);
		typename std::remove_reference<decltype(Face_B(*cell.data))>::type face_db{
			cell_length[2]*cedge_e[2] - cell_length[1]*cedge_e[1],
			cell_length[0]*cedge_e[0] - cell_length[2]*cedge_e[2],
			cell_length[1]*cedge_e[1] - cell_length[0]*cedge_e[0]
		};

		for (const auto& neighbor: cell.neighbors_of) {
			const auto& nedge_e = Edge_E(*neighbor.data);
			if (neighbor.x == -1 and neighbor.y == 0 and neighbor.z == 0) {
				face_db[1] += cell_length[2] * nedge_e[2];
				face_db[2] -= cell_length[1] * nedge_e[1];
			}
			if (neighbor.x == 0 and neighbor.y == -1 and neighbor.z == 0) {
				face_db[0] -= cell_length[2] * nedge_e[2];
				face_db[2] += cell_length[0] * nedge_e[0];
			}
			if (neighbor.x == 0 and neighbor.y == 0 and neighbor.z == -1) {
				face_db[0] += cell_length[1] * nedge_e[1];
				face_db[1] -= cell_length[0] * nedge_e[0];
			}
		}
		Face_B(*cell.data)[0] -= dt*face_db[0]/cell_area[0];
		Face_B(*cell.data)[1] -= dt*face_db[1]/cell_area[1];
		Face_B(*cell.data)[2] -= dt*face_db[2]/cell_area[2];
	}
}


/*!
Applies the MHD solution to normal cells of \p grid.
*/
template <
	class Solver_Info,
	class Cell_Iter,
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Volume_Magnetic_Field_Getter,
	class Face_Magnetic_Field_Getter,
	class Solver_Info_Getter
> void average_magnetic_field(
	const Cell_Iter& cells,
	const Mass_Density_Getter Mas,
	const Momentum_Density_Getter Mom,
	const Total_Energy_Density_Getter Nrj,
	const Volume_Magnetic_Field_Getter VMag,
	const Face_Magnetic_Field_Getter FMag,
	const Solver_Info_Getter Sol_Info,
	const double adiabatic_index,
	const double vacuum_permeability,
	const bool constant_thermal_pressure
) {
for (const auto& cell: cells) {
		if ((Sol_Info(*cell.data) & pamhd::mhd::Solver_Info::dont_solve) > 0) {
			continue;
		}
		if (constant_thermal_pressure and Mas(*cell.data) <= 0) {
			continue;
		}

		const auto old_pressure = [&](){
			if (constant_thermal_pressure) {
				return pamhd::mhd::get_pressure(
					Mas(*cell.data), Mom(*cell.data), Nrj(*cell.data), VMag(*cell.data),
					adiabatic_index, vacuum_permeability
				);
			} else {
				return 0.0;
			}
		}();

		// value in case neighbor isn't available
		VMag(*cell.data) = FMag(*cell.data);

		for (const auto& neighbor: cell.neighbors_of) {
			if ((Sol_Info(*neighbor.data) & pamhd::mhd::Solver_Info::dont_solve) > 0) {
				continue;
			}
			if (neighbor.x == -1 and neighbor.y == 0 and neighbor.z == 0) {
				VMag(*cell.data)[0] = 0.5 * (FMag(*cell.data)[0] + FMag(*neighbor.data)[0]);
			}
			if (neighbor.x == 0 and neighbor.y == -1 and neighbor.z == 0) {
				VMag(*cell.data)[1] = 0.5 * (FMag(*cell.data)[1] + FMag(*neighbor.data)[1]);
			}
			if (neighbor.x == 0 and neighbor.y == 0 and neighbor.z == -1) {
				VMag(*cell.data)[2] = 0.5 * (FMag(*cell.data)[2] + FMag(*neighbor.data)[2]);
			}
		}

		if (constant_thermal_pressure) {
			const auto vel = (Mom(*cell.data)/Mas(*cell.data)).eval();
			Nrj(*cell.data) = pamhd::mhd::get_total_energy_density(
				Mas(*cell.data), vel, old_pressure, VMag(*cell.data),
				adiabatic_index, vacuum_permeability
			);
		}
	}
}


}} // namespaces


#endif // ifndef PAMHD_MHD_SOLVE_STAGGERED_HPP
