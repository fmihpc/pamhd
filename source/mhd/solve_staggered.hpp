/*
Solves the MHD part of PAMHD using an external flux function.

Copyright 2014, 2015, 2016, 2017 Ilja Honkonen
Copyright 2018, 2019, 2022,
          2023, 2024, 2025 Finnish Meteorological Institute
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

#ifndef PAMHD_MHD_SOLVE_STAGGERED_HPP
#define PAMHD_MHD_SOLVE_STAGGERED_HPP


#include "array"
#include "cmath"
#include "limits"
#include "string"
#include "tuple"
#include "vector"

#include "dccrg.hpp"
#include "gensimcell.hpp"
#include "prettyprint.hpp"

#include "grid/amr.hpp"
#include "mhd/rusanov.hpp"
#include "mhd/hll_athena.hpp"
#include "mhd/hlld_athena.hpp"
#include "mhd/options.hpp"
#include "mhd/roe_athena.hpp"
#include "mhd/variables.hpp"
#include "substepping.hpp"
#include "variables.hpp"


namespace pamhd {
namespace mhd {


template <
	class Grid,
	class Cell_Iter,
	class Neighbor_Iter,
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Magnetic_Field_Getter,
	class Background_Magnetic_Field_Getter,
	class Solver_Info_Getter,
	class Substepping_Period_Getter,
	class Solver
> std::tuple<double, detail::MHD> get_flux(
	const Grid& grid,
	const Cell_Iter& cell,
	const Neighbor_Iter& neighbor,
	const int& dir,
	const Mass_Density_Getter& Mas,
	const Momentum_Density_Getter& Mom,
	const Total_Energy_Density_Getter& Nrj,
	const Magnetic_Field_Getter& Vol_B,
	const Background_Magnetic_Field_Getter& Bg_B,
	const Solver_Info_Getter& SInfo,
	const Substepping_Period_Getter& Substep,
	const Solver& solver,
	const double& adiabatic_index,
	const double& vacuum_permeability,
	const double& dt
) {
	using std::tie;
	using std::tuple;

	// shorthand for variables of internal MHD data type
	const pamhd::mhd::Mass_Density mas_int{};
	const pamhd::mhd::Momentum_Density mom_int{};
	const pamhd::mhd::Total_Energy_Density nrj_int{};
	const pamhd::Magnetic_Field mag_int{};

	detail::MHD state_neg, state_pos;
	state_neg[mas_int] = Mas.data(*cell.data);
	state_neg[mom_int] = get_rotated_vector(Mom.data(*cell.data), abs(dir));
	state_neg[nrj_int] = Nrj.data(*cell.data);
	state_neg[mag_int] = get_rotated_vector(Vol_B.data(*cell.data), abs(dir));

	state_pos[mas_int] = Mas.data(*neighbor.data);
	state_pos[mom_int] = get_rotated_vector(Mom.data(*neighbor.data), abs(dir));
	state_pos[nrj_int] = Nrj.data(*neighbor.data);
	state_pos[mag_int] = get_rotated_vector(Vol_B.data(*neighbor.data), abs(dir));

	const Magnetic_Field::data_type bg_face_b
		= get_rotated_vector(Bg_B.data(*cell.data)(dir), abs(dir));
	detail::MHD flux;
	double max_vel = -1;
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
			tie(flux, max_vel) = SOLVER(pamhd::mhd::get_flux_rusanov);
			break;
		case pamhd::mhd::Solver::hll_athena:
			tie(flux, max_vel) = SOLVER(pamhd::mhd::athena::get_flux_hll);
			break;
		case pamhd::mhd::Solver::hlld_athena:
			tie(flux, max_vel) = SOLVER(pamhd::mhd::athena::get_flux_hlld);
			break;
		case pamhd::mhd::Solver::roe_athena:
			tie(flux, max_vel) = SOLVER(pamhd::mhd::athena::get_flux_roe);
			break;
		case pamhd::mhd::Solver::hybrid:
			break;
		default:
			std::cerr <<  __FILE__ << "(" << __LINE__ << ") "
				<< "Invalid solver" << std::endl;
			abort();
		}
		#undef SOLVER
	} catch (const std::domain_error& error) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << ") "
			<< "Rank " << grid.get_rank()
			<< " solution failed with dt " << dt
			<< " between cells " << cell.id
			<< " and " << neighbor.id
			<< " of cell type " << SInfo.data(*cell.data)
			<< " and " << SInfo.data(*neighbor.data)
			<< " with substeps " << Substep.data(*cell.data)
			<< " and " << Substep.data(*neighbor.data)
			<< " at " << grid.geometry.get_center(cell.id)
			<< " and " << grid.geometry.get_center(neighbor.id)
			<< " in direction " << dir
			<< " with states (mass, momentum, total energy, magnetic field): "
			<< Mas.data(*cell.data) << ", "
			<< Mom.data(*cell.data) << ", "
			<< Nrj.data(*cell.data) << ", "
			<< Vol_B.data(*cell.data) << " and "
			<< Mas.data(*neighbor.data) << ", "
			<< Mom.data(*neighbor.data) << ", "
			<< Nrj.data(*neighbor.data) << ", "
			<< Vol_B.data(*neighbor.data)
			<< " because: " << error.what()
			<< std::endl;
		abort();
	}

	// rotate flux back
	flux[mom_int] = get_rotated_vector(flux[mom_int], -abs(dir));
	flux[mag_int] = get_rotated_vector(flux[mag_int], -abs(dir));

	return std::make_tuple(max_vel, flux);
}


/*!
Calculates MHD fluxes in/out of given cells.

Flux getters with array indices [0..5] should return variables in this order:
-x,+x,-y,+y,-z,+z.

Saves fluxes of cells with SInfo.data(*cell_data) == 1,
ignores cells with SInfo < 0.
*/
template <
	class Cell_Iter,
	class Grid,
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Volume_Magnetic_Field_Getter,
	class Face_Magnetic_Field_Getter,
	class Face_dB_Getter,
	class Background_Magnetic_Field_Getter,
	class Mass_Density_Flux_Getters,
	class Momentum_Density_Flux_Getters,
	class Total_Energy_Density_Flux_Getters,
	class Magnetic_Field_Flux_Getters,
	class Solver_Info_Getter,
	class Substepping_Period_Getter,
	class Max_Velocity_Getter
> void get_fluxes(
	const Solver& solver,
	const Cell_Iter& cells,
	Grid& grid,
	const int& current_substep,
	const double& adiabatic_index,
	const double& vacuum_permeability,
	const double& dt,
	const Mass_Density_Getter& Mas,
	const Momentum_Density_Getter& Mom,
	const Total_Energy_Density_Getter& Nrj,
	const Volume_Magnetic_Field_Getter& Vol_B,
	const Face_Magnetic_Field_Getter& Face_B,
	const Face_dB_Getter& Face_dB,
	const Background_Magnetic_Field_Getter& Bg_B,
	const Mass_Density_Flux_Getters& Mas_f,
	const Momentum_Density_Flux_Getters& Mom_f,
	const Total_Energy_Density_Flux_Getters& Nrj_f,
	const Magnetic_Field_Flux_Getters& Mag_f,
	const Solver_Info_Getter& SInfo,
	const Substepping_Period_Getter& Substep,
	const Max_Velocity_Getter& Max_v
) try {
	using std::abs;
	using std::array;
	using std::get;
	using std::min;
	using std::runtime_error;
	using std::tie;
	using std::to_string;

	const pamhd::mhd::Mass_Density mas_int{};
	const pamhd::mhd::Momentum_Density mom_int{};
	const pamhd::mhd::Total_Energy_Density nrj_int{};
	const pamhd::Magnetic_Field mag_int{};

	for (const auto& cell: cells) {
		if (cell.data == nullptr) continue;
		if (SInfo.data(*cell.data) < 0) continue;

		// skip if no local effect
		if (not cell.is_local) {
			bool skip = true;
			for (const auto& neighbor: cell.neighbors_of) {
				if (
					neighbor.face_neighbor == 0
					and neighbor.edge_neighbor[0] < 0
				) continue;
				if (neighbor.is_local) {
					skip = false;
					break;
				}
			}
			if (skip) continue;
		}

		const int csub = Substep.data(*cell.data);
		int min_min_sub = csub;

		const auto [cell_dx, cell_dy, cell_dz]
			= grid.geometry.get_length(cell.id);

		array<bool, 3> missing_flux{false, false, false};
		for (const auto& neighbor: cell.neighbors_of) {
			const auto& fn = neighbor.face_neighbor;
			if (fn == 0) continue;
			if (neighbor.data == nullptr) continue;
			if (SInfo.data(*neighbor.data) < 0) continue;

			const int nsub = Substep.data(*neighbor.data);
			min_min_sub = min(nsub, min_min_sub);
			if (
				current_substep % csub != 0
				and current_substep % nsub != 0
			) {
				continue;
			}

			if (neighbor.relative_size > 0) {
				const size_t dim = abs(fn) - 1;
				missing_flux[(dim + 1) % 3] = true;
				missing_flux[(dim + 2) % 3] = true;
			}

			// solve flux only in +dir
			if (fn < 0) continue;

			const auto [max_vel, flux] = get_flux(
				grid, cell, neighbor, fn, Mas, Mom, Nrj,
				Vol_B, Bg_B, SInfo, Substep, solver,
				adiabatic_index, vacuum_permeability, dt);

			if (solver != Solver::hybrid) {
				Max_v.data(*cell.data)(fn) = max_vel;
				Max_v.data(*neighbor.data)(-fn) = max_vel;
			}

			// cell size and substep factors for fluxes
			const auto min_sub = min(csub, nsub);
			const auto min_dt = dt * min_sub;
			double cfac = min_dt, nfac = min_dt;
			// average smaller fluxes through large face
			if (neighbor.relative_size > 0) {
				cfac /= 4;
			} else if (neighbor.relative_size < 0) {
				nfac /= 4;
			}

			const auto [neigh_dx, neigh_dy, neigh_dz]
				= grid.geometry.get_length(neighbor.id);
			const auto
				min_dx = min(cell_dx, neigh_dx),
				min_dy = min(cell_dy, neigh_dy),
				min_dz = min(cell_dz, neigh_dz);
			const auto
				Vx_avg
					= Mom.data(*cell.data)[0]/Mas.data(*cell.data)/2
					+ Mom.data(*neighbor.data)[0]/Mas.data(*neighbor.data)/2,
				Vy_avg
					= Mom.data(*cell.data)[1]/Mas.data(*cell.data)/2
					+ Mom.data(*neighbor.data)[1]/Mas.data(*neighbor.data)/2,
				Vz_avg
					= Mom.data(*cell.data)[2]/Mas.data(*cell.data)/2
					+ Mom.data(*neighbor.data)[2]/Mas.data(*neighbor.data)/2,
				Bx_avg
					= Vol_B.data(*cell.data)[0]/2
					+ Vol_B.data(*neighbor.data)[0]/2
					+ Bg_B.data(*cell.data)(-1)[0]/4
					+ Bg_B.data(*cell.data)(+1)[0]/4
					+ Bg_B.data(*neighbor.data)(-1)[0]/4
					+ Bg_B.data(*neighbor.data)(+1)[0]/4,
				By_avg
					= Vol_B.data(*cell.data)[1]/2
					+ Vol_B.data(*neighbor.data)[1]/2
					+ Bg_B.data(*cell.data)(-2)[1]/4
					+ Bg_B.data(*cell.data)(+2)[1]/4
					+ Bg_B.data(*neighbor.data)(-2)[1]/4
					+ Bg_B.data(*neighbor.data)(+2)[1]/4,
				Bz_avg
					= Vol_B.data(*cell.data)[2]/2
					+ Vol_B.data(*neighbor.data)[2]/2
					+ Bg_B.data(*cell.data)(-3)[2]/4
					+ Bg_B.data(*cell.data)(+3)[2]/4
					+ Bg_B.data(*neighbor.data)(-3)[2]/4
					+ Bg_B.data(*neighbor.data)(+3)[2]/4;
			const double inv_min_dt = [&](){
				if (min_dt == 0) return 0.0;
				else return 1.0 / min_dt;}();

			if (fn == +1) {
				if (solver != Solver::hybrid and cell.is_local) {
					Mas_f(*cell.data, +1) += cfac * flux[mas_int];
					Mom_f(*cell.data, +1) += cfac * flux[mom_int];
					Nrj_f(*cell.data, +1) += cfac * flux[nrj_int];
					Mag_f(*cell.data, +1) += cfac * flux[mag_int];
				}

				if (solver != Solver::hybrid and neighbor.is_local) {
					Mas_f(*neighbor.data, -1) += nfac * flux[mas_int];
					Mom_f(*neighbor.data, -1) += nfac * flux[mom_int];
					Nrj_f(*neighbor.data, -1) += nfac * flux[nrj_int];
					Mag_f(*neighbor.data, -1) += nfac * flux[mag_int];
				}

				const auto E_source = [&]()->array<double, 3> {
					if (solver != Solver::hybrid) {
						return {
							flux[mag_int][0],
							flux[mag_int][1],
							flux[mag_int][2]};
					} else {
						const auto Bx = [&]()->double {
							if (neighbor.relative_size < 0) {
								return Face_B.data(*neighbor.data)(-1)
									+ Bg_B.data(*neighbor.data)(-1)[0];
							} else {
								return Face_B.data(*cell.data)(+1)
									+ Bg_B.data(*cell.data)(+1)[0];
							}
						}();
						return {
							// not used in assign_face_dBs_fx
							0,
							// remove +dt*dz term applied
							// in assign_face_dBs_fx
							(Vy_avg*Bx - Vx_avg*By_avg),
								//* inv_min_dt / min_dz,
							// remove -dt*dy term
							-(Vx_avg*Bz_avg - Vz_avg*Bx)};
								//* inv_min_dt / min_dy};
					}
				}();
				assign_face_dBs_fx(
					grid, cell, neighbor, E_source, Face_dB,
					min_dy, min_dz, min_dt);
			}

			if (fn == +2) {
				if (solver != Solver::hybrid and cell.is_local) {
					Mas_f(*cell.data, +2) += cfac * flux[mas_int];
					Mom_f(*cell.data, +2) += cfac * flux[mom_int];
					Nrj_f(*cell.data, +2) += cfac * flux[nrj_int];
					Mag_f(*cell.data, +2) += cfac * flux[mag_int];
				}

				if (solver != Solver::hybrid and neighbor.is_local) {
					Mas_f(*neighbor.data, -2) += nfac * flux[mas_int];
					Mom_f(*neighbor.data, -2) += nfac * flux[mom_int];
					Nrj_f(*neighbor.data, -2) += nfac * flux[nrj_int];
					Mag_f(*neighbor.data, -2) += nfac * flux[mag_int];
				}

				const auto E_source = [&]()->array<double, 3> {
					if (solver != Solver::hybrid) {
						return {
							flux[mag_int][0],
							flux[mag_int][1],
							flux[mag_int][2]};
					} else {
						const auto By = [&]()->double {
							if (neighbor.relative_size < 0) {
								return Face_B.data(*neighbor.data)(-2)
									+ Bg_B.data(*neighbor.data)(-2)[1];
							} else {
								return Face_B.data(*cell.data)(+2)
									+ Bg_B.data(*cell.data)(+2)[1];
							}
						}();
						return {
							// remove -dt*dz use in assign...fy
							-(Vy_avg*Bx_avg - Vx_avg*By),
								//* inv_min_dt / min_dz,
							0, // not used
							// remove +dt*dx
							(Vz_avg*By - Vy_avg*Bz_avg)};
								//* inv_min_dt / min_dx};
					}
				}();
				assign_face_dBs_fy(
					grid, cell, neighbor, E_source, Face_dB,
					min_dx, min_dz, min_dt);
			}

			if (fn == +3) {
				if (solver != Solver::hybrid and cell.is_local) {
					Mas_f(*cell.data, +3) += cfac * flux[mas_int];
					Mom_f(*cell.data, +3) += cfac * flux[mom_int];
					Nrj_f(*cell.data, +3) += cfac * flux[nrj_int];
					Mag_f(*cell.data, +3) += cfac * flux[mag_int];
				}

				if (solver != Solver::hybrid and neighbor.is_local) {
					Mas_f(*neighbor.data, -3) += nfac * flux[mas_int];
					Mom_f(*neighbor.data, -3) += nfac * flux[mom_int];
					Nrj_f(*neighbor.data, -3) += nfac * flux[nrj_int];
					Mag_f(*neighbor.data, -3) += nfac * flux[mag_int];
				}

				const auto E_source = [&]()->array<double, 3> {
					if (solver != Solver::hybrid) {
						return {
							flux[mag_int][0],
							flux[mag_int][1],
							flux[mag_int][2]};
					} else {
						const auto Bz = [&]()->double {
							if (neighbor.relative_size < 0) {
								return Face_B.data(*neighbor.data)(-3)
									+ Bg_B.data(*neighbor.data)(-3)[2];
							} else {
								return Face_B.data(*cell.data)(+3)
									+ Bg_B.data(*cell.data)(+3)[2];
							}
						}();
						return {
							// remove +dt*dy term
							(Vx_avg*Bz - Vz_avg*Bx_avg),
								//* inv_min_dt / min_dy,
							// remove -dt*dx
							-(Vz_avg*By_avg - Vy_avg*Bz),
								//* inv_min_dt / min_dx,
							0};
					}
				}();
				assign_face_dBs_fz(
					grid, cell, neighbor, E_source, Face_dB,
					min_dx, min_dy, min_dt);
			}
		}

		if (missing_flux[0]) {
			const auto [max_vel, flux] = get_flux(
				grid, cell, cell, +1, Mas, Mom, Nrj, Vol_B, Bg_B,
				SInfo, Substep, solver, adiabatic_index,
				vacuum_permeability, dt);

			assign_missing_face_dBs_fx(
				grid, cell, flux[mag_int], Face_dB,
				cell_dy, cell_dz, dt*min_min_sub);
		}
		if (missing_flux[1]) {
			const auto [max_vel, flux] = get_flux(
				grid, cell, cell, +2, Mas, Mom, Nrj, Vol_B, Bg_B,
				SInfo, Substep, solver, adiabatic_index,
				vacuum_permeability, dt);

			assign_missing_face_dBs_fy(
				grid, cell, flux[mag_int], Face_dB,
				cell_dx, cell_dz, dt*min_min_sub);
		}
		if (missing_flux[2]) {
			const auto [max_vel, flux] = get_flux(
				grid, cell, cell, +3, Mas, Mom, Nrj, Vol_B, Bg_B,
				SInfo, Substep, solver, adiabatic_index,
				vacuum_permeability, dt);

			assign_missing_face_dBs_fz(
				grid, cell, flux[mag_int], Face_dB,
				cell_dx, cell_dy, dt*min_min_sub);
		}
	}
	Max_v.type().is_stale = true;

} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


/*
Edge directed along x:
	Flux through touching face with normal in y direction:
		+Bz flux is added to E
	Flux through face with normal in z direction:
		-By flux is added to E
Edge y:
	Flux x: -Bz flux added
	Flux z: +Bx flux added
Edge z:
	Flux x: +By
	Flux y: -Bx
*/

// changes in face magnetic fields from x-directed flux
template <
	class Grid,
	class Cell_Iter,
	class Neighbor_Iter,
	class Flux_Vec,
	class Face_dB_Getter
> void assign_face_dBs_fx(
	const Grid& grid,
	const Cell_Iter& cell, // flux out of this cell
	const Neighbor_Iter& flux_neigh, // flux into this cell
	const Flux_Vec& flux_mag, // only magnetic part of flux needed
	const Face_dB_Getter& Face_dB,
	const double& dy,
	const double& dz,
	const double& dt
) {
	using std::runtime_error;
	using std::to_string;

	if (flux_neigh.face_neighbor != +1) {
		throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
	}

	const double
		// average E from 4 edge-touching faces
		Ey = -flux_mag[2] * dt * dy / 4,
		Ez = +flux_mag[1] * dt * dz / 4;

	const auto
		cleni = grid.mapping.get_cell_length_in_indices(cell.id),
		flux_nleni = grid.mapping.get_cell_length_in_indices(flux_neigh.id);

	auto& cfdb = Face_dB.data(*cell.data);
	{ // handle cell and flux neighbor separately
		auto
			&nfdb = Face_dB.data(*flux_neigh.data);
		// cfdb(+1) and nfdb(-1) cancel out on shared face
		// cfdb(-1) and nfdb(+1) not affected

		// nfdb(-1) is affected if neighbor larger than cell
		if (flux_neigh.relative_size < 0) {
			if (flux_neigh.y == 0) {
				nfdb(-1) -= Ez;
			} else if (cleni == flux_neigh.y + flux_nleni) {
				nfdb(-1) += Ez;
			} else {
				throw runtime_error(__FILE__ "("
					+ to_string(__LINE__) + ")");
			}
			if (flux_neigh.z == 0) {
				nfdb(-1) += Ey;
			} else if (cleni == flux_neigh.z + flux_nleni) {
				nfdb(-1) -= Ey;
			} else {
				throw runtime_error(__FILE__ "("
					+ to_string(__LINE__) + ")");
			}
		}
		// no effect if neighbor smaller and -2 faces don't align
		if (flux_neigh.relative_size <= 0 or flux_neigh.y == 0) {
			cfdb(-2) -= Ez;
		}
		// no effect if cell smaller and -2 faces don't align
		if (flux_neigh.relative_size >= 0 or flux_neigh.y == 0) {
			nfdb(-2) += Ez;
		}
		if (flux_neigh.relative_size <= 0 or cleni == flux_neigh.y + flux_nleni) {
			cfdb(+2) -= Ez;
		}
		if (flux_neigh.relative_size >= 0 or cleni == flux_neigh.y + flux_nleni) {
			nfdb(+2) += Ez;
		}
		if (flux_neigh.relative_size <= 0 or flux_neigh.z == 0) {
			cfdb(-3) += Ey;
		}
		if (flux_neigh.relative_size >= 0 or flux_neigh.z == 0) {
			nfdb(-3) -= Ey;
		}
		if (flux_neigh.relative_size <= 0 or cleni == flux_neigh.z + flux_nleni) {
			cfdb(+3) += Ey;
		}
		if (flux_neigh.relative_size >= 0 or cleni == flux_neigh.z + flux_nleni) {
			nfdb(+3) -= Ey;
		}
	}

	// assumes cells of same ref lvl are of same physical size
	for (const auto& neighbor: cell.neighbors_of) {
		if (neighbor.data == nullptr) continue;
		//if (neighbor.id == flux_neigh.id) continue;

		const auto& fn = neighbor.face_neighbor;
		const auto& en = neighbor.edge_neighbor;
		if (fn == 0 and en[0] < 0) continue;

		const auto nleni = grid.mapping.get_cell_length_in_indices(neighbor.id);
		auto& nfdb = Face_dB.data(*neighbor.data);

		// flux neighbor handled above
		if (fn == +1 and neighbor.id != flux_neigh.id) {
			if (neighbor.y != flux_neigh.y and neighbor.z != flux_neigh.z) {
				// no effect
			} else if (neighbor.y < flux_neigh.y) {
				// corresponding change to larger cell
				cfdb(+1) += Ez;
				nfdb(-1) += Ez;
				nfdb(+2) += Ez;
			} else if (neighbor.y > flux_neigh.y) {
				cfdb(+1) -= Ez;
				nfdb(-1) -= Ez;
				nfdb(-2) += Ez;
			} else if (neighbor.z < flux_neigh.z) {
				cfdb(+1) -= Ey;
				nfdb(-1) -= Ey;
				nfdb(+3) -= Ey;
			} else if (neighbor.z > flux_neigh.z) {
				cfdb(+1) += Ey;
				nfdb(-1) += Ey;
				nfdb(-3) -= Ey;
			} else {
				throw std::runtime_error(__FILE__ "("
					+ std::to_string(__LINE__) + ")");
			}
		}

		// neighbor edge might be shorter than given dx/dy/dz
		const auto [nEy, nEz] = [&](){
			if (
				neighbor.relative_size > 0
				and flux_neigh.relative_size == 0
			) {
				return std::make_tuple(Ey/2, Ez/2);
			} else {
				return std::make_tuple(Ey, Ez);
			}
		}();

		// no effect if neighbor's and cell's +1 faces don't align
		if (cleni == neighbor.x + nleni) {
			if (fn == -2
				// no effect if smaller flux neighbor...
				and (flux_neigh.relative_size <= 0
					// ...doesn't touch cell's -2 face...
					or (flux_neigh.y == 0
						// or isn't aligned with cell's smaller neighbor
						and (neighbor.relative_size <= 0
							or neighbor.z == flux_neigh.z)))
			) {
				nfdb(+1) += nEz;
				nfdb(+2) -= nEz;
			}
			if (fn == +2
				and (flux_neigh.relative_size <= 0
					or (cleni == flux_neigh.y + flux_nleni
						and (neighbor.relative_size <= 0
							or neighbor.z == flux_neigh.z)))
			) {
				nfdb(+1) -= nEz;
				nfdb(-2) -= nEz;
			}
			if (fn == -3
				and (flux_neigh.relative_size <= 0
					or (flux_neigh.z == 0
						and (neighbor.relative_size <= 0
							or neighbor.y == flux_neigh.y)))
			) {
				nfdb(+1) -= nEy;
				nfdb(+3) += nEy;
			}
			if (fn == +3
				and (flux_neigh.relative_size <= 0
					or (cleni == flux_neigh.z + flux_nleni
						and (neighbor.relative_size <= 0
							or neighbor.y == flux_neigh.y)))
			) {
				nfdb(+1) += nEy;
				nfdb(-3) += nEy;
			}
		}
		if (
			en[0] == 1
			and en[1] == +1
			and en[2] == -1
			// no effect if flux neighbor doesn't touch this neighbor
			and flux_neigh.z == 0
			and (neighbor.relative_size <= 0
				or flux_neigh.relative_size <= 0
				// nor if smaller neighbors not aligned
				or neighbor.y == flux_neigh.y)
		) {
			nfdb(-1) -= nEy;
			nfdb(+3) -= nEy;
		}
		if (
			en[0] == 1
			and en[1] == +1
			and en[2] == +1
			and cleni == flux_neigh.z + flux_nleni
			and (neighbor.relative_size <= 0
				or flux_neigh.relative_size <= 0
				or neighbor.y == flux_neigh.y)
		) {
			nfdb(-1) += nEy;
			nfdb(-3) -= nEy;
		}
		if (
			en[0] == 2
			and en[1] == +1
			and en[2] == -1
			and flux_neigh.y == 0
			and (neighbor.relative_size <= 0
				or flux_neigh.relative_size <= 0
				or neighbor.z == flux_neigh.z)
		) {
			nfdb(-1) += nEz;
			nfdb(+2) += nEz;
		}
		if (
			en[0] == 2
			and en[1] == +1
			and en[2] == +1
			and cleni == flux_neigh.y + flux_nleni
			and (neighbor.relative_size <= 0
				or flux_neigh.relative_size <= 0
				or neighbor.z == flux_neigh.z)
		) {
			nfdb(-1) -= nEz;
			nfdb(-2) += nEz;
		}
	}
}

template <
	class Grid,
	class Cell_Iter,
	class Flux_Vec,
	class Face_dB_Getter
> void assign_missing_face_dBs_fx(
	const Grid& grid,
	const Cell_Iter& cell,
	const Flux_Vec& flux_mag,
	const Face_dB_Getter& Face_dB,
	double dy,
	double dz,
	const double& dt
) {
	using std::runtime_error;
	using std::to_string;

	// only smaller neighbors affected
	const double
		nEy = -flux_mag[2] * dt * dy / 8,
		nEz = +flux_mag[1] * dt * dz / 8;

	const auto cleni = grid.mapping.get_cell_length_in_indices(cell.id);

	for (const auto& neighbor: cell.neighbors_of) {
		if (neighbor.data == nullptr) continue;
		if (neighbor.relative_size <= 0) continue;

		const auto& fn = neighbor.face_neighbor;
		if (fn == 0 or abs(fn) == 1) continue;

		const auto nleni = grid.mapping.get_cell_length_in_indices(neighbor.id);
		auto& nfdb = Face_dB.data(*neighbor.data);

		if (fn == -2) {
			if (neighbor.x == 0) {
				nfdb(+1) += nEz;
				nfdb(+2) -= nEz;
			} else if (cleni == neighbor.x + nleni) {
				nfdb(-1) += nEz;
				nfdb(+2) += nEz;
			} else {
				throw runtime_error(__FILE__ "("
					+ to_string(__LINE__) + ")");
			}
		}
		if (fn == +2) {
			if (neighbor.x == 0) {
				nfdb(+1) -= nEz;
				nfdb(-2) -= nEz;
			} else if (cleni == neighbor.x + nleni) {
				nfdb(-1) -= nEz;
				nfdb(-2) += nEz;
			} else {
				throw runtime_error(__FILE__ "("
					+ to_string(__LINE__) + ")");
			}
		}
		if (fn == -3) {
			if (neighbor.x == 0) {
				nfdb(+1) -= nEy;
				nfdb(+3) += nEy;
			} else if (cleni == neighbor.x + nleni) {
				nfdb(-1) -= nEy;
				nfdb(+3) -= nEy;
			} else {
				throw runtime_error(__FILE__ "("
					+ to_string(__LINE__) + ")");
			}
		}
		if (fn == +3) {
			if (neighbor.x == 0) {
				nfdb(+1) += nEy;
				nfdb(-3) += nEy;
			} else if (cleni == neighbor.x + nleni) {
				nfdb(-1) += nEy;
				nfdb(-3) -= nEy;
			} else {
				throw runtime_error(__FILE__ "("
					+ to_string(__LINE__) + ")");
			}
		}
	}
}

// changes in face magnetic fields from y-directed flux
template <
	class Grid,
	class Cell_Iter,
	class Neighbor_Iter,
	class Flux_Vec,
	class Face_dB_Getter
> void assign_face_dBs_fy(
	const Grid& grid,
	const Cell_Iter& cell,
	const Neighbor_Iter& flux_neigh,
	const Flux_Vec& flux_mag,
	const Face_dB_Getter& Face_dB,
	const double& dx,
	const double& dz,
	const double& dt
) {
	using std::runtime_error;
	using std::to_string;

	const double
		Ex = +flux_mag[2] * dt * dx / 4,
		Ez = -flux_mag[0] * dt * dz / 4;

	if (flux_neigh.face_neighbor != +2) {
		throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
	}

	const auto
		cleni = grid.mapping.get_cell_length_in_indices(cell.id),
		flux_nleni = grid.mapping.get_cell_length_in_indices(flux_neigh.id);

	auto& cfdb = Face_dB.data(*cell.data);
	{
		auto& nfdb = Face_dB.data(*flux_neigh.data);

		if (flux_neigh.relative_size < 0) {
			if (flux_neigh.x == 0) {
				nfdb(-2) += Ez;
			} else if (cleni == flux_neigh.x + flux_nleni) {
				nfdb(-2) -= Ez;
			} else {
				throw runtime_error(__FILE__ "("
					+ to_string(__LINE__) + ")");
			}
			if (flux_neigh.z == 0) {
				nfdb(-2) -= Ex;
			} else if (cleni == flux_neigh.z + flux_nleni) {
				nfdb(-2) += Ex;
			} else {
				throw runtime_error(__FILE__ "("
					+ to_string(__LINE__) + ")");
			}
		}
		if (flux_neigh.relative_size <= 0 or flux_neigh.x == 0) {
			cfdb(-1) += Ez;
		}
		if (flux_neigh.relative_size >= 0 or flux_neigh.x == 0) {
			nfdb(-1) -= Ez;
		}
		if (flux_neigh.relative_size <= 0 or cleni == flux_neigh.x + flux_nleni) {
			cfdb(+1) += Ez;
		}
		if (flux_neigh.relative_size >= 0 or cleni == flux_neigh.x + flux_nleni) {
			nfdb(+1) -= Ez;
		}
		if (flux_neigh.relative_size <= 0 or flux_neigh.z == 0) {
			cfdb(-3) -= Ex;
		}
		if (flux_neigh.relative_size >= 0 or flux_neigh.z == 0) {
			nfdb(-3) += Ex;
		}
		if (flux_neigh.relative_size <= 0 or cleni == flux_neigh.z + flux_nleni) {
			cfdb(+3) -= Ex;
		}
		if (flux_neigh.relative_size >= 0 or cleni == flux_neigh.z + flux_nleni) {
			nfdb(+3) += Ex;
		}
	}

	for (const auto& neighbor: cell.neighbors_of) {
		if (neighbor.data == nullptr) continue;
		const auto& fn = neighbor.face_neighbor;
		const auto& en = neighbor.edge_neighbor;
		if (fn == 0 and en[0] < 0) continue;

		const auto nleni = grid.mapping.get_cell_length_in_indices(neighbor.id);
		auto& nfdb = Face_dB.data(*neighbor.data);

		if (fn == +2 and neighbor.id != flux_neigh.id) {
			if (neighbor.x != flux_neigh.x and neighbor.z != flux_neigh.z) {
				// no effect
			} else if (neighbor.x < flux_neigh.x) {
				cfdb(+2) -= Ez;
				nfdb(-2) -= Ez;
				nfdb(+1) -= Ez;
			} else if (neighbor.x > flux_neigh.x) {
				cfdb(+2) += Ez;
				nfdb(-2) += Ez;
				nfdb(-1) -= Ez;
			} else if (neighbor.z < flux_neigh.z) {
				cfdb(+2) += Ex;
				nfdb(-2) += Ex;
				nfdb(+3) += Ex;
			} else if (neighbor.z > flux_neigh.z) {
				cfdb(+2) -= Ex;
				nfdb(-2) -= Ex;
				nfdb(-3) += Ex;
			} else {
				throw std::runtime_error(__FILE__ "("
					+ std::to_string(__LINE__) + ")");
			}
		}

		const auto [nEx, nEz] = [&](){
			if (
				neighbor.relative_size > 0
				and flux_neigh.relative_size == 0
			) {
				return std::make_tuple(Ex/2, Ez/2);
			} else {
				return std::make_tuple(Ex, Ez);
			}
		}();

		if (cleni == neighbor.y + nleni) {
			if (fn == -1
				and (flux_neigh.relative_size <= 0
					or (flux_neigh.x == 0
						and (neighbor.relative_size <= 0
							or neighbor.z == flux_neigh.z)))
			) {
				nfdb(+1) += nEz;
				nfdb(+2) -= nEz;
			}
			if (fn == +1
				and (flux_neigh.relative_size <= 0
					or (cleni == flux_neigh.x + flux_nleni
						and (neighbor.relative_size <= 0
							or neighbor.z == flux_neigh.z)))
			) {
				nfdb(-1) += nEz;
				nfdb(+2) += nEz;
			}
			if (fn == -3
				and (flux_neigh.relative_size <= 0
					or (flux_neigh.z == 0
						and (neighbor.relative_size <= 0
							or neighbor.x == flux_neigh.x)))
			) {
				nfdb(+2) += nEx;
				nfdb(+3) -= nEx;
			}
			if (fn == +3
				and (flux_neigh.relative_size <= 0
					or (cleni == flux_neigh.z + flux_nleni
						and (neighbor.relative_size <= 0
							or neighbor.x == flux_neigh.x)))
			) {
				nfdb(+2) -= nEx;
				nfdb(-3) -= nEx;
			}
		}
		if (
			en[0] == 0
			and en[1] == +1
			and en[2] == -1
			and flux_neigh.z == 0
			and (neighbor.relative_size <= 0
				or flux_neigh.relative_size <= 0
				or neighbor.x == flux_neigh.x)
		) {
			nfdb(-2) += nEx;
			nfdb(+3) += nEx;
		}
		if (
			en[0] == 0
			and en[1] == +1
			and en[2] == +1
			and cleni == flux_neigh.z + flux_nleni
			and (neighbor.relative_size <= 0
				or flux_neigh.relative_size <= 0
				or neighbor.x == flux_neigh.x)
		) {
			nfdb(-2) -= nEx;
			nfdb(-3) += nEx;
		}
		if (
			en[0] == 2
			and en[1] == -1
			and en[2] == +1
			and flux_neigh.x == 0
			and (neighbor.relative_size <= 0
				or flux_neigh.relative_size <= 0
				or neighbor.z == flux_neigh.z)
		) {
			nfdb(+1) -= nEz;
			nfdb(-2) -= nEz;
		}
		if (
			en[0] == 2
			and en[1] == +1
			and en[2] == +1
			and cleni == flux_neigh.x + flux_nleni
			and (neighbor.relative_size <= 0
				or flux_neigh.relative_size <= 0
				or neighbor.z == flux_neigh.z)
		) {
			nfdb(-1) -= nEz;
			nfdb(-2) += nEz;
		}
	}
}

template <
	class Grid,
	class Cell_Iter,
	class Flux_Vec,
	class Face_dB_Getter
> void assign_missing_face_dBs_fy(
	const Grid& grid,
	const Cell_Iter& cell,
	const Flux_Vec& flux_mag,
	const Face_dB_Getter& Face_dB,
	const double& dx,
	const double& dz,
	const double& dt
) {
	using std::runtime_error;
	using std::to_string;

	const double
		nEx = +flux_mag[2] * dt * dx / 8,
		nEz = -flux_mag[0] * dt * dz / 8;

	const auto cleni = grid.mapping.get_cell_length_in_indices(cell.id);

	for (const auto& neighbor: cell.neighbors_of) {
		if (neighbor.data == nullptr) continue;
		if (neighbor.relative_size <= 0) continue;

		const auto& fn = neighbor.face_neighbor;
		if (fn == 0 or abs(fn) == 2) continue;

		const auto nleni = grid.mapping.get_cell_length_in_indices(neighbor.id);
		auto& nfdb = Face_dB.data(*neighbor.data);

		if (fn == -1) {
			if (neighbor.y == 0) {
				nfdb(+1) += nEz;
				nfdb(+2) -= nEz;
			} else if (cleni == neighbor.y + nleni) {
				nfdb(+1) -= nEz;
				nfdb(-2) -= nEz;
			} else {
				throw runtime_error(__FILE__ "("
					+ to_string(__LINE__) + ")");
			}
		}
		if (fn == +1) {
			if (neighbor.y == 0) {
				nfdb(-1) += nEz;
				nfdb(+2) += nEz;
			} else if (cleni == neighbor.y + nleni) {
				nfdb(-1) -= nEz;
				nfdb(-2) += nEz;
			} else {
				throw runtime_error(__FILE__ "("
					+ to_string(__LINE__) + ")");
			}
		}
		if (fn == -3) {
			if (neighbor.y == 0) {
				nfdb(+2) += nEx;
				nfdb(+3) -= nEx;
			} else if (cleni == neighbor.y + nleni) {
				nfdb(-2) += nEx;
				nfdb(+3) += nEx;
			} else {
				throw runtime_error(__FILE__ "("
					+ to_string(__LINE__) + ")");
			}
		}
		if (fn == +3) {
			if (neighbor.y == 0) {
				nfdb(+2) -= nEx;
				nfdb(-3) -= nEx;
			} else if (cleni == neighbor.y + nleni) {
				nfdb(-2) -= nEx;
				nfdb(-3) += nEx;
			} else {
				throw runtime_error(__FILE__ "("
					+ to_string(__LINE__) + ")");
			}
		}
	}
}

template <
	class Grid,
	class Cell_Iter,
	class Neighbor_Iter,
	class Flux_Vec,
	class Face_dB_Getter
> void assign_face_dBs_fz(
	const Grid& grid,
	const Cell_Iter& cell,
	const Neighbor_Iter& flux_neigh,
	const Flux_Vec& flux_mag,
	const Face_dB_Getter& Face_dB,
	const double& dx,
	const double& dy,
	const double& dt
) {
	using std::runtime_error;
	using std::to_string;

	const double
		Ex = -flux_mag[1] * dt * dx / 4,
		Ey = +flux_mag[0] * dt * dy / 4;

	if (flux_neigh.face_neighbor != +3) {
		throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
	}

	const auto
		cleni = grid.mapping.get_cell_length_in_indices(cell.id),
		flux_nleni = grid.mapping.get_cell_length_in_indices(flux_neigh.id);
	auto& cfdb = Face_dB.data(*cell.data);
	{
		auto& nfdb = Face_dB.data(*flux_neigh.data);

		if (flux_neigh.relative_size < 0) {
			if (flux_neigh.x == 0) {
				nfdb(-3) -= Ey;
			} else if (cleni == flux_neigh.x + flux_nleni) {
				nfdb(-3) += Ey;
			} else {
				throw runtime_error(__FILE__ "("
					+ to_string(__LINE__) + ")");
			}
			if (flux_neigh.y == 0) {
				nfdb(-3) += Ex;
			} else if (cleni == flux_neigh.y + flux_nleni) {
				nfdb(-3) -= Ex;
			} else {
				throw runtime_error(__FILE__ "("
					+ to_string(__LINE__) + ")");
			}
		}
		if (flux_neigh.relative_size <= 0 or flux_neigh.x == 0) {
			cfdb(-1) -= Ey;
		}
		if (flux_neigh.relative_size >= 0 or flux_neigh.x == 0) {
			nfdb(-1) += Ey;
		}
		if (flux_neigh.relative_size <= 0 or cleni == flux_neigh.x + flux_nleni) {
			cfdb(+1) -= Ey;
		}
		if (flux_neigh.relative_size >= 0 or cleni == flux_neigh.x + flux_nleni) {
			nfdb(+1) += Ey;
		}
		if (flux_neigh.relative_size <= 0 or flux_neigh.y == 0) {
			cfdb(-2) += Ex;
		}
		if (flux_neigh.relative_size >= 0 or flux_neigh.y == 0) {
			nfdb(-2) -= Ex;
		}
		if (flux_neigh.relative_size <= 0 or cleni == flux_neigh.y + flux_nleni) {
			cfdb(+2) += Ex;
		}
		if (flux_neigh.relative_size >= 0 or cleni == flux_neigh.y + flux_nleni) {
			nfdb(+2) -= Ex;
		}
	}

	for (const auto& neighbor: cell.neighbors_of) {
		if (neighbor.data == nullptr) continue;
		const auto& fn = neighbor.face_neighbor;
		const auto& en = neighbor.edge_neighbor;
		if (fn == 0 and en[0] < 0) continue;

		const auto nleni = grid.mapping.get_cell_length_in_indices(neighbor.id);
		auto& nfdb = Face_dB.data(*neighbor.data);

		if (fn == +3 and neighbor.id != flux_neigh.id) {
			if (neighbor.x != flux_neigh.x and neighbor.y != flux_neigh.y) {
				// no effect
			} else if (neighbor.x < flux_neigh.x) {
				cfdb(+3) += Ey;
				nfdb(-3) += Ey;
				nfdb(+1) += Ey;
			} else if (neighbor.x > flux_neigh.x) {
				cfdb(+3) -= Ey;
				nfdb(-3) -= Ey;
				nfdb(-1) += Ey;
			} else if (neighbor.y < flux_neigh.y) {
				cfdb(+3) -= Ex;
				nfdb(-3) -= Ex;
				nfdb(+2) -= Ex;
			} else if (neighbor.y > flux_neigh.y) {
				cfdb(+3) += Ex;
				nfdb(-3) += Ex;
				nfdb(-2) -= Ex;
			} else {
				throw std::runtime_error(__FILE__ "("
					+ std::to_string(__LINE__) + ")");
			}
		}

		const auto [nEx, nEy] = [&](){
			if (
				neighbor.relative_size > 0
				and flux_neigh.relative_size == 0
			) {
				return std::make_tuple(Ex/2, Ey/2);
			} else {
				return std::make_tuple(Ex, Ey);
			}
		}();

		if (cleni == neighbor.z + nleni) {
			if (fn == -1
				and (flux_neigh.relative_size <= 0
					or (flux_neigh.x == 0
						and (neighbor.relative_size <= 0
							or neighbor.y == flux_neigh.y)))
			) {
				nfdb(+1) -= nEy;
				nfdb(+3) += nEy;
			}
			if (fn == +1
				and (flux_neigh.relative_size <= 0
					or (cleni == flux_neigh.x + flux_nleni
						and (neighbor.relative_size <= 0
							or neighbor.y == flux_neigh.y)))
			) {
				nfdb(-1) -= nEy;
				nfdb(+3) -= nEy;
			}
			if (fn == -2
				and (flux_neigh.relative_size <= 0
					or (flux_neigh.y == 0
						and (neighbor.relative_size <= 0
							or neighbor.x == flux_neigh.x)))
			) {
				nfdb(+2) += nEx;
				nfdb(+3) -= nEx;
			}
			if (fn == +2
				and (flux_neigh.relative_size <= 0
					or (cleni == flux_neigh.y + flux_nleni
						and (neighbor.relative_size <= 0
							or neighbor.x == flux_neigh.x)))
			) {
				nfdb(-2) += nEx;
				nfdb(+3) += nEx;
			}
		}
		if (
			en[0] == 0
			and en[1] == -1
			and en[2] == +1
			and flux_neigh.y == 0
			and (neighbor.relative_size <= 0
				or flux_neigh.relative_size <= 0
				or neighbor.x == flux_neigh.x)
		) {
			nfdb(+2) -= nEx;
			nfdb(-3) -= nEx;
		}
		if (
			en[0] == 0
			and en[1] == +1
			and en[2] == +1
			and cleni == flux_neigh.y + flux_nleni
			and (neighbor.relative_size <= 0
				or flux_neigh.relative_size <= 0
				or neighbor.x == flux_neigh.x)
		) {
			nfdb(-2) -= nEx;
			nfdb(-3) += nEx;
		}
		if (
			en[0] == 1
			and en[1] == -1
			and en[2] == +1
			and flux_neigh.x == 0
			and (neighbor.relative_size <= 0
				or flux_neigh.relative_size <= 0
				or neighbor.y == flux_neigh.y)
		) {
			nfdb(+1) += nEy;
			nfdb(-3) += nEy;
		}
		if (
			en[0] == 1
			and en[1] == +1
			and en[2] == +1
			and cleni == flux_neigh.x + flux_nleni
			and (neighbor.relative_size <= 0
				or flux_neigh.relative_size <= 0
				or neighbor.y == flux_neigh.y)
		) {
			nfdb(-1) += nEy;
			nfdb(-3) -= nEy;
		}
	}
}

template <
	class Grid,
	class Cell_Iter,
	class Flux_Vec,
	class Face_dB_Getter
> void assign_missing_face_dBs_fz(
	const Grid& grid,
	const Cell_Iter& cell,
	const Flux_Vec& flux_mag,
	const Face_dB_Getter& Face_dB,
	const double& dx,
	const double& dy,
	const double& dt
) {
	using std::runtime_error;
	using std::to_string;

	const double
		nEx = -flux_mag[1] * dt * dx / 8,
		nEy = +flux_mag[0] * dt * dy / 8;

	const auto cleni = grid.mapping.get_cell_length_in_indices(cell.id);

	for (const auto& neighbor: cell.neighbors_of) {
		if (neighbor.data == nullptr) continue;
		if (neighbor.relative_size <= 0) continue;

		const auto& fn = neighbor.face_neighbor;
		if (fn == 0 or abs(fn) == 3) continue;

		const auto nleni = grid.mapping.get_cell_length_in_indices(neighbor.id);
		auto& nfdb = Face_dB.data(*neighbor.data);

		if (fn == -1) {
			if (neighbor.z == 0) {
				nfdb(+1) -= nEy;
				nfdb(+3) += nEy;
			} else if (cleni == neighbor.z + nleni) {
				nfdb(+1) += nEy;
				nfdb(-3) += nEy;
			} else {
				throw runtime_error(__FILE__ "("
					+ to_string(__LINE__) + ")");
			}
		}
		if (fn == +1) {
			if (neighbor.z == 0) {
				nfdb(-1) -= nEy;
				nfdb(+3) -= nEy;
			} else if (cleni == neighbor.z + nleni) {
				nfdb(-1) += nEy;
				nfdb(-3) -= nEy;
			} else {
				throw runtime_error(__FILE__ "("
					+ to_string(__LINE__) + ")");
			}
		}
		if (fn == -2) {
			if (neighbor.z == 0) {
				nfdb(+2) += nEx;
				nfdb(+3) -= nEx;
			} else if (cleni == neighbor.z + nleni) {
				nfdb(+2) -= nEx;
				nfdb(-3) -= nEx;
			} else {
				throw runtime_error(__FILE__ "("
					+ to_string(__LINE__) + ")");
			}
		}
		if (fn == +2) {
			if (neighbor.z == 0) {
				nfdb(-2) += nEx;
				nfdb(+3) += nEx;
			} else if (cleni == neighbor.z + nleni) {
				nfdb(-2) -= nEx;
				nfdb(-3) += nEx;
			} else {
				throw runtime_error(__FILE__ "("
					+ to_string(__LINE__) + ")");
			}
		}
	}
}


/*!
Set new MHD state in cells with SInfo.data(*cell.data) == 1.

Sets new MHD state based on fluxes and face B changes.
Zeroes fluxes and face B changes afterwards.
*/
template <
	class Cells,
	class Grid,
	class Face_Magnetic_Field_Getter,
	class Face_dB_Getter,
	class Solver_Info_Getter,
	class Substepping_Period_Getter,
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Volume_Magnetic_Field_Getter,
	class Mass_Density_Flux_Getters,
	class Momentum_Density_Flux_Getters,
	class Total_Energy_Density_Flux_Getters,
	class Magnetic_Field_Flux_Getters
> void update_mhd_state(
	const Cells& cells,
	Grid& grid,
	const int& current_substep,
	const Face_Magnetic_Field_Getter& Face_B,
	const Face_dB_Getter& Face_dB,
	const Solver_Info_Getter& SInfo,
	const Substepping_Period_Getter& Substep,
	const Mass_Density_Getter& Mas,
	const Momentum_Density_Getter& Mom,
	const Total_Energy_Density_Getter& Nrj,
	const Volume_Magnetic_Field_Getter& Vol_B,
	const Mass_Density_Flux_Getters& Mas_f,
	const Momentum_Density_Flux_Getters& Mom_f,
	const Total_Energy_Density_Flux_Getters& Nrj_f,
	const Magnetic_Field_Flux_Getters& Mag_f
) try {
	using std::get;
	using std::runtime_error;
	using std::to_string;

	for (const auto& cell: cells) {
		if (SInfo.data(*cell.data) > 0) {
			if (current_substep % Substep.data(*cell.data) != 0) {
				continue;
			}

			const auto [dx, dy, dz] = grid.geometry.get_length(cell.id);

			Mas.data(*cell.data) += (Mas_f(*cell.data, -1) - Mas_f(*cell.data, +1)) / dx;
			Mas.data(*cell.data) += (Mas_f(*cell.data, -2) - Mas_f(*cell.data, +2)) / dy;
			Mas.data(*cell.data) += (Mas_f(*cell.data, -3) - Mas_f(*cell.data, +3)) / dz;
			Mom.data(*cell.data) += (Mom_f(*cell.data, -1) - Mom_f(*cell.data, +1)) / dx;
			Nrj.data(*cell.data) += (Nrj_f(*cell.data, -1) - Nrj_f(*cell.data, +1)) / dx;
			Vol_B.data(*cell.data) += (Mag_f(*cell.data, -1) - Mag_f(*cell.data, +1)) / dx;
			Mom.data(*cell.data) += (Mom_f(*cell.data, -2) - Mom_f(*cell.data, +2)) / dy;
			Nrj.data(*cell.data) += (Nrj_f(*cell.data, -2) - Nrj_f(*cell.data, +2)) / dy;
			Vol_B.data(*cell.data) += (Mag_f(*cell.data, -2) - Mag_f(*cell.data, +2)) / dy;
			Mom.data(*cell.data) += (Mom_f(*cell.data, -3) - Mom_f(*cell.data, +3)) / dz;
			Nrj.data(*cell.data) += (Nrj_f(*cell.data, -3) - Nrj_f(*cell.data, +3)) / dz;
			Vol_B.data(*cell.data) += (Mag_f(*cell.data, -3) - Mag_f(*cell.data, +3)) / dz;

			const std::array<double, 3> area{dy*dz, dx*dz, dx*dy};
			for (size_t dim: {0, 1, 2}) {
				Face_B.data(*cell.data)(dim, -1) += Face_dB.data(*cell.data)(dim, -1) / area[dim];
				Face_B.data(*cell.data)(dim, +1) += Face_dB.data(*cell.data)(dim, +1) / area[dim];
			}
		}

		for (int dir: {-3,-2,-1,+1,+2,+3}) {
			Mas_f(*cell.data, dir) =
			Nrj_f(*cell.data, dir) = 0;
			Mom_f(*cell.data, dir) =
			Mag_f(*cell.data, dir) = {0, 0, 0};
		}
		Face_dB.data(*cell.data) = {0, 0, 0, 0, 0, 0};
	}
	Mas.type().is_stale = true;
	Mom.type().is_stale = true;
	Nrj.type().is_stale = true;
	Vol_B.type().is_stale = true;
	Face_B.type().is_stale = true;

} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


/*!
Makes B consistent in given cells.

Face B can have two values at any location since
cells store it at every face. FMagP is used by default
so FMagN in neighbor on positive side should be equal.
With adaptive mesh refinement, in inreasing order of
importance:

1) Vol_B = FMagN = FmagP if neighbors missing or type < 0
2) Face B of smaller neighbors overrides any face B
3) FMagN overrides FMagP of larger neighbor on negative side
4) Vol_B = 0.5*(FmagN+FMagP)

After above decisions remaining face B values which were
overriden are copied/averaged from other cell(s).

If constant_thermal_pressure == true total energy is
adjusted after averaging volume B.
*/
template <
	class Cell_Iter,
	class Grid,
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Volume_Magnetic_Field_Getter,
	class Face_Magnetic_Field_Getter,
	class Solver_Info_Getter,
	class Substepping_Period_Getter
> void update_B_consistency(
	const int current_substep,
	const Cell_Iter& cells,
	Grid& grid,
	const Mass_Density_Getter& Mas,
	const Momentum_Density_Getter& Mom,
	const Total_Energy_Density_Getter& Nrj,
	const Volume_Magnetic_Field_Getter& Vol_B,
	const Face_Magnetic_Field_Getter& Face_B,
	const Solver_Info_Getter& SInfo,
	const Substepping_Period_Getter& Substep,
	const double& adiabatic_index,
	const double& vacuum_permeability,
	const bool& constant_thermal_pressure
) try {
	using std::runtime_error;
	using std::to_string;

	using Cell = Grid::cell_data_type;

	bool update_copies = false;
	if (Face_B.type().is_stale) {
		update_copies = true;
		Cell::set_transfer_all(true, Face_B.type());
		Face_B.type().is_stale = false;
	}
	if (SInfo.type().is_stale) {
		update_copies = true;
		Cell::set_transfer_all(true, SInfo.type());
		SInfo.type().is_stale = false;
	}
	if (update_copies) {
		grid.update_copies_of_remote_neighbors();
	}
	Cell::set_transfer_all(false, Face_B.type(), SInfo.type());

	for (const auto& cell: cells) {
		if (SInfo.data(*cell.data) < 0) continue;
		if (current_substep % Substep.data(*cell.data) != 0) continue;
		if (constant_thermal_pressure and Mas.data(*cell.data) <= 0) continue;

		const auto old_pressure = [&](){
			if (constant_thermal_pressure) {
				return pamhd::mhd::get_pressure(
					Mas.data(*cell.data), Mom.data(*cell.data), Nrj.data(*cell.data), Vol_B.data(*cell.data),
					adiabatic_index, vacuum_permeability);
			} else {
				return 0.0;
			}
		}();

		const auto& c_face_b = Face_B.data(*cell.data);
		for (const auto& neighbor: cell.neighbors_of) {
			if (SInfo.data(*neighbor.data) < 0) continue;
			const auto& fn = neighbor.face_neighbor;
			if (fn <= 0) continue;

			const auto dim = abs(fn) - 1;
			const auto& n_face_b = Face_B.data(*neighbor.data);
			if (c_face_b(dim, +1) != n_face_b(dim, -1)) {
//FIXME				abort();
			}
		}
		for (size_t dim: {0, 1, 2}) {
			Vol_B.data(*cell.data)[dim] = 0.5 * (c_face_b(dim, -1) + c_face_b(dim, +1));
		}
		if (constant_thermal_pressure) {
			const auto vel = (Mom.data(*cell.data)/Mas.data(*cell.data)).eval();
			Nrj.data(*cell.data) = pamhd::mhd::get_total_energy_density(
				Mas.data(*cell.data), vel, old_pressure, Vol_B.data(*cell.data),
				adiabatic_index, vacuum_permeability
			);
		}
	}
	if (constant_thermal_pressure) {
		Nrj.type().is_stale = true;
	}
	Vol_B.type().is_stale = true;

} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}

template <
	class Solver,
	class Grid,
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Volume_Magnetic_Field_Getter,
	class Face_Magnetic_Field_Getter,
	class Face_dB_Getter,
	class Background_Magnetic_Field_Getter,
	class Mass_Density_Flux_Getters,
	class Momentum_Density_Flux_Getters,
	class Total_Energy_Density_Flux_Getters,
	class Magnetic_Field_Flux_Getters,
	class Solver_Info_Getter,
	class Substepping_Period_Getter,
	class Max_Velocity_Getter
> void get_all_fluxes(
	const double& sub_dt,
	const Solver& solver,
	Grid& grid,
	const int& substep,
	const double& adiabatic_index,
	const double& vacuum_permeability,
	const Mass_Density_Getter& Mas,
	const Momentum_Density_Getter& Mom,
	const Total_Energy_Density_Getter& Nrj,
	const Volume_Magnetic_Field_Getter& Vol_B,
	const Face_Magnetic_Field_Getter& Face_B,
	const Face_dB_Getter& Face_dB,
	const Background_Magnetic_Field_Getter& Bg_B,
	const Mass_Density_Flux_Getters& Mas_f,
	const Momentum_Density_Flux_Getters& Mom_f,
	const Total_Energy_Density_Flux_Getters& Nrj_f,
	const Magnetic_Field_Flux_Getters& Mag_f,
	const Solver_Info_Getter& SInfo,
	const Substepping_Period_Getter& Substep,
	const Max_Velocity_Getter& Max_v
) try {
	using Cell = Grid::cell_data_type;

	bool update_copies = false;
	if (
		Mas.type().is_stale or Mom.type().is_stale
		or Nrj.type().is_stale or Vol_B.type().is_stale
	) {
		update_copies = true;
		Cell::set_transfer_all(true,
			Mas.type(), Mom.type(), Nrj.type(), Vol_B.type());
		Mas.type().is_stale = false;
		Mom.type().is_stale = false;
		Nrj.type().is_stale = false;
		Vol_B.type().is_stale = false;
	}
	if (SInfo.type().is_stale) {
		update_copies = true;
		Cell::set_transfer_all(true, SInfo.type());
		SInfo.type().is_stale = false;
	}
	if (Substep.type().is_stale) {
		update_copies = true;
		Cell::set_transfer_all(true, Substep.type());
		Substep.type().is_stale = false;
	}
	if (Bg_B.type().is_stale) {
		update_copies = true;
		Cell::set_transfer_all(true, Bg_B.type());
		Bg_B.type().is_stale = false;
	}
	if (update_copies) {
		grid.start_remote_neighbor_copy_updates();
	}

	pamhd::mhd::get_fluxes(
		solver, grid.inner_cells(), grid, substep,
		adiabatic_index, vacuum_permeability, sub_dt,
		Mas, Mom, Nrj, Vol_B, Face_B, Face_dB, Bg_B,
		Mas_f, Mom_f, Nrj_f, Mag_f,
		SInfo, Substep, Max_v
	);

	if (update_copies) {
		grid.wait_remote_neighbor_copy_update_receives();
	}

	pamhd::mhd::get_fluxes(
		solver, grid.outer_cells(), grid, substep,
		adiabatic_index, vacuum_permeability, sub_dt,
		Mas, Mom, Nrj, Vol_B, Face_B, Face_dB, Bg_B,
		Mas_f, Mom_f, Nrj_f, Mag_f,
		SInfo, Substep, Max_v
	);

	if (update_copies) {
		grid.wait_remote_neighbor_copy_update_sends();
	}
	Cell::set_transfer_all(false,
		Mas.type(), Mom.type(), Nrj.type(), Vol_B.type(),
		SInfo.type(), Substep.type(), Bg_B.type());

	pamhd::mhd::get_fluxes(
		solver, grid.remote_cells(), grid, substep,
		adiabatic_index, vacuum_permeability, sub_dt,
		Mas, Mom, Nrj, Vol_B, Face_B, Face_dB, Bg_B,
		Mas_f, Mom_f, Nrj_f, Mag_f,
		SInfo, Substep, Max_v
	);
	Max_v.type().is_stale = true;

} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


//! Returns length of timestep taken.
template <
	class Solver,
	class Grid,
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Magnetic_Field_Getter,
	class Face_Magnetic_Field_Getter,
	class Face_dB_Getter,
	class Background_Magnetic_Field_Getter,
	class Mass_Density_Flux_Getters,
	class Momentum_Density_Flux_Getters,
	class Total_Energy_Density_Flux_Getters,
	class Magnetic_Field_Flux_Getters,
	class Solver_Info_Getter,
	class Timestep_Getter,
	class Substepping_Period_Getter,
	class Substep_Min_Getter,
	class Substep_Max_Getter,
	class Max_Velocity_Getter
> double timestep(
	const Solver solver,
	Grid& grid,
	Options& options_mhd,
	const double simulation_time,
	double max_time_step,
	const double time_step_factor,
	const double adiabatic_index,
	const double vacuum_permeability,
	const Mass_Density_Getter& Mas,
	const Momentum_Density_Getter& Mom,
	const Total_Energy_Density_Getter& Nrj,
	const Magnetic_Field_Getter& Vol_B,
	const Face_Magnetic_Field_Getter& Face_B,
	const Face_dB_Getter& Face_dB,
	const Background_Magnetic_Field_Getter& Bg_B,
	const Mass_Density_Flux_Getters& Mas_f,
	const Momentum_Density_Flux_Getters& Mom_f,
	const Total_Energy_Density_Flux_Getters& Nrj_f,
	const Magnetic_Field_Flux_Getters& Mag_f,
	const Solver_Info_Getter& SInfo,
	const Timestep_Getter& Timestep,
	const Substepping_Period_Getter& Substep,
	const Substep_Min_Getter& Substep_Min,
	const Substep_Max_Getter& Substep_Max,
	const Max_Velocity_Getter& Max_v
) try {
	set_minmax_substepping_period(
		simulation_time, grid, options_mhd,
		Substep_Min, Substep_Max
	);

	minimize_timestep(
		grid, time_step_factor, SInfo, Timestep, Max_v
	);

	const double sub_dt = pamhd::set_minmax_substepping_period(
		grid, max_time_step, SInfo,
		Timestep, Substep_Min, Substep_Max
	);

	restrict_substepping_period(
		grid, Substep, Substep_Max, SInfo
	);

	const int max_substep = update_substeps(grid, SInfo, Substep);
	if (grid.get_rank() == 0) {std::cout
		<< "Substep: " << sub_dt << ", largest substep period: "
		<< max_substep << std::endl;
	}

	double total_dt = 0;
	for (int substep = 1; substep <= max_substep; substep += 1) {
		total_dt += sub_dt;

		get_all_fluxes(
			sub_dt, solver, grid, substep, adiabatic_index,
			vacuum_permeability, Mas, Mom, Nrj, Vol_B,
			Face_B, Face_dB, Bg_B, Mas_f, Mom_f, Nrj_f,
			Mag_f, SInfo, Substep, Max_v
		);

		update_mhd_state(
			grid.local_cells(), grid, substep,
			Face_B, Face_dB, SInfo, Substep,
			Mas, Mom, Nrj, Vol_B,
			Mas_f, Mom_f, Nrj_f, Mag_f
		);

		// constant thermal pressure when updating vol B after solution
		update_B_consistency(
			substep, grid.local_cells(), grid,
			Mas, Mom, Nrj, Vol_B, Face_B,
			SInfo, Substep,
			adiabatic_index,
			vacuum_permeability,
			true
		);
	}

	return total_dt;

} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


//! Restricts Timestep based on MHD physics
template <
	class Grid,
	class Solver_Info_Getter,
	class Timestep_Getter,
	class Max_Velocity_Getter
> void minimize_timestep(
	Grid& grid,
	const double& dt_factor,
	const Solver_Info_Getter& SInfo,
	const Timestep_Getter& Timestep,
	const Max_Velocity_Getter& Max_v
) try {
	using std::min;
	using std::numeric_limits;

	using Cell = Grid::cell_data_type;

	bool update_copies = false;
	if (Max_v.type().is_stale) {
		update_copies = true;
		Cell::set_transfer_all(true, Max_v.type());
		Max_v.type().is_stale = false;
	}
	if (SInfo.type().is_stale) {
		update_copies = true;
		Cell::set_transfer_all(true, SInfo.type());
		SInfo.type().is_stale = false;
	}
	if (update_copies) {
		grid.update_copies_of_remote_neighbors();
	}
	Cell::set_transfer_all(false, Max_v.type(), SInfo.type());
	Timestep.type().is_stale = true;

	for (const auto& cell: grid.local_cells()) {
		Timestep.data(*cell.data) = numeric_limits<double>::max();
	}

	for (const auto& cell: grid.all_cells()) {
		if (SInfo.data(*cell.data) < 0) {
			continue;
		}

		const auto cell_len = grid.geometry.get_length(cell.id);
		for (const auto& neighbor: cell.neighbors_of) {
			const auto& fn = neighbor.face_neighbor;
			if (
				fn <= 0
				or neighbor.data == nullptr
				or SInfo.data(*neighbor.data) < 0
			) {
				continue;
			}

			const auto neigh_dim = size_t(abs(fn) - 1);
			const auto neigh_len = grid.geometry.get_length(neighbor.id);
			const auto max_v = std::max(
				Max_v.data(*cell.data)(fn),
				Max_v.data(*neighbor.data)(-fn));

			if (cell.is_local) {
				if (max_v <= 0) {
					Timestep.data(*cell.data) = 0;
				} else {
					Timestep.data(*cell.data) = min(
						dt_factor * cell_len[neigh_dim] / max_v,
						Timestep.data(*cell.data));
				}
			}
			if (neighbor.is_local) {
				if (max_v <= 0) {
					Timestep.data(*neighbor.data) = 0;
				} else {
					Timestep.data(*neighbor.data) = min(
						dt_factor * neigh_len[neigh_dim] / max_v,
						Timestep.data(*neighbor.data));
				}
			}
		}
	}

} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


}} // namespaces


#endif // ifndef PAMHD_MHD_SOLVE_STAGGERED_HPP
