/*
Solves the MHD part of PAMHD using an external flux function.

Copyright 2014, 2015, 2016, 2017 Ilja Honkonen
Copyright 2018, 2019, 2022, 2023, 2024 Finnish Meteorological Institute
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
	const Magnetic_Field_Getter& Mag,
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
	state_neg[mas_int] = Mas(*cell.data);
	state_neg[mom_int] = get_rotated_vector(Mom(*cell.data), abs(dir));
	state_neg[nrj_int] = Nrj(*cell.data);
	state_neg[mag_int] = get_rotated_vector(Mag(*cell.data), abs(dir));

	state_pos[mas_int] = Mas(*neighbor.data);
	state_pos[mom_int] = get_rotated_vector(Mom(*neighbor.data), abs(dir));
	state_pos[nrj_int] = Nrj(*neighbor.data);
	state_pos[mag_int] = get_rotated_vector(Mag(*neighbor.data), abs(dir));

	const Magnetic_Field::data_type bg_face_b
		= get_rotated_vector(Bg_B(*cell.data)(dir), abs(dir));
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
			<< " of cell type " << SInfo(*cell.data)
			<< " and " << SInfo(*neighbor.data)
			<< " with substeps " << Substep(*cell.data)
			<< " and " << Substep(*neighbor.data)
			<< " at " << grid.geometry.get_center(cell.id)
			<< " and " << grid.geometry.get_center(neighbor.id)
			<< " in direction " << dir
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

	// rotate flux back
	flux[mom_int] = get_rotated_vector(flux[mom_int], -abs(dir));
	flux[mag_int] = get_rotated_vector(flux[mag_int], -abs(dir));

	return std::make_tuple(max_vel, flux);
}


/*!
Calculates MHD fluxes in/out of given cells.

Flux getters with array indices [0..5] should return variables in this order:
-x,+x,-y,+y,-z,+z.

Saves fluxes of cells with SInfo(*cell_data) == 1,
ignores cells with SInfo < 0.
*/
template <
	class Cell_Iter,
	class Grid,
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Magnetic_Field_Getter,
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
	const Solver solver,
	const Cell_Iter& cells,
	Grid& grid,
	const int current_substep,
	const double adiabatic_index,
	const double vacuum_permeability,
	const double dt,
	const Mass_Density_Getter Mas,
	const Momentum_Density_Getter Mom,
	const Total_Energy_Density_Getter Nrj,
	const Magnetic_Field_Getter Mag,
	const Face_dB_Getter Face_dB,
	const Background_Magnetic_Field_Getter Bg_B,
	const Mass_Density_Flux_Getters Mas_f,
	const Momentum_Density_Flux_Getters Mom_f,
	const Total_Energy_Density_Flux_Getters Nrj_f,
	const Magnetic_Field_Flux_Getters Mag_f,
	const Solver_Info_Getter SInfo,
	const Substepping_Period_Getter Substep,
	const Max_Velocity_Getter Max_v
) try {
	using std::abs;
	using std::get;
	using std::min;
	using std::runtime_error;
	using std::tie;
	using std::to_string;

	const auto
		Mas_fnx = get<0>(Mas_f), Mas_fpx = get<1>(Mas_f),
		Mas_fny = get<2>(Mas_f), Mas_fpy = get<3>(Mas_f),
		Mas_fnz = get<4>(Mas_f), Mas_fpz = get<5>(Mas_f),
		Mom_fnx = get<0>(Mom_f), Mom_fpx = get<1>(Mom_f),
		Mom_fny = get<2>(Mom_f), Mom_fpy = get<3>(Mom_f),
		Mom_fnz = get<4>(Mom_f), Mom_fpz = get<5>(Mom_f),
		Nrj_fnx = get<0>(Nrj_f), Nrj_fpx = get<1>(Nrj_f),
		Nrj_fny = get<2>(Nrj_f), Nrj_fpy = get<3>(Nrj_f),
		Nrj_fnz = get<4>(Nrj_f), Nrj_fpz = get<5>(Nrj_f),
		Mag_fnx = get<0>(Mag_f), Mag_fpx = get<1>(Mag_f),
		Mag_fny = get<2>(Mag_f), Mag_fpy = get<3>(Mag_f),
		Mag_fnz = get<4>(Mag_f), Mag_fpz = get<5>(Mag_f);

	const pamhd::mhd::Mass_Density mas_int{};
	const pamhd::mhd::Momentum_Density mom_int{};
	const pamhd::mhd::Total_Energy_Density nrj_int{};
	const pamhd::Magnetic_Field mag_int{};

	for (const auto& cell: cells) {
		if (cell.data == nullptr) continue;
		if (SInfo(*cell.data) < 0) continue;

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

		const int csub = Substep(*cell.data);
		int min_min_sub = csub;

		const auto [cell_dx, cell_dy, cell_dz]
			= grid.geometry.get_length(cell.id);

		std::array<bool, 3> missing_flux{false, false, false};
		for (const auto& neighbor: cell.neighbors_of) {
			const auto& fn = neighbor.face_neighbor;
			if (fn == 0) continue;
			if (neighbor.data == nullptr) continue;
			if (SInfo(*neighbor.data) < 0) continue;

			const int nsub = Substep(*neighbor.data);
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
				Mag, Bg_B, SInfo, Substep, solver,
				adiabatic_index, vacuum_permeability, dt);

			Max_v(*cell.data)(fn) = max_vel;
			Max_v(*neighbor.data)(-fn) = max_vel;

			// cell size and substep factors for fluxes
			const auto min_sub = min(csub, nsub);
			double cfac = dt*min_sub, nfac = dt*min_sub;
			// average smaller fluxes through large face
			if (neighbor.relative_size > 0) {
				cfac /= 4;
			} else if (neighbor.relative_size < 0) {
				nfac /= 4;
			}

			const auto [neigh_dx, neigh_dy, neigh_dz]
				= grid.geometry.get_length(neighbor.id);

			if (fn == +1) {
				if (cell.is_local) {
					Mas_fpx(*cell.data) += cfac * flux[mas_int];
					Mom_fpx(*cell.data) += cfac * flux[mom_int];
					Nrj_fpx(*cell.data) += cfac * flux[nrj_int];
					Mag_fpx(*cell.data) += cfac * flux[mag_int];
				}

				if (neighbor.is_local) {
					Mas_fnx(*neighbor.data) += nfac * flux[mas_int];
					Mom_fnx(*neighbor.data) += nfac * flux[mom_int];
					Nrj_fnx(*neighbor.data) += nfac * flux[nrj_int];
					Mag_fnx(*neighbor.data) += nfac * flux[mag_int];
				}

				assign_face_dBs_fx(
					grid, cell, neighbor, flux[mag_int], Face_dB,
					min(cell_dy, neigh_dy), min(cell_dz, neigh_dz),
					dt*min_sub);
			}

			if (fn == +2) {
				if (cell.is_local) {
					Mas_fpy(*cell.data) += cfac * flux[mas_int];
					Mom_fpy(*cell.data) += cfac * flux[mom_int];
					Nrj_fpy(*cell.data) += cfac * flux[nrj_int];
					Mag_fpy(*cell.data) += cfac * flux[mag_int];
				}

				if (neighbor.is_local) {
					Mas_fny(*neighbor.data) += nfac * flux[mas_int];
					Mom_fny(*neighbor.data) += nfac * flux[mom_int];
					Nrj_fny(*neighbor.data) += nfac * flux[nrj_int];
					Mag_fny(*neighbor.data) += nfac * flux[mag_int];
				}

				assign_face_dBs_fy(
					grid, cell, neighbor, flux[mag_int], Face_dB,
					min(cell_dx, neigh_dx), min(cell_dz, neigh_dz),
					dt*min_sub);
			}

			if (fn == +3) {
				if (cell.is_local) {
					Mas_fpz(*cell.data) += cfac * flux[mas_int];
					Mom_fpz(*cell.data) += cfac * flux[mom_int];
					Nrj_fpz(*cell.data) += cfac * flux[nrj_int];
					Mag_fpz(*cell.data) += cfac * flux[mag_int];
				}

				if (neighbor.is_local) {
					Mas_fnz(*neighbor.data) += nfac * flux[mas_int];
					Mom_fnz(*neighbor.data) += nfac * flux[mom_int];
					Nrj_fnz(*neighbor.data) += nfac * flux[nrj_int];
					Mag_fnz(*neighbor.data) += nfac * flux[mag_int];
				}

				assign_face_dBs_fz(
					grid, cell, neighbor, flux[mag_int], Face_dB,
					min(cell_dx, neigh_dx), min(cell_dy, neigh_dy),
					dt*min_sub);
			}
		}

		if (missing_flux[0]) {
			const auto [max_vel, flux] = get_flux(
				grid, cell, cell, +1, Mas, Mom, Nrj, Mag, Bg_B,
				SInfo, Substep, solver, adiabatic_index,
				vacuum_permeability, dt);

			assign_missing_face_dBs_fx(
				grid, cell, flux[mag_int], Face_dB,
				cell_dy, cell_dz, dt*min_min_sub);
		}
		if (missing_flux[1]) {
			const auto [max_vel, flux] = get_flux(
				grid, cell, cell, +2, Mas, Mom, Nrj, Mag, Bg_B,
				SInfo, Substep, solver, adiabatic_index,
				vacuum_permeability, dt);

			assign_missing_face_dBs_fy(
				grid, cell, flux[mag_int], Face_dB,
				cell_dx, cell_dz, dt*min_min_sub);
		}
		if (missing_flux[2]) {
			const auto [max_vel, flux] = get_flux(
				grid, cell, cell, +3, Mas, Mom, Nrj, Mag, Bg_B,
				SInfo, Substep, solver, adiabatic_index,
				vacuum_permeability, dt);

			assign_missing_face_dBs_fz(
				grid, cell, flux[mag_int], Face_dB,
				cell_dx, cell_dy, dt*min_min_sub);
		}
	}
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

	auto& cfdb = Face_dB(*cell.data);
	{ // handle cell and flux neighbor separately
		auto
			&nfdb = Face_dB(*flux_neigh.data);
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
		auto& nfdb = Face_dB(*neighbor.data);

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
		auto& nfdb = Face_dB(*neighbor.data);

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

	auto& cfdb = Face_dB(*cell.data);
	{
		auto& nfdb = Face_dB(*flux_neigh.data);

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
		auto& nfdb = Face_dB(*neighbor.data);

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
		auto& nfdb = Face_dB(*neighbor.data);

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
	auto& cfdb = Face_dB(*cell.data);
	{
		auto& nfdb = Face_dB(*flux_neigh.data);

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
		auto& nfdb = Face_dB(*neighbor.data);

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
		auto& nfdb = Face_dB(*neighbor.data);

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
Set new MHD state in cells with SInfo(*cell.data) == 1.

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
	class Magnetic_Field_Getter,
	class Mass_Density_Flux_Getters,
	class Momentum_Density_Flux_Getters,
	class Total_Energy_Density_Flux_Getters,
	class Magnetic_Field_Flux_Getters
> void update_mhd_state(
	const Cells& cells,
	Grid& grid,
	const int current_substep,
	const Face_Magnetic_Field_Getter Face_B,
	const Face_dB_Getter Face_dB,
	const Solver_Info_Getter SInfo,
	const Substepping_Period_Getter Substep,
	const Mass_Density_Getter Mas,
	const Momentum_Density_Getter Mom,
	const Total_Energy_Density_Getter Nrj,
	const Magnetic_Field_Getter Mag,
	const Mass_Density_Flux_Getters Mas_f,
	const Momentum_Density_Flux_Getters Mom_f,
	const Total_Energy_Density_Flux_Getters Nrj_f,
	const Magnetic_Field_Flux_Getters Mag_f
) try {
	using std::get;
	using std::runtime_error;
	using std::to_string;

	const auto
		Mas_fnx = get<0>(Mas_f), Mas_fpx = get<1>(Mas_f),
		Mas_fny = get<2>(Mas_f), Mas_fpy = get<3>(Mas_f),
		Mas_fnz = get<4>(Mas_f), Mas_fpz = get<5>(Mas_f),
		Mom_fnx = get<0>(Mom_f), Mom_fpx = get<1>(Mom_f),
		Mom_fny = get<2>(Mom_f), Mom_fpy = get<3>(Mom_f),
		Mom_fnz = get<4>(Mom_f), Mom_fpz = get<5>(Mom_f),
		Nrj_fnx = get<0>(Nrj_f), Nrj_fpx = get<1>(Nrj_f),
		Nrj_fny = get<2>(Nrj_f), Nrj_fpy = get<3>(Nrj_f),
		Nrj_fnz = get<4>(Nrj_f), Nrj_fpz = get<5>(Nrj_f),
		Mag_fnx = get<0>(Mag_f), Mag_fpx = get<1>(Mag_f),
		Mag_fny = get<2>(Mag_f), Mag_fpy = get<3>(Mag_f),
		Mag_fnz = get<4>(Mag_f), Mag_fpz = get<5>(Mag_f);

	for (const auto& cell: cells) {
		if (SInfo(*cell.data) > 0) {
			if (current_substep % Substep(*cell.data) != 0) {
				continue;
			}

			const auto [dx, dy, dz] = grid.geometry.get_length(cell.id);

			Mas(*cell.data) += (Mas_fnx(*cell.data) - Mas_fpx(*cell.data)) / dx;
			Mas(*cell.data) += (Mas_fny(*cell.data) - Mas_fpy(*cell.data)) / dy;
			Mas(*cell.data) += (Mas_fnz(*cell.data) - Mas_fpz(*cell.data)) / dz;
			Mom(*cell.data) += (Mom_fnx(*cell.data) - Mom_fpx(*cell.data)) / dx;
			Nrj(*cell.data) += (Nrj_fnx(*cell.data) - Nrj_fpx(*cell.data)) / dx;
			Mag(*cell.data) += (Mag_fnx(*cell.data) - Mag_fpx(*cell.data)) / dx;
			Mom(*cell.data) += (Mom_fny(*cell.data) - Mom_fpy(*cell.data)) / dy;
			Nrj(*cell.data) += (Nrj_fny(*cell.data) - Nrj_fpy(*cell.data)) / dy;
			Mag(*cell.data) += (Mag_fny(*cell.data) - Mag_fpy(*cell.data)) / dy;
			Mom(*cell.data) += (Mom_fnz(*cell.data) - Mom_fpz(*cell.data)) / dz;
			Nrj(*cell.data) += (Nrj_fnz(*cell.data) - Nrj_fpz(*cell.data)) / dz;
			Mag(*cell.data) += (Mag_fnz(*cell.data) - Mag_fpz(*cell.data)) / dz;

			const std::array<double, 3> area{dy*dz, dx*dz, dx*dy};
			for (size_t dim: {0, 1, 2}) {
				Face_B(*cell.data)(dim, -1) += Face_dB(*cell.data)(dim, -1) / area[dim];
				Face_B(*cell.data)(dim, +1) += Face_dB(*cell.data)(dim, +1) / area[dim];
			}
		}

		Mas_fnx(*cell.data) =
		Mas_fny(*cell.data) =
		Mas_fnz(*cell.data) =
		Mas_fpx(*cell.data) =
		Mas_fpy(*cell.data) =
		Mas_fpz(*cell.data) =
		Nrj_fnx(*cell.data) =
		Nrj_fny(*cell.data) =
		Nrj_fnz(*cell.data) =
		Nrj_fpx(*cell.data) =
		Nrj_fpy(*cell.data) =
		Nrj_fpz(*cell.data) = 0;

		Mom_fnx(*cell.data) =
		Mom_fny(*cell.data) =
		Mom_fnz(*cell.data) =
		Mom_fpx(*cell.data) =
		Mom_fpy(*cell.data) =
		Mom_fpz(*cell.data) =
		Mag_fnx(*cell.data) =
		Mag_fny(*cell.data) =
		Mag_fnz(*cell.data) =
		Mag_fpx(*cell.data) =
		Mag_fpy(*cell.data) =
		Mag_fpz(*cell.data) = {0, 0, 0};

		Face_dB(*cell.data) = {0, 0, 0, 0, 0, 0};
	}
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

1) VMag = FMagN = FmagP if neighbors missing or type < 0
2) Face B of smaller neighbors overrides any face B
3) FMagN overrides FMagP of larger neighbor on negative side
4) VMag = 0.5*(FmagN+FMagP)

After above decisions remaining face B values which were
overriden are copied/averaged from other cell(s).

If constant_thermal_pressure == true total energy is
adjusted after averaging volume B.
*/
template <
	class Cell_Iter,
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
	const Mass_Density_Getter Mas,
	const Momentum_Density_Getter Mom,
	const Total_Energy_Density_Getter Nrj,
	const Volume_Magnetic_Field_Getter VMag,
	const Face_Magnetic_Field_Getter Face_B,
	const Solver_Info_Getter SInfo,
	const Substepping_Period_Getter Substep,
	const double adiabatic_index,
	const double vacuum_permeability,
	const bool constant_thermal_pressure
) try {
	using std::runtime_error;
	using std::to_string;

	for (const auto& cell: cells) {
		if (SInfo(*cell.data) < 0) continue;
		if (current_substep % Substep(*cell.data) != 0) continue;
		if (constant_thermal_pressure and Mas(*cell.data) <= 0) continue;

		const auto old_pressure = [&](){
			if (constant_thermal_pressure) {
				return pamhd::mhd::get_pressure(
					Mas(*cell.data), Mom(*cell.data), Nrj(*cell.data), VMag(*cell.data),
					adiabatic_index, vacuum_permeability);
			} else {
				return 0.0;
			}
		}();

		const auto& c_face_b = Face_B(*cell.data);
		for (const auto& neighbor: cell.neighbors_of) {
			if (SInfo(*neighbor.data) < 0) continue;
			const auto& fn = neighbor.face_neighbor;
			if (fn <= 0) continue;

			const auto dim = abs(fn) - 1;
			const auto& n_face_b = Face_B(*neighbor.data);
			if (c_face_b(dim, +1) != n_face_b(dim, -1)) {
//FIXME				abort();
			}
		}
		for (size_t dim: {0, 1, 2}) {
			VMag(*cell.data)[dim] = 0.5 * (c_face_b(dim, -1) + c_face_b(dim, +1));
		}
		if (constant_thermal_pressure) {
			const auto vel = (Mom(*cell.data)/Mas(*cell.data)).eval();
			Nrj(*cell.data) = pamhd::mhd::get_total_energy_density(
				Mas(*cell.data), vel, old_pressure, VMag(*cell.data),
				adiabatic_index, vacuum_permeability
			);
		}
	}
} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}

/*! Converts substepping periods from 2^N to N format.

Returns new largest N.
*/
template <
	class Grid,
	class Solver_Info_Getter,
	class Substepping_Period_Getter
> int update_substeps(
	Grid& grid,
	const Solver_Info_Getter SInfo,
	const Substepping_Period_Getter Substep
) try {
	int max_local = -1;
	for (const auto& cell: grid.local_cells()) {
		if (SInfo(*cell.data) < 0) {
			continue;
		}
		Substep(*cell.data) = 1 << Substep(*cell.data);
		max_local = std::max(Substep(*cell.data), max_local);
	}
	int max_global = -1;
	auto comm = grid.get_communicator();
	if (
		MPI_Allreduce(
			&max_local,
			&max_global,
			1,
			MPI_INT,
			MPI_MAX,
			comm
		) != MPI_SUCCESS
	) {
		std::cerr << __FILE__ "(" << __LINE__
			<< "): Couldn't reduce max substep." << std::endl;
		abort();
	}
	MPI_Comm_free(&comm);

	return max_global;

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
	const Mass_Density_Getter Mas,
	const Momentum_Density_Getter Mom,
	const Total_Energy_Density_Getter Nrj,
	const Magnetic_Field_Getter Mag,
	const Face_Magnetic_Field_Getter Face_B,
	const Face_dB_Getter Face_dB,
	const Background_Magnetic_Field_Getter Bg_B,
	const Mass_Density_Flux_Getters Mas_fs,
	const Momentum_Density_Flux_Getters Mom_fs,
	const Total_Energy_Density_Flux_Getters Nrj_fs,
	const Magnetic_Field_Flux_Getters Mag_fs,
	const Solver_Info_Getter SInfo,
	const Timestep_Getter Timestep,
	const Substepping_Period_Getter Substep,
	const Substep_Min_Getter Substep_Min,
	const Substep_Max_Getter Substep_Max,
	const Max_Velocity_Getter Max_v
) try {
	using std::max;
	using std::min;

	set_minmax_substepping_period(
		simulation_time,
		grid,
		options_mhd,
		Substep_Min,
		Substep_Max
	);

	Grid::cell_data_type::set_transfer_all(true, pamhd::mhd::Max_Velocity());
	grid.update_copies_of_remote_neighbors();
	Grid::cell_data_type::set_transfer_all(false, pamhd::mhd::Max_Velocity());

	const double sub_dt = set_minmax_substepping_period(
		grid,
		max_time_step,
		time_step_factor,
		SInfo,
		Timestep,
		Substep_Min,
		Substep_Max,
		Max_v
	);

	Grid::cell_data_type::set_transfer_all(true, pamhd::mhd::Substepping_Period());
	restrict_substepping_period(
		grid,
		Substep,
		Substep_Max,
		SInfo
	);

	const int max_substep = update_substeps(grid, SInfo, Substep);
	grid.update_copies_of_remote_neighbors();
	Grid::cell_data_type::set_transfer_all(false, pamhd::mhd::Substepping_Period());
	if (grid.get_rank() == 0) {
		std::cout
			<< "Substep: " << sub_dt << ", largest substep period: "
			<< max_substep << std::endl;
	}

	double total_dt = 0;
	for (int substep = 1; substep <= max_substep; substep += 1) {
		total_dt += sub_dt;

		Grid::cell_data_type::set_transfer_all(true, pamhd::mhd::MHD_State_Conservative());
		grid.start_remote_neighbor_copy_updates();
		pamhd::mhd::get_fluxes(
			solver, grid.inner_cells(), grid, substep,
			adiabatic_index, vacuum_permeability, sub_dt,
			Mas, Mom, Nrj, Mag, Face_dB, Bg_B,
			Mas_fs, Mom_fs, Nrj_fs, Mag_fs,
			SInfo, Substep, Max_v
		);

		grid.wait_remote_neighbor_copy_update_receives();

		pamhd::mhd::get_fluxes(
			solver, grid.outer_cells(), grid, substep,
			adiabatic_index, vacuum_permeability, sub_dt,
			Mas, Mom, Nrj, Mag, Face_dB, Bg_B,
			Mas_fs, Mom_fs, Nrj_fs, Mag_fs,
			SInfo, Substep, Max_v
		);

		grid.wait_remote_neighbor_copy_update_sends();

		pamhd::mhd::get_fluxes(
			solver, grid.remote_cells(), grid, substep,
			adiabatic_index, vacuum_permeability, sub_dt,
			Mas, Mom, Nrj, Mag, Face_dB, Bg_B,
			Mas_fs, Mom_fs, Nrj_fs, Mag_fs,
			SInfo, Substep, Max_v
		);
		Grid::cell_data_type::set_transfer_all(false, pamhd::mhd::MHD_State_Conservative());

		pamhd::mhd::update_mhd_state(
			grid.local_cells(), grid,
			substep, Face_B, Face_dB, SInfo,
			Substep, Mas, Mom, Nrj, Mag,
			Mas_fs, Mom_fs, Nrj_fs, Mag_fs
		);

		Grid::cell_data_type::set_transfer_all(true,
			pamhd::Face_Magnetic_Field(),
			// update pressure for B consistency calculation
			pamhd::mhd::MHD_State_Conservative()
		);
		grid.update_copies_of_remote_neighbors();

		// constant thermal pressure when updating vol B after solution
		pamhd::mhd::update_B_consistency(
			substep, grid.local_cells(),
			Mas, Mom, Nrj, Mag, Face_B,
			SInfo, Substep,
			adiabatic_index,
			vacuum_permeability,
			true
		);

		grid.update_copies_of_remote_neighbors();
		Grid::cell_data_type::set_transfer_all(false,
			pamhd::Face_Magnetic_Field(),
			pamhd::mhd::MHD_State_Conservative()
		);

	}

	return total_dt;

} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


/*!
Reduces difference in max substep between neighbors to <= 1.

Assumes 2^N format for substep periods.
*/
template <
	class Grid,
	class Substepping_Period_Getter,
	class Substep_Max_Getter,
	class Solver_Info_Getter
> void restrict_substepping_period(
	Grid& grid,
	const Substepping_Period_Getter Substep,
	const Substep_Max_Getter Substep_Max,
	const Solver_Info_Getter SInfo
) try {
	using std::max;
	using std::min;

	for (const auto& cell: grid.local_cells()) {
		if (SInfo(*cell.data) < 0) {
			continue;
		}
		Substep(*cell.data) = Substep_Max(*cell.data);
	}

	uint64_t modified_cells = 0;
	auto comm = grid.get_communicator();
	do {
		grid.update_copies_of_remote_neighbors();
		uint64_t modified_cells_local = 0;
		for (const auto& cell: grid.local_cells()) {
			if (SInfo(*cell.data) < 0) {
				continue;
			}
			for (const auto& neighbor: cell.neighbors_of) {
				if (SInfo(*neighbor.data) < 0) {
					continue;
				}
				if (neighbor.face_neighbor == 0 and neighbor.edge_neighbor[0] < 0) {
					continue;
				}

				if (Substep(*cell.data) > Substep(*neighbor.data) + 1) {
					Substep(*cell.data) = Substep(*neighbor.data) + 1;
					modified_cells_local++;
				}
			}
		}
		if (
			MPI_Allreduce(
				&modified_cells_local,
				&modified_cells,
				1,
				MPI_UINT64_T,
				MPI_SUM,
				comm
			) != MPI_SUCCESS
		) {
			std::cerr << __FILE__ "(" << __LINE__
				<< "): Couldn't reduce modified cell count."
				<< std::endl;
			abort();
		}
	} while (modified_cells > 0);
	grid.update_copies_of_remote_neighbors();
	MPI_Comm_free(&comm);

} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


/*!
Sets min and max substepping periods based on geometry.

Offsets minimum substep so that at it will be 0 in at
least one cell.

Assumes that substep is stored in 2^N format.
*/
template <
	class Grid,
	class Substep_Min_Getter,
	class Substep_Max_Getter
> void set_minmax_substepping_period(
	double time,
	Grid& grid,
	Options& options_mhd,
	const Substep_Min_Getter Substep_Min,
	const Substep_Max_Getter Substep_Max
) try {
	int min_substep_min_local = std::numeric_limits<int>::max();
	for (const auto& cell: grid.local_cells()) {
		const auto c = grid.geometry.get_center(cell.id);
		const auto
			r = std::sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]),
			lat = std::asin(c[2] / r),
			lon = std::atan2(c[1], c[0]);
		if (options_mhd.substep_min_i >= 0) {
			Substep_Min(*cell.data) = options_mhd.substep_min_i;
		} else {
			Substep_Min(*cell.data) = options_mhd.substep_min_e.evaluate(time, c[0], c[1], c[2], r, lat, lon);
		}
		if (options_mhd.substep_max_i >= 0) {
			Substep_Max(*cell.data) = options_mhd.substep_max_i;
		} else {
			Substep_Max(*cell.data) = options_mhd.substep_max_e.evaluate(time, c[0], c[1], c[2], r, lat, lon);
		}
		min_substep_min_local = std::min(Substep_Min(*cell.data), min_substep_min_local);
	}

	int min_substep_min_global = std::numeric_limits<int>::max();
	auto comm = grid.get_communicator();
	if (
		MPI_Allreduce(
			&min_substep_min_local,
			&min_substep_min_global,
			1,
			MPI_INT,
			MPI_MIN,
			comm
		) != MPI_SUCCESS
	) {
		std::cerr << __FILE__ "(" << __LINE__
			<< "): Couldn't reduce min_substep_min."
			<< std::endl;
		abort();
	}
	MPI_Comm_free(&comm);

	if (min_substep_min_global > 0) {
		for (const auto& cell: grid.local_cells()) {
			Substep_Min(*cell.data) -= min_substep_min_global;
		}
	}
} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


/*! Restricts min and max substepping periods based on physics.

Assumes substep is in 2^N format indicating period of solving
a cell, or in other words, indicating length of cell's timestep
in substeps.

If max_dt is smaller than smallest allowed timestep (including
dt_factor) in any cell, returns max_dt and sets substep range
to [0, 0] in all cells.
Otherwise returns largest substep length satisfying each cell's
combination of timestep and minimum substep period.
For example if one cell's timestep is dt and minimum substep is
0 and another cell's timestep is 1.5*dt and minimum substep is
1, returns 1.5*dt/2.

Should be called after geometry version of this function.
*/
template <
	class Grid,
	class Solver_Info_Getter,
	class Timestep_Getter,
	class Substep_Min_Getter,
	class Substep_Max_Getter,
	class Max_Velocity_Getter
> double set_minmax_substepping_period(
	Grid& grid,
	const double max_dt,
	const double dt_factor,
	const Solver_Info_Getter SInfo,
	const Timestep_Getter Timestep,
	const Substep_Min_Getter Substep_Min,
	const Substep_Max_Getter Substep_Max,
	const Max_Velocity_Getter Max_v
) try {
	using std::cerr;
	using std::endl;
	using std::max;
	using std::min;
	using std::numeric_limits;

	// calculate range of dts allowed in cells
	double
		smallest_dt_local = numeric_limits<double>::max(),
		largest_dt_local = 0;

	for (const auto& cell: grid.local_cells()) {
		Timestep(*cell.data) = numeric_limits<double>::max();
	}
	for (const auto& cell: grid.all_cells()) {
		if (SInfo(*cell.data) < 0) {
			continue;
		}

		const auto cell_len = grid.geometry.get_length(cell.id);
		for (const auto& neighbor: cell.neighbors_of) {
			const auto& fn = neighbor.face_neighbor;
			if (
				fn <= 0
				or neighbor.data == nullptr
				or SInfo(*neighbor.data) < 0
			) {
				continue;
			}

			const auto neigh_dim = size_t(abs(fn) - 1);
			const auto neigh_len = grid.geometry.get_length(neighbor.id);
			const auto max_v = max(
				Max_v(*cell.data)(fn),
				Max_v(*neighbor.data)(-fn));

			if (cell.is_local) {
				if (max_v <= 0) {
					Timestep(*cell.data) = 0;
				} else {
					Timestep(*cell.data) = min(
						cell_len[neigh_dim] / max_v,
						Timestep(*cell.data));
				}
			}
			if (neighbor.is_local) {
				if (max_v <= 0) {
					Timestep(*neighbor.data) = 0;
				} else {
					Timestep(*neighbor.data) = min(
						neigh_len[neigh_dim] / max_v,
						Timestep(*neighbor.data));
				}
			}
		}
	}
	for (const auto& cell: grid.local_cells()) {
		if (SInfo(*cell.data) < 0) {
			continue;
		}
		Timestep(*cell.data) *= dt_factor;
		smallest_dt_local = min(Timestep(*cell.data), smallest_dt_local);
		largest_dt_local = max(Timestep(*cell.data), largest_dt_local);
	}
	auto comm = grid.get_communicator();
	double smallest_dt_global = -1, largest_dt_global = -1;
	if (
		MPI_Allreduce(
			&smallest_dt_local,
			&smallest_dt_global,
			1,
			MPI_DOUBLE,
			MPI_MIN,
			comm
		) != MPI_SUCCESS
	) {
		cerr << __FILE__ "(" << __LINE__
			<< "): Couldn't reduce smallest_dt." << endl;
		abort();
	}
	if (
		MPI_Allreduce(
			&largest_dt_local,
			&largest_dt_global,
			1,
			MPI_DOUBLE,
			MPI_MIN,
			comm
		) != MPI_SUCCESS
	) {
		cerr << __FILE__ "(" << __LINE__
			<< "): Couldn't reduce largest_dt." << endl;
		abort();
	}

	if (max_dt <= smallest_dt_global) {
		for (const auto& cell: grid.local_cells()) {
			Substep_Min(*cell.data) =
			Substep_Max(*cell.data) = 0;
		}
		MPI_Comm_free(&comm);
		return max_dt;
	}

	// minimum substep length
	double ret_val_local = numeric_limits<double>::max();
	for (const auto& cell: grid.local_cells()) {
		if (SInfo(*cell.data) < 0) {
			continue;
		}
		ret_val_local = min(ret_val_local,
			Timestep(*cell.data) / (1 << Substep_Min(*cell.data)));
	}

	double ret_val_global = numeric_limits<double>::max();
	if (
		MPI_Allreduce(
			&ret_val_local,
			&ret_val_global,
			1,
			MPI_DOUBLE,
			MPI_MIN,
			comm
		) != MPI_SUCCESS
	) {
		cerr << __FILE__ "(" << __LINE__
			<< "): Couldn't reduce ret_val." << endl;
		abort();
	}
	MPI_Comm_free(&comm);

	// decrease too large max substep periods
	for (const auto& cell: grid.local_cells()) {
		if (SInfo(*cell.data) < 0) {
			continue;
		}

		const auto& dt = Timestep(*cell.data);
		if (dt > 0) {
			int step = max(0, (int)std::floor(std::log2(dt / ret_val_global)));
			if (Substep_Max(*cell.data) > step) {
				Substep_Max(*cell.data) = step;
			}
			if (Substep_Min(*cell.data) > Substep_Max(*cell.data)) {
				cerr << "Unexpected substeps in cell "
					<< cell.id << ": "
					<< Substep_Min(*cell.data) << ", "
					<< Substep_Max(*cell.data) << ", "
					<< dt << ", " << ret_val_global << endl;
				abort();
			}
		} else {
			Substep_Min(*cell.data) =
			Substep_Max(*cell.data) = 0;
		}
	}

	return ret_val_global;

} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


}} // namespaces


#endif // ifndef PAMHD_MHD_SOLVE_STAGGERED_HPP
