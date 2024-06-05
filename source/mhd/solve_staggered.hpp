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
#include "mhd/roe_athena.hpp"
#include "mhd/variables.hpp"
#include "variables.hpp"


namespace pamhd {
namespace mhd {


/*!
Calculates MHD fluxes in/out of given cells.

Returns the maximum allowed length of time step for the next step on this process.

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
	class Background_Magnetic_Field_Getter,
	class Mass_Density_Flux_Getters,
	class Momentum_Density_Flux_Getters,
	class Total_Energy_Density_Flux_Getters,
	class Magnetic_Field_Flux_Getters,
	class Primary_Face_Getter,
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
	const Mass_Density_Getter Mas,
	const Momentum_Density_Getter Mom,
	const Total_Energy_Density_Getter Nrj,
	const Magnetic_Field_Getter Mag,
	const Background_Magnetic_Field_Getter Bg_B,
	const Mass_Density_Flux_Getters Mas_f,
	const Momentum_Density_Flux_Getters Mom_f,
	const Total_Energy_Density_Flux_Getters Nrj_f,
	const Magnetic_Field_Flux_Getters Mag_f,
	const Primary_Face_Getter PFace,
	const Solver_Info_Getter SInfo,
	const Substepping_Period_Getter Substep,
	const Max_Velocity_Getter Max_v
) try {
	using std::abs;
	using std::get;
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

	// shorthand for referring to variables of internal MHD data type
	const pamhd::mhd::Mass_Density mas_int{};
	const pamhd::mhd::Momentum_Density mom_int{};
	const pamhd::mhd::Total_Energy_Density nrj_int{};
	const pamhd::Magnetic_Field mag_int{};

	for (const auto& cell: cells) {
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

		// solve flux also if boundary cell on negative side of face
		if (SInfo(*cell.data) < 0) {
			continue;
		}

		if (current_substep % Substep(*cell.data) != 0) {
			continue;
		}

		const auto& cpface = PFace(*cell.data);

		// maximum substep period in neighborhood
		int max_substep = Substep(*cell.data);
		for (const auto& neighbor: cell.neighbors_of) {
			if (neighbor.face_neighbor == 0 and neighbor.edge_neighbor[0] < 0) {
				continue;
			}
			if (SInfo(*neighbor.data) < 0) {
				continue;
			}
			max_substep = std::max(Substep(*neighbor.data), max_substep);
		}

		for (const auto& neighbor: cell.neighbors_of) {
			const auto& n = neighbor.face_neighbor;
			// flux stored in neighbor and calculated by its owner
			if (n == 0 or not cpface(n)) continue;

			if (SInfo(*neighbor.data) < 0) {
				continue;
			}

			if (current_substep % Substep(*neighbor.data) != 0) {
				continue;
			}

			detail::MHD state_neg, state_pos;
			// take into account direction of neighbor from cell
			state_neg[mas_int] = Mas(*cell.data);
			state_neg[mom_int] = get_rotated_vector(Mom(*cell.data), abs(n));
			state_neg[nrj_int] = Nrj(*cell.data);
			state_neg[mag_int] = get_rotated_vector(Mag(*cell.data), abs(n));

			state_pos[mas_int] = Mas(*neighbor.data);
			state_pos[mom_int] = get_rotated_vector(Mom(*neighbor.data), abs(n));
			state_pos[nrj_int] = Nrj(*neighbor.data);
			state_pos[mag_int] = get_rotated_vector(Mag(*neighbor.data), abs(n));
			if (n < 0) {
				const auto temp = state_neg;
				state_neg = state_pos;
				state_pos = temp;
			}

			const Magnetic_Field::data_type bg_face_b
				= get_rotated_vector(Bg_B(*cell.data)(n), abs(n));
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
					<< "Solution failed between cells " << cell.id
					<< " and " << neighbor.id
					<< " of cell type " << SInfo(*cell.data)
					<< " and " << SInfo(*neighbor.data)
					<< " at " << grid.geometry.get_center(cell.id)
					<< " and " << grid.geometry.get_center(neighbor.id)
					<< " in direction " << n
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

			Max_v(*cell.data)(n) = max_vel;

			// rotate flux back
			flux[mom_int] = get_rotated_vector(flux[mom_int], -abs(n));
			flux[mag_int] = get_rotated_vector(flux[mag_int], -abs(n));

			if (n == -1) {
				Mas_fnx(*cell.data) = flux[mas_int];
				Mom_fnx(*cell.data) = flux[mom_int];
				Nrj_fnx(*cell.data) = flux[nrj_int];
				Mag_fnx(*cell.data) = flux[mag_int];
			}
			if (n == +1) {
				Mas_fpx(*cell.data) = flux[mas_int];
				Mom_fpx(*cell.data) = flux[mom_int];
				Nrj_fpx(*cell.data) = flux[nrj_int];
				Mag_fpx(*cell.data) = flux[mag_int];
			}
			if (n == -2) {
				Mas_fny(*cell.data) = flux[mas_int];
				Mom_fny(*cell.data) = flux[mom_int];
				Nrj_fny(*cell.data) = flux[nrj_int];
				Mag_fny(*cell.data) = flux[mag_int];
			}
			if (n == +2) {
				Mas_fpy(*cell.data) = flux[mas_int];
				Mom_fpy(*cell.data) = flux[mom_int];
				Nrj_fpy(*cell.data) = flux[nrj_int];
				Mag_fpy(*cell.data) = flux[mag_int];
			}
			if (n == -3) {
				Mas_fnz(*cell.data) = flux[mas_int];
				Mom_fnz(*cell.data) = flux[mom_int];
				Nrj_fnz(*cell.data) = flux[nrj_int];
				Mag_fnz(*cell.data) = flux[mag_int];
			}
			if (n == +3) {
				Mas_fpz(*cell.data) = flux[mas_int];
				Mom_fpz(*cell.data) = flux[mom_int];
				Nrj_fpz(*cell.data) = flux[nrj_int];
				Mag_fpz(*cell.data) = flux[mag_int];
			}
		}
	}
} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


/*!
Calculates edge electric fields in given cells.

Also updates MHD state with fluxes.

Saves edge E of cells with SInfo(*cell_data) == 1,
ignores cells with SInfo < 0.
*/
template <
	class Grid,
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Magnetic_Field_Getter,
	class Edge_Electric_Field_Getter,
	class Mass_Density_Flux_Getters,
	class Momentum_Density_Flux_Getters,
	class Total_Energy_Density_Flux_Getters,
	class Magnetic_Field_Flux_Getters,
	class Primary_Face_Getter,
	class Primary_Edge_Getter,
	class Face_Info_Getter,
	class Substepping_Period_Getter
> void get_edge_electric_field(
	Grid& grid,
	const double sub_dt,
	const int current_substep,
	const Mass_Density_Getter Mas,
	const Momentum_Density_Getter Mom,
	const Total_Energy_Density_Getter Nrj,
	const Magnetic_Field_Getter Mag,
	const Edge_Electric_Field_Getter Edge_E,
	const Mass_Density_Flux_Getters Mas_f,
	const Momentum_Density_Flux_Getters Mom_f,
	const Total_Energy_Density_Flux_Getters Nrj_f,
	const Magnetic_Field_Flux_Getters Mag_f,
	const Primary_Face_Getter PFace,
	const Primary_Edge_Getter PEdge,
	const Face_Info_Getter FInfo,
	const Substepping_Period_Getter Substep
) try {
	using std::array;
	using std::get;
	using std::make_tuple;
	using std::min;
	using std::pow;
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

	for (const auto& cell: grid.local_cells()) {
		if (current_substep % Substep(*cell.data) == 0) {
			continue;
		}

		int max_substep = Substep(*cell.data);
		for (const auto& neighbor: cell.neighbors_of) {
			const auto& fn = neighbor.face_neighbor;
			if (fn == 0 and neighbor.edge_neighbor[0] < 0) {
				continue;
			}
			if (fn != 0 and PFace(*neighbor.data)(-fn) and FInfo(*neighbor.data)(-fn) < 0) {
				continue;
			}
			max_substep = std::max(Substep(*neighbor.data), max_substep);
		}

		auto& edge_e = Edge_E(*cell.data);
		// track number of contributions to edge_e()s
		grid::Edge_Type<int> e_items;
		for (int d1: {0, 1, 2}) // parallel dim of edge
		for (int d2: {-1, +1}) // side in 1st perpendicular dim of edge
		for (int d3: {-1, +1}) { // side of cell in 2nd perp dim
			edge_e(d1,d2,d3) = 0.0;
			e_items(d1,d2,d3) = 0;
		}

		const auto& cpface = PFace(*cell.data); // primary face true/false
		const auto& cfinfo = FInfo(*cell.data); // face type -1/0/+1
		const auto& cpedge = PEdge(*cell.data); // primary edge true/false

		/*! Contributions to edge electric fields

		Edge directed along x:
			Flux through touching face with normal in y direction:
				-Bz flux is added to E
			Flux through face with normal in z direction:
				+By flux is added to E
		Edge y:
			Flux x: +Bz flux added
			Flux z: -Bx flux added
		Edge z:
			Flux x: -By
			Flux y: +Bx
		*/

		// cell's own face can contribute to all touching edges
		size_t Bd = 9;
		int d1 = 9, d2 = 9, d3 = 9;
		if (cpface(-1) and cfinfo(-1) >= 0) {
			if (d1 = 1, d2 = -1, d3 = -1, Bd = d1 + 1;
				// non-primary edges are handled by another cell
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) += Mag_fnx(*cell.data)[Bd];
				e_items(d1,d2,d3)++;
			}
			if (d1 = 1, d2 = -1, d3 = +1, Bd = d1 + 1;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) += Mag_fnx(*cell.data)[Bd];
				e_items(d1,d2,d3)++;
			}
			if (d1 = 2, d2 = -1, d3 = -1, Bd = d1 - 1;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) -= Mag_fnx(*cell.data)[Bd];
				e_items(d1,d2,d3)++;
			}
			if (d1 = 2, d2 = -1, d3 = +1, Bd = d1 - 1;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) -= Mag_fnx(*cell.data)[Bd];
				e_items(d1,d2,d3)++;
			}
		}

		if (cpface(+1) and cfinfo(+1) >= 0) {
			if (d1 = 1, d2 = +1, d3 = -1, Bd = d1 + 1;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) += Mag_fpx(*cell.data)[Bd];
				e_items(d1,d2,d3)++;
			}
			if (d1 = 1, d2 = +1, d3 = +1, Bd = d1 + 1;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) += Mag_fpx(*cell.data)[Bd];
				e_items(d1,d2,d3)++;
			}
			if (d1 = 2, d2 = +1, d3 = -1, Bd = d1 - 1;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) -= Mag_fpx(*cell.data)[Bd];
				e_items(d1,d2,d3)++;
			}
			if (d1 = 2, d2 = +1, d3 = +1, Bd = d1 - 1;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) -= Mag_fpx(*cell.data)[Bd];
				e_items(d1,d2,d3)++;
			}
		}

		if (cpface(-2) and cfinfo(-2) >= 0) {
			if (d1 = 0, d2 = -1, d3 = -1, Bd = d1 + 2;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) -= Mag_fny(*cell.data)[Bd];
				e_items(d1,d2,d3)++;
			}
			if (d1 = 0, d2 = -1, d3 = +1, Bd = d1 + 2;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) -= Mag_fny(*cell.data)[Bd];
				e_items(d1,d2,d3)++;
			}
			if (d1 = 2, d2 = -1, d3 = -1, Bd = d1 - 2;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) += Mag_fny(*cell.data)[Bd];
				e_items(d1,d2,d3)++;
			}
			if (d1 = 2, d2 = +1, d3 = -1, Bd = d1 - 2;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) += Mag_fny(*cell.data)[Bd];
				e_items(d1,d2,d3)++;
			}
		}

		if (cpface(+2) and cfinfo(+2) >= 0) {
			if (d1 = 0, d2 = +1, d3 = -1, Bd = d1 + 2;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) -= Mag_fpy(*cell.data)[Bd];
				e_items(d1,d2,d3)++;
			}
			if (d1 = 0, d2 = +1, d3 = +1, Bd = d1 + 2;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) -= Mag_fpy(*cell.data)[Bd];
				e_items(d1,d2,d3)++;
			}
			if (d1 = 2, d2 = -1, d3 = +1, Bd = d1 - 2;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) += Mag_fpy(*cell.data)[Bd];
				e_items(d1,d2,d3)++;
			}
			if (d1 = 2, d2 = +1, d3 = +1, Bd = d1 - 2;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) += Mag_fpy(*cell.data)[Bd];
				e_items(d1,d2,d3)++;
			}
		}

		if (cpface(-3) and cfinfo(-3) >= 0) {
			if (d1 = 0, d2 = -1, d3 = -1, Bd = d1 + 1;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) += Mag_fnz(*cell.data)[Bd];
				e_items(d1,d2,d3)++;
			}
			if (d1 = 0, d2 = +1, d3 = -1, Bd = d1 + 1;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) += Mag_fnz(*cell.data)[Bd];
				e_items(d1,d2,d3)++;
			}
			if (d1 = 1, d2 = -1, d3 = -1, Bd = d1 - 1;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) -= Mag_fnz(*cell.data)[Bd];
				e_items(d1,d2,d3)++;
			}
			if (d1 = 1, d2 = +1, d3 = -1, Bd = d1 - 1;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) -= Mag_fnz(*cell.data)[Bd];
				e_items(d1,d2,d3)++;
			}
		}

		if (cpface(+3) and cfinfo(+3) >= 0) {
			if (d1 = 0, d2 = -1, d3 = +1, Bd = d1 + 1;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) += Mag_fpz(*cell.data)[Bd];
				e_items(d1,d2,d3)++;
			}
			if (d1 = 0, d2 = +1, d3 = +1, Bd = d1 + 1;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) += Mag_fpz(*cell.data)[Bd];
				e_items(d1,d2,d3)++;
			}
			if (d1 = 1, d2 = -1, d3 = +1, Bd = d1 - 1;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) -= Mag_fpz(*cell.data)[Bd];
				e_items(d1,d2,d3)++;
			}
			if (d1 = 1, d2 = +1, d3 = +1, Bd = d1 - 1;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) -= Mag_fpz(*cell.data)[Bd];
				e_items(d1,d2,d3)++;
			}
		}

		const auto [dx, dy, dz] = grid.geometry.get_length(cell.id);
		const auto cleni = grid.mapping.get_cell_length_in_indices(cell.id);

		for (const auto& neighbor: cell.neighbors_of) {
			const auto& fn = neighbor.face_neighbor;
			const auto& en = neighbor.edge_neighbor;
			if (fn == 0 and en[0] < 0) {
				continue;
			}

			if (current_substep % Substep(*neighbor.data) != 0) {
				continue;
			}

			const auto dt = sub_dt * std::max(Substep(*cell.data), Substep(*neighbor.data));
			const auto [dtdx, dtdy, dtdz] = [&]{
				if (neighbor.relative_size > 0) {
					return make_tuple(dt/dx/4, dt/dy/4, dt/dz/4);
				} else {
					return make_tuple(dt/dx, dt/dy, dt/dz);
				}
			}();

			/*
			Apply fluxes
			*/

			if (fn == -1) {
				if (cpface(fn)) {
					Mas(*cell.data) += Mas_fnx(*cell.data)*dtdx;
					Mom(*cell.data) += Mom_fnx(*cell.data)*dtdx;
					Nrj(*cell.data) += Nrj_fnx(*cell.data)*dtdx;
					Mag(*cell.data) += Mag_fnx(*cell.data)*dtdx;
				} else {
					Mas(*cell.data) += Mas_fpx(*neighbor.data)*dtdx;
					Mom(*cell.data) += Mom_fpx(*neighbor.data)*dtdx;
					Nrj(*cell.data) += Nrj_fpx(*neighbor.data)*dtdx;
					Mag(*cell.data) += Mag_fpx(*neighbor.data)*dtdx;
				}
			}

			if (fn == +1) {
				if (cpface(fn)) {
					Mas(*cell.data) -= Mas_fpx(*cell.data)*dtdx;
					Mom(*cell.data) -= Mom_fpx(*cell.data)*dtdx;
					Nrj(*cell.data) -= Nrj_fpx(*cell.data)*dtdx;
					Mag(*cell.data) -= Mag_fpx(*cell.data)*dtdx;
				} else {
					Mas(*cell.data) -= Mas_fnx(*neighbor.data)*dtdx;
					Mom(*cell.data) -= Mom_fnx(*neighbor.data)*dtdx;
					Nrj(*cell.data) -= Nrj_fnx(*neighbor.data)*dtdx;
					Mag(*cell.data) -= Mag_fnx(*neighbor.data)*dtdx;
				}
			}

			if (fn == -2) {
				if (cpface(fn)) {
					Mas(*cell.data) += Mas_fny(*cell.data)*dtdy;
					Mom(*cell.data) += Mom_fny(*cell.data)*dtdy;
					Nrj(*cell.data) += Nrj_fny(*cell.data)*dtdy;
					Mag(*cell.data) += Mag_fny(*cell.data)*dtdy;
				} else {
					Mas(*cell.data) += Mas_fpy(*neighbor.data)*dtdy;
					Mom(*cell.data) += Mom_fpy(*neighbor.data)*dtdy;
					Nrj(*cell.data) += Nrj_fpy(*neighbor.data)*dtdy;
					Mag(*cell.data) += Mag_fpy(*neighbor.data)*dtdy;
				}
			}

			if (fn == +2) {
				if (cpface(fn)) {
					Mas(*cell.data) -= Mas_fpy(*cell.data)*dtdy;
					Mom(*cell.data) -= Mom_fpy(*cell.data)*dtdy;
					Nrj(*cell.data) -= Nrj_fpy(*cell.data)*dtdy;
					Mag(*cell.data) -= Mag_fpy(*cell.data)*dtdy;
				} else {
					Mas(*cell.data) -= Mas_fny(*neighbor.data)*dtdy;
					Mom(*cell.data) -= Mom_fny(*neighbor.data)*dtdy;
					Nrj(*cell.data) -= Nrj_fny(*neighbor.data)*dtdy;
					Mag(*cell.data) -= Mag_fny(*neighbor.data)*dtdy;
				}
			}

			if (fn == -3) {
				if (cpface(fn)) {
					Mas(*cell.data) += Mas_fnz(*cell.data)*dtdz;
					Mom(*cell.data) += Mom_fnz(*cell.data)*dtdz;
					Nrj(*cell.data) += Nrj_fnz(*cell.data)*dtdz;
					Mag(*cell.data) += Mag_fnz(*cell.data)*dtdz;
				} else {
					Mas(*cell.data) += Mas_fpz(*neighbor.data)*dtdz;
					Mom(*cell.data) += Mom_fpz(*neighbor.data)*dtdz;
					Nrj(*cell.data) += Nrj_fpz(*neighbor.data)*dtdz;
					Mag(*cell.data) += Mag_fpz(*neighbor.data)*dtdz;
				}
			}

			if (fn == +3) {
				if (cpface(fn)) {
					Mas(*cell.data) -= Mas_fpz(*cell.data)*dtdz;
					Mom(*cell.data) -= Mom_fpz(*cell.data)*dtdz;
					Nrj(*cell.data) -= Nrj_fpz(*cell.data)*dtdz;
					Mag(*cell.data) -= Mag_fpz(*cell.data)*dtdz;
				} else {
					Mas(*cell.data) -= Mas_fnz(*neighbor.data)*dtdz;
					Mom(*cell.data) -= Mom_fnz(*neighbor.data)*dtdz;
					Nrj(*cell.data) -= Nrj_fnz(*neighbor.data)*dtdz;
					Mag(*cell.data) -= Mag_fnz(*neighbor.data)*dtdz;
				}
			}

			/*
			Contributions to edge electric fields
			*/
			const auto& npface = PFace(*neighbor.data);
			const auto& nfinfo = FInfo(*neighbor.data);
			const auto nleni = grid.mapping.get_cell_length_in_indices(neighbor.id);

			if (fn == -1) {
				if (d1 = 1, d2 = -1, d3 = -1, Bd = 0; cpedge(d1,d2,d3)) {
					if (npface(-3) and nfinfo(-3) >= 0 and neighbor.z == 0) {
						edge_e(d1,d2,d3) -= Mag_fnz(*neighbor.data)[Bd];
						e_items(d1,d2,d3)++;
					}
					if (cleni == neighbor.z + nleni) {
						throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					}
				}

				if (d1 = 1, d2 = -1, d3 = +1, Bd = 0; cpedge(d1,d2,d3)) {
					if (npface(+3) and nfinfo(+3) >= 0 and cleni == neighbor.z + nleni) {
						edge_e(d1,d2,d3) -= Mag_fpz(*neighbor.data)[Bd];
						e_items(d1,d2,d3)++;
					}
					/*
					Substitute flux from missing face in larger neighbor
					with one on opposite side of edge, in this cell.
					*/
					if (
						neighbor.relative_size < 0
						and cpface(+3)
						and cfinfo(+3) >= 0
						and neighbor.z == 0
					) {
						edge_e(d1,d2,d3) -= Mag_fpz(*cell.data)[Bd];
						e_items(d1,d2,d3)++;
					}
					if (neighbor.relative_size >= 0 and neighbor.z == 0) {
						throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					}
				}

				if (d1 = 2, d2 = -1, d3 = -1, Bd = 0; cpedge(d1,d2,d3)) {
					if (npface(-2) and nfinfo(-2) >= 0 and neighbor.y == 0) {
						edge_e(d1,d2,d3) += Mag_fny(*neighbor.data)[Bd];
						e_items(d1,d2,d3)++;
					}
					if (cleni == neighbor.y + nleni) {
						throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					}
				}

				if (d1 = 2, d2 = -1, d3 = +1, Bd = 0; cpedge(d1,d2,d3)) {
					if (npface(+2) and nfinfo(+2) >= 0 and cleni == neighbor.y + nleni) {
						edge_e(d1,d2,d3) += Mag_fpy(*neighbor.data)[Bd];
						e_items(d1,d2,d3)++;
					}
					if (
						neighbor.relative_size < 0
						and cfinfo(+2) >= 0
						and cpface(+2)
						and neighbor.y == 0
					) {
						edge_e(d1,d2,d3) += Mag_fpy(*cell.data)[Bd];
						e_items(d1,d2,d3)++;
					}
					if (neighbor.relative_size >= 0 and neighbor.y == 0) {
						throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					}
				}
			}

			if (fn == +1) {
				if (d1 = 2, d2 = +1, d3 = -1, Bd = 0; cpedge(d1,d2,d3)) {
					if (npface(-2) and nfinfo(-2) >= 0 and neighbor.y == 0) {
						edge_e(d1,d2,d3) += Mag_fny(*neighbor.data)[Bd];
						e_items(d1,d2,d3)++;
					}
				}

				if (d1 = 2, d2 = +1, d3 = +1, Bd = 0; cpedge(d1,d2,d3)) {
					if (npface(+2) and nfinfo(+2) >= 0 and cleni == neighbor.y + nleni) {
						edge_e(d1,d2,d3) += Mag_fpy(*neighbor.data)[Bd];
						e_items(d1,d2,d3)++;
					}
					if (
						neighbor.relative_size < 0
						and cpface(+2)
						and cfinfo(+2) >= 0
						and neighbor.y == 0
					) {
						edge_e(d1,d2,d3) += Mag_fpy(*cell.data)[Bd];
						e_items(d1,d2,d3)++;
					}
				}

				if (d1 = 1, d2 = +1, d3 = -1, Bd = 0; cpedge(d1,d2,d3)) {
					if (npface(-3) and nfinfo(-3) >= 0 and neighbor.z == 0) {
						edge_e(d1,d2,d3) -= Mag_fnz(*neighbor.data)[Bd];
						e_items(d1,d2,d3)++;
					}
				}

				if (d1 = 1, d2 = +1, d3 = +1, Bd = 0; cpedge(d1,d2,d3)) {
					if (npface(+3) and nfinfo(+3) >= 0 and cleni == neighbor.z + nleni) {
						edge_e(d1,d2,d3) -= Mag_fpz(*neighbor.data)[Bd];
						e_items(d1,d2,d3)++;
					}
					if (
						neighbor.relative_size < 0
						and cpface(+3)
						and cfinfo(+3) >= 0
						and neighbor.z == 0
					) {
						edge_e(d1,d2,d3) -= Mag_fpz(*cell.data)[Bd];
						e_items(d1,d2,d3)++;
					}
				}
			}

			if (fn == -2) {
				if (d1 = 2, d2 = -1, d3 = -1, Bd = 1; cpedge(d1,d2,d3)) {
					if (npface(-1) and nfinfo(-1) >= 0 and neighbor.x == 0) {
						edge_e(d1,d2,d3) -= Mag_fnx(*neighbor.data)[Bd];
						e_items(d1,d2,d3)++;
					}
				}

				if (d1 = 2, d2 = +1, d3 = -1, Bd = 1; cpedge(d1,d2,d3)) {
					if (nfinfo(+1) >= 0 and npface(+1) and cleni == neighbor.x + nleni) {
						edge_e(d1,d2,d3) -= Mag_fpx(*neighbor.data)[Bd];
						e_items(d1,d2,d3)++;
					}
					if (
						neighbor.relative_size < 0
						and cpface(+1)
						and cfinfo(+1) >= 0
						and neighbor.x == 0
					) {
						edge_e(d1,d2,d3) -= Mag_fpx(*cell.data)[Bd];
						e_items(d1,d2,d3)++;
					}
				}

				if (d1 = 0, d2 = -1, d3 = -1, Bd = 1; cpedge(d1,d2,d3)) {
					if (npface(-3) and nfinfo(-3) >= 0 and neighbor.z == 0) {
						edge_e(d1,d2,d3) += Mag_fnz(*neighbor.data)[Bd];
						e_items(d1,d2,d3)++;
					}
				}

				if (d1 = 0, d2 = -1, d3 = +1, Bd = 1; cpedge(d1,d2,d3)) {
					if (npface(+3) and nfinfo(+3) >= 0 and cleni == neighbor.z + nleni) {
						edge_e(d1,d2,d3) += Mag_fpz(*neighbor.data)[Bd];
						e_items(d1,d2,d3)++;
					}
					if (
						neighbor.relative_size < 0
						and cpface(+3)
						and cfinfo(+3) >= 0
						and neighbor.z == 0
					) {
						edge_e(d1,d2,d3) += Mag_fpz(*cell.data)[Bd];
						e_items(d1,d2,d3)++;
					}
				}
			}

			if (fn == +2) {
				if (d1 = 2, d2 = -1, d3 = +1, Bd = 1; cpedge(d1,d2,d3)) {
					if (nfinfo(-1) >= 0 and npface(-1) and neighbor.x == 0) {
						edge_e(d1,d2,d3) -= Mag_fnx(*neighbor.data)[Bd];
						e_items(d1,d2,d3)++;
					}
				}

				if (d1 = 2, d2 = +1, d3 = +1, Bd = 1; cpedge(d1,d2,d3)) {
					if (npface(+1) and nfinfo(+1) >= 0 and cleni == neighbor.x + nleni) {
						edge_e(d1,d2,d3) -= Mag_fpx(*neighbor.data)[Bd];
						e_items(d1,d2,d3)++;
					}
					if (
						neighbor.relative_size < 0
						and cpface(+1)
						and cfinfo(+1) >= 0
						and neighbor.x == 0
					) {
						edge_e(d1,d2,d3) -= Mag_fpx(*cell.data)[Bd];
						e_items(d1,d2,d3)++;
					}
				}

				if (d1 = 0, d2 = +1, d3 = -1, Bd = 1; cpedge(d1,d2,d3)) {
					if (nfinfo(-3) >= 0 and npface(-3) and neighbor.z == 0) {
						edge_e(d1,d2,d3) += Mag_fnz(*neighbor.data)[Bd];
						e_items(d1,d2,d3)++;
					}
				}

				if (d1 = 0, d2 = +1, d3 = +1, Bd = 1; cpedge(d1,d2,d3)) {
					if (npface(+3) and nfinfo(+3) >= 0 and cleni == neighbor.z + nleni) {
						edge_e(d1,d2,d3) += Mag_fpz(*neighbor.data)[Bd];
						e_items(d1,d2,d3)++;
					}
					if (
						neighbor.relative_size < 0
						and cpface(+3)
						and cfinfo(+3) >= 0
						and neighbor.z == 0
					) {
						edge_e(d1,d2,d3) += Mag_fpz(*cell.data)[Bd];
						e_items(d1,d2,d3)++;
					}
				}
			}

			if (fn == -3) {
				if (d1 = 1, d2 = -1, d3 = -1, Bd = 2; cpedge(d1,d2,d3)) {
					if (npface(-1) and nfinfo(-1) >= 0 and neighbor.x == 0) {
						edge_e(d1,d2,d3) += Mag_fnx(*neighbor.data)[Bd];
						e_items(d1,d2,d3)++;
					}
				}

				if (d1 = 1, d2 = +1, d3 = -1, Bd = 2; cpedge(d1,d2,d3)) {
					if (npface(+1) and nfinfo(+1) >= 0 and cleni == neighbor.x + nleni) {
						edge_e(d1,d2,d3) += Mag_fpx(*neighbor.data)[Bd];
						e_items(d1,d2,d3)++;
					}
					if (
						neighbor.relative_size < 0
						and cpface(+1)
						and cfinfo(+1) >= 0
						and neighbor.x == 0
					) {
						edge_e(d1,d2,d3) += Mag_fpx(*cell.data)[Bd];
						e_items(d1,d2,d3)++;
					}
				}

				if (d1 = 0, d2 = -1, d3 = -1, Bd = 2; cpedge(d1,d2,d3)) {
					if (npface(-2) and nfinfo(-2) >= 0 and neighbor.y == 0) {
						edge_e(d1,d2,d3) -= Mag_fny(*neighbor.data)[Bd];
						e_items(d1,d2,d3)++;
					}
				}

				if (d1 = 0, d2 = +1, d3 = -1, Bd = 2; cpedge(d1,d2,d3)) {
					if (npface(+2) and nfinfo(+2) >= 0 and cleni == neighbor.y + nleni) {
						edge_e(d1,d2,d3) -= Mag_fpy(*neighbor.data)[Bd];
						e_items(d1,d2,d3)++;
					}
					if (
						neighbor.relative_size < 0
						and cpface(+2)
						and cfinfo(+2) >= 0
						and neighbor.y == 0
					) {
						edge_e(d1,d2,d3) -= Mag_fpy(*cell.data)[Bd];
						e_items(d1,d2,d3)++;
					}
				}
			}

			if (fn == +3) {
				if (d1 = 1, d2 = -1, d3 = +1, Bd = 2; cpedge(d1,d2,d3)) {
					if (npface(-1) and nfinfo(-1) >= 0 and neighbor.x == 0) {
						edge_e(d1,d2,d3) += Mag_fnx(*neighbor.data)[Bd];
						e_items(d1,d2,d3)++;
					}
				}

				if (d1 = 1, d2 = +1, d3 = +1, Bd = 2; cpedge(d1,d2,d3)) {
					if (npface(+1) and nfinfo(+1) >= 0 and cleni == neighbor.x + nleni) {
						edge_e(d1,d2,d3) += Mag_fpx(*neighbor.data)[Bd];
						e_items(d1,d2,d3)++;
					}
					if (
						neighbor.relative_size < 0
						and cpface(+1)
						and cfinfo(+1) >= 0
						and neighbor.x == 0
					) {
						edge_e(d1,d2,d3) += Mag_fpx(*cell.data)[Bd];
						e_items(d1,d2,d3)++;
					}
				}

				if (d1 = 0, d2 = -1, d3 = +1, Bd = 2; cpedge(d1,d2,d3)) {
					if (npface(-2) and nfinfo(-2) >= 0 and neighbor.y == 0) {
						edge_e(d1,d2,d3) -= Mag_fny(*neighbor.data)[Bd];
						e_items(d1,d2,d3)++;
					}
				}

				if (d1 = 0, d2 = +1, d3 = +1, Bd = 2; cpedge(d1,d2,d3)) {
					if (npface(+2) and nfinfo(+2) >= 0 and cleni == neighbor.y + nleni) {
						edge_e(d1,d2,d3) -= Mag_fpy(*neighbor.data)[Bd];
						e_items(d1,d2,d3)++;
					}
					if (
						neighbor.relative_size < 0
						and cpface(+2)
						and cfinfo(+2) >= 0
						and neighbor.y == 0
					) {
						edge_e(d1,d2,d3) -= Mag_fpy(*cell.data)[Bd];
						e_items(d1,d2,d3)++;
					}
				}
			}

			if (d1 = 0, d2 = -1, d3 = -1; cpedge(d1,d2,d3) and
				en[0] == d1 and en[1] == d2 and en[2] == d3
			) {
				if (Bd = 2; nfinfo(+2) >= 0 and npface(+2)) {
					edge_e(d1,d2,d3) -= Mag_fpy(*neighbor.data)[Bd];
					e_items(d1,d2,d3)++;
				}
				if (Bd = 1; nfinfo(+3) >= 0 and npface(+3)) {
					edge_e(d1,d2,d3) += Mag_fpz(*neighbor.data)[Bd];
					e_items(d1,d2,d3)++;
				}
			}

			if (d1 = 0, d2 = -1, d3 = +1; cpedge(d1,d2,d3) and
				en[0] == d1 and en[1] == d2 and en[2] == d3
			) {
				if (Bd = 2; nfinfo(+2) >= 0 and npface(+2)) {
					edge_e(d1,d2,d3) -= Mag_fpy(*neighbor.data)[Bd];
					e_items(d1,d2,d3)++;
				}
				if (Bd = 1; nfinfo(-3) >= 0 and npface(-3)) {
					edge_e(d1,d2,d3) += Mag_fnz(*neighbor.data)[Bd];
					e_items(d1,d2,d3)++;
				}
			}

			if (d1 = 0, d2 = +1, d3 = -1; cpedge(d1,d2,d3) and
				en[0] == d1 and en[1] == d2 and en[2] == d3
			) {
				if (Bd = 2; nfinfo(-2) >= 0 and npface(-2)) {
					edge_e(d1,d2,d3) -= Mag_fny(*neighbor.data)[Bd];
					e_items(d1,d2,d3)++;
				}
				if (Bd = 1; nfinfo(+3) >= 0 and npface(+3)) {
					edge_e(d1,d2,d3) += Mag_fpz(*neighbor.data)[Bd];
					e_items(d1,d2,d3)++;
				}
			}

			if (d1 = 0, d2 = +1, d3 = +1; cpedge(d1,d2,d3) and
				en[0] == d1 and en[1] == d2 and en[2] == d3
			) {
				if (Bd = 2; nfinfo(-2) >= 0 and npface(-2)) {
					edge_e(d1,d2,d3) -= Mag_fny(*neighbor.data)[Bd];
					e_items(d1,d2,d3)++;
				}
				if (Bd = 1; nfinfo(-3) >= 0 and npface(-3)) {
					edge_e(d1,d2,d3) += Mag_fnz(*neighbor.data)[Bd];
					e_items(d1,d2,d3)++;
				}
			}

			if (d1 = 1, d2 = -1, d3 = -1; cpedge(d1,d2,d3) and
				en[0] == d1 and en[1] == d2 and en[2] == d3
			) {
				if (Bd = 2; nfinfo(+1) >= 0 and npface(+1)) {
					edge_e(d1,d2,d3) += Mag_fpx(*neighbor.data)[Bd];
					e_items(d1,d2,d3)++;
				}
				if (Bd = 0; nfinfo(+3) >= 0 and npface(+3)) {
					edge_e(d1,d2,d3) -= Mag_fpz(*neighbor.data)[Bd];
					e_items(d1,d2,d3)++;
				}
			}

			if (d1 = 1, d2 = -1, d3 = +1; cpedge(d1,d2,d3) and
				en[0] == d1 and en[1] == d2 and en[2] == d3
			) {
				if (Bd = 2; nfinfo(+1) >= 0 and npface(+1)) {
					edge_e(d1,d2,d3) += Mag_fpx(*neighbor.data)[Bd];
					e_items(d1,d2,d3)++;
				}
				if (Bd = 0; nfinfo(-3) >= 0 and npface(-3)) {
					edge_e(d1,d2,d3) -= Mag_fnz(*neighbor.data)[Bd];
					e_items(d1,d2,d3)++;
				}
			}

			if (d1 = 1, d2 = +1, d3 = -1; cpedge(d1,d2,d3) and
				en[0] == d1 and en[1] == d2 and en[2] == d3
			) {
				if (Bd = 2; nfinfo(-1) >= 0 and npface(-1)) {
					edge_e(d1,d2,d3) += Mag_fnx(*neighbor.data)[Bd];
					e_items(d1,d2,d3)++;
				}
				if (Bd = 0; nfinfo(+3) >= 0 and npface(+3)) {
					edge_e(d1,d2,d3) -= Mag_fpz(*neighbor.data)[Bd];
					e_items(d1,d2,d3)++;
				}
			}

			if (d1 = 1, d2 = +1, d3 = +1; cpedge(d1,d2,d3) and
				en[0] == d1 and en[1] == d2 and en[2] == d3
			) {
				if (Bd = 2; nfinfo(-1) >= 0 and npface(-1)) {
					edge_e(d1,d2,d3) += Mag_fnx(*neighbor.data)[Bd];
					e_items(d1,d2,d3)++;
				}
				if (Bd = 0; nfinfo(-3) >= 0 and npface(-3)) {
					edge_e(d1,d2,d3) -= Mag_fnz(*neighbor.data)[Bd];
					e_items(d1,d2,d3)++;
				}
			}

			if (d1 = 2, d2 = -1, d3 = -1; cpedge(d1,d2,d3) and
				en[0] == d1 and en[1] == d2 and en[2] == d3
			) {
				if (Bd = 1; nfinfo(+1) >= 0 and npface(+1)) {
					edge_e(d1,d2,d3) -= Mag_fpx(*neighbor.data)[Bd];
					e_items(d1,d2,d3)++;
				}
				if (Bd = 0; nfinfo(+2) >= 0 and npface(+2)) {
					edge_e(d1,d2,d3) += Mag_fpy(*neighbor.data)[Bd];
					e_items(d1,d2,d3)++;
				}
			}

			if (d1 = 2, d2 = -1, d3 = +1; cpedge(d1,d2,d3) and
				en[0] == d1 and en[1] == d2 and en[2] == d3
			) {
				if (Bd = 1; nfinfo(+1) >= 0 and npface(+1)) {
					edge_e(d1,d2,d3) -= Mag_fpx(*neighbor.data)[Bd];
					e_items(d1,d2,d3)++;
				}
				if (Bd = 0; nfinfo(-2) >= 0 and npface(-2)) {
					edge_e(d1,d2,d3) += Mag_fny(*neighbor.data)[Bd];
					e_items(d1,d2,d3)++;
				}
			}

			if (d1 = 2, d2 = +1, d3 = -1; cpedge(d1,d2,d3) and
				en[0] == d1 and en[1] == d2 and en[2] == d3
			) {
				if (Bd = 1; nfinfo(-1) >= 0 and npface(-1)) {
					edge_e(d1,d2,d3) -= Mag_fnx(*neighbor.data)[Bd];
					e_items(d1,d2,d3)++;
				}
				if (Bd = 0; nfinfo(+2) >= 0 and npface(+2)) {
					edge_e(d1,d2,d3) += Mag_fpy(*neighbor.data)[Bd];
					e_items(d1,d2,d3)++;
				}
			}

			if (d1 = 2, d2 = +1, d3 = +1; cpedge(d1,d2,d3) and
				en[0] == d1 and en[1] == d2 and en[2] == d3
			) {
				if (Bd = 1; nfinfo(-1) >= 0 and npface(-1)) {
					edge_e(d1,d2,d3) -= Mag_fnx(*neighbor.data)[Bd];
					e_items(d1,d2,d3)++;
				}
				if (Bd = 0; nfinfo(-2) >= 0 and npface(-2)) {
					edge_e(d1,d2,d3) += Mag_fny(*neighbor.data)[Bd];
					e_items(d1,d2,d3)++;
				}
			}

			// these shouldn't be possible but kept for reference
			if (fn == -1 and npface(+1)) {
				if (cpedge(1,-1,-1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(1,-1,-1) += get<1>(Mag_f)(*neighbor.data)[2];
					//e_items[1][0][0]++;
				}
				if (cpedge(1,-1,+1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(1,-1,+1) += get<1>(Mag_f)(*neighbor.data)[2];
					//e_items[1][0][1]++;
				}
				if (cpedge(2,-1,-1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(2,-1,-1) -= get<1>(Mag_f)(*neighbor.data)[1];
					//e_items[2][0][0]++;
				}
				if (cpedge(2,-1,+1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(2,-1,+1) -= get<1>(Mag_f)(*neighbor.data)[1];
					//e_items[2][0][1]++;
				}
			}
			if (fn == +1 and npface(-1)) {
				if (cpedge(1,+1,-1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(1,+1,-1) += get<0>(Mag_f)(*neighbor.data)[2];
					//e_items[1][1][0]++;
				}
				if (cpedge(1,+1,+1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(1,+1,+1) += get<0>(Mag_f)(*neighbor.data)[2];
					//e_items[1][1][1]++;
				}
				if (cpedge(2,+1,-1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(2,+1,-1) -= get<0>(Mag_f)(*neighbor.data)[1];
					//e_items[2][1][0]++;
				}
				if (cpedge(2,+1,+1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(2,+1,+1) -= get<0>(Mag_f)(*neighbor.data)[1];
					//e_items[2][1][1]++;
				}
			}
			if (fn == -2 and npface(+2)) {
				if (cpedge(0,-1,-1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(0,-1,-1) -= get<3>(Mag_f)(*neighbor.data)[2];
					//e_items[0][0][0]++;
				}
				if (cpedge(0,-1,+1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(0,-1,+1) -= get<3>(Mag_f)(*neighbor.data)[2];
					//e_items[0][0][1]++;
				}
				if (cpedge(2,-1,-1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(2,-1,-1) += get<3>(Mag_f)(*neighbor.data)[0];
					//e_items[2][0][0]++;
				}
				if (cpedge(2,+1,-1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(2,+1,-1) += get<3>(Mag_f)(*neighbor.data)[0];
					//e_items[2][1][0]++;
				}
			}
			if (fn == +2 and npface(-2)) {
				if (cpedge(0,+1,-1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(0,+1,-1) -= get<2>(Mag_f)(*neighbor.data)[2];
					//e_items[0][1][0]++;
				}
				if (cpedge(0,+1,+1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(0,+1,+1) -= get<2>(Mag_f)(*neighbor.data)[2];
					//e_items[0][1][1]++;
				}
				if (cpedge(2,-1,+1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(2,-1,+1) += get<2>(Mag_f)(*neighbor.data)[0];
					//e_items[2][0][1]++;
				}
				if (cpedge(2,+1,+1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(2,+1,+1) += get<2>(Mag_f)(*neighbor.data)[0];
					//e_items[2][1][1]++;
				}
			}
			if (fn == -3 and npface(+3)) {
				if (cpedge(0,-1,-1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(0,-1,-1) += get<5>(Mag_f)(*neighbor.data)[1];
					//e_items[0][0][0]++;
				}
				if (cpedge(0,+1,-1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(0,+1,-1) += get<5>(Mag_f)(*neighbor.data)[1];
					//e_items[0][1][0]++;
				}
				if (cpedge(1,-1,-1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(1,-1,-1) -= get<5>(Mag_f)(*neighbor.data)[0];
					//e_items[1][0][0]++;
				}
				if (cpedge(1,+1,-1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(1,+1,-1) -= get<5>(Mag_f)(*neighbor.data)[0];
					//e_items[1][1][0]++;
				}
			}
			if (fn == +3 and npface(-3)) {
				if (cpedge(0,-1,+1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(0,-1,+1) += get<4>(Mag_f)(*neighbor.data)[1];
					//e_items[0][0][1]++;
				}
				if (cpedge(0,+1,+1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(0,+1,+1) += get<4>(Mag_f)(*neighbor.data)[1];
					//e_items[0][1][1]++;
				}
				if (cpedge(1,-1,+1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(1,-1,+1) -= get<4>(Mag_f)(*neighbor.data)[0];
					//e_items[1][0][1]++;
				}
				if (cpedge(1,+1,+1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(1,+1,+1) -= get<4>(Mag_f)(*neighbor.data)[0];
					//e_items[1][1][1]++;
				}
			}
		}

		for (int d1: {0, 1, 2})
		for (int d2: {-1, +1})
		for (int d3: {-1, +1}) {
			if (not cpedge(d1,d2,d3)) {
				if (e_items(d1,d2,d3) != 0) throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
				continue;
			}
			if (e_items(d1,d2,d3) > 0) edge_e(d1,d2,d3) /= e_items(d1,d2,d3);
		}

		const auto [rx, ry, rz] = grid.geometry.get_center(cell.id);
		if (Mas(*cell.data) < 0) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << ") "
				<< "Negative mass density in cell " << cell.id
				<< " at " << rx << ", " << ry << ", " << rz << std::endl;
			abort();
		}

		if (Nrj(*cell.data) < 0) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << ") "
				<< "Negative total energy density in cell " << cell.id
				<< " at " << rx << ", " << ry << ", " << rz << std::endl;
			abort();
		}
	}
} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


/*! Solves new face magnetic fields from edge electric fields.

Equations 13-15 of https://doi.org/10.1006/jcph.1998.6153

Saves new face B of cells with SInfo(*cell_data) == 1,
ignores cells with SInfo < 0.
*/
template <
	class Cells,
	class Grid,
	class Face_Magnetic_Field_Getter,
	class Edge_Electric_Field_Getter,
	class Primary_Face_Getter,
	class Primary_Edge_Getter,
	class Face_Info_Getter,
	class Substepping_Period_Getter
> void get_face_magnetic_field(
	const Cells& cells,
	Grid& grid,
	const double sub_dt,
	const int current_substep,
	const Face_Magnetic_Field_Getter Face_B,
	const Edge_Electric_Field_Getter Edge_E,
	const Primary_Face_Getter PFace,
	const Primary_Edge_Getter PEdge,
	const Face_Info_Getter FInfo,
	const Substepping_Period_Getter Substep
) try {
	using std::runtime_error;
	using std::to_string;

	for (const auto& cell: cells) {
		if (current_substep % Substep(*cell.data) != 0) {
			continue;
		}

		int max_substep = Substep(*cell.data);
		for (const auto& neighbor: cell.neighbors_of) {
			const auto& fn = neighbor.face_neighbor;
			if (fn == 0 and neighbor.edge_neighbor[0] < 0) {
				continue;
			}
			if (fn != 0 and PFace(*neighbor.data)(-fn) and FInfo(*neighbor.data)(-fn) < 0) {
				continue;
			}
			max_substep = std::max(Substep(*neighbor.data), max_substep);
		}

		const auto cell_length = grid.geometry.get_length(cell.id);

		typename std::remove_reference<decltype(Face_B(*cell.data))>::type face_db{0, 0, 0, 0, 0, 0};

		const auto& cpface = PFace(*cell.data);
		const auto& cfinfo = FInfo(*cell.data);
		const auto& cpedge = PEdge(*cell.data);
		const auto& cedge_e = Edge_E(*cell.data);

		/*! Contributions from edge electric fields to face Bs

		Right-hand line integral over edges touching each face:
		dBx: E(2, *,+1)dz - E(2, *,-1)dz - E(1, *,+1)dx + E(1, *,-1)dx
		dBy: E(0, *,+1)dx - E(0, *,-1)dx - E(2,+1, *)dz + E(2,-1, *)dz
		dBz: E(1,+1, *)dy - E(1,-1, *)dy + E(0,-1, *)dx - E(0,+1, *)dx
		*/
		if (cfinfo(0,-1) >= 0 and cpface(0,-1)) {
			if (cpedge(1,-1,-1)) face_db(0, -1) += cell_length[1]*cedge_e(1,-1,-1);
			if (cpedge(1,-1,+1)) face_db(0, -1) -= cell_length[1]*cedge_e(1,-1,+1);
			if (cpedge(2,-1,-1)) face_db(0, -1) -= cell_length[2]*cedge_e(2,-1,-1);
			if (cpedge(2,-1,+1)) face_db(0, -1) += cell_length[2]*cedge_e(2,-1,+1);
		}
		if (cfinfo(0,+1) >= 0 and cpface(0,+1)) {
			if (cpedge(1,+1,-1)) face_db(0, +1) += cell_length[1]*cedge_e(1,+1,-1);
			if (cpedge(1,+1,+1)) face_db(0, +1) -= cell_length[1]*cedge_e(1,+1,+1);
			if (cpedge(2,+1,-1)) face_db(0, +1) -= cell_length[2]*cedge_e(2,+1,-1);
			if (cpedge(2,+1,+1)) face_db(0, +1) += cell_length[2]*cedge_e(2,+1,+1);
		}
		if (cfinfo(1,-1) >= 0 and cpface(1,-1)) {
			if (cpedge(0,-1,-1)) face_db(1, -1) -= cell_length[0]*cedge_e(0,-1,-1);
			if (cpedge(0,-1,+1)) face_db(1, -1) += cell_length[0]*cedge_e(0,-1,+1);
			if (cpedge(2,-1,-1)) face_db(1, -1) += cell_length[2]*cedge_e(2,-1,-1);
			if (cpedge(2,+1,-1)) face_db(1, -1) -= cell_length[2]*cedge_e(2,+1,-1);
		}
		if (cfinfo(1,+1) >= 0 and cpface(1,+1)) {
			if (cpedge(0,+1,-1)) face_db(1, +1) -= cell_length[0]*cedge_e(0,+1,-1);
			if (cpedge(0,+1,+1)) face_db(1, +1) += cell_length[0]*cedge_e(0,+1,+1);
			if (cpedge(2,-1,+1)) face_db(1, +1) += cell_length[2]*cedge_e(2,-1,+1);
			if (cpedge(2,+1,+1)) face_db(1, +1) -= cell_length[2]*cedge_e(2,+1,+1);
		}
		if (cfinfo(2,-1) >= 0 and cpface(2,-1)) {
			if (cpedge(0,-1,-1)) face_db(2, -1) += cell_length[0]*cedge_e(0,-1,-1);
			if (cpedge(0,+1,-1)) face_db(2, -1) -= cell_length[0]*cedge_e(0,+1,-1);
			if (cpedge(1,-1,-1)) face_db(2, -1) -= cell_length[1]*cedge_e(1,-1,-1);
			if (cpedge(1,+1,-1)) face_db(2, -1) += cell_length[1]*cedge_e(1,+1,-1);
		}
		if (cfinfo(2,+1) >= 0 and cpface(2,+1)) {
			if (cpedge(0,-1,+1)) face_db(2, +1) += cell_length[0]*cedge_e(0,-1,+1);
			if (cpedge(0,+1,+1)) face_db(2, +1) -= cell_length[0]*cedge_e(0,+1,+1);
			if (cpedge(1,-1,+1)) face_db(2, +1) -= cell_length[1]*cedge_e(1,-1,+1);
			if (cpedge(1,+1,+1)) face_db(2, +1) += cell_length[1]*cedge_e(1,+1,+1);
		}

		const auto cleni = grid.mapping.get_cell_length_in_indices(cell.id);

		for (const auto& neighbor: cell.neighbors_of) {
			const auto& fn = neighbor.face_neighbor;
			const auto& en = neighbor.edge_neighbor;

			if (fn == 0 and en[0] < 0) {
				continue;
			}

			if (current_substep % Substep(*neighbor.data) != 0) {
				continue;
			}

			const auto neigh_length = grid.geometry.get_length(neighbor.id);
			const auto nleni = grid.mapping.get_cell_length_in_indices(neighbor.id);

			const auto& npedge = PEdge(*neighbor.data);
			const auto& nedge_e = Edge_E(*neighbor.data);

			if (cfinfo(-1) >= 0 and cpface(-1) and neighbor.x == 0) {
				if (fn == -2 and npedge(2,-1,+1)) face_db(0,-1) -= neigh_length[2]*nedge_e(2,-1,+1);
				if (fn == +2 and npedge(2,-1,-1)) face_db(0,-1) += neigh_length[2]*nedge_e(2,-1,-1);
				if (fn == -3 and npedge(1,-1,+1)) face_db(0,-1) += neigh_length[1]*nedge_e(1,-1,+1);
				if (fn == +3 and npedge(1,-1,-1)) face_db(0,-1) -= neigh_length[1]*nedge_e(1,-1,-1);
			}
			if (cfinfo(+1) >= 0 and cpface(+1) and cleni == neighbor.x + nleni) {
				if (fn == -2 and npedge(2,+1,+1)) face_db(0,+1) -= neigh_length[2]*nedge_e(2,+1,+1);
				if (fn == +2 and npedge(2,+1,-1)) face_db(0,+1) += neigh_length[2]*nedge_e(2,+1,-1);
				if (fn == -3 and npedge(1,+1,+1)) face_db(0,+1) += neigh_length[1]*nedge_e(1,+1,+1);
				if (fn == +3 and npedge(1,+1,-1)) face_db(0,+1) -= neigh_length[1]*nedge_e(1,+1,-1);
			}
			if (cfinfo(-2) >= 0 and cpface(-2) and neighbor.y == 0) {
				if (fn == -1 and npedge(2,+1,-1)) face_db(1,-1) += neigh_length[2]*nedge_e(2,+1,-1);
				if (fn == +1 and npedge(2,-1,-1)) face_db(1,-1) -= neigh_length[2]*nedge_e(2,-1,-1);
				if (fn == -3 and npedge(0,-1,+1)) face_db(1,-1) -= neigh_length[0]*nedge_e(0,-1,+1);
				if (fn == +3 and npedge(0,-1,-1)) face_db(1,-1) += neigh_length[0]*nedge_e(0,-1,-1);
			}
			if (cfinfo(+2) >= 0 and cpface(+2) and cleni == neighbor.y + nleni) {
				if (fn == -1 and npedge(2,+1,+1)) face_db(1,+1) += neigh_length[2]*nedge_e(2,+1,+1);
				if (fn == +1 and npedge(2,-1,+1)) face_db(1,+1) -= neigh_length[2]*nedge_e(2,-1,+1);
				if (fn == -3 and npedge(0,+1,+1)) face_db(1,+1) -= neigh_length[0]*nedge_e(0,+1,+1);
				if (fn == +3 and npedge(0,+1,-1)) face_db(1,+1) += neigh_length[0]*nedge_e(0,+1,-1);
			}
			if (cfinfo(-3) >= 0 and cpface(-3) and neighbor.z == 0) {
				if (fn == -1 and npedge(1,+1,-1)) face_db(2,-1) -= neigh_length[1]*nedge_e(1,+1,-1);
				if (fn == +1 and npedge(1,-1,-1)) face_db(2,-1) += neigh_length[1]*nedge_e(1,-1,-1);
				if (fn == -2 and npedge(0,+1,-1)) face_db(2,-1) += neigh_length[0]*nedge_e(0,+1,-1);
				if (fn == +2 and npedge(0,-1,-1)) face_db(2,-1) -= neigh_length[0]*nedge_e(0,-1,-1);
			}
			if (cfinfo(+3) >= 0 and cpface(+3) and cleni == neighbor.z + nleni) {
				if (fn == -1 and npedge(1,+1,+1)) face_db(2,+1) -= neigh_length[1]*nedge_e(1,+1,+1);
				if (fn == +1 and npedge(1,-1,+1)) face_db(2,+1) += neigh_length[1]*nedge_e(1,-1,+1);
				if (fn == -2 and npedge(0,+1,+1)) face_db(2,+1) += neigh_length[0]*nedge_e(0,+1,+1);
				if (fn == +2 and npedge(0,-1,+1)) face_db(2,+1) -= neigh_length[0]*nedge_e(0,-1,+1);
			}
			if (en[0] == 0 and en[1] == -1 and en[2] == -1 and npedge(0,+1,+1)) {
				if (cpedge(0,-1,-1)) throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
				if (cfinfo(-2) >= 0 and cpface(-2)) face_db(1,-1) -= neigh_length[0]*nedge_e(0,+1,+1);
				if (cfinfo(-3) >= 0 and cpface(-3)) face_db(2,-1) += neigh_length[0]*nedge_e(0,+1,+1);
			}
			if (en[0] == 0 and en[1] == -1 and en[2] == +1 and npedge(0,+1,-1)) {
				if (cpedge(0,-1,+1)) throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
				if (cfinfo(-2) >= 0 and cpface(-2)) face_db(1,-1) += neigh_length[0]*nedge_e(0,+1,-1);
				if (cfinfo(+3) >= 0 and cpface(+3)) face_db(2,+1) += neigh_length[0]*nedge_e(0,+1,-1);
			}
			if (en[0] == 0 and en[1] == +1 and en[2] == -1 and npedge(0,-1,+1)) {
				if (cpedge(0,+1,-1)) throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
				if (cfinfo(+2) >= 0 and cpface(+2)) face_db(1,+1) -= neigh_length[0]*nedge_e(0,-1,+1);
				if (cfinfo(-3) >= 0 and cpface(-3)) face_db(2,-1) -= neigh_length[0]*nedge_e(0,-1,+1);
			}
			if (en[0] == 0 and en[1] == +1 and en[2] == +1 and npedge(0,-1,-1)) {
				if (cpedge(0,+1,+1)) throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
				if (cfinfo(+2) >= 0 and cpface(+2)) face_db(1,+1) += neigh_length[0]*nedge_e(0,-1,-1);
				if (cfinfo(+3) >= 0 and cpface(+3)) face_db(2,+1) -= neigh_length[0]*nedge_e(0,-1,-1);
			}
			if (en[0] == 1 and en[1] == -1 and en[2] == -1 and npedge(1,+1,+1)) {
				if (cpedge(1,-1,-1)) throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
				if (cfinfo(-1) >= 0 and cpface(-1)) face_db(0,-1) += neigh_length[1]*nedge_e(1,+1,+1);
				if (cfinfo(-3) >= 0 and cpface(-3)) face_db(2,-1) -= neigh_length[1]*nedge_e(1,+1,+1);
			}
			if (en[0] == 1 and en[1] == -1 and en[2] == +1 and npedge(1,+1,-1)) {
				if (cpedge(1,-1,+1)) throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
				if (cfinfo(-1) >= 0 and cpface(-1)) face_db(0,-1) -= neigh_length[1]*nedge_e(1,+1,-1);
				if (cfinfo(+3) >= 0 and cpface(+3)) face_db(2,+1) -= neigh_length[1]*nedge_e(1,+1,-1);
			}
			if (en[0] == 1 and en[1] == +1 and en[2] == -1 and npedge(1,-1,+1)) {
				if (cpedge(1,+1,-1)) throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
				if (cfinfo(+1) >= 0 and cpface(+1)) face_db(0,+1) += neigh_length[1]*nedge_e(1,-1,+1);
				if (cfinfo(-3) >= 0 and cpface(-3)) face_db(2,-1) += neigh_length[1]*nedge_e(1,-1,+1);
			}
			if (en[0] == 1 and en[1] == +1 and en[2] == +1 and npedge(1,-1,-1)) {
				if (cpedge(1,+1,+1)) throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
				if (cfinfo(+1) >= 0 and cpface(+1)) face_db(0,+1) -= neigh_length[1]*nedge_e(1,-1,-1);
				if (cfinfo(+3) >= 0 and cpface(+3)) face_db(2,+1) += neigh_length[1]*nedge_e(1,-1,-1);
			}
			if (en[0] == 2 and en[1] == -1 and en[2] == -1 and npedge(2,+1,+1)) {
				if (cpedge(2,-1,-1)) throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
				if (cfinfo(-1) >= 0 and cpface(-1)) face_db(0,-1) -= neigh_length[2]*nedge_e(2,+1,+1);
				if (cfinfo(-2) >= 0 and cpface(-2)) face_db(1,-1) += neigh_length[2]*nedge_e(2,+1,+1);
			}
			if (en[0] == 2 and en[1] == -1 and en[2] == +1 and npedge(2,+1,-1)) {
				if (cpedge(2,-1,+1)) throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
				if (cfinfo(-1) >= 0 and cpface(-1)) face_db(0,-1) += neigh_length[2]*nedge_e(2,+1,-1);
				if (cfinfo(+2) >= 0 and cpface(+2)) face_db(1,+1) += neigh_length[2]*nedge_e(2,+1,-1);
			}
			if (en[0] == 2 and en[1] == +1 and en[2] == -1 and npedge(2,-1,+1)) {
				if (cpedge(2,+1,-1)) throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
				if (cfinfo(+1) >= 0 and cpface(+1)) face_db(0,+1) -= neigh_length[2]*nedge_e(2,-1,+1);
				if (cfinfo(-2) >= 0 and cpface(-2)) face_db(1,-1) -= neigh_length[2]*nedge_e(2,-1,+1);
			}
			if (en[0] == 2 and en[1] == +1 and en[2] == +1 and npedge(2,-1,-1)) {
				if (cpedge(2,+1,+1)) throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
				if (cfinfo(+1) >= 0 and cpface(+1)) face_db(0,+1) += neigh_length[2]*nedge_e(2,-1,-1);
				if (cfinfo(+2) >= 0 and cpface(+2)) face_db(1,+1) -= neigh_length[2]*nedge_e(2,-1,-1);
			}
			// these shouldn't be possible but kept for reference
			if (cpface(-1) and fn == -1) {
				if (npedge(1,+1,-1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db(0, -1) += neigh_length[1]*nedge_e(1,+1,-1);
				}
				if (npedge(1,+1,+1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db(0, -1) -= neigh_length[1]*nedge_e(1,+1,+1);
				}
				if (npedge(2,+1,-1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db(0, -1) -= neigh_length[2]*nedge_e(2,+1,-1);
				}
				if (npedge(2,+1,+1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db(0, -1) += neigh_length[2]*nedge_e(2,+1,+1);
				}
			}
			if (cpface(+1) and fn == +1) {
				if (npedge(1,-1,-1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db(0, +1) += neigh_length[1]*nedge_e(1,-1,-1);
				}
				if (npedge(1,-1,+1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db(0, +1) -= neigh_length[1]*nedge_e(1,-1,+1);
				}
				if (npedge(2,-1,-1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db(0, +1) -= neigh_length[2]*nedge_e(2,-1,-1);
				}
				if (npedge(2,-1,+1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db(0, +1) += neigh_length[2]*nedge_e(2,-1,+1);
				}
			}
			if (cpface(-2) and fn == -2) {
				if (npedge(0,+1,-1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db(1, -1) -= neigh_length[0]*nedge_e(0,+1,-1);
				}
				if (npedge(0,+1,+1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db(1, -1) += neigh_length[0]*nedge_e(0,+1,+1);
				}
				if (npedge(2,-1,+1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db(1, -1) += neigh_length[2]*nedge_e(2,-1,+1);
				}
				if (npedge(2,+1,+1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db(1, -1) -= neigh_length[2]*nedge_e(2,+1,+1);
				}
			}
			if (cpface(+2) and fn == +2) {
				if (npedge(0,-1,-1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db(1, +1) -= neigh_length[0]*nedge_e(0,-1,-1);
				}
				if (npedge(0,-1,+1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db(1, +1) += neigh_length[0]*nedge_e(0,-1,+1);
				}
				if (npedge(2,-1,-1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db(1, +1) += neigh_length[2]*nedge_e(2,-1,-1);
				}
				if (npedge(2,+1,-1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db(1, +1) -= neigh_length[2]*nedge_e(2,+1,-1);
				}
			}
			if (cpface(-3) and fn == -3) {
				if (npedge(0,-1,+1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db(2, -1) += neigh_length[0]*nedge_e(0,-1,+1);
				}
				if (npedge(0,+1,+1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db(2, -1) -= neigh_length[0]*nedge_e(0,+1,+1);
				}
				if (npedge(1,-1,+1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db(2, -1) -= neigh_length[1]*nedge_e(1,-1,+1);
				}
				if (npedge(1,+1,+1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db(2, -1) += neigh_length[1]*nedge_e(1,+1,+1);
				}
			}
			if (cpface(+3) and fn == +3) {
				if (npedge(0,-1,-1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db(2, +1) += neigh_length[0]*nedge_e(0,-1,-1);
				}
				if (npedge(0,+1,-1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db(2, +1) -= neigh_length[0]*nedge_e(0,+1,-1);
				}
				if (npedge(1,-1,-1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db(2, +1) -= neigh_length[1]*nedge_e(1,-1,-1);
				}
				if (npedge(1,+1,-1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db(2, +1) += neigh_length[1]*nedge_e(1,+1,-1);
				}
			}
		}
		const std::array<double, 3> area{
			cell_length[1]*cell_length[2],
			cell_length[0]*cell_length[2],
			cell_length[0]*cell_length[1]
		};
		const auto dt = sub_dt * Substep(*cell.data);
		for (auto dir: {-3,-2,-1,+1,+2,+3}) {
			if (cfinfo(dir) >= 0 and cpface(dir)) {
				Face_B(*cell.data)(dir) -= dt*face_db(dir) / area[abs(dir) - 1];
			}
		}
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
	class Primary_Face_Getter,
	class Face_Type_Getter,
	class Substepping_Period_Getter
> void update_B_consistency(
	const Cell_Iter& cells,
	const Mass_Density_Getter Mas,
	const Momentum_Density_Getter Mom,
	const Total_Energy_Density_Getter Nrj,
	const Volume_Magnetic_Field_Getter VMag,
	const Face_Magnetic_Field_Getter Face_B,
	const Primary_Face_Getter PFace,
	const Face_Type_Getter FInfo,
	const Substepping_Period_Getter Substep,
	const double adiabatic_index,
	const double vacuum_permeability,
	const bool constant_thermal_pressure
) try {
	using std::runtime_error;
	using std::to_string;

	for (const auto& cell: cells) {
		int max_substep = Substep(*cell.data);
		for (const auto& neighbor: cell.neighbors_of) {
			const auto& fn = neighbor.face_neighbor;
			if (fn == 0 and neighbor.edge_neighbor[0] < 0) {
				continue;
			}
			if (fn != 0 and PFace(*neighbor.data)(-fn) and FInfo(*neighbor.data)(-fn) < 0) {
				continue;
			}
			max_substep = std::max(Substep(*neighbor.data), max_substep);
		}

		if (constant_thermal_pressure and Mas(*cell.data) <= 0) {
			continue;
		}

		const auto old_pressure = [&](){
			if (constant_thermal_pressure) {
				return pamhd::mhd::get_pressure(
					Mas(*cell.data), Mom(*cell.data), Nrj(*cell.data), VMag(*cell.data),
					adiabatic_index, vacuum_permeability);
			} else {
				return 0.0;
			}
		}();

		auto& c_face_b = Face_B(*cell.data);
		const auto& cpface = PFace(*cell.data);
		const auto& cfinfo = FInfo(*cell.data);

		for (auto dim: {0, 1, 2})
		for (auto side: {-1, +1}) {
			if (not cpface(dim, side)) {
				c_face_b(dim, side) = 0;
			} else if (cfinfo(dim, side) < 0) { // handle dont_solve face
				c_face_b(dim, side) = 0;
				if (cpface(dim, -side) and cfinfo(dim, -side) >= 0) {
					c_face_b(dim, side) = c_face_b(dim, -side);
				}
			}
		}

		for (const auto& neighbor: cell.neighbors_of) {
			const auto& fn = neighbor.face_neighbor;
			if (fn == 0) continue;

			const double factor = [&]{
				if (neighbor.relative_size > 0) {
					return 0.25;
				} else {
					return 1.0;
				}
			}();

			const auto& n_face_b = Face_B(*neighbor.data);
			const auto& npface = PFace(*neighbor.data);
			const auto& nfinfo = FInfo(*neighbor.data);

			if (not cpface(fn)) {
				if (not npface(-fn)) throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
				if (nfinfo(-fn) < 0) continue;

				c_face_b(fn) += n_face_b(-fn) * factor;
			}
			// handle dont_solve face
			if (cpface(-fn) and cfinfo(-fn) < 0 and npface(-fn) and nfinfo(-fn) >= 0) {
				c_face_b(-fn) += n_face_b(-fn) * factor;
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


//! Returns length of substep and max substep interval among all cells
template <
	class Grid,
	class Primary_Face_Getter,
	class Solver_Info_Getter,
	class Substepping_Period_Getter,
	class Max_Velocity_Getter
> std::tuple<double, int> update_substeps(
	Grid& grid,
	const double max_dt,
	const double dt_factor,
	const Primary_Face_Getter PFace,
	const Solver_Info_Getter SInfo,
	const Substepping_Period_Getter Substep,
	const Max_Velocity_Getter Max_v
) try {
	using std::make_tuple;
	using std::max;
	using std::min;

	if (max_dt <= 0) {
		return make_tuple(0.0, 1);
	}

	double
		min_dt_local = std::numeric_limits<double>::max(),
		min_dt_global = std::numeric_limits<double>::max();
	for (const auto& cell: grid.local_cells()) {
		if (SInfo(*cell.data) < 0) {
			continue;
		}

		const auto clen = grid.geometry.get_length(cell.id);
		const auto& cpface = PFace(*cell.data);
		for (const auto& neighbor: cell.neighbors_of) {
			const auto& fn = neighbor.face_neighbor;
			if (fn == 0 or SInfo(*neighbor.data) < 0) {
				continue;
			}

			const auto neigh_dim = size_t(abs(fn) - 1);
			if (cpface(fn)) {
				min_dt_local = min(clen[neigh_dim] / Max_v(*cell.data)(fn), min_dt_local);
			} else {
				min_dt_local = min(clen[neigh_dim] / Max_v(*neighbor.data)(-fn), min_dt_local);
				const auto& npface = PFace(*neighbor.data);
				if (not npface(-fn)) {
					throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
				}
			}
		}
	}
	auto comm = grid.get_communicator();
	if (
		MPI_Allreduce(
			&min_dt_local,
			&min_dt_global,
			1,
			MPI_DOUBLE,
			MPI_MIN,
			comm
		) != MPI_SUCCESS
	) {
		std::cerr << __FILE__ "(" << __LINE__
			<< "): Couldn't reduce min_dt_local."
			<< std::endl;
		abort();
	}
	min_dt_global *= dt_factor;
	if (min_dt_global > max_dt) {
		min_dt_global = max_dt;
	}

	for (const auto& cell: grid.local_cells()) {
		if (SInfo(*cell.data) < 0) {
			continue;
		}

		double min_dt_cell = std::numeric_limits<double>::max();
		const auto clen = grid.geometry.get_length(cell.id);
		const auto& cpface = PFace(*cell.data);

		for (const auto& neighbor: cell.neighbors_of) {
			const auto& fn = neighbor.face_neighbor;
			if (fn == 0 or SInfo(*neighbor.data) < 0) {
				continue;
			}

			const auto neigh_dim = size_t(abs(fn) - 1);
			if (cpface(fn)) {
				min_dt_cell = min(clen[neigh_dim] / Max_v(*cell.data)(fn), min_dt_cell);
			} else {
				min_dt_cell = min(clen[neigh_dim] / Max_v(*neighbor.data)(-fn), min_dt_cell);
			}
		}
		min_dt_cell *= dt_factor;
		// 2^N form
		if (min_dt_cell > 0) {
			Substep(*cell.data) = max(0, (int)std::floor(std::log2(min_dt_cell / min_dt_global)));
		} else {
			Substep(*cell.data) = 0;
		}
	}

	// reduce substep difference between neighbors to 2x
	Grid::cell_data_type::set_transfer_all(true, pamhd::mhd::Substepping_Period());
	uint64_t modified_cells = 0;
	int max_substep = 0;
	do {
		grid.update_copies_of_remote_neighbors();
		uint64_t modified_cells_local = 0;
		int max_substep_local = 0;
		for (const auto& cell: grid.local_cells()) {
			if (SInfo(*cell.data) < 0) {
				continue;
			}
			for (const auto& neighbor: cell.neighbors_of) {
				if (neighbor.face_neighbor == 0 and neighbor.edge_neighbor[0] < 0) {
					continue;
				}
				if (SInfo(*neighbor.data) < 0) {
					continue;
				}

				if (Substep(*cell.data) > Substep(*neighbor.data) + 1) {
					Substep(*cell.data) = Substep(*neighbor.data) + 1;
					modified_cells_local++;
				}
				max_substep_local = max(Substep(*cell.data), max_substep_local);
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
		if (modified_cells == 0) {
			if (
				MPI_Allreduce(
					&max_substep_local,
					&max_substep,
					1,
					MPI_INT,
					MPI_MAX,
					comm
				) != MPI_SUCCESS
			) {
				std::cerr << __FILE__ "(" << __LINE__
					<< "): Couldn't reduce modified cell count."
					<< std::endl;
				abort();
			}
		}
	} while (modified_cells > 0);
	Grid::cell_data_type::set_transfer_all(false, pamhd::mhd::Substepping_Period());
	MPI_Comm_free(&comm);

	// convert substep from 2^N to N
	max_substep = 1 << max_substep;
	for (const auto& cell: grid.local_cells()) {
		if (SInfo(*cell.data) < 0) {
			continue;
		}
		Substep(*cell.data) = 1 << Substep(*cell.data);
	}

	return make_tuple(min_dt_global, max_substep);

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
	class Edge_Electric_Field_Getter,
	class Background_Magnetic_Field_Getter,
	class Mass_Density_Flux_Getters,
	class Momentum_Density_Flux_Getters,
	class Total_Energy_Density_Flux_Getters,
	class Magnetic_Field_Flux_Getters,
	class Primary_Face_Getter,
	class Primary_Edge_Getter,
	class Solver_Info_Getter,
	class Face_Info_Getter,
	class Substepping_Period_Getter,
	class Max_Velocity_Getter
> double timestep(
	const Solver solver,
	Grid& grid,
	const double max_time_step,
	const double time_step_factor,
	const double adiabatic_index,
	const double vacuum_permeability,
	const Mass_Density_Getter Mas,
	const Momentum_Density_Getter Mom,
	const Total_Energy_Density_Getter Nrj,
	const Magnetic_Field_Getter Mag,
	const Face_Magnetic_Field_Getter Face_B,
	const Edge_Electric_Field_Getter Edge_E,
	const Background_Magnetic_Field_Getter Bg_B,
	const Mass_Density_Flux_Getters Mas_fs,
	const Momentum_Density_Flux_Getters Mom_fs,
	const Total_Energy_Density_Flux_Getters Nrj_fs,
	const Magnetic_Field_Flux_Getters Mag_fs,
	const Primary_Face_Getter PFace,
	const Primary_Edge_Getter PEdge,
	const Solver_Info_Getter SInfo,
	const Face_Info_Getter FInfo,
	const Substepping_Period_Getter Substep,
	const Max_Velocity_Getter Max_v
) try {
	using std::max;
	using std::min;

	const auto [sub_dt, max_substep] = update_substeps(
		grid, max_time_step, time_step_factor,
		PFace, SInfo, Substep, Max_v);

std::cout << "\n" __FILE__ "(" << __LINE__ << ") Entering substep loop: " << sub_dt << " " << max_substep << std::endl;
	double total_dt = 0;
	for (int substep = 1; substep <= max_substep; substep++) {
		total_dt += sub_dt;

	Grid::cell_data_type::set_transfer_all(true, pamhd::mhd::MHD_State_Conservative());
	grid.start_remote_neighbor_copy_updates();

	pamhd::mhd::get_fluxes(
		solver, grid.inner_cells(), grid, substep,
		adiabatic_index, vacuum_permeability,
		Mas, Mom, Nrj, Mag, Bg_B,
		Mas_fs, Mom_fs, Nrj_fs, Mag_fs,
		PFace, SInfo, Substep, Max_v
	);

	grid.wait_remote_neighbor_copy_update_receives();

	pamhd::mhd::get_fluxes(
		solver, grid.outer_cells(), grid, substep,
		adiabatic_index, vacuum_permeability,
		Mas, Mom, Nrj, Mag, Bg_B,
		Mas_fs, Mom_fs, Nrj_fs, Mag_fs,
		PFace, SInfo, Substep, Max_v
	);

	grid.wait_remote_neighbor_copy_update_sends();
	Grid::cell_data_type::set_transfer_all(false, pamhd::mhd::MHD_State_Conservative());

	Grid::cell_data_type::set_transfer_all(true,
		pamhd::mhd::MHD_Flux(),
		pamhd::mhd::Max_Velocity());
	grid.update_copies_of_remote_neighbors();
	Grid::cell_data_type::set_transfer_all(false,
		pamhd::mhd::MHD_Flux(),
		pamhd::mhd::Max_Velocity());

	// TODO: split into inner and outer cells
	pamhd::mhd::get_edge_electric_field(
		grid, sub_dt, substep,
		Mas, Mom, Nrj, Mag, Edge_E,
		Mas_fs, Mom_fs, Nrj_fs, Mag_fs,
		PFace, PEdge, FInfo, Substep
	);

	Grid::cell_data_type::set_transfer_all(true,
		pamhd::Edge_Electric_Field(),
		// update pressure for B consistency calculation
		pamhd::mhd::MHD_State_Conservative()
	);
	grid.update_copies_of_remote_neighbors();
	Grid::cell_data_type::set_transfer_all(false,
		pamhd::Edge_Electric_Field(),
		pamhd::mhd::MHD_State_Conservative()
	);

	pamhd::mhd::get_face_magnetic_field(
		grid.local_cells(), grid,
		sub_dt, substep, Face_B, Edge_E,
		PFace, PEdge, FInfo, Substep
	);
	Grid::cell_data_type::set_transfer_all(true, pamhd::Face_Magnetic_Field());
	grid.update_copies_of_remote_neighbors();
	Grid::cell_data_type::set_transfer_all(false, pamhd::Face_Magnetic_Field());

	// constant thermal pressure when updating vol B after solution
	pamhd::mhd::update_B_consistency(
		grid.local_cells(),
		Mas, Mom, Nrj, Mag, Face_B,
		PFace, FInfo, Substep,
		adiabatic_index,
		vacuum_permeability,
		true
	);
	Grid::cell_data_type::set_transfer_all(true,
		pamhd::Face_Magnetic_Field(),
		pamhd::mhd::MHD_State_Conservative()
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


}} // namespaces


#endif // ifndef PAMHD_MHD_SOLVE_STAGGERED_HPP
