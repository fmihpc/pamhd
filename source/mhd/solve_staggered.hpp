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
	class Solver_Info_Getter
> double get_fluxes(
	const Solver solver,
	const Cell_Iter& cells,
	Grid& grid,
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
	const Solver_Info_Getter SInfo
) try {
	using std::abs;
	using std::get;
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

	// maximum allowed next time step for cells of this process
	double max_dt = std::numeric_limits<double>::max();

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
		const auto cell_length = grid.geometry.get_length(cell.id);

		const auto& cpface = PFace(*cell.data);

		for (const auto& neighbor: cell.neighbors_of) {
			const auto n = neighbor.face_neighbor;
			// flux stored in neighbor and calculated by its owner
			if (n == 0 or not cpface(n)) continue;

			if (SInfo(*neighbor.data) < 0) {
				continue;
			}

			detail::MHD state_neg, state_pos;
			// take into account direction of neighbor from cell
			state_neg[mas_int] = Mas(*cell.data);
			state_neg[mom_int] = get_rotated_vector(Mom(*cell.data), abs(neighbor.face_neighbor));
			state_neg[nrj_int] = Nrj(*cell.data);
			state_neg[mag_int] = get_rotated_vector(Mag(*cell.data), abs(neighbor.face_neighbor));

			state_pos[mas_int] = Mas(*neighbor.data);
			state_pos[mom_int] = get_rotated_vector(Mom(*neighbor.data), abs(neighbor.face_neighbor));
			state_pos[nrj_int] = Nrj(*neighbor.data);
			state_pos[mag_int] = get_rotated_vector(Mag(*neighbor.data), abs(neighbor.face_neighbor));
			if (neighbor.face_neighbor < 0) {
				const auto temp = state_neg;
				state_neg = state_pos;
				state_pos = temp;
			}

			const size_t neighbor_dim = size_t(abs(neighbor.face_neighbor) - 1);
			Magnetic_Field::data_type bg_face_b;
			if (neighbor.face_neighbor < 0) {
				bg_face_b = Bg_B(*neighbor.data)(neighbor_dim, 1);
			} else {
				bg_face_b = Bg_B(*cell.data)(neighbor_dim, 1);
			}
			bg_face_b = get_rotated_vector(bg_face_b, abs(neighbor.face_neighbor));

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
					<< " of cell type " << SInfo(*cell.data)
					<< " and " << SInfo(*neighbor.data)
					<< " at " << grid.geometry.get_center(cell.id)
					<< " and " << grid.geometry.get_center(neighbor.id)
					<< " in direction " << neighbor.face_neighbor
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
			flux[mom_int] = get_rotated_vector(flux[mom_int], -abs(neighbor.face_neighbor));
			flux[mag_int] = get_rotated_vector(flux[mag_int], -abs(neighbor.face_neighbor));

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

	return max_dt;

} catch (const std::exception& e) {
	throw std::runtime_error(__func__ + std::string(": ") + e.what());
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
	class Face_Info_Getter
> void get_edge_electric_field(
	Grid& grid,
	const double dt,
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
	const Face_Info_Getter FInfo
) try {
	using std::array;
	using std::get;
	using std::make_tuple;
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
		auto& edge_e = Edge_E(*cell.data);
		// track number of contributions to edge_e()s
		array<array<array<int, 2>, 2>, 3> e_items;
		for (size_t d1: {0, 1, 2}) // parallel dim of edge
		for (size_t d2: {0, 1}) // side in 1st perpendicular dim of edge
		for (size_t d3: {0, 1}) { // side of cell in 2nd perp dim
			edge_e(d1,d2,d3) = 0.0;
			e_items[d1][d2][d3] = 0;
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
		if (cpface(-1) and cfinfo(-1) >= 0) {
			{const size_t d1 = 1, d2 = 0, Bd = 2;
			// y directed edge, -x side of cell center, -z side of cell center
			if (const size_t d3 = 0;
				// non-primary edges are handled by another cell
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) += Mag_fnx(*cell.data)[Bd];
				e_items[d1][d2][d3]++;
			}
			// y directed edge, -x side of cell center, +z side of cell center
			if (const size_t d3 = 1;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) += Mag_fnx(*cell.data)[Bd];
				e_items[d1][d2][d3]++;
			}}

			{const size_t d1 = 2, d2 = 0, Bd = 1;
			if (const size_t d3 = 0;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) -= Mag_fnx(*cell.data)[Bd];
				e_items[d1][d2][d3]++;
			}
			// z directed edge, -x side of cell center, +y side of cell center
			if (const size_t d3 = 1;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) -= Mag_fnx(*cell.data)[Bd];
				e_items[d1][d2][d3]++;
			}}
		}

		if (cpface(+1) and cfinfo(+1) >= 0) {
			{const size_t d1 = 1, d2 = 1, Bd = 2;
			if (const size_t d3 = 0;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) += Mag_fpx(*cell.data)[Bd];
				e_items[d1][d2][d3]++;
			}
			if (const size_t d3 = 1;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) += Mag_fpx(*cell.data)[Bd];
				e_items[d1][d2][d3]++;
			}}

			{const size_t d1 = 2, d2 = 1, Bd = 1;
			if (const size_t d3 = 0;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) -= Mag_fpx(*cell.data)[Bd];
				e_items[d1][d2][d3]++;
			}
			if (const size_t d3 = 1;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) -= Mag_fpx(*cell.data)[Bd];
				e_items[d1][d2][d3]++;
			}}
		}

		if (cpface(-2) and cfinfo(-2) >= 0) {
			{const size_t d1 = 0, d2 = 0, Bd = 2;
			if (const size_t d3 = 0;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) -= Mag_fny(*cell.data)[Bd];
				e_items[d1][d2][d3]++;
			}
			if (const size_t d3 = 1;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) -= Mag_fny(*cell.data)[Bd];
				e_items[d1][d2][d3]++;
			}}

			{const size_t d1 = 2, d3 = 0, Bd = 0;
			if (const size_t d2 = 0;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) += Mag_fny(*cell.data)[Bd];
				e_items[d1][d2][d3]++;
			}
			if (const size_t d2 = 1;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) += Mag_fny(*cell.data)[Bd];
				e_items[d1][d2][d3]++;
			}}
		}

		if (cpface(+2) and cfinfo(+2) >= 0) {
			{const size_t d1 = 0, d2 = 1, Bd = 2;
			if (const size_t d3 = 0;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) -= Mag_fpy(*cell.data)[Bd];
				e_items[d1][d2][d3]++;
			}
			if (const size_t d3 = 1;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) -= Mag_fpy(*cell.data)[Bd];
				e_items[d1][d2][d3]++;
			}}

			{const size_t d1 = 2, d3 = 1, Bd = 0;
			if (const size_t d2 = 0;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) += Mag_fpy(*cell.data)[Bd];
				e_items[d1][d2][d3]++;
			}
			if (const size_t d2 = 1;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) += Mag_fpy(*cell.data)[Bd];
				e_items[d1][d2][d3]++;
			}}
		}

		if (cpface(-3) and cfinfo(-3) >= 0) {
			{const size_t d1 = 0, d3 = 0, Bd = 1;
			if (const size_t d2 = 0;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) += Mag_fnz(*cell.data)[Bd];
				e_items[d1][d2][d3]++;
			}
			if (const size_t d2 = 1;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) += Mag_fnz(*cell.data)[Bd];
				e_items[d1][d2][d3]++;
			}}

			{const size_t d1 = 1, d3 = 0, Bd = 0;
			if (const size_t d2 = 0;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) -= Mag_fnz(*cell.data)[Bd];
				e_items[d1][d2][d3]++;
			}
			if (const size_t d2 = 1;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) -= Mag_fnz(*cell.data)[Bd];
				e_items[d1][d2][d3]++;
			}}
		}

		if (cpface(+3) and cfinfo(+3) >= 0) {
			{const size_t d1 = 0, d3 = 1, Bd = 1;
			if (const size_t d2 = 0;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) += Mag_fpz(*cell.data)[Bd];
				e_items[d1][d2][d3]++;
			}
			if (const size_t d2 = 1;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) += Mag_fpz(*cell.data)[Bd];
				e_items[d1][d2][d3]++;
			}}

			{const size_t d1 = 1, d3 = 1, Bd = 0;
			if (const size_t d2 = 0;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) -= Mag_fpz(*cell.data)[Bd];
				e_items[d1][d2][d3]++;
			}
			if (const size_t d2 = 1;
				cpedge(d1,d2,d3)
			) {
				edge_e(d1,d2,d3) -= Mag_fpz(*cell.data)[Bd];
				e_items[d1][d2][d3]++;
			}}
		}

		const auto [dx, dy, dz] = grid.geometry.get_length(cell.id);
		const auto cleni = grid.mapping.get_cell_length_in_indices(cell.id);

		for (const auto& neighbor: cell.neighbors_of) {
			const auto& fn = neighbor.face_neighbor;
			const auto& en = neighbor.edge_neighbor;
			if (fn == 0 and en[0] < 0) {
				continue;
			}

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
				constexpr size_t d2 = 0, Bd = 0;

				if (constexpr size_t d1 = 1, d3 = 0; cpedge(d1,d2,d3)) {
					if (npface(-3) and nfinfo(-3) >= 0 and neighbor.z == 0) {
						edge_e(d1,d2,d3) -= Mag_fnz(*neighbor.data)[Bd];
						e_items[d1][d2][d3]++;
					}
					if (cleni == neighbor.z + nleni) {
						throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					}
				}

				if (constexpr size_t d1 = 1, d3 = 1; cpedge(d1,d2,d3)) {
					if (npface(+3) and nfinfo(+3) >= 0 and cleni == neighbor.z + nleni) {
						edge_e(d1,d2,d3) -= Mag_fpz(*neighbor.data)[Bd];
						e_items[d1][d2][d3]++;
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
						e_items[d1][d2][d3]++;
					}
					if (neighbor.relative_size >= 0 and neighbor.z == 0) {
						throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					}
				}

				if (constexpr size_t d1 = 2, d3 = 0; cpedge(d1,d2,d3)) {
					if (npface(-2) and nfinfo(-2) >= 0 and neighbor.y == 0) {
						edge_e(d1,d2,d3) += Mag_fny(*neighbor.data)[Bd];
						e_items[d1][d2][d3]++;
					}
					if (cleni == neighbor.y + nleni) {
						throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					}
				}

				if (constexpr size_t d1 = 2, d3 = 1; cpedge(d1,d2,d3)) {
					if (npface(+2) and nfinfo(+2) >= 0 and cleni == neighbor.y + nleni) {
						edge_e(d1,d2,d3) += Mag_fpy(*neighbor.data)[Bd];
						e_items[d1][d2][d3]++;
					}
					if (
						neighbor.relative_size < 0
						and cfinfo(+2) >= 0
						and cpface(+2)
						and neighbor.y == 0
					) {
						edge_e(d1,d2,d3) += Mag_fpy(*cell.data)[Bd];
						e_items[d1][d2][d3]++;
					}
					if (neighbor.relative_size >= 0 and neighbor.y == 0) {
						throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					}
				}
			}

			if (fn == +1) {
				constexpr size_t d2 = 1, Bd = 0;

				if (constexpr size_t d1 = 2, d3 = 0; cpedge(d1,d2,d3)) {
					if (npface(-2) and nfinfo(-2) >= 0 and neighbor.y == 0) {
						edge_e(d1,d2,d3) += Mag_fny(*neighbor.data)[Bd];
						e_items[d1][d2][d3]++;
					}
				}

				if (constexpr size_t d1 = 2, d3 = 1; cpedge(d1,d2,d3)) {
					if (npface(+2) and nfinfo(+2) >= 0 and cleni == neighbor.y + nleni) {
						edge_e(d1,d2,d3) += Mag_fpy(*neighbor.data)[Bd];
						e_items[d1][d2][d3]++;
					}
					if (
						neighbor.relative_size < 0
						and cpface(+2)
						and cfinfo(+2) >= 0
						and neighbor.y == 0
					) {
						edge_e(d1,d2,d3) += Mag_fpy(*cell.data)[Bd];
						e_items[d1][d2][d3]++;
					}
				}

				if (constexpr size_t d1 = 1, d3 = 0; cpedge(d1,d2,d3)) {
					if (npface(-3) and nfinfo(-3) >= 0 and neighbor.z == 0) {
						edge_e(d1,d2,d3) -= Mag_fnz(*neighbor.data)[Bd];
						e_items[d1][d2][d3]++;
					}
				}

				if (constexpr size_t d1 = 1, d3 = 1; cpedge(d1,d2,d3)) {
					if (npface(+3) and nfinfo(+3) >= 0 and cleni == neighbor.z + nleni) {
						edge_e(d1,d2,d3) -= Mag_fpz(*neighbor.data)[Bd];
						e_items[d1][d2][d3]++;
					}
					if (
						neighbor.relative_size < 0
						and cpface(+3)
						and cfinfo(+3) >= 0
						and neighbor.z == 0
					) {
						edge_e(d1,d2,d3) -= Mag_fpz(*cell.data)[Bd];
						e_items[d1][d2][d3]++;
					}
				}
			}

			if (fn == -2) {
				constexpr size_t Bd = 1;

				if (constexpr size_t d1 = 2, d2 = 0, d3 = 0; cpedge(d1,d2,d3)) {
					if (npface(-1) and nfinfo(-1) >= 0 and neighbor.x == 0) {
						edge_e(d1,d2,d3) -= Mag_fnx(*neighbor.data)[Bd];
						e_items[d1][d2][d3]++;
					}
				}

				if (constexpr size_t d1 = 2, d2 = 1, d3 = 0; cpedge(d1,d2,d3)) {
					if (nfinfo(+1) >= 0 and npface(+1) and cleni == neighbor.x + nleni) {
						edge_e(d1,d2,d3) -= Mag_fpx(*neighbor.data)[Bd];
						e_items[d1][d2][d3]++;
					}
					if (
						neighbor.relative_size < 0
						and cpface(+1)
						and cfinfo(+1) >= 0
						and neighbor.x == 0
					) {
						edge_e(d1,d2,d3) -= Mag_fpx(*cell.data)[Bd];
						e_items[d1][d2][d3]++;
					}
				}

				if (constexpr size_t d1 = 0, d2 = 0, d3 = 0; cpedge(d1,d2,d3)) {
					if (npface(-3) and nfinfo(-3) >= 0 and neighbor.z == 0) {
						edge_e(d1,d2,d3) += Mag_fnz(*neighbor.data)[Bd];
						e_items[d1][d2][d3]++;
					}
				}

				if (constexpr size_t d1 = 0, d2 = 0, d3 = 1; cpedge(d1,d2,d3)) {
					if (npface(+3) and nfinfo(+3) >= 0 and cleni == neighbor.z + nleni) {
						edge_e(d1,d2,d3) += Mag_fpz(*neighbor.data)[Bd];
						e_items[d1][d2][d3]++;
					}
					if (
						neighbor.relative_size < 0
						and cpface(+3)
						and cfinfo(+3) >= 0
						and neighbor.z == 0
					) {
						edge_e(d1,d2,d3) += Mag_fpz(*cell.data)[Bd];
						e_items[d1][d2][d3]++;
					}
				}
			}

			if (fn == +2) {
				constexpr size_t Bd = 1;

				if (constexpr size_t d1 = 2, d2 = 0, d3 = 1; cpedge(d1,d2,d3)) {
					if (nfinfo(-1) >= 0 and npface(-1) and neighbor.x == 0) {
						edge_e(d1,d2,d3) -= Mag_fnx(*neighbor.data)[Bd];
						e_items[d1][d2][d3]++;
					}
				}

				if (constexpr size_t d1 = 2, d2 = 1, d3 = 1; cpedge(d1,d2,d3)) {
					if (npface(+1) and nfinfo(+1) >= 0 and cleni == neighbor.x + nleni) {
						edge_e(d1,d2,d3) -= Mag_fpx(*neighbor.data)[Bd];
						e_items[d1][d2][d3]++;
					}
					if (
						neighbor.relative_size < 0
						and cpface(+1)
						and cfinfo(+1) >= 0
						and neighbor.x == 0
					) {
						edge_e(d1,d2,d3) -= Mag_fpx(*cell.data)[Bd];
						e_items[d1][d2][d3]++;
					}
				}

				if (constexpr size_t d1 = 0, d2 = 1, d3 = 0; cpedge(d1,d2,d3)) {
					if (nfinfo(-3) >= 0 and npface(-3) and neighbor.z == 0) {
						edge_e(d1,d2,d3) += Mag_fnz(*neighbor.data)[Bd];
						e_items[d1][d2][d3]++;
					}
				}

				if (constexpr size_t d1 = 0, d2 = 1, d3 = 1; cpedge(d1,d2,d3)) {
					if (npface(+3) and nfinfo(+3) >= 0 and cleni == neighbor.z + nleni) {
						edge_e(d1,d2,d3) += Mag_fpz(*neighbor.data)[Bd];
						e_items[d1][d2][d3]++;
					}
					if (
						neighbor.relative_size < 0
						and cpface(+3)
						and cfinfo(+3) >= 0
						and neighbor.z == 0
					) {
						edge_e(d1,d2,d3) += Mag_fpz(*cell.data)[Bd];
						e_items[d1][d2][d3]++;
					}
				}
			}

			if (fn == -3) {
				constexpr size_t d3 = 0, Bd = 2;

				if (constexpr size_t d1 = 1, d2 = 0; cpedge(d1,d2,d3)) {
					if (npface(-1) and nfinfo(-1) >= 0 and neighbor.x == 0) {
						edge_e(d1,d2,d3) += Mag_fnx(*neighbor.data)[Bd];
						e_items[d1][d2][d3]++;
					}
				}

				if (constexpr size_t d1 = 1, d2 = 1; cpedge(d1,d2,d3)) {
					if (npface(+1) and nfinfo(+1) >= 0 and cleni == neighbor.x + nleni) {
						edge_e(d1,d2,d3) += Mag_fpx(*neighbor.data)[Bd];
						e_items[d1][d2][d3]++;
					}
					if (
						neighbor.relative_size < 0
						and cpface(+1)
						and cfinfo(+1) >= 0
						and neighbor.x == 0
					) {
						edge_e(d1,d2,d3) += Mag_fpx(*cell.data)[Bd];
						e_items[d1][d2][d3]++;
					}
				}

				if (constexpr size_t d1 = 0, d2 = 0; cpedge(d1,d2,d3)) {
					if (npface(-2) and nfinfo(-2) >= 0 and neighbor.y == 0) {
						edge_e(d1,d2,d3) -= Mag_fny(*neighbor.data)[Bd];
						e_items[d1][d2][d3]++;
					}
				}

				if (constexpr size_t d1 = 0, d2 = 1; cpedge(d1,d2,d3)) {
					if (npface(+2) and nfinfo(+2) >= 0 and cleni == neighbor.y + nleni) {
						edge_e(d1,d2,d3) -= Mag_fpy(*neighbor.data)[Bd];
						e_items[d1][d2][d3]++;
					}
					if (
						neighbor.relative_size < 0
						and cpface(+2)
						and cfinfo(+2) >= 0
						and neighbor.y == 0
					) {
						edge_e(d1,d2,d3) -= Mag_fpy(*cell.data)[Bd];
						e_items[d1][d2][d3]++;
					}
				}
			}

			if (fn == +3) {
				constexpr size_t d3 = 1, Bd = 2;

				if (constexpr size_t d1 = 1, d2 = 0; cpedge(d1,d2,d3)) {
					if (npface(-1) and nfinfo(-1) >= 0 and neighbor.x == 0) {
						edge_e(d1,d2,d3) += Mag_fnx(*neighbor.data)[Bd];
						e_items[d1][d2][d3]++;
					}
				}

				if (constexpr size_t d1 = 1, d2 = 1; cpedge(d1,d2,d3)) {
					if (npface(+1) and nfinfo(+1) >= 0 and cleni == neighbor.x + nleni) {
						edge_e(d1,d2,d3) += Mag_fpx(*neighbor.data)[Bd];
						e_items[d1][d2][d3]++;
					}
					if (
						neighbor.relative_size < 0
						and cpface(+1)
						and cfinfo(+1) >= 0
						and neighbor.x == 0
					) {
						edge_e(d1,d2,d3) += Mag_fpx(*cell.data)[Bd];
						e_items[d1][d2][d3]++;
					}
				}

				if (constexpr size_t d1 = 0, d2 = 0; cpedge(d1,d2,d3)) {
					if (npface(-2) and nfinfo(-2) >= 0 and neighbor.y == 0) {
						edge_e(d1,d2,d3) -= Mag_fny(*neighbor.data)[Bd];
						e_items[d1][d2][d3]++;
					}
				}

				if (constexpr size_t d1 = 0, d2 = 1; cpedge(d1,d2,d3)) {
					if (npface(+2) and nfinfo(+2) >= 0 and cleni == neighbor.y + nleni) {
						edge_e(d1,d2,d3) -= Mag_fpy(*neighbor.data)[Bd];
						e_items[d1][d2][d3]++;
					}
					if (
						neighbor.relative_size < 0
						and cpface(+2)
						and cfinfo(+2) >= 0
						and neighbor.y == 0
					) {
						edge_e(d1,d2,d3) -= Mag_fpy(*cell.data)[Bd];
						e_items[d1][d2][d3]++;
					}
				}
			}

			if (constexpr size_t d1 = 0, d2 = 0, d3 = 0;
				cpedge(d1,d2,d3) and en[0] == d1 and en[1] == d2 and en[2] == d3
			) {
				if (nfinfo(+2) >= 0 and npface(+2)) {
					edge_e(d1,d2,d3) -= Mag_fpy(*neighbor.data)[2];
					e_items[d1][d2][d3]++;
				}
				if (nfinfo(+3) >= 0 and npface(+3)) {
					edge_e(d1,d2,d3) += Mag_fpz(*neighbor.data)[1];
					e_items[d1][d2][d3]++;
				}
			}

			if (constexpr size_t d1 = 0, d2 = 0, d3 = 1;
				cpedge(d1,d2,d3) and en[0] == d1 and en[1] == d2
			) {
				if (en[2] == d3) {
					if (nfinfo(+2) >= 0 and npface(+2)) {
						edge_e(d1,d2,d3) -= Mag_fpy(*neighbor.data)[2];
						e_items[d1][d2][d3]++;
					}
					if (nfinfo(-3) >= 0 and npface(-3)) {
						edge_e(d1,d2,d3) += Mag_fnz(*neighbor.data)[1];
						e_items[d1][d2][d3]++;
					}
				}
			}

			if (constexpr size_t d1 = 0, d2 = 1, d3 = 0;
				cpedge(d1,d2,d3) and en[0] == d1 and en[2] == d3
			) {
				if (en[1] == d2) {
					if (nfinfo(-2) >= 0 and npface(-2)) {
						edge_e(d1,d2,d3) -= Mag_fny(*neighbor.data)[2];
						e_items[d1][d2][d3]++;
					}
					if (nfinfo(+3) >= 0 and npface(+3)) {
						edge_e(d1,d2,d3) += Mag_fpz(*neighbor.data)[1];
						e_items[d1][d2][d3]++;
					}
				}
			}

			if (constexpr size_t d1 = 0, d2 = 1, d3 = 1;
				cpedge(d1,d2,d3) and en[0] == d1
			) {
				if (en[1] == d2 and en[2] == d3) {
					if (nfinfo(-2) >= 0 and npface(-2)) {
						edge_e(d1,d2,d3) -= Mag_fny(*neighbor.data)[2];
						e_items[d1][d2][d3]++;
					}
					if (nfinfo(-3) >= 0 and npface(-3)) {
						edge_e(d1,d2,d3) += Mag_fnz(*neighbor.data)[1];
						e_items[d1][d2][d3]++;
					}
				}
			}

			if (constexpr size_t d1 = 1, d2 = 0, d3 = 0;
				cpedge(d1,d2,d3) and en[0] == d1 and en[1] == d2 and en[2] == d3
			) {
				if (nfinfo(+1) >= 0 and npface(+1)) {
					edge_e(d1,d2,d3) += Mag_fpx(*neighbor.data)[2];
					e_items[d1][d2][d3]++;
				}
				if (nfinfo(+3) >= 0 and npface(+3)) {
					edge_e(d1,d2,d3) -= Mag_fpz(*neighbor.data)[0];
					e_items[d1][d2][d3]++;
				}
			}

			if (constexpr size_t d1 = 1, d2 = 0, d3 = 1;
				cpedge(d1,d2,d3) and en[0] == d1 and en[1] == d2
			) {
				if (en[2] == d3) {
					if (nfinfo(+1) >= 0 and npface(+1)) {
						edge_e(d1,d2,d3) += Mag_fpx(*neighbor.data)[2];
						e_items[d1][d2][d3]++;
					}
					if (nfinfo(-3) >= 0 and npface(-3)) {
						edge_e(d1,d2,d3) -= Mag_fnz(*neighbor.data)[0];
						e_items[d1][d2][d3]++;
					}
				}
			}

			if (constexpr size_t d1 = 1, d2 = 1, d3 = 0;
				cpedge(d1,d2,d3) and en[0] == d1 and en[2] == d3
			) {
				if (en[1] == d2) {
					if (nfinfo(-1) >= 0 and npface(-1)) {
						edge_e(d1,d2,d3) += Mag_fnx(*neighbor.data)[2];
						e_items[d1][d2][d3]++;
					}
					if (nfinfo(+3) >= 0 and npface(+3)) {
						edge_e(d1,d2,d3) -= Mag_fpz(*neighbor.data)[0];
						e_items[d1][d2][d3]++;
					}
				}
			}

			if (constexpr size_t d1 = 1, d2 = 1, d3 = 1;
				cpedge(d1,d2,d3) and en[0] == d1
			) {
				if (en[1] == d2 and en[2] == d3) {
					if (nfinfo(-1) >= 0 and npface(-1)) {
						edge_e(d1,d2,d3) += Mag_fnx(*neighbor.data)[2];
						e_items[d1][d2][d3]++;
					}
					if (nfinfo(-3) >= 0 and npface(-3)) {
						edge_e(d1,d2,d3) -= Mag_fnz(*neighbor.data)[0];
						e_items[d1][d2][d3]++;
					}
				}
			}

			if (constexpr size_t d1 = 2, d2 = 0, d3 = 0;
				cpedge(d1,d2,d3) and en[0] == d1 and en[1] == d2 and en[2] == d3
			) {
				if (nfinfo(+1) >= 0 and npface(+1)) {
					edge_e(d1,d2,d3) -= Mag_fpx(*neighbor.data)[1];
					e_items[d1][d2][d3]++;
				}
				if (nfinfo(+2) >= 0 and npface(+2)) {
					edge_e(d1,d2,d3) += Mag_fpy(*neighbor.data)[0];
					e_items[d1][d2][d3]++;
				}
			}

			if (constexpr size_t d1 = 2, d2 = 0, d3 = 1;
				cpedge(d1,d2,d3) and en[0] == d1 and en[1] == d2
			) {
				if (en[2] == d3) {
					if (nfinfo(+1) >= 0 and npface(+1)) {
						edge_e(d1,d2,d3) -= Mag_fpx(*neighbor.data)[1];
						e_items[d1][d2][d3]++;
					}
					if (nfinfo(-2) >= 0 and npface(-2)) {
						edge_e(d1,d2,d3) += Mag_fny(*neighbor.data)[0];
						e_items[d1][d2][d3]++;
					}
				}
			}

			if (constexpr size_t d1 = 2, d2 = 1, d3 = 0;
				cpedge(d1,d2,d3) and en[0] == d1 and en[2] == d3
			) {
				if (en[1] == d2) {
					if (nfinfo(-1) >= 0 and npface(-1)) {
						edge_e(d1,d2,d3) -= Mag_fnx(*neighbor.data)[1];
						e_items[d1][d2][d3]++;
					}
					if (nfinfo(+2) >= 0 and npface(+2)) {
						edge_e(d1,d2,d3) += Mag_fpy(*neighbor.data)[0];
						e_items[d1][d2][d3]++;
					}
				}
			}

			if (constexpr size_t d1 = 2, d2 = 1, d3 = 1;
				cpedge(d1,d2,d3) and en[0] == d1
			) {
				if (en[1] == d2 and en[2] == d3) {
					if (nfinfo(-1) >= 0 and npface(-1)) {
						edge_e(d1,d2,d3) -= Mag_fnx(*neighbor.data)[1];
						e_items[d1][d2][d3]++;
					}
					if (nfinfo(-2) >= 0 and npface(-2)) {
						edge_e(d1,d2,d3) += Mag_fny(*neighbor.data)[0];
						e_items[d1][d2][d3]++;
					}
				}
			}

			// these shouldn't be possible but kept for reference
			if (fn == -1 and npface(+1)) {
				if (cpedge(1,0,0)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(1,0,0) += get<1>(Mag_f)(*neighbor.data)[2];
					//e_items[1][0][0]++;
				}
				if (cpedge(1,0,1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(1,0,1) += get<1>(Mag_f)(*neighbor.data)[2];
					//e_items[1][0][1]++;
				}
				if (cpedge(2,0,0)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(2,0,0) -= get<1>(Mag_f)(*neighbor.data)[1];
					//e_items[2][0][0]++;
				}
				if (cpedge(2,0,1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(2,0,1) -= get<1>(Mag_f)(*neighbor.data)[1];
					//e_items[2][0][1]++;
				}
			}
			if (fn == +1 and npface(-1)) {
				if (cpedge(1,1,0)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(1,1,0) += get<0>(Mag_f)(*neighbor.data)[2];
					//e_items[1][1][0]++;
				}
				if (cpedge(1,1,1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(1,1,1) += get<0>(Mag_f)(*neighbor.data)[2];
					//e_items[1][1][1]++;
				}
				if (cpedge(2,1,0)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(2,1,0) -= get<0>(Mag_f)(*neighbor.data)[1];
					//e_items[2][1][0]++;
				}
				if (cpedge(2,1,1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(2,1,1) -= get<0>(Mag_f)(*neighbor.data)[1];
					//e_items[2][1][1]++;
				}
			}
			if (fn == -2 and npface(+2)) {
				if (cpedge(0,0,0)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(0,0,0) -= get<3>(Mag_f)(*neighbor.data)[2];
					//e_items[0][0][0]++;
				}
				if (cpedge(0,0,1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(0,0,1) -= get<3>(Mag_f)(*neighbor.data)[2];
					//e_items[0][0][1]++;
				}
				if (cpedge(2,0,0)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(2,0,0) += get<3>(Mag_f)(*neighbor.data)[0];
					//e_items[2][0][0]++;
				}
				if (cpedge(2,1,0)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(2,1,0) += get<3>(Mag_f)(*neighbor.data)[0];
					//e_items[2][1][0]++;
				}
			}
			if (fn == +2 and npface(-2)) {
				if (cpedge(0,1,0)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(0,1,0) -= get<2>(Mag_f)(*neighbor.data)[2];
					//e_items[0][1][0]++;
				}
				if (cpedge(0,1,1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(0,1,1) -= get<2>(Mag_f)(*neighbor.data)[2];
					//e_items[0][1][1]++;
				}
				if (cpedge(2,0,1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(2,0,1) += get<2>(Mag_f)(*neighbor.data)[0];
					//e_items[2][0][1]++;
				}
				if (cpedge(2,1,1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(2,1,1) += get<2>(Mag_f)(*neighbor.data)[0];
					//e_items[2][1][1]++;
				}
			}
			if (fn == -3 and npface(+3)) {
				if (cpedge(0,0,0)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(0,0,0) += get<5>(Mag_f)(*neighbor.data)[1];
					//e_items[0][0][0]++;
				}
				if (cpedge(0,1,0)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(0,1,0) += get<5>(Mag_f)(*neighbor.data)[1];
					//e_items[0][1][0]++;
				}
				if (cpedge(1,0,0)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(1,0,0) -= get<5>(Mag_f)(*neighbor.data)[0];
					//e_items[1][0][0]++;
				}
				if (cpedge(1,1,0)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(1,1,0) -= get<5>(Mag_f)(*neighbor.data)[0];
					//e_items[1][1][0]++;
				}
			}
			if (fn == +3 and npface(-3)) {
				if (cpedge(0,0,1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(0,0,1) += get<4>(Mag_f)(*neighbor.data)[1];
					//e_items[0][0][1]++;
				}
				if (cpedge(0,1,1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(0,1,1) += get<4>(Mag_f)(*neighbor.data)[1];
					//e_items[0][1][1]++;
				}
				if (cpedge(1,0,1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(1,0,1) -= get<4>(Mag_f)(*neighbor.data)[0];
					//e_items[1][0][1]++;
				}
				if (cpedge(1,1,1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//edge_e(1,1,1) -= get<4>(Mag_f)(*neighbor.data)[0];
					//e_items[1][1][1]++;
				}
			}
		}

		for (size_t d1: {0, 1, 2})
		for (size_t d2: {0, 1})
		for (size_t d3: {0, 1}) {
			if (not cpedge(d1,d2,d3)) {
				if (e_items[d1][d2][d3] != 0) throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
				continue;
			}
			if (e_items[d1][d2][d3] > 0) edge_e(d1,d2,d3) /= e_items[d1][d2][d3];
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
	throw std::runtime_error(__func__ + std::string(": ") + e.what());
}


/*! Solves new face magnetic fields from edge electric fields.

Equations 13-15 of https://doi.org/10.1006/jcph.1998.6153

Saves new face B of cells with SInfo(*cell_data) == 1,
ignores cells with SInfo < 0.
*/
template <
	class Cells,
	class Grid,
	class Face_Magnetic_Field_Pos_Getter,
	class Face_Magnetic_Field_Neg_Getter,
	class Edge_Electric_Field_Getter,
	class Primary_Face_Getter,
	class Primary_Edge_Getter,
	class Face_Info_Getter
> void get_face_magnetic_field(
	const Cells& cells,
	Grid& grid,
	const double dt,
	const Face_Magnetic_Field_Pos_Getter Face_B_pos,
	const Face_Magnetic_Field_Neg_Getter Face_B_neg,
	const Edge_Electric_Field_Getter Edge_E,
	const Primary_Face_Getter PFace,
	const Primary_Edge_Getter PEdge,
	const Face_Info_Getter FInfo
) try {
	using std::runtime_error;
	using std::to_string;

	for (const auto& cell: cells) {
		const auto cell_length = grid.geometry.get_length(cell.id);

		typename std::remove_reference<decltype(Face_B_pos(*cell.data))>::type face_db_p{0, 0, 0}, face_db_n{0, 0, 0};

		const auto& cpface = PFace(*cell.data);
		const auto& cfinfo = FInfo(*cell.data);
		const auto& cpedge = PEdge(*cell.data);
		const auto& cedge_e = Edge_E(*cell.data);

		/*! Contributions from edge electric fields to face Bs

		Right-hand line integral over edges touching each face:
		dBx: E(2,*,1)dz - E(2,*,0)dz - E(1,*,1)dx + E(1,*,0)dx
		dBy: E(0,*,1)dx - E(0,*,0)dx - E(2,1,*)dz + E(2,0,*)dz
		dBz: E(1,1,*)dy - E(1,0,*)dy + E(0,0,*)dx - E(0,1,*)dx
		*/
		if (cfinfo(0,-1) >= 0 and cpface(0,-1)) {
			if (cpedge(1,0,0)) face_db_n[0] += cell_length[1]*cedge_e(1,0,0);
			if (cpedge(1,0,1)) face_db_n[0] -= cell_length[1]*cedge_e(1,0,1);
			if (cpedge(2,0,0)) face_db_n[0] -= cell_length[2]*cedge_e(2,0,0);
			if (cpedge(2,0,1)) face_db_n[0] += cell_length[2]*cedge_e(2,0,1);
		}
		if (cfinfo(0,+1) >= 0 and cpface(0,+1)) {
			if (cpedge(1,1,0)) face_db_p[0] += cell_length[1]*cedge_e(1,1,0);
			if (cpedge(1,1,1)) face_db_p[0] -= cell_length[1]*cedge_e(1,1,1);
			if (cpedge(2,1,0)) face_db_p[0] -= cell_length[2]*cedge_e(2,1,0);
			if (cpedge(2,1,1)) face_db_p[0] += cell_length[2]*cedge_e(2,1,1);
		}
		if (cfinfo(1,-1) >= 0 and cpface(1,-1)) {
			if (cpedge(0,0,0)) face_db_n[1] -= cell_length[0]*cedge_e(0,0,0);
			if (cpedge(0,0,1)) face_db_n[1] += cell_length[0]*cedge_e(0,0,1);
			if (cpedge(2,0,0)) face_db_n[1] += cell_length[2]*cedge_e(2,0,0);
			if (cpedge(2,1,0)) face_db_n[1] -= cell_length[2]*cedge_e(2,1,0);
		}
		if (cfinfo(1,+1) >= 0 and cpface(1,+1)) {
			if (cpedge(0,1,0)) face_db_p[1] -= cell_length[0]*cedge_e(0,1,0);
			if (cpedge(0,1,1)) face_db_p[1] += cell_length[0]*cedge_e(0,1,1);
			if (cpedge(2,0,1)) face_db_p[1] += cell_length[2]*cedge_e(2,0,1);
			if (cpedge(2,1,1)) face_db_p[1] -= cell_length[2]*cedge_e(2,1,1);
		}
		if (cfinfo(2,-1) >= 0 and cpface(2,-1)) {
			if (cpedge(0,0,0)) face_db_n[2] += cell_length[0]*cedge_e(0,0,0);
			if (cpedge(0,1,0)) face_db_n[2] -= cell_length[0]*cedge_e(0,1,0);
			if (cpedge(1,0,0)) face_db_n[2] -= cell_length[1]*cedge_e(1,0,0);
			if (cpedge(1,1,0)) face_db_n[2] += cell_length[1]*cedge_e(1,1,0);
		}
		if (cfinfo(2,+1) >= 0 and cpface(2,+1)) {
			if (cpedge(0,0,1)) face_db_p[2] += cell_length[0]*cedge_e(0,0,1);
			if (cpedge(0,1,1)) face_db_p[2] -= cell_length[0]*cedge_e(0,1,1);
			if (cpedge(1,0,1)) face_db_p[2] -= cell_length[1]*cedge_e(1,0,1);
			if (cpedge(1,1,1)) face_db_p[2] += cell_length[1]*cedge_e(1,1,1);
		}

		const auto cleni = grid.mapping.get_cell_length_in_indices(cell.id);

		for (const auto& neighbor: cell.neighbors_of) {
			const auto& fn = neighbor.face_neighbor;
			const auto& en = neighbor.edge_neighbor;

			if (fn == 0 and en[0] < 0) {
				continue;
			}

			const auto neigh_length = grid.geometry.get_length(neighbor.id);
			const auto nleni = grid.mapping.get_cell_length_in_indices(neighbor.id);

			const auto& npedge = PEdge(*neighbor.data);
			const auto& nedge_e = Edge_E(*neighbor.data);

			if (cfinfo(-1) >= 0 and cpface(-1) and neighbor.x == 0) {
				if (fn == -2 and npedge(2,0,1)) face_db_n[0] -= neigh_length[2]*nedge_e(2,0,1);
				if (fn == +2 and npedge(2,0,0)) face_db_n[0] += neigh_length[2]*nedge_e(2,0,0);
				if (fn == -3 and npedge(1,0,1)) face_db_n[0] += neigh_length[1]*nedge_e(1,0,1);
				if (fn == +3 and npedge(1,0,0)) face_db_n[0] -= neigh_length[1]*nedge_e(1,0,0);
			}
			if (cfinfo(+1) >= 0 and cpface(+1) and cleni == neighbor.x + nleni) {
				if (fn == -2 and npedge(2,1,1)) face_db_p[0] -= neigh_length[2]*nedge_e(2,1,1);
				if (fn == +2 and npedge(2,1,0)) face_db_p[0] += neigh_length[2]*nedge_e(2,1,0);
				if (fn == -3 and npedge(1,1,1)) face_db_p[0] += neigh_length[1]*nedge_e(1,1,1);
				if (fn == +3 and npedge(1,1,0)) face_db_p[0] -= neigh_length[1]*nedge_e(1,1,0);
			}
			if (cfinfo(-2) >= 0 and cpface(-2) and neighbor.y == 0) {
				if (fn == -1 and npedge(2,1,0)) face_db_n[1] += neigh_length[2]*nedge_e(2,1,0);
				if (fn == +1 and npedge(2,0,0)) face_db_n[1] -= neigh_length[2]*nedge_e(2,0,0);
				if (fn == -3 and npedge(0,0,1)) face_db_n[1] -= neigh_length[0]*nedge_e(0,0,1);
				if (fn == +3 and npedge(0,0,0)) face_db_n[1] += neigh_length[0]*nedge_e(0,0,0);
			}
			if (cfinfo(+2) >= 0 and cpface(+2) and cleni == neighbor.y + nleni) {
				if (fn == -1 and npedge(2,1,1)) face_db_p[1] += neigh_length[2]*nedge_e(2,1,1);
				if (fn == +1 and npedge(2,0,1)) face_db_p[1] -= neigh_length[2]*nedge_e(2,0,1);
				if (fn == -3 and npedge(0,1,1)) face_db_p[1] -= neigh_length[0]*nedge_e(0,1,1);
				if (fn == +3 and npedge(0,1,0)) face_db_p[1] += neigh_length[0]*nedge_e(0,1,0);
			}
			if (cfinfo(-3) >= 0 and cpface(-3) and neighbor.z == 0) {
				if (fn == -1 and npedge(1,1,0)) face_db_n[2] -= neigh_length[1]*nedge_e(1,1,0);
				if (fn == +1 and npedge(1,0,0)) face_db_n[2] += neigh_length[1]*nedge_e(1,0,0);
				if (fn == -2 and npedge(0,1,0)) face_db_n[2] += neigh_length[0]*nedge_e(0,1,0);
				if (fn == +2 and npedge(0,0,0)) face_db_n[2] -= neigh_length[0]*nedge_e(0,0,0);
			}
			if (cfinfo(+3) >= 0 and cpface(+3) and cleni == neighbor.z + nleni) {
				if (fn == -1 and npedge(1,1,1)) face_db_p[2] -= neigh_length[1]*nedge_e(1,1,1);
				if (fn == +1 and npedge(1,0,1)) face_db_p[2] += neigh_length[1]*nedge_e(1,0,1);
				if (fn == -2 and npedge(0,1,1)) face_db_p[2] += neigh_length[0]*nedge_e(0,1,1);
				if (fn == +2 and npedge(0,0,1)) face_db_p[2] -= neigh_length[0]*nedge_e(0,0,1);
			}
			if (en[0] == 0 and en[1] == 0 and en[2] == 0 and npedge(0,1,1)) {
				if (cpedge(0,0,0)) throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
				if (cfinfo(-2) >= 0 and cpface(-2)) face_db_n[1] -= neigh_length[0]*nedge_e(0,1,1);
				if (cfinfo(-3) >= 0 and cpface(-3)) face_db_n[2] += neigh_length[0]*nedge_e(0,1,1);
			}
			if (en[0] == 0 and en[1] == 0 and en[2] == 1 and npedge(0,1,0)) {
				if (cpedge(0,0,1)) throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
				if (cfinfo(-2) >= 0 and cpface(-2)) face_db_n[1] += neigh_length[0]*nedge_e(0,1,0);
				if (cfinfo(+3) >= 0 and cpface(+3)) face_db_p[2] += neigh_length[0]*nedge_e(0,1,0);
			}
			if (en[0] == 0 and en[1] == 1 and en[2] == 0 and npedge(0,0,1)) {
				if (cpedge(0,1,0)) throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
				if (cfinfo(+2) >= 0 and cpface(+2)) face_db_p[1] -= neigh_length[0]*nedge_e(0,0,1);
				if (cfinfo(-3) >= 0 and cpface(-3)) face_db_n[2] -= neigh_length[0]*nedge_e(0,0,1);
			}
			if (en[0] == 0 and en[1] == 1 and en[2] == 1 and npedge(0,0,0)) {
				if (cpedge(0,1,1)) throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
				if (cfinfo(+2) >= 0 and cpface(+2)) face_db_p[1] += neigh_length[0]*nedge_e(0,0,0);
				if (cfinfo(+3) >= 0 and cpface(+3)) face_db_p[2] -= neigh_length[0]*nedge_e(0,0,0);
			}
			if (en[0] == 1 and en[1] == 0 and en[2] == 0 and npedge(1,1,1)) {
				if (cpedge(1,0,0)) throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
				if (cfinfo(-1) >= 0 and cpface(-1)) face_db_n[0] += neigh_length[1]*nedge_e(1,1,1);
				if (cfinfo(-3) >= 0 and cpface(-3)) face_db_n[2] -= neigh_length[1]*nedge_e(1,1,1);
			}
			if (en[0] == 1 and en[1] == 0 and en[2] == 1 and npedge(1,1,0)) {
				if (cpedge(1,0,1)) throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
				if (cfinfo(-1) >= 0 and cpface(-1)) face_db_n[0] -= neigh_length[1]*nedge_e(1,1,0);
				if (cfinfo(+3) >= 0 and cpface(+3)) face_db_p[2] -= neigh_length[1]*nedge_e(1,1,0);
			}
			if (en[0] == 1 and en[1] == 1 and en[2] == 0 and npedge(1,0,1)) {
				if (cpedge(1,1,0)) throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
				if (cfinfo(+1) >= 0 and cpface(+1)) face_db_p[0] += neigh_length[1]*nedge_e(1,0,1);
				if (cfinfo(-3) >= 0 and cpface(-3)) face_db_n[2] += neigh_length[1]*nedge_e(1,0,1);
			}
			if (en[0] == 1 and en[1] == 1 and en[2] == 1 and npedge(1,0,0)) {
				if (cpedge(1,1,1)) throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
				if (cfinfo(+1) >= 0 and cpface(+1)) face_db_p[0] -= neigh_length[1]*nedge_e(1,0,0);
				if (cfinfo(+3) >= 0 and cpface(+3)) face_db_p[2] += neigh_length[1]*nedge_e(1,0,0);
			}
			if (en[0] == 2 and en[1] == 0 and en[2] == 0 and npedge(2,1,1)) {
				if (cpedge(2,0,0)) throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
				if (cfinfo(-1) >= 0 and cpface(-1)) face_db_n[0] -= neigh_length[2]*nedge_e(2,1,1);
				if (cfinfo(-2) >= 0 and cpface(-2)) face_db_n[1] += neigh_length[2]*nedge_e(2,1,1);
			}
			if (en[0] == 2 and en[1] == 0 and en[2] == 1 and npedge(2,1,0)) {
				if (cpedge(2,0,1)) throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
				if (cfinfo(-1) >= 0 and cpface(-1)) face_db_n[0] += neigh_length[2]*nedge_e(2,1,0);
				if (cfinfo(+2) >= 0 and cpface(+2)) face_db_p[1] += neigh_length[2]*nedge_e(2,1,0);
			}
			if (en[0] == 2 and en[1] == 1 and en[2] == 0 and npedge(2,0,1)) {
				if (cpedge(2,1,0)) throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
				if (cfinfo(+1) >= 0 and cpface(+1)) face_db_p[0] -= neigh_length[2]*nedge_e(2,0,1);
				if (cfinfo(-2) >= 0 and cpface(-2)) face_db_n[1] -= neigh_length[2]*nedge_e(2,0,1);
			}
			if (en[0] == 2 and en[1] == 1 and en[2] == 1 and npedge(2,0,0)) {
				if (cpedge(2,1,1)) throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
				if (cfinfo(+1) >= 0 and cpface(+1)) face_db_p[0] += neigh_length[2]*nedge_e(2,0,0);
				if (cfinfo(+2) >= 0 and cpface(+2)) face_db_p[1] -= neigh_length[2]*nedge_e(2,0,0);
			}
			// these shouldn't be possible but kept for reference
			if (cpface(-1) and fn == -1) {
				if (npedge(1,1,0)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db_n[0] += neigh_length[1]*nedge_e(1,1,0);
				}
				if (npedge(1,1,1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db_n[0] -= neigh_length[1]*nedge_e(1,1,1);
				}
				if (npedge(2,1,0)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db_n[0] -= neigh_length[2]*nedge_e(2,1,0);
				}
				if (npedge(2,1,1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db_n[0] += neigh_length[2]*nedge_e(2,1,1);
				}
			}
			if (cpface(+1) and fn == +1) {
				if (npedge(1,0,0)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db_p[0] += neigh_length[1]*nedge_e(1,0,0);
				}
				if (npedge(1,0,1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db_p[0] -= neigh_length[1]*nedge_e(1,0,1);
				}
				if (npedge(2,0,0)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db_p[0] -= neigh_length[2]*nedge_e(2,0,0);
				}
				if (npedge(2,0,1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db_p[0] += neigh_length[2]*nedge_e(2,0,1);
				}
			}
			if (cpface(-2) and fn == -2) {
				if (npedge(0,1,0)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db_n[1] -= neigh_length[0]*nedge_e(0,1,0);
				}
				if (npedge(0,1,1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db_n[1] += neigh_length[0]*nedge_e(0,1,1);
				}
				if (npedge(2,0,1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db_n[1] += neigh_length[2]*nedge_e(2,0,1);
				}
				if (npedge(2,1,1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db_n[1] -= neigh_length[2]*nedge_e(2,1,1);
				}
			}
			if (cpface(+2) and fn == +2) {
				if (npedge(0,0,0)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db_p[1] -= neigh_length[0]*nedge_e(0,0,0);
				}
				if (npedge(0,0,1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db_p[1] += neigh_length[0]*nedge_e(0,0,1);
				}
				if (npedge(2,0,0)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db_p[1] += neigh_length[2]*nedge_e(2,0,0);
				}
				if (npedge(2,1,0)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db_p[1] -= neigh_length[2]*nedge_e(2,1,0);
				}
			}
			if (cpface(-3) and fn == -3) {
				if (npedge(0,0,1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db_n[2] += neigh_length[0]*nedge_e(0,0,1);
				}
				if (npedge(0,1,1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db_n[2] -= neigh_length[0]*nedge_e(0,1,1);
				}
				if (npedge(1,0,1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db_n[2] -= neigh_length[1]*nedge_e(1,0,1);
				}
				if (npedge(1,1,1)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db_n[2] += neigh_length[1]*nedge_e(1,1,1);
				}
			}
			if (cpface(+3) and fn == +3) {
				if (npedge(0,0,0)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db_p[2] += neigh_length[0]*nedge_e(0,0,0);
				}
				if (npedge(0,1,0)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db_p[2] -= neigh_length[0]*nedge_e(0,1,0);
				}
				if (npedge(1,0,0)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db_p[2] -= neigh_length[1]*nedge_e(1,0,0);
				}
				if (npedge(1,1,0)) {
					throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
					//face_db_p[2] += neigh_length[1]*nedge_e(1,1,0);
				}
			}
		}
		const std::array<double, 3> area{
			cell_length[1]*cell_length[2],
			cell_length[0]*cell_length[2],
			cell_length[0]*cell_length[1]
		};
		if (cfinfo(-1) >= 0 and cpface(-1)) Face_B_neg(*cell.data)[0] -= dt*face_db_n[0]/area[0];
		if (cfinfo(+1) >= 0 and cpface(+1)) Face_B_pos(*cell.data)[0] -= dt*face_db_p[0]/area[0];
		if (cfinfo(-2) >= 0 and cpface(-2)) Face_B_neg(*cell.data)[1] -= dt*face_db_n[1]/area[1];
		if (cfinfo(+2) >= 0 and cpface(+2)) Face_B_pos(*cell.data)[1] -= dt*face_db_p[1]/area[1];
		if (cfinfo(-3) >= 0 and cpface(-3)) Face_B_neg(*cell.data)[2] -= dt*face_db_n[2]/area[2];
		if (cfinfo(+3) >= 0 and cpface(+3)) Face_B_pos(*cell.data)[2] -= dt*face_db_p[2]/area[2];
	}
} catch (const std::exception& e) {
	throw std::runtime_error(__func__ + std::string(": ") + e.what());
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
	class Face_Magnetic_Field_Pos_Getter,
	class Face_Magnetic_Field_Neg_Getter,
	class Primary_Face_Getter,
	class Face_Type_Getter
> void update_B_consistency(
	const Cell_Iter& cells,
	const Mass_Density_Getter Mas,
	const Momentum_Density_Getter Mom,
	const Total_Energy_Density_Getter Nrj,
	const Volume_Magnetic_Field_Getter VMag,
	const Face_Magnetic_Field_Pos_Getter Face_B_pos,
	const Face_Magnetic_Field_Neg_Getter Face_B_neg,
	const Primary_Face_Getter PFace,
	const Face_Type_Getter FInfo,
	const double adiabatic_index,
	const double vacuum_permeability,
	const bool constant_thermal_pressure
) try {
	using std::runtime_error;
	using std::to_string;

	for (const auto& cell: cells) {
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

		auto
			&c_face_b_neg = Face_B_neg(*cell.data),
			&c_face_b_pos = Face_B_pos(*cell.data);

		const auto& cpface = PFace(*cell.data);
		const auto& cfinfo = FInfo(*cell.data);

		for (const size_t dim: {0, 1, 2})
		for (const int side: {-1, +1}) {
			if (not cpface(dim, side)) {
				if (side < 0) {
					c_face_b_neg[dim] = 0;
				} else {
					c_face_b_pos[dim] = 0;
				}
			} else if (cfinfo(dim, side) < 0) { // handle dont_solve face
				if (side < 0) {
					c_face_b_neg[dim] = 0;
				} else {
					c_face_b_pos[dim] = 0;
				}
				if (cpface(dim, -side) and cfinfo(dim, -side) >= 0) {
					if (side < 0) {
						c_face_b_neg[dim] = c_face_b_pos[dim];
					} else {
						c_face_b_pos[dim] = c_face_b_neg[dim];
					}
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

			const auto
				&n_face_b_neg = Face_B_neg(*neighbor.data),
				&n_face_b_pos = Face_B_pos(*neighbor.data);

			const auto& npface = PFace(*neighbor.data);
			const auto& nfinfo = FInfo(*neighbor.data);

			const size_t dim = std::abs(fn) - 1;
			if (not cpface(fn)) {
				if (not npface(-fn)) throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
				if (nfinfo(-fn) < 0) continue;

				if (fn < 0) {
					c_face_b_neg[dim] += n_face_b_pos[dim] * factor;
				} else {
					c_face_b_pos[dim] += n_face_b_neg[dim] * factor;
				}
			}
			// handle dont_solve face
			if (cpface(-fn) and cfinfo(-fn) < 0 and npface(-fn) and nfinfo(-fn) >= 0) {
				if (fn < 0) {
					c_face_b_pos[dim] += n_face_b_pos[dim] * factor;
				} else {
					c_face_b_neg[dim] += n_face_b_neg[dim] * factor;
				}
			}
		}

		VMag(*cell.data) = 0.5 * (c_face_b_neg + c_face_b_pos);

		if (constant_thermal_pressure) {
			const auto vel = (Mom(*cell.data)/Mas(*cell.data)).eval();
			Nrj(*cell.data) = pamhd::mhd::get_total_energy_density(
				Mas(*cell.data), vel, old_pressure, VMag(*cell.data),
				adiabatic_index, vacuum_permeability
			);
		}
	}
} catch (const std::exception& e) {
	throw std::runtime_error(__func__ + std::string(": ") + e.what());
}


}} // namespaces


#endif // ifndef PAMHD_MHD_SOLVE_STAGGERED_HPP
