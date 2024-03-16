/*
Handles boundary logic of MHD part of PAMHD.

Copyright 2015, 2016, 2017 Ilja Honkonen
Copyright 2023, 2024 Finnish Meteorological Institute
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

#ifndef PAMHD_MHD_BOUNDARIES_HPP
#define PAMHD_MHD_BOUNDARIES_HPP


#include "algorithm"
#include "cmath"
#include "iterator"
#include "limits"
#include "map"
#include "set"
#include "string"
#include "utility"
#include "vector"

#include "dccrg.hpp"

#include "grid/amr.hpp"
#include "mhd/common.hpp"
#include "mhd/variables.hpp"


namespace pamhd {
namespace mhd {


// as set_solver_info() below but for magnetic field only
// TODO: remove magnetic field from set_solver_info()
template<
	class Solver_Info,
	class Grid,
	class Boundaries,
	class Boundary_Geometries,
	class Solver_Info_Getter
> void set_solver_info_magnetic(
	Grid& grid,
	Boundaries& boundaries,
	const Boundary_Geometries& geometries,
	const Solver_Info_Getter& Sol_Info
) {
	using std::runtime_error;
	using std::to_string;

	boundaries.classify(grid, geometries, Sol_Info);

	for (const auto& cell: grid.cells) {
		Sol_Info(*cell.data) = 0;
	}

	// magnetic field
	constexpr pamhd::Magnetic_Field B{};
	for (const auto& cell: boundaries.get_value_boundary_cells(B)) {
		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
		}
		Sol_Info(*cell_data) |= Solver_Info::magnetic_field_bdy;
	}
	for (const auto& item: boundaries.get_copy_boundary_cells(B)) {
		auto* const cell_data = grid[item[0]];
		if (cell_data == nullptr) {
			throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
		}
		Sol_Info(*cell_data) |= Solver_Info::magnetic_field_bdy;
	}
	const std::set<uint64_t> dont_solve_mag(
		boundaries.get_dont_solve_cells(B).cbegin(),
		boundaries.get_dont_solve_cells(B).cend()
	);

	// don't solve cells in which no variable is solved
	for (auto& cell: dont_solve_mag) {
		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
		}
		Sol_Info(*cell_data) |= Solver_Info::dont_solve;
	}

	grid.update_copies_of_remote_neighbors();
}


/*!
Prepares boundary information needed for MHD solver about each simulation cell.

MPI transfer of Sol_Info variable must be switched on before calling this.
*/
template<
	class Solver_Info,
	class Grid,
	class Boundaries,
	class Boundary_Geometries,
	class Solver_Info_Getter
> void set_solver_info(
	Grid& grid,
	Boundaries& boundaries,
	const Boundary_Geometries& geometries,
	const Solver_Info_Getter& Sol_Info
) {
	using std::runtime_error;
	using std::to_string;

	boundaries.classify(grid, geometries, Sol_Info);

	for (const auto& cell: grid.cells) {
		Sol_Info(*cell.data) = 0;
	}

	// number density
	constexpr pamhd::mhd::Number_Density N{};
	for (const auto& cell: boundaries.get_value_boundary_cells(N)) {
		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
		}
		Sol_Info(*cell_data) |= Solver_Info::mass_density_bdy;
	}
	for (const auto& item: boundaries.get_copy_boundary_cells(N)) {
		auto* const cell_data = grid[item[0]];
		if (cell_data == nullptr) {
			throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
		}
		Sol_Info(*cell_data) |= Solver_Info::mass_density_bdy;
	}
	const std::set<uint64_t> dont_solve_mass(
		boundaries.get_dont_solve_cells(N).cbegin(),
		boundaries.get_dont_solve_cells(N).cend()
	);

	// velocity
	constexpr pamhd::mhd::Velocity V{};
	for (const auto& cell: boundaries.get_value_boundary_cells(V)) {
		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
		}
		Sol_Info(*cell_data) |= Solver_Info::velocity_bdy;
	}
	for (const auto& item: boundaries.get_copy_boundary_cells(V)) {
		auto* const cell_data = grid[item[0]];
		if (cell_data == nullptr) {
			throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
		}
		Sol_Info(*cell_data) |= Solver_Info::velocity_bdy;
	}
	const std::set<uint64_t> dont_solve_velocity(
		boundaries.get_dont_solve_cells(V).cbegin(),
		boundaries.get_dont_solve_cells(V).cend()
	);

	// pressure
	constexpr pamhd::mhd::Pressure P{};
	for (const auto& cell: boundaries.get_value_boundary_cells(P)) {
		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
		}
		Sol_Info(*cell_data) |= Solver_Info::pressure_bdy;
	}
	for (const auto& item: boundaries.get_copy_boundary_cells(P)) {
		auto* const cell_data = grid[item[0]];
		if (cell_data == nullptr) {
			throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
		}
		Sol_Info(*cell_data) |= Solver_Info::pressure_bdy;
	}
	const std::set<uint64_t> dont_solve_pressure(
		boundaries.get_dont_solve_cells(P).cbegin(),
		boundaries.get_dont_solve_cells(P).cend()
	);

	// magnetic field
	constexpr pamhd::Magnetic_Field B{};
	for (const auto& cell: boundaries.get_value_boundary_cells(B)) {
		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
		}
		Sol_Info(*cell_data) |= Solver_Info::magnetic_field_bdy;
	}
	for (const auto& item: boundaries.get_copy_boundary_cells(B)) {
		auto* const cell_data = grid[item[0]];
		if (cell_data == nullptr) {
			throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
		}
		Sol_Info(*cell_data) |= Solver_Info::magnetic_field_bdy;
	}
	const std::set<uint64_t> dont_solve_mag(
		boundaries.get_dont_solve_cells(B).cbegin(),
		boundaries.get_dont_solve_cells(B).cend()
	);

	if (
		not std::equal(
			dont_solve_mass.cbegin(),
			dont_solve_mass.cend(),
			dont_solve_velocity.cbegin()
		)
	) {
		throw std::invalid_argument(
			std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
			+ "Mass density and velocity dont_solves aren't equal."
		);
	}
	if (
		not std::equal(
			dont_solve_mass.cbegin(),
			dont_solve_mass.cend(),
			dont_solve_pressure.cbegin()
		)
	) {
		throw std::invalid_argument(
			std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
			+ "Mass density and pressure dont_solves aren't equal."
		);
	}
	if (
		not std::equal(
			dont_solve_mass.cbegin(),
			dont_solve_mass.cend(),
			dont_solve_mag.cbegin()
		)
	) {
		throw std::invalid_argument(
			std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
			+ "Mass density and magnetic field dont_solves aren't equal."
		);
	}

	// don't solve cells in which no variable is solved
	for (auto& cell: dont_solve_mass) {
		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
		}
		Sol_Info(*cell_data) |= Solver_Info::dont_solve;
	}

	grid.update_copies_of_remote_neighbors();
}


// as apply_magnetic_field_boundaries() below but for staggered solver.
template<
	class Grid,
	class Boundaries,
	class Boundary_Geometries,
	class Magnetic_Field_Getter,
	class Primary_Face_Getter,
	class Face_Info_Getter
> void apply_magnetic_field_boundaries_staggered(
	Grid& grid,
	Boundaries& boundaries,
	const Boundary_Geometries& bdy_geoms,
	const double t,
	const Magnetic_Field_Getter& Face_B,
	const Primary_Face_Getter& PFace,
	const Face_Info_Getter& FInfo
) try {
	using std::array;
	using std::asin;
	using std::atan2;
	using std::runtime_error;
	using std::sqrt;
	using std::to_string;

	constexpr pamhd::Magnetic_Field B{};
	for (
		size_t i = 0;
		i < boundaries.get_number_of_value_boundaries(B);
		i++
	) {
		auto& value_bdy = boundaries.get_value_boundary(B, i);
		const auto& geometry_id = value_bdy.get_geometry_id();
		const auto& cells = bdy_geoms.get_cells(geometry_id);
		for (const auto& cell: cells) {
			auto* const cell_data = grid[cell];
			if (cell_data == nullptr) {
				throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
			}

			const auto
				c = grid.geometry.get_center(cell),
				len = grid.geometry.get_length(cell);

			const auto get_B = [&value_bdy, &t](const array<double, 3>& c) -> pamhd::Magnetic_Field::data_type {
				const auto
					r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]),
					lat = asin(c[2] / c[0]),
					lon = atan2(c[1], c[0]);
				return value_bdy.get_data(t, c[0], c[1], c[2], r, lat, lon);
			};

			const auto& pface = PFace(*cell_data);
			const auto& finfo = FInfo(*cell_data);
			for (auto dim: {0, 1, 2})
			for (auto side: {-1, +1}) {
				if (pface(dim, side) and finfo(dim, side) == 0) {
					auto r = c;
					r[dim] += side * len[dim]/2;
					Face_B(*cell_data)(dim, side) = get_B(r)[dim];
				}
			}
		}
	}

	std::set<uint64_t> cp_bdy_cells;
	for (const auto& items: boundaries.get_copy_boundary_cells(B)) {
		if (items.size() == 0) throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
		cp_bdy_cells.insert(*items.begin());
	}

	for (const auto& cell: grid.local_cells()) {
		if (cp_bdy_cells.count(cell.id) == 0) continue;

		auto& c_face_b = Face_B(*cell.data);
		const auto& cpface = PFace(*cell.data);
		const auto& cfinfo = FInfo(*cell.data);

		pamhd::grid::Face_Type<int> nr_values{0, 0, 0, 0, 0, 0};
		for (auto dim: {0, 1, 2}) {
			for (auto side: {-1, +1}) {
				if (cpface(dim, side) and cfinfo(dim, side) < 1) c_face_b(dim, side) = 0;
			}

			if (cpface(dim, -1) and cpface(dim, +1)) {
				if (cfinfo(dim, -1) == 0 and cfinfo(dim, +1) == 0) throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
				if (cfinfo(dim, -1) == -1 and cfinfo(dim, +1) == -1) throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");

				if (cfinfo(dim, -1) == 0 and cfinfo(dim, +1) == 1) {
					c_face_b(dim, -1) += c_face_b(dim, +1);
					nr_values(dim, -1)++;
				}
				if (cfinfo(dim, -1) == 1 and cfinfo(dim, +1) == 0) {
					c_face_b(dim, +1) += c_face_b(dim, -1);
					nr_values(dim, +1)++;
				}
			}
		}

		// relative face neighbor sizes
		int nxs = 9, pxs = 9, nys = 9, pys = 9, nzs = 9, pzs = 9;
		for (const auto& neighbor: cell.neighbors_of) {
			switch (neighbor.face_neighbor) {
			case -1:
				nxs = neighbor.relative_size;
				break;
			case +1:
				pxs = neighbor.relative_size;
				break;
			case -2:
				nys = neighbor.relative_size;
				break;
			case +2:
				pys = neighbor.relative_size;
				break;
			case -3:
				nzs = neighbor.relative_size;
				break;
			case +3:
				pzs = neighbor.relative_size;
				break;
			default: break;
			}
		}

		for (const auto& neighbor: cell.neighbors_of) {
			auto& n_face_b = Face_B(*neighbor.data);

			const auto& fn = neighbor.face_neighbor;
			const auto& en = neighbor.edge_neighbor;
			if (fn == 0 and en[0] < 0) continue;

			const auto& npface = PFace(*neighbor.data);
			const auto& nfinfo = FInfo(*neighbor.data);

			if (fn != 0) {
				for (auto dir: {-3,-2,-1,+1,+2,+3}) {
					if (
						cpface(dir) and cfinfo(dir) == 0 and
						npface(dir) and nfinfo(dir) == 1
					) {
						c_face_b(dir) += n_face_b(dir);
						nr_values(dir)++;
					}
				}
			}

			if (en[0] == 0 and en[1] == -1 and en[2] == -1) {
				for (auto dir: {-2, -3}) {
					if (
						cpface(+dir) and cfinfo(+dir) == 0 and
						npface(-dir) and nfinfo(-dir) == 1
					) {
						c_face_b(+dir) += n_face_b(-dir);
						nr_values(+dir)++;
					}
				}
			}

			if (en[0] == 0 and en[1] == -1) {
				const int d3 = +1;
				if (en[2] == d3) {
					for (auto dir: {-2, +3}) {
						if (
							cpface(+dir) and cfinfo(+dir) == 0 and
							npface(-dir) and nfinfo(-dir) == 1
						) {
							c_face_b(+dir) += n_face_b(-dir);
							nr_values(+dir)++;
						}
					}
				}
				if (en[2] != d3 and nys < 0) {
					if (const int dir = +3;
						cpface(dir) and cfinfo(dir) == 0 and
						npface(dir) and nfinfo(dir) == 1
					) {
						c_face_b(dir) += n_face_b(dir);
						nr_values(dir)++;
					}
				}
			}

			if (en[0] == 0 and en[2] == -1) {
				const int d2 = +1;
				if (en[1] == d2) {
					for (auto dir: {+2, -3}) {
						if (
							cpface(+dir) and cfinfo(+dir) == 0 and
							npface(-dir) and nfinfo(-dir) == 1
						) {
							c_face_b(+dir) += n_face_b(-dir);
							nr_values(+dir)++;
						}
					}
				}
				if (en[1] != d2 and nzs < 0) {
					if (const int dir = +2;
						cpface(dir) and cfinfo(dir) == 0 and
						npface(dir) and nfinfo(dir) == 1
					) {
						c_face_b(dir) += n_face_b(dir);
						nr_values(dir)++;
					}
				}
			}

			if (en[0] == 0) {
				const int d2 = +1, d3 = +1;
				if (en[1] == d2 and en[2] == d3) {
					for (auto dir: {+2, -3}) {
						if (
							cpface(+dir) and cfinfo(+dir) == 0 and
							npface(-dir) and nfinfo(-dir) == 1
						) {
							c_face_b(+dir) += n_face_b(-dir);
							nr_values(+dir)++;
						}
					}
				}
				if (en[1] != d2 and en[2] == d3 and pzs < 0) {
					if (const int dir = +2;
						cpface(dir) and cfinfo(dir) == 0 and
						npface(dir) and nfinfo(dir) == 1
					) {
						c_face_b(dir) += n_face_b(dir);
						nr_values(dir)++;
					}
				}
				if (en[1] == d2 and en[2] != d3 and pys < 0) {
					if (const int dir = +3;
						cpface(dir) and cfinfo(dir) == 0 and
						npface(dir) and nfinfo(dir) == 1
					) {
						c_face_b(dir) += n_face_b(dir);
						nr_values(dir)++;
					}
				}
			}

			if (en[0] == 1 and en[1] == -1 and en[2] == -1) {
				for (auto dir: {-1, -3}) {
					if (
						cpface(+dir) and cfinfo(+dir) == 0 and
						npface(-dir) and nfinfo(-dir) == 1
					) {
						c_face_b(+dir) += n_face_b(-dir);
						nr_values(+dir)++;
					}
				}
			}

			if (en[0] == 1 and en[1] == -1) {
				const int d3 = +1;
				if (en[2] == d3) {
					for (auto dir: {-1, +3}) {
						if (
							cpface(+dir) and cfinfo(+dir) == 0 and
							npface(-dir) and nfinfo(-dir) == 1
						) {
							c_face_b(+dir) += n_face_b(-dir);
							nr_values(+dir)++;
						}
					}
				}
				if (en[2] != d3 and nxs < 0) {
					if (const int dir = +3;
						cpface(dir) and cfinfo(dir) == 0 and
						npface(dir) and nfinfo(dir) == 1
					) {
						c_face_b(dir) += n_face_b(dir);
						nr_values(dir)++;
					}
				}
			}

			if (en[0] == 1 and en[2] == -1) {
				const int d2 = +1;
				if (en[1] == d2) {
					for (auto dir: {+1, -3}) {
						if (
							cpface(+dir) and cfinfo(+dir) == 0 and
							npface(-dir) and nfinfo(-dir) == 1
						) {
							c_face_b(+dir) += n_face_b(-dir);
							nr_values(+dir)++;
						}
					}
				}
				if (en[1] != d2 and nzs < 0) {
					if (const int dir = +1;
						cpface(dir) and cfinfo(dir) == 0 and
						npface(dir) and nfinfo(dir) == 1
					) {
						c_face_b(dir) += n_face_b(dir);
						nr_values(dir)++;
					}
				}
			}

			if (en[0] == 1) {
				const int d2 = +1, d3 = +1;
				if (en[1] == d2 and en[2] == d3) {
					for (auto dir: {+1, +3}) {
						if (
							cpface(+dir) and cfinfo(+dir) == 0 and
							npface(-dir) and nfinfo(-dir) == 1
						) {
							c_face_b(+dir) += n_face_b(-dir);
							nr_values(+dir)++;
						}
					}
				}
				if (en[1] != d2 and en[2] == d3 and pzs < 0) {
					if (const int dir = +1;
						cpface(dir) and cfinfo(dir) == 0 and
						npface(dir) and nfinfo(dir) == 1
					) {
						c_face_b(dir) += n_face_b(dir);
						nr_values(dir)++;
					}
				}
				if (en[1] == d2 and en[2] != d3 and pxs < 0) {
					if (const int dir = +3;
						cpface(dir) and cfinfo(dir) == 0 and
						npface(dir) and nfinfo(dir) == 1
					) {
						c_face_b(dir) += n_face_b(dir);
						nr_values(dir)++;
					}
				}
			}

			if (en[0] == 2 and en[1] == -1 and en[2] == -1) {
				for (auto dir: {-1, -2}) {
					if (
						cpface(+dir) and cfinfo(+dir) == 0 and
						npface(-dir) and nfinfo(-dir) == 1
					) {
						c_face_b(+dir) += n_face_b(-dir);
						nr_values(+dir)++;
					}
				}
			}

			if (en[0] == 2 and en[1] == -1) {
				const int d3 = +1;
				if (en[2] == d3) {
					for (auto dir: {-1, +2}) {
						if (
							cpface(+dir) and cfinfo(+dir) == 0 and
							npface(-dir) and nfinfo(-dir) == 1
						) {
							c_face_b(+dir) += n_face_b(-dir);
							nr_values(+dir)++;
						}
					}
				}
				if (en[2] != d3 and nxs < 0) {
					if (const int dir = +2;
						cpface(dir) and cfinfo(dir) == 0 and
						npface(dir) and nfinfo(dir) == 1
					) {
						c_face_b(dir) += n_face_b(dir);
						nr_values(dir)++;
					}
				}
			}

			if (en[0] == 2 and en[2] == -1) {
				const int d2 = +1;
				if (en[1] == d2) {
					for (auto dir: {+1, -2}) {
						if (
							cpface(+dir) and cfinfo(+dir) == 0 and
							npface(-dir) and nfinfo(-dir) == 1
						) {
							c_face_b(+dir) += n_face_b(-dir);
							nr_values(+dir)++;
						}
					}
				}
				if (en[1] != d2 and nys < 0) {
					if (const int dir = +1;
						cpface(dir) and cfinfo(dir) == 0 and
						npface(dir) and nfinfo(dir) == 1
					) {
						c_face_b(dir) += n_face_b(dir);
						nr_values(dir)++;
					}
				}
			}

			if (en[0] == 2) {
				const int d2 = +1, d3 = +1;
				if (en[1] == d2 and en[2] == d3) {
					for (auto dir: {+1, +2}) {
						if (
							cpface(+dir) and cfinfo(+dir) == 0 and
							npface(-dir) and nfinfo(-dir) == 1
						) {
							c_face_b(+dir) += n_face_b(-dir);
							nr_values(+dir)++;
						}
					}
				}
				if (en[1] != d2 and en[2] == d3 and pys < 0) {
					if (const int dir = +1;
						cpface(dir) and cfinfo(dir) == 0 and
						npface(dir) and nfinfo(dir) == 1
					) {
						c_face_b(dir) += n_face_b(dir);
						nr_values(dir)++;
					}
				}
				if (en[1] == d2 and en[2] != d3 and pxs < 0) {
					if (const int dir = +2;
						cpface(dir) and cfinfo(dir) == 0 and
						npface(dir) and nfinfo(dir) == 1
					) {
						c_face_b(dir) += n_face_b(dir);
						nr_values(dir)++;
					}
				}
			}
		}

		for (auto dir: {-3,-2,-1,+1,+2,+3}) {
			if (cpface(dir) and cfinfo(dir) == 0) {
				if (nr_values(dir) == 0) {
					throw runtime_error(
						__FILE__ "(" + to_string(__LINE__) + "): "
						+ to_string(cell.id) + " " + to_string(dir));
				}
				c_face_b(dir) /= nr_values(dir);
			}
		}
	}
} catch (const std::exception& e) {
	throw std::runtime_error(__func__ + std::string(": ") + e.what());
}


// as apply_boundaries() below but for magnetic field only.
template<
	class Grid,
	class Boundaries,
	class Boundary_Geometries,
	class Magnetic_Field_Getter
> void apply_magnetic_field_boundaries(
	Grid& grid,
	Boundaries& boundaries,
	const Boundary_Geometries& bdy_geoms,
	const double simulation_time,
	const Magnetic_Field_Getter& Mag
) try {
	using std::runtime_error;
	using std::to_string;

	// magnetic field
	constexpr pamhd::Magnetic_Field B{};
	for (
		size_t i = 0;
		i < boundaries.get_number_of_value_boundaries(B);
		i++
	) {
		auto& value_bdy = boundaries.get_value_boundary(B, i);
		const auto& geometry_id = value_bdy.get_geometry_id();
		const auto& cells = bdy_geoms.get_cells(geometry_id);
		for (const auto& cell: cells) {
			const auto c = grid.geometry.get_center(cell);
			const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
			const auto
				lat = asin(c[2] / r),
				lon = atan2(c[1], c[0]);

			const auto magnetic_field = value_bdy.get_data(
				simulation_time,
				c[0], c[1], c[2],
				r, lat, lon
			);

			auto* const cell_data = grid[cell];
			if (cell_data == nullptr) {
				throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
			}

			Mag(*cell_data) = magnetic_field;
		}
	}
	for (const auto& item: boundaries.get_copy_boundary_cells(B)) {
		if (item.size() < 2) {
			throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
		}

		pamhd::Magnetic_Field::data_type source_value{0, 0, 0};
		for (size_t i = 1; i < item.size(); i++) {
			auto* source_data = grid[item[i]];
			if (source_data == nullptr) {
				throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
			}

			source_value += Mag(*source_data);
		}
		source_value /= item.size() - 1;

		auto *target_data = grid[item[0]];
		if (target_data == nullptr) {
			throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
		}

		Mag(*target_data) = source_value;
	}
} catch (const std::exception& e) {
	throw std::runtime_error(__func__ + std::string(": ") + e.what());
}


/*!
Applies boundaries of all simulation variables.

Value boundaries are applied to all cells
within that bundary's geometry. Value
boundaries are applied in order given in json
data, in case several overlap in geometry so
last one remains in effect.

Copy boundaries are applied after all value
boundaries. In case of more than one normal
neighbor their average is copied, vector
variables are processed by component.
*/
template<
	class Grid,
	class Boundaries,
	class Boundary_Geometries,
	class Mass_Getter,
	class Momentum_Getter,
	class Energy_Getter,
	class Magnetic_Field_Getter,
	class Cell_Info_Getter
> void apply_fluid_boundaries(
	Grid& grid,
	Boundaries& boundaries,
	const Boundary_Geometries& bdy_geoms,
	const double simulation_time,
	const Mass_Getter& Mas,
	const Momentum_Getter& Mom,
	const Energy_Getter& Nrj,
	const Magnetic_Field_Getter& Mag,
	const Cell_Info_Getter& CInfo,
	const double proton_mass,
	const double adiabatic_index,
	const double vacuum_permeability
) try {
	using std::runtime_error;
	using std::to_string;

	std::set<uint64_t> cp_bdy_cells;

	// number density
	constexpr pamhd::mhd::Number_Density N{};
	for (
		size_t i = 0;
		i < boundaries.get_number_of_value_boundaries(N);
		i++
	) {
		auto& value_bdy = boundaries.get_value_boundary(N, i);
		const auto& geometry_id = value_bdy.get_geometry_id();
		const auto& cells = bdy_geoms.get_cells(geometry_id);
		for (const auto& cell: cells) {
			const auto c = grid.geometry.get_center(cell);
			const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
			const auto
				lat = asin(c[2] / r),
				lon = atan2(c[1], c[0]);

			const auto mass_density
				= proton_mass
				* value_bdy.get_data(
					simulation_time,
					c[0], c[1], c[2],
					r, lat, lon
				);

			auto* const cell_data = grid[cell];
			if (cell_data == nullptr) {
				throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
			}

			Mas(*cell_data) = mass_density;
		}
	}

	cp_bdy_cells.clear();
	for (const auto& item: boundaries.get_copy_boundary_cells(N)) {
		if (item.size() < 1) throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
		cp_bdy_cells.insert(*item.begin());
	}

	for (const auto& cell: grid.local_cells()) {
		if (cp_bdy_cells.count(cell.id) == 0) continue;
		size_t nr_face_values = 0, nr_edge_values = 0;
		pamhd::mhd::Number_Density::data_type
			source_f_value = 0, source_e_value = 0;
		for (const auto& neighbor: cell.neighbors_of) {
			if (CInfo(*neighbor.data) < 1) continue;
			if (neighbor.face_neighbor != 0) {
				source_f_value += Mas(*neighbor.data);
				nr_face_values++;
			}
			if (neighbor.edge_neighbor[0] >= 0) {
				source_e_value += Mas(*neighbor.data);
				nr_edge_values++;
			}
		}
		if (nr_face_values > 0) {
			Mas(*cell.data) = source_f_value / nr_face_values;
		} else if (nr_edge_values > 0) {
			Mas(*cell.data) = source_e_value / nr_edge_values;
		}
	}

	// velocity
	constexpr pamhd::mhd::Velocity V{};
	for (
		size_t i = 0;
		i < boundaries.get_number_of_value_boundaries(V);
		i++
	) {
		auto& value_bdy = boundaries.get_value_boundary(V, i);
		const auto& geometry_id = value_bdy.get_geometry_id();
		const auto& cells = bdy_geoms.get_cells(geometry_id);
		for (const auto& cell: cells) {
			const auto c = grid.geometry.get_center(cell);
			const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
			const auto
				lat = asin(c[2] / r),
				lon = atan2(c[1], c[0]);

			const auto velocity = value_bdy.get_data(
				simulation_time,
				c[0], c[1], c[2],
				r, lat, lon
			);

			auto* const cell_data = grid[cell];
			if (cell_data == nullptr) {
				throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
			}

			Mom(*cell_data) = Mas(*cell_data) * velocity;
		}
	}

	cp_bdy_cells.clear();
	for (const auto& item: boundaries.get_copy_boundary_cells(V)) {
		if (item.size() < 1) throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
		cp_bdy_cells.insert(*item.begin());
	}

	for (const auto& cell: grid.local_cells()) {
		if (cp_bdy_cells.count(cell.id) == 0) continue;
		size_t nr_face_values = 0, nr_edge_values = 0;
		pamhd::mhd::Velocity::data_type
			source_f_value{0, 0, 0}, source_e_value{0, 0, 0};
		for (const auto& neighbor: cell.neighbors_of) {
			if (CInfo(*neighbor.data) < 1) continue;
			if (neighbor.face_neighbor != 0) {
				source_f_value += pamhd::mhd::get_velocity(
					Mom(*neighbor.data),
					Mas(*neighbor.data)
				);
				nr_face_values++;
			}
			if (neighbor.edge_neighbor[0] >= 0) {
				source_e_value += pamhd::mhd::get_velocity(
					Mom(*neighbor.data),
					Mas(*neighbor.data)
				);
				nr_edge_values++;
			}
		}
		if (nr_face_values > 0) {
			Mom(*cell.data) = Mas(*cell.data) * source_f_value / nr_face_values;
		} else if (nr_edge_values > 0) {
			Mom(*cell.data) = Mas(*cell.data) * source_e_value / nr_edge_values;
		}
	}

	// pressure
	constexpr pamhd::mhd::Pressure P{};
	for (
		size_t i = 0;
		i < boundaries.get_number_of_value_boundaries(P);
		i++
	) {
		auto& value_bdy = boundaries.get_value_boundary(P, i);
		const auto& geometry_id = value_bdy.get_geometry_id();
		const auto& cells = bdy_geoms.get_cells(geometry_id);
		for (const auto& cell: cells) {
			const auto c = grid.geometry.get_center(cell);
			const auto r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
			const auto
				lat = asin(c[2] / r),
				lon = atan2(c[1], c[0]);

			const auto pressure = value_bdy.get_data(
				simulation_time,
				c[0], c[1], c[2],
				r, lat, lon
			);

			auto* const cell_data = grid[cell];
			if (cell_data == nullptr) {
				throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
			}

			Nrj(*cell_data) = pamhd::mhd::get_total_energy_density(
				Mas(*cell_data),
				pamhd::mhd::get_velocity(
					Mom(*cell_data),
					Mas(*cell_data)
				),
				pressure,
				Mag(*cell_data),
				adiabatic_index,
				vacuum_permeability
			);
		}
	}

	cp_bdy_cells.clear();
	for (const auto& item: boundaries.get_copy_boundary_cells(P)) {
		if (item.size() < 1) throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
		cp_bdy_cells.insert(*item.begin());
	}

	for (const auto& cell: grid.local_cells()) {
		if (cp_bdy_cells.count(cell.id) == 0) continue;
		size_t nr_face_values = 0, nr_edge_values = 0;
		pamhd::mhd::Pressure::data_type
			source_f_value = 0, source_e_value = 0;
		for (const auto& neighbor: cell.neighbors_of) {
			if (CInfo(*neighbor.data) < 1) continue;
			if (neighbor.face_neighbor != 0) {
				source_f_value += pamhd::mhd::get_pressure(
					Mas(*neighbor.data),
					Mom(*neighbor.data),
					Nrj(*neighbor.data),
					Mag(*neighbor.data),
					adiabatic_index,
					vacuum_permeability
				);
				nr_face_values++;
			}
			if (neighbor.edge_neighbor[0] >= 0) {
				source_e_value += pamhd::mhd::get_pressure(
					Mas(*neighbor.data),
					Mom(*neighbor.data),
					Nrj(*neighbor.data),
					Mag(*neighbor.data),
					adiabatic_index,
					vacuum_permeability
				);
				nr_edge_values++;
			}
		}
		if (nr_face_values > 0) {
			Nrj(*cell.data) = pamhd::mhd::get_total_energy_density(
				Mas(*cell.data),
				pamhd::mhd::get_velocity(
					Mom(*cell.data),
					Mas(*cell.data)
				),
				source_f_value / nr_face_values,
				Mag(*cell.data),
				adiabatic_index,
				vacuum_permeability
			);
		} else if (nr_edge_values > 0) {
			Nrj(*cell.data) = pamhd::mhd::get_total_energy_density(
				Mas(*cell.data),
				pamhd::mhd::get_velocity(
					Mom(*cell.data),
					Mas(*cell.data)
				),
				source_e_value / nr_edge_values,
				Mag(*cell.data),
				adiabatic_index,
				vacuum_permeability
			);
		}
	}
} catch (const std::exception& e) {
	throw std::runtime_error(__func__ + std::string(": ") + e.what());
}


/*! Classifies primary faces of cells.

Types: normal cell/face == 1, boundary == 0, dont_solve == -1.

CInfo must be set for cells.

Update of FInfo variable between processes must be enabled.
*/
template <
	class Grid,
	class Primary_Face_Getter,
	class Cell_Info_Getter,
	class Face_Info_Getter
> void classify_faces(
	Grid& grid,
	const Primary_Face_Getter& PFace,
	const Cell_Info_Getter& CInfo,
	const Face_Info_Getter& FInfo
) try {
	const std::array<int, 6> all_dirs{-3,-2,-1,+1,+2,+3};
	// face(s) sharing edge(s) with normal cell(s) is normal
	for (const auto& cell: grid.local_cells()) {
		for (const int dir: all_dirs) {
			FInfo(*cell.data)(dir) = -99;
		}

		const auto& ccinfo = CInfo(*cell.data);
		const auto& cpface = PFace(*cell.data);
		if (ccinfo == 1) {
			for (const int dir: all_dirs) {
				if (PFace(*cell.data)(dir)) {
					FInfo(*cell.data)(dir) = 1;
				}
			}
			continue;
		}

		std::set<int> dirs; // faces to set as normal
		for (const auto& neighbor: cell.neighbors_of) {
			const auto& fn = neighbor.face_neighbor;
			if (
				fn != 0
				and cpface(fn)
				and CInfo(*neighbor.data) == 1
			) FInfo(*cell.data)(fn) = 1;
		}
	}

	grid.update_copies_of_remote_neighbors();

	/*
	Face(s) sharing edge(s) with normal face(s) or face(s)
	on opposite side of cell from normal face(s) are boundary
	*/
	for (const auto& cell: grid.local_cells()) {
		std::set<int> dirs; // faces to set as boundary

		const auto& cpface = PFace(*cell.data);
		const auto& cfinfo = FInfo(*cell.data);
		for (const int dir: all_dirs) {
			if (cpface(dir) and cfinfo(dir) == 1) {
				dirs.insert(all_dirs.cbegin(), all_dirs.cend());
				break;
			}
		}
		for (const auto& neighbor: cell.neighbors_of) {
			if (dirs.size() == all_dirs.size()) break;

			const auto& fn = neighbor.face_neighbor;
			const auto& en = neighbor.edge_neighbor;
			if (fn == 0 and en[0] < 0) continue;

			const auto& npface = PFace(*neighbor.data);
			const auto& nfinfo = FInfo(*neighbor.data);

			if (fn != 0) {
				if (npface(fn) and nfinfo(fn) == 1) dirs.insert(fn);
				if (npface(-fn) and nfinfo(-fn) == 1) {
					dirs.insert(all_dirs.cbegin(), all_dirs.cend());
					break;
				}
				for (const int dir: all_dirs) {
					if (dir == fn or dir == -fn) continue;
					if (npface(dir) and nfinfo(dir) == 1) dirs.insert(dir);
				}
			}

			int dir1 = 9, dir2 = 9;
			if (en[0] == 0) {
				if (en[1] == -1 and en[2] == -1) {
					if (dir1 = +2, dir2 = +3;
						   (npface(dir1) and nfinfo(dir1) == 1)
						or (npface(dir2) and nfinfo(dir2) == 1)
					) dirs.insert({-dir1, -dir2});
				}
				if (en[1] == -1 and en[2] == +1) {
					if (dir1 = +2, dir2 = -3;
						   (npface(dir1) and nfinfo(dir1) == 1)
						or (npface(dir2) and nfinfo(dir2) == 1)
					) dirs.insert({-dir1, -dir2});
				}
				if (en[1] == +1 and en[2] == -1) {
					if (dir1 = -2, dir2 = +3;
						   (npface(dir1) and nfinfo(dir1) == 1)
						or (npface(dir2) and nfinfo(dir2) == 1)
					) dirs.insert({-dir1, -dir2});
				}
				if (en[1] == +1 and en[2] == +1) {
					if (dir1 = -2, dir2 = -3;
						   (npface(dir1) and nfinfo(dir1) == 1)
						or (npface(dir2) and nfinfo(dir2) == 1)
					) dirs.insert({-dir1, -dir2});
				}
			}
			if (en[0] == 1) {
				if (en[1] == -1 and en[2] == -1) {
					if (dir1 = +1, dir2 = +3;
						   (npface(dir1) and nfinfo(dir1) == 1)
						or (npface(dir2) and nfinfo(dir2) == 1)
					) dirs.insert({-dir1, -dir2});
				}
				if (en[1] == -1 and en[2] == +1) {
					if (dir1 = +1, dir2 = -3;
						   (npface(dir1) and nfinfo(dir1) == 1)
						or (npface(dir2) and nfinfo(dir2) == 1)
					) dirs.insert({-dir1, -dir2});
				}
				if (en[1] == +1 and en[2] == -1) {
					if (dir1 = -1, dir2 = +3;
						   (npface(dir1) and nfinfo(dir1) == 1)
						or (npface(dir2) and nfinfo(dir2) == 1)
					) dirs.insert({-dir1, -dir2});
				}
				if (en[1] == +1 and en[2] == +1) {
					if (dir1 = -1, dir2 = -3;
						   (npface(dir1) and nfinfo(dir1) == 1)
						or (npface(dir2) and nfinfo(dir2) == 1)
					) dirs.insert({-dir1, -dir2});
				}
			}
			if (en[0] == 2) {
				if (en[1] == -1 and en[2] == -1) {
					if (dir1 = +1, dir2 = +2;
						   (npface(dir1) and nfinfo(dir1) == 1)
						or (npface(dir2) and nfinfo(dir2) == 1)
					) dirs.insert({-dir1, -dir2});
				}
				if (en[1] == -1 and en[2] == +1) {
					if (dir1 = +1, dir2 = -2;
						   (npface(dir1) and nfinfo(dir1) == 1)
						or (npface(dir2) and nfinfo(dir2) == 1)
					) dirs.insert({-dir1, -dir2});
				}
				if (en[1] == +1 and en[2] == -1) {
					if (dir1 = -1, dir2 = +2;
						   (npface(dir1) and nfinfo(dir1) == 1)
						or (npface(dir2) and nfinfo(dir2) == 1)
					) dirs.insert({-dir1, -dir2});
				}
				if (en[1] == +1 and en[2] == +1) {
					if (dir1 = -1, dir2 = -2;
						   (npface(dir1) and nfinfo(dir1) == 1)
						or (npface(dir2) and nfinfo(dir2) == 1)
					) dirs.insert({-dir1, -dir2});
				}
			}
		}

		for (int dir: dirs) {
			if (cpface(dir) and cfinfo(dir) == -99) FInfo(*cell.data)(dir) = 0;
		}
	}

	grid.update_copies_of_remote_neighbors();

	// face(s) sharing edge(s) with boundary face(s) is dont_solve
	for (const auto& cell: grid.local_cells()) {
		std::set<int> dirs;

		const auto& cpface = PFace(*cell.data);
		const auto& cfinfo = FInfo(*cell.data);
		for (const int dir: all_dirs) {
			if (cpface(dir) and cfinfo(dir) == 0) {
				dirs.insert(all_dirs.cbegin(), all_dirs.cend());
				break;
			}
		}

		for (const auto& neighbor: cell.neighbors_of) {
			if (dirs.size() == all_dirs.size()) break;

			const auto& fn = neighbor.face_neighbor;
			const auto& en = neighbor.edge_neighbor;
			if (fn == 0 and en[0] < 0) continue;

			const auto& npface = PFace(*neighbor.data);
			const auto& nfinfo = FInfo(*neighbor.data);

			if (fn != 0) {
				if (npface(fn) and nfinfo(fn) == 0) dirs.insert(fn);
				if (npface(-fn) and nfinfo(-fn) == 0) {
					dirs.insert(all_dirs.cbegin(), all_dirs.cend());
					break;
				}
				for (const int dir: all_dirs) {
					if (dir == fn or dir == -fn) continue;
					if (npface(dir) and nfinfo(dir) == 0) dirs.insert(dir);
				}
			}

			int dir1 = 9, dir2 = 9;
			if (en[0] == 0) {
				if (en[1] == -1 and en[2] == -1) {
					if (dir1 = +2, dir2 = +3;
						   (npface(dir1) and nfinfo(dir1) == 0)
						or (npface(dir2) and nfinfo(dir2) == 0)
					) dirs.insert({-dir1, -dir2});
				}
				if (en[1] == -1 and en[2] == +1) {
					if (dir1 = +2, dir2 = -3;
						   (npface(dir1) and nfinfo(dir1) == 0)
						or (npface(dir2) and nfinfo(dir2) == 0)
					) dirs.insert({-dir1, -dir2});
				}
				if (en[1] == +1 and en[2] == -1) {
					if (dir1 = -2, dir2 = +3;
						   (npface(dir1) and nfinfo(dir1) == 0)
						or (npface(dir2) and nfinfo(dir2) == 0)
					) dirs.insert({-dir1, -dir2});
				}
				if (en[1] == +1 and en[2] == +1) {
					if (dir1 = -2, dir2 = -3;
						   (npface(dir1) and nfinfo(dir1) == 0)
						or (npface(dir2) and nfinfo(dir2) == 0)
					) dirs.insert({-dir1, -dir2});
				}
			}
			if (en[0] == 1) {
				if (en[1] == -1 and en[2] == -1) {
					if (dir1 = +1, dir2 = +3;
						   (npface(dir1) and nfinfo(dir1) == 0)
						or (npface(dir2) and nfinfo(dir2) == 0)
					) dirs.insert({-dir1, -dir2});
				}
				if (en[1] == -1 and en[2] == +1) {
					if (dir1 = +1, dir2 = -3;
						   (npface(dir1) and nfinfo(dir1) == 0)
						or (npface(dir2) and nfinfo(dir2) == 0)
					) dirs.insert({-dir1, -dir2});
				}
				if (en[1] == +1 and en[2] == -1) {
					if (dir1 = -1, dir2 = +3;
						   (npface(dir1) and nfinfo(dir1) == 0)
						or (npface(dir2) and nfinfo(dir2) == 0)
					) dirs.insert({-dir1, -dir2});
				}
				if (en[1] == +1 and en[2] == +1) {
					if (dir1 = -1, dir2 = -3;
						   (npface(dir1) and nfinfo(dir1) == 0)
						or (npface(dir2) and nfinfo(dir2) == 0)
					) dirs.insert({-dir1, -dir2});
				}
			}
			if (en[0] == 2) {
				if (en[1] == -1 and en[2] == -1) {
					if (dir1 = +1, dir2 = +2;
						   (npface(dir1) and nfinfo(dir1) == 0)
						or (npface(dir2) and nfinfo(dir2) == 0)
					) dirs.insert({-dir1, -dir2});
				}
				if (en[1] == -1 and en[2] == +1) {
					if (dir1 = +1, dir2 = -2;
						   (npface(dir1) and nfinfo(dir1) == 0)
						or (npface(dir2) and nfinfo(dir2) == 0)
					) dirs.insert({-dir1, -dir2});
				}
				if (en[1] == +1 and en[2] == -1) {
					if (dir1 = -1, dir2 = +2;
						   (npface(dir1) and nfinfo(dir1) == 0)
						or (npface(dir2) and nfinfo(dir2) == 0)
					) dirs.insert({-dir1, -dir2});
				}
				if (en[1] == +1 and en[2] == +1) {
					if (dir1 = -1, dir2 = -2;
						   (npface(dir1) and nfinfo(dir1) == 0)
						or (npface(dir2) and nfinfo(dir2) == 0)
					) dirs.insert({-dir1, -dir2});
				}
			}
		}

		for (const int dir: dirs) {
			if (cpface(dir) and cfinfo(dir) == -99) FInfo(*cell.data)(dir) = -1;
		}
	}
} catch (const std::exception& e) {
	throw std::runtime_error(__func__ + std::string(": ") + e.what());
}


/*! Classifies primary edges of given cells.

Types: normal face/edge == 1, boundary == 0, dont_solve == -1.
Edge type is maximum of touching faces' types.

FInfo must be set for given cells and their edge neighbors'
faces before calling this function.
*/
template <
	class Cells,
	//class Grid,
	class Primary_Face_Getter,
	class Primary_Edge_Getter,
	class Face_Info_Getter,
	class Edge_Info_Getter
> void classify_edges(
	const Cells& cells,
	const Primary_Face_Getter& PFace,
	const Primary_Edge_Getter& PEdge,
	const Face_Info_Getter& FInfo,
	const Edge_Info_Getter& EInfo
) try {
	using std::max;

	for (const auto& cell: cells) {
		auto& ceinfo = EInfo(*cell.data);
		for (int d1: {0, 1, 2})
		for (int d2: {-1, +1})
		for (int d3: {-1, +1}) {
			ceinfo(d1,d2,d3) = -2;
		}

		const auto& cpface = PFace(*cell.data);
		const auto& cpedge = PEdge(*cell.data);
		const auto& cfinfo = FInfo(*cell.data);

		if (cpface(-1)) {
			for (int d1: {1, 2})
			for (int d3: {-1, +1}) {
				if (constexpr int d2 = -1; cpedge(d1,d2,d3)) {
					ceinfo(d1,d2,d3) = max(cfinfo(d2), ceinfo(d1,d2,d3));
				}
			}
		}

		size_t d1 = 9;
		int d2 = 9, d3 = 9;
		if (cpface(+1)) {
			{d1 = 1, d2 = +1;
			if (d3 = -1; cpedge(d1,d2,d3)) {
				ceinfo(d1,d2,d3) = max(cfinfo(+1), ceinfo(d1,d2,d3));
			}
			if (d3 = +1; cpedge(d1,d2,d3)) {
				ceinfo(d1,d2,d3) = max(cfinfo(+1), ceinfo(d1,d2,d3));
			}}

			{d1 = 2, d2 = +1;
			if (d3 = -1; cpedge(d1,d2,d3)) {
				ceinfo(d1,d2,d3) = max(cfinfo(+1), ceinfo(d1,d2,d3));
			}
			if (d3 = +1; cpedge(d1,d2,d3)) {
				ceinfo(d1,d2,d3) = max(cfinfo(+1), ceinfo(d1,d2,d3));
			}}
		}

		if (cpface(-2)) {
			{d1 = 0, d2 = -1;
			if (d3 = -1; cpedge(d1,d2,d3)) {
				ceinfo(d1,d2,d3) = max(cfinfo(-2), ceinfo(d1,d2,d3));
			}
			if (d3 = +1; cpedge(d1,d2,d3)) {
				ceinfo(d1,d2,d3) = max(cfinfo(-2), ceinfo(d1,d2,d3));
			}}

			{d1 = 2, d3 = -1;
			if (d2 = -1; cpedge(d1,d2,d3)) {
				ceinfo(d1,d2,d3) = max(cfinfo(-2), ceinfo(d1,d2,d3));
			}
			if (d2 = +1; cpedge(d1,d2,d3)) {
				ceinfo(d1,d2,d3) = max(cfinfo(-2), ceinfo(d1,d2,d3));
			}}
		}

		if (cpface(+2)) {
			{d1 = 0, d2 = +1;
			if (d3 = -1; cpedge(d1,d2,d3)) {
				ceinfo(d1,d2,d3) = max(cfinfo(+2), ceinfo(d1,d2,d3));
			}
			if (d3 = +1; cpedge(d1,d2,d3)) {
				ceinfo(d1,d2,d3) = max(cfinfo(+2), ceinfo(d1,d2,d3));
			}}

			{d1 = 2, d3 = +1;
			if (d2 = -1; cpedge(d1,d2,d3)) {
				ceinfo(d1,d2,d3) = max(cfinfo(+2), ceinfo(d1,d2,d3));
			}
			if (d2 = +1; cpedge(d1,d2,d3)) {
				ceinfo(d1,d2,d3) = max(cfinfo(+2), ceinfo(d1,d2,d3));
			}}
		}

		if (cpface(-3)) {
			{d1 = 0, d3 = -1;
			if (d2 = -1; cpedge(d1,d2,d3)) {
				ceinfo(d1,d2,d3) = max(cfinfo(-3), ceinfo(d1,d2,d3));
			}
			if (d2 = +1; cpedge(d1,d2,d3)) {
				ceinfo(d1,d2,d3) = max(cfinfo(-3), ceinfo(d1,d2,d3));
			}}

			{d1 = 1, d3 = -1;
			if (d2 = -1; cpedge(d1,d2,d3)) {
				ceinfo(d1,d2,d3) = max(cfinfo(-3), ceinfo(d1,d2,d3));
			}
			if (d2 = +1; cpedge(d1,d2,d3)) {
				ceinfo(d1,d2,d3) = max(cfinfo(-3), ceinfo(d1,d2,d3));
			}}
		}

		if (cpface(+3)) {
			{d1 = 0, d3 = +1;
			if (d2 = -1; cpedge(d1,d2,d3)) {
				ceinfo(d1,d2,d3) = max(cfinfo(+3), ceinfo(d1,d2,d3));
			}
			if (d2 = +1; cpedge(d1,d2,d3)) {
				ceinfo(d1,d2,d3) = max(cfinfo(+3), ceinfo(d1,d2,d3));
			}}

			{d1 = 1, d3 = +1;
			if (d2 = -1; cpedge(d1,d2,d3)) {
				ceinfo(d1,d2,d3) = max(cfinfo(+3), ceinfo(d1,d2,d3));
			}
			if (d2 = +1; cpedge(d1,d2,d3)) {
				ceinfo(d1,d2,d3) = max(cfinfo(+3), ceinfo(d1,d2,d3));
			}}
		}

		for (const auto& neighbor: cell.neighbors_of) {
			const auto& fn = neighbor.face_neighbor;
			const auto& en = neighbor.edge_neighbor;
			if (fn == 0 and en[0] < 0) {
				continue;
			}

			const auto& npface = PFace(*neighbor.data);
			const auto& nfinfo = FInfo(*neighbor.data);

			int d1 = 9, d2 = 9, d3 = 9;
			if (fn == -1) {
				d2 = -1;
				if (d1 = 1, d3 = -1;
					cpedge(d1,d2,d3) and npface(-3)// and neighbor.z == 0
				) {
					ceinfo(d1,d2,d3) = max(nfinfo(-3), ceinfo(d1,d2,d3));
				}

				if (d1 = 1, d3 = +1;
					cpedge(d1,d2,d3) and npface(+3)// and cleni == neighbor.z + nleni
				) {
					ceinfo(d1,d2,d3) = max(nfinfo(+3), ceinfo(d1,d2,d3));
				}

				if (d1 = 2, d3 = -1;
					cpedge(d1,d2,d3) and npface(-2)// and neighbor.y == 0
				) {
					ceinfo(d1,d2,d3) = max(nfinfo(-2), ceinfo(d1,d2,d3));
				}

				if (d1 = 2, d3 = +1;
					cpedge(d1,d2,d3) and npface(+2)// and cleni == neighbor.y + nleni
				) {
					ceinfo(d1,d2,d3) = max(nfinfo(+2), ceinfo(d1,d2,d3));
				}
			}

			if (fn == +1) {
				d2 = +1;
				if (d1 = 2, d3 = -1;
					cpedge(d1,d2,d3) and npface(-2)// and neighbor.y == 0
				) {
					ceinfo(d1,d2,d3) = max(nfinfo(-2), ceinfo(d1,d2,d3));
				}

				if (d1 = 2, d3 = +1;
					cpedge(d1,d2,d3) and npface(+2)// and cleni == neighbor.y + nleni
				) {
					ceinfo(d1,d2,d3) = max(nfinfo(+2), ceinfo(d1,d2,d3));
				}

				if (d1 = 1, d3 = -1;
					cpedge(d1,d2,d3) and npface(-3)// and neighbor.z == 0
				) {
					ceinfo(d1,d2,d3) = max(nfinfo(-3), ceinfo(d1,d2,d3));
				}

				if (d1 = 1, d3 = +1;
					cpedge(d1,d2,d3) and npface(+3)// and cleni == neighbor.z + nleni
				) {
					ceinfo(d1,d2,d3) = max(nfinfo(+3), ceinfo(d1,d2,d3));
				}
			}

			if (fn == -2) {
				if (d1 = 2, d2 = -1, d3 = -1;
					cpedge(d1,d2,d3) and npface(-1)// and neighbor.x == 0
				) {
					ceinfo(d1,d2,d3) = max(nfinfo(-1), ceinfo(d1,d2,d3));
				}

				if (d1 = 2, d2 = +1, d3 = -1;
					cpedge(d1,d2,d3) and npface(+1)// and cleni == neighbor.x + nleni
				) {
					ceinfo(d1,d2,d3) = max(nfinfo(+1), ceinfo(d1,d2,d3));
				}

				if (d1 = 0, d2 = -1, d3 = -1;
					cpedge(d1,d2,d3) and npface(-3)// and neighbor.z == 0
				) {
					ceinfo(d1,d2,d3) = max(nfinfo(-3), ceinfo(d1,d2,d3));
				}

				if (d1 = 0, d2 = -1, d3 = +1;
					cpedge(d1,d2,d3) and npface(+3)// and cleni == neighbor.z + nleni
				) {
					ceinfo(d1,d2,d3) = max(nfinfo(+3), ceinfo(d1,d2,d3));
				}
			}

			if (fn == +2) {
				if (d1 = 2, d2 = -1, d3 = +1;
					cpedge(d1,d2,d3) and npface(-1)// and neighbor.x == 0
				) {
					ceinfo(d1,d2,d3) = max(nfinfo(-1), ceinfo(d1,d2,d3));
				}

				if (d1 = 2, d2 = +1, d3 = +1;
					cpedge(d1,d2,d3) and npface(+1)// and cleni == neighbor.x + nleni
				) {
					ceinfo(d1,d2,d3) = max(nfinfo(+1), ceinfo(d1,d2,d3));
				}

				if (d1 = 0, d2 = +1, d3 = -1;
					cpedge(d1,d2,d3) and npface(-3)// and neighbor.z == 0
				) {
					ceinfo(d1,d2,d3) = max(nfinfo(-3), ceinfo(d1,d2,d3));
				}

				if (d1 = 0, d2 = +1, d3 = +1;
					cpedge(d1,d2,d3) and npface(+3)// and cleni == neighbor.z + nleni
				) {
					ceinfo(d1,d2,d3) = max(nfinfo(+3), ceinfo(d1,d2,d3));
				}
			}

			if (fn == -3) {
				d3 = -1;
				if (d1 = 1, d2 = -1;
					cpedge(d1,d2,d3) and npface(-1)// and neighbor.x == 0
				) {
					ceinfo(d1,d2,d3) = max(nfinfo(-1), ceinfo(d1,d2,d3));
				}

				if (d1 = 1, d2 = +1;
					cpedge(d1,d2,d3) and npface(+1) //and cleni == neighbor.x + nleni
				) {
					ceinfo(d1,d2,d3) = max(nfinfo(+1), ceinfo(d1,d2,d3));
				}

				if (d1 = 0, d2 = -1;
					cpedge(d1,d2,d3) and npface(-2) //and neighbor.y == 0
				) {
					ceinfo(d1,d2,d3) = max(nfinfo(-2), ceinfo(d1,d2,d3));
				}

				if (d1 = 0, d2 = +1;
					cpedge(d1,d2,d3) and npface(+2) //and cleni == neighbor.y + nleni
				) {
					ceinfo(d1,d2,d3) = max(nfinfo(+2), ceinfo(d1,d2,d3));
				}
			}

			if (fn == +3) {
				d3 = +1;
				if (d1 = 1, d2 = -1;
					cpedge(d1,d2,d3) and npface(-1) //and neighbor.x == 0
				) {
					ceinfo(d1,d2,d3) = max(nfinfo(-1), ceinfo(d1,d2,d3));
				}

				if (d1 = 1, d2 = +1;
					cpedge(d1,d2,d3) and npface(+1) //and cleni == neighbor.x + nleni
				) {
					ceinfo(d1,d2,d3) = max(nfinfo(+1), ceinfo(d1,d2,d3));
				}

				if (d1 = 0, d2 = -1;
					cpedge(d1,d2,d3) and npface(-2) //and neighbor.y == 0
				) {
					ceinfo(d1,d2,d3) = max(nfinfo(-2), ceinfo(d1,d2,d3));
				}

				if (d1 = 0, d2 = +1;
					cpedge(d1,d2,d3) and npface(+2) //and cleni == neighbor.y + nleni
				) {
					ceinfo(d1,d2,d3) = max(nfinfo(+2), ceinfo(d1,d2,d3));
				}
			}

			if (d1 = 0, d2 = -1, d3 = -1; cpedge(d1,d2,d3) and
				en[0] == d1 and en[1] == d2 and en[2] == d3
			) {
				if (npface(+2)) {
					ceinfo(d1,d2,d3) = max(nfinfo(+2), ceinfo(d1,d2,d3));
				}
				if (npface(+3)) {
					ceinfo(d1,d2,d3) = max(nfinfo(+3), ceinfo(d1,d2,d3));
				}
			}

			if (d1 = 0, d2 = -1, d3 = +1; cpedge(d1,d2,d3) and
				en[0] == d1 and en[1] == d2 and en[2] == d3
			) {
				if (npface(+2)) {
					ceinfo(d1,d2,d3) = max(nfinfo(+2), ceinfo(d1,d2,d3));
				}
				if (npface(-3)) {
					ceinfo(d1,d2,d3) = max(nfinfo(-3), ceinfo(d1,d2,d3));
				}
			}

			if (d1 = 0, d2 = +1, d3 = -1; cpedge(d1,d2,d3) and
				en[0] == d1 and en[1] == d2 and en[2] == d3
			) {
				if (npface(-2)) {
					ceinfo(d1,d2,d3) = max(nfinfo(-2), ceinfo(d1,d2,d3));
				}
				if (npface(+3)) {
					ceinfo(d1,d2,d3) = max(nfinfo(+3), ceinfo(d1,d2,d3));
				}
			}

			if (d1 = 0, d2 = +1, d3 = +1; cpedge(d1,d2,d3) and
				en[0] == d1 and en[1] == d2 and en[2] == d3
			) {
				if (npface(-2)) {
					ceinfo(d1,d2,d3) = max(nfinfo(-2), ceinfo(d1,d2,d3));
				}
				if (npface(-3)) {
					ceinfo(d1,d2,d3) = max(nfinfo(-3), ceinfo(d1,d2,d3));
				}
			}

			if (d1 = 1, d2 = -1, d3 = -1; cpedge(d1,d2,d3) and
				en[0] == d1 and en[1] == d2 and en[2] == d3
			) {
				if (npface(+1)) {
					ceinfo(d1,d2,d3) = max(nfinfo(+1), ceinfo(d1,d2,d3));
				}
				if (npface(+3)) {
					ceinfo(d1,d2,d3) = max(nfinfo(+3), ceinfo(d1,d2,d3));
				}
			}

			if (d1 = 1, d2 = -1, d3 = +1; cpedge(d1,d2,d3) and
				en[0] == d1 and en[1] == d2 and en[2] == d3
			) {
				if (npface(+1)) {
					ceinfo(d1,d2,d3) = max(nfinfo(+1), ceinfo(d1,d2,d3));
				}
				if (npface(-3)) {
					ceinfo(d1,d2,d3) = max(nfinfo(-3), ceinfo(d1,d2,d3));
				}
			}

			if (d1 = 1, d2 = +1, d3 = -1; cpedge(d1,d2,d3) and
				en[0] == d1 and en[1] == d2 and en[2] == d3
			) {
				if (npface(-1)) {
					ceinfo(d1,d2,d3) = max(nfinfo(-1), ceinfo(d1,d2,d3));
				}
				if (npface(+3)) {
					ceinfo(d1,d2,d3) = max(nfinfo(+3), ceinfo(d1,d2,d3));
				}
			}

			if (d1 = 1, d2 = +1, d3 = +1; cpedge(d1,d2,d3) and
				en[0] == d1 and en[1] == d2 and en[2] == d3
			) {
				if (npface(-1)) {
					ceinfo(d1,d2,d3) = max(nfinfo(-1), ceinfo(d1,d2,d3));
				}
				if (npface(-3)) {
					ceinfo(d1,d2,d3) = max(nfinfo(-3), ceinfo(d1,d2,d3));
				}
			}

			if (d1 = 2, d2 = -1, d3 = -1;
				cpedge(d1,d2,d3) and en[0] == d1 and en[1] == d2 and en[2] == d3
			) {
				if (npface(+1)) {
					ceinfo(d1,d2,d3) = max(nfinfo(+1), ceinfo(d1,d2,d3));
				}
				if (npface(+2)) {
					ceinfo(d1,d2,d3) = max(nfinfo(+2), ceinfo(d1,d2,d3));
				}
			}

			if (d1 = 2, d2 = -1, d3 = +1; cpedge(d1,d2,d3) and
				en[0] == d1 and en[1] == d2 and en[2] == d3
			) {
				if (npface(+1)) {
					ceinfo(d1,d2,d3) = max(nfinfo(+1), ceinfo(d1,d2,d3));
				}
				if (npface(-2)) {
					ceinfo(d1,d2,d3) = max(nfinfo(-2), ceinfo(d1,d2,d3));
				}
			}

			if (d1 = 2, d2 = +1, d3 = -1; cpedge(d1,d2,d3) and
				en[0] == d1 and en[1] == d2 and en[2] == d3
			) {
				if (npface(-1)) {
					ceinfo(d1,d2,d3) = max(nfinfo(-1), ceinfo(d1,d2,d3));
				}
				if (npface(+2)) {
					ceinfo(d1,d2,d3) = max(nfinfo(+2), ceinfo(d1,d2,d3));
				}
			}

			if (d1 = 2, d2 = +1, d3 = +1; cpedge(d1,d2,d3) and
				en[0] == d1 and en[1] == d2 and en[2] == d3
			) {
				if (npface(-1)) {
					ceinfo(d1,d2,d3) = max(nfinfo(-1), ceinfo(d1,d2,d3));
				}
				if (npface(-2)) {
					ceinfo(d1,d2,d3) = max(nfinfo(-2), ceinfo(d1,d2,d3));
				}
			}
		}
	}
} catch (const std::exception& e) {
	throw std::runtime_error(__func__ + std::string(": ") + e.what());
}


//! Makes sure boundary cells and their neighbors dont' change size
template <
	class Grid,
	class Cell_Info_Getter
> void enforce_boundary_cell_sizes(
	Grid& grid,
	const Cell_Info_Getter& CInfo
) try {
	for (const auto& cell: grid.local_cells()) {
		const auto& ccinfo = CInfo(*cell.data);
		if (ccinfo < 1) {
			grid.dont_refine(cell.id);
			grid.dont_unrefine(cell.id);
			continue;
		}
		for (const auto& neighbor: cell.neighbors_of) {
			const auto& ncinfo = CInfo(*neighbor.data);
			if (ncinfo < 1) {
				grid.dont_refine(cell.id);
				grid.dont_unrefine(cell.id);
				break;
			}
		}
	}
} catch (const std::exception& e) {
	throw std::runtime_error(__func__ + std::string(": ") + e.what());
}

}} // namespaces


#endif // ifndef PAMHD_MHD_BOUNDARIES_HPP
