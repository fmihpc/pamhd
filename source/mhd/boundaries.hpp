/*
Handles boundary logic of MHD part of PAMHD.

Copyright 2015, 2016, 2017 Ilja Honkonen
Copyright 2023, 2024, 2025 Finnish Meteorological Institute
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
#include "mhd/solve.hpp"
#include "mhd/variables.hpp"


namespace pamhd {
namespace mhd {


/*!
Prepares boundary information needed for solver(s)
about each simulation cell.

MPI transfer of variable returned by SInfo must be
switched on before calling this.

1 means normal read-write cell
0 means boundary read-only cell
-1 means dont_solve cell that's not to be read nor written
*/
template<
	class Grid,
	class Boundaries,
	class Boundary_Geometries,
	class Solver_Info_Getter
> void set_solver_info(
	Grid& grid,
	Boundaries& boundaries,
	const Boundary_Geometries& geometries,
	const Solver_Info_Getter& SInfo
) try {
	using std::runtime_error;
	using std::to_string;

	using Cell = Grid::cell_data_type;
	Cell::set_transfer_all(true, SInfo.type());

	boundaries.classify(grid, geometries, SInfo);

	for (const auto& cell: grid.local_cells()) {
		SInfo.data(*cell.data) = 1;
	}

	// number density
	constexpr pamhd::mhd::Number_Density N{};
	for (const auto& cell: boundaries.get_value_boundary_cells(N)) {
		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
		}
		SInfo.data(*cell_data) = 0;
	}
	for (const auto& item: boundaries.get_copy_boundary_cells(N)) {
		auto* const cell_data = grid[item[0]];
		if (cell_data == nullptr) {
			throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
		}
		SInfo.data(*cell_data) = 0;
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
		SInfo.data(*cell_data) = 0;
	}
	for (const auto& item: boundaries.get_copy_boundary_cells(V)) {
		auto* const cell_data = grid[item[0]];
		if (cell_data == nullptr) {
			throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
		}
		SInfo.data(*cell_data) = 0;
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
		SInfo.data(*cell_data) = 0;
	}
	for (const auto& item: boundaries.get_copy_boundary_cells(P)) {
		auto* const cell_data = grid[item[0]];
		if (cell_data == nullptr) {
			throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
		}
		SInfo.data(*cell_data) = 0;
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
		SInfo.data(*cell_data) = 0;
	}
	for (const auto& item: boundaries.get_copy_boundary_cells(B)) {
		auto* const cell_data = grid[item[0]];
		if (cell_data == nullptr) {
			throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
		}
		SInfo.data(*cell_data) = 0;
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

	// don't solve cells
	for (auto& cell: dont_solve_mass) {
		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
		}
		SInfo.data(*cell_data) = -1;
	}

	// also don't solve cells too far away from local
	for (const auto& cell: grid.remote_cells()) {
		bool solve = false;
		for (const auto& neighbor: cell.neighbors_to) {
			if (neighbor.is_local) {
				solve = true;
				break;
			}
		}
		if (not solve) SInfo.data(*cell.data) = -1;
	}

	grid.update_copies_of_remote_neighbors();
	Cell::set_transfer_all(false, SInfo.type());

} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


// as apply_magnetic_field_boundaries() below but for staggered solver.
template<
	class Grid,
	class Boundaries,
	class Boundary_Geometries,
	class Magnetic_Field_Getter,
	class Face_Info_Getter
> void apply_magnetic_field_boundaries_staggered(
	Grid& grid,
	Boundaries& boundaries,
	const Boundary_Geometries& bdy_geoms,
	const double t,
	const Magnetic_Field_Getter& Face_B,
	const Face_Info_Getter& FInfo
) try {
	using std::array;
	using std::asin;
	using std::atan2;
	using std::runtime_error;
	using std::sqrt;
	using std::to_string;

	using Cell = Grid::cell_data_type;
	Cell::set_transfer_all(true, Face_B.type(), FInfo.type());
	grid.update_copies_of_remote_neighbors();
	Cell::set_transfer_all(false, Face_B.type(), FInfo.type());

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

			const auto& finfo = FInfo.data(*cell_data);
			for (auto dim: {0, 1, 2})
			for (auto side: {-1, +1}) {
				if (finfo(dim, side) == 0) {
					auto r = c;
					r[dim] += side * len[dim]/2;
					Face_B.data(*cell_data)(dim, side) = get_B(r)[dim];
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

		auto& c_face_b = Face_B.data(*cell.data);
		const auto& cfinfo = FInfo.data(*cell.data);

		pamhd::Face_Type<int> nr_values{0, 0, 0, 0, 0, 0};
		for (auto dim: {0, 1, 2}) {
			for (auto side: {-1, +1}) {
				if (cfinfo(dim, side) < 1) c_face_b(dim, side) = 0;
			}

			if (cfinfo(dim, -1) == 0 and cfinfo(dim, +1) == 0) throw runtime_error(
				__FILE__ "(" + to_string(__LINE__) + ") "
				+ to_string(dim) + " " + to_string(cell.id)
			);
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
			auto& n_face_b = Face_B.data(*neighbor.data);

			const auto& fn = neighbor.face_neighbor;
			const auto& en = neighbor.edge_neighbor;
			if (fn == 0 and en[0] < 0) continue;

			const auto& nfinfo = FInfo.data(*neighbor.data);

			if (fn != 0) {
				for (auto dir: {-3,-2,-1,+1,+2,+3}) {
					if (cfinfo(dir) == 0 and nfinfo(dir) == 1) {
						c_face_b(dir) += n_face_b(dir);
						nr_values(dir)++;
					}
				}
			}

			if (en[0] == 0 and en[1] == -1 and en[2] == -1) {
				for (auto dir: {-2, -3}) {
					if (cfinfo(+dir) == 0 and nfinfo(-dir) == 1) {
						c_face_b(+dir) += n_face_b(-dir);
						nr_values(+dir)++;
					}
				}
			}

			if (en[0] == 0 and en[1] == -1) {
				const int d3 = +1;
				if (en[2] == d3) {
					for (auto dir: {-2, +3}) {
						if (cfinfo(+dir) == 0 and nfinfo(-dir) == 1) {
							c_face_b(+dir) += n_face_b(-dir);
							nr_values(+dir)++;
						}
					}
				}
				if (en[2] != d3 and nys < 0) {
					if (const int dir = +3;
						cfinfo(dir) == 0 and nfinfo(dir) == 1
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
						if (cfinfo(+dir) == 0 and nfinfo(-dir) == 1) {
							c_face_b(+dir) += n_face_b(-dir);
							nr_values(+dir)++;
						}
					}
				}
				if (en[1] != d2 and nzs < 0) {
					if (const int dir = +2;
						cfinfo(dir) == 0 and nfinfo(dir) == 1
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
						if (cfinfo(+dir) == 0 and nfinfo(-dir) == 1) {
							c_face_b(+dir) += n_face_b(-dir);
							nr_values(+dir)++;
						}
					}
				}
				if (en[1] != d2 and en[2] == d3 and pzs < 0) {
					if (const int dir = +2;
						cfinfo(dir) == 0 and nfinfo(dir) == 1
					) {
						c_face_b(dir) += n_face_b(dir);
						nr_values(dir)++;
					}
				}
				if (en[1] == d2 and en[2] != d3 and pys < 0) {
					if (const int dir = +3;
						cfinfo(dir) == 0 and nfinfo(dir) == 1
					) {
						c_face_b(dir) += n_face_b(dir);
						nr_values(dir)++;
					}
				}
			}

			if (en[0] == 1 and en[1] == -1 and en[2] == -1) {
				for (auto dir: {-1, -3}) {
					if (cfinfo(+dir) == 0 and nfinfo(-dir) == 1) {
						c_face_b(+dir) += n_face_b(-dir);
						nr_values(+dir)++;
					}
				}
			}

			if (en[0] == 1 and en[1] == -1) {
				const int d3 = +1;
				if (en[2] == d3) {
					for (auto dir: {-1, +3}) {
						if (cfinfo(+dir) == 0 and nfinfo(-dir) == 1) {
							c_face_b(+dir) += n_face_b(-dir);
							nr_values(+dir)++;
						}
					}
				}
				if (en[2] != d3 and nxs < 0) {
					if (const int dir = +3;
						cfinfo(dir) == 0 and nfinfo(dir) == 1
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
						if (cfinfo(+dir) == 0 and nfinfo(-dir) == 1) {
							c_face_b(+dir) += n_face_b(-dir);
							nr_values(+dir)++;
						}
					}
				}
				if (en[1] != d2 and nzs < 0) {
					if (const int dir = +1;
						cfinfo(dir) == 0 and nfinfo(dir) == 1
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
						if (cfinfo(+dir) == 0 and nfinfo(-dir) == 1) {
							c_face_b(+dir) += n_face_b(-dir);
							nr_values(+dir)++;
						}
					}
				}
				if (en[1] != d2 and en[2] == d3 and pzs < 0) {
					if (const int dir = +1;
						cfinfo(dir) == 0 and nfinfo(dir) == 1
					) {
						c_face_b(dir) += n_face_b(dir);
						nr_values(dir)++;
					}
				}
				if (en[1] == d2 and en[2] != d3 and pxs < 0) {
					if (const int dir = +3;
						cfinfo(dir) == 0 and nfinfo(dir) == 1
					) {
						c_face_b(dir) += n_face_b(dir);
						nr_values(dir)++;
					}
				}
			}

			if (en[0] == 2 and en[1] == -1 and en[2] == -1) {
				for (auto dir: {-1, -2}) {
					if (cfinfo(+dir) == 0 and nfinfo(-dir) == 1) {
						c_face_b(+dir) += n_face_b(-dir);
						nr_values(+dir)++;
					}
				}
			}

			if (en[0] == 2 and en[1] == -1) {
				const int d3 = +1;
				if (en[2] == d3) {
					for (auto dir: {-1, +2}) {
						if (cfinfo(+dir) == 0 and nfinfo(-dir) == 1) {
							c_face_b(+dir) += n_face_b(-dir);
							nr_values(+dir)++;
						}
					}
				}
				if (en[2] != d3 and nxs < 0) {
					if (const int dir = +2;
						cfinfo(dir) == 0 and nfinfo(dir) == 1
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
						if (cfinfo(+dir) == 0 and nfinfo(-dir) == 1) {
							c_face_b(+dir) += n_face_b(-dir);
							nr_values(+dir)++;
						}
					}
				}
				if (en[1] != d2 and nys < 0) {
					if (const int dir = +1;
						cfinfo(dir) == 0 and nfinfo(dir) == 1
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
						if (cfinfo(+dir) == 0 and nfinfo(-dir) == 1) {
							c_face_b(+dir) += n_face_b(-dir);
							nr_values(+dir)++;
						}
					}
				}
				if (en[1] != d2 and en[2] == d3 and pys < 0) {
					if (const int dir = +1;
						cfinfo(dir) == 0 and nfinfo(dir) == 1
					) {
						c_face_b(dir) += n_face_b(dir);
						nr_values(dir)++;
					}
				}
				if (en[1] == d2 and en[2] != d3 and pxs < 0) {
					if (const int dir = +2;
						cfinfo(dir) == 0 and nfinfo(dir) == 1
					) {
						c_face_b(dir) += n_face_b(dir);
						nr_values(dir)++;
					}
				}
			}
		}

		for (auto dir: {-3,-2,-1,+1,+2,+3}) {
			if (cfinfo(dir) == 0) {
				if (nr_values(dir) == 0) {
					throw runtime_error(
						__FILE__ "(" + to_string(__LINE__) + "): "
						+ to_string(cell.id) + " " + to_string(dir));
				}
				c_face_b(dir) /= nr_values(dir);
			}
		}
	}
	Cell::set_transfer_all(true, Face_B.type());
	grid.update_copies_of_remote_neighbors();
	Cell::set_transfer_all(false, Face_B.type());

} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


// as apply_boundaries() below but for magnetic field only.
template<
	class Grid,
	class Boundaries,
	class Boundary_Geometries,
	class Volume_Magnetic_Field_Getter
> void apply_magnetic_field_boundaries(
	Grid& grid,
	Boundaries& boundaries,
	const Boundary_Geometries& bdy_geoms,
	const double& simulation_time,
	const Volume_Magnetic_Field_Getter& Vol_B
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

			Vol_B.data(*cell_data) = magnetic_field;
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

			source_value += Vol_B.data(*source_data);
		}
		source_value /= item.size() - 1;

		auto *target_data = grid[item[0]];
		if (target_data == nullptr) {
			throw runtime_error(__FILE__ "(" + to_string(__LINE__) + ")");
		}

		Vol_B.data(*target_data) = source_value;
	}
} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
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
	class Volume_Magnetic_Field_Getter,
	class Cell_Info_Getter
> void apply_fluid_boundaries(
	Grid& grid,
	Boundaries& boundaries,
	const Boundary_Geometries& bdy_geoms,
	const double simulation_time,
	const Mass_Getter& Mas,
	const Momentum_Getter& Mom,
	const Energy_Getter& Nrj,
	const Volume_Magnetic_Field_Getter& Vol_B,
	const Cell_Info_Getter& CInfo,
	const double& proton_mass,
	const double& adiabatic_index,
	const double& vacuum_permeability
) try {
	using std::runtime_error;
	using std::to_string;

	using Cell = Grid::cell_data_type;

	bool update_copies = false;
	if (Nrj.type().is_stale) {
		update_copies = true;
		// these become stale again below
		Cell::set_transfer_all(true,
			Mas.type(), Mom.type(), Nrj.type(), Vol_B.type());
	}
	if (CInfo.type().is_stale) {
		update_copies = true;
		Cell::set_transfer_all(true, CInfo.type());
		CInfo.type().is_stale = false;
	}
	if (update_copies) {
		grid.update_copies_of_remote_neighbors();
	}
	Cell::set_transfer_all(false,
		Mas.type(), Mom.type(), Nrj.type(),
		Vol_B.type(), CInfo.type());

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

			Mas.data(*cell_data) = mass_density;
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
			if (CInfo.data(*neighbor.data) < 1) continue;
			if (neighbor.face_neighbor != 0) {
				source_f_value += Mas.data(*neighbor.data);
				nr_face_values++;
			}
			if (neighbor.edge_neighbor[0] >= 0) {
				source_e_value += Mas.data(*neighbor.data);
				nr_edge_values++;
			}
		}
		if (nr_face_values > 0) {
			Mas.data(*cell.data) = source_f_value / nr_face_values;
		} else if (nr_edge_values > 0) {
			Mas.data(*cell.data) = source_e_value / nr_edge_values;
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

			Mom.data(*cell_data) = pamhd::mul(Mas.data(*cell_data), velocity);
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
			if (CInfo.data(*neighbor.data) < 1) continue;
			if (neighbor.face_neighbor != 0) {
				source_f_value = pamhd::add(source_f_value, pamhd::mhd::get_velocity(
					Mom.data(*neighbor.data),
					Mas.data(*neighbor.data)
				));
				nr_face_values++;
			}
			if (neighbor.edge_neighbor[0] >= 0) {
				source_e_value = pamhd::add(source_e_value, pamhd::mhd::get_velocity(
					Mom.data(*neighbor.data),
					Mas.data(*neighbor.data)
				));
				nr_edge_values++;
			}
		}
		if (nr_face_values > 0) {
			Mom.data(*cell.data) = pamhd::mul(
				pamhd::mul(Mas.data(*cell.data), source_f_value),
				1 / nr_face_values);
		} else if (nr_edge_values > 0) {
			Mom.data(*cell.data) = pamhd::mul(
				pamhd::mul(Mas.data(*cell.data), source_e_value),
				1 / nr_edge_values);
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

			Nrj.data(*cell_data) = pamhd::mhd::get_total_energy_density(
				Mas.data(*cell_data),
				pamhd::mhd::get_velocity(
					Mom.data(*cell_data),
					Mas.data(*cell_data)
				),
				pressure,
				Vol_B.data(*cell_data),
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
			if (CInfo.data(*neighbor.data) < 1) continue;
			if (neighbor.face_neighbor != 0) {
				source_f_value += pamhd::mhd::get_pressure(
					Mas.data(*neighbor.data),
					Mom.data(*neighbor.data),
					Nrj.data(*neighbor.data),
					Vol_B.data(*neighbor.data),
					adiabatic_index,
					vacuum_permeability
				);
				nr_face_values++;
			}
			if (neighbor.edge_neighbor[0] >= 0) {
				source_e_value += pamhd::mhd::get_pressure(
					Mas.data(*neighbor.data),
					Mom.data(*neighbor.data),
					Nrj.data(*neighbor.data),
					Vol_B.data(*neighbor.data),
					adiabatic_index,
					vacuum_permeability
				);
				nr_edge_values++;
			}
		}
		if (nr_face_values > 0) {
			Nrj.data(*cell.data) = pamhd::mhd::get_total_energy_density(
				Mas.data(*cell.data),
				pamhd::mhd::get_velocity(
					Mom.data(*cell.data),
					Mas.data(*cell.data)
				),
				source_f_value / nr_face_values,
				Vol_B.data(*cell.data),
				adiabatic_index,
				vacuum_permeability
			);
		} else if (nr_edge_values > 0) {
			Nrj.data(*cell.data) = pamhd::mhd::get_total_energy_density(
				Mas.data(*cell.data),
				pamhd::mhd::get_velocity(
					Mom.data(*cell.data),
					Mas.data(*cell.data)
				),
				source_e_value / nr_edge_values,
				Vol_B.data(*cell.data),
				adiabatic_index,
				vacuum_permeability
			);
		}
	}

} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


/*! Classifies cell faces.

Types: normal cell/face == 1, boundary == 0, dont_solve == -1.

CInfo must be set for cells.
*/
template <
	class Grid,
	class Cell_Info_Getter,
	class Face_Info_Getter
> void classify_faces(
	Grid& grid,
	const Cell_Info_Getter& CInfo,
	const Face_Info_Getter& FInfo
) try {
	using std::max;
	using std::min;
	using std::runtime_error;
	using std::to_string;

	using Cell = Grid::cell_data_type;
	Cell::set_transfer_all(true, CInfo.type(), FInfo.type());
	grid.update_copies_of_remote_neighbors();

	const std::array<int, 6> all_dirs{-3,-2,-1,+1,+2,+3};

	// face(s) sharing edge(s) with normal cell(s) is normal
	for (const auto& cell: grid.local_cells()) {
		auto& cfinfo = FInfo.data(*cell.data);

		for (auto dir: all_dirs) {
			cfinfo(dir) = -99;
		}

		const auto& ccinfo = CInfo.data(*cell.data);
		if (ccinfo == 1) {
			for (auto dir: all_dirs) {
				cfinfo(dir) = 1;
			}
			continue;
		}

		for (const auto& neighbor: cell.neighbors_of) {
			const auto& ncinfo = CInfo.data(*neighbor.data);
			if (ncinfo < 1) continue;

			const auto& fn = neighbor.face_neighbor;
			if (fn != 0) {
				for (auto dir: all_dirs) {
					if (dir == -fn) continue;
					cfinfo(dir) = 1;
				}
			}
			const auto& en = neighbor.edge_neighbor;
			if (en[0] >= 0) {
				const int
					dir1 = 1 + std::min((en[0] + 1) % 3, (en[0] + 2) % 3),
					dir2 = 1 + std::max((en[0] + 1) % 3, (en[0] + 2) % 3);
				if (en[1] < 0) {
					cfinfo(-dir1) = 1;
				} else {
					cfinfo(+dir1) = 1;
				}
				if (en[2] < 0) {
					cfinfo(-dir2) = 1;
				} else {
					cfinfo(+dir2) = 1;
				}
			}
		}
	}

	grid.update_copies_of_remote_neighbors();

	for (const auto& cell: grid.local_cells()) {
		const auto& cfinfo = FInfo.data(*cell.data);
		for (const auto& neighbor: cell.neighbors_of) {
			const auto& fn = neighbor.face_neighbor;
			const auto& en = neighbor.edge_neighbor;
			if (fn != 0 and en[0] < 0) continue;
			const auto& nfinfo = FInfo.data(*neighbor.data);
			if (fn != 0) {
				if (cfinfo(fn) != nfinfo(-fn)) {
					throw runtime_error(
						__FILE__ "(" + to_string(__LINE__) + "): "
						+ to_string(cell.id) + " "
						+ to_string(neighbor.id) + " "
						+ to_string(fn) + " "
						+ to_string(cfinfo(fn)) + " "
						+ to_string(nfinfo(-fn))
					);
				}
			}
		}
	}

	/*
	Face(s) sharing edge(s) with normal face(s) or face(s)
	on opposite side of cell from normal face(s) are boundary
	*/
	for (const auto& cell: grid.local_cells()) {
		// if even one cell face is normal other faces are boundary
		bool has_normal_face = false;
		auto& cfinfo = FInfo.data(*cell.data);
		for (auto dir: all_dirs) {
			if (cfinfo(dir) == 1) {
				has_normal_face = true;
				break;
			}
		}
		if (not has_normal_face) continue;

		for (auto dir: all_dirs) {
			if (cfinfo(dir) < 0) cfinfo(dir) = 0;
		}
	}

	grid.update_copies_of_remote_neighbors();

	// as above but for neighboring cells
	for (const auto& cell: grid.local_cells()) {
		auto& cfinfo = FInfo.data(*cell.data);
		for (const auto& neighbor: cell.neighbors_of) {
			const auto& nfinfo = FInfo.data(*neighbor.data);
			const auto& fn = neighbor.face_neighbor;
			if (fn != 0) {
				for (auto dir: all_dirs) {
					if (nfinfo(dir) < 1) continue;
					if (dir == fn and cfinfo(dir) < 0) {
						cfinfo(dir) = 0;
						continue;
					}
					if (dir == -fn and cfinfo(fn) < 1) throw runtime_error(
						__FILE__ "(" + to_string(__LINE__) + "): "
						+ to_string(cell.id) + " "
						+ to_string(neighbor.id) + " "
						+ to_string(dir) + " "
						+ to_string(cfinfo(fn)) + " "
						+ to_string(nfinfo(-fn))
					);

					if (nfinfo(dir) == 1 and cfinfo(dir) < 0) {
						cfinfo(dir) = 0;
					}
				}
			}
			const auto& en = neighbor.edge_neighbor;
			if (en[0] >= 0) {
				const int
					dir1 = 1 + std::min((en[0] + 1) % 3, (en[0] + 2) % 3),
					dir2 = 1 + std::max((en[0] + 1) % 3, (en[0] + 2) % 3),
					sign1 = [&en](){
						if (en[1] < 0) return -1;
						else return +1;}(),
					sign2 = [&en](){
						if (en[2] < 0) return -1;
						else return +1;}();

				if (nfinfo(-sign1*dir1) == 1 or nfinfo(-sign2*dir2) == 1) {
					if (cfinfo(sign1*dir1) < 0) cfinfo(sign1*dir1) = 0;
					if (cfinfo(sign2*dir2) < 0) cfinfo(sign2*dir2) = 0;
				}
			}
		}
	}

	grid.update_copies_of_remote_neighbors();

	// rest are dont_solve
	for (const auto& cell: grid.local_cells()) {
		auto& cfinfo = FInfo.data(*cell.data);
		for (auto dir: all_dirs) {
			if (cfinfo(dir) < -1) cfinfo(dir) = -1;
		}
	}
	grid.update_copies_of_remote_neighbors();
	Cell::set_transfer_all(false, CInfo.type(), FInfo.type());
	CInfo.type().is_stale = false;
	FInfo.type().is_stale = false;

} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
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
		const auto& ccinfo = CInfo.data(*cell.data);
		if (ccinfo < 1) {
			grid.dont_refine(cell.id);
			grid.dont_unrefine(cell.id);
			continue;
		}
		for (const auto& neighbor: cell.neighbors_of) {
			const auto& ncinfo = CInfo.data(*neighbor.data);
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


template<
	class Grid,
	class Geometries,
	class Boundaries,
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Magnetic_Field_Getter,
	class Face_Magnetic_Field_Getter,
	class Solver_Info_Getter,
	class Face_Info_Getter,
	class Substepping_Period_Getter
> void apply_boundaries(
	Grid& grid,
	Geometries& geometries,
	Boundaries& boundaries,
	const double& simulation_time,
	const double& proton_mass,
	const double& adiabatic_index,
	const double& vacuum_permeability,
	const Mass_Density_Getter& Mas,
	const Momentum_Density_Getter& Mom,
	const Total_Energy_Density_Getter& Nrj,
	const Magnetic_Field_Getter& Mag,
	const Face_Magnetic_Field_Getter& Face_B,
	const Solver_Info_Getter& SInfo,
	const Face_Info_Getter& FInfo,
	const Substepping_Period_Getter& Substep
) try {
	pamhd::mhd::apply_fluid_boundaries(
		grid, boundaries, geometries, simulation_time,
		Mas, Mom, Nrj, Mag, SInfo, proton_mass,
		adiabatic_index, vacuum_permeability
	);

	pamhd::mhd::apply_magnetic_field_boundaries_staggered(
		grid, boundaries, geometries,
		simulation_time, Face_B, FInfo
	);

	pamhd::mhd::update_B_consistency(
		0, grid.local_cells(), grid, Mas, Mom, Nrj,
		Mag, Face_B, SInfo, Substep, adiabatic_index,
		vacuum_permeability, true
	);
} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}

}} // namespaces


#endif // ifndef PAMHD_MHD_BOUNDARIES_HPP
