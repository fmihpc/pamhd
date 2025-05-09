/*
Interpolation functions between cells, faces, edges, vertices.

Copyright 2025 Finnish Meteorological Institute
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

#ifndef PAMHD_MATH_INTERPOLATION_HPP
#define PAMHD_MATH_INTERPOLATION_HPP


#include "string"

#include "grid/amr.hpp"


namespace pamhd {
namespace math {


/*! Interpolates cell volume average to cell vertices.

Assumes that both source and target are same type.
Assumes target variable is of type grid::Vertex_Type.

Interpolation is performed for cells with type > 0,
cells with type < 0 are ignored.
*/
template<
	class Grid,
	class Source_Getter,
	class Target_Getter,
	class Cell_Type_Getter
> void volume2vertex(
	Grid& grid,
	const Source_Getter& Src,
	const Target_Getter& Tgt,
	const Cell_Type_Getter& Type
) try {
	using std::abs;

	using Cell = std::remove_reference_t<decltype(grid)>::cell_data_type;

	Tgt.type().is_stale = true;

	bool update_copies = false;
	if (Src.type().is_stale) {
		update_copies = true;
		Cell::set_transfer_all(true, Src.type());
		Src.type().is_stale = false;
	}
	if (Type.type().is_stale) {
		update_copies = true;
		Cell::set_transfer_all(true, Type.type());
		Type.type().is_stale = false;
	}
	if (update_copies) {
		grid.update_copies_of_remote_neighbors();
		Cell::set_transfer_all(false, Src.type(), Type.type());
	}

	for (const auto& cell: grid.local_cells()) {
		if (Type.data(*cell.data) <= 0) continue;

		grid::Vertex_Type<int> nr_values;
		auto& tgt = Tgt.data(*cell.data);

		for (int x: {-1, +1})
		for (int y: {-1, +1})
		for (int z: {-1, +1}) {
			tgt(x, y, z) = Src.data(*cell.data);
			nr_values(x, y, z) = 1;
		}

		for (const auto& neighbor: cell.neighbors_of) {
			if (Type.data(*neighbor.data) < 0) continue;
			const auto& src = Src.data(*neighbor.data);

			const auto& fn = neighbor.face_neighbor;
			if (abs(fn) == 1) {
				for (int d1: {-1, +1} )
				for (int d2: {-1, +1} ) {
					tgt(fn, d1, d2) += src;
					nr_values(fn, d1, d2) += 1;
				}
			}
			if (abs(fn) == 2) {
				for (int d1: {-1, +1} )
				for (int d2: {-1, +1} ) {
					tgt(d1, fn, d2) += src;
					nr_values(d1, fn, d2) += 1;
				}
			}
			if (abs(fn) == 3) {
				for (int d1: {-1, +1} )
				for (int d2: {-1, +1} ) {
					tgt(d1, d2, fn) += src;
					nr_values(d1, d2, fn) += 1;
				}
			}

			const auto& en = neighbor.edge_neighbor;
			if (en[0] == 0) {
				for (int side: {-1, +1}) {
					tgt(side, en[1], en[2]) += src;
					nr_values(side, en[1], en[2]) += 1;
				}
			}
			if (en[0] == 1) {
				for (int side: {-1, +1}) {
					tgt(en[1], side, en[2]) += src;
					nr_values(en[1], side, en[2]) += 1;
				}
			}
			if (en[0] == 2) {
				for (int side: {-1, +1}) {
					tgt(en[1], en[2], side) += src;
					nr_values(en[1], en[2], side) += 1;
				}
			}

			const auto& vn = neighbor.vertex_neighbor;
			if (vn[0] != 0 and vn[1] != 0 and vn[2] != 0) {
				tgt(vn[0], vn[1], vn[2]) += src;
				nr_values(vn[0], vn[1], vn[2]) += 1;
			}
		}

		for (int x: {-1, +1})
		for (int y: {-1, +1})
		for (int z: {-1, +1}) {
			if (nr_values(x, y, z) < 1 or nr_values(x, y, z) > 8) {
				throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
			}
			tgt(x, y, z) /= nr_values(x, y, z);
		}
	}
} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


/*! Interpolates cell vertex data to cell center.

Assumes that both source and target are same type.
Assumes source variable is of type grid::Vertex_Type.

Interpolation is performed for cells with type > 0.
*/
template<
	class Grid,
	class Source_Getter,
	class Target_Getter,
	class Cell_Type_Getter
> void vertex2volume(
	Grid& grid,
	const Source_Getter& Src,
	const Target_Getter& Tgt,
	const Cell_Type_Getter& Type
) try {
	Tgt.type().is_stale = true;
	for (const auto& cell: grid.local_cells()) {
		if (Type.data(*cell.data) <= 0) continue;
		Tgt.data(*cell.data)
			=(Src.data(*cell.data)(-1, -1, -1)
			+ Src.data(*cell.data)(-1, -1, +1)
			+ Src.data(*cell.data)(-1, +1, -1)
			+ Src.data(*cell.data)(-1, +1, +1)
			+ Src.data(*cell.data)(+1, -1, -1)
			+ Src.data(*cell.data)(+1, -1, +1)
			+ Src.data(*cell.data)(+1, +1, -1)
			+ Src.data(*cell.data)(+1, +1, +1)) / 8;
	}
} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


/*! Interpolates cell faces to cell center.

Assumes source is grid::Face_Type and target is 3d.
Data from faces with normal in N dimension is interpolated
to component N of target.

Interpolation is performed for cells with type > 0.
*/
template<
	class Grid,
	class Source_Getter,
	class Target_Getter,
	class Cell_Type_Getter
> void face2volume(
	Grid& grid,
	const Source_Getter& Src,
	const Target_Getter& Tgt,
	const Cell_Type_Getter& Type
) try {
	Tgt.type().is_stale = true;
	for (const auto& cell: grid.local_cells()) {
		if (Type.data(*cell.data) <= 0) continue;
		for (size_t dim: {0, 1, 2}) {
			Tgt.data(*cell.data)[dim]
				=(Src.data(*cell.data)(dim, -1)
				+ Src.data(*cell.data)(dim, +1)) / 2;
		}
	}
} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


template<
	class Grid,
	class Source_Getter,
	class Target_Getter,
	class Cell_Type_Getter
> void vertex2face(
	Grid& grid,
	const Source_Getter& Src,
	const Target_Getter& Tgt,
	const Cell_Type_Getter& Type
) try {
} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


/*! Interpolates cell edges to cell vertices.

Assumes source is 1d Edge_Type and target is 3d Vertex_Type.
Data from edges along N dimension is interpolated
to component N of targets.

Interpolation is performed for cells with type > 0,
cells with type < 0 are ignored.
*/
template<
	class Grid,
	class Source_Getter,
	class Target_Getter,
	class Cell_Type_Getter
> void edge2vertex(
	Grid& grid,
	const Source_Getter& Src,
	const Target_Getter& Tgt,
	const Cell_Type_Getter& Type
) try {
	using std::abs;

	using Cell = std::remove_reference_t<decltype(grid)>::cell_data_type;

	Tgt.type().is_stale = true;

	bool update_copies = false;
	if (Src.type().is_stale) {
		update_copies = true;
		Cell::set_transfer_all(true, Src.type());
		Src.type().is_stale = false;
	}
	if (Type.type().is_stale) {
		update_copies = true;
		Cell::set_transfer_all(true, Type.type());
		Type.type().is_stale = false;
	}
	if (update_copies) {
		grid.update_copies_of_remote_neighbors();
		Cell::set_transfer_all(false, Src.type(), Type.type());
	}

	for (const auto& cell: grid.local_cells()) {
		if (Type.data(*cell.data) <= 0) continue;

		grid::Vertex_Type<int> nr_values;
		auto& tgt = Tgt.data(*cell.data);

		tgt(-1, -1, -1) =
		tgt(-1, -1, +1) =
		tgt(-1, +1, -1) =
		tgt(-1, +1, +1) =
		tgt(+1, -1, -1) =
		tgt(+1, -1, +1) =
		tgt(+1, +1, -1) =
		tgt(+1, +1, +1) = {0, 0, 0};
		nr_values(-1, -1, -1) =
		nr_values(-1, -1, +1) =
		nr_values(-1, +1, -1) =
		nr_values(-1, +1, +1) =
		nr_values(+1, -1, -1) =
		nr_values(+1, -1, +1) =
		nr_values(+1, +1, -1) =
		nr_values(+1, +1, +1) = 0;

		for (int side: {-1, +1})
		for (int perp1: {-1, +1})
		for (int perp2: {-1, +1}) {
			tgt(side, perp1, perp2)[0] += Src.data(*cell.data)(0, perp1, perp2);
			nr_values(side, perp1, perp2) += 1;

			tgt(perp1, side, perp2)[1] += Src.data(*cell.data)(1, perp1, perp2);
			nr_values(perp1, side, perp2) += 1;

			tgt(perp1, perp2, side)[2] += Src.data(*cell.data)(2, perp1, perp2);
			nr_values(perp1, perp2, side) += 1;
		}

		for (const auto& neighbor: cell.neighbors_of) {
			const auto& fn = neighbor.face_neighbor;
			if (fn == 0) continue;
			if (Type.data(*neighbor.data) < 0) continue;

			const size_t dim = abs(fn) - 1;
			const auto& src = Src.data(*neighbor.data);
			if (abs(fn) == 1) {
				for (int d1: {-1, +1})
				for (int d2: {-1, +1}) {
					tgt(fn, d1, d2)[dim] += src(dim, d1, d2);
					nr_values(fn, d1, d2) += 1;
				}
			}
			if (abs(fn) == 2) {
				for (int d1: {-1, +1})
				for (int d2: {-1, +1}) {
					tgt(d1, fn, d2)[dim] += src(dim, d1, d2);
					nr_values(d1, fn, d2) += 1;
				}
			}
			if (abs(fn) == 3) {
				for (int d1: {-1, +1})
				for (int d2: {-1, +1}) {
					tgt(d1, d2, fn)[dim] += src(dim, d1, d2);
					nr_values(d1, d2, fn) += 1;
				}
			}
		}

		for (int x: {-1, +1})
		for (int y: {-1, +1})
		for (int z: {-1, +1}) {
			if (nr_values(x, y, z) < 1 or nr_values(x, y, z) > 8) {
				throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
			}
			for (size_t dim: {0, 1, 2}) {
				tgt(x, y, z)[dim] /= nr_values(x, y, z);
			}
		}
	}
} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


}} // namespaces

#endif // ifndef PAMHD_MATH_INTERPOLATION_HPP
