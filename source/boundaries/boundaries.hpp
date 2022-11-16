/*
Class for handling boundaries of one simulation variable.

Copyright 2016, 2017 Ilja Honkonen
Copyright 2019 Finnish Meteorological Institute
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

#ifndef PAMHD_BOUNDARIES_BOUNDARIES_HPP
#define PAMHD_BOUNDARIES_BOUNDARIES_HPP


#include "stdexcept"
#include "string"
#include "type_traits"
#include "utility"
#include "vector"

#include "dccrg.hpp"

#include "boundaries/copy_boundaries.hpp"
#include "boundaries/geometries.hpp"
#include "boundaries/value_boundaries.hpp"


namespace pamhd {
namespace boundaries {


template<
	class Cell_Id,
	class Geometry_Id,
	class Variable
> class Boundaries
{
public:

	/*!
	Creates boundaries for one simulation variable from given json data.

	\see Value_Boundaries::set and Copy_Boundaries::set for requirements
	on json data.
	*/
	void set(const rapidjson::Value& object)
	{
		try {
			this->value_boundaries.set(object);
		} catch (const std::invalid_argument& error) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Couldn't set value boundaries of variable "
				+ Variable::get_option_name() + ": "
				+ error.what()
			);
		}
		try {
			this->copy_boundaries.set(object);
		} catch (const std::invalid_argument& error) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Couldn't set copy boundaries of variable "
				+ Variable::get_option_name() + ": "
				+ error.what()
			);
		}
	}


	/*!
	Classifies cells into boundaries and dont_solve cells.

	Transfer of Cell_Type variable between processes
	must be enabled during call to this function.
	*/
	template<
		class Grid,
		class Geometries,
		class Cell_Type_Getter
	> void classify(
		Grid& grid,
		const Geometries& geometries,
		const Cell_Type_Getter& Cell_Type
	) {
		using std::get;
		using std::to_string;

		// cell types for internal use
		constexpr typename std::remove_reference<
			decltype(Cell_Type(*grid[0]))
		>::type
			normal_cell{1},
			dont_solve_cell{2},
			value_bdy_cell{3},
			copy_bdy_cell{4};

		// all cells are normal by default
		for (const auto& cell: grid.local_cells()) {
			Cell_Type(*cell.data) = normal_cell;
		}

		// value boundary cells
		for (const auto& value_bdy: this->value_boundaries.boundaries) {
			for (const auto& cell_id: geometries.get_cells(value_bdy.get_geometry_id())) {
				auto* const cell_data = grid[cell_id];
				if (cell_data == nullptr) {
					throw std::invalid_argument(
						std::string(__FILE__ "(") + to_string(__LINE__) + "): "
						"No data for cell " + to_string(cell_id)
					);
				}
				Cell_Type(*cell_data) = value_bdy_cell;
			}
		}

		// copy boundary cells
		for (const auto& copy_geom_id: this->copy_boundaries.geometry_ids) {
			for (const auto& cell_id: geometries.get_cells(copy_geom_id)) {
				auto* const cell_data = grid[cell_id];
				if (cell_data == nullptr) {
					throw std::invalid_argument(
						std::string(__FILE__ "(") + to_string(__LINE__) + "): "
						"No data for cell " + to_string(cell_id)
					);
				}
				if (Cell_Type(*cell_data) != normal_cell) {
					continue;
				}
				Cell_Type(*cell_data) = copy_bdy_cell;
			}
		}

		// tell other processes about this proc's cell types
		grid.update_copies_of_remote_neighbors();

		// turn boundary cells without normal neighbors into dont_solve
		for (const auto& cell: grid.local_cells()) {
			if (Cell_Type(*cell.data) == normal_cell) {
				continue;
			}

			bool has_normal_neighbor = false;
			for (const auto& neighbor: cell.neighbors_of) {
				// cell itself doesn't count
				if (cell.id == neighbor.id) {
					continue;
				}

				if (Cell_Type(*neighbor.data) == normal_cell) {
					has_normal_neighbor = true;
					break;
				}
			}

			if (not has_normal_neighbor) {
				Cell_Type(*cell.data) = dont_solve_cell;
			}
		}

		// tell other processes about dont solve cells
		grid.update_copies_of_remote_neighbors();

		// get sources of copy boundary cells
		for (const auto& cell: grid.local_cells()) {
			if (Cell_Type(*cell.data) != copy_bdy_cell) {
				continue;
			}

			std::vector<decltype(cell.id)> source{cell.id};
			for (const auto& neighbor: cell.neighbors_of) {
				if (Cell_Type(*neighbor.data) == normal_cell) {
					source.push_back(neighbor.id);
				}
			}

			// turn copy cell without sources into dont solve
			if (source.size() < 2) {
				Cell_Type(*cell.data) = dont_solve_cell;
			} else {
				this->copy_boundaries.push_back_source(source);
			}
		}

		// tell other processes about newest dont solve cells
		grid.update_copies_of_remote_neighbors();

		for (const auto& cell: grid.local_cells()) {
			switch (Cell_Type(*cell.data)) {
				case dont_solve_cell:
					this->dont_solve_cells.push_back(cell.id);
					break;
				case value_bdy_cell:
					this->value_bdy_cells.push_back(cell.id);
					break;
				default:
					break;
			}
		}
	}


	const std::vector<Cell_Id>& get_dont_solve_cells() const
	{
		return this->dont_solve_cells;
	}

	const std::vector<Cell_Id>& get_value_boundary_cells() const
	{
		return this->value_bdy_cells;
	}

	const std::vector<std::vector<Cell_Id>>& get_copy_boundary_cells() const
	{
		return this->copy_boundaries.copy_sources;
	}


	size_t get_number_of_value_boundaries() const
	{
		return this->value_boundaries.get_number_of_boundaries();
	}

	const Value_Boundary<Geometry_Id, Variable>& get_value_boundary(const size_t& boundary_id) const
	{
		return this->value_boundaries.get_value_boundary(boundary_id);
	}



private:

	Value_Boundaries<Geometry_Id, Variable> value_boundaries;
	Copy_Boundaries<Cell_Id, Geometry_Id, Variable> copy_boundaries;

	std::vector<Cell_Id> dont_solve_cells, value_bdy_cells;
};

}} // namespaces

#endif // ifndef PAMHD_BOUNDARIES_BOUNDARIES_HPP
