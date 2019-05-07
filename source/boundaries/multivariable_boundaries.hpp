/*
Class for handling boundaries of all simulation variables.

Copyright 2016, 2017 Ilja Honkonen
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

#ifndef PAMHD_BOUNDARIES_MULTIVARIABLE_BOUNDARIES_HPP
#define PAMHD_BOUNDARIES_MULTIVARIABLE_BOUNDARIES_HPP


#include "stdexcept"
#include "string"
#include "type_traits"
#include "utility"
#include "vector"

#include "dccrg.hpp"

#include "boundaries/geometries.hpp"
#include "boundaries/boundaries.hpp"


namespace pamhd {
namespace boundaries {


template<
	class Cell_Id,
	class Geometry_Id,
	class... Variables
> class Multivariable_Boundaries {};

// stops recursion over simulation variables, does nothing
template<
	class Cell_Id,
	class Geometry_Id
> class Multivariable_Boundaries<Cell_Id, Geometry_Id>
{
public:
	void get_normal_cells() const {}
	void get_dont_solve_cells() const {}
	void get_value_boundary_cells() const {}
	void get_copy_boundary_cells() const {}
	void set(const rapidjson::Value&) const {}
	template<
		class Grid,
		class Geoms,
		class Cell_Type_Getter
	> void classify(
		Grid&,
		const Geoms&,
		const Cell_Type_Getter&
	) const {}
	void get_number_of_value_boundaries() const {}
	void get_value_boundary() const {}
};

template<
	class Cell_Id,
	class Geometry_Id,
	class Current_Variable,
	class... Rest
> class Multivariable_Boundaries<
	Cell_Id,
	Geometry_Id,
	Current_Variable,
	Rest...
> :
	public Multivariable_Boundaries<Cell_Id, Geometry_Id, Rest...>
{
public:

	using Multivariable_Boundaries<Cell_Id, Geometry_Id, Rest...>::get_normal_cells;
	using Multivariable_Boundaries<Cell_Id, Geometry_Id, Rest...>::get_dont_solve_cells;
	using Multivariable_Boundaries<Cell_Id, Geometry_Id, Rest...>::get_value_boundary_cells;
	using Multivariable_Boundaries<Cell_Id, Geometry_Id, Rest...>::get_copy_boundary_cells;
	using Multivariable_Boundaries<Cell_Id, Geometry_Id, Rest...>::get_number_of_value_boundaries;
	using Multivariable_Boundaries<Cell_Id, Geometry_Id, Rest...>::get_value_boundary;


	void set(const rapidjson::Value& object)
	{
		const std::string name = Current_Variable::get_option_name();
		if (not object.HasMember(name.c_str())) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Object doesn't have a " + name + " key."
			);
		}

		try {
			this->boundaries.set(object[name.c_str()]);
		} catch (const std::invalid_argument& error) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Couldn't set boundaries of variable "
				+ Current_Variable::get_option_name() + ": "
				+ error.what()
			);
		}

		Multivariable_Boundaries<Cell_Id, Geometry_Id, Rest...>::set(object);
	}


	/*!
	Classifies cell in grid into normal, boundary and dont_solve cells.

	Transfer of Cell_Type variable between processes must have been
	enabled before calling this function.
	*/
	template<
		class Grid,
		class Geoms,
		class Cell_Type_Getter
	> void classify(
		Grid& grid,
		const Geoms& geometries,
		const Cell_Type_Getter& Cell_Type
	) {
		this->boundaries.classify(grid, geometries, Cell_Type);

		Multivariable_Boundaries<
			Cell_Id,
			Geometry_Id,
			Rest...
		>::classify(grid, geometries, Cell_Type);
	}


	const std::vector<Cell_Id>& get_dont_solve_cells(
		const Current_Variable&
	) const {
		return this->boundaries.get_dont_solve_cells();
	}

	const std::vector<Cell_Id>& get_value_boundary_cells(
		const Current_Variable&
	) const {
		return this->boundaries.get_value_boundary_cells();
	}

	const std::vector<std::vector<Cell_Id>>& get_copy_boundary_cells(
		const Current_Variable&
	) const {
		return this->boundaries.get_copy_boundary_cells();
	}


	size_t get_number_of_value_boundaries(const Current_Variable&) const
	{
		return this->boundaries.get_number_of_value_boundaries();
	}

	const Value_Boundary<Geometry_Id, Current_Variable>& get_value_boundary(
		const Current_Variable&,
		const size_t& boundary_id
	) const {
		return this->boundaries.get_value_boundary(boundary_id);
	}
	Value_Boundary<Geometry_Id, Current_Variable>& get_value_boundary(
		const Current_Variable&,
		const size_t& boundary_id
	) {
		return const_cast<Value_Boundary<Geometry_Id, Current_Variable>&>(
			static_cast<
				const Multivariable_Boundaries<
					Cell_Id,
					Geometry_Id,
					Current_Variable,
					Rest...
				>&
			>(*this).boundaries.get_value_boundary(boundary_id)
		);
	}



private:

	Boundaries<Cell_Id, Geometry_Id, Current_Variable> boundaries;
};

}} // namespaces

#endif // ifndef PAMHD_BOUNDARIES_MULTIVARIABLE_BOUNDARIES_HPP
