/*
Class for handling all value boundaries of one simulation variable.

Copyright 2014, 2015, 2016, 2017 Ilja Honkonen
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

#ifndef PAMHD_BOUNDARIES_VALUE_BOUNDARIES_HPP
#define PAMHD_BOUNDARIES_VALUE_BOUNDARIES_HPP


#include "algorithm"
#include "array"
#include "stdexcept"
#include "string"
#include "type_traits"
#include "utility"

#include "rapidjson/document.h"
#include "mpParser.h"

#include "boundaries/geometries.hpp"
#include "boundaries/value_boundary.hpp"


namespace pamhd {
namespace boundaries {


/*!
Collection of value boundaries creatable json data.
*/
template<
	class Geometry_Id,
	class Variable
> class Value_Boundaries
{
public:

	Value_Boundaries() : boundaries(boundaries_rw) {}

	Value_Boundaries(
		const Value_Boundaries<Geometry_Id, Variable>& other
	) : boundaries(boundaries_rw), boundaries_rw(other.boundaries_rw) {}

	Value_Boundaries(
		Value_Boundaries<Geometry_Id, Variable>&& other
	) : boundaries(boundaries_rw), boundaries_rw(std::move(other.boundaries_rw)) {}


	/*!
	Prepares value boundaries from given rapidjson object.

	Example json object:
	\verbatim
	{
		"value-boundaries": [
			{"geometry-id": 0, "time-stamps": [1], "values": [123]}
			{"geometry-id": 1, "time-stamps": [-1, 2], "values": ["sin(t)", "cos(x)"]}
		]
	}
	\endverbatim
	*/
	void set(const rapidjson::Value& object)
	{
		using std::to_string;

		if (not object.HasMember("value-boundaries")) {
			return;
		}
		const auto& json_bdys = object["value-boundaries"];

		if (not json_bdys.IsArray()) {
			throw std::invalid_argument(__FILE__ ": value-boundaries is not an array.");
		}

		this->boundaries_rw.clear();
		this->boundaries_rw.resize(json_bdys.Size());

		for (size_t i = 0; i < json_bdys.Size(); i++) {
			try {
				this->boundaries_rw[i].set(json_bdys[i]);
			} catch (const std::invalid_argument& error) {
				throw std::invalid_argument(
					std::string(__FILE__ "(") + to_string(__LINE__) + ")"
					+ ": Couldn't set value boundary number " + to_string(i + 1) + ": "
					+ error.what()
				);
			}
		}
	}


	size_t get_number_of_boundaries() const
	{
		return this->boundaries.size();
	}

	/*!
	Boundary id is in range 0..number of boundaries - 1
	*/
	const Value_Boundary<Geometry_Id, Variable>& get_value_boundary(const size_t& boundary_id) const
	{
		return this->boundaries[boundary_id];
	}


	typename Variable::data_type get_data(
		const size_t& boundary_index,
		const double& t,
		const double& x,
		const double& y,
		const double& z,
		const double& radius,
		const double& latitude,
		const double& longitude
	) {
		if (boundary_index >= this->boundaries.size()) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + ") "
				+ "Invalid boundary index: " + std::to_string(boundary_index)
				+ ", should be < " + std::to_string(this->boundaries.size())
			);
		}

		return this->boundaries_rw[boundary_index].get_data(
			t, x, y, z, radius, latitude, longitude
		);
	}


	//! all value boundaries of one simulation variable
	const std::vector<Value_Boundary<Geometry_Id, Variable>>& boundaries;


private:

	std::vector<Value_Boundary<Geometry_Id, Variable>> boundaries_rw;

};

}} // namespaces

#endif // ifndef PAMHD_BOUNDARIES_VALUE_BOUNDARIES_HPP
