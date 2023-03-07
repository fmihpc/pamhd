/*
Class for handling all initial conditions of one simulation variable.

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

#ifndef PAMHD_BOUNDARIES_INITIAL_CONDITIONS_HPP
#define PAMHD_BOUNDARIES_INITIAL_CONDITIONS_HPP


#include "algorithm"
#include "array"
#include "stdexcept"
#include "string"
#include "type_traits"
#include "utility"

#include "rapidjson/document.h"
#include "mpParser.h"

#include "boundaries/initial_condition.hpp"
#include "math/expression.hpp"


namespace pamhd {
namespace boundaries {


/*!
Collection of simulation initial conditions creatable from json data.
*/
template<
	class Geometry_Id,
	class Variable
> class Initial_Conditions
{
public:

	/*!
	Prepares initial conditions from given rapidjson object.

	Given object must have a "default" key with corresponding
	value as a number, string (math expression), or object
	(3d point cloud data), \see Initial_Condition::set() for
	details.

	Example json object:
	\verbatim
	{
		"default": "sin(2*pi*x/123)",
		"initial-conditions": [
			{"geometry-id": 0, "value": 123}
			{"geometry-id": 1, "value": "sin(t)"}
		]
	}
	\endverbatim
	*/
	void set(const rapidjson::Value& object)
	{
		if (not object.HasMember("default")) {
			throw std::invalid_argument(__FILE__ ": object doesn't have a default key.");
		}
		const auto& default_ = object["default"];

		try {
			fill_variable_value_from_json(
				default_,
				this->default_init_cond_type,
				this->default_number_value,
				this->default_math_expression,
				this->default_coordinates,
				this->default_coordinate_type,
				this->default_data
			);
		} catch (const std::invalid_argument& error) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + ") "
				+ " Couldn't set variable from json object: "
				+ error.what()
			);
		}

		if (not object.HasMember("initial-conditions")) {
			return;
		}
		const auto& json_regions = object["initial-conditions"];

		if (not json_regions.IsArray()) {
			throw std::invalid_argument(__FILE__ ": initial-conditions is not an array.");
		}

		this->regions.clear();
		this->regions.resize(json_regions.Size());

		for (size_t i = 0; i < json_regions.Size(); i++) {
			this->regions[i].set(json_regions[i]);
		}
	}


	typename Variable::data_type get_default_data(
		const double& t,
		const double& x,
		const double& y,
		const double& z,
		const double& radius,
		const double& latitude,
		const double& longitude
	) {
		switch (this->default_init_cond_type) {
		case 1:
			return this->default_number_value;

		case 2:
			return this->default_math_expression.evaluate(
				t, x, y, z, radius, latitude, longitude
			);

		// nearest neighbor interpolation in point cloud
		case 3:
			return this->default_data[
				find_data(
					[&]()
						-> std::array<double, 3>
					{
						switch (this->default_coordinate_type) {
						case 1:
							return {x, y, z};
						case 2:
							return {radius, latitude, longitude};
						default:
							throw std::out_of_range(__FILE__ ": Invalid coordinate type.");
						}
					}(),
					this->default_coordinates
				)
			];

		default:
			throw std::out_of_range(__FILE__ ": Invalid initial condition type.");
		}
	}


	size_t get_number_of_regions() const
	{
		return this->regions.size();
	}

	/*!
	Region id is in range 0..number of regions - 1
	*/
	const Initial_Condition<Geometry_Id, Variable>& get_initial_condition(const size_t& region_id) const
	{
		return this->regions[region_id];
	}


	typename Variable::data_type get_data(
		const size_t& region_index,
		const double& t,
		const double& x,
		const double& y,
		const double& z,
		const double& radius,
		const double& latitude,
		const double& longitude
	) {
		return this->regions[region_index].get_data(
			t, x, y, z, radius, latitude, longitude
		);
	}


private:

	/*
	Variables for handling default initial condition, see
	corresponding names in initial_condition.hpp for more info
	*/
	typename Variable::data_type default_number_value;
	math::Expression<Variable> default_math_expression;
	std::array<std::vector<double>, 3> default_coordinates;
	int default_coordinate_type = -1;
	std::vector<typename Variable::data_type> default_data;
	int default_init_cond_type = -1;

	// regions of non-default initial condition loaded from json data.
	std::vector<Initial_Condition<Geometry_Id, Variable>> regions;

};

}} // namespaces

#endif // ifndef PAMHD_BOUNDARIES_INITIAL_CONDITIONS_HPP
