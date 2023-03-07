/*
Class for handling the initial condition of a simulation in one region.

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

#ifndef PAMHD_BOUNDARIES_INITIAL_CONDITION_HPP
#define PAMHD_BOUNDARIES_INITIAL_CONDITION_HPP


#include "algorithm"
#include "array"
#include "stdexcept"
#include "string"
#include "type_traits"
#include "utility"

#include "rapidjson/document.h"
#include "mpParser.h"

#include "boundaries/common.hpp"
#include "math/expression.hpp"


namespace pamhd {
namespace boundaries {


/*!
Simulation initial condition creatable from json data.
*/
template<
	class Geometry_Id,
	class Variable
> class Initial_Condition
{
public:


	/*!
	Prepares initial condition from given rapidjson object.

	Given object must have "geometry-id" key with integer
	value and "value" key with either integer, string or
	object as value. If object then it must have keys
	denoting coordinate planes where data points are
	located and a "data" key.
	Coordinates and data are arrays of numbers.
	Coordinate keys are either "x", "y" and "z" for
	cartesian geometry or "radius", "lat" and "lon" for
	spherical data. "data" is in order
	(x_1, y_1, z_1), (x_2, y_1, z_1), ...

	Geometry ids refer to geometries of
	source/boundaries/geometries.hpp and must be >= 0.

	Example jsons which could be given as object:
	\verbatim
	{"geometry-id": 0, "value": 123}
	\endverbatim
	\verbatim
	{"geometry-id": 1, "value": "sin(t)"}
	\endverbatim
	\verbatim
	{
		"geometry-id": 0,
		"value": {
			"x": [1, 2, 3, 4],
			"y": [-3, 2, 5, 8],
			"z": [0.1, 0.2, 40, 123],
			"data": [1, 2, 3, ..., 64]
		}
	}
	\endverbatim
	*/
	void set(const rapidjson::Value& object)
	{
		if (not object.HasMember("geometry-id")) {
			throw std::invalid_argument(__FILE__ ": object doesn't have a geometry-id key.");
		}
		const auto& geometry_id = object["geometry-id"];

		if (not geometry_id.IsUint()) {
			throw std::invalid_argument(__FILE__ ": geometry-id is not unsigned int.");
		}

		this->geometry_id = geometry_id.GetUint();


		if (not object.HasMember("value")) {
			throw std::invalid_argument(__FILE__ ": object doesn't have a value key.");
		}
		const auto& value = object["value"];

		try {
			fill_variable_value_from_json(
				value,
				this->init_cond_type,
				this->number_value,
				this->math_expression,
				this->coordinates,
				this->coordinate_type,
				this->data
			);
		} catch (const std::invalid_argument& error) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + ") "
				+ " Couldn't set variable from json object: "
				+ error.what()
			);
		}
	}


	typename Variable::data_type get_data(
		const double& t,
		const double& x,
		const double& y,
		const double& z,
		const double& radius,
		const double& latitude,
		const double& longitude
	) {
		switch (this->init_cond_type) {
		case 1:
			return this->number_value;

		case 2:
			return this->math_expression.evaluate(
				t, x, y, z, radius, latitude, longitude
			);

		// nearest neighbor interpolation in point cloud
		case 3:
			return this->data[
				find_data(
					[&]()
						-> std::array<double, 3>
					{
						switch (this->coordinate_type) {
						case 1:
							return {x, y, z};
						case 2:
							return {radius, latitude, longitude};
						default:
							throw std::out_of_range(__FILE__ ": Invalid coordinate type.");
						}
					}(),
					this->coordinates
				)
			];

		default:
			throw std::out_of_range(__FILE__ ": Invalid initial condition type.");
		}
	}


	Geometry_Id get_geometry_id() const
	{
		return this->geometry_id;
	}

	void set_geometry_id(const Geometry_Id& id)
	{
		this->geometry_id = id;
	}


private:

	// region where this initial condition is applied
	Geometry_Id geometry_id;

	// initial condition if number was given in json file
	typename Variable::data_type number_value;

	// initial condition if string was given in json file
	math::Expression<Variable> math_expression;

	// initial condition if object was given in json file
	std::array<std::vector<double>, 3> coordinates; // e.g. x,y,z coordinates of points
	int coordinate_type = -1; // cartesian == 1, geographic == 2

	// data in order x[0],y[0],z[0]; x[1],y[0],z[0]; ...
	std::vector<typename Variable::data_type> data;

	int init_cond_type = -1; // number == 1, string == 2, object == 3
};

}} // namespaces

#endif // ifndef PAMHD_BOUNDARIES_INITIAL_CONDITION_HPP
