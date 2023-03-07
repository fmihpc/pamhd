/*
Functionality used by different parts of PAMHD boundaries.

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

#ifndef PAMHD_BOUNDARIES_COMMON_HPP
#define PAMHD_BOUNDARIES_COMMON_HPP


#include "iterator"
#include "string"
#include "type_traits"

#include "rapidjson/document.h"

#include "math/expression.hpp"


namespace pamhd {
namespace boundaries {


/*!
Returns value of type Variable::data_type when it's signed integral.

Reads a signed integral number from given rapidjson object, throws
if given object contains a value of different type.

Object is interpreted as signed 64 bit integer, if it doesn't fit
into Variable::data_type this function should be specialized for
smaller integers (use https://github.com/nasailja/pamhd/pulls or
submit https://github.com/nasailja/pamhd/issues).

Example:
\verbatim
// Variable = struct {using data_type = int};
const char* json = "[123]"
rapidjson::Document doc; doc.Parse(json);
typename Variable::data_type var = get_json_value(doc[0]);
\endverbatim
*/
template<class Variable> typename std::enable_if<
	std::is_integral<typename Variable::data_type>::value
		&& std::is_signed<typename Variable::data_type>::value,
	typename Variable::data_type
>::type get_json_value(const rapidjson::Value& object)
{
	if (not object.IsInt64()) {
		throw std::invalid_argument(
			std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
			+ "Object isn't of signed integral type."
		);
	}
	return object.GetInt64();
}

//! returns value of type Variable::data_type when it's unsigned integral.
template<class Variable> typename std::enable_if<
	std::is_integral<typename Variable::data_type>::value
		&& !std::is_signed<typename Variable::data_type>::value,
	typename Variable::data_type
>::type get_json_value(const rapidjson::Value& object)
{
	if (not object.IsUint64()) {
		throw std::invalid_argument(
			std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
			+ "Object isn't of unsigned integral type."
		);
	}
	return object.GetUint64();
}

//! returns value of type Variable::data_type when it's floating point.
template<class Variable> typename std::enable_if<
	std::is_floating_point<typename Variable::data_type>::value,
	typename Variable::data_type
>::type get_json_value(const rapidjson::Value& object)
{
	if (not object.IsNumber()) {
		throw std::invalid_argument(
			std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
			+ "Object isn't of floating point type."
		);
	}
	return object.GetDouble();
}

//! returns value of type Variable::data_type when it's array of signed integral values.
template<class Variable> typename std::enable_if<
	std::tuple_size<typename Variable::data_type>::value >= 0
		&& std::is_integral<typename Variable::data_type::value_type>::value
		&& std::is_signed<typename Variable::data_type::value_type>::value,
	typename Variable::data_type
>::type get_json_value(const rapidjson::Value& object)
{
	constexpr auto size = std::tuple_size<typename Variable::data_type>::value;

	if (not object.IsArray()) {
		throw std::invalid_argument(
			std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
			+ "Object isn't an array."
		);
	}
	if (object.Size() != size) {
		throw std::invalid_argument(
			std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
			+ "Array has wrong size: " + std::to_string(object.Size())
			+ ", should be " + std::to_string(size)
		);
	}

	typename Variable::data_type ret_val;
	for (size_t i = 0; i < size; i++) {
		if (not object[i].IsInt64()) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Element of array isn't signed integral number."
			);
		}
		ret_val[i] = object[i].GetInt64();
	}
	return ret_val;
}

//! returns value of type Variable::data_type when it's array of unsigned integral values.
template<class Variable> typename std::enable_if<
	std::tuple_size<typename Variable::data_type>::value >= 0
		&& std::is_integral<typename Variable::data_type::value_type>::value
		&& !std::is_signed<typename Variable::data_type::value_type>::value,
	typename Variable::data_type
>::type get_json_value(const rapidjson::Value& object)
{
	constexpr auto size = std::tuple_size<typename Variable::data_type>::value;

	if (not object.IsArray()) {
		throw std::invalid_argument(
			std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
			+ "Object isn't an array."
		);
	}
	if (object.Size() != size) {
		throw std::invalid_argument(
			std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
			+ "Array has wrong size: " + std::to_string(object.Size())
			+ ", should be " + std::to_string(size)
		);
	}

	typename Variable::data_type ret_val;
	for (size_t i = 0; i < size; i++) {
		if (not object[i].IsUint64()) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Element of array isn't unsigned integral number."
			);
		}
		ret_val[i] = object[i].GetUint64();
	}
	return ret_val;
}

//! returns value of type Variable::data_type when it's array of floating point values.
template<class Variable> typename std::enable_if<
	std::tuple_size<typename Variable::data_type>::value >= 0
		&& std::is_floating_point<typename Variable::data_type::value_type>::value,
	typename Variable::data_type
>::type get_json_value(const rapidjson::Value& object)
{
	constexpr auto size = std::tuple_size<typename Variable::data_type>::value;

	if (not object.IsArray()) {
		throw std::invalid_argument(
			std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
			+ "Object isn't an array."
		);
	}
	if (object.Size() != size) {
		throw std::invalid_argument(
			std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
			+ "Array has wrong size: " + std::to_string(object.Size())
			+ ", should be " + std::to_string(size)
		);
	}

	typename Variable::data_type ret_val;
	for (size_t i = 0; i < size; i++) {
		if (not object[i].IsNumber()) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Element of array isn't a number."
			);
		}
		ret_val[i] = object[i].GetDouble();
	}
	return ret_val;
}


#ifdef EIGEN_WORLD_VERSION

//! returns value of type Variable::data_type when it's Eigen::Vector of floating point values with compile time size.
template<class Variable> typename std::enable_if<
	Variable::data_type::ColsAtCompileTime == 1
		&& Variable::data_type::RowsAtCompileTime >= 1
		&& std::is_floating_point<typename Variable::data_type::value_type>::value,
	typename Variable::data_type
>::type get_json_value(const rapidjson::Value& object)
{
	constexpr auto size = Variable::data_type::RowsAtCompileTime;

	if (not object.IsArray()) {
		throw std::invalid_argument(
			std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
			+ "Object isn't an array."
		);
	}
	if (object.Size() != size) {
		throw std::invalid_argument(
			std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
			+ "Array has wrong size: " + std::to_string(object.Size())
			+ ", should be " + std::to_string(size)
		);
	}

	typename Variable::data_type ret_val;
	for (size_t i = 0; i < size; i++) {
		if (not object[i].IsNumber()) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Element of array isn't a number."
			);
		}
		ret_val(i) = object[i].GetDouble();
	}
	return ret_val;
}

#endif


/*!
Fills given variables from given json data, throws in case of error.

\see fill_variable_value_from_json().
*/
template<class Variable> void fill_variable_from_json_grid_data(
	const rapidjson::Value& object,
	std::array<std::vector<double>, 3>& coordinates,
	int& coordinate_type,
	std::vector<typename Variable::data_type>& data
) {
	if (object.HasMember("x")) {
		coordinate_type = 1;
	} else if (object.HasMember("radius")) {
		coordinate_type = 2;
	} else {
		throw std::invalid_argument(__FILE__ ": object doesn't have either x or radius key.");
	}

	std::array<const char*, 3> components{"x", "y", "z"};
	if (coordinate_type == 2) {
		components[0] = "radius";
		components[1] = "lat";
		components[2] = "lon";
	}

	// load coordinates of data points
	size_t component_i = 0;
	for (const auto& component: components) {
		if (not object.HasMember(component)) {
			throw std::invalid_argument(__FILE__ ": object object doesn't have required key.");
		}

		const auto& coord = object[component];
		if (not coord.IsArray()) {
			throw std::invalid_argument(__FILE__ ": coordinate isn't array.");
		}

		coordinates[component_i].clear();
		for (auto i = coord.Begin(); i != coord.End(); i++) {
			coordinates[component_i].push_back(i->GetDouble());
		}

		if (
			not std::is_sorted(
				coordinates[component_i].cbegin(),
				coordinates[component_i].cend()
			)
		) {
			throw std::invalid_argument(__FILE__ ": coordinates aren't sorted in non-descending order.");
		}

		component_i++;
	}

	// load data
	if (object.HasMember("data")) {
		const auto& json_data = object["data"];

		if (not json_data.IsArray()) {
			throw std::invalid_argument(__FILE__ ": data isn't array.");
		}

		data.clear();
		for (auto i = json_data.Begin(); i != json_data.End(); i++) {
			data.push_back(get_json_value<Variable>(*i));
		}
	} else {
		throw std::invalid_argument(__FILE__ ": object doesn't have a x or radius key.");
	}
}


/*!
Fills one of given arguments based on given rapidjson object.

If given object is a number value_type is set to 1 and
number_value is set to the number.

If given object is a string value_type is set to 2 and
math_expression is initialized to the string.

If given object is an object value_type is set to 3 and
coordinates, coordinate_type and data are set. In this case
given object must have a "data" key and either "x", "y", "z"
or "radius", "lat", "lon" keys. All values behind keys must be
arrays and length of "data" must equal the product of lengths
of either "x", "y", "z" or "radius", "lat", "lon". All arrays
must consist of numbers. data is filled in order where "x"
or "radius" coordinate increases fastest and "z" or "lon"
slowest, e.g. in last example below data[0] represents
(1, -3, 0.1), data[1] represents (2, -3, 0.1), data[4]
represents (1, -3, 0.1), data[16] represents (1, -3, 0.2), etc.

Example jsons which could be given as object:
\verbatim
123
\endverbatim
\verbatim
"sin(t)"
\endverbatim
\verbatim
{
	"x": [1, 2, 3, 4],
	"y": [-3, 2, 5, 8],
	"z": [0.1, 0.2, 40, 123],
	"data": [1, 2, 3, ..., 64]
}
\endverbatim
*/
template<class Variable> void fill_variable_value_from_json(
	const rapidjson::Value& object,
	int& value_type,
	typename Variable::data_type& number_value,
	math::Expression<Variable>& math_expression,
	std::array<std::vector<double>, 3>& coordinates,
	int& coordinate_type,
	std::vector<typename Variable::data_type>& data
) {
	// one scalar or vector value
	if (object.IsNumber() or object.IsArray()) {

		value_type = 1;
		number_value = get_json_value<Variable>(object);

	// math expression
	} else if (object.IsString()) {

		value_type = 2;
		math_expression.set_expression(object.GetString());

	// 3d grid of scalar values
	} else if (object.IsObject()) {

		value_type = 3;
		fill_variable_from_json_grid_data<Variable>(
			object,
			coordinates,
			coordinate_type,
			data
		);

	} else {
		throw std::invalid_argument(
			std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
			+ "Invalid object type."
		);
	}
}


/*!
Time series version of fill_variable_value_from_json()
*/
template<class Variable> void fill_variable_value_from_json(
	const rapidjson::Value& object,
	int& value_type,
	std::vector<typename Variable::data_type>& number_values,
	std::vector<math::Expression<Variable>>& math_expressions,
	std::array<std::vector<double>, 3>& coordinates,
	int& coordinate_type,
	std::vector<typename Variable::data_type>& data
) {
	if (not object.IsArray() and not object.IsObject()) {
		throw std::invalid_argument(
			std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
			+ "Not given an array or object."
		);
	}

	if (object.IsArray()) {
		if (object.Size() == 0) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Array is empty."
			);
		}

		if (object[0].IsNumber() or object[0].IsArray()) {
			value_type = 1;

			number_values.clear();
			number_values.reserve(object.Size());
			for (size_t i = 0; i < object.Size(); i++) {
				number_values.push_back(get_json_value<Variable>(object[i]));
			}

		} else if (object[0].IsString()) {
			value_type = 2;

			math_expressions.clear();
			math_expressions.resize(object.Size());
			for (size_t i = 0; i < object.Size(); i++) {
				math_expressions[i].set_expression(object[i].GetString());
			}
		} else {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "First object has invalid type."
			);
		}

	} else {
		value_type = 3;
		fill_variable_from_json_grid_data<Variable>(
			object,
			coordinates,
			coordinate_type,
			data
		);
	}
}


/*!
Returns closest value from given 3d cartesian grid.

Finds coordinate in all_coordinates that is closest to
request_cordinate and returns corresponding item from data.

Assumes that data is on a cartesian 3d grid and
in each dimension separately the location of all
points are given by all_coordinate[0], [1] and [2]
respectively.

For example if x == all_coordinates[0], etc. then
data[0] corresponds to point (x[0], y[0], z[0]) and if
x.size() == N then data[2*N] corresponds to point
(x[0], y[2], z[0]), etc.
*/
size_t find_data(
	const std::array<double, 3>& request_coordinate,
	const std::array<std::vector<double>, 3>& all_coordinates
) {
	std::array<size_t, 3> index{0, 0, 0};
	for (size_t i = 0; i < all_coordinates.size(); i++) {
		// location of first value not smaller than data point coordinate
		const auto after = std::lower_bound(
			all_coordinates[i].cbegin(),
			all_coordinates[i].cend(),
			request_coordinate[i]
		);

		if (after == all_coordinates[i].cbegin()) {
			continue;
		} else if (after == all_coordinates[i].cend()) {
			index[i] = all_coordinates[i].size() - 1;
		} else {
			const auto before = after - 1;
			if (
				std::abs(request_coordinate[i] - *before)
				<= std::abs(request_coordinate[i] - *after)
			) {
				index[i] = std::distance(all_coordinates[i].cbegin(), before);
			} else {
				index[i] = std::distance(all_coordinates[i].cbegin(), after);
			}
		}
	}

	const size_t data_i
		= index[0]
		+ index[1] * all_coordinates[0].size()
		+ index[2] * all_coordinates[0].size() * all_coordinates[1].size();

	return data_i;
}


}} // namespaces

#endif // ifndef PAMHD_BOUNDARIES_COMMON_HPP
