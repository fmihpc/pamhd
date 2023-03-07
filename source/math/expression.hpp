/*
Class for evaluating math expression of a simulation variable.

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

#ifndef PAMHD_BOUNDARIES_MATH_EXPRESSION_HPP
#define PAMHD_BOUNDARIES_MATH_EXPRESSION_HPP


#include "stdexcept"
#include "string"
#include "type_traits"
#include "utility"

#include "mpParser.h"


namespace pamhd {
namespace boundaries {


/*!
Variables available in expression, all are of floating point type:
	- t: simulation time
	- x, y, z: center of simulation cell in cartesian coordinates
	- radius, lat, lon: center in spherical coordinates
*/
template<class Variable> class Math_Expression
{
public:

	Math_Expression() :
		t_var(&t_val),
		x_var(&x_val), y_var(&y_val), z_var(&z_val),
		radius_var(&radius_val), lat_var(&lat_val), lon_var(&lon_val)
	{
		this->parser.DefineVar("t", this->t_var);
		this->parser.DefineVar("x", this->x_var);
		this->parser.DefineVar("y", this->y_var);
		this->parser.DefineVar("z", this->z_var);
		this->parser.DefineVar("radius", this->radius_var);
		this->parser.DefineVar("lat", this->lat_var);
		this->parser.DefineVar("lon", this->lon_var);
	}

	Math_Expression(const Math_Expression& other) :
		t_val(other.t_val),
		x_val(other.x_val), y_val(other.y_val), z_val(other.z_val),
		radius_val(other.radius_val), lat_val(other.lat_val), lon_val(other.lon_val),
		t_var(&t_val),
		x_var(&x_val), y_var(&y_val), z_var(&z_val),
		radius_var(&radius_val), lat_var(&lat_val), lon_var(&lon_val)
	{
		this->parser.DefineVar("t", this->t_var);
		this->parser.DefineVar("x", this->x_var);
		this->parser.DefineVar("y", this->y_var);
		this->parser.DefineVar("z", this->z_var);
		this->parser.DefineVar("radius", this->radius_var);
		this->parser.DefineVar("lat", this->lat_var);
		this->parser.DefineVar("lon", this->lon_var);
	}

	Math_Expression(Math_Expression&& other) = delete;


	void add_expression_variable(
		const std::string& name,
		mup::Variable& variable
	) {
		this->parser.DefineVar(name, variable);
	}


	void clear_expression_variables()
	{
		this->parser.ClearVar();
	}


	void set_expression(const std::string& expression)
	{
		this->parser.SetExpr(expression);
	}

	std::string get_expression() const
	{
		return this->parser.GetExpr();
	}


	typename Variable::data_type evaluate(
		const double& t,
		const double& x,
		const double& y,
		const double& z,
		const double& radius,
		const double& latitude,
		const double& longitude
	) {
		this->t_val = t;
		this->x_val = x;
		this->y_val = y;
		this->z_val = z;
		this->radius_val = radius;
		this->lat_val = latitude;
		this->lon_val = longitude;

		try {
			return this->evaluate_impl();
		} catch(const mup::ParserError& error) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Couldn't evaluate expression \"" + this->parser.GetExpr()
				+ "\" for variable " + Variable::get_name()
				+ ": " + error.GetMsg()
			);
		}
	}


private:

	mup::ParserX parser = mup::ParserX(mup::pckCOMMON | mup::pckNON_COMPLEX | mup::pckMATRIX | mup::pckUNIT);
	mup::Value t_val, x_val, y_val, z_val, radius_val, lat_val, lon_val;
	mup::Variable t_var, x_var, y_var, z_var, radius_var, lat_var, lon_var;


	//! used to evaluate math expression if Variable::data_type is integral
	template<class V = Variable> typename std::enable_if<
		std::is_integral<typename V::data_type>::value,
		typename V::data_type
	>::type evaluate_impl()
	{
		const auto& evaluated = this->parser.Eval();
		if (evaluated.GetType() != 'i' and evaluated.GetType() != 'f') {
			throw std::invalid_argument(
				std::string("Expression \"")
				+ this->parser.GetExpr()
				+ std::string("\" is neither integer nor floating point.")
			);
		}
		return evaluated.GetInteger();
	}


	//! used to evaluate math expression if Variable::data_type is floating point
	template<class V = Variable> typename std::enable_if<
		std::is_floating_point<typename V::data_type>::value,
		typename V::data_type
	>::type evaluate_impl()
	{
		const auto& evaluated = this->parser.Eval();
		if (evaluated.GetType() != 'i' and evaluated.GetType() != 'f') {
			throw std::invalid_argument(
				std::string("Expression \"")
				+ this->parser.GetExpr()
				+ std::string("\" is neither integer nor floating point.")
			);
		}
		return evaluated.GetFloat();
	}


	//! used to evaluate math expression if Variable::data_type is array of floating point numbers
	template<class V = Variable> typename std::enable_if<
		std::tuple_size<typename V::data_type>::value >= 0
			&& std::is_floating_point<typename V::data_type::value_type>::value,
		typename V::data_type
	>::type evaluate_impl()
	{
		const auto& evaluated
			= [this](){
				const auto& temp = this->parser.Eval();
				if (temp.GetType() != 'm') {
					throw std::invalid_argument(
						std::string("Expression \"")
						+ this->parser.GetExpr()
						+ std::string("\" is not an array.")
					);
				}
				return temp.GetArray();
			}();

		// TODO: merge with identical code in int version
		if (evaluated.GetRows() != 1) {
			throw std::invalid_argument(
				std::string("Invalid number of rows in expression \"")
				+ parser.GetExpr()
				+ std::string("\" for variable ")
				+ Variable::get_name()
				+ std::string(", should be 1\n")
			);
		}

		constexpr size_t size = std::tuple_size<typename V::data_type>::value;

		if (evaluated.GetCols() != size) {
			throw std::invalid_argument(
				std::string("Invalid number of columns in expression \"")
				+ parser.GetExpr()
				+ std::string("\" for variable ")
				+ Variable::get_name()
				+ std::string("\n")
			);
		}

		typename Variable::data_type ret_val;
		for (size_t i = 0; i < size; i++) {
			ret_val[i] = evaluated.At(int(i)).GetFloat();
		}

		return ret_val;
	}


	//! used to evaluate math expression if Variable::data_type is array of integral numbers
	template<class V = Variable> typename std::enable_if<
		std::tuple_size<typename V::data_type>::value >= 0
			&& std::is_integral<typename V::data_type::value_type>::value,
		typename V::data_type
	>::type evaluate_impl()
	{
		const auto& evaluated
			= [this](){
				const auto& temp = this->parser.Eval();
				if (temp.GetType() != 'm') {
					throw std::invalid_argument(
						std::string("Expression \"")
						+ this->parser.GetExpr()
						+ std::string("\" is not an array.")
					);
				}
				return temp.GetArray();
			}();

		if (evaluated.GetRows() != 1) {
			throw std::invalid_argument(
				std::string("Invalid number of rows in expression \"")
				+ parser.GetExpr()
				+ std::string("\" for variable ")
				+ Variable::get_name()
				+ std::string(", should be 1\n")
			);
		}

		constexpr size_t size = std::tuple_size<typename Variable::data_type>::value;

		if (evaluated.GetCols() != size) {
			throw std::invalid_argument(
				std::string("Invalid number of columns in expression \"")
				+ parser.GetExpr()
				+ std::string("\" for variable ")
				+ Variable::get_name()
				+ std::string("\n")
			);
		}

		typename Variable::data_type ret_val;
		for (size_t i = 0; i < size; i++) {
			ret_val[i] = evaluated.At(int(i)).GetInteger();
		}

		return ret_val;
	}


	#ifdef EIGEN_WORLD_VERSION

	//! used to evaluate math expression if Variable::data_type is Eigen::Matrix<double or float, N, 1>
	template<class V = Variable> typename std::enable_if<
		V::data_type::ColsAtCompileTime == 1
			&& V::data_type::RowsAtCompileTime >= 1
			&& std::is_floating_point<typename V::data_type::value_type>::value,
		typename V::data_type
	>::type evaluate_impl()
	{
		const auto& evaluated
			= [this](){
				const auto& temp = this->parser.Eval();
				if (temp.GetType() != 'm') {
					throw std::invalid_argument(
						std::string("Expression ")
						+ this->parser.GetExpr()
						+ std::string(" is not an array.")
					);
				}
				return temp.GetArray();
			}();

		if (evaluated.GetRows() != 1) {
			throw std::invalid_argument(
				std::string("Invalid number of rows in expression ")
				+ parser.GetExpr()
				+ std::string(" for variable ")
				+ Variable::get_name()
				+ std::string(", should be 1\n")
			);
		}

		constexpr size_t size = V::data_type::RowsAtCompileTime;

		if (evaluated.GetCols() != size) {
			throw std::invalid_argument(
				std::string("Invalid number of columns in expression ")
				+ parser.GetExpr()
				+ std::string(" for variable ")
				+ Variable::get_name()
				+ std::string("\n")
			);
		}

		typename Variable::data_type ret_val;
		for (size_t i = 0; i < size; i++) {
			ret_val(i) = evaluated.At(int(i)).GetFloat();
		}

		return ret_val;
	}

	#endif // ifdef EIGEN_WORLD_VERSION
};

}} // namespaces

#endif // ifndef PAMHD_BOUNDARIES_MATH_EXPRESSION_HPP
