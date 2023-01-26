/*
Class for handling grid options.

Copyright 2014, 2015, 2016, 2017 Ilja Honkonen
Copyright 2022, 2023 Finnish Meteorological Institute
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

#ifndef PAMHD_GRID_OPTIONS_HPP
#define PAMHD_GRID_OPTIONS_HPP


#include "array"
#include "cstdint"
#include "stdexcept"
#include "string"

#include "mpParser.h"
#include "rapidjson/document.h"


namespace pamhd {
namespace grid {

struct Number_Of_Cells {
	using data_type = std::array<uint64_t, 3>;
	static const std::string get_name() { return {"number of cells"}; }
	static const std::string get_option_name() { return {"cells"}; }
};

struct Volume {
	using data_type = std::array<double, 3>;
	static const std::string get_name() { return {"volume"}; }
	static const std::string get_option_name() { return {"volume"}; }
};

struct Start {
	using data_type = std::array<double, 3>;
	static const std::string get_name() { return {"start"}; }
	static const std::string get_option_name() { return {"start"}; }
};

struct Periodic {
	using data_type = std::array<bool, 3>;
	static const std::string get_name() { return {"periodic"}; }
	static const std::string get_option_name() { return {"periodic"}; }
};

class Options
{
public:
	Options() = default;
	Options(const Options& other) = default;
	Options(Options&& other) = default;

	Options(const rapidjson::Value& object)
	{
		this->set(object);
	};


	void set(const rapidjson::Value& object)
	{
		if (not object.HasMember("grid-options")) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Object doesn't have a grid-options key."
			);
		}
		const auto& grid_options = object["grid-options"];


		// grid periodicity
		const auto periodic_name = Periodic::get_option_name();
		if (not grid_options.HasMember(periodic_name.c_str())) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "grid-options doesn't have a " + periodic_name + " key."
			);
		}
		const auto& periodic_json = grid_options[periodic_name.c_str()];

		if (not periodic_json.IsString()) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Expression for grid periodicity isn't a string '"
				+ parser.GetExpr()
			);
		}

		this->parser.SetExpr(periodic_json.GetString());
		auto evaluated = this->parser.Eval().GetArray();
		if (evaluated.GetRows() != 1) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Invalid number of rows in expression '" + parser.GetExpr()
				+ "': " + std::to_string(evaluated.GetRows()) + ", should be 1"
			);
		}
		if (evaluated.GetCols() != 3) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Invalid number of columns in expression '" + parser.GetExpr()
				+ "': " + std::to_string(evaluated.GetRows()) + ", should be 3"
			);
		}

		this->periodic[0] = evaluated.At(0).GetBool();
		this->periodic[1] = evaluated.At(1).GetBool();
		this->periodic[2] = evaluated.At(2).GetBool();


		// number of grid cells
		const auto nr_cells_name = Number_Of_Cells::get_option_name();
		if (not grid_options.HasMember(nr_cells_name.c_str())) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "grid-options doesn't have a " + nr_cells_name + " key."
			);
		}
		const auto& nr_cells_json = grid_options[nr_cells_name.c_str()];

		if (not nr_cells_json.IsString()) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Expression for number of cells isn't a string."
			);
		}

		this->parser.SetExpr(nr_cells_json.GetString());
		evaluated = this->parser.Eval().GetArray();
		if (evaluated.GetRows() != 1) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Invalid number of rows in expression '" + parser.GetExpr()
				+ "': " + std::to_string(evaluated.GetRows()) + ", should be 1"
			);
		}
		if (evaluated.GetCols() != 3) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Invalid number of columns in expression '" + parser.GetExpr()
				+ "': " + std::to_string(evaluated.GetRows()) + ", should be 3"
			);
		}

		this->cells[0] = evaluated.At(0).GetFloat();
		this->cells[1] = evaluated.At(1).GetFloat();
		this->cells[2] = evaluated.At(2).GetFloat();


		// grid volume
		mup::Value cells_val(3, 0);
		cells_val.At(0) = mup::int_type(this->cells[0]);
		cells_val.At(1) = mup::int_type(this->cells[1]);
		cells_val.At(2) = mup::int_type(this->cells[2]);
		mup::Variable cells_var{&cells_val};
		this->parser.DefineVar("cells", cells_var);

		const auto volume_name = Volume::get_option_name();
		if (not grid_options.HasMember(volume_name.c_str())) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "grid-options doesn't have a " + volume_name + " key."
			);
		}
		const auto& volume_json = grid_options[volume_name.c_str()];

		if (not volume_json.IsString()) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Expression for grid volume isn't a string '"
				+ parser.GetExpr()
			);
		}

		this->parser.SetExpr(volume_json.GetString());
		evaluated = this->parser.Eval().GetArray();
		if (evaluated.GetRows() != 1) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Invalid number of rows in expression '" + parser.GetExpr()
				+ "': " + std::to_string(evaluated.GetRows()) + ", should be 1"
			);
		}
		if (evaluated.GetCols() != 3) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Invalid number of columns in expression '" + parser.GetExpr()
				+ "': " + std::to_string(evaluated.GetRows()) + ", should be 3"
			);
		}

		this->volume[0] = evaluated.At(0).GetFloat();
		this->volume[1] = evaluated.At(1).GetFloat();
		this->volume[2] = evaluated.At(2).GetFloat();


		// grid starting coordinate
		mup::Value volume_val(3, 0);
		volume_val.At(0) = this->volume[0];
		volume_val.At(1) = this->volume[1];
		volume_val.At(2) = this->volume[2];
		mup::Variable volume_var{&volume_val};
		this->parser.DefineVar("volume", volume_var);

		const auto start_name = Start::get_option_name();
		if (not grid_options.HasMember(start_name.c_str())) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "grid-options doesn't have a " + start_name + " key."
			);
		}
		const auto& start_json = grid_options[start_name.c_str()];

		if (not start_json.IsString()) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Expression for grid start isn't a string '"
				+ parser.GetExpr()
			);
		}

		this->parser.SetExpr(start_json.GetString());
		evaluated = this->parser.Eval().GetArray();
		if (evaluated.GetRows() != 1) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Invalid number of rows in expression '" + parser.GetExpr()
				+ "': " + std::to_string(evaluated.GetRows()) + ", should be 1"
			);
		}
		if (evaluated.GetCols() != 3) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "Invalid number of columns in expression '" + parser.GetExpr()
				+ "': " + std::to_string(evaluated.GetRows()) + ", should be 3"
			);
		}

		this->start[0] = evaluated.At(0).GetFloat();
		this->start[1] = evaluated.At(1).GetFloat();
		this->start[2] = evaluated.At(2).GetFloat();


		// maximum refinement level
		if (not object.HasMember("max-ref-lvl")) {
			return;
		}
		const auto& max_ref_lvl_json = object["max-ref-lvl"];
		if (not max_ref_lvl_json.IsInt()) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "JSON item max-ref-lvl must be integer."
			);
		}
		max_ref_lvl = object["max-ref-lvl"].GetInt();
		if (max_ref_lvl < 0) {
			throw std::invalid_argument(
				std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
				+ "JSON item max-ref-lvl cannot be negative."
			);
		}
	}


	const typename Periodic::data_type& get_periodic() const
	{
		return this->periodic;
	}

	const typename Number_Of_Cells::data_type& get_number_of_cells() const
	{
		return this->cells;
	}

	const typename Volume::data_type& get_volume() const
	{
		return this->volume;
	}

	const typename Start::data_type& get_start() const
	{
		return this->start;
	}

	const int& get_max_ref_lvl() const
	{
		return this->max_ref_lvl;
	}


private:

	typename Periodic::data_type periodic;
	typename Number_Of_Cells::data_type cells;
	typename Volume::data_type volume;
	typename Start::data_type start;
	int max_ref_lvl = 0;

	mup::ParserX parser = mup::ParserX(mup::pckCOMMON | mup::pckNON_COMPLEX | mup::pckMATRIX | mup::pckUNIT);
};

}} // namespaces

#endif // ifndef PAMHD_GRID_OPTIONS_HPP
