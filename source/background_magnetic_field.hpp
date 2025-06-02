/*
Handles background magnetic field of MHD part of PAMHD.

Copyright 2014, 2016, 2017 Ilja Honkonen
Copyright 2025 Finnish Meteorological Institute
All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.


Author(s): Ilja Honkonen
*/

#ifndef PAMHD_BACKGROUND_MAGNETIC_FIELD_HPP
#define PAMHD_BACKGROUND_MAGNETIC_FIELD_HPP


#include "cmath"
#include "limits"
#include "string"

#include "rapidjson/document.h"

#include "common_functions.hpp"


namespace pamhd {


template<class Scalar> struct Line_Dipole
{
	/*
	Line dipole, dimension and direction start from 1.

	If dimension == direction then gives a dipole field
	in plane perpendicular to dimension when position along
	dimension == 0.

	Example:
		dimension = 2, direction = -2 gives:
			r4 = (x^2+z^2)
			Bx = magnitude * +2 * (z^2 - x^2) / r4
			Bz = magnitude * -2*x*z / r4
	*/
	std::array<Scalar, 2> position = {0, 0};
	unsigned int line_dimension = 0;
	int field_direction_through_line = 0;
	Scalar magnitude = 0;
};


template<class Scalar, class Vector> class Background_Magnetic_Field
{

private:
	Vector constant = {0, 0, 0};
	std::vector<std::pair<Vector, Vector>> dipole_moments_positions;
	std::vector<Line_Dipole<Scalar>> line_dipoles;
	Scalar min_distance = std::numeric_limits<Scalar>::epsilon();
	bool exists_ = false;


public:
	Background_Magnetic_Field() = default;
	Background_Magnetic_Field(const Background_Magnetic_Field& other) = default;
	Background_Magnetic_Field(Background_Magnetic_Field&& other) = default;
	Background_Magnetic_Field(const rapidjson::Value& object)
	{
		this->set(object);
	}


	bool exists() const {
		return this->exists_;
	}


	void set(const rapidjson::Value& object) {
		using std::abs;
		using std::isnormal;

		if (not object.HasMember("background-magnetic-field")) {
			return;
		}
		const auto& bg_B = object["background-magnetic-field"];

		if (bg_B.HasMember("value")) {
			const auto& value = bg_B["value"];
			if (not value.IsArray()) {
				throw std::invalid_argument(
					std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
					+ "Background magnetic field value is not an array."
				);
			}

			if (value.Size() != 3) {
				throw std::invalid_argument(
					std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
					+ "Background magnetic field value must have 3 components."
				);
			}

			this->constant[0] = value[0].GetDouble();
			this->constant[1] = value[1].GetDouble();
			this->constant[2] = value[2].GetDouble();

			if (
				this->constant[0] != 0
				or this->constant[1] != 0
				or this->constant[2] != 0
			) {
				this->exists_ = true;
			}
		}

		if (bg_B.HasMember("min-distance")) {
			const auto& min_distance_json = bg_B["min-distance"];
			if (not min_distance_json.IsNumber()) {
				throw std::invalid_argument(
					std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
					+ "Background magnetic field minimum distance is not a number."
				);
			}

			this->min_distance = min_distance_json.GetDouble();
		}

		if (bg_B.HasMember("line-dipoles")) {
			const auto& dipoles = bg_B["line-dipoles"];
			if (not dipoles.IsArray()) {
				throw std::invalid_argument(
					std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
					+ "Dipoles item must be an array."
				);
			}

			for (rapidjson::SizeType i = 0; i < dipoles.Size(); i++) {
				Line_Dipole<Scalar> dipole;

				if (not dipoles[i].HasMember("line-dimension")) {
					throw std::invalid_argument(
						std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
						+ "Line dipole missing line-dimension member."
					);
				}
				const auto& line_dim_json = dipoles[i]["line-dimension"];
				if (not line_dim_json.IsUint()) {
					throw std::invalid_argument(
						std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
						+ "line-dimension must be an unsigned integer."
					);
				}

				dipole.line_dimension = line_dim_json.GetUint();
				if (dipole.line_dimension == 0 or dipole.line_dimension > 3) {
					throw std::invalid_argument(
						std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
						+ "line-dimension must be > 0 and < 4."
					);
				}

				if (not dipoles[i].HasMember("field-direction-through-line")) {
					throw std::invalid_argument(
						std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
						+ "Line dipole missing field-direction-through-line member."
					);
				}
				const auto& field_dir_json = dipoles[i]["field-direction-through-line"];
				if (not field_dir_json.IsInt()) {
					throw std::invalid_argument(
						std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
						+ "field-direction-through-line must be an integer."
					);
				}

				dipole.field_direction_through_line = field_dir_json.GetInt();
				if (dipole.field_direction_through_line == 0 or abs(dipole.field_direction_through_line) > 3) {
					throw std::invalid_argument(
						std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
						+ "must be 0 < |field-direction-through-line| < 4."
					);
				}

				if (not dipoles[i].HasMember("magnitude")) {
					throw std::invalid_argument(
						std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
						+ "Line dipole missing magnitude member."
					);
				}
				const auto& magnitude_json = dipoles[i]["magnitude"];
				if (not magnitude_json.IsDouble()) {
					throw std::invalid_argument(
						std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
						+ "magnitude must be a floating point number."
					);
				}
				dipole.magnitude = magnitude_json.GetDouble();

				if (not dipoles[i].HasMember("position")) {
					throw std::invalid_argument(
						std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
						+ std::to_string(i + 1) + "th line dipole missing a position member."
					);
				}
				const auto& position_json = dipoles[i]["position"];
				if (not position_json.IsArray()) {
					throw std::invalid_argument(
						std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
						+ std::to_string(i + 1) + "th line dipole position must be an array."
					);
				}
				if (position_json.Size() != 3) {
					throw std::invalid_argument(
						std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
						+ std::to_string(i + 1) + "th line dipole position must have 3 components."
					);
				}

				switch (dipole.line_dimension) {
				case 1:
					dipole.position[0] = position_json[1].GetDouble();
					dipole.position[1] = position_json[2].GetDouble();
					break;
				case 2:
					dipole.position[0] = position_json[0].GetDouble();
					dipole.position[1] = position_json[2].GetDouble();
					break;
				case 3:
					dipole.position[0] = position_json[0].GetDouble();
					dipole.position[1] = position_json[1].GetDouble();
					break;
				default:
					std::cerr << __FILE__ "(" << __LINE__ << ")" << std::endl;
					abort();
				}

				this->line_dipoles.push_back(dipole);
			}
		}

		if (bg_B.HasMember("dipoles")) {
			const auto& dipoles = bg_B["dipoles"];
			if (not dipoles.IsArray()) {
				throw std::invalid_argument(
					std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
					+ "Dipoles item must be an array."
				);
			}

			for (rapidjson::SizeType i = 0; i < dipoles.Size(); i++) {
				if (not dipoles[i].HasMember("moment")) {
					throw std::invalid_argument(
						std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
						+ std::to_string(i + 1) + "th dipole moment missing a moment member."
					);
				}
				const auto& moment_json = dipoles[i]["moment"];
				if (not moment_json.IsArray()) {
					throw std::invalid_argument(
						std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
						+ std::to_string(i + 1) + "th dipole moment must be an array."
					);
				}
				if (moment_json.Size() != 3) {
					throw std::invalid_argument(
						std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
						+ std::to_string(i + 1) + "th dipole moment array must have 3 components."
					);
				}

				const Vector moment{
					moment_json[0].GetDouble(),
					moment_json[1].GetDouble(),
					moment_json[2].GetDouble(),
				};
				if (moment[0] != 0 or moment[1] != 0 or moment[2] != 0) {
					this->exists_ = true;
				}

				if (not dipoles[i].HasMember("position")) {
					throw std::invalid_argument(
						std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
						+ std::to_string(i + 1) + "th dipole moment missing a position member."
					);
				}

				const auto& position_json = dipoles[i]["position"];
				if (not position_json.IsArray()) {
					throw std::invalid_argument(
						std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
						+ std::to_string(i + 1) + "th dipole position must be an array."
					);
				}
				if (position_json.Size() != 3) {
					throw std::invalid_argument(
						std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
						+ std::to_string(i + 1) + "th dipole position array must have 3 components."
					);
				}

				const Vector position{
					position_json[0].GetDouble(),
					position_json[1].GetDouble(),
					position_json[2].GetDouble(),
				};

				this->dipole_moments_positions.emplace_back(moment, position);
			}
		}
	}


	/*
	Modified from https://github.com/nasailja/background_B/blob/master/source/dipole.hpp

	Returns magnetic field from given sources at given position.
	If given position is closer to a dipole than this->min_distance,
	that dipole contributes mu0/2/pi*moment/min_distance^3 instead.
	*/
	Vector get_background_field(
		const Vector& field_position,
		const Scalar& vacuum_permeability
	) const {
		using std::abs;
		using std::pow;

		Vector ret_val{this->constant};

		for (const auto& line_dipole: this->line_dipoles) {
			const Vector
				dip_pos = [&](){
					Vector pos{0, 0, 0};
					switch (line_dipole.line_dimension) {
					case 1:
						pos[0] = field_position[0];
						pos[1] = line_dipole.position[0];
						pos[2] = line_dipole.position[1];
						break;
					case 2:
						pos[0] = line_dipole.position[0];
						pos[1] = field_position[1];
						pos[2] = line_dipole.position[1];
						break;
					case 3:
						pos[0] = line_dipole.position[0];
						pos[1] = line_dipole.position[1];
						pos[2] = field_position[2];
						break;
					default:
						std::cerr << __FILE__ "(" << __LINE__ << ")" << std::endl;
						abort();
					}
					return pos;
				}(),
				r = pamhd::add(field_position, pamhd::neg(dip_pos));

			const Scalar r1 = pamhd::norm(r);
			if (r1 < this->min_distance) {
				if (line_dipole.line_dimension == abs(line_dipole.field_direction_through_line)) {
					ret_val[line_dipole.line_dimension - 1]
						+= vacuum_permeability * line_dipole.magnitude
						/ (2 * M_PI * pow(this->min_distance, 3));
				} else {
					if (line_dipole.field_direction_through_line > 0) {
						ret_val[line_dipole.field_direction_through_line - 1]
							+= vacuum_permeability * line_dipole.magnitude
							/ (2 * M_PI * pow(this->min_distance, 3));
					} else {
						ret_val[abs(line_dipole.field_direction_through_line) - 1]
							-= vacuum_permeability * line_dipole.magnitude
							/ (2 * M_PI * pow(this->min_distance, 3));
					}
				}
			} else {
				Vector value{
					line_dipole.magnitude * vacuum_permeability / (2 * M_PI * r1*r1*r1*r1),
					line_dipole.magnitude * vacuum_permeability / (2 * M_PI * r1*r1*r1*r1),
					line_dipole.magnitude * vacuum_permeability / (2 * M_PI * r1*r1*r1*r1)
				};
				size_t value_dim1 = 0, value_dim2 = 0, coord_dim1 = 0, coord_dim2 = 0;

				if (line_dipole.line_dimension == 1) {
					if (abs(line_dipole.field_direction_through_line) == 1) {

						value[1] = value[2] = 0;
						value[0] *= -r1 / 2;

					} else {

						value[0] = 0;
						coord_dim1 = 1;
						coord_dim2 = 2;

						if (abs(line_dipole.field_direction_through_line) == 2) {
							value_dim1 = 2;
							value_dim2 = 1;
						} else {
							value_dim1 = 1;
							value_dim2 = 2;
						}
					}
				} else if (line_dipole.line_dimension == 2) {
					if (abs(line_dipole.field_direction_through_line) == 2) {

						value[0] = value[2] = 0;
						value[1] *= -r1 / 2;

					} else {

						value[1] = 0;
						coord_dim1 = 0;
						coord_dim2 = 2;

						if (abs(line_dipole.field_direction_through_line) == 1) {
							value_dim1 = 2;
							value_dim2 = 0;
						} else {
							value_dim1 = 0;
							value_dim2 = 2;
						}
					}
				} else if (line_dipole.line_dimension == 3) {
					if (abs(line_dipole.field_direction_through_line) == 3) {

						value[0] = value[1] = 0;
						value[2] *= -r1 / 2;

					} else {

						value[2] = 0;
						coord_dim1 = 0;
						coord_dim2 = 1;

						if (abs(line_dipole.field_direction_through_line) == 1) {
							value_dim1 = 1;
							value_dim2 = 0;
						} else {
							value_dim1 = 0;
							value_dim2 = 1;
						}
					}
				}

				value[value_dim1] *= r[coord_dim1] * r[coord_dim2];
				value[value_dim2] *= (r[coord_dim2]*r[coord_dim2] - r[coord_dim1]*r[coord_dim1]);
				if (line_dipole.field_direction_through_line < 0) {
					value[0] *= -1;
					value[1] *= -1;
					value[2] *= -1;
				}

				ret_val[0] += value[0];
				ret_val[1] += value[1];
				ret_val[2] += value[2];
			}
		}

		for (const auto& dip_mom_pos: this->dipole_moments_positions) {

			const Vector r = pamhd::add(field_position, pamhd::neg(dip_mom_pos.second));
			const Scalar r1 = pamhd::norm(r);
			if (r1 < this->min_distance) {
				ret_val = pamhd::add(ret_val,
					pamhd::mul(
						dip_mom_pos.first,
						vacuum_permeability / (2 * M_PI * std::pow(this->min_distance, 3))));
			} else {
				const Vector
					r_unit = pamhd::mul(r, 1 / r1),
					projected_dip_mom = pamhd::mul(
						pamhd::dot(dip_mom_pos.first, r_unit),
						r_unit);
				ret_val = pamhd::add(ret_val,
					pamhd::mul(
						pamhd::add(
							pamhd::mul(3, projected_dip_mom),
							pamhd::neg(dip_mom_pos.first)),
						vacuum_permeability / (4 * M_PI) / (r1*r1*r1)));
			}
		}

		return ret_val;
	}

};

} // namespace


#endif // ifndef PAMHD_BACKGROUND_MAGNETIC_FIELD_HPP
