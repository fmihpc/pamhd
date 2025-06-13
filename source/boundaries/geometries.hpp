/*
Boundary geometries class of PAMHD

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


#ifndef PAMHD_GEOMETRIES_HPP
#define PAMHD_GEOMETRIES_HPP


#include "algorithm"
#include "cstdlib"
#include "exception"
#include "iostream"
#include "limits"
#include "map"
#include "string"
#include "utility"

#include "rapidjson/document.h"

#include "boundaries/box.hpp"
#include "boundaries/sphere.hpp"


namespace pamhd {
namespace boundaries {

/*!
Vector_T is assumed to be std::array or similar.
*/
template<
	class Geometry_Id,
	class Vector,
	class Scalar,
	class Cell_Id
> class Geometries
{
public:

	Geometries() = default;

	Geometries(
		const Geometries<Geometry_Id, Vector, Scalar, Cell_Id>&
	) = default;

	Geometries(
		Geometries<Geometry_Id, Vector, Scalar, Cell_Id>&&
	) = default;

	Geometries(const rapidjson::Value& object)
	{
		this->set(object);
	}


	/*!
	Creates geometries from given rapidjson object.

	Given object must be an array of objects each of which
	contains a supported geometry object.

	Example json object:
	\verbatim
	{
		"geometries": [
			{"box": {
				"start": [0, 0, 0],
				"end": [1, 1, 1]
			}},
			{"sphere": {
				"center": [-1, -2, -3],
				"radius": 0.01
			}}
		]
	}
	\endverbatim
	*/
	void set(const rapidjson::Value& object)
	{
		if (not object.HasMember("geometries")) {
			return;
		}
		const auto& json_geometries = object["geometries"];

		if (not json_geometries.IsArray()) {
			throw std::invalid_argument("geometries item is not an array (\"geometries\": [...]).");
		}

		for (rapidjson::SizeType i = 0; i < json_geometries.Size(); i++) {
			this->boxes.erase(i);
			this->spheres.erase(i);

			if (json_geometries[i].HasMember("box")) {
				this->boxes[i]
					= Box<Vector, Cell_Id>(
						json_geometries[i]["box"]
					);
			} else if (json_geometries[i].HasMember("sphere")) {
				this->spheres[i]
					= Sphere<Vector, Scalar, Cell_Id>(
						json_geometries[i]["sphere"]
					);
			} else {
				throw std::invalid_argument(
					__FILE__ ": Item in geometries array is not of a supported type."
				);
			}
		}
	}


	/*!
	Calls overlaps function of each existing geometry.

	Returns ids of all geometries that overlap given volume.
	*/
	std::vector<Geometry_Id> overlaps(
		const Vector& start,
		const Vector& end,
		const Cell_Id& cell_id
	) {
		std::vector<Geometry_Id> geometry_ids;
		geometry_ids.reserve(this->boxes.size() + this->spheres.size());
		for (const auto& box: this->boxes) {
			geometry_ids.push_back(box.first);
		}
		for (const auto& sphere: this->spheres) {
			geometry_ids.push_back(sphere.first);
		}

		std::vector<Geometry_Id> overlapping;
		if (geometry_ids.size() == 0) {
			return overlapping;
		}

		std::sort(geometry_ids.begin(), geometry_ids.end());

		for (const auto& geometry_id: geometry_ids) {

			if (this->boxes.count(geometry_id) > 0) {
				if (this->boxes.at(geometry_id).overlaps(start, end, cell_id)) {
					overlapping.push_back(geometry_id);
				}
			} else if (this->spheres.count(geometry_id) > 0) {
				if (this->spheres.at(geometry_id).overlaps(start, end, cell_id)) {
					overlapping.push_back(geometry_id);
				}
			} else {
				throw std::out_of_range(
					std::string(__FILE__ "(") + std::to_string(__LINE__) + "): "
					+ "Internal error"
				);
			}
		}

		return overlapping;
	}


	/*!
	Returns cells belonging to given geometry id.

	Throws an exception if given geometry id doesn't exist.
	*/
	const std::vector<Cell_Id>& get_cells(const Geometry_Id& geometry) const
	{
		if (this->boxes.count(geometry) > 0) {
			return this->boxes.at(geometry).cells;
		} else if (this->spheres.count(geometry) > 0) {
			return this->spheres.at(geometry).cells;
		}

		throw std::invalid_argument(__FILE__ ": Given geometry id doesn't exist.");
	}


	void clear_cells(const Geometry_Id& geometry) {
		if (this->boxes.count(geometry) > 0) {
			this->boxes.at(geometry).clear_cells();
		} else if (this->spheres.count(geometry) > 0) {
			this->spheres.at(geometry).clear_cells();
		}
	}


	std::vector<Geometry_Id> get_geometry_ids() const
	{
		std::vector<Geometry_Id> ids;
		ids.reserve(this->boxes.size() + this->spheres.size());

		for (const auto& box: this->boxes) {
			ids.push_back(box.first);
		}

		for (const auto& sphere: this->spheres) {
			ids.push_back(sphere.first);
		}

		return ids;
	}


private:

	std::map<
		Geometry_Id, // geometry id shared between all geometries
		Box<Vector, Cell_Id>
	> boxes;

	std::map<
		Geometry_Id,
		Sphere<Vector, Scalar, Cell_Id>
	> spheres;

};

}} // namespaces

#endif // ifndef PAMHD_GEOMETRIES_HPP
