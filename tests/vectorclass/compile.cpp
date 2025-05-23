/*
Tests compilation with vectorclass (www.agner.org/optimize/#vectorclass).

Copyright 2014, 2015, 2016, 2017 Ilja Honkonen
Copyright 2025 Finnish Meteorological Institute
All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.


Author(s): Ilja Honkonen
*/

#include "array"
#include "cstdlib"

#include "vectorclass.h"

int main()
{
	const Vec4i
		a(1, 2, 3, 4),
		b(-4, -3, -2, -1),
		c = a + b;

	const std::array<Vec4i, 2>
		e{{
			{1, 2, 3, 4},
			{-4, -3, -2, -1}
		}},
		f{{
			e[0] + e[1],
			e[0] * e[1]
		}};

	if (
		c[0] == -3
		and f[0][0] == -3
		and f[1][0] == -4
	) {
		return EXIT_SUCCESS;
	} else {
		return EXIT_FAILURE;
	}
}
