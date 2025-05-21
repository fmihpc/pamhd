/*
Tests mul(a, b).

Copyright 2025 Finnish Meteorological Institute
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


Author(s): Ilja Honkonen
*/


#include "array"
#include "deque"
#include "cstdlib"
#include "iostream"
#include "vector"

#include "prettyprint.hpp"

#include "common_functions.hpp"


int main()
{
	using std::array;
	using std::cerr;
	using std::deque;
	using std::endl;
	using std::vector;

	{
	const double v1{2}, v2{4};
	const auto result = pamhd::mul(v1, v2);
	if (result != 8) { cerr << __FILE__"("<<__LINE__<<") "
		" Incorrect result: " << result << endl;}
	}

	{
	const int v1{-3}, v2{2};
	const auto result = pamhd::mul(v1, v2);
	if (result != -6) { cerr << __FILE__"("<<__LINE__<<") "
		" Incorrect result: " << result << endl;}
	}

	{
	const array<float, 4> v1{1,2,3,4}, v2{-5,-6,-7,-8};
	const auto result = pamhd::mul(v1, v2);
	if (result != array<float, 4>{-5,-12,-21,-32}) { cerr << __FILE__"("<<__LINE__<<") "
		" Incorrect result: " << result << endl;}
	}

	{
	const array<int, 2> v1{-1,2};
	const int v2{4};
	const auto result = pamhd::mul(v1, v2);
	if (result != array<int, 2>{-4,8}) { cerr << __FILE__"("<<__LINE__<<") "
		" Incorrect result: " << result << endl;}
	}

	{
	const int v1{-4};
	const array<int, 1> v2{1};
	const auto result = pamhd::mul(v1, v2);
	if (result != array<int, 1>{-4}) { cerr << __FILE__"("<<__LINE__<<") "
		" Incorrect result: " << result << endl;}
	}

	{
	const deque<int> v1{9,8}, v2{-4,-3};
	const auto result = pamhd::mul(v1, v2);
	if (result != deque<int>{-36,-24}) { cerr << __FILE__"("<<__LINE__<<") "
		" Incorrect result: " << result << endl;}
	}

	{
	const vector<int> v1{3,2,1}, v2{2,3,4};
	const auto result = pamhd::mul(v1, v2);
	if (result != vector<int>{6,6,4}) { cerr<<__FILE__"("<<__LINE__<<") "
		" Incorrect result: " << result << endl;}
	}

	return EXIT_SUCCESS;
}
