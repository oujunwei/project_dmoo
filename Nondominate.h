/*!
Copyright (c) 2013, 申瑞珉 (Alden Ruimin Shen)
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#pragma once

#include <cassert>
#include <list>

namespace otl
{
namespace utility
{
template <typename _TIterator, typename _TIndividual, typename _TDominate>
bool IdentifyElite(_TIterator individual, std::list<_TIndividual> &population, std::list<_TIndividual> &nondominate, _TDominate dominate)
{
	for (_TIterator elite = nondominate.begin(); elite != nondominate.end();)
	{
		if (dominate(*individual, *elite))
		{
			typename std::list<_TIndividual>::iterator move = elite;
			++elite;
			population.splice(population.begin(), nondominate, move);
		}
		else if (dominate(*elite, *individual))
			return false;
		else
			++elite;
	}
	return true;
}

template <typename _TIndividual, typename _TDominate>
std::list<_TIndividual> ExtractNondominate(std::list<_TIndividual> &population, _TDominate dominate)
{
	typedef typename std::list<_TIndividual>::iterator _TIterator;
	assert(!population.empty());
	std::list<_TIndividual> nondominate;
	//c1.splice(c1.beg,c2,c2.beg)    将c2的beg位置的元素连接到c1的beg位置，并且在c2中施放掉beg位置的元素
	nondominate.splice(nondominate.end(), population, population.begin());
	for (_TIterator individual = population.begin(); individual != population.end();)
	{
		assert(!nondominate.empty());
		if (IdentifyElite(individual, population, nondominate, dominate))
		{
			typename std::list<_TIndividual>::iterator move = individual;
			++individual;
			nondominate.splice(nondominate.begin(), population, move);
		}
		else
			++individual;
	}
	return nondominate;
}
}
}
