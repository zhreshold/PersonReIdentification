/***********************************************************************/
/*
/*   Script File: utility.cpp
/*
/*   Description:
/*
/*   Helpful utility functions
/*
/*
/*   Author: Joshua Zhang
/*   Date:  Sep-2014
/*
/*   Copyright(C) 2014  Joshua Zhang	 - All Rights Reserved.
/*
/*   This software is available for non-commercial use only. It must
/*   not be modified and distributed without prior permission of the
/*   author. The author is not responsible for implications from the
/*   use of this software.
/*
/***********************************************************************/

#include "utility.h"
#include <vector>
#include <algorithm>


void rgb2hsv(float r, float g, float b, float &h, float &s, float &v)
{
	// fast version converter

	float K = 0.f;

	if (g < b)
	{
		std::swap(g, b);
		K = -1.f;
	}

	if (r < g)
	{
		std::swap(r, g);
		K = -2.f / 6.f - K;
	}

	float chroma = r - std::min(g, b);
	h = fabs(K + (g - b) / (6.f * chroma + 1e-20f));
	s = chroma / (r + 1e-20f);
	v = r;
}

int get_combination(int m, int n)
{
	if (m < 0 || m > n)
	{
		return 0;
	}

	int iComb = 1;
	int i = 0;
	while (i < n)
	{
		++i;
		iComb *= m - i + 1;
		iComb /= i;
	}

	return iComb;
}
