/***********************************************************************/
/*
/*   Script File: hkmeans.cpp
/*
/*   Description:
/*
/*   Hierachical kmeans class
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

#include "hkmeans.h"

using namespace cv;

hkmeans::hkmeans()
{
}

hkmeans::~hkmeans()
{
}

void hkmeans::extractFeature(hkPatch &patch, Mat input)
{
	if (input.rows < 1 || input.cols < 1)
	{
		printf("Input patch is empty!\n");
		exit(-1);
	}

	



}