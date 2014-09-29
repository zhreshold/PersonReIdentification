/***********************************************************************/
/*
/*   Script File: hkmeans.h
/*
/*   Description:
/*
/*   Hierachical kmeans class declaration
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

#ifndef _HKMEANS_H_
#define _HKMEANS_H_

#include <vector>
#include <opencv\cv.h>

using namespace std;

// ! whether or not store patch data for illustration
#define	STORE_PATCH_DATA	1

typedef float	hk_type;
typedef unsigned char	uchar;

struct hkPatch
{
	int		flag;

	int		width;
	int		height;
	int		channel;
	int		step;

	// ! store patch data if turned on
	vector<uchar> data;

	// ! patch feature
	vector<hk_type> feature;
};

struct hkNode
{
	hkNode*	parent;
	vector<hk_type> center;
	vector<hkNode*> childs;
};



class hkmeans
{
public:
	hkmeans();
	~hkmeans();

	// ! partition K in each level
	int	K;
	// ! root node
	hkNode root;
	// ! training data
	vector<hkPatch>	patches;

	// static functions
	void	extractFeature(hkPatch &patch, cv::Mat input);

private:

	

};





#endif