/***********************************************************************/
/*
/*   Script File: feature.h
/*
/*   Description:
/*
/*   Feature class header
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

#ifndef _FEATURE_H_
#define	_FEATURE_H_

#include "define.h"
#include "dataset.h"
#include <vector>
#include <opencv\cv.h>
#include <opencv\highgui.h>
extern "C"
{
#include <vl/gmm.h>
#include <vl/fisher.h>
}

// color labels definition. 
// note: do not change the starting value 0 or add unsorted value in enum!!
#define	COLOR_LABEL_METHOD		1
#if COLOR_LABEL_METHOD == 0
enum COLOR_LABEL_ENUM
{
	COLOR_BLACK = 0,
	COLOR_DEEPBLUE,	
	COLOR_GRAY,		
	COLOR_GREEN,	
	COLOR_LIGHTBLUE,	
	COLOR_MIXED,		
	COLOR_ORANGE,	
	COLOR_PINK,		
	COLOR_PURPLE,	
	COLOR_RED,		
	COLOR_WHITE,		
	COLOR_YELLOW
	COLOR_LABEL_COUNT,
};	// COLOR_LABEL_METHOD == 0
#elif COLOR_LABEL_METHOD == 1
enum COLOR_LABEL_ENUM
{
	COLOR_BLACK = 0,
	COLOR_GRAY,
	COLOR_WHITEGRAY,
	COLOR_WHITE,
	COLOR_RED,
	COLOR_ORANGE,
	COLOR_YELLOW,
	COLOR_CHARTREUSE_GREEN,
	COLOR_GREEN,
	COLOR_SPRING_GREEN,
	COLOR_CYAN,
	COLOR_AZURE,
	COLOR_BLUE,
	COLOR_VIOLET,
	COLOR_MAGENTA,
	COLOR_ROSE,
	COLOR_LABEL_COUNT,
};
#endif // COLOR_LABEL_METHOD == 1

using namespace std;
using namespace cv;

class pri_feat
{
public:
	pri_feat();
	pri_feat(int indicator);
	~pri_feat();


	// local functions
	static void		get_hsv_color_labels(vector<UChar> &hsvColorLabels);
	static void		rgb2hsv(float r, float g, float b, float &h, float &s, float &v);
	static void		get_local_descriptor(int row, int col, vector<LDType> &ldBuffPixel, Mat hsvImage);
	

	// public functions
	void	init(pri_dataset &dataset);
	void	init_new_gmm(pri_dataset &dataset);
	void	extract_feature();

private:
	vector<vector<FeatureType>> imgFeat;					// image feature buffer
	vector<int>					imgFeatID;					// id vector of image feature
	vector<vector<FeatureType>> pairFeat;					// pairwise feature buffer
	vector<UChar>				m_hsvColorLabels;			// color labels look-up table in HSV space, 16x16x16 = 4096 dim
	vector<LDType>				ldBuffer;					// local descriptor buffer
	VlGMM*						m_gmm;						// gaussian mixture models 

	int							numPersons;					// number of individuals for current experiment
	int							numShots;					// number of shots of each individual
	vector<string>				filenames;					// image lists
	vector<vector<int>>			queryIdx;					// index of images used for experiment, DIM = (num_person) x ( num_shots)
	Rect						gROI;						// global ROI

	// image buffer
	Mat							image;						// original rgb image
	Mat							mask;						// image foreground mask
	Mat							hsvImage;					// hsv image
	Mat							labImage;					// Lab image
	Mat							grayImage;					// grayscale image


	// private functions
	void	extract_feature_image(vector<FeatureType> &feat);
	void	extract_feature_block(vector<FeatureType> &blockFeat, Rect blkROI);
	int		collect_local_descriptors();
	void	load_gmm();
	
};









#endif