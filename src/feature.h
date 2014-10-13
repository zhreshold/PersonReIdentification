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
#include "hkmeans.h"
#include <vector>
#include <opencv\cv.h>
#include <opencv\highgui.h>
extern "C"
{
#include <vl/gmm.h>
#include <vl/fisher.h>
#include <vl/hog.h>
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

struct sort_descend
{
	float	score;
	int		id;

	// overload comparison function to descend order
	bool operator < (sort_descend &s) const { return score > s.score; }
};

class pri_feat
{
public:
	pri_feat();
	pri_feat(int indicator);
	~pri_feat();

	// public functions
	void	init(pri_dataset &dataset);
	void	init_new_gmm(pri_dataset &dataset);
	void	extract_feature();
	void	save_pairwise_feature_block();						// save pair-wise feature to files
	void	save_pairwise_feature_image();						// save combined block feature as image feature
	void	train_block_models();
	void	load_block_weights();
	void	train_image_model();
	void	load_image_weights();
	void	rank_cmc();
	float	get_rank_n(int n);									// return the n-th rank
	void	partition_sort();
	void	init_new_kmeans(hkmeans &km, pri_dataset &dataset);


#if DEV_DEBUG
	// debug functions
	void	debug_show_top_n(int n);					// show top n match results
	void	debug_show_strip_score();					// show pairwise strip score and image for observation
	void	partition_sort_show();
#endif

private:
	vector<vector<FeatureType>> imgFeat;					// image feature buffer
	vector<int>					imgFeatID;					// id vector of image feature
	vector<vector<FeatureType>> pairFeat;					// pairwise feature buffer
	vector<UChar>				m_hsvColorLabels;			// color labels look-up table in HSV space, 16x16x16 = 4096 dim
	vector<LDType>				ldBuffer;					// local descriptor buffer
	VlGMM*						m_gmm;						// gaussian mixture models
	VlHog*						m_hog;						// hog
	vector<vector<float>>		blockWeights;				// block-wise SVM weights, DIM = (numBlocks x dim) or ( 1 x dim) if use unified model
	vector<float>				imageWeights;				// image-wise SVM weights
	hkmeans						m_km;

	int							numPersons;					// number of individuals for current experiment
	int							numShots;					// number of shots of each individual
	vector<string>				filenames;					// image lists
	vector<vector<int>>			queryIdx;					// index of images used for experiment, DIM = (num_person) x ( num_shots)
	Rect						gROI;						// global ROI
	vector<vector<int>>			pairIdxIntra;				// intra pairs index table, DIM = (numIntraPairs) x 2
	vector<vector<int>>			pairIdxInter;				// inter pairs index table

	// image buffer
	Mat							image;						// original rgb image
	Mat							mask;						// image foreground mask
	Mat							hsvImage;					// hsv image
	Mat							labImage;					// Lab image
	Mat							grayImage;					// grayscale image

	// sort results for matching
	vector<vector<sort_descend>>	results;					// stores sorted results for matching
	vector<float>					ranks;						// store ranks;
	vector<vector<sort_descend>>	partitionSort;				// sort for each partition


	// private functions
	void	get_hsv_color_labels();					// compute HSV color label look-up table
	void	extract_feature_image(vector<FeatureType> &feat);						// image-wise feature
	void	extract_feature_block(vector<FeatureType> &blockFeat, Rect blkROI);		// block-wise feature
	int		collect_local_descriptors();											// collect local descriptor for GMMs
	void	load_gmm();																// load trained GMMs
	void	create_pair_index();													// create pairs index table
	void	get_local_descriptor(int row, int col, vector<LDType> &ldBuffPixel);	// collect local descriptors

	// similarity functions
	float	similarity_score(float f1, float f2);
	float	similarity_score_2(float f1, float f2);
	float	dist_score(float f1, float f2);
	float	hist_similartity_score(vector<FeatureType> f1, vector<FeatureType> f2);
	float	hist_similartity_score_2(vector<FeatureType> f1, vector<FeatureType> f2);
	float	hist_dist_score(vector<FeatureType> f1, vector<FeatureType> f2);
	void	write_similarity_to_file(ofstream &file, vector<FeatureType> f1, vector<FeatureType> f2, int featLen, int k);	// write similarity to file
	void	extract_feature_image_kmeans(hkmeans &km);
	void	get_combine_image_feature(vector<FeatureType> & combFeat, vector<FeatureType> f1, vector<FeatureType> f2);
	void    get_combine_image_feature_no_weight(vector<FeatureType> & combFeat, vector<FeatureType> f1, vector<FeatureType> f2);
	float	image_pairwise_score(int idx1, int idx2);
	
	
	
};









#endif