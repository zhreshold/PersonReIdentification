/***********************************************************************/
/*
/*   Script File: dataset.h
/*
/*   Description:
/*
/*   Dataset import class declaration
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


#ifndef _DATASET_H_
#define	_DATASET_H_

#include "define.h"
#include <string>
#include <vector>
#include <opencv\cv.h>


class pri_dataset
{
public:
	pri_dataset();
	pri_dataset(const int phase);
	~pri_dataset();

	int					num_person();					// return number of individuals
	int					num_shot();						// return number of shots of each individual
	std::vector<std::string>		get_filenames();				// return filename list
	std::vector<std::vector<int>> get_query_index();				// return query index
	cv::Rect				get_roi();						// return global ROI

private:
	int					m_phase;						// phase = 1 if training,  0 if testing
	int					datasetType;					// dataset = VIPER, CUHK01...
	int					numPersons;						// number of individuals for current experiment
	int					numShots;						// number of shots of each individual
	int					numPersonTotal;					// number of individuals in the entire dataset

	cv::Rect						gROI;							// global ROI for specific dataset
	std::vector<std::string>		filenames;						// image lists
	std::vector<std::vector<int>>	queryIdx;						// index of images used for experiment, DIM = (num_person) x ( num_shots)

	// private functions
	void init_viper();
};





#endif