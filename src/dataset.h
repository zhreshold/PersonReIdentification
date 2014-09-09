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

using namespace std;


class pri_dataset
{
public:
	pri_dataset();
	pri_dataset(int	phase);
	~pri_dataset();

	int				num_person();
	int				num_shot();

private:
	int				m_phase;						// phase = 1 if training,  0 if testing
	int				dataset_type;					// dataset = VIPER, CUHK01...
	int				num_persons;					// number of individuals for current experiment
	int				num_shots;						// number of shots of each individual
	int				num_person_total;				// number of individuals in the entire dataset
	vector<string>	filenames;						// image lists
	vector< vector<int> > query_idx;				// index of images used for experiment, DIM = (num_person) x ( num_shots)

	// private functions
	void init_viper();
};





#endif