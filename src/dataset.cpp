/***********************************************************************/
/*
/*   Script File: dataset.cpp
/*
/*   Description:
/*
/*   Dataset import class
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


#include "dataset.h"
#include <fstream>
#include <algorithm>


pri_dataset::pri_dataset()
{
	dataset_type = DATASET_TYPE;
	m_phase = PHASE_TRAINING;

	if (dataset_type == VIPER)
		init_viper();
}

pri_dataset::pri_dataset(int phase)
{
	dataset_type = DATASET_TYPE;
	m_phase = phase;

	if (dataset_type == VIPER)
		init_viper();
}

pri_dataset::~pri_dataset()
{
	dataset_type = NULL;
	m_phase = PHASE_NO_STATUS;
	num_persons = 0;
	num_person_total = 0;
	num_shots = 0;
	filenames.clear();
	query_idx.clear();
}

void pri_dataset::init_viper()
{
	// config
	num_person_total = 632;
	num_shots = 2;


	// path
	string		path = ROOT_PATH + string("data/VIPeR/");
	char		filename[256];
	int			index = 0;
	vector<int>	perm;
	//fstream		file;

	
	if (m_phase == PHASE_TRAINING)		
	{
		// init for training
		for (int i = 0; i < num_person_total; i++)
			perm.push_back(i);

		num_persons = NUM_PERSON_TRAIN;

		// get random permutation by shuffling
		random_shuffle(perm.begin(), perm.end());

		for (int i = 0; i < num_persons; i++)
		{
			vector<int>		personQuery;

			// image 0 for person i
			sprintf_s(filename, "%03d_00.bmp", perm[i]);
			filenames.push_back(filename);
			personQuery.push_back(index);
			index++;

			// image 1 for person i
			sprintf_s(filename, "%03d_01.bmp", perm[i]);
			filenames.push_back(filename);
			personQuery.push_back(index);
			index++;

			// push back query indexes for person i
			query_idx.push_back(personQuery);
		}

		// write the rest to file for test
		ofstream file(ROOT_PATH + string("cache/testPartition.txt") , ios::out | ios::trunc);
		if (!file.is_open())
		{
			exit(ERR_FILE_UNABLE_TO_OPEN);
		}

		for (int i = num_persons; i < num_person_total; i++)
		{
			file << perm[i] << endl;
		}

		file.close();

	}
	else
	{
		// init for testing
		num_persons = num_person_total - NUM_PERSON_TRAIN;

		// load partition for test from file
		ifstream file(ROOT_PATH + string("cache/testPartition.txt"), ios::in);
		if (!file.is_open())
		{
			exit(ERR_FILE_NOT_EXIST);
		}

		for (int i = 0; i < num_persons; i++)
		{
			int		id;
			file >> id;
			perm.push_back(id);
		}

		file.close();

		// load image names according to the partition
		for (int i = 0; i < num_persons; i++)
		{
			vector<int>		personQuery;

			// image 0 for person i
			sprintf_s(filename, "%03d_00.bmp", perm[i]);
			filenames.push_back(filename);
			personQuery.push_back(index);
			index++;

			// image 1 for person i
			sprintf_s(filename, "%03d_01.bmp", perm[i]);
			filenames.push_back(filename);
			personQuery.push_back(index);
			index++;

			// push back query indexes for person i
			query_idx.push_back(personQuery);
		}

	}
	
	

}

int	pri_dataset::num_person()
{
	return num_persons;
}

int	pri_dataset::num_shot()
{
	return num_shots;
}