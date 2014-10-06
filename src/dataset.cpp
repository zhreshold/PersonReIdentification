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

using namespace std;
using namespace cv;


pri_dataset::pri_dataset()
{
	datasetType = DATASET_TYPE;
	m_phase = PHASE_TRAINING;

	if (datasetType == VIPER)
		init_viper();
}

pri_dataset::pri_dataset(const int phase)
{
	datasetType = DATASET_TYPE;
	m_phase = phase;

	if (datasetType == VIPER)
		init_viper();
}

pri_dataset::~pri_dataset()
{
	datasetType = NULL;
	m_phase = PHASE_NO_STATUS;
	numPersons = 0;
	numPersonTotal = 0;
	numShots = 0;
	filenames.clear();
	queryIdx.clear();
}

void pri_dataset::init_viper()
{
	// config
	numPersonTotal = 632;
	numShots = 2;
	gROI = Rect(4, 16, 40, 106);


	// path
	string		path = ROOT_PATH + string("data/VIPeR/");
	char		filename[256];
	int			index = 0;
	vector<int>	perm;
	//fstream		file;

	
	if (m_phase == PHASE_TRAINING)		
	{
		// init for training
		for (int i = 0; i < numPersonTotal; i++)
			perm.push_back(i);

		numPersons = NUM_PERSON_TRAIN;

		// get random permutation by shuffling
		srand(RANDOM_SEED);
		random_shuffle(perm.begin(), perm.end());

		for (int i = 0; i < numPersons; i++)
		{
			vector<int>		personQuery;

			// image 0 for person i
			sprintf_s(filename, "%03d_00.bmp", perm[i]);
			filenames.push_back(path + string(filename));
			personQuery.push_back(index);
			index++;

			// image 1 for person i
			sprintf_s(filename, "%03d_01.bmp", perm[i]);
			filenames.push_back(path + string(filename));
			personQuery.push_back(index);
			index++;

			// push back query indexes for person i, append person id
			personQuery.push_back(perm[i]);
			queryIdx.push_back(personQuery);
		}

		// write the rest to file for test
		ofstream file(ROOT_PATH + string("cache/testPartition.txt") , ios::out | ios::trunc);
		if (!file.is_open())
		{
			exit(ERR_FILE_UNABLE_TO_OPEN);
		}

		for (int i = numPersons; i < numPersonTotal; i++)
		{
			file << perm[i] << endl;
		}

		file.close();

	}
	else if (m_phase == PHASE_TESTING)
	{
		// init for testing
		numPersons = numPersonTotal - NUM_PERSON_TRAIN;

		// load partition for test from file
		ifstream file(ROOT_PATH + string("cache/testPartition.txt"), ios::in);
		if (!file.is_open())
		{
			exit(ERR_FILE_NOT_EXIST);
		}

		for (int i = 0; i < numPersons; i++)
		{
			int		id;
			file >> id;
			perm.push_back(id);
		}

		file.close();

		// load image names according to the partition
		for (int i = 0; i < numPersons; i++)
		{
			vector<int>		personQuery;

			// image 0 for person i
			sprintf_s(filename, "%03d_00.bmp", perm[i]);
			filenames.push_back(path + string(filename));
			personQuery.push_back(index);
			index++;

			// image 1 for person i
			sprintf_s(filename, "%03d_01.bmp", perm[i]);
			filenames.push_back(path + string(filename));
			personQuery.push_back(index);
			index++;

			// push back query indexes for person i, append person id
			personQuery.push_back(perm[i]);
			queryIdx.push_back(personQuery);
		}

	}
	else
	{
		// init for collecting all the feature points
		numPersons = numPersonTotal;

		// load image names
		for (int i = 0; i < numPersons; i++)
		{
			// image 0 for person i
			sprintf_s(filename, "%03d_00.bmp", i);
			filenames.push_back(path + string(filename));

			// image 1 for person i
			sprintf_s(filename, "%03d_01.bmp", i);
			filenames.push_back(path + string(filename));
		}
	}
	
	

}

int	pri_dataset::num_person()
{
	return numPersons;
}

int	pri_dataset::num_shot()
{
	return numShots;
}

vector<string> pri_dataset::get_filenames()
{
	return filenames;
}

vector<vector<int>> pri_dataset::get_query_index()
{
	return queryIdx;
}

Rect pri_dataset::get_roi()
{
	return gROI;
}