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
	{
		init_viper();
	}
	else
	{

	}

	obtain_rand_blocks();
}

pri_dataset::pri_dataset(const int phase)
{
	datasetType = DATASET_TYPE;
	m_phase = phase;

	if (datasetType == VIPER)
	{
		init_viper();
	}
	else
	{

	}

	obtain_rand_blocks();
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

vector<Rect> pri_dataset::get_rand_blocks()
{
	return randBlocks;
}

void pri_dataset::obtain_rand_blocks()
{
	randBlocks.clear();

	if (NEW_RAND_BLOCKS)
	{
		srand(time(NULL));
		while (randBlocks.size() < MAX_NUM_RAND_BLOCKS)
		{
			const int bw = MAX_RAND_BLOCK_WIDTH - MIN_RAND_BLOCK_WIDTH + 1;
			const int bh = MAX_RAND_BLOCK_HEIGHT - MIN_RAND_BLOCK_HEIGHT + 1;
			int x = rand() % gROI.width + gROI.x;
			int y = rand() % gROI.height + gROI.y;
			int	width = rand() % bw + MIN_RAND_BLOCK_WIDTH;
			int height = rand() % bh + MIN_RAND_BLOCK_HEIGHT;

			if ((x + width < gROI.x + gROI.width) && (y + height < gROI.y + gROI.height))
			{
				randBlocks.push_back(Rect(x, y, width, height));
			}
		}

		// erase duplicates if any
		for (int i = 0; i < randBlocks.size(); i++)
		{
			for (int j = i + 1; j < randBlocks.size();)
			{
				if (randBlocks[i] == randBlocks[j])
				{
					randBlocks.erase(randBlocks.begin() + j);
				}
				else
				{
					j++;
				}
			}
		}

		//save to file
		string filename;
		if (datasetType == VIPER)
		{
			filename = ROOT_PATH + string("cache/rand_blocks_viper.txt");
		}
		else
		{
			filename = ROOT_PATH + string("cache/rand_blocks_other.txt");
		}
		ofstream fp;
		fp.open(filename.c_str(), ios::out | ios::trunc);
		if (!fp.is_open())
		{
			printf("Error: can't open file to write blocks!\n");
			getchar();
			exit(-1);
		}

		fp << randBlocks.size() << endl;
		for (int i = 0; i < randBlocks.size(); i++)
		{
			fp << randBlocks[i].x << " ";
			fp << randBlocks[i].y << " ";
			fp << randBlocks[i].width << " ";
			fp << randBlocks[i].height << endl;
		}

		fp.close();
	}
	else
	{
		//load from file
		string filename;
		if (datasetType == VIPER)
		{
			filename = ROOT_PATH + string("cache/rand_blocks_viper.txt");
		}
		else
		{
			filename = ROOT_PATH + string("cache/rand_blocks_other.txt");
		}
		ifstream fp;
		fp.open(filename.c_str(), ios::in);
		if (!fp.is_open())
		{
			printf("Error: can't open file to read blocks!\n");
			getchar();
			exit(-1);
		}

		int size, x, y, width, height;

		fp >> size;

		if (size > MAX_NUM_RAND_BLOCKS)
		{
			size = MAX_NUM_RAND_BLOCKS;
		}

		for (int i = 0; i < size; i++)
		{
			fp >> x;
			fp >> y;
			fp >> width;
			fp >> height;
			randBlocks.push_back(Rect(x, y, width, height));
		}

		fp.close();

	}
}