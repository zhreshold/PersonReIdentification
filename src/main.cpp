/***********************************************************************/
/*
/*   Script File: main.cpp
/*
/*   Description:
/*
/*   Main entrance
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

#include "define.h"
#include "feature.h"
#include "dataset.h"

// training 1 or testing 0
int	programPhase = 0;

void test()
{
	pri_dataset		dataset(PHASE_TESTING);
	pri_feat		feature;

	feature.init(dataset);

	cout << "Extracting image features..." << endl;
	feature.extract_feature();

	cout << "Matching query images in gallery..." << endl;
	feature.rank_cmc();

	// display some results
	cout << "Rank 1 : " << feature.get_rank_n(1) * 100 << endl;
	cout << "Rank 5 : " << feature.get_rank_n(5) * 100 << endl;
	cout << "Rank 10 : " << feature.get_rank_n(10) * 100 << endl;
	cout << "Rank 20 : " << feature.get_rank_n(20) * 100 << endl;
	cout << "Rank 50 : " << feature.get_rank_n(50) * 100 << endl;
	cout << "Rank 100 : " << feature.get_rank_n(100) * 100 << endl;
	cout << "Rank 316 : " << feature.get_rank_n(316) * 100 << endl;

}

void train()
{
	pri_dataset	dataset(PHASE_TRAINING);
	pri_feat	feature;

	feature.init(dataset);

	cout << "Extracting image features..." << endl;
	feature.extract_feature();

	cout << "Save pairwise image feature to file..." << endl;
	feature.save_pairwise_feature_image();

	cout << "Training image SVM model..." << endl;
	feature.train_image_model();
}


void main()
{
	// train new GMM
	if (NEW_GMM_MODELS)
	{
		pri_dataset	wholeDataset(PHASE_COLLECTING);
		pri_feat	featNewGMM(0);

		featNewGMM.init_new_gmm(wholeDataset);
	}

	if (programPhase == 1)
	{
		train();
	}
	else
	{
		test();
	}

	system("pause");
}