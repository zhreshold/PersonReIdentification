/***********************************************************************/
/*
/*   Script File: train.cpp
/*
/*   Description:
/*
/*   Main for training module
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
#include "dataset.h"
#include "feature.h"
#include <iostream>


int	main()
{
	// train new model
	if (NEW_GMM_MODELS)
	{
		pri_dataset	wholeDataset(PHASE_COLLECTING);
		pri_feat	featNewGMM(0);

		featNewGMM.init_new_gmm(wholeDataset);
	}
	//system("pause");
	//return(0);

	pri_dataset		dataset(PHASE_TRAINING);
	pri_feat		feature;


	feature.init(dataset);

	cout << "Extracting image features..." << endl;
	feature.extract_feature();

	cout << "Save pairwise block feature to files..." << endl;
	feature.save_pairwise_feature_block();

	cout << "Training block-wise SVM models..." << endl;
	feature.train_block_models();

	

	//cout << "Save pairwise image feature to file..." << endl;
	//feature.save_pairwise_feature_image();

	//cout << "Training image SVM model..." << endl;
	//feature.train_image_model();



	system("pause");
}