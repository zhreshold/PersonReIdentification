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


int	main()
{
	// train new model
	if (NEW_GMM_MODELS)
	{
		pri_dataset	wholeDataset(PHASE_COLLECTING);
		pri_feat	featNewGMM(0);

		featNewGMM.init_new_gmm(wholeDataset);
	}
	

	pri_dataset		dataset(PHASE_TRAINING);
	pri_feat		feature;

	feature.init(dataset);
	feature.extract_feature();



	system("pause");
}