/***********************************************************************/
/*
/*   Script File: test.cpp
/*
/*   Description:
/*
/*   Main for testing module
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
	pri_dataset		dataset(PHASE_TESTING);
	pri_feat		feature;

	feature.init(dataset);

	cout << "Extracting image features..." << endl;
	feature.extract_feature();





	//cout << "Saving partition sorts..." << endl;
	//feature.partition_sort();

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

#if DEV_DEBUG
	cout << "Debug partition sort..." << endl;
	feature.partition_sort_show();

	//cout << "Debug strip scores..." << endl;
	//feature.debug_show_strip_score();

	//feature.debug_show_top_n(10);
#endif

	system("pause");
}