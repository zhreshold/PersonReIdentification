/***********************************************************************/
/*
/*   Script File: feature.cpp
/*
/*   Description:
/*
/*   Feature class
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

#include "feature.h"
#include "utility.h"
#include <fstream>
#include <numeric>


pri_feat::pri_feat()
{
	// init
	numPersons = 0;
	numShots = 0;

	// compute hsv color label look-up table
	get_hsv_color_labels();

	// load GMM
	load_gmm();
}

pri_feat::pri_feat(int	indicator)
{
	if (indicator == 1)
	{
		// normal init
		numPersons = 0;
		numShots = 0;

		// compute hsv color label look-up table
		get_hsv_color_labels();

		// load GMM
		load_gmm();
	}
	else
	{
		// init
		numPersons = 0;
		numShots = 0;
	}

}

pri_feat::~pri_feat()
{
	numPersons = NULL;
	numShots = NULL;
	imgFeat.clear();
	imgFeatID.clear();
	pairFeat.clear();
	m_hsvColorLabels.clear();
	ldBuffer.clear();
	filenames.clear();
	queryIdx.clear();
	vl_gmm_delete(m_gmm);
	m_gmm = NULL;
	blockWeights.clear();
	imageWeights.clear();
	results.clear();
}

void pri_feat::init(pri_dataset &dataset)
{
	// copy infos from dataset class
	numPersons = dataset.num_person();
	numShots = dataset.num_shot();
	filenames = dataset.get_filenames();
	queryIdx = dataset.get_query_index();
	gROI = dataset.get_roi();

	// validation check
	if (queryIdx.size() != numPersons)
		exit(ERR_SIZE_MISMATCH);

	if (numPersons < 1)
		exit(ERR_EMPTY_SET);

	if (queryIdx[0].size() != numShots+1)
		exit(ERR_SIZE_MISMATCH);

	printf("Current # individuals: %d\n", numPersons);
	printf("Current # image per individual: %d\n", numShots);

}

void pri_feat::init_new_gmm(pri_dataset &dataset)
{
	// copy infos from dataset
	numPersons = dataset.num_person();
	numShots = dataset.num_shot();
	filenames = dataset.get_filenames();

	// validation check
	if (numPersons < 1)
		exit(ERR_EMPTY_SET);
	printf("Current # individuals: %d\n", numPersons);
	printf("Current # image per individual: %d\n", numShots);


	// collect local descriptors and train GMM
	ofstream	file(ROOT_PATH + string("cache/GMM.dat"), ios::out|ios::trunc);
	if (!file.is_open())
		exit(ERR_FILE_UNABLE_TO_OPEN);

	// get local descriptors
	printf("Extracting local descriptors...\n");
	for (int i = 0; i < filenames.size(); i++)
	{
		image = imread(filenames[i], 1);
		int status = collect_local_descriptors();

		if (status == 1)
			break;
	}

	printf("Total # pixels extracted: %d\n", ldBuffer.size() / LD_DIM);

	// train GMM
	printf("Training GMM, be patient...\n");
	m_gmm = vl_gmm_new(VL_TYPE_FLOAT, LD_DIM, LD_NUM_CLUSTERS);
	vl_gmm_cluster(m_gmm, ldBuffer.data(), ldBuffer.size()/LD_DIM);

	// save GMM
	printf("Saving GMM to file...\n");
	int		gSize = LD_DIM * LD_NUM_CLUSTERS;
	LDType*	gBuffer;

	// means
	gBuffer = (LDType*)vl_gmm_get_means(m_gmm);
	for (int i = 0; i < gSize; i++)
	{
		file << gBuffer[i] << "\t";
	}
	file << endl;

	// covariances
	gBuffer = (LDType*)vl_gmm_get_covariances(m_gmm);
	for (int i = 0; i < gSize; i++)
	{
		file << gBuffer[i] << "\t";
	}
	file << endl;

	// priors
	gBuffer = (LDType*)vl_gmm_get_priors(m_gmm);
	for (int i = 0; i < gSize; i++)
	{
		file << gBuffer[i] << "\t";
	}
	file << endl;


	file.close();

}

void pri_feat::load_gmm()
{
	ifstream	file(ROOT_PATH + string("cache/GMM.dat"), ios::in);
	if (!file.is_open())
		exit(ERR_FILE_NOT_EXIST);

	vector<LDType>	gBuffer;
	int				gSize = LD_DIM * LD_NUM_CLUSTERS;
	gBuffer.resize(gSize);
	//fill(gBuffer.begin(), gBuffer.end(), 0);

	m_gmm = vl_gmm_new(VL_TYPE_FLOAT, LD_DIM, LD_NUM_CLUSTERS);

	// load with means, covariances and priors
	for (int j = 0; j < gSize; j++)
	{
		file >> gBuffer[j];
	}
	vl_gmm_set_means(m_gmm, gBuffer.data());
	for (int j = 0; j < gSize; j++)
	{
		file >> gBuffer[j];
	}
	vl_gmm_set_covariances(m_gmm, gBuffer.data());
	for (int j = 0; j < gSize; j++)
	{
		file >> gBuffer[j];
	}
	vl_gmm_set_priors(m_gmm, gBuffer.data());
	

	file.close();
}


void pri_feat::get_hsv_color_labels()
{
	int			vr, vg, vb;
	float		r, g, b, h, s, v;
	int			CL;		

	// init look-up  table
	m_hsvColorLabels.clear();
	m_hsvColorLabels.resize(4096);
	fill(m_hsvColorLabels.begin(), m_hsvColorLabels.end(), 0);

	for (r = 0; r < 256; r += 16)
	for (g = 0; g < 256; g += 16)
	for (b = 0; b < 256; b += 16)
	{
		rgb2hsv(r, g, b, h, s, v);
		v /= 256.0;

#if (COLOR_LABEL_METHOD == 0)
		{
			if(v < 0.45)
			{
			CL = COLOR_BLACK;
			}
			else
			{
			if(s < 0.25)
			{
				//no color, pale, do Black, white and gray
				if(v < 0.5) CL = COLOR_BLACK;
				if((v>=0.5) && (v<0.9)) CL = COLOR_GRAY;
				if(v>=0.9) CL = COLOR_WHITE;
			}
			else
			{
				//it has colors
				if(h < 0.05) CL = COLOR_RED;
				if((h>=0.05) && (h<0.13)) CL = COLOR_ORANGE;
				if((h>=0.13) && (h<0.3)) CL = COLOR_YELLOW;
				if((h>=0.3) && (h<0.5)) CL = COLOR_GREEN;
				if((h>=0.5) && (h<0.7)) 
				{
					if(s < 0.8) CL = COLOR_LIGHTBLUE;
					else CL = COLOR_DEEPBLUE;
				}
				if((h>=0.7) && (h<0.8)) CL = COLOR_PURPLE;
				if((h>=0.8) && (h<1.0)) CL = COLOR_RED;
			}
		}
#elif (COLOR_LABEL_METHOD == 1)
		//new color label definitions, see http://en.wikipedia.org/wiki/Wikipedia:WikiProject_Color/Normalized_Color_Coordinates
		if ( v < 0.3)
		{
			CL = COLOR_BLACK;
		}
		else
		{
			if (s < 0.2)
			{
				//pale colors
				if (v < 0.3) CL = COLOR_BLACK;
				if ((v >= 0.3) && (v<0.55)) CL = COLOR_GRAY;
				if ((v >= 0.55) && (v < 0.8)) CL = COLOR_WHITEGRAY;
				if (v >= 0.8) CL = COLOR_WHITE;
			}
			else if (1)
			{
				if (h < 0.041) CL = COLOR_RED;
				if (h >= 0.041 && h < 0.125) CL = COLOR_ORANGE;
				if (h >= 0.125 && h < 0.208) CL = COLOR_YELLOW;
				if (h >= 0.208 && h < 0.292) CL = COLOR_CHARTREUSE_GREEN;
				if (h >= 0.292 && h < 0.375) CL = COLOR_GREEN;
				if (h >= 0.375 && h < 0.458) CL = COLOR_GREEN;
				if (h >= 0.458 && h < 0.542) CL = COLOR_CYAN;
				if (h >= 0.542 && h < 0.625) CL = COLOR_BLUE;
				if (h >= 0.625 && h < 0.708) CL = COLOR_BLUE;
				if (h >= 0.708 && h < 0.792) CL = COLOR_VIOLET;
				if (h >= 0.792 && h < 0.875) CL = COLOR_MAGENTA;
				if (h >= 0.875 && h < 0.958) CL = COLOR_ROSE;
				if (h >= 0.958) CL = COLOR_RED;
			}
			else
			{
				//brief version
				if (h < 0.083) CL = COLOR_RED;
				if (h >= 0.083 && h < 0.25) CL = COLOR_YELLOW;
				if (h >= 0.25 && h < 0.417) CL = COLOR_GREEN;
				if (h >= 0.417 && h < 0.583) CL = COLOR_CYAN;
				if (h >= 0.583 && h < 0.75) CL = COLOR_BLUE;
				if (h >= 0.75 && h < 0.917) CL = COLOR_MAGENTA;
				if (h >= 0.917) CL = COLOR_RED;
			}
		}

#endif

		// fill in the table
		vr = (int)r / 16;
		vg = (int)g / 16;
		vb = (int)b / 16;
		m_hsvColorLabels[vr * 256 + vg * 16 + vb] = CL;

	}
}



void pri_feat::create_pair_index()
{
	int		numIntraPerPerson = get_combination(2, numShots);
	int		numInterPerPerson = numIntraPerPerson * SVM_NEG_POS_RATIO;
	if (numInterPerPerson > numPersons - 1)
		numInterPerPerson = numPersons - 1;

	// clean up
	pairIdxIntra.clear();
	pairIdxInter.clear();

	// intra pairs first
	for (int i = 0; i < numPersons; i++)
	{
		for (int p1 = 0; p1 < numShots; p1++)
		for (int p2 = p1 + 1; p2 < numShots; p2++)
		{
			int	idx1 = queryIdx[i][p1];
			int	idx2 = queryIdx[i][p2];
			vector<int> idxIntra;
			idxIntra.push_back(idx1);
			idxIntra.push_back(idx2);
			pairIdxIntra.push_back(idxIntra);
		}
	}

	// inter pairs
	vector<int>		permIn, permOut;
	for (int i = 0; i < numPersons; i++)
	{
		permOut.push_back(i);
	}
	for (int i = 0; i < numShots; i++)
	{
		permIn.push_back(i);
	}

	srand(RANDOM_SEED);
	for (int i = 0; i < numPersons; i++)
	{
		random_shuffle(permOut.begin(), permOut.end());

		vector<int>	tmpPermOut = permOut;
		// excluse self
		for (vector<int>::iterator j = tmpPermOut.begin(); j != tmpPermOut.end(); j++)
		{
			if (*j == i)
			{
				tmpPermOut.erase(j);
				break;
			}
		}

		for (int k = 0; k < numInterPerPerson; k++)
		{

			// pick up randomly inside class
			random_shuffle(permIn.begin(), permIn.end());
			int		idx1 = queryIdx[i][permIn[0]];

			// shuffle again
			random_shuffle(permIn.begin(), permIn.end());
			int		idx2 = queryIdx[tmpPermOut[k]][permIn[0]];

			// push back
			vector<int> idxInter;
			idxInter.push_back(idx1);
			idxInter.push_back(idx2);
			pairIdxInter.push_back(idxInter);
		}

	}
}



void pri_feat::extract_feature()
{								
	int					numel = numPersons * numShots;			// number of images to process
	vector<FeatureType>	tmpFeat;								// temporal feature buffer

	// reset
	imgFeat.clear();
	imgFeatID.clear();

	for (int i = 0; i < numPersons; i++)
	{
		vector<int> query_person = queryIdx[i];
		for (int j = 0; j < numShots; j++)
		{
			// retrieve the index in filelist
			int		idx = query_person[j];
			// load image
			image = imread(filenames[idx], 1);
			// use default mask in this version
			mask.create(image.rows, image.cols, CV_8UC1);
			mask.setTo(1);

			// extract feature
			extract_feature_image(tmpFeat);
			// push back feature and id
			imgFeat.push_back(tmpFeat);
			imgFeatID.push_back(query_person[numShots]);

		}
	}
	
}


void pri_feat::extract_feature_image(vector<FeatureType> &feat)
{
	
	// copy from global ROI
	Rect				roi = gROI;
	int					blockHeight, blockWidth;
	vector<FeatureType>	blockFeatBuffer, imageFeatBuffer;

	// calculate block width and height
	//blockHeight = roi.height / IMAGE_PARTITION_Y;
	//blockWidth = roi.width / IMAGE_PARTITION_X;
	double	divy = IMAGE_PARTITION_Y - IMAGE_PARTITION_Y * PARTITION_OVERLAP + PARTITION_OVERLAP;
	double	divx = IMAGE_PARTITION_X - IMAGE_PARTITION_X * PARTITION_OVERLAP + PARTITION_OVERLAP;
	blockHeight = (int)floor(roi.height / divy);
	blockWidth = (int)floor(roi.width / divx);
	int	stepHeight = blockHeight - (int)ceil(blockHeight * PARTITION_OVERLAP);
	int stepWidth = blockWidth - (int)ceil(blockWidth * PARTITION_OVERLAP);

	// validation check
	if (image.rows != mask.rows || image.cols != mask.cols)
		exit(ERR_SIZE_MISMATCH);

	if (roi.x < 1 || roi.x + roi.width > image.cols 
		|| roi.y < 1 || roi.y + roi.height > image.rows)
	{
		printf("Invalid ROI: out of boarder!\n");
		system("pause");
		exit(ERR_SEE_NOTICE);
	}

	if (blockHeight < 5 || blockWidth < 5)
	{
		printf("Block too small, check paramters!\n");
		system("pause");
		exit(ERR_SEE_NOTICE);
	}

	// convert to various color space
	cvtColor(image, hsvImage, CV_BGR2HSV);
	cvtColor(image, labImage, CV_BGR2Lab);
	cvtColor(image, grayImage, CV_BGR2GRAY);

	feat.clear();

	// extract feature in each block
	for (int i = 0; i < IMAGE_PARTITION_Y; i++)
	{
		for (int j = 0; j < IMAGE_PARTITION_X; j++)
		{
			int		row = roi.y + i * stepHeight;								// block start y
			int		col = roi.x + j * stepWidth;								// block start x
			Rect	blockROI = Rect(col, row, blockWidth, blockHeight);			// block ROI

			extract_feature_block(blockFeatBuffer, blockROI);
			feat.insert(feat.end(), blockFeatBuffer.begin(), blockFeatBuffer.end());
		}
	}

}

void pri_feat::get_local_descriptor(int row, int col, vector<LDType> &ldBuffPixel)
{
	// boundary check, disabled
	/*if (row < 1 || row > hsvImage.rows || col < 1 || col > hsvImage.cols)
	{
		exit(ERR_OUT_OF_BOUNDARY);
	}*/

	int			h, s, v;
	int			hl, ht, hb, hr, sl, st, sb, sr, vl, vt, vb, vr;
	LDType		val;

	ldBuffPixel.clear();

	h = hsvImage.at<Vec3b>(row, col)[0];
	s = hsvImage.at<Vec3b>(row, col)[1];
	v = hsvImage.at<Vec3b>(row, col)[2];
	hl = hsvImage.at<Vec3b>(row, col - 1)[0];
	sl = hsvImage.at<Vec3b>(row, col - 1)[1];
	vl = hsvImage.at<Vec3b>(row, col - 1)[2];
	ht = hsvImage.at<Vec3b>(row - 1, col)[0];
	st = hsvImage.at<Vec3b>(row - 1, col)[1];
	vt = hsvImage.at<Vec3b>(row - 1, col)[2];
	hr = hsvImage.at<Vec3b>(row, col + 1)[0];
	sr = hsvImage.at<Vec3b>(row, col + 1)[1];
	vr = hsvImage.at<Vec3b>(row, col + 1)[2];
	hb = hsvImage.at<Vec3b>(row + 1, col)[0];
	sb = hsvImage.at<Vec3b>(row + 1, col)[1];
	vb = hsvImage.at<Vec3b>(row + 1, col)[2];

	// hue, 1st order, 2nd order
	val = h / 180.f;
	ldBuffPixel.push_back(val);
	val = (hr - h) / 180.f;
	ldBuffPixel.push_back(val);
	val = (hb - h) / 180.f;
	ldBuffPixel.push_back(val);
	val = (hl + hr - (h << 1)) / 180.f;
	ldBuffPixel.push_back(val);
	val = (hb + ht - (h << 1)) / 180.f;
	ldBuffPixel.push_back(val);

	// saturation, 1st order, 2nd order
	val = s / 255.f;
	ldBuffPixel.push_back(val);
	val = (sr - s) / 255.f;
	ldBuffPixel.push_back(val);
	val = (sb - s) / 255.f;
	ldBuffPixel.push_back(val);
	val = (sl + sr - (s << 1)) / 255.f;
	ldBuffPixel.push_back(val);
	val = (sb + st - (s << 1)) / 255.f;
	ldBuffPixel.push_back(val);

	// value, 1st order, 2nd order
	val = v / 255.f;
	ldBuffPixel.push_back(val);
	val = (vr - v) / 255.f;
	ldBuffPixel.push_back(val);
	val = (vb - v) / 255.f;
	ldBuffPixel.push_back(val);
	val = (vl + vr - (v << 1)) / 255.f;
	ldBuffPixel.push_back(val);
	val = (vb + vt - (v << 1)) / 255.f;
	ldBuffPixel.push_back(val);

}

void pri_feat::extract_feature_block(vector<FeatureType> &blockFeat, Rect blkROI)
{
	int			r, g, b;
	int			h, s, v;
	int			l, a, bb;
	//int			hl, ht, hb, hr, sl, st, sb, sr, vl, vt, vb, vr;
	//LDType		val;
	int			bin, offset;
	int			pixelCount = 0;

	// reset buffer for histogram
	ldBuffer.clear();
	blockFeat.clear();
	blockFeat.resize(128 + COLOR_LABEL_COUNT);
	std::fill(blockFeat.begin(), blockFeat.end(), 0.f);

	for (int row = blkROI.y; row < blkROI.y + blkROI.height; row++)
	{
		for (int col = blkROI.x; col < blkROI.x + blkROI.width; col++)
		{
			// check mask at this pixel and its 4 neighbors
			if (mask.at<UChar>(row, col)
				&& mask.at<UChar>(row, col - 1)
				&& mask.at<UChar>(row, col + 1)
				&& mask.at<UChar>(row - 1, col)
				&& mask.at<UChar>(row + 1, col))
			{
				// retrieve values
				b = image.at<Vec3b>(row, col)[0];
				g = image.at<Vec3b>(row, col)[1];
				r = image.at<Vec3b>(row, col)[2];
				h = hsvImage.at<Vec3b>(row, col)[0];
				s = hsvImage.at<Vec3b>(row, col)[1];
				v = hsvImage.at<Vec3b>(row, col)[2];
				l = labImage.at<Vec3b>(row, col)[0];
				a = labImage.at<Vec3b>(row, col)[1];
				bb = labImage.at<Vec3b>(row, col)[2];

				// --------------------------quantize to histogram----------------------------//
				// [0-136] --> RGB, HSV, Lab
				// [0-15] --> red
				bin = r >> 4;
				offset = 0;
				blockFeat.at(bin + offset)++;
				// [16-31] --> green
				bin = g >> 4;
				offset += 16;
				blockFeat.at(bin + offset)++;
				// [32 - 47] --> blue
				bin = b >> 4;
				offset += 16;
				blockFeat.at(bin + offset)++;
				// [48 - 63] --> hue
				bin = (h << 4) / 180;			// original h ~ [0 180]
				//assert(bin >= 0 && bin < 16);
				offset += 16;
				blockFeat.at(bin + offset)++;
				// [64 - 79] --> saturation
				bin = s >> 4;
				offset += 16;
				blockFeat.at(bin + offset)++;
				// [80 - 95] --> value
				bin = v >> 4;
				offset += 16;
				blockFeat.at(bin + offset)++;
				// [96 - 111] --> a
				bin = a >> 4;
				offset += 16;
				blockFeat.at(bin + offset)++;
				// [112 - 127] --> b
				bin = bb >> 4;
				offset += 16;
				blockFeat.at(bin + offset)++;
				// [128 - 143] or [128 - 139] --> color label
				r >>= 4;
				g >>= 4;
				b >>= 4;
				bin = m_hsvColorLabels.at((r << 8) + (g << 4) + b);
				offset += 16;
				blockFeat.at(bin + offset)++;
				//------------------------------end histogram------------------------------//

				// ------------------------compute local descriptor------------------------//
		
				vector<LDType>	ldBuffPixel;
				get_local_descriptor(row, col, ldBuffPixel);
				ldBuffer.insert(ldBuffer.end(), ldBuffPixel.begin(), ldBuffPixel.end());


				// --------------------------end local descriptor--------------------------//

				pixelCount++;
			}
		}
	}

	// L1-norm color histogram
	if (pixelCount > 0)
	{
		for (vector<FeatureType>::iterator iter = blockFeat.begin(); iter != blockFeat.end(); iter++)
			*iter /= pixelCount;
	}

	// add color percent feature
	int		percentPartition = 5;
	offset = 128;

	for (int i = 0; i < COLOR_LABEL_COUNT; i++)
	for (int j = 0; j < 5; j++)
	{
		bin = i * percentPartition + j + offset;
		if (blockFeat.at(bin) > j / (FeatureType)percentPartition)
			blockFeat.push_back(1);
		else
			blockFeat.push_back(0);
	}

	// LBP feature


	// fisher encoder
	vector<FeatureType> fisherBuffer;
	fisherBuffer.resize(LD_DIM * LD_NUM_CLUSTERS * 2);
	fill(fisherBuffer.begin(), fisherBuffer.end(), 0.f);
	float*		means = (float*)vl_gmm_get_means(m_gmm);
	vl_fisher_encode(fisherBuffer.data(), VL_TYPE_FLOAT, vl_gmm_get_means(m_gmm), LD_DIM, LD_NUM_CLUSTERS,
		vl_gmm_get_covariances(m_gmm), vl_gmm_get_priors(m_gmm), ldBuffer.data(), ldBuffer.size() / LD_DIM, VL_FISHER_FLAG_IMPROVED);

	blockFeat.insert(blockFeat.end(), fisherBuffer.begin(), fisherBuffer.end());
	
}


int	pri_feat::collect_local_descriptors()
{
	size_t	maxSize = MAX_VECTOR_SIZE;

	if (image.rows < 1 || image.cols < 1)
	{
		exit(ERR_SEE_NOTICE);
		printf("Invalid access to image, empty or not exist!\n");
	}

	cvtColor(image, hsvImage, CV_BGR2HSV);
	vector<LDType> ldBuffPixel;

	for (int row = 1; row < image.rows - 1; row++)
	for (int col = 1; col < image.cols - 1; col++)
	{
		get_local_descriptor(row, col, ldBuffPixel);
		if (ldBuffer.size() + ldBuffPixel.size() < maxSize)
			ldBuffer.insert(ldBuffer.end(), ldBuffPixel.begin(), ldBuffPixel.end());
		else
			return 1;
	}

	return 0;
}

void pri_feat::save_pairwise_feature_block()
{
	if (imgFeat.size() < 1)
	{
		return;
	}

	create_pair_index();

	const int	numPartitions = IMAGE_PARTITION_Y * IMAGE_PARTITION_X;
	size_t		featLen = imgFeat[0].size() / numPartitions;
	string		path = ROOT_PATH + string("cache/");

	if (USE_UNIFIED_MODEL)
	{
		ofstream	file(path + "blocks.trdat", ios::out | ios::trunc);
		if (!file.is_open())
			exit(ERR_FILE_UNABLE_TO_OPEN);

		// intra pairs first
		for (int i = 0; i < pairIdxIntra.size(); i++)
		{
			// switch regions to operate
			for (int k = 0; k < numPartitions; k++)
			{
				// feature range
				vector<FeatureType>::iterator	startIter = imgFeat[pairIdxIntra[i][0]].begin() + featLen * k;
				vector<FeatureType>::iterator	endIter = startIter + featLen;

				// partial feature
				vector<FeatureType> f1(startIter, endIter);
				startIter = imgFeat[pairIdxIntra[i][1]].begin() + featLen * k;
				endIter = startIter + featLen;
				vector<FeatureType> f2(startIter, endIter);


				// write pairwise feature to file
				file << 1 << "\t";									// intra pair label = 1
				write_similarity_to_file(file, f1, f2);
			}
		}

		// then inter pairs
		for (int i = 0; i < pairIdxInter.size(); i++)
		{
			// switch regions to operate
			for (int k = 0; k < numPartitions; k++)
			{
				// feature range
				vector<FeatureType>::iterator	startIter = imgFeat[pairIdxInter[i][0]].begin() + featLen * k;
				vector<FeatureType>::iterator	endIter = startIter + featLen;

				// partial feature
				vector<FeatureType> f1(startIter, endIter);
				startIter = imgFeat[pairIdxInter[i][1]].begin() + featLen * k;
				endIter = startIter + featLen;
				vector<FeatureType> f2(startIter, endIter);


				// write pairwise feature to file
				file << -1 << "\t";									// inter pair label = 0
				write_similarity_to_file(file, f1, f2);
			}
		}

		file.close();
	}
	else
	{
		vector<ofstream> fileList;
		fileList.resize(numPartitions);
		char filename[260];

		// open files to write
		for (int i = 0; i < numPartitions; i++)
		{
			sprintf(filename, "block_%d.trdat", i);
			fileList[i].open(path + filename, ios::out | ios::trunc);
			if (!fileList[i].is_open())
				exit(ERR_FILE_UNABLE_TO_OPEN);
		}


		// intra pairs first
		for (int i = 0; i < pairIdxIntra.size(); i++)
		{
			// switch regions to operate
			for (int k = 0; k < numPartitions; k++)
			{
				// feature range
				vector<FeatureType>::iterator	startIter = imgFeat[pairIdxIntra[i][0]].begin() + featLen * k;
				vector<FeatureType>::iterator	endIter = startIter + featLen;

				// partial feature
				vector<FeatureType> f1(startIter, endIter);
				startIter = imgFeat[pairIdxIntra[i][1]].begin() + featLen * k;
				endIter = startIter + featLen;
				vector<FeatureType> f2(startIter, endIter);


				// write pairwise feature to file
				fileList[k] << 1 << "\t";									// intra pair label = 1
				write_similarity_to_file(fileList[k], f1, f2);
			}
		}

		// then inter pairs
		for (int i = 0; i < pairIdxInter.size(); i++)
		{
			// switch regions to operate
			for (int k = 0; k < numPartitions; k++)
			{
				// feature range
				vector<FeatureType>::iterator	startIter = imgFeat[pairIdxInter[i][0]].begin() + featLen * k;
				vector<FeatureType>::iterator	endIter = startIter + featLen;

				// partial feature
				vector<FeatureType> f1(startIter, endIter);
				startIter = imgFeat[pairIdxInter[i][1]].begin() + featLen * k;
				endIter = startIter + featLen;
				vector<FeatureType> f2(startIter, endIter);


				// write pairwise feature to file
				fileList[k] << -1 << "\t";									// inter pair label = 0
				write_similarity_to_file(fileList[k], f1, f2);
			}
		}

		// close files
		for (int i = 0; i < numPartitions; i++)
		{
			fileList[i].close();
		}
	}
}

float pri_feat::similarity_score(float f1, float f2)
{
	float	val, v1, v2, minVal, maxVal;

	v1 = abs(f1);
	v2 = abs(f2);

	if (v1 < v2)
	{
		minVal = v1;
		maxVal = v2;
	}
	else
	{
		minVal = v2;
		maxVal = v1;
	}

	if (maxVal == 0 || minVal == 0)
		val = 0;
	else
	{
		val = minVal / maxVal;
	}

	return val;
}


void pri_feat::write_similarity_to_file(ofstream &file, vector<FeatureType> f1, vector<FeatureType> f2)
{
	if (f1.size() != f2.size())
		return;

	float		val;
	for (int i = 0; i < f1.size(); i++)
	{
		val = similarity_score(f1[i], f2[i]);

		if ( val != 0)
		{
			file << i+1 << ":" << val << " ";
		}
	}

	// if the last element is 0, write it to maintain the proper size
	if (val == 0)
		file << f1.size() << ":" << val << " ";

	file << endl;		// write end of line
}


void pri_feat::train_block_models()
{
	// call external exe to train SVM models
	char	filename[260];
	char	modelname[260];
	string	exePath = ROOT_PATH + string("resource/");
	string	exeFile = string("ranksvm64.exe");

	if (USE_UNIFIED_MODEL)
	{
		sprintf(filename, "../cache/blocks.trdat");
		sprintf(modelname, "../cache/blocks.model");
		char	systemcall[2048];
		sprintf(systemcall, "cd /d %s && %s %s %s", exePath.c_str(), exeFile.c_str(), filename, modelname);
		system(systemcall);
	}
	else
	{
		const int	numPartitions = IMAGE_PARTITION_Y * IMAGE_PARTITION_X;
		for (int i = 0; i < numPartitions; i++)
		{
			sprintf(filename, "../cache/block_%d.trdat", i);
			sprintf(modelname, "../cache/block_%d.model", i);
			char	systemcall[2048];
			sprintf(systemcall, "cd /d %s && %s %s %s", exePath.c_str(), exeFile.c_str(), filename, modelname);
			system(systemcall);
		}
	}
}

void pri_feat::train_image_model()
{
	// call external exe to train SVM model
	char	filename[260];
	char	modelname[260];
	string	exePath = ROOT_PATH + string("resource/");
	string	exeFile = string("ranksvm64.exe");

	sprintf(filename, "../cache/image.trdat");
	sprintf(modelname, "../cache/image.model");
	char	systemcall[2048];
	sprintf(systemcall, "cd /d %s && %s %s %s", exePath.c_str(), exeFile.c_str(), filename, modelname);
	system(systemcall);
}


void pri_feat::load_block_weights()
{
	blockWeights.clear();

	if (USE_UNIFIED_MODEL)
	{
		char		filename[260];
		sprintf(filename, "cache/blocks.model");
		ifstream	file(ROOT_PATH + string(filename), ios::in);
		if (!file.is_open())
			exit(ERR_FILE_UNABLE_TO_OPEN);

		string	line = "";
		int		featLen = 0;

		while (line != "w")
		{
			getline(file, line);

			// try to get the feature length
			if (line.substr(0, 10) == "nr_feature")
				featLen = stoi(line.substr(11, line.size() - 11));
		}

		if (featLen < 1)
		{
			printf("Error feature length detected!\n");
			exit(ERR_SEE_NOTICE);
		}

		// temporal weights vector
		vector<float>	weights;
		for (int j = 0; j < featLen; j++)
		{
			float	buf;
			file >> buf;
			weights.push_back(buf);
		}
		// push back to 2-D vector
		blockWeights.push_back(weights);

		file.close();
	}
	else
	{
		const int	numPartitions = IMAGE_PARTITION_Y * IMAGE_PARTITION_X;
		for (int i = 0; i < numPartitions; i++)
		{
			char		filename[260];
			sprintf(filename, "cache/block_%d.model", i);
			ifstream	file(ROOT_PATH + string(filename), ios::in);
			if (!file.is_open())
				exit(ERR_FILE_UNABLE_TO_OPEN);

			string	line = "";
			int		featLen = 0;

			while (line != "w")
			{
				getline(file, line);

				// try to get the feature length
				if (line.substr(0, 10) == "nr_feature")
					featLen = stoi(line.substr(11, line.size() - 11));
			}

			if (featLen < 1)
			{
				printf("Error feature length detected!\n");
				exit(ERR_SEE_NOTICE);
			}

			// temporal weights vector
			vector<float>	weights;
			for (int j = 0; j < featLen; j++)
			{
				float	buf;
				file >> buf;
				weights.push_back(buf);
			}

			file.close();

			// push back to 2-D vector
			blockWeights.push_back(weights);
		}
	}
}

void pri_feat::load_image_weights()
{
	imageWeights.clear();

	char		filename[260];
	sprintf(filename, "cache/image.model");
	ifstream	file(ROOT_PATH + string(filename), ios::in);
	if (!file.is_open())
		exit(ERR_FILE_UNABLE_TO_OPEN);

	string	line = "";
	int		featLen = 0;

	while (line != "w")
	{
		getline(file, line);

		// try to get the feature length
		if (line.substr(0, 10) == "nr_feature")
			featLen = stoi(line.substr(11, line.size() - 11));
	}

	if (featLen < 1)
	{
		printf("Error feature length detected!\n");
		exit(ERR_SEE_NOTICE);
	}

	// temporal weights vector
	for (int j = 0; j < featLen; j++)
	{
		float	buf;
		file >> buf;
		imageWeights.push_back(buf);
	}

	file.close();

}

void pri_feat::save_pairwise_feature_image()
{
	// validation check
	if (imgFeat.size() < 1)
	{
		return;
	}

	load_block_weights();
	create_pair_index();

	if (blockWeights.size() < 1)
	{
		return;
	}

	const int	numPartitions = IMAGE_PARTITION_Y * IMAGE_PARTITION_X;
	size_t		featLen = imgFeat[0].size() / numPartitions;

	if (blockWeights[0].size() != featLen)
	{
		printf("Load weights length mismatch!\n");
		exit(ERR_SEE_NOTICE);
	}

	if (pairIdxIntra.size() < 1 || pairIdxInter.size() < 1)
		create_pair_index();

	// calculate similarities and write to file
	ofstream	file(ROOT_PATH + string("cache/image.trdat"), ios::out | ios::trunc);
	if (!file.is_open())
		exit(ERR_FILE_UNABLE_TO_OPEN);

	// intra pairs first
	for (int i = 0; i < pairIdxIntra.size(); i++)
	{
		vector<float>	combFeat;
		get_combine_image_feature(combFeat, imgFeat[pairIdxIntra[i][0]], imgFeat[pairIdxIntra[i][1]]);
		file << 1 << "\t";		// intra label 1
		for (int j = 0; j < combFeat.size(); j++)
			file << j + 1 << ":" << combFeat[j] << " ";
		file << endl;
	}

	// then inter pairs
	for (int i = 0; i < pairIdxInter.size(); i++)
	{
		vector<float>	combFeat;
		get_combine_image_feature(combFeat, imgFeat[pairIdxInter[i][0]], imgFeat[pairIdxInter[i][1]]);
		file << -1 << "\t";		// inter label -1
		for (int j = 0; j < combFeat.size(); j++)
			file << j + 1 << ":" << combFeat[j] << " ";
		file << endl;
	}

	file.close();
}

void pri_feat::get_combine_image_feature(vector<FeatureType> & combFeat, vector<FeatureType> f1, vector<FeatureType> f2)
{
	const int numBlocks = IMAGE_PARTITION_Y * IMAGE_PARTITION_X;

	// size check
	if (f1.size() != f2.size())
		return;

	// reset
	combFeat.clear();

	
	if (blockWeights.size() == 1 && USE_UNIFIED_MODEL)
	{
		// duplicate weights
		vector<float>	weights = blockWeights[0];
		for (int i = 1; i < numBlocks; i++)
			blockWeights.push_back(weights);
	}
	else if (blockWeights.size() > 1)
	{
		if (blockWeights.size() != numBlocks)
			exit(ERR_SIZE_MISMATCH);
	}
	else
	{
		printf("Error: block weights size invalid!\n");
		exit(ERR_SEE_NOTICE);
	}

	size_t	featLen = blockWeights[0].size();
	if (featLen * numBlocks != f1.size())
		exit(ERR_SIZE_MISMATCH);

	// basic block similarity
	int		featPtr = 0;
	for (int i = 0; i < numBlocks; i++)
	{
		float	sim = 0;
		for (int j = 0; j < featLen; j++)
		{
			sim += blockWeights[i][j] * similarity_score((float)f1[featPtr], (float)f2[featPtr]);
			featPtr++;
		}
		combFeat.push_back(sim);
	}

	// calculate mean and standard deviation
	float	sum = std::accumulate(combFeat.begin(), combFeat.end(), 0.f);
	float	mean = sum / numBlocks;
	float	accum = 0.0;
	std::for_each(combFeat.begin(), combFeat.end(), 
	[&](const float d)
	{ 
		accum += (d - mean)*(d - mean); 
	}
	);

	float stdev = sqrt(accum / numBlocks);
	combFeat.push_back(mean);
	combFeat.push_back(stdev);
	

	// try some combination
	// combination of two blocks, mean and standard dev
	for (int i = 0; i < numBlocks; i++)
	for (int j = i + IMAGE_PARTITION_X; j < numBlocks; j += IMAGE_PARTITION_X)
	{
		float	sim = (combFeat[i] + combFeat[j]) / 2;
		float	stdev2 = (combFeat[i] - sim) * (combFeat[i] - sim) + (combFeat[j] - sim) * (combFeat[j] - sim);
		stdev2 = sqrt(stdev2 / 2);
		combFeat.push_back(sim);
		combFeat.push_back(stdev2);
	}
	// combination of three blocks
	for (int i = 0; i < numBlocks; i++)
	for (int j = i + IMAGE_PARTITION_X; j < numBlocks; j += IMAGE_PARTITION_X)
	for (int k = j + IMAGE_PARTITION_X; k < numBlocks; k += IMAGE_PARTITION_X)
	{
		float	sim = (combFeat[i] + combFeat[j] + combFeat[k]) / 3;
		float	stdev3 = (combFeat[i] - sim) * (combFeat[i] - sim) 
						+ (combFeat[j] - sim) * (combFeat[j] - sim) 
						+ (combFeat[k] - sim) * (combFeat[k] - sim);
		stdev3 = sqrt(stdev3 / 3);
		combFeat.push_back(sim);
		combFeat.push_back(stdev3);
	}
	
}


float pri_feat::image_pairwise_score(int idx1, int idx2)
{

	// init
	if (blockWeights.size() < 1)
		load_block_weights();

	if (imageWeights.size() < 1)
		load_image_weights();

	if (blockWeights.size() < 1 || imageWeights.size() < 1)
	{
		printf("Error: SVM model weights invalid!\n");
		exit(ERR_SEE_NOTICE);
	}

	
	vector<float>	combFeat;

	// get lower level score
	get_combine_image_feature(combFeat, imgFeat[idx1], imgFeat[idx2]);

	float	score = 0;
	for (int i = 0; i < imageWeights.size(); i++)
		score += combFeat[i] * imageWeights[i];

	return score;

}

void pri_feat::rank_cmc()
{
	results.clear();
	ranks.clear();
	ranks.resize(numPersons);
	fill(ranks.begin(), ranks.end(), 0.f);

	if (imgFeat.size() < 1)
		return;

	for (int i = 0; i < numPersons; i++)
	{
		vector<sort_descend> search;
		int	qIdx = queryIdx[i][0];		// use the first image in each query
		for (int j = 0; j < numPersons; j++)
		for (int k = 1; k < numShots; k++)
		{
			int		gIdx = queryIdx[j][k];	// gallery index
			float	score = image_pairwise_score(qIdx, gIdx);
			sort_descend sort;
			sort.score = score;
			sort.id = imgFeatID[gIdx];
			search.push_back(sort);
		}

		// sort according to the scores
		sort(search.begin(), search.end());
		results.push_back(search);
	}

	// statistic rank info
	for (int i = 0; i < numPersons; i++)
	{
		int	qId = queryIdx[i][numShots];	// query image ID
		int	rank = numPersons-1;
		for (int j = 0; j < results[i].size(); j++)
		{
			if (qId == results[i][j].id)
			{
				rank = j;
				break;
			}
		}
		ranks[rank]++;
	}

	// get percentage
	for (int i = 0; i < ranks.size(); i++)
		ranks[i] /= numPersons;

	// get cumulative sum
	vector<float> tmp = ranks;
	std::partial_sum(tmp.begin(), tmp.end(), ranks.begin(), plus<float>());
}

float pri_feat::get_rank_n(int n)
{
	// n = n - 1 to fit the structure
	return ranks[n - 1];
}

#if DEV_DEBUG

void pri_feat::debug_show_top_n(int n)
{
	// try to put the top n images side by side in horizontal direction
	Mat		qImage = imread(filenames[0], 1);
	// create a canvas that is height x ((n+2) * width) : one query, one target, n top
	Mat		canvas;
	canvas.create(qImage.rows, qImage.cols * (n + 2), qImage.type());

	Rect	roi;
	char	filename[260];

	for (int idx = 0; idx < results.size(); idx++)
	{
		cout << "--------------------------------------------------" << endl;
		cout << "Debug top n for: " << queryIdx[idx][numShots] << endl;
		// first copy query image
		roi = Rect(0, 0, qImage.cols, qImage.rows);
		sprintf(filename, "%sdata/VIPeR/%03d_00.bmp", ROOT_PATH, queryIdx[idx][numShots]);
		qImage = imread(filename, 1);
		qImage.copyTo(canvas(roi));

		// top n images
		for (int i = 0; i < n; i++)
		{
			roi = Rect(qImage.cols * (i + 1), 0, qImage.cols, qImage.rows);
			sprintf(filename, "%sdata/VIPeR/%03d_01.bmp", ROOT_PATH, results[idx][i].id);
			cout << "Top score #" << i + 1 << ": " << results[idx][i].score << endl;
			qImage = imread(filename, 1);
			qImage.copyTo(canvas(roi));
		}

		// target image
		cout << "----------" << endl;
		int	rank = numPersons;
		for (int i = 0; i < numPersons; i++)
		{
			if (results[idx][i].id == queryIdx[idx][numShots])
			{
				rank = i;
				break;
			}
		}
		cout << "self rank: " << rank + 1 << endl;
		cout << "self target score: " << results[idx][rank].score << endl;
		roi = Rect(canvas.cols - qImage.cols, 0, qImage.cols, qImage.rows);
		sprintf(filename, "%sdata/VIPeR/%03d_01.bmp", ROOT_PATH, queryIdx[idx][numShots]);
		qImage = imread(filename, 1);
		qImage.copyTo(canvas(roi));

		// show image
		destroyAllWindows();
		imshow("Top n debug window", canvas);
		waitKey(0);
	}
	
}

#endif