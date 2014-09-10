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
#include <algorithm>
#include <fstream>


pri_feat::pri_feat()
{
	// init
	numPersons = 0;
	numShots = 0;

	// compute hsv color label look-up table
	get_hsv_color_labels(m_hsvColorLabels);

	// load GMM
	load_gmm(m_gmm);
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

	for (int i = 0; i < m_gmm.size(); i++)
	{
		vl_gmm_delete(m_gmm[i]);
	}
	m_gmm.clear();
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

void pri_feat::load_gmm(vector<VlGMM*> &gmm)
{
	ifstream	file(ROOT_PATH + string("cache/GMM.dat"), ios::in);
	if (!file.is_open())
		exit(ERR_FILE_NOT_EXIST);

	vector<LDType>	gBuffer;
	VlGMM*			gmmPtr;
	int				numPartition = IMAGE_PARTITION_Y * IMAGE_PARTITION_X;
	int				gSize = LD_DIM * LD_NUM_CLUSTERS;
	gBuffer.resize(gSize);
	//fill(gBuffer.begin(), gBuffer.end(), 0);

	// load with means, covariances and priors
	for (int i = 0; i < numPartition; i++)
	{
		gmmPtr = vl_gmm_new(VL_TYPE_FLOAT, LD_DIM, LD_NUM_CLUSTERS);
		gmm.push_back(gmmPtr);
		for (int j = 0; j < gSize; j++)
		{
			file >> gBuffer[j];
		}
		vl_gmm_set_means(gmm[i], gBuffer.data());
		for (int j = 0; j < gSize; j++)
		{
			file >> gBuffer[j];
		}
		vl_gmm_set_covariances(gmm[i], gBuffer.data());
		for (int j = 0; j < gSize; j++)
		{
			file >> gBuffer[j];
		}
		vl_gmm_set_priors(gmm[i], gBuffer.data());
	}

	file.close();
}


void pri_feat::get_hsv_color_labels(vector<UChar> &hsvColorLabels)
{
	int			vr, vg, vb;
	float		r, g, b, h, s, v;
	int			CL;		

	// init look-up  table
	hsvColorLabels.clear();
	hsvColorLabels.resize(4096);
	fill(hsvColorLabels.begin(), hsvColorLabels.end(), 0);

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
		hsvColorLabels[vr * 256 + vg * 16 + vb] = CL;

	}
}

void pri_feat::rgb2hsv(float r, float g, float b, float &h, float &s, float &v)
{
	// fast version converter

	float K = 0.f;

	if (g < b)
	{
		std::swap(g, b);
		K = -1.f;
	}

	if (r < g)
	{
		std::swap(r, g);
		K = -2.f / 6.f - K;
	}

	float chroma = r - std::min(g, b);
	h = fabs(K + (g - b) / (6.f * chroma + 1e-20f));
	s = chroma / (r + 1e-20f);
	v = r;
}


void pri_feat::extract_feature()
{
	Mat					image;
	Mat					mask;									// foreground mask
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
			extract_feature_image(tmpFeat, image, mask);
			// push back feature and id
			imgFeat.push_back(tmpFeat);
			imgFeatID.push_back(query_person[numShots]);

		}
	}
	
}


void pri_feat::extract_feature_image(vector<FeatureType> &feat, Mat image, Mat mask)
{
	
	// copy from global ROI
	Rect				roi = gROI;
	int					blockHeight, blockWidth;
	vector<FeatureType>	blockFeatBuffer, imageFeatBuffer;
	Mat					hsvImage, labImage;

	// calculate block width and height
	blockHeight = roi.height / IMAGE_PARTITION_Y;
	blockWidth = roi.width / IMAGE_PARTITION_X;

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

	// convert to HSV color space
	cvtColor(image, hsvImage, CV_BGR2HSV);
	cvtColor(image, labImage, CV_BGR2Lab);

	// extract feature in each block
	for (int i = 0; i < IMAGE_PARTITION_Y; i++)
	{
		for (int j = 0; j < IMAGE_PARTITION_X; j++)
		{
			int		row = roi.y + i * blockHeight;								// block start y
			int		col = roi.x + j * blockWidth;								// block start x
			Rect	blockROI = Rect(col, row, blockWidth, blockHeight);			// block ROI

			extract_feature_block(blockFeatBuffer, image, hsvImage, labImage, mask, blockROI);
		}
	}

}

void pri_feat::extract_feature_block(vector<FeatureType> &blockFeat, Mat image, Mat hsvImage, Mat labImage, Mat mask, Rect blkROI)
{
	int			r, g, b;
	int			h, s, v;
	int			l, a, bb;
	int			hl, ht, hb, hr, sl, st, sb, sr, vl, vt, vb, vr;
	LDType		val;
	int			bin, offset;
	int			pixelCount = 0;

	// reset buffer for histogram
	ldBuffer.clear();
	blockFeat.clear();
	if (COLOR_LABEL_METHOD == 0)
		blockFeat.resize(140);
	else if (COLOR_LABEL_METHOD == 1)
		blockFeat.resize(144);
	fill(blockFeat.begin(), blockFeat.end(), 0);

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
				bin = m_hsvColorLabels.at(r << 8 + g << 4 + b);
				offset += 16;
				blockFeat.at(bin + offset)++;
				//------------------------------end histogram------------------------------//

				// ------------------------compute local descriptor------------------------//
		
				//retrieve values
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
				val = h / 180.0;
				ldBuffer.push_back(val);
				val = (hr - h) / 180.0;
				ldBuffer.push_back(val);
				val = (hb - h) / 180.0;
				ldBuffer.push_back(val);
				val = (hl + hr - h << 1) / 180.0;
				ldBuffer.push_back(val);
				val = (hb + ht - h << 1) / 180.0;
				ldBuffer.push_back(val);

				// saturation, 1st order, 2nd order



				// --------------------------end local descriptor--------------------------//

				pixelCount++;
			}
		}
	}
	
}