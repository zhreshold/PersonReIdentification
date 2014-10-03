/***********************************************************************/
/*
/*   Script File: hkmeans.h
/*
/*   Description:
/*
/*   Hierachical kmeans class declaration
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

#ifndef _HKMEANS_H_
#define _HKMEANS_H_

#include "define.h"
#include <fstream>
#include <vector>
extern "C"
{
#include <vl\kmeans.h>
}



using namespace std;

// ! whether or not store patch data for illustration
#define	STORE_PATCH_DATA	1
#define MIN_NODE_SIZE		100

typedef float	hk_type;
typedef unsigned char	uchar;

struct hkPatch
{
	int		flag;

	int		width;
	int		height;
	int		channel;
	int		step;

	// ! store patch data if turned on
	vector<uchar> data;

	// ! patch feature
	vector<hk_type> feature;
};

struct hkNode
{
	hkNode*	parent;
	int		level;
	int		dim;
	vector<hk_type> center;
	vector<hkNode> childs;
};

void destroyHkTree(hkNode &root);



class hkmeans
{
public:
	hkmeans();
	~hkmeans();

	// ! partition K in each level
	int	K;
	// ! depth of tree
	int depth;
	// ! root node
	hkNode root;
	// ! training data
	vector<hkPatch>	patches;


	// public functions
	// ! clear all data
	void clear();

	// ! train hierachical kmeans
	void train();
	// ! refine tree
	void refine(hkNode* node, int totalData);
	// ! save tree
	void write_tree_to_file(const char* filename);
	// ! load tree
	void read_tree_from_file(const char* filename);
	// ! debug centers
	void debug_output_tree();

	// static functions
	void run_with_vlfeat_kmeans(float* data, float* centers, vector<vl_uint32> &assignments, int numData, int dim, int numCenters, int maxIter);
	void train_recursive_tree(hkNode* node, vector<float> &data, int dim, int numCenters, int maxIter);
	void quantize(hkNode* node, vector<float> &data, vector<int> &assignment, int dim);
	void encode(hkNode* node, vector<float> &data, vector<int> &hist, int dim, int flag);

private:
	VlKMeans*		kmeans;
	
	// private functions
	
	void write_node_to_file(ofstream &fp, hkNode* node, vector<int> code);
	void read_node_from_file(ifstream &fp, hkNode* node);
	
};





#endif