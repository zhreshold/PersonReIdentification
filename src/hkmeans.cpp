/***********************************************************************/
/*
/*   Script File: hkmeans.cpp
/*
/*   Description:
/*
/*   Hierachical kmeans class
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

#include "hkmeans.h"
#include <chrono>
#include <fstream>
#include <string>
#include <Windows.h>
#include <opencv\cv.h>
#include <opencv\highgui.h>

void destroyHkTree(hkNode &root)
{
	while (root.childs.size() > 0)
	{
		destroyHkTree(root.childs.back());
		root.childs.pop_back();
	}

	root.center.clear();

	root.parent = NULL;
	
}

hkmeans::hkmeans()
{
	K = 0;
	depth = 0;
	root.parent = NULL;
	root.level = 0;

	kmeans = vl_kmeans_new(VL_TYPE_FLOAT, VlDistanceL2);
	vl_kmeans_set_algorithm(kmeans, VlKMeansLloyd);
}

hkmeans::~hkmeans()
{
	clear();
	vl_kmeans_delete(kmeans);
}

void hkmeans::clear()
{
	patches.clear();
	destroyHkTree(root);
}

void hkmeans::train()
{
	if (K < 2)
	{
		printf("Error: K is too small!\n");
		exit(-1);
	}

	if (patches.size() < 1)
	{
		printf("Error: no training data!\n");
		exit(-1);
	}



	vector<float>	feat;
	int				dim = patches[0].feature.size();
	int				numData;
	vector<float>	centers;
	int				maxIter = 50;
	vector<vl_uint32>		assignments;
	centers.resize(K);


	// load data
	for (int i = 0; i < patches.size(); i++)
	{
		feat.insert(feat.end(), patches[i].feature.begin(), patches[i].feature.end());
	}
	numData = patches.size();

	train_recursive_tree(&root, feat, dim, K, maxIter);

	//refine(&root, numData);
	
	
}

void hkmeans::refine(hkNode* node, int totalData)
{
	if (node->childs.size() == 0)
	{
		// leaf node
		
	}

	for (int i = 0; i < node->childs.size(); i++)
	{
		refine(&node->childs[i], totalData);
	}
}


void hkmeans::run_with_vlfeat_kmeans(float* data, float* centers, vector<vl_uint32> &assignments, int numData, int dim, int numCenters, int maxIter)
{
	
	vl_kmeans_init_centers_with_rand_data(kmeans, data, dim, numData, numCenters);
	printf("Start...\n");
	//vl_kmeans_set_max_num_iterations(kmeans, 1);
	//for (int i = 0; i < maxIter; i++)
	//{
	//	printf("Iter # %d...", i);
	//	auto begin = std::chrono::high_resolution_clock::now();
	//	vl_kmeans_refine_centers(kmeans, data, numData);
	//	auto end = std::chrono::high_resolution_clock::now();
	//	auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
	//	printf("Elapsed time: %.2f seconds.\n", ms / 1000.0);

	//}

	vl_kmeans_set_max_num_iterations(kmeans, maxIter);
	auto begin = std::chrono::high_resolution_clock::now();
	vl_kmeans_refine_centers(kmeans, data, numData);
	auto end = std::chrono::high_resolution_clock::now();
	auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
	printf("Elapsed time: %.2f seconds.\n", ms / 1000.0);

	float*	c = (float*)vl_kmeans_get_centers(kmeans);

	for (int i = 0; i < numCenters * dim; i++)
		centers[i] = c[i];

	// assignment
	assignments.resize(numData);
	vector<float>	distances;
	distances.resize(numData);
	vl_kmeans_quantize(kmeans, assignments.data(), distances.data(), data, numData);

	
}

void hkmeans::train_recursive_tree(hkNode* node, vector<float> &data, int dim, int numCenters, int maxIter)
{
	int				numData;
	vector<float>	centers;
	vector<vl_uint32>		assignments;
	centers.resize(numCenters*dim);

	numData = data.size() / dim;
	printf("Level %d, #data: %d...\n", node->level, numData);
	run_with_vlfeat_kmeans(data.data(), centers.data(), assignments, numData, dim, numCenters, maxIter);
	node->center = centers;
	node->dim = dim;

	if (node->level < depth-1)
	{
		for (int i = 0; i < numCenters; i++)
		{
			vector<float> nextData;
			for (int j = 0; j < numData; j++)
			{
				if (assignments[j] == i)
				{
					nextData.insert(nextData.end(), data.begin() + j *dim, data.begin() + j * dim + dim);
				}
			}

			if (nextData.size()/dim > MIN_NODE_SIZE)
			{
				hkNode cNode;
				cNode.parent = node;
				cNode.level = node->level + 1;
				train_recursive_tree(&cNode, nextData, dim, numCenters, maxIter);
				node->childs.push_back(cNode);
			}
		}
	}
	else
	{
		// reached bottom of the tree
		
	}


}

void hkmeans::quantize(hkNode* node, vector<float> &data, vector<int> &assignment, int dim)
{
	if (dim != node->dim)
		return;

	float	minDist = 1e06;
	int		a = -1;
	for (int i = 0; i < K; i++)
	{
		float tmpDist = 0;
		for (int j = 0; j < dim; j++)
		{
			float dist = node->center[i*dim + j] - data[j];
			tmpDist += dist * dist;
		}
		if (tmpDist < minDist)
		{
			a = i;
			minDist = tmpDist;
		}
	}

	assignment.push_back(a);

	if (node->childs.size() > 0)
	{
		quantize(&node->childs[a], data, assignment, dim);
	}
	
	return;

}

void hkmeans::encode(hkNode* node, vector<float> &data, vector<int> &hist, int dim, int flag)
{
	if (dim != node->dim)
		return;

	if (flag == 1)
	{
		float	minDist = 1e06;
		int		a = -1;
		for (int i = 0; i < K; i++)
		{
			float tmpDist = 0;
			for (int j = 0; j < dim; j++)
			{
				float dist = node->center[i*dim + j] - data[j];
				tmpDist += dist * dist;
			}
			if (tmpDist < minDist)
			{
				a = i;
				minDist = tmpDist;
			}
		}

		for (int i = 0; i < K; i++)
		{
			if (i == a)
				hist.push_back(1);
			else
				hist.push_back(0);
		}

		for (int i = 0; i < node->childs.size(); i++)
		{
			if (i == a)
				encode(&node->childs[i], data, hist, dim, 1);
			else
				encode(&node->childs[i], data, hist, dim, 0);
		}
	}
	else
	{
		for (int i = 0; i < K; i++)
			hist.push_back(0);

		for (int i = 0; i < node->childs.size(); i++)
			encode(&node->childs[i], data, hist, dim, 0);
	}
}

void hkmeans::debug_output_tree()
{
#if STORE_PATCH_DATA
	cout << "Saving patch images..." << endl;
	for (int i = 0; i < patches.size(); i++)
	{
		cout << 100.0 * i / patches.size() << "%..." << endl;
		cv::Mat thisPatch(patches[i].height, patches[i].width, CV_8UC3, patches[i].data.data());
		vector<int> assign;
		quantize(&root, patches[i].feature, assign, patches[i].feature.size());
		string	folder("../../../debug/hkmeans/");
		for (int j = 0; j < assign.size(); j++)
		{
			folder = folder + std::to_string(assign[j]) + "/";
			if (CreateDirectory(folder.c_str(), NULL) ||
				ERROR_ALREADY_EXISTS == GetLastError())
			{
				// CopyFile(...)
				if ( j > 3)
					cv::imwrite(folder + std::to_string(i) + ".jpg", thisPatch);
			}
			else
			{
				// Failed to create directory.
				exit(-2);
			}

		}

	}

#endif
}

void hkmeans::write_tree_to_file(const char* filename)
{
	ofstream	fp;
	fp.open(filename, ios::out | ios::trunc);
	if (!fp.is_open())
	{
		cout << "Error: can't open file to write tree!" << endl;
		exit(-1);
	}

	// first K and depth
	fp << K << " " << depth << endl;
	vector<int> code;
	code.push_back(0);
	write_node_to_file(fp, &root, code);



	fp.close();

}

void hkmeans::write_node_to_file(ofstream &fp, hkNode* node, vector<int> code)
{
	fp << (int)(code.size()) << " ";
	for (int i = 0; i < code.size(); i++)
	{
		fp << code[i] << " ";
	}

	fp << node->level << " ";
	fp << node->dim << " ";
	for (int i = 0; i < node->center.size(); i++)
	{
		fp << node->center[i] << " ";
	}

	fp << endl;

	for (int i = 0; i < node->childs.size(); i++)
	{
		vector<int> childCode = code;
		childCode.push_back(i);
		write_node_to_file(fp, &node->childs[i], childCode);
	}
}


void hkmeans::read_tree_from_file(const char* filename)
{
	ifstream fp(filename, ios::in);
	if (!fp.is_open())
	{
		cout << "Error: can't open file to read tree!" << endl;
		exit(-1);
	}

	int		tmpK, tmpDepth;

	fp >> tmpK;
	fp >> tmpDepth;

	K = tmpK;
	depth = tmpDepth;

	int line = 1;
	while (!fp.eof())
	{
		//cout << line++ << endl;
		read_node_from_file(fp, &root);
	}


	fp.close();
}

void hkmeans::read_node_from_file(ifstream &fp, hkNode* node)
{
	int nCode = 0;
	int tmp;

	fp >> nCode;

	if (nCode < 1 || nCode > depth)
		return;

	vector<int> code;
	for (int i = 0; i < nCode; i++)
	{
		fp >> tmp;
		code.push_back(tmp);
	}

	hkNode tmpNode;

	fp >> tmp;
	tmpNode.level = tmp;
	fp >> tmp;
	tmpNode.dim = tmp;

	float c;
	for (int i = 0; i < tmpNode.dim * K; i++)
	{
		fp >> c;
		tmpNode.center.push_back(c);
	}

	// insert node
	hkNode* insertNode = node;

	if (code.size() == 1)
	{
		insertNode->center = tmpNode.center;
		insertNode->level = tmpNode.level;;
		insertNode->dim = tmpNode.dim;
	}
	else
	{
		for (int i = 1; i < code.size() - 1; i++)
		{
			insertNode = &insertNode->childs[code[i]];
		}

		insertNode->childs.push_back(tmpNode);
	}

	

}