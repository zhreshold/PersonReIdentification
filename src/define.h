/***********************************************************************/
/*                                                                    
/*   Script File: define.h                                                          
/*                        
/*   Description:
/*                                               
/*   Definitions and typedef      
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



#ifndef __DEFINE_H__
#define __DEFINE_H__

#define	_CRT_SECURE_NO_WARNINGS

/*-------------Do not change between dashed lines----------------*/

// macro
#define	Malloc(type, size)		(type*)malloc((size)*sizeof(type))

// typedef
typedef unsigned int			UInt;
typedef	unsigned char			UChar;
typedef unsigned long			ULong;
typedef	double					Double;
typedef bool					Bool;

enum ERROR_CODES
{
	ERR_OK,
	ERR_SEE_NOTICE,
	ERR_FILE_NOT_EXIST,
	ERR_FILE_UNABLE_TO_OPEN,
	ERR_OUT_OF_BOUNDARY,
	ERR_SIZE_MISMATCH,
	ERR_EMPTY_SET,
};

enum MODULE_PHASE
{
	PHASE_NO_STATUS,
	PHASE_TRAINING,
	PHASE_TESTING,
	PHASE_COLLECTING,
	
};


/*-------------Do not change between dashed lines----------------*/

typedef	float	FeatureType;						// final feature type
typedef	float	LDType;								// local descriptor type

#define	MAX_VECTOR_SIZE			16589934592
#define	ROOT_PATH				"../../../"
#define	RANDOM_SEED				7000				// set to fixed number for debug only


// datasets to test
enum DATASET_TYPES
{
	VIPER = 1, CUHK01,
};
#define	DATASET_TYPE			VIPER				// VIPER, CUHK01...
#define	NUM_PERSON_TRAIN		316					// number of individuals used for train
#define	IMAGE_PARTITION_X		1					// number of partitions along x
#define	IMAGE_PARTITION_Y		6					// number of partitions along y


// local descriptor & GMM model
#define	LD_DIM					15					// local descriptor dimension
#define	LD_NUM_CLUSTERS			30					// number of GMM clusters
#define	NEW_GMM_MODELS			0					// 1: train new models, 0: use existed models

// SVM
#define	USE_UNIFIED_MODEL		0					// 1: use single model for all regions, 0: use seperate models
#define	SVM_NEG_POS_RATIO		5					// maximum ratio of negtive/positive samples

/****************************************************************/

#endif // _DEFINE_H_