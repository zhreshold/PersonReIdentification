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
	ERR_FILE_NOT_EXIST,
	ERR_FILE_UNABLE_TO_OPEN,
	ERR_OUT_OF_BOUNDARY,
};

enum MODULE_PHASE
{
	PHASE_NO_STATUS,
	PHASE_TRAINING,
	PHASE_TESTING,
	
};


/*-------------Do not change between dashed lines----------------*/

typedef	float	FEATURE_TYPE;

#define	ROOT_PATH				"../../../"


// datasets to test
enum DATASET_TYPES
{
	VIPER = 1, CUHK01,
};
#define	DATASET_TYPE			VIPER				// VIPER, CUHK01...
#define	NUM_PERSON_TRAIN		316					// number of individuals used for train

/****************************************************************/

#endif // _DEFINE_H_