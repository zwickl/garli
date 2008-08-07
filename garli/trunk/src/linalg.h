/*	linalg.h
|
|	Prototypes for matrix-inversion and eigensystem functions
|
|	Copyright (c) 1998 by David L. Swofford, Smithsonian Institution.
|	All rights reserved.
|
|	NOTE: if ANSI function prototypes are not supported, define NO_PROTOTYPES
|		  before including this file.
*/

#define RC_COMPLEX_EVAL 2	/* code that complex eigenvalue obtained */

extern int  InvertMatrix (MODEL_FLOAT **a, int n, MODEL_FLOAT *col, int *indx, MODEL_FLOAT **a_inv);
extern int  LUDecompose (MODEL_FLOAT **a, int n, MODEL_FLOAT *vv, int *indx, MODEL_FLOAT *pd);
int  EigenRealGeneral (int n, MODEL_FLOAT **a, MODEL_FLOAT *v, MODEL_FLOAT *vi, MODEL_FLOAT **u, int *iwork, MODEL_FLOAT *work);


//these are actually from John's MCMC.c file
void CalcCijk (MODEL_FLOAT *c_ijk, int n, const MODEL_FLOAT **u,  const MODEL_FLOAT **v);
void CalcPij (const MODEL_FLOAT *c_ijk, int n, const MODEL_FLOAT *eigenValues, MODEL_FLOAT v, MODEL_FLOAT r, MODEL_FLOAT **p, MODEL_FLOAT *EigValexp);

