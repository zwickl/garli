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

extern int  InvertMatrix (FLOAT_TYPE **a, int n, FLOAT_TYPE *col, int *indx, FLOAT_TYPE **a_inv);
extern int  LUDecompose (FLOAT_TYPE **a, int n, FLOAT_TYPE *vv, int *indx, FLOAT_TYPE *pd);
int  EigenRealGeneral (int n, FLOAT_TYPE **a, FLOAT_TYPE *v, FLOAT_TYPE *vi, FLOAT_TYPE **u, int *iwork, FLOAT_TYPE *work);


//these are actually from John's MCMC.c file
void CalcCijk (FLOAT_TYPE *c_ijk, int n, const FLOAT_TYPE **u,  const FLOAT_TYPE **v);
void CalcPij (const FLOAT_TYPE *c_ijk, int n, const FLOAT_TYPE *eigenValues, FLOAT_TYPE v, FLOAT_TYPE r, FLOAT_TYPE **p, FLOAT_TYPE *EigValexp);

