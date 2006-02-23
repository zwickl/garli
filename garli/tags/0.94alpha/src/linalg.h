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

extern int  InvertMatrix (double **a, int n, double *col, int *indx, double **a_inv);
extern int  LUDecompose (double **a, int n, double *vv, int *indx, double *pd);
int  EigenRealGeneral (int n, double **a, double *v, double *vi, double **u, int *iwork, double *work);


//these are actually from John's MCMC.c file
void CalcCijk (double *c_ijk, int n, const double **u,  const double **v);
void CalcPij (const double *c_ijk, int n, const double *eigenValues, double v, double r, double **p, double *EigValexp);

