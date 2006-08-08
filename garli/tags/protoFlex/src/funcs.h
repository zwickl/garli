// GARLI version 0.93 source code
// Copyright  2005 by Derrick J. Zwickl
// All rights reserved.
//
// This code may be used and modified for non-commercial purposes
// but redistribution in any form requires written permission.
// Please contact:
//
//  Derrick Zwickl
//	Integrative Biology, UT
//	1 University Station, C0930
//	Austin, TX  78712
//  email: zwickl@mail.utexas.edu
//
//	Note: In 2006  moving to NESCENT (The National
//	Evolutionary Synthesis Center) for a postdoc

#ifndef FUNCS_H
#define FUNCS_H

//a variety of functions that don't belong to any class

#include <stdlib.h>

#include "parameters.h"
#include "population.h"
#include "mlhky.h"
#ifdef UNIX
#include <sys/mman.h>
#endif

extern rng rnd;

class InternalState{
	private:
	char best;
	double probs[4];
	
	public:
	InternalState(double *tots){
		char bases[4]={'A', 'C', 'G', 'T'};

		double tot=0.0;
		tot = tots[0] + tots[1] + tots[2] + tots[3];
		int max1=(tots[0] > tots[1] ? 0:1);
		int max2=(tots[2] > tots[3] ? 2:3);	
		best=bases[(tots[max1] > tots[max2] ? max1 : max2)];
		for(int i=0;i<4;i++)
			probs[i]=tots[i]/tot;
		}
	void Output(ofstream &out){
		out << best << "\t" << probs[0] << "\t" << probs[1] <<  "\t" << probs[2] << "\t" <<  probs[3] << "\n";
		
		}
	};

int FileExists(const char* s);
int ReadData(const Parameters& params, HKYData* data);
int ReadData(const char* filename, HKYData* data);
void GetRestartParams(Parameters& params);
int RandomInt(int lb, int ub);
double RandomFrac();
double RandomDouble(double lb, double ub);
int mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double (*func)(TreeNode*, Tree*, double), TreeNode *thisnode, Tree *thistree);
int DZbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double (*func)(TreeNode*, Tree*, double), TreeNode *thisnode, Tree *thistree);
double brent(double ax, double bx, double cx, double (*f)(TreeNode *, Tree*, double), double tol, double *xmin, TreeNode *thisnode, Tree *thistree);
double DZbrent(double ax, double bx, double cx, double fa, double fb, double fc, double (*f)(TreeNode *, Tree*, double), double tol, double *xmin, TreeNode *thisnode, Tree *thistree);
void DirichletRandomVariable (double *alp, double *z, int n);
void InferStatesFromCla(char *states, double *cla, int nchar);
vector<InternalState *> *InferStatesFromCla(double *cla, int nchar);
double CalculatePDistance(const char *str1, const char *str2, int nchar);
#ifndef GANESH
double CalculateHammingDistance(const char *str1, const char *str2, int nchar);
#else
double CalculateHammingDistance(const char *str1, const char *str2, const int *col_count, int nchar)
#endif
void SampleBranchLengthCurve(double (*func)(TreeNode*, Tree*, double, bool), TreeNode *thisnode, Tree *thistree);
void CalcDerivCLASPartialInternalRateHet(CondLikeArray *destCLA, CondLikeArray *destD1CLA, CondLikeArray *destD2CLA, const CondLikeArray *partialCLA, const CondLikeArray *childCLA,
	const double *prmat, const double *d1mat, const double *d2mat, int nchar);
void CalcDerivCLASPartialTerminalRateHet(CondLikeArray *destCLA, CondLikeArray *destD1CLA, CondLikeArray *destD2CLA, const CondLikeArray *partialCLA,
	const double *prmat, const double *d1mat, const double *d2mat, const char *Ldata, int nchar);

void CalcFullCLAInternalInternalRateHet(CondLikeArray *destCLA, const CondLikeArray *LCLA, const CondLikeArray *RCLA, const double *Lpr, const double *Rpr, int nchar, int nRateCats=4);
void CalcFullCLAInternalInternal(CondLikeArray *destCLA, const CondLikeArray *LCLA, const CondLikeArray *RCLA, const double *Lpr, const double *Rpr, int nchar);
void CalcFullCLATerminalTerminalRateHet(CondLikeArray *destCLA, const double *Lpr, const double *Rpr, char *Ldata, char *Rdata, int nchar, int nRateCats=4);
void CalcFullCLATerminalTerminal(CondLikeArray *destCLA, const double *Lpr, const double *Rpr, char *Ldata, char *Rdata, int nchar);
void CalcFullCLAInternalTerminalRateHet(CondLikeArray *destCLA, const CondLikeArray *LCLA, const double *pr1, const double *pr2, char *data2, int nchar, int nRateCats=4);
void CalcFullCLAInternalTerminal(CondLikeArray *destCLA, const CondLikeArray *LCLA, const double *pr1, const double *pr2, char *data2, int nchar);
void CalcFullCLAPartialInternalRateHet(CondLikeArray *destCLA, const CondLikeArray *LCLA, const double *pr1, CondLikeArray *partialCLA, int nchar, int nRateCats=4);
void CalcFullCLAPartialInternal(CondLikeArray *destCLA, const CondLikeArray *LCLA, const double *pr1, CondLikeArray *partialCLA, int nchar);
void CalcFullCLAPartialTerminalRateHet(CondLikeArray *destCLA, const CondLikeArray *partialCLA, const double *Lpr, char *Ldata, int nchar, int nRateCats=4);
void CalcFullCLAPartialTerminal(CondLikeArray *destCLA, const CondLikeArray *partialCLA, const double *Lpr, char *Ldata, int nchar);

int DZbrak(double *worstOuter, double *mid, double *bestOuter, double *worstOuterL, double *midL, double *bestOuterL, double (*func)(TreeNode*, Tree*, double, bool), TreeNode *thisnode, Tree *thistree);
double DZbrent(double ax, double bx, double cx, double fa, double fx, double fc, double (*f)(TreeNode *, Tree*, double, bool), double tol, double *xmin, TreeNode *thisnode, Tree *thistree);
/*
void CalcFullCLAInternalInternalRateHet(double *dest, const double *LCL, const double *RCL, const double *Lpr, const double *Rpr, int nchar);
void CalcFullCLATerminalTerminalRateHet(double *dest, const double *LCL, const double *RCL, const double *Lpr, const double *Rpr, unsigned char *Ldata, unsigned char *Rdata, int nchar);
void CalcFullCLAInternalTerminalRateHet(double *dest, const double *CL1, const double *pr1, const double *pr2, unsigned char *data2, int nchar);
void CalcFullCLAInternalInternal(double *dest, const double *LCL, const double *RCL, const double *Lpr, const double *Rpr, int nchar);
void CalcFullCLATerminalTerminal(double *dest, const double *LCL, const double *RCL, const double *Lpr, const double *Rpr, unsigned char *Ldata, unsigned char *Rdata, int nchar);
void CalcFullCLAInternalTerminal(double *dest, const double *CL1, const double *pr1, const double *pr2, unsigned char *data2, int nchar);
*/
int gsl_min_find_bracket(double (*f)(TreeNode *, Tree*, double),double *x_minimum,double * f_minimum,double * x_lower, double * f_lower, double * x_upper, double * f_upper, size_t eval_max, TreeNode *thisnode, Tree *thistree);

inline void ArrayMultiply(double *dest, const double *source, int num){
	//simply multiplies each element in dest by the corresponding element in source, up to num
	for(register int i=0;i<num;i++)
		*(dest++) *= *(source++);
	}
	
inline void CalcSiteCLARateHetEquals(double *dest, const double *tCL, const double *tp){
	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.
	for(int i=0;i<4;i++){
		*(dest++)=tp[0]*tCL[0]+tp[1]*tCL[1]+tp[2]*tCL[2]+tp[3]*tCL[3];
		*(dest++)=tp[4]*tCL[0]+tp[5]*tCL[1]+tp[6]*tCL[2]+tp[7]*tCL[3];
		*(dest++)=tp[8]*tCL[0]+tp[9]*tCL[1]+tp[10]*tCL[2]+tp[11]*tCL[3];
		*(dest++)=tp[12]*tCL[0]+tp[13]*tCL[1]+tp[14]*tCL[2]+tp[15]*tCL[3];
		tp+=16;
		tCL+=4;
		}
	}

inline void CalcSiteCLARateHetTimes(double *dest, const double *tCL, const double *tp){
	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.
	for(int i=0;i<4;i++){
		*(dest++)*=tp[0]*tCL[0]+tp[1]*tCL[1]+tp[2]*tCL[2]+tp[3]*tCL[3];
		*(dest++)*=tp[4]*tCL[0]+tp[5]*tCL[1]+tp[6]*tCL[2]+tp[7]*tCL[3];
		*(dest++)*=tp[8]*tCL[0]+tp[9]*tCL[1]+tp[10]*tCL[2]+tp[11]*tCL[3];
		*(dest++)*=tp[12]*tCL[0]+tp[13]*tCL[1]+tp[14]*tCL[2]+tp[15]*tCL[3];
		tp+=16;
		tCL+=4;
		}
					}
			 



template<class T>
void ScrambleArray(int n, T ar[])	{
	int times = n*2;

	int x, y;
	T temp;
	for (int i = 0; i < times; ++i)	{
//		x = rand() % n;
//		y = rand() % n;
		x = rnd.random_int(RAND_MAX) % n;
		y = rnd.random_int(RAND_MAX) % n;

		if (x != y)	{
			temp = ar[x];
			ar[x] = ar[y];
			ar[y] = temp;
		}
	}

};

#endif
