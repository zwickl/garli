// GARLI version 0.96b8 source code
// Copyright 2005-2008 Derrick J. Zwickl
// email: zwickl@nescent.org
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef FUNCS_H
#define FUNCS_H

//a variety of functions that don't belong to any class

#include <stdlib.h>

#include "population.h"
#include "sequencedata.h"
#ifdef UNIX
#include <sys/mman.h>
#endif

extern rng rnd;

class StateSet{
	protected:
		vector<string> states;
		int numStates;
	public:
		StateSet(int ns){
			numStates = ns;
			assert(numStates == 4 || numStates == 20);
			if(numStates == 4){
				states.push_back("A");
				states.push_back("C");
				states.push_back("G");
				states.push_back("T");
				}
			else if(numStates == 20){
				states.push_back("A");
				states.push_back("C");
				states.push_back("D");
				states.push_back("E");
				states.push_back("F");
				states.push_back("G");
				states.push_back("H");
				states.push_back("I");
				states.push_back("K");
				states.push_back("L");
				states.push_back("M");
				states.push_back("N");
				states.push_back("P");
				states.push_back("Q");
				states.push_back("R");
				states.push_back("S");
				states.push_back("T");
				states.push_back("V");
				states.push_back("W");
				states.push_back("Y");
				}
			}
		StateSet(const GeneticCode *code){
			numStates = code->NumStates();
			for(int s = 0;s < numStates;s++)
				states.push_back(code->LookupCodonDisplayFromIndex(s));
			}
		void OutputInternalStateHeader(ofstream &out) const{
			out << "site\tbestState(prob)\t";
			for(int s = 0;s < numStates;s++)
				out << "prob(" << states[s] << ")\t";
			out << endl;
			}
		const string GetState(int s) const{
			return states[s];
			}
	};

class InternalState{
	protected:
		int best;
		int numStates;
		vector<FLOAT_TYPE> probs;

	public:
		InternalState(int ns){
			numStates = ns;
			probs.resize(numStates);
			}
		void CalcProbs(const FLOAT_TYPE *tots){
			FLOAT_TYPE tot=0.0;
			best = 0;
			FLOAT_TYPE bestVal = ZERO_POINT_ZERO;

			for(int s = 0;s < numStates;s++)
				tot += tots[s];
			for(int i=0;i<numStates;i++){
				probs[i]=tots[i]/tot;
				if(probs[i] > bestVal){
					bestVal = probs[i];
					best = i;
					}
				}
			}
	void Output(ofstream &out, const StateSet &states) const{
		out << states.GetState(best) << "(" << probs[best] << ")\t";
		for(int s = 0;s < numStates;s++)
			out << probs[s] << "\t";
		out << endl;
		}
	};

bool FloatingPointEquals(const FLOAT_TYPE first, const FLOAT_TYPE sec, const FLOAT_TYPE epsilon);

#if defined(SINGLE_PRECISION_FLOATS) && (!defined(_MSC_VER)) || (defined(BOINC) && defined (_WIN32))
//Overloaded versions of min and max that take different types for the two arguments
//This should not be used in hot code when possible, and conditional comp should
//be used to make two different versions of the code
float min(const double first, const float second);
float min(const float first, const double second);
float max(const double first, const float second);
float max(const float first, const double second);
#endif

void OutputImportantDefines();

int FileExists(const char* s);
bool FileIsFasta(const char *name);
bool FileIsNexus(const char *name);
int ReadData(GeneralGamlConfig *, SequenceData* data);
bool ReadData(const char* filename, SequenceData* data);
//void GetRestartParams(Parameters& params);
int RandomInt(int lb, int ub);
FLOAT_TYPE RandomFrac();
FLOAT_TYPE RandomDouble(FLOAT_TYPE lb, FLOAT_TYPE ub);
int mnbrak(FLOAT_TYPE *ax, FLOAT_TYPE *bx, FLOAT_TYPE *cx, FLOAT_TYPE *fa, FLOAT_TYPE *fb, FLOAT_TYPE *fc, FLOAT_TYPE (*func)(TreeNode*, Tree*, FLOAT_TYPE), TreeNode *thisnode, Tree *thistree);
int DZbrak(FLOAT_TYPE *ax, FLOAT_TYPE *bx, FLOAT_TYPE *cx, FLOAT_TYPE *fa, FLOAT_TYPE *fb, FLOAT_TYPE *fc, FLOAT_TYPE (*func)(TreeNode*, Tree*, FLOAT_TYPE), TreeNode *thisnode, Tree *thistree);
FLOAT_TYPE brent(FLOAT_TYPE ax, FLOAT_TYPE bx, FLOAT_TYPE cx, FLOAT_TYPE (*f)(TreeNode *, Tree*, FLOAT_TYPE), FLOAT_TYPE tol, FLOAT_TYPE *xmin, TreeNode *thisnode, Tree *thistree);
FLOAT_TYPE DZbrent(FLOAT_TYPE ax, FLOAT_TYPE bx, FLOAT_TYPE cx, FLOAT_TYPE fa, FLOAT_TYPE fb, FLOAT_TYPE fc, FLOAT_TYPE (*f)(TreeNode *, Tree*, FLOAT_TYPE), FLOAT_TYPE tol, FLOAT_TYPE *xmin, TreeNode *thisnode, Tree *thistree);
void DirichletRandomVariable (FLOAT_TYPE *alp, FLOAT_TYPE *z, int n);

void InferStatesFromCla(vector<InternalState> &stateVec, const FLOAT_TYPE *cla, int nchar, int nstates);
FLOAT_TYPE CalculateHammingDistance(const char *str1, const char *str2, const int *counts, int nchar, int nstates);

void SampleBranchLengthCurve(FLOAT_TYPE (*func)(TreeNode*, Tree*, FLOAT_TYPE, bool), TreeNode *thisnode, Tree *thistree);

int DZbrak(FLOAT_TYPE *worstOuter, FLOAT_TYPE *mid, FLOAT_TYPE *bestOuter, FLOAT_TYPE *worstOuterL, FLOAT_TYPE *midL, FLOAT_TYPE *bestOuterL, FLOAT_TYPE (*func)(TreeNode*, Tree*, FLOAT_TYPE, bool), TreeNode *thisnode, Tree *thistree);
FLOAT_TYPE DZbrent(FLOAT_TYPE ax, FLOAT_TYPE bx, FLOAT_TYPE cx, FLOAT_TYPE fa, FLOAT_TYPE fx, FLOAT_TYPE fc, FLOAT_TYPE (*f)(TreeNode *, Tree*, FLOAT_TYPE, bool), FLOAT_TYPE tol, FLOAT_TYPE *xmin, TreeNode *thisnode, Tree *thistree);
/*
void CalcFullCLAInternalInternalRateHet(FLOAT_TYPE *dest, const FLOAT_TYPE *LCL, const FLOAT_TYPE *RCL, const FLOAT_TYPE *Lpr, const FLOAT_TYPE *Rpr, int nchar);
void CalcFullCLATerminalTerminalRateHet(FLOAT_TYPE *dest, const FLOAT_TYPE *LCL, const FLOAT_TYPE *RCL, const FLOAT_TYPE *Lpr, const FLOAT_TYPE *Rpr, unsigned char *Ldata, unsigned char *Rdata, int nchar);
void CalcFullCLAInternalTerminalRateHet(FLOAT_TYPE *dest, const FLOAT_TYPE *CL1, const FLOAT_TYPE *pr1, const FLOAT_TYPE *pr2, unsigned char *data2, int nchar);
void CalcFullCLAInternalInternal(FLOAT_TYPE *dest, const FLOAT_TYPE *LCL, const FLOAT_TYPE *RCL, const FLOAT_TYPE *Lpr, const FLOAT_TYPE *Rpr, int nchar);
void CalcFullCLATerminalTerminal(FLOAT_TYPE *dest, const FLOAT_TYPE *LCL, const FLOAT_TYPE *RCL, const FLOAT_TYPE *Lpr, const FLOAT_TYPE *Rpr, unsigned char *Ldata, unsigned char *Rdata, int nchar);
void CalcFullCLAInternalTerminal(FLOAT_TYPE *dest, const FLOAT_TYPE *CL1, const FLOAT_TYPE *pr1, const FLOAT_TYPE *pr2, unsigned char *data2, int nchar);
*/
int gsl_min_find_bracket(FLOAT_TYPE (*f)(TreeNode *, Tree*, FLOAT_TYPE),FLOAT_TYPE *x_minimum,FLOAT_TYPE * f_minimum,FLOAT_TYPE * x_lower, FLOAT_TYPE * f_lower, FLOAT_TYPE * x_upper, FLOAT_TYPE * f_upper, size_t eval_max, TreeNode *thisnode, Tree *thistree);

inline void ArrayMultiply(FLOAT_TYPE *dest, const FLOAT_TYPE *source, int num){
	//simply multiplies each element in dest by the corresponding element in source, up to num
	for(register int i=0;i<num;i++)
		*(dest++) *= *(source++);
	}

inline void CalcSiteCLARateHetEquals(FLOAT_TYPE *dest, const FLOAT_TYPE *tCL, const FLOAT_TYPE *tp){
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

inline void CalcSiteCLARateHetTimes(FLOAT_TYPE *dest, const FLOAT_TYPE *tCL, const FLOAT_TYPE *tp){
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
