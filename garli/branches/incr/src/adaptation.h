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


#ifndef ADAPTION_H
#define ADAPTION_H

#include <fstream>
#include <vector>

#include "configoptions.h"
#include "hashdefines.h"

class MFILE;

class Adaptation{
	public:
	//here are all of the scalars:
	//4 ints, 1 bool, 23 FLOAT_TYPEs
	unsigned intervalLength; 
	unsigned intervalsToStore;
	FLOAT_TYPE lastgenscore;
	FLOAT_TYPE laststepscore;

	FLOAT_TYPE improveOverStoredIntervals;
	
	bool reset;

	FLOAT_TYPE startOptPrecision;
	FLOAT_TYPE branchOptPrecision;
	FLOAT_TYPE minOptPrecision;
	FLOAT_TYPE precReductionFactor;
	int numPrecReductions;
	
	FLOAT_TYPE topoWeight;
	FLOAT_TYPE modWeight;
	FLOAT_TYPE brlenWeight;
	
	FLOAT_TYPE randNNIweight;
	FLOAT_TYPE origRandNNIweight;
	FLOAT_TYPE randSPRweight;
	FLOAT_TYPE limSPRweight;
	
	FLOAT_TYPE recTopImproveSize;
	FLOAT_TYPE exNNIprob;
	FLOAT_TYPE exlimSPRprob;

	FLOAT_TYPE topoMutateProb;
	FLOAT_TYPE randNNIprob;
	FLOAT_TYPE randSPRprob;
	FLOAT_TYPE limSPRprob;
	
	FLOAT_TYPE modelMutateProb;

	unsigned limSPRrange;

	//the vectors. All will be of length intervalsToStore
	std::vector<FLOAT_TYPE> improvetotal;

	std::vector<FLOAT_TYPE> randNNI;
	std::vector<int> randNNInum;
	
	std::vector<FLOAT_TYPE> exNNI;
	std::vector<int> exNNInum;

	std::vector<FLOAT_TYPE> randSPR;
	std::vector<int> randSPRnum;
	
	std::vector<FLOAT_TYPE> limSPR;
	std::vector<int> limSPRnum;

	std::vector<FLOAT_TYPE> exlimSPR;
	std::vector<int> exlimSPRnum;

	std::vector<FLOAT_TYPE> randRecom;
	std::vector<int> randRecomnum;
	
	std::vector<FLOAT_TYPE> bipartRecom;
	std::vector<int> bipartRecomnum;

	std::vector<FLOAT_TYPE> onlyBrlen;
	std::vector<int> onlyBrlennum;

	std::vector<FLOAT_TYPE> anyModel;
	std::vector<int> anyModelnum;
	
#ifdef MPI_VERSION
	std::vector<FLOAT_TYPE> fromRemoteSubtree;
	std::vector<FLOAT_TYPE> fromRemoteNonSubtree;
	std::vector<int> bestFromRemoteNum;
	std::vector<FLOAT_TYPE> bestFromRemote;
#endif

	Adaptation(const GeneralGamlConfig *gc);
	public:
	void SetChangeableVariablesFromConfAfterReadingCheckpoint(const GeneralGamlConfig *gc);
	void PrepareForNextInterval();
	void UpdateProbs();
	void OutputProbs(std::ostream &plog, int gen);
	void BeginProbLog(std::ostream &plot, int gen);
	bool ReducePrecision(){
		if(FloatingPointEquals(branchOptPrecision, minOptPrecision, 1e-10) || numPrecReductions == 0) return false;
		if(topoMutateProb > .01 || FloatingPointEquals(topoWeight, 0.0, 1e-10)){
			//changing this to a linear reduction in prec.  Geometric was too fast
			//branchOptPrecision*=precReductionFactor;
			branchOptPrecision -= precReductionFactor;
			if((branchOptPrecision < minOptPrecision) || (branchOptPrecision - minOptPrecision) < (precReductionFactor / 2.0)) branchOptPrecision=minOptPrecision;
			return true;
			}
		else return false;
		}

	void WriteToCheckpoint(OUTPUT_CLASS &) const;
	void ReadFromCheckpoint(FILE *);
	void NormalizeMutateProbs();
};

#endif
