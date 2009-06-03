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

#include <iostream>
#include <fstream>
#include <cstdio>

using namespace std;

#include "defs.h"
#include "funcs.h"
#include "adaptation.h"
#include "math.h"
#include "configoptions.h"
#include "individual.h"


void Adaptation::NormalizeMutateProbs() {
	const FLOAT_TYPE generalClassTot = topoWeight+modWeight+brlenWeight;
	topoMutateProb = topoWeight / generalClassTot;
	modelMutateProb = modWeight / generalClassTot; 	

#	ifdef GANESH
		const FLOAT_TYPE topoTypeTot = randNNIweight + randSPRweight + limSPRweight + randPECRweight;
		randNNIprob = randNNIweight / topoTypeTot;
		randSPRprob = randSPRweight / topoTypeTot;
		limSPRprob = limSPRweight / topoTypeTot;
		randPECRprob = randPECRweight / topoTypeTot;
#	else
		const FLOAT_TYPE topoTypeTot = randNNIweight + randSPRweight + limSPRweight;
		randNNIprob = randNNIweight / topoTypeTot;
		randSPRprob = randSPRweight / topoTypeTot;
		limSPRprob = limSPRweight / topoTypeTot;
#	endif
}

Adaptation::Adaptation(const GeneralGamlConfig *gc) {
	intervalsToStore = gc->intervalsToStore;
	intervalLength = gc->intervalLength;
	startOptPrecision = branchOptPrecision = gc->startOptPrec;
	minOptPrecision = gc->minOptPrec;

	numPrecReductions = gc->numPrecReductions;
	// old geometric precReduction:
	//precReductionFactor = pow((minOptPrecision/startOptPrecision), 1.0/numPrecReductions);

	precReductionFactor = (gc->numPrecReductions > 0 ? (startOptPrecision - minOptPrecision)/FLOAT_TYPE(numPrecReductions) : gc->precReductionFactor);

	reset = false;

	topoWeight = gc->topoWeight;
	modWeight = gc->modWeight;
	brlenWeight = gc->brlenWeight;
	randNNIweight = gc->randNNIweight;
	origRandNNIweight = randNNIweight;
	randSPRweight = gc->randSPRweight;
	limSPRweight = gc->limSPRweight;
	limSPRrange = gc->limSPRrange;
#	ifdef GANESH
		randPECRweight = gc->randPECRweight;
#	endif

	NormalizeMutateProbs();

	exNNIprob = 0.0;
	exlimSPRprob = 0.0;

	lastgenscore = 0.0;
	laststepscore = 0.0;
	improveOverStoredIntervals = 0.0;
	recTopImproveSize = 1.0;

	randNNI.assign(intervalsToStore, 0.0);
	randNNInum.assign(intervalsToStore, 0);
	exNNI.assign(intervalsToStore, 0.0);
	exNNInum.assign(intervalsToStore, 0);
	randSPR.assign(intervalsToStore, 0.0);
	randSPRnum.assign(intervalsToStore, 0);
#	ifdef GANESH
		randPECR.assign(intervalsToStore, 0.0);
		randPECRnum.assign(intervalsToStore, 0);
#	endif
	limSPR.assign(intervalsToStore, 0.0);
	limSPRnum.assign(intervalsToStore, 0);
	exlimSPR.assign(intervalsToStore, 0.0);
	exlimSPRnum.assign(intervalsToStore, 0);
	randRecom.assign(intervalsToStore, 0.0);
	randRecomnum.assign(intervalsToStore, 0);
	bipartRecom.assign(intervalsToStore, 0.0);
	bipartRecomnum.assign(intervalsToStore, 0);
	onlyBrlen.assign(intervalsToStore, 0.0);
	onlyBrlennum.assign(intervalsToStore, 0);
	improvetotal.assign(intervalsToStore, 0.0);
	anyModel.assign(intervalsToStore, 0.0);
	anyModelnum.assign(intervalsToStore, 0);

#	ifdef MPI_VERSION
		bestFromRemote.assign(intervalsToStore, 0.0);
		bestFromRemoteNum.assign(intervalsToStore, 0);
#	endif
}


void Adaptation::SetChangeableVariablesFromConfAfterReadingCheckpoint(const GeneralGamlConfig *gc) {
	//this just sets a few adaptation members to their values in the conf file
	//potentially overriding some values read from a checkpoint.  This would
	//be useful if multiple runs with different settings were to be restarted from
	//a single checkpoint
	minOptPrecision = gc->minOptPrec;

	numPrecReductions = gc->numPrecReductions;
	if (gc->numPrecReductions > 0)
		precReductionFactor = (gc->startOptPrec- minOptPrecision)/FLOAT_TYPE(numPrecReductions);
	else
		precReductionFactor = gc->precReductionFactor;

	topoWeight = gc->topoWeight;
	modWeight = gc->modWeight;
	brlenWeight = gc->brlenWeight;
	
	origRandNNIweight = randNNIweight = gc->randNNIweight;
	randSPRweight = gc->randSPRweight;
	limSPRweight = gc->limSPRweight;
	limSPRrange = gc->limSPRrange;
	}


void Adaptation::WriteToCheckpoint(OUTPUT_CLASS &out) const{
	//this function assumes that it has been passed an MFILE that is already open for 
	//binary writing

	//7/13/07 changing this to calculate the actual size of the chunk of scalars
	//(the number of bytes between the start of the object and the first nonscalar
	//data member) rather than counting the number of each type and adding it up 
	//manually.  This should make it work irrespective of things like memory padding
	//for data member alignment, which could vary between platforms and compilers
	intptr_t scalarSize = (intptr_t) &improvetotal - (intptr_t) this;
	out.WRITE_TO_FILE(this, scalarSize, 1);

	//now the arrays, which should be of length intervalsToStore
	out.WRITE_TO_FILE(&improvetotal[0], sizeof(FLOAT_TYPE), intervalsToStore);

	out.WRITE_TO_FILE(&randNNI[0], sizeof(FLOAT_TYPE), intervalsToStore);
	out.WRITE_TO_FILE(&randNNInum[0], sizeof(int), intervalsToStore);
	
	out.WRITE_TO_FILE(&exNNI[0], sizeof(FLOAT_TYPE), intervalsToStore);
	out.WRITE_TO_FILE(&exNNInum[0], sizeof(int), intervalsToStore);

	out.WRITE_TO_FILE(&randSPR[0], sizeof(FLOAT_TYPE), intervalsToStore);
	out.WRITE_TO_FILE(&randSPRnum[0], sizeof(int), intervalsToStore);
	
	out.WRITE_TO_FILE(&limSPR[0], sizeof(FLOAT_TYPE), intervalsToStore);
	out.WRITE_TO_FILE(&limSPRnum[0], sizeof(int), intervalsToStore);

	out.WRITE_TO_FILE(&exlimSPR[0], sizeof(FLOAT_TYPE), intervalsToStore);
	out.WRITE_TO_FILE(&exlimSPRnum[0], sizeof(int), intervalsToStore);

	out.WRITE_TO_FILE(&randRecom[0], sizeof(FLOAT_TYPE), intervalsToStore);
	out.WRITE_TO_FILE(&randRecomnum[0], sizeof(int), intervalsToStore);
	
	out.WRITE_TO_FILE(&bipartRecom[0], sizeof(FLOAT_TYPE), intervalsToStore);
	out.WRITE_TO_FILE(&bipartRecomnum[0], sizeof(int), intervalsToStore);

	out.WRITE_TO_FILE(&onlyBrlen[0], sizeof(FLOAT_TYPE), intervalsToStore);
	out.WRITE_TO_FILE(&onlyBrlennum[0], sizeof(int), intervalsToStore);

	out.WRITE_TO_FILE(&anyModel[0], sizeof(FLOAT_TYPE), intervalsToStore);
	out.WRITE_TO_FILE(&anyModelnum[0], sizeof(int), intervalsToStore);
	}

void Adaptation::ReadFromCheckpoint(FILE *in) {
	//this function assumes that it has been passed a FILE* that is already open for 
	//binary reading

	//7/13/07 changing this to calculate the actual size of the chunk of scalars
	//(the number of bytes between the start of the object and the first nonscalar
	//data member) rather than counting the number of each type and adding it up 
	//manually.  This should make it work irrespective of things like memory padding
	//for data member alignment, which could vary between platforms and compilers
	intptr_t scalarSize = (intptr_t) &improvetotal - (intptr_t) this;

	fread((char *) this, 1, scalarSize, in);

	//now the arrays, which should be of length intervalsToStore
	fread((char *) &improvetotal[0], sizeof(FLOAT_TYPE), intervalsToStore, in);

	fread((char *) &randNNI[0], sizeof(FLOAT_TYPE), intervalsToStore, in);
	fread((char *) &randNNInum[0], sizeof(int), intervalsToStore, in);
	
	fread((char *) &exNNI[0], sizeof(FLOAT_TYPE), intervalsToStore, in);
	fread((char *) &exNNInum[0], sizeof(int), intervalsToStore, in);

	fread((char *) &randSPR[0], sizeof(FLOAT_TYPE), intervalsToStore, in);
	fread((char *) &randSPRnum[0], sizeof(int), intervalsToStore, in);
	
	fread((char *) &limSPR[0], sizeof(FLOAT_TYPE), intervalsToStore, in);
	fread((char *) &limSPRnum[0], sizeof(int), intervalsToStore, in);

	fread((char *) &exlimSPR[0], sizeof(FLOAT_TYPE), intervalsToStore, in);
	fread((char *) &exlimSPRnum[0], sizeof(int), intervalsToStore, in);

	fread((char *) &randRecom[0], sizeof(FLOAT_TYPE), intervalsToStore, in);
	fread((char *) &randRecomnum[0], sizeof(int), intervalsToStore, in);
	
	fread((char *) &bipartRecom[0], sizeof(FLOAT_TYPE), intervalsToStore, in);
	fread((char *) &bipartRecomnum[0], sizeof(int), intervalsToStore, in);

	fread((char *) &onlyBrlen[0], sizeof(FLOAT_TYPE), intervalsToStore, in);
	fread((char *) &onlyBrlennum[0], sizeof(int), intervalsToStore, in);

	fread((char *) &anyModel[0], sizeof(FLOAT_TYPE), intervalsToStore, in);
	fread((char *) &anyModelnum[0], sizeof(int), intervalsToStore, in);
	}


void Adaptation::PrepareForNextInterval() {
 //if we're on the first generation of a new recording period, shift everything over
	for(int i = intervalsToStore-1; i > 0; i--) {
		improvetotal[i] = improvetotal[i-1];
		randNNI[i] = randNNI[i-1];
		randNNInum[i] = randNNInum[i-1];
		exNNI[i] = exNNI[i-1];
		exNNInum[i] = exNNInum[i-1];
		 
		randSPR[i] = randSPR[i-1];
		randSPRnum[i] = randSPRnum[i-1];
		limSPR[i] = limSPR[i-1];
		limSPRnum[i] = limSPRnum[i-1];
		exlimSPR[i] = exlimSPR[i-1];
		exlimSPRnum[i] = exlimSPRnum[i-1];
#		ifdef GANESH
			randPECR[i] = randPECR[i-1];
			randPECRnum[i] = randPECRnum[i-1];			
#		endif			
		randRecom[i] = randRecom[i-1];
		randRecomnum[i] = randRecomnum[i-1];
		bipartRecom[i]=bipartRecom[i-1];
		bipartRecomnum[i] = bipartRecomnum[i-1];

		onlyBrlen[i] = onlyBrlen[i-1];
		onlyBrlennum[i] = onlyBrlennum[i-1];
		anyModel[i] = anyModel[i-1];
		anyModelnum[i] = anyModelnum[i-1];
#		ifdef MPI_VERSION
			bestFromRemote[i] = bestFromRemote[i-1];
			bestFromRemoteNum[i] = bestFromRemoteNum[i-1];
#		endif
	}

	// clean up for next entry
	improvetotal[0] = 0.0;
	randNNI[0] = 0.0;
	randNNInum[0] = 0;
	exNNI[0] = 0.0;
	exNNInum[0] = 0;
		 
	randSPR[0] = 0.0;
	randSPRnum[0] = 0;
	limSPR[0] = 0.0;
	limSPRnum[0] = 0;
	exlimSPR[0] = 0.0;
	exlimSPRnum[0] = 0;
#	ifdef GANESH
		randPECR[0] = 0.0;
		randPECRnum[0] = 0; 
#	endif
	randRecom[0] = 0.0;
	randRecomnum[0] = 0;
	bipartRecom[0]=0.0;
	bipartRecomnum[0]=0;

	onlyBrlen[0] = 0.0;
	onlyBrlennum[0] = 0;
	anyModel[0] = 0.0;
	anyModelnum[0] = 0;
#	ifdef MPI_VERSION
		bestFromRemote[0]=0.0;bestFromRemoteNum[0]=0;
#	endif
	}

void Adaptation::BeginProbLog(std::ostream &plog, int gen) {
	plog << "gen\tmod\ttopo\tbrlen\tNNI\trandSPR\tlimSPR\n";
	//gen could be non-zero if we restarted
	OutputProbs(plog, gen);
	}

void Adaptation::OutputProbs(std::ostream &plog, int gen) {
	plog << gen << "\t" << modelMutateProb << "\t" << topoMutateProb << "\t" << (1.0-modelMutateProb-topoMutateProb) << "\t";
	plog << randNNIprob << "\t" << randSPRprob << "\t" << limSPRprob << endl;
	}

void Adaptation::UpdateProbs() {
	FLOAT_TYPE topoTot = 0.0;
	FLOAT_TYPE modTot = 0.0;
	FLOAT_TYPE onlyBrlenTot = 0.0;
	int numTopos = 0, numMod = 0, numOnlyBrlen = 0;
	
	FLOAT_TYPE totRandNNI = 0.0, totLimSPR = 0.0, totRandSPR = 0.0;
	int totNumRandNNI = 0, totNumLimSPR = 0, totNumRandSPR = 0;
	FLOAT_TYPE totBipartRecom = 0.0;
	int totNumBipartRecom = 0;

#	ifdef MPI_VERSION
		FLOAT_TYPE totalFromRemote;
#	endif


#	ifdef GANESH
		FLOAT_TYPE totRandPECR = 0.0;
		int totNumRandPECR = 0;
#	endif
	
	for(unsigned i = 0;i < intervalsToStore;i++) {
#		ifdef GANESH
			topoTot += randNNI[i] + randSPR[i] + limSPR[i] + randPECR[i];
			numTopos += randNNInum[i] + randSPRnum[i] + limSPRnum[i] +
			randPECRnum[i];
#		else
			topoTot += randNNI[i] + randSPR[i] + limSPR[i];
			numTopos += randNNInum[i] + randSPRnum[i] + limSPRnum[i];
#		endif
		modTot += anyModel[i];
		numMod += anyModelnum[i];
		onlyBrlenTot += onlyBrlen[i];
		numOnlyBrlen += onlyBrlennum[i];
		totRandNNI += randNNI[i];
		totNumRandNNI += randNNInum[i];
		totLimSPR += limSPR[i];
		totNumLimSPR += limSPRnum[i];
		totRandSPR += randSPR[i];
		totBipartRecom += bipartRecom[i];
		totNumBipartRecom += bipartRecomnum[i];
		
#		ifdef GANESH
			totRandPECR += randPECR[i];
			totNumRandPECR += randPECRnum[i];
#		endif
		totNumRandSPR += randSPRnum[i];	
#		ifdef	MPI_VERSION
			totalFromRemote += bestFromRemote[i];
#		endif
	}
	
	//version 0.95b3 - The reduction of precision that used to appear here has been
	//moved to Adaptation::ReducePrecision, which is called from Run, MasterMaster and
	//RemoteSubtreeWorker when lastTopoImprove is > that #int * intLength generations ago
	const FLOAT_TYPE perTopo = topoWeight + (numTopos!=0 ? topoTot/numTopos : 0.0);
	const FLOAT_TYPE perModel = modWeight + (numMod != 0 ? modTot/numMod : 0.0);
	const FLOAT_TYPE perBrlen = brlenWeight + (numOnlyBrlen > 0 ? onlyBrlenTot/numOnlyBrlen : 0.0);
	const FLOAT_TYPE perBipartRecom (totNumBipartRecom > 0 ? totBipartRecom/totNumBipartRecom : 0.0);
	const FLOAT_TYPE generalClassTot = perTopo + perModel + perBrlen;
	
	FLOAT_TYPE brlenOnlyMut;

	//only update these probs if model mutations are turned off completely
	//or if some model mutations have been done (ie not in subtree mode)
	if (anyModelnum[0] != 0 || FloatingPointEquals(modWeight, 0.0, 1e-10)) {
		brlenOnlyMut = perBrlen/generalClassTot;
		modelMutateProb = perModel/generalClassTot;
		topoMutateProb = perTopo/generalClassTot;
	}

	//enforce a minimum probability
	if (modWeight != 0.0 && topoWeight != 0.0) {
		FLOAT_TYPE minProb= (FLOAT_TYPE) 0.02;	
		if (topoMutateProb < minProb) {
			modelMutateProb -= minProb - topoMutateProb;
			topoMutateProb = minProb;
			}
		if (modelMutateProb < minProb) {
			topoMutateProb -= minProb - modelMutateProb;
			modelMutateProb = minProb;
			}
		if (1.0 - (modelMutateProb + topoMutateProb) < minProb) {
			const FLOAT_TYPE diff = minProb - (FLOAT_TYPE)(1.0 - (modelMutateProb + topoMutateProb));
			if (modelMutateProb - diff/2.0 > .02 && topoMutateProb - diff/2.0 > .02) {
				modelMutateProb -= diff/(FLOAT_TYPE)2.0;
				topoMutateProb -= diff/(FLOAT_TYPE)2.0;
				}
			else{
				if (modelMutateProb - diff/2.0 < .02)
					topoMutateProb -= diff;
				else 
					modelMutateProb -= diff;
				}		

			brlenOnlyMut = minProb;
			}
		}
		
	if (totNumRandNNI == 0)
		totNumRandNNI = 1;
	if (totNumLimSPR == 0)
		totNumLimSPR = 1;
	if (totNumRandSPR == 0)
		totNumRandSPR = 1;
#	ifdef GANESH
		if (totNumRandPECR == 0)
			totNumRandPECR = 1;
#	endif
	//Because NNI's chosen by an SPR mutator are marked as NNI's, this needs to be done to keep from 
	//giving NNI's some prob even when the weight was 0.0
	FLOAT_TYPE perRandNNI= (randNNIweight == ZERO_POINT_ZERO ? ZERO_POINT_ZERO : totRandNNI/totNumRandNNI + randNNIweight);
	FLOAT_TYPE perLimSPR=  (limSPRweight  == ZERO_POINT_ZERO ? ZERO_POINT_ZERO : totLimSPR/totNumLimSPR + limSPRweight);
	FLOAT_TYPE perRandSPR= (limSPRweight  == ZERO_POINT_ZERO ? ZERO_POINT_ZERO : totRandSPR/totNumRandSPR + randSPRweight);
	
#	ifdef GANESH
		FLOAT_TYPE perRandPECR = totRandPECR/totNumRandPECR + randPECRweight;
		const FLOAT_TYPE typeOfTopoTot = perRandNNI+perLimSPR+perRandSPR+perRandPECR;
		randNNIprob = perRandNNI/typeOfTopoTot;
		randSPRprob = perRandSPR/typeOfTopoTot;
		limSPRprob = perLimSPR/typeOfTopoTot;
		randPECRprob = perRandPECR/typeOfTopoTot;
#	else
		const FLOAT_TYPE typeOfTopoTot = perRandNNI+perLimSPR+perRandSPR;
		if (typeOfTopoTot > ZERO_POINT_ZERO) {
			randNNIprob = perRandNNI/typeOfTopoTot;
			randSPRprob = perRandSPR/typeOfTopoTot;
			limSPRprob = perLimSPR/typeOfTopoTot;
		}
#	endif
}
