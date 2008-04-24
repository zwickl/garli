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

Adaptation::Adaptation(const GeneralGamlConfig *gc){

	intervalsToStore=gc->intervalsToStore;
	intervalLength=gc->intervalLength;

	startOptPrecision = branchOptPrecision = gc->startOptPrec;
	minOptPrecision = gc->minOptPrec;

	numPrecReductions=gc->numPrecReductions;
	if(gc->numPrecReductions > 0)//changing prec reduction to linear rather than geometric
		//precReductionFactor = pow((minOptPrecision/startOptPrecision), 1.0/numPrecReductions);
		precReductionFactor = (startOptPrecision - minOptPrecision)/FLOAT_TYPE(numPrecReductions);
	else
		precReductionFactor = gc->precReductionFactor;

	reset=false;

	topoWeight = gc->topoWeight;
	modWeight = gc->modWeight;
	brlenWeight = gc->brlenWeight;
	
	origRandNNIweight = randNNIweight = gc->randNNIweight;
	randSPRweight = gc->randSPRweight;
	limSPRweight = gc->limSPRweight;
	limSPRrange = gc->limSPRrange;
#ifdef GANESH
	randPECRweight = gc->randPECRweight;
#endif

	FLOAT_TYPE tot = topoWeight+modWeight+brlenWeight;
	topoMutateProb = topoWeight / tot;
	modelMutateProb = modWeight / tot; 	

#ifdef GANESH
	tot = randNNIweight + randSPRweight + limSPRweight + randPECRweight;
	randNNIprob = randNNIweight / tot;
	randSPRprob = randSPRweight / tot;
	limSPRprob = limSPRweight / tot;
	randPECRprob = randPECRweight / tot;
#else
	tot = randNNIweight + randSPRweight + limSPRweight;
	randNNIprob = randNNIweight / tot;
	randSPRprob = randSPRweight / tot;
	limSPRprob = limSPRweight / tot;
#endif

	exNNIprob = 0.0;
	exlimSPRprob = 0.0;

	lastgenscore = 0.0;
	laststepscore=0.0;
	improveOverStoredIntervals=0.0;
	recTopImproveSize = 1.0;

	randNNI = new FLOAT_TYPE[intervalsToStore]; randNNInum = new int[intervalsToStore];
	exNNI = new FLOAT_TYPE[intervalsToStore]; exNNInum = new int[intervalsToStore];
	randSPR = new FLOAT_TYPE[intervalsToStore]; randSPRnum = new int[intervalsToStore];
#ifdef GANESH
	randPECR = new FLOAT_TYPE[intervalsToStore]; randPECRnum = new int[intervalsToStore];
#endif
	limSPR = new FLOAT_TYPE[intervalsToStore]; limSPRnum = new int[intervalsToStore];
	exlimSPR = new FLOAT_TYPE[intervalsToStore]; exlimSPRnum = new int[intervalsToStore];
	randRecom = new FLOAT_TYPE[intervalsToStore];	 randRecomnum = new int[intervalsToStore];
	bipartRecom = new FLOAT_TYPE[intervalsToStore];	 bipartRecomnum = new int[intervalsToStore];
	onlyBrlen = new FLOAT_TYPE[intervalsToStore]; onlyBrlennum = new int[intervalsToStore];
	improvetotal = new FLOAT_TYPE[intervalsToStore];
	anyModel = new FLOAT_TYPE[intervalsToStore]; anyModelnum = new int[intervalsToStore];

#ifdef MPI_VERSION
	bestFromRemote=new FLOAT_TYPE[intervalsToStore];bestFromRemoteNum=new int[intervalsToStore];
#endif
	for(unsigned i=0;i<intervalsToStore;i++){
	randNNI[i] = 0.0; randNNInum[i] = 0;
	exNNI[i] = 0.0; exNNInum[i] = 0;
	randSPR[i] = 0.0; randSPRnum[i] = 0;
#ifdef GANESH
	randPECR[i] = 0.0; randPECRnum[i] = 0;
#endif
	limSPR[i] = 0.0; limSPRnum[i] = 0;
	exlimSPR[i] = 0.0; exlimSPRnum[i] = 0;
	randRecom[i] = 0.0;	 randRecomnum[i] = 0;
	bipartRecom[i]= 0.0;	 bipartRecomnum[i]= 0;
	onlyBrlen[i] = 0.0; onlyBrlennum[i] = 0;
	improvetotal[i] = 0.0;
	anyModel[i] = 0.0; anyModelnum[i] = 0;
#ifdef MPI_VERSION
	bestFromRemote[i]=0.0;	bestFromRemoteNum[i]=0;
#endif

	}
}

Adaptation::~Adaptation(){
	delete []randNNI; delete []randNNInum;
	delete []exNNI; delete []exNNInum;
	delete []randSPR; delete []randSPRnum ;
#ifdef GANESH
	delete []randPECR ; delete []randPECRnum ;
#endif
	delete []limSPR ; delete []limSPRnum ;
	delete []exlimSPR ; delete []exlimSPRnum ;
	delete []randRecom ;	 delete []randRecomnum ;
	delete []bipartRecom ;	 delete []bipartRecomnum ;
	delete []onlyBrlen ; delete []onlyBrlennum ;
	delete []improvetotal ;
	delete []anyModel ; delete []anyModelnum ;
	}

void Adaptation::SetChangeableVariablesFromConfAfterReadingCheckpoint(const GeneralGamlConfig *gc){
	//this just sets a few adaptation members to their values in the conf file
	//potentially overriding some values read from a checkpoint.  This would
	//be useful if multiple runs with different settings were to be restarted from
	//a single checkpoint
	minOptPrecision = gc->minOptPrec;

	numPrecReductions=gc->numPrecReductions;
	if(gc->numPrecReductions > 0)
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
	out.WRITE_TO_FILE(improvetotal, sizeof(FLOAT_TYPE), intervalsToStore);

	out.WRITE_TO_FILE(randNNI, sizeof(FLOAT_TYPE), intervalsToStore);
	out.WRITE_TO_FILE(randNNInum, sizeof(int), intervalsToStore);
	
	out.WRITE_TO_FILE(exNNI, sizeof(FLOAT_TYPE), intervalsToStore);
	out.WRITE_TO_FILE(exNNInum, sizeof(int), intervalsToStore);

	out.WRITE_TO_FILE(randSPR, sizeof(FLOAT_TYPE), intervalsToStore);
	out.WRITE_TO_FILE(randSPRnum, sizeof(int), intervalsToStore);
	
	out.WRITE_TO_FILE(limSPR, sizeof(FLOAT_TYPE), intervalsToStore);
	out.WRITE_TO_FILE(limSPRnum, sizeof(int), intervalsToStore);

	out.WRITE_TO_FILE(exlimSPR, sizeof(FLOAT_TYPE), intervalsToStore);
	out.WRITE_TO_FILE(exlimSPRnum, sizeof(int), intervalsToStore);

	out.WRITE_TO_FILE(randRecom, sizeof(FLOAT_TYPE), intervalsToStore);
	out.WRITE_TO_FILE(randRecomnum, sizeof(int), intervalsToStore);
	
	out.WRITE_TO_FILE(bipartRecom, sizeof(FLOAT_TYPE), intervalsToStore);
	out.WRITE_TO_FILE(bipartRecomnum, sizeof(int), intervalsToStore);

	out.WRITE_TO_FILE(onlyBrlen, sizeof(FLOAT_TYPE), intervalsToStore);
	out.WRITE_TO_FILE(onlyBrlennum, sizeof(int), intervalsToStore);

	out.WRITE_TO_FILE(anyModel, sizeof(FLOAT_TYPE), intervalsToStore);
	out.WRITE_TO_FILE(anyModelnum, sizeof(int), intervalsToStore);
	}

/*
void Adaptation::WriteToCheckpoint(ofstream &out) const{
	//this function assumes that it has been passed a stream that is already open for 
	//binary writing
	assert(out.good());

	//7/13/07 changing this to calculate the actual size of the chunk of scalars
	//(the number of bytes between the start of the object and the first nonscalar
	//data member) rather than counting the number of each type and adding it up 
	//manually.  This should make it work irrespective of things like memory padding
	//for data member alignment, which could vary between platforms and compilers
	intptr_t scalarSize = (intptr_t) &improvetotal - (intptr_t) this;

	out.write((char *) this, (streamsize) scalarSize);

	//now the arrays, which should be of length intervalsToStore
	out.write((char *) improvetotal, sizeof(FLOAT_TYPE)*intervalsToStore);

	out.write((char *) randNNI, sizeof(FLOAT_TYPE)*intervalsToStore);
	out.write((char *) randNNInum, sizeof(int)*intervalsToStore);
	
	out.write((char *) exNNI, sizeof(FLOAT_TYPE)*intervalsToStore);
	out.write((char *) exNNInum, sizeof(int)*intervalsToStore);

	out.write((char *) randSPR, sizeof(FLOAT_TYPE)*intervalsToStore);
	out.write((char *) randSPRnum, sizeof(int)*intervalsToStore);
	
	out.write((char *) limSPR, sizeof(FLOAT_TYPE)*intervalsToStore);
	out.write((char *) limSPRnum, sizeof(int)*intervalsToStore);

	out.write((char *) exlimSPR, sizeof(FLOAT_TYPE)*intervalsToStore);
	out.write((char *) exlimSPRnum, sizeof(int)*intervalsToStore);

	out.write((char *) randRecom, sizeof(FLOAT_TYPE)*intervalsToStore);
	out.write((char *) randRecomnum, sizeof(int)*intervalsToStore);
	
	out.write((char *) bipartRecom, sizeof(FLOAT_TYPE)*intervalsToStore);
	out.write((char *) bipartRecomnum, sizeof(int)*intervalsToStore);

	out.write((char *) onlyBrlen, sizeof(FLOAT_TYPE)*intervalsToStore);
	out.write((char *) onlyBrlennum, sizeof(int)*intervalsToStore);

	out.write((char *) anyModel, sizeof(FLOAT_TYPE)*intervalsToStore);
	out.write((char *) anyModelnum, sizeof(int)*intervalsToStore);
	}
*/

void Adaptation::ReadFromCheckpoint(FILE *in){
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
	fread((char *) improvetotal, sizeof(FLOAT_TYPE), intervalsToStore, in);

	fread((char *) randNNI, sizeof(FLOAT_TYPE), intervalsToStore, in);
	fread((char *) randNNInum, sizeof(int), intervalsToStore, in);
	
	fread((char *) exNNI, sizeof(FLOAT_TYPE), intervalsToStore, in);
	fread((char *) exNNInum, sizeof(int), intervalsToStore, in);

	fread((char *) randSPR, sizeof(FLOAT_TYPE), intervalsToStore, in);
	fread((char *) randSPRnum, sizeof(int), intervalsToStore, in);
	
	fread((char *) limSPR, sizeof(FLOAT_TYPE), intervalsToStore, in);
	fread((char *) limSPRnum, sizeof(int), intervalsToStore, in);

	fread((char *) exlimSPR, sizeof(FLOAT_TYPE), intervalsToStore, in);
	fread((char *) exlimSPRnum, sizeof(int), intervalsToStore, in);

	fread((char *) randRecom, sizeof(FLOAT_TYPE), intervalsToStore, in);
	fread((char *) randRecomnum, sizeof(int), intervalsToStore, in);
	
	fread((char *) bipartRecom, sizeof(FLOAT_TYPE), intervalsToStore, in);
	fread((char *) bipartRecomnum, sizeof(int), intervalsToStore, in);

	fread((char *) onlyBrlen, sizeof(FLOAT_TYPE), intervalsToStore, in);
	fread((char *) onlyBrlennum, sizeof(int), intervalsToStore, in);

	fread((char *) anyModel, sizeof(FLOAT_TYPE), intervalsToStore, in);
	fread((char *) anyModelnum, sizeof(int), intervalsToStore, in);
	}


void Adaptation::PrepareForNextInterval(){
 //if we're on the first generation of a new recording period, shift everything over
	for(int i=intervalsToStore-1; i>0; i--){
		improvetotal[i] = improvetotal[i-1];
		randNNI[i] = randNNI[i-1];randNNInum[i] = randNNInum[i-1];
		exNNI[i] = exNNI[i-1]; exNNInum[i] = exNNInum[i-1];
		 
		randSPR[i] = randSPR[i-1];randSPRnum[i] = randSPRnum[i-1];
		limSPR[i] = limSPR[i-1]; limSPRnum[i] = limSPRnum[i-1];
		exlimSPR[i] = exlimSPR[i-1]; exlimSPRnum[i] = exlimSPRnum[i-1];
#ifdef GANESH
		randPECR[i] = randPECR[i-1]; randPECRnum[i] = randPECRnum[i-1];			
#endif			
//		taxonSwap[i] = taxonSwap[i-1];taxonSwapnum[i] = taxonSwapnum[i-1];
		randRecom[i] = randRecom[i-1];randRecomnum[i] = randRecomnum[i-1];
		bipartRecom[i]=bipartRecom[i-1];bipartRecomnum[i] = bipartRecomnum[i-1];

		onlyBrlen[i] = onlyBrlen[i-1];onlyBrlennum[i] = onlyBrlennum[i-1];
		anyModel[i] = anyModel[i-1]; anyModelnum[i] = anyModelnum[i-1];
//		slopes[i] = slopes[i-1];
#ifdef MPI_VERSION
		bestFromRemote[i] = bestFromRemote[i-1];
		bestFromRemoteNum[i] = bestFromRemoteNum[i-1];
#endif
//		fromRemoteSubtree[i] = fromRemoteSubtree[i-1];
//		fromRemoteNonSubtree[i] = fromRemoteNonSubtree[i-1];
		}// end of for loop

	// clean up for next entry
	improvetotal[0] = 0.0;
	randNNI[0] = 0.0;randNNInum[0] = 0;
	exNNI[0] = 0.0;exNNInum[0] = 0;
		 
	randSPR[0] = 0.0;randSPRnum[0] = 0;
	limSPR[0] = 0.0;limSPRnum[0] = 0;
	exlimSPR[0] = 0.0;exlimSPRnum[0] = 0;
#ifdef GANESH
	randPECR[0] = 0.0;randPECRnum[0] = 0; 
#endif
//	taxonSwap[0] = 0.0;taxonSwapnum[0] = 0;
	randRecom[0] = 0.0;randRecomnum[0] = 0;
	bipartRecom[0]=0.0;bipartRecomnum[0]=0;

	onlyBrlen[0] = 0.0;onlyBrlennum[0] = 0;
	anyModel[0] = 0.0;anyModelnum[0] = 0;
#ifdef MPI_VERSION
	bestFromRemote[0]=0.0;bestFromRemoteNum[0]=0;
#endif
//	fromRemoteSubtree[0] = 0.0;
//	fromRemoteNonSubtree[0] = 0.0;
	}

void Adaptation::BeginProbLog(ofstream &plog, int gen){
	plog << "gen\tmod\ttopo\tbrlen\tNNI\trandSPR\tlimSPR\n";
	//gen could be non-zero if we restarted
	OutputProbs(plog, gen);
	}

void Adaptation::OutputProbs(ofstream &plog, int gen){
	plog << gen << "\t" << modelMutateProb << "\t" << topoMutateProb << "\t" << (1.0-modelMutateProb-topoMutateProb) << "\t";
	plog << randNNIprob << "\t" << randSPRprob << "\t" << limSPRprob << endl;
	}

void Adaptation::UpdateProbs(){
	FLOAT_TYPE topoTot=0.0, modTot=0.0, onlyBrlenTot=0.0;
	int numTopos=0, numMod=0, numOnlyBrlen=0;
	
	FLOAT_TYPE totRandNNI=0.0, totLimSPR=0.0, totRandSPR=0.0;
	int totNumRandNNI=0, totNumLimSPR=0, totNumRandSPR=0;
	FLOAT_TYPE totBipartRecom=0.0;
	int totNumBipartRecom=0;

#ifdef MPI_VERSION
	FLOAT_TYPE totalFromRemote;
#endif


#ifdef GANESH
	FLOAT_TYPE totRandPECR=0.0;
	int totNumRandPECR=0;
#endif
	
	for(unsigned i=0;i<intervalsToStore;i++){
#ifdef GANESH
		topoTot += randNNI[i] + randSPR[i] + limSPR[i] + randPECR[i];
		numTopos += randNNInum[i] + randSPRnum[i] + limSPRnum[i] +
		randPECRnum[i];
#else
		topoTot += randNNI[i] + randSPR[i] + limSPR[i];
		numTopos += randNNInum[i] + randSPRnum[i] + limSPRnum[i];
#endif
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
		
#ifdef GANESH
		totRandPECR += randPECR[i];
		totNumRandPECR += randPECRnum[i];
#endif
		totNumRandSPR += randSPRnum[i];	

#ifdef	MPI_VERSION
	totalFromRemote += bestFromRemote[i];
#endif
		}
	
	FLOAT_TYPE perTopo, perModel, perBrlen;
	if(numTopos!=0) perTopo=(topoTot/numTopos);
	else perTopo=0.0;
	if(numMod!=0) perModel=(modTot/numMod);
	else perModel=0.0;
	if(numOnlyBrlen>0) perBrlen=(onlyBrlenTot/numOnlyBrlen);
	else perBrlen= 0.0;
	FLOAT_TYPE perBipartRecom;
	if(totNumBipartRecom > 0) perBipartRecom=(totBipartRecom/totNumBipartRecom);
	else perBipartRecom=0.0;
	
	//version 0.95b3 - The reduction of precision that used to appear here has been
	//moved to Adaptation::ReducePrecision, which is called from Run, MasterMaster and
	//RemoteSubtreeWorker when lastTopoImprove is > that #int * intLength generations ago
	
	perTopo += topoWeight;
	perModel += modWeight;
	perBrlen += brlenWeight;
	
	FLOAT_TYPE tot=perTopo+perModel+perBrlen;

	FLOAT_TYPE brlenOnlyMut;

	//only update these probs if model mutations are turned off completely
	//or if some model mutations have been done (ie not in subtree mode)
	if(anyModelnum[0]!=0 || FloatingPointEquals(modWeight, 0.0, 1e-10)){
		brlenOnlyMut=perBrlen/tot;
		modelMutateProb = perModel/tot;
		topoMutateProb = perTopo/tot;
		}

	//enforce a minimum probability
	if(modWeight != 0.0 && topoWeight != 0.0){
		FLOAT_TYPE minProb= (FLOAT_TYPE) 0.02;	
		if(topoMutateProb < minProb){
			modelMutateProb -= minProb - topoMutateProb;
			topoMutateProb=minProb;
			}
		if(modelMutateProb < minProb){
			topoMutateProb -= minProb - modelMutateProb;
			modelMutateProb=minProb;
			}
		if(1.0 - (modelMutateProb + topoMutateProb) < minProb){
			FLOAT_TYPE diff=minProb - (FLOAT_TYPE)(1.0 - (modelMutateProb + topoMutateProb));
			if(modelMutateProb - diff/2.0 > .02 && topoMutateProb - diff/2.0 > .02){
				modelMutateProb -= diff/(FLOAT_TYPE)2.0;
				topoMutateProb -= diff/(FLOAT_TYPE)2.0;
				}
			else{
				if(modelMutateProb - diff/2.0 < .02){
					topoMutateProb -= diff;
					}
				else modelMutateProb -= diff;
				}		

			brlenOnlyMut=minProb;
			}
		}
		
//		brlenOnlyMut = 1.0 - topoMutateProb - modelMutateProb;
//		}
/*	else{
		scaler=(1-brlenOnlyMut) / (topoMutateProb + modelMutateProb);
		modelMutateProb *= scaler;
		topoMutateProb *= scaler;
		}
*/
	if(totNumRandNNI==0) totNumRandNNI=1;
	if(totNumLimSPR==0) totNumLimSPR=1;
	if(totNumRandSPR==0) totNumRandSPR=1;
#ifdef GANESH
	if(totNumRandPECR==0) totNumRandPECR=1;
#endif
	//Because NNI's chosen by an SPR mutator are marked as NNI's, this needs to be done to keep from 
	//giving NNI's some prob even when the weight was 0.0
	FLOAT_TYPE perRandNNI= (randNNIweight == ZERO_POINT_ZERO ? ZERO_POINT_ZERO : totRandNNI/totNumRandNNI + randNNIweight);
	FLOAT_TYPE perLimSPR=  (limSPRweight  == ZERO_POINT_ZERO ? ZERO_POINT_ZERO : totLimSPR/totNumLimSPR + limSPRweight);
	FLOAT_TYPE perRandSPR= (limSPRweight  == ZERO_POINT_ZERO ? ZERO_POINT_ZERO : totRandSPR/totNumRandSPR + randSPRweight);
	
#ifdef GANESH
	FLOAT_TYPE perRandPECR=totRandPECR/totNumRandPECR + randPECRweight;
	tot=perRandNNI+perLimSPR+perRandSPR+perRandPECR;
	randNNIprob=perRandNNI/tot;
	randSPRprob=perRandSPR/tot;
	limSPRprob=perLimSPR/tot;
	randPECRprob=perRandPECR/tot;
#else

	tot=perRandNNI+perLimSPR+perRandSPR;
	if(tot > ZERO_POINT_ZERO){
		randNNIprob=perRandNNI/tot;
		randSPRprob=perRandSPR/tot;
		limSPRprob=perLimSPR/tot;
		}
#endif
	}
