

// GARLI version 0.951 source code
// Copyright  2005-2006 by Derrick J. Zwickl
// All rights reserved.
//
// This code may be used and modified for non-commercial purposes
// but redistribution in any form requires written permission.
// Please contact:
//
//  Derrick Zwickl
//	National Evolutionary Synthesis Center
//	2024 W. Main Street, Suite A200
//	Durham, NC 27705
//  email: zwickl@nescent.org
//


#include <iostream>
#include <fstream>
#include <cstdio>

using namespace std;

#include "defs.h"
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
		precReductionFactor = (startOptPrecision - minOptPrecision)/double(numPrecReductions);
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

	double tot = topoWeight+modWeight+brlenWeight;
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

	randNNI = new double[intervalsToStore]; randNNInum = new int[intervalsToStore];
	exNNI = new double[intervalsToStore]; exNNInum = new int[intervalsToStore];
	randSPR = new double[intervalsToStore]; randSPRnum = new int[intervalsToStore];
#ifdef GANESH
	randPECR = new double[intervalsToStore]; randPECRnum = new int[intervalsToStore];
#endif
	limSPR = new double[intervalsToStore]; limSPRnum = new int[intervalsToStore];
	exlimSPR = new double[intervalsToStore]; exlimSPRnum = new int[intervalsToStore];
	randRecom = new double[intervalsToStore];	 randRecomnum = new int[intervalsToStore];
	bipartRecom = new double[intervalsToStore];	 bipartRecomnum = new int[intervalsToStore];
	onlyBrlen = new double[intervalsToStore]; onlyBrlennum = new int[intervalsToStore];
	improvetotal = new double[intervalsToStore];
	anyModel = new double[intervalsToStore]; anyModelnum = new int[intervalsToStore];

#ifdef MPI_VERSION
	bestFromRemote=new double[intervalsToStore];bestFromRemoteNum=new int[intervalsToStore];
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
		precReductionFactor = (gc->startOptPrec- minOptPrecision)/double(numPrecReductions);
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

void Adaptation::WriteToCheckpoint(ofstream &out){
	//this function assumes that it has been passed a stream that is already open for 
	//binary writing
	assert(out.good());

	//first take care of the scalars, which all come first in the class
	int scalarSize = sizeof(int)*4 + sizeof(bool) + sizeof(double)*22;
	out.write((char *) this, scalarSize);

	//now the arrays, which should be of length intervalsToStore
	out.write((char *) improvetotal, sizeof(double)*intervalsToStore);

	out.write((char *) randNNI, sizeof(double)*intervalsToStore);
	out.write((char *) randNNInum, sizeof(int)*intervalsToStore);
	
	out.write((char *) exNNI, sizeof(double)*intervalsToStore);
	out.write((char *) exNNInum, sizeof(int)*intervalsToStore);

	out.write((char *) randSPR, sizeof(double)*intervalsToStore);
	out.write((char *) randSPRnum, sizeof(int)*intervalsToStore);
	
	out.write((char *) limSPR, sizeof(double)*intervalsToStore);
	out.write((char *) limSPRnum, sizeof(int)*intervalsToStore);

	out.write((char *) exlimSPR, sizeof(double)*intervalsToStore);
	out.write((char *) exlimSPRnum, sizeof(int)*intervalsToStore);

	out.write((char *) randRecom, sizeof(double)*intervalsToStore);
	out.write((char *) randRecomnum, sizeof(int)*intervalsToStore);
	
	out.write((char *) bipartRecom, sizeof(double)*intervalsToStore);
	out.write((char *) bipartRecomnum, sizeof(int)*intervalsToStore);

	out.write((char *) onlyBrlen, sizeof(double)*intervalsToStore);
	out.write((char *) onlyBrlennum, sizeof(int)*intervalsToStore);

	out.write((char *) anyModel, sizeof(double)*intervalsToStore);
	out.write((char *) anyModelnum, sizeof(int)*intervalsToStore);
	}

void Adaptation::ReadFromCheckpoint(ifstream &in){
	//this function assumes that it has been passed a stream that is already open for 
	//binary reading
	assert(in.good());
	int scalarSize = sizeof(int)*4 + sizeof(bool) + sizeof(double)*22;

	in.read((char *) this, scalarSize);

	//now the arrays, which should be of length intervalsToStore
	in.read((char *) improvetotal, sizeof(double)*intervalsToStore);

	in.read((char *) randNNI, sizeof(double)*intervalsToStore);
	in.read((char *) randNNInum, sizeof(int)*intervalsToStore);
	
	in.read((char *) exNNI, sizeof(double)*intervalsToStore);
	in.read((char *) exNNInum, sizeof(int)*intervalsToStore);

	in.read((char *) randSPR, sizeof(double)*intervalsToStore);
	in.read((char *) randSPRnum, sizeof(int)*intervalsToStore);
	
	in.read((char *) limSPR, sizeof(double)*intervalsToStore);
	in.read((char *) limSPRnum, sizeof(int)*intervalsToStore);

	in.read((char *) exlimSPR, sizeof(double)*intervalsToStore);
	in.read((char *) exlimSPRnum, sizeof(int)*intervalsToStore);

	in.read((char *) randRecom, sizeof(double)*intervalsToStore);
	in.read((char *) randRecomnum, sizeof(int)*intervalsToStore);
	
	in.read((char *) bipartRecom, sizeof(double)*intervalsToStore);
	in.read((char *) bipartRecomnum, sizeof(int)*intervalsToStore);

	in.read((char *) onlyBrlen, sizeof(double)*intervalsToStore);
	in.read((char *) onlyBrlennum, sizeof(int)*intervalsToStore);

	in.read((char *) anyModel, sizeof(double)*intervalsToStore);
	in.read((char *) anyModelnum, sizeof(int)*intervalsToStore);
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
	double topoTot=0.0, modTot=0.0, onlyBrlenTot=0.0;
	int numTopos=0, numMod=0, numOnlyBrlen=0;
	
	double totRandNNI=0.0, totLimSPR=0.0, totRandSPR=0.0;
	int totNumRandNNI=0, totNumLimSPR=0, totNumRandSPR=0;
	double totBipartRecom=0.0;
	int totNumBipartRecom=0;

#ifdef MPI_VERSION
	double totalFromRemote;
#endif


#ifdef GANESH
	double totRandPECR=0.0;
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
	
	double perTopo, perModel, perBrlen;
	if(numTopos!=0) perTopo=(topoTot/numTopos);
	else perTopo=0.0;
	if(numMod!=0) perModel=(modTot/numMod);
	else perModel=0.0;
	if(numOnlyBrlen>0) perBrlen=(onlyBrlenTot/numOnlyBrlen);
	else perBrlen= 0.0;
	double perBipartRecom;
	if(totNumBipartRecom > 0) perBipartRecom=(totBipartRecom/totNumBipartRecom);
	else perBipartRecom=0.0;
	
	//version 0.95b3 - The reduction of precision that used to appear here has been
	//moved to Adaptation::ReducePrecision, which is called from Run, MasterMaster and
	//RemoteSubtreeWorker when lastTopoImprove is > that #int * intLength generations ago
	
	perTopo += topoWeight;
	perModel += modWeight;
	perBrlen += brlenWeight;
	
	double tot=perTopo+perModel+perBrlen;

	double brlenOnlyMut;

	//only update these probs if model mutations are turned off completely
	//or if some model mutations have been done (ie not in subtree mode)
	if(anyModelnum[0]!=0 || modWeight == 0.0){
		brlenOnlyMut=perBrlen/tot;
		modelMutateProb = perModel/tot;
		topoMutateProb = perTopo/tot;
		}

	//enforce a minimum probability
	if(modWeight != 0.0 && topoWeight != 0.0){
		double minProb=.02;	
		if(topoMutateProb < minProb){
			modelMutateProb -= minProb - topoMutateProb;
			topoMutateProb=minProb;
			}
		if(modelMutateProb < minProb){
			topoMutateProb -= minProb - modelMutateProb;
			modelMutateProb=minProb;
			}
		if(1.0 - (modelMutateProb + topoMutateProb) < minProb){
			double diff=minProb - (1.0 - (modelMutateProb + topoMutateProb));
			if(modelMutateProb - diff/2.0 > .02 && topoMutateProb - diff/2.0 > .02){
				modelMutateProb -= diff/2.0;
				topoMutateProb -= diff/2.0;
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
		
	double perRandNNI=totRandNNI/totNumRandNNI + randNNIweight;
	double perLimSPR=totLimSPR/totNumLimSPR + limSPRweight;
	double perRandSPR=totRandSPR/totNumRandSPR + randSPRweight;
#ifdef GANESH
	double perRandPECR=totRandPECR/totNumRandPECR + randPECRweight;
	tot=perRandNNI+perLimSPR+perRandSPR+perRandPECR;
	randNNIprob=perRandNNI/tot;
	randSPRprob=perRandSPR/tot;
	limSPRprob=perLimSPR/tot;
	randPECRprob=perRandPECR/tot;
#else
	tot=perRandNNI+perLimSPR+perRandSPR;
	randNNIprob=perRandNNI/tot;
	randSPRprob=perRandSPR/tot;
	limSPRprob=perLimSPR/tot;
#endif
	}
