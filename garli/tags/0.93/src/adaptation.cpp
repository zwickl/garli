

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

#include <iostream>
#include <fstream>
#include <cstdio>

using namespace std;

#include "adaptation.h"
#include "math.h"
#include "configoptions.h"
#include "individual.h"

Adaptation::Adaptation(const GeneralGamlConfig *gc){

	intervalsToStore=gc->intervalsToStore;
	intervalLength=gc->intervalLength;

	startOptPrecision = branchOptPrecision = gc->startOptPrec;
	minOptPrecision = gc->minOptPrec;
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
	taxonSwapprob=0.0;

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
	for(int i=0;i<intervalsToStore;i++){
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

void Adaptation::BeginProbLog(ofstream &plog){
	plog << "gen\tmod\ttopo\tbrlen\tNNI\trandSPR\tlimSPR\n";
	OutputProbs(plog, 0);
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
	
	for(int i=0;i<intervalsToStore;i++){
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
	
	//if there have been no beneficial topo mutations recently, decrease the precision
	//9/12/05 changing this to require no _significant_ topo improvements, rather than any
	//at all.  ie, istead of perTopo==0.0, perTopo < 0.0001
#ifndef MPI_VERSION
	if(perTopo<0.0001 && topoMutateProb > .1) branchOptPrecision *= precReductionFactor;
#else
	if(perTopo<0.0001 && topoMutateProb > .1 && perBipartRecom<0.0001) branchOptPrecision *= precReductionFactor;
#endif
	if(branchOptPrecision < minOptPrecision) branchOptPrecision=minOptPrecision;
	
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
