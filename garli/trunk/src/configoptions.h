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


#ifndef CONFIGOPTIONS_H
#define CONFIGOPTIONS_H

#include <string>

using std::string;

#include "hashdefines.h"

class GeneralGamlConfig{
	public:
	//these options will be the same regardless of whether a population is master or remote
	int logevery;
	int saveevery;
	double megsClaMemory;
	double availableMemory;
	int randseed;

	string datafname;
	string method;
	string ofprefix;
	string streefname;
	string constraintfile;
	bool refineStart;
	bool outputTreelog;
	bool outputMostlyUselessFiles;
	bool outputPhylipTree;
	bool dontInferProportionInvariant;
	bool useflexrates;
	int numratecats;

	bool enforceTermConditions;
	int lastTopoImproveThresh;
	double improveOverStoredIntervalsThresh;
	double significantTopoChange;

	//all of the following options can vary between master and remote
	//general population stuff
	int max_nindivs, min_nindivs;
	int holdover;
	double selectionIntensity;
	double holdoverPenalty;
	int stopgen;
	int stoptime;

	double startOptPrec;
	double minOptPrec;
	double precReductionFactor;
	int numPrecReductions;

	double treeRejectionThreshold;

	//parameters affecting proportion of mutations
	double topoWeight;
		double randNNIweight;
		double randSPRweight;
		double limSPRweight;
        double randPECRweight;
	double modWeight;
	double brlenWeight;

	int intervalLength;
	int intervalsToStore;

	int bootstrapReps;
	bool inferInternalStateProbs;

	//parameters affecting other details of mutations				
	double maxBrlenMuts, minBrlenMuts;
	int gammaShapeBrlen;
	int gammaShapeModel;
	int limSPRrange;		

	//perturbation parameters
	int pertType;			
	double pertThresh;
	int minPertInterval;
	int maxPertsNoImprove;
	bool restartAfterAbandon;
	int gensBeforeRestart;
	
	double ratchetProportion;
	double ratchetOffThresh;
	int ratchetMaxGen;
	
	int nniTargetAccepts;
	int nniMaxAttempts;
	
	int numSprCycles;
	int sprPertRange;

	//the number of seconds between remote tree sends
	double sendInterval;
	
 	// methods
	GeneralGamlConfig();
	int Read(const char*, bool isMaster=false);
	int Serialize(char**, int*) const;
	int Deserialize(char*, int);
	bool operator==(const GeneralGamlConfig&) const;
	};
		
class MasterGamlConfig: public GeneralGamlConfig{
	public:
	//parallel behavior parameters-stored in pop->paraMan on master only
	double startUpdateThresh;
	double minUpdateThresh;
	double updateReductionFactor;
		
	int subtreeInterval;
	double subtreeStartThresh;
	int minSubtreeSize;
	int targetSubtreeSize;
	double orphanFactor;
	
	int maxRecomIndivs;
/*	
	int pertType;			
	double pertThresh;
	double pertAmount;
	int maxPertsNoImprove;
	
	double ratchetProportion;
	double ratchetOffThresh;
	int ratchetMaxGen;
	
	double nniAcceptThresh;
	int numSprCycles;
	int sprPertRange;
*/
	int bootstrapReps;
	double bootTermThresh;

 	// methods
	MasterGamlConfig();
	int Read(const char*, bool isMaster=false);
	};

#endif

