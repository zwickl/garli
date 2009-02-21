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


#ifndef CONFIGOPTIONS_H
#define CONFIGOPTIONS_H

#include <string.h>
#include <string>

using std::string;
#include "hashdefines.h"
#include "defs.h"

enum GarliRunMode {
	NORMAL_RUN_MODE = 0,
	NSWAPS_RUN_MODE = 1,
	SWAPPER_BY_DIST_NOT_FURTHEST_RUN_MODE = 2,
	SWAPPER_BY_CUT_NOT_FURTHEST_RUN_MODE = 3,
	SWAPPER_RANDOM_RUN_MODE = 4,
	SWAPPER_BY_DIST_FURTHEST_RUN_MODE = 5,
	SWAPPER_BY_CUT_FURTHEST_RUN_MODE = 6,
	NONREDUCING_VAR_START_RUN_MODE = 7,
	SITE_RATES_RUN_MODE = 8,
	REDUCING_VAR_START_RUN_MODE = 9,
	ADD_TAXON_RUN_MODE = 10,
	GENERATE_TREES_RUNMODE = 21
	};
typedef std::pair<std::string, std::string> KeyValueStringPair;
typedef std::pair<bool, KeyValueStringPair> AttemptedParseResult;
class GeneralGamlConfig {
	public:
	//these options will be the same regardless of whether a population is master or remote

	//output related
	string ofprefix;
	unsigned logevery;
	unsigned saveevery;
	bool outputTreelog;
	bool outputMostlyUselessFiles;
	bool outputPhylipTree;
	bool outputCurrentBestTopology;

	bool collapseBranches;
	//this is just a string that I can use for whatever I want in special runmodes
	string arbitraryString;

	//starting the run
	int randseed;

	public:
		enum StartingTree {
			RANDOM_START,
			STEPWISE_ADDITION_START,
			FROM_FILE_START,
			INCOMPLETE_FROM_FILE_START
		};
		StartingTree GetStartMode() const;
		const char * GetTreeFilename() const;
		bool ReadTreeFromFile() const {
			StartingTree m = GetStartMode();
			return (m == FROM_FILE_START || m == INCOMPLETE_FROM_FILE_START);
		}
		bool		  IsValid() const;
		
		AttemptedParseResult ParseLineIntoConfigObject(const std::string line);		
		bool  		* GetBoolOptReference(const char * n) ;
		double 		* GetDoubleOptReference(const char * n) ;
		int			* GetIntOptReference(const char * n);
		std::string * GetStringOptReference(const char * n);
		unsigned 	* GetUnsignedOptReference(const char * n);		
	private:
		string streefname;
		string incompletetreefname;
	public:
	bool refineStart;

	//general run details
	string datafname;
	string constraintfile;
	FLOAT_TYPE megsClaMemory;
	FLOAT_TYPE availableMemory;
	bool restart;
	bool checkpoint;
	FLOAT_TYPE significantTopoChange;
	string outgroupString;
	unsigned searchReps;
	public:
		GarliRunMode GetRunMode() const {return GarliRunMode(runmode);}
	private:
		unsigned runmode;
	public:
	//finishing the run
	bool enforceTermConditions;
	unsigned lastTopoImproveThresh;
	FLOAT_TYPE improveOverStoredIntervalsThresh;
	unsigned stopgen;
	unsigned stoptime;

	unsigned attachmentsPerTaxon;

	//model settings
	string datatype;
	string geneticCode;
	string stateFrequencies; //equal, estimate, emprical, fixed
	string rateMatrix;		 //6rate, 2rate, 1rate, fixed, custom(
	string proportionInvariant; //none, fixed, estimate
	string rateHetModel;			//gamma, gammafixed, flex, none
	unsigned numRateCats;	

	//all of the following options can vary between master and remote
	//general population stuff
	unsigned nindivs;
	unsigned holdover;
	FLOAT_TYPE selectionIntensity;
	FLOAT_TYPE holdoverPenalty;

	FLOAT_TYPE startOptPrec;
	FLOAT_TYPE minOptPrec;
	int numPrecReductions;
	FLOAT_TYPE precReductionFactor; //deprecated
	FLOAT_TYPE treeRejectionThreshold;

	//parameters affecting proportion of mutations
	FLOAT_TYPE topoWeight;
		FLOAT_TYPE randNNIweight;
		FLOAT_TYPE randSPRweight;
		FLOAT_TYPE limSPRweight;
//      FLOAT_TYPE randPECRweight;
	FLOAT_TYPE modWeight;
	FLOAT_TYPE brlenWeight;

	unsigned intervalLength;
	unsigned intervalsToStore;

	//parameters affecting other details of mutations				
	FLOAT_TYPE meanBrlenMuts;
	FLOAT_TYPE gammaShapeBrlen;
	FLOAT_TYPE gammaShapeModel;
	unsigned limSPRrange;		
	FLOAT_TYPE uniqueSwapBias;
	FLOAT_TYPE distanceSwapBias;
	
	//optional analyses
	unsigned bootstrapReps;
	FLOAT_TYPE resampleProportion;
	bool inferInternalStateProbs;

#ifdef INCLUDE_PERTURBATION
	//perturbation parameters
	int pertType;			
	FLOAT_TYPE pertThresh;
	int minPertInterval;
	int maxPertsNoImprove;
	bool restartAfterAbandon;
	int gensBeforeRestart;
	
	FLOAT_TYPE ratchetProportion;
	FLOAT_TYPE ratchetOffThresh;
	int ratchetMaxGen;
	
	int nniTargetAccepts;
	int nniMaxAttempts;
	
	int numSprCycles;
	int sprPertRange;
#endif

	//the number of seconds between remote tree sends (parallel only)
	FLOAT_TYPE sendInterval;
	
	//by default these come from the defs.h file, but could be overriden
	FLOAT_TYPE minBrlen;
	FLOAT_TYPE maxBrlen;
	FLOAT_TYPE startingBrlen;

 	// methods
	GeneralGamlConfig();
	GeneralGamlConfig(const GeneralGamlConfig &);
	int Read(const char*, bool isMaster=false);
	int Serialize(char**, int*) const;
	int Deserialize(char*, int);
	bool operator==(const GeneralGamlConfig&) const;
	};
		
class MasterGamlConfig : public GeneralGamlConfig{
	public:
	//parallel behavior parameters-stored in pop->paraMan on master only
	FLOAT_TYPE startUpdateThresh;
	FLOAT_TYPE minUpdateThresh;
	FLOAT_TYPE updateReductionFactor;
		
	int subtreeInterval;
	FLOAT_TYPE subtreeStartThresh;
	int minSubtreeSize;
	int targetSubtreeSize;
	FLOAT_TYPE orphanFactor;
	
	int maxRecomIndivs;
/*	
	int pertType;			
	FLOAT_TYPE pertThresh;
	FLOAT_TYPE pertAmount;
	int maxPertsNoImprove;
	
	FLOAT_TYPE ratchetProportion;
	FLOAT_TYPE ratchetOffThresh;
	int ratchetMaxGen;
	
	FLOAT_TYPE nniAcceptThresh;
	int numSprCycles;
	int sprPertRange;
*/
	int bootstrapReps;
	FLOAT_TYPE bootTermThresh;

 	// methods
	MasterGamlConfig();
	MasterGamlConfig(const MasterGamlConfig & other)
		:GeneralGamlConfig(other) {
		*this = other;
	}
	int Read(const char*, bool isMaster=false);
	};

#endif

