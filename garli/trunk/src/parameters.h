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

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "configoptions.h"
#include "stricl.h"

class rng;
class DNAData;
class HKYData;

class Parameters{
	public:
	HKYData* data;

	// parameters that CANNOT be set in the configuration file
	int restart;
	int stopnow;
	int showProgress;
	char logfname[85];
	char startfname[85];
	string constraintFile;
	char treefname[85];
	char gnufname[85];
	NxsString starting_tree;
	double prev_time;
	long prev_generations;

	// parameters that CAN be set in the configuration file
	int nindivs;//, min_nindivs, max_nindivs;
	int stopgen;
	int stoptime;
	int holdover;
	int logEvery;
	int fatlog;
	int saveEvery;
	long randomSeed;
	int myrank;
	int nprocesses;
	char datafname[81];
	char ofprefix[81];
	char plottitle[81];
	char statefname[81];
	double gammaShapeBrlen;
	double gammaShapeModel;
	double selectionIntensity;
	double holdoverPenalty;
	double meanBrlenMuts;//, min_brlen_muts, max_brlen_muts;
	double treeRejectionThreshold;

	Parameters() : data(0), restart(0), stopnow(0), showProgress(0)
		{ FactorySettings(); }
	~Parameters() {}

	void FactorySettings();
	static void ShowCorrespondence( ostream& out );
	static void ShowDefinitions( ostream& out );

	void ReadTestIndivs( istream& inf );
	int CheckValidity( char* msg, int msglen );
	void BriefReport( ostream& out );
	void ReadParams(istream&);
	int SetParams(const GeneralGamlConfig&, const HKYData&);

	int Serialize(char*& buf, int* size_);  // cjb
	int Deserialize(const char* buf);  // cjb
	bool operator==(const Parameters& rhs) const; // cjb - to test serialization code
	//int Randomize(); // cjb - semi randomly determine running parameters
	Parameters& operator=(const Parameters& rhs);  // cjb


};
// istream& operator>>( istream&, Parameters& p );
 ostream& operator<<( ostream&, Parameters& p );
	
#endif

