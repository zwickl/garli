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

//	NOTE: Portions of this source adapted from GAML source, written by Paul O. Lewis

#include <iomanip>

#include "parameters.h"
#include "datamatr.h"
#include "rng.h"
#include "funcs.h"

#define MAX_NAME_LENGTH  80

//
//
// Methods for class Parameters
//
//

void Parameters::FactorySettings()
{
/*	// parameters that CAN be set in the configuration file
	strcpy( ofprefix, "sample" );
	strcpy( plottitle, "Title" );
	nindivs			= DEF_NINDIVS;
	holdover		= DEF_HOLDOVER;
	logEvery		= DEF_LOGEVERY;
	fatlog		= 0;
	saveEvery		= DEF_SAVEEVERY;
	randomSeed		= 0L;
	rnd.randomize();
	gammaShape		= DEF_BRLEN_SHAPE;
	mean_brlen_muts	= DEF_BRLEN_MUT_PROB;

	// parameters that CANNOT be set in the configuration file
	myrank			= 0;
	nprocesses		= 1;
	prev_time		= 0.0;
	prev_generations	= 0L;
	restart			= 0;
	stopnow			= 0;
	showProgress		= 0;
	starting_tree           = "";
	strcpy( datafname, "gaml.dat" );
	strcpy( logfname, "sample.log" );
	strcpy( treefname, "sample.tre" );
	strcpy( gnufname, "sample.gnu" );
	strcpy( statefname, "sample.gml" );

	strcpy( startfname, "random" );
*/
}

// Check params whose validity does not depend on the data having been read
// returns 1 if parameters are valid, 0 if not valid
int Parameters::CheckValidity( char* msg, int msglen )
{
	if( nprocesses > nindivs+1 ) {
		strncpy( msg, "Sorry, population size must be at least one less than number of processors.", msglen );
		return 0;
	}
/*	if( topoMutProb < 0.0 || topoMutProb > 1.0 ) {
		strncpy( msg, "Sorry, mutation probabilities must be between 0.0 and 1.0.", msglen );
		return 0;
	}

	if( brlenMutProb < 0.0 || brlenMutProb > 1.0 ) {
		strncpy( msg, "Sorry, mutation probabilities must be between 0.0 and 1.0.", msglen );
		return 0;
	}
	if( crossoverProb < 0.0 || crossoverProb > 1.0 ) {
		strncpy( msg, "Sorry, crossover probabilities must be between 0.0 and 1.0.", msglen );
		return 0;
	}
	if( gammaShapeBrlen <= 0.0 ) {
		strncpy( msg, "Sorry, gamma shape must be greater than 0.0.", msglen );
		return 0;
	}
*/
	return 1;
}

int Parameters::SetParams(const GeneralGamlConfig& conf, const HKYData& data_in)	{
	// general config
	strncpy(datafname, conf.datafname.c_str(), 80);
	strncpy(ofprefix, conf.ofprefix.c_str(), 80);
	if (conf.streefname.length() == 0)
		strcpy(startfname, "random");
	else
		strncpy(startfname, conf.streefname.c_str(), 84);
	logEvery = conf.logevery;
	saveEvery = conf.saveevery;

	nindivs = RandomInt(conf.min_nindivs, conf.max_nindivs);
	stopgen = conf.stopgen;
	stoptime = conf.stoptime;
	meanBrlenMuts = RandomDouble(conf.minBrlenMuts, conf.maxBrlenMuts);
	holdover = conf.holdover;
	gammaShapeBrlen = conf.gammaShapeBrlen;
	gammaShapeModel = conf.gammaShapeBrlen;
	selectionIntensity=conf.selectionIntensity;
	holdoverPenalty=conf.holdoverPenalty;
	treeRejectionThreshold = conf.treeRejectionThreshold;

	// construct other settings from the previous ones
	sprintf(logfname, "%s.log", ofprefix);
	sprintf(gnufname, "%s.gnu", ofprefix);
	sprintf(treefname, "%s.best.tre", ofprefix);
	sprintf(statefname, "%s.gml", ofprefix);

	// now the random seedSetParams
	if(conf.randseed < 1)
		randomSeed = RandomInt(1, 1000);
	else randomSeed=conf.randseed;
	rnd.set_seed(randomSeed);

	// finally set the data pointer
	data = const_cast<HKYData*>(&data_in);
	
	return 0;
}

void Parameters::BriefReport( ostream& out )
{
#if defined( ALLOW_USERTREES )
	if( strcmp( testfname, "none" ) != 0 ) {
		out << "  Optimizing parameters of trees in file: " << testfname << endl;
		return;
	}
#endif

	if( restart )
		out << "  Program has been restarted using population recorded in state file " << endl;

	if( randomSeed == 0L ) {
		out << "  Random number seed not specified in the configuration file" << endl;
		out << "    Seed chosen using system clock: " << rnd.init_seed() << endl;
	}
	else {
		out << "  Random number seed specified was " << rnd.init_seed() << endl;
	}

	out << "  Number of individuals: " << nindivs << endl;
	if( strcmp( startfname, "random" ) == 0 )
		out << "  Random trees used to create initial population" << endl;
	else
		out << "  Initial population created by mutating tree in file: " << startfname << endl;
	out << "  Number of generations between log records " << logEvery << endl;
	out << "  Save state of the GA every " << saveEvery << " generations" << endl;
	out << "  Gamma shape used for modifying branch lengths: " << gammaShapeBrlen << endl;

	if( showProgress )
		out << "  Verbose output requested." << endl;
	else
		out << "  Terse output requested." << endl;

	out << "  File names created using the prefix " << ofprefix << endl;
	out << "    data   : " << datafname << endl;
	out << "    log    : " << logfname << endl;
	out << "    gnuplot: " << gnufname << endl;
	out << "    tree   : " << treefname << endl;
	out << "    state  : " << statefname << endl;

	out << "  Title used in GnuPlot command file: " << plottitle << endl;
}

void Parameters::ShowDefinitions( ostream& out )
{
	out << "  datafname     = file name of the data file containing sequences" << endl;
	out << "  ofprefix      = most output files will begin with this prefix" << endl;
	out << "  plottitle     = title to be specified in GnuPlot command file" << endl;
	out << "  randomseed    = seed for pseudorandom number generator" << endl;
	out << "                  if 0, seed will be obtained from system clock" << endl;
	out << "  startwith     = file containing tree definition to use in creating" << endl;
	out << "                  starting population (leave out to get random trees)" << endl;
	out << "  logevery      = number of generations between log file entries" << endl;
	out << "  fatlog        = specifies that extra info be saved in log file entries" << endl;
	out << "                  if 0 (the default) only gen, best score, and time saved" << endl;
	out << "                  specify 1 for extra information in each log file entry" << endl;
	out << "  saveevery     = GAML saves the current population in a 'state' file" << endl;
	out << "                  this specifies the number of generations between such saves" << endl;
	out << "  nindivs       = number of individuals in the population" << endl;
	out << "  holdover      = number of unmodified copies of best individual" << endl;
	out << "                  made each generation" << endl;
	out << "  crossoverprob = crossover probability" << endl;
	out << "  topomutprob   = topological mutation probability" << endl;
	out << "  brlenmutprob  = branch length mutation probability" << endl;
	out << "  gammashape    = shape used for gamma distribution used in" << endl;
	out << "                  modifying branch lengths and other continuous parameters" << endl;
	out << "  kappaprob     = kappa mutation probability" << endl;

#if defined( ALLOW_USERTREES )
	out << "  usertrees     = if specified, GAML will evaluate user trees in this file" << endl;
	out << "                  and then terminate immediately" << endl;
#endif
}

void Parameters::ShowCorrespondence( ostream& out )
{
	out << "  nindivs       = n" << endl;
	out << "  holdover      = k" << endl;
	out << "  topomutprob   = mu" << endl;
	out << "  brlenmutprob  = lambda" << endl;
	out << "  gammashape    = alpha" << endl;
	out << "  kappaprob     = pi" << endl;
	out << "  crossoverprob = r" << endl;
}

ostream& operator<<( ostream& outf, Parameters& p )
{
	outf << setw(20) << "datafname " << setw(12) << p.datafname << endl;
#if defined( ALLOW_USERTREES )
	outf << setw(20) << "usertrees " << setw(12) << p.testfname << endl;
#endif
	if( strcmp( p.startfname, "random" ) != 0 )
		outf << setw(20) << "startwith " << setw(12) << p.startfname << endl;
	outf << setw(20) << "ofprefix " << setw(12) << p.ofprefix << endl;
	outf << setw(20) << "plottitle " << setw(12) << p.plottitle << endl;
	outf << setw(20) << "randomseed " << setw(12) << p.randomSeed << endl;
	outf << setw(20) << "logevery " << setw(12) << p.logEvery << endl;
	outf << setw(20) << "fatlog " << setw(12) << p.fatlog << endl;
	outf << setw(20) << "saveevery " << setw(12) << p.saveEvery << endl;
	outf << setw(20) << "nindivs " << setw(12) << p.nindivs << endl;
	outf << setw(20) << "holdover " << setw(12) << p.holdover << endl;
	outf << setw(20) << "gammashape " << setw(12) << setprecision(DEF_PRECISION) << p.gammaShapeBrlen << endl;
	return outf;
}

/*	note, data is just a pointer to some already allocated HKYData (i.e. it isn't "owned" by an instance of Parameter) so we don't need to serialize it.
	rng has no dynamic data, so there is no need to serialize it.
*/
int Parameters::Serialize(char*& buf, int* size_)	{
	int& size = *size_;
	int len = starting_tree.stringlen() + 1;  // don't forget the null terminator
	size = sizeof(Parameters) + sizeof(len) + len;
	char* p = buf = new char[size];

	memcpy(p, this, sizeof(Parameters));
	p += sizeof(Parameters);

	memcpy(p, &len, sizeof(len));
	p += sizeof(len);

	memcpy(p, starting_tree.getstr(), len);
	p += len;

	assert(p - buf == size);
	return 0;
}

int Parameters::Deserialize(const char* buf)	{
	int size;
	const char* p = buf;

	memcpy(this, p, sizeof(Parameters));
	p += sizeof(Parameters);

	memcpy(&size, p, sizeof(int));
	p += sizeof(int);

	starting_tree.reallocate();  // cjb - here's our cheap hack!
	starting_tree = p;
	p += size;

	assert(size == starting_tree.stringlen() + 1);
	assert(p == buf + size + sizeof(Parameters) + sizeof(int));

	return 0;
}

Parameters& Parameters::operator=(const Parameters& rhs)	{
	memcpy(this, &rhs, sizeof(Parameters));
	strncpy(logfname, rhs.logfname, 85);
	strncpy(startfname, rhs.startfname, 85);
#if defined( ALLOW_USERTREES )
	strncpy(testfname, rhs.testfname, 85);
#endif
	strncpy(treefname, rhs.treefname, 85);
	strncpy(gnufname, rhs.gnufname, 85);
	strncpy(datafname, rhs.datafname, 81);
	strncpy(ofprefix, rhs.ofprefix, 81);
	strncpy(plottitle, rhs.plottitle, 81);
	strncpy(statefname, rhs.statefname, 81);
	return *this;
}

bool Parameters::operator==(const Parameters& rhs) const	{

	if (this == &rhs)
		return true;

	if (starting_tree != rhs.starting_tree)
		return false;

	// cheap hack

	NxsString temp1 = starting_tree;
	NxsString temp2 = rhs.starting_tree;

	const_cast<NxsString&>(starting_tree) = rhs.starting_tree;

	bool rv = (memcmp(this, &rhs, sizeof(NxsString)) == 0) ? true : false;

	const_cast<NxsString&>(starting_tree) = temp1;
	const_cast<NxsString&>(rhs.starting_tree) = temp2;

	return rv;
}



