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
//
//	NOTE: Portions of this source adapted from:
//	Press, W. H., B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling.  1992.
//	Numerical Recipes in C : The Art of Scientific Computing.  Cambridge University Press, Cambridge.

#if defined(_MSC_VER)
//POL 23-Feb-2006 VC doesn't have this header, and it was not needed to compile
//#	include <unistd.h>
#else
#	include <unistd.h>
#endif

#include "defs.h"
#include "funcs.h"
#include "population.h"
#include "tree.h"
#include "outputman.h"
#include "garlireader.h"

extern OutputManager outman;

#undef ROOT_OPT
#define FOURTH_ROOT

#define LOG_MIN_BRLEN log(min_brlen)

bool FloatingPointEquals(const FLOAT_TYPE first, const FLOAT_TYPE sec, const FLOAT_TYPE epsilon){
	FLOAT_TYPE diff = fabs(first - sec);
	return (diff < epsilon);
	}

//this is for sticking info about what is defined into log files, for later checking
void OutputImportantDefines(){
	outman.DebugMessage("#####\nThe following are/are not defined:");

#ifdef RESCALE_ARRAY_LENGTH
	outman.DebugMessage("RESCALE_ARRAY_LENGTH = %d", RESCALE_ARRAY_LENGTH);
#endif

	outman.DebugMessageNoCR("LUMP_LIKES : ");
#ifdef LUMP_LIKES
	outman.DebugMessage("%d", LUMP_FREQ);
#else
	outman.DebugMessage("no");
#endif

#ifdef DEBUG_SCORES
	outman.DebugMessage("DEBUG_SCORES");
#endif

#ifdef OPT_DEBUG
	outman.DebugMessage("OPT_DEBUG");
#endif

#ifdef VARIABLE_OPTIMIZATION
	outman.DebugMessage("VARIABLE_OPTIMIZATION");
#endif

#ifdef NO_EVOLUTION
	outman.DebugMessage("NO_EVOLUTION");
#endif

#ifdef SWAP_BASED_TERMINATION
	outman.DebugMessage("SWAP_BASED_TERMINATION");
#endif

#ifdef MORE_DETERM_PARAM_OPT
	outman.DebugMessage("MORE_DETERM_PARAM_OPT = yes");
#else
	outman.DebugMessage("MORE_DETERM_PARAM_OPT = no");
#endif

#ifdef ADAPTIVE_BOUNDED_OPT
	outman.DebugMessage("ADAPTIVE_BOUNDED_OPT = yes");
#else
	outman.DebugMessage("ADAPTIVE_BOUNDED_OPT = no");
#endif

#ifdef PUSH_TO_MIN_BLEN
	outman.DebugMessage("PUSH_TO_MIN_BLEN = yes");
#else
	outman.DebugMessage("PUSH_TO_MIN_BLEN = no");
#endif

#ifdef DEBUG_MESSAGES
	outman.DebugMessage("DEBUG_MESSAGES = yes");
#else
	outman.DebugMessage("DEBUG_MESSAGES = no");
#endif

#ifdef BOUND_DIGITS
	outman.DebugMessage("BOUND_DIGITS = %d", BOUND_DIGITS);
#else
	outman.DebugMessage("BOUND_DIGITS = default");
#endif

#ifdef OPT_BOUNDED_RESTORE
	outman.DebugMessage("OPT_BOUNDED_RESTORE = yes");
#else
	outman.DebugMessage("OPT_BOUNDED_RESTORE = no");
#endif

#ifdef FINAL_RESTORE_BLENS
	outman.DebugMessage("FINAL_RESTORE_BLENS = yes");
#else
	outman.DebugMessage("FINAL_RESTORE_BLENS = no");
#endif

	outman.DebugMessage("#####\n");
	}

#ifdef BROOK_GPU
#include <brook/brook.hpp>
//#include <brook/profiler.hpp>
//using namespace brook::internal;

void  BranchLike2 (::brook::stream des,
		::brook::stream res,
		::brook::stream pmat);

void  SecondBranchLike (::brook::stream des,
		::brook::stream part,
		::brook::stream res,
		::brook::stream pmat);

void  Product4 (::brook::stream des1,
		::brook::stream des2,
		::brook::stream res);

		brook::stream LCLstream(brook::getStreamType(( float4  *)0), 890, -1);
		brook::stream RCLstream(brook::getStreamType(( float4  *)0), 890, -1);
		brook::stream deststream(brook::getStreamType(( float4  *)0), 890, -1);
//		brook::stream tempstream(brook::getStreamType(( float4  *)0), 890, -1);
//		brook::stream tempstream2(brook::getStreamType(( float4  *)0), 890, -1);
		brook::stream Lprstream(brook::getStreamType(( float  *)0), 16, -1);
		brook::stream Rprstream(brook::getStreamType(( float  *)0), 16, -1);
#endif


//a variety of functions that don't belong to any class

#if defined(SINGLE_PRECISION_FLOATS) && (!defined(_MSC_VER)) || (defined(BOINC) && defined (_WIN32))
//Overloaded versions of min and max that take different types for the two arguments
//This should not be used in hot code when possible, and conditional comp should
//be used to make two different versions of the code
float min(const double first, const float second) {return min((float) first, second);}
float min(const float first, const double second) {return min(first, (float) second);}
float max(const double first, const float second) {return max((float) first, second);}
float max(const float first, const double second) {return max(first, (float) second);}
#endif

int FileExists( const char* s )
{

#ifdef POWERMAC_VERSION	// cjb

	if( s && access( s, 0 ) == 0 )
		return 1;
#else
	if (s)	{
		ifstream test(s);
		if(test.good())	{
			test.close();
			return 1;
		}
	}
#endif

	return 0;
}
/*
// GetRestartParams extracts the following pieces of information from
// the first line of the state file:
//	1. the last generation done on previous run (prev_generations)
//		- we want to start this run with generation prev_generations
//	2. the time recorded after the last generation was completed (prev_time)
//		- we'll start counting seconds with prev_time rather than zero
//	3. the random number seed with which to begin the run (randomSeed)
//		- we thus start exactly where we left off before
//
void GetRestartParams( Parameters& params )
{
	if( !FileExists( params.statefname ) ) throw ErrorException("Error opening state file: %s", params.statefname);

	ifstream sf( params.statefname );
	sf >> params.prev_generations >> params.prev_time >> params.randomSeed;
	sf.close();

	rnd.set_seed( params.randomSeed );
}
*/

int GetToken( FILE *in, char* tokenbuf, int maxlen){
	int ok = 1;

	int i;
	char ch = ' ';

	// skip leading whitespace
	while( in && ( isspace(ch) || ch == '[' ) ){
		ch = getc(in);
		}
	if( !in ) return 0;

	tokenbuf[0] = ch;
	tokenbuf[1] = '\0';
	tokenbuf[maxlen-1] = '\0';

	for( i = 1; i < maxlen-1; i++ ) {
		ch = getc(in);
		if( isspace(ch) || ch == ']' )
			break;
		tokenbuf[i] = ch;
		tokenbuf[i+1] = '\0';
	}

	if( i >= maxlen-1 )
		ok = 0;

	return ok;
}

bool FileIsNexus(const char *name){
	if (!FileExists(name))	{
		throw ErrorException("could not open file: %s!", name);
		}

	bool nexus = false;
	FILE *inf;
#ifdef BOINC
	inf = boinc_fopen(name, "r");
#else
	inf = fopen(name, "r");
#endif
	char buf[1024];
	GetToken(inf, buf, 1024);
	if(!(_stricmp(buf, "#NEXUS"))) nexus = true;

	fclose(inf);
	return nexus;
	}

bool FileIsFasta(const char *name){
	if (!FileExists(name))	{
		throw ErrorException("could not open file: %s!", name);
		}

	bool fasta = false;
	FILE *inf;
#ifdef BOINC
	inf = boinc_fopen(name, "r");
#else
	inf = fopen(name, "r");
#endif
	char buf[1024];
	GetToken(inf, buf, 1024);
	if(buf[0] == '>') fasta = true;

	fclose(inf);
	return fasta;
	}

/* //the ReadData within the GarliReader should now be used to read all data files
bool ReadData(const char* filename, SequenceData* data)	{
	bool usedNCL = false;

	if (!FileExists(filename))	{
		throw ErrorException("data file not found: %s!", filename);
		}

	if(FileIsNexus(filename)){
		outman.UserMessage("Attempting to read data file in Nexus format (using NCL): %s ...", filename);
		GarliReader &reader = GarliReader::GetInstance();
		int err = reader.HandleExecute(filename, true);
		if(err) throw ErrorException("Problem reading nexus datafile");
		//moving error checking and finding of correct char block into individual CreateMatrix functions
	//	NxsCharactersBlock *chars = reader.GetCharactersBlock();
	//	if(modSpec.IsAminoAcid() && modSpec.IsCodonAminoAcid()==false && chars->GetDataType() != NxsCharactersBlock::protein)
	//		throw ErrorException("protein data specified, but nexus file does not contain protein data!");
		data->CreateMatrixFromNCL(reader);
		usedNCL = true;
		}
	else if(FileIsFasta(filename)){
		outman.UserMessage("Attempting to read data file in Fasta format: %s ...", filename);
		data->ReadFasta(filename);
		}
	else{
		outman.UserMessage("Attempting to read data file in Phylip format: %s ...", filename);
		data->ReadPhylip(filename);
		}
*/
/*
	if(modSpec.IsCodon()){
		assert(0);
		if(modSpec.IsVertMitoCode()){
			static_cast<CodonData*>(data)->SetVertMitoCode();
			}
		static_cast<CodonData*>(data)->FillCodonMatrix(false);
		}
	else if(modSpec.IsCodonAminoAcid()){
		assert(0);
		if(modSpec.IsVertMitoCode()){
			static_cast<CodonData*>(data)->SetVertMitoCode();
			}
		static_cast<CodonData*>(data)->SetAminoAcid();
		static_cast<CodonData*>(data)->FillCodonMatrix(true);
		}
*/
/*
	// report summary statistics about the data
	data->Summarize();
	outman.UserMessage("\nData summary:");
	outman.UserMessage(" %d taxa", data->NTax());
	outman.UserMessage(" %d total characters.", data->NChar());
	outman.UserMessage("   %d constant characters.", data->NConstant());
	outman.UserMessage("   %d parsimony-informative characters.", data->NInformative());
	outman.UserMessage("   %d autapomorphic characters.", data->NAutapomorphic());
	//int total = data->NConstant() + data->NInformative() + data->NAutapomorphic();
	outman.flush();

	//if(modSpec.IsNucleotide()){
	if(1){
		// try to compress
		if (!data->Dense())	{
			outman.UserMessage("Compressing data matrix...");
			data->Collapse();
			outman.UserMessage("%d columns in data matrix after compression.", data->NChar());
			}
		else {
			outman.UserMessage("Datafile already compressed.");
			outman.UserMessage("%d columns in compressed data matrix.\n", data->NChar());
			}

		if(modSpec.IsNucleotide()){
			data->DetermineConstantSites();
//			if(!data->Dense()) data->Save(filename, "new");
			}
		else if(modSpec.IsAminoAcid()){
//			static_cast<CodonData*>(data)->DetermineConstantAASites();
			}
		}
*/
/*	return usedNCL;
	}
*/

int ReadData(GeneralGamlConfig *conf, NucleotideData* data)	{
	assert(0);
	// regurgitate params specified
/*	if( params.restart ) {
		outman.UserMessage("Restarting using state file \"%s\"", params.statefname);
		GetRestartParams( const_cast<Parameters&>(params) );
		outman.UserMessage("random number seed set to %d", params.randomSeed);
		outman.UserMessage("last generation from previous run was %d", params.prev_generations);
		outman.UserMessage("starting with previous elapsed time, which was %d seconds", params.prev_time);
		}
*/
//	const_cast<Parameters&>(params).BriefReport( cout );
//	outman.UserMessage("");

	// Check to be sure data file exists
	//
/*	if( !FileExists( conf->datafname.c_str() ) ) throw ErrorException("data file does not exist: %s", conf->datafname.c_str());

	// Read in the data matrix
	outman.flush();
	outman.UserMessage("Reading data file %s...", conf->datafname.c_str());
	data->Read( conf->datafname.c_str() );

	// report summary statistics about data
	data->Summarize();
	outman.UserMessage(" %d constant characters.", data->NConstant());
	outman.UserMessage(" %d parsimony-informative characters.", data->NInformative());
	outman.UserMessage(" %d autapomorphic characters.", data->NAutapomorphic());
	int total = data->NConstant() + data->NInformative() + data->NAutapomorphic();
	outman.UserMessage(" %d total characters.", total);
	outman.flush();

	//DZ Only compress and write data to file if dense=0 (data is not already compressed)
	if(!(data->Dense())){
		outman.UserMessage("Compressing data file...");
		data->Collapse();
		data->Save("compdata.nex", "new");
		outman.UserMessage("  %d columns in data matrix after compression", data->NChar());
		}
	else outman.UserMessage("Datafile already compressed.\n %d columns in compressed data matrix", data->NChar());

	data->DetermineConstantSites();
	outman.UserMessage("");
	outman.flush();
*/	return 0;

}

int RandomInt(int lb, int ub)	{
	return lb + rand() % (ub-lb+1);
}

FLOAT_TYPE RandomFrac()	{
	return (FLOAT_TYPE) (rand() / RAND_MAX);
}

FLOAT_TYPE RandomDouble(FLOAT_TYPE lb, FLOAT_TYPE ub)	{
	return lb + RandomFrac() * (ub - lb);
}


//Bracketing func from Numerical Recipies
//attempts to use parabolic fit to bracket, otherwise
//uses golden section
#define GOLD 1.618034
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define GLIMIT 100.0
//#define GLIMIT 1.6
#define ITMAX 50
#define TINY 1.0e-20
#define SHFT(a,b,c,d) a=b;b=c;c=d;
#define SIGN(a,b) ((b)>ZERO_POINT_ZERO ? fabs(a) : -fabs(a))
#define FMAX(a,b) ((a)>(b) ? (a):(b))

//This version takes a node pointer and optimizes blens
int mnbrak(FLOAT_TYPE *ax, FLOAT_TYPE *bx, FLOAT_TYPE *cx, FLOAT_TYPE *fa, FLOAT_TYPE *fb, FLOAT_TYPE *fc, FLOAT_TYPE (*func)(TreeNode*, Tree*, FLOAT_TYPE), TreeNode *thisnode, Tree *thistree){
	FLOAT_TYPE ulim, u, r, q, fu;

//	ofstream brak("brakdebug.log", ios::app);
//	brak.precision(10);
//	brak << "node " << thisnode->nodeNum << "\n";

	*fa=(*func)(thisnode, thistree, *ax);
	*fb=(*func)(thisnode, thistree, *bx);
	*fc=(*func)(thisnode, thistree, *cx);

	//hopefully we passsed in a good bracket.  If so, get out.
	if(*fb < *fa && *fb < *fc) return 0;

	/*

	if(*fb > *fa){
		SHFT(dum, *ax, *bx, dum)
		SHFT(dum, *fb, *fa, dum)
		}

	*cx=(*bx)+GOLD*(*bx-*ax);
	*cx = (*cx > min_brlen ? (*cx < DEF_MAX_BRLEN ? *cx : DEF_MAX_BRLEN) : min_brlen);
	*fc=(*func)(thisnode, thistree, *cx);
*/
	while(*fb>*fc){

		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/(FLOAT_TYPE)(2.0*SIGN(FMAX(fabs(q-r),TINY), q-r));
		u = (FLOAT_TYPE)(u > DEF_MIN_BRLEN ? (u < DEF_MAX_BRLEN ? u : DEF_MAX_BRLEN) : DEF_MIN_BRLEN);
		ulim=(FLOAT_TYPE)((*bx)+GLIMIT*(*cx-*bx));
		if((*bx-u)*(u-*cx)>ZERO_POINT_ZERO){
			fu=(*func)(thisnode, thistree, u);
			if(fu < *fc){
				*ax=*bx;
				*bx=u;
				*fa=*fb;
				*fb=fu;
				return 0;
				}
			else if(fu > *fb){
				*cx=u;
				*fc=fu;
				return 0;
				}
			u=(FLOAT_TYPE)((*cx)+GOLD*(*cx-*bx));
			//DZ 10/27/03 don't let this evaluate totally insane blens
/*			if(u>=.69){ //=ln(2)
				if(max(*ax, max(*bx, *cx)) < .69) u=.69;
				else{
					 u=max(*ax, max(*bx, *cx)) + .69;
					 }
				limited=true;
				}
*/			fu=(*func)(thisnode, thistree, u);
			}
		else if((*cx-u)*(u-ulim)>ZERO_POINT_ZERO){
			//DZ 10/27/03 don't let this evaluate totally insane blens
	/*		if(u>=.69){ //=ln(2)
				if(max(*ax, max(*bx, *cx)) < .69) u=.69;
				else{
					u=max(*ax, max(*bx, *cx)) + .69;
					}
				limited=true;
				}
	*/		fu=(*func)(thisnode, thistree, u);
			if(fu <*fc){
				SHFT(*bx, *cx, u, *cx+(FLOAT_TYPE)GOLD*(*cx-*bx));
				SHFT(*fb, *fc, fu, (*func)(thisnode, thistree, u));
				}
			}
		else if((u-ulim)*(ulim-*cx) >ZERO_POINT_ZERO){
			u=ulim;
			//DZ 10/27/03 don't let this evaluate totally insane blens
	/*		if(u>=.69){ //=ln(2)
				if(max(*ax, max(*bx, *cx)) < .69) u=.69;
				else{
					u=max(*ax, max(*bx, *cx)) + .69;
					}
				limited=true;
				}
	*/		fu=(*func)(thisnode, thistree, u);
			}
		else{
			u=(*cx)+(FLOAT_TYPE)GOLD*(*cx-*bx);
			fu=(*func)(thisnode, thistree, u);
			}
		SHFT(*ax, *bx, *cx, u)
		SHFT(*fa, *fb, *fc, fu)
	/*	if(((*ax < -18.42) && (*bx < -18.42)) || ((*ax<-10) && (*bx<-10) && (*cx>1))){
			//DZ 12-18-03 if our three best points are all < ln(1e-8), just give up and take that as a blen
			//the MLE is probably essentially 0.  Note that sometimes when ax and bx are very small this
			//func tries very large values for cx, which I think is a bug.  This hack also avoids that
			return 1;
			}
	*/	}
	return 0;
	}


//This version takes a node pointer and optimizes blens
FLOAT_TYPE brent(FLOAT_TYPE ax, FLOAT_TYPE bx, FLOAT_TYPE cx, FLOAT_TYPE (*f)(TreeNode *, Tree*, FLOAT_TYPE), FLOAT_TYPE tol, FLOAT_TYPE *xmin, TreeNode *thisnode, Tree *thistree){
	 int iter;
	 FLOAT_TYPE a, b, d, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
	 FLOAT_TYPE e=ZERO_POINT_ZERO;

	 a=(ax < cx ? ax : cx); //make a the smallest of the three bracket points
	 b=(ax > cx ? ax : cx); //and b the largest
	 x=w=v=bx;				//make x the current minimum, as well as w and v

	 fw=fv=fx=(*f)(thisnode, thistree, x);

	 for(iter=1;iter<=ITMAX;iter++){
	 	xm=ZERO_POINT_FIVE*(a+b);		//xm is the midpoint of the bracket (of a and b)

	 	tol2=(FLOAT_TYPE)(2.0*(tol1=(FLOAT_TYPE)(tol*fabs(x)+ZEPS)));

	 	if (fabs(x-xm) <= (tol2-ZERO_POINT_FIVE*(b-a))){ //termination condition
	 		*xmin=x;							//if the distance between x and bracket mean is <
	 		return fx;
	 		}
	 	if (fabs(e) > tol1){	//construct a trial parabolic fit
	 		r=(x-w)*(fx-fv);
	 		q=(x-v)*(fx-fw);
	 		p=(x-v)*q-(x-w)*r;
	 		q=(FLOAT_TYPE)(2.0*(q-r));
	 		if(q>ZERO_POINT_ZERO) p=-p;
	 		q=fabs(q);
	 		etemp=e;
	 		e=d;
	 		if(fabs(p) >= fabs(ZERO_POINT_FIVE*q*etemp)||p<=q*(a-x) || p>=q*(b-x)) //determine if the parabolic fit is good
	 			d=(FLOAT_TYPE)(CGOLD*(e=(x>=xm?a-x:b-x)));  //if not

	 		else{				//if so, take the parabolic step
	 			d=p/q;
	 			u=x+d;
	 			if(u-a < tol2||b-u<tol2)
	 				d=SIGN(tol1,xm-x);
	 			}
	 		}
	 	else{
	 		d=(FLOAT_TYPE)(CGOLD*(e=(x>=xm?a-x:b-x))); //e is the distance moved in the step before last
	 		}							 //d is golden section of that (.38.... times)

	 	u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));//u is the next point to be evaluated
	 												//it is x+d or ?
	 	fu=(*f)(thisnode, thistree, u);
	 	if(fu<=fx){						//if our new try at u is better than the previous min at x

	 		if(u>=x) a=x; else b=x;		//if u is > x, x becomes a, otherwise it becomes b
	 		SHFT(v,w,x,u);				//w becomes v, x becomes w and u becomes x
	 		SHFT(fv,fw,fx,fu);
	 		}

	 	else{							//if our new try at u is worse than the old min at x

	 		if(u<x) a=u; else b=u;		//if u is < x, u becomes a, otherwise it becomes b

	 		if(fu<=fw||w==x){			//if the score at u is < the score of the last attempt or x was the last attempt
	 			v=w;					//w is the second best point
	 			w=u;					//v is the third best point
	 			fv=fw;
	 			fw=fu;
	 			}
	 		else if(fu<=fv||v==x||v==w){
	 			v=u;
	 			fv=fu;
	 			}
	 		}
	 	}
	 *xmin=x;
	 return fx;
	 }

#ifdef FOURTH_ROOT
#define effectiveMin 0.01
#define effectiveMax 1.77827941
#define sweetspot 0.08
#define smallShift 0.02
#elif ROOT_OPT
#define effectiveMin 0.0001
#define effectiveMax 3.16227766
#define sweetspot 0.0016
#define smallShift 0.0004
#else
#define effectiveMin min_brlen
#define effectiveMax max_brlen
#define sweetspot 0.00000256
#define smallShift 0.00000016
#endif


//My version of the bracketing fuction that can abort and stop evaluating under certain conditions
//it assumes that the raw branch lengths are being passed in (not logs) so values are bounded
//by the minimum branch length
int DZbrak(FLOAT_TYPE *worstOuter, FLOAT_TYPE *mid, FLOAT_TYPE *bestOuter, FLOAT_TYPE *worstOuterL, FLOAT_TYPE *midL, FLOAT_TYPE *bestOuterL, FLOAT_TYPE (*func)(TreeNode*, Tree*, FLOAT_TYPE, bool), TreeNode *thisnode, Tree *thistree){
	//points are always passed in such that worstOuter < mid < bestOuter
	FLOAT_TYPE nextTry, r, q, nextTryL, dum;
	bool possibleZeroMLE=false;

	*worstOuterL=(*func)(thisnode, thistree, *worstOuter, true);
	if(*midL<0)
		*midL=(*func)(thisnode, thistree, *mid, true);

	if(*midL>*worstOuterL){//the min must be to the left (or maybe between) of our evals, so don't bother
		//evaluating the bestOuter we passed in, which we know is to the right. Either evaluate the min,
		//or if worstOuter already is the the min evaluate a point between the current evals
		if(*worstOuter==effectiveMin){
			SHFT(dum, *bestOuter, *worstOuter, *mid)
			SHFT(dum, *bestOuterL, *worstOuterL, *midL)
			*mid=(FLOAT_TYPE)((*worstOuter+*bestOuter)*ZERO_POINT_FIVE);
			*midL=(*func)(thisnode, thistree, *mid, true);
			if(*bestOuterL < *midL && !(*mid > sweetspot)) return 1;
			}
		else if(!(*worstOuter > sweetspot)){
			SHFT(dum, *mid, *worstOuter, dum)
			SHFT(dum, *midL, *worstOuterL, dum)
			*bestOuter=(FLOAT_TYPE)effectiveMin;
			*bestOuterL=(*func)(thisnode, thistree, *bestOuter, true);
			if(*bestOuterL < *midL && !(*mid > sweetspot)) return 1;
			}
		else{
			*bestOuter=(FLOAT_TYPE)(sweetspot-.02);
			possibleZeroMLE=true;
			SHFT(dum, *worstOuter, *mid, dum)
			SHFT(dum, *worstOuterL, *midL, dum)
			*bestOuterL=(*func)(thisnode, thistree, *bestOuter, true);
			}
		}
	else
		*bestOuterL=(*func)(thisnode, thistree, *bestOuter, true);

	/*There are a pretty limited number of cases for each loop here
	case 1: We have three points that define a bracket -> exit with 0
	case 2: We have three points with sucessively better scores to the right
			2a: The curvature implied by the three points is convex or only slightly concave which makes
				the parabolic estimate very poor. -> Do a GOLD step.
			2b: The parabolic estimate is between mid and bestOuter -> take parabolic step
			2c: The parabolic estimate is to the right of bestOuter -> ?

	case 3: We have three points with sucessively better scores to the left
			3a: The best score is at the minimum allowed value, suggesting a possible zero MLE
				3a1: Parabolic estimate is between mid and minimum value.  ??
				3a2: Parabolic estimate is < minimum -> return zeroMLE=true
			3b: The parabolic estimate is between mid and bestOuter -> take parabolic step
			3c: The parabolic estimate is to the left of bestOuter  -> take parabolic step
	*/

	if(*worstOuterL < *bestOuterL){
		SHFT(dum, *worstOuter, *bestOuter, dum)
		SHFT(dum, *worstOuterL, *bestOuterL, dum)
		}

	do{
		if(*midL < *worstOuterL && *midL < *bestOuterL){//case 1, got a bracket
			if(*bestOuter==effectiveMin && (*worstOuter - *mid)> .2){
				//nextTry=(*mid+*worstOuter)*.5;
				nextTry=(FLOAT_TYPE)0.16;
				nextTryL=(*func)(thisnode, thistree, nextTry, true);
				assert(nextTryL < *worstOuterL);
				if(nextTryL < *midL){
					SHFT(*bestOuter, *mid, nextTry, dum)
					SHFT(*bestOuterL, *midL, nextTryL, dum)
					}
				else if(nextTryL < *bestOuterL){
					SHFT(*worstOuter, *bestOuter, nextTry, dum)
					SHFT(*worstOuterL, *bestOuterL, nextTryL, dum)
					}
				else{
					*worstOuter=nextTry;
					*worstOuterL=nextTryL;
					}
				}
			return 0;
			}
		else{
			FLOAT_TYPE diffMidBestL=(*midL-*bestOuterL);
			FLOAT_TYPE diffMidBest=(*mid-*bestOuter);
			FLOAT_TYPE diffMidWorstL=(*worstOuterL-*midL);
			FLOAT_TYPE diffMidWorst=(*worstOuter-*mid);
			if(*worstOuter < *bestOuter){ //case 2
				//check the curvature

				FLOAT_TYPE slopeRatio=(diffMidBestL/diffMidBest) / (diffMidWorstL/diffMidWorst);
				if(slopeRatio > 0.9){
					nextTry=(FLOAT_TYPE)((*bestOuter)+2.0*(*bestOuter-*mid));//case 2a
					}
				else{ //case 2b and 2c
					r=diffMidWorst*diffMidBestL;
					q=diffMidBest*diffMidWorstL;
					nextTry=(FLOAT_TYPE)((*mid)-(diffMidBest*q-(*mid-*worstOuter)*r)/(2.0*SIGN(FMAX(fabs(q-r),TINY), q-r)));
					if(/*nextTry > *bestOuter && */fabs(nextTry-*bestOuter) < smallShift){
						//if the parabolic estimate is very near our current best it tends to take
						//a while to get the bracket, so just push it a little further to the right
						nextTry += (FLOAT_TYPE)smallShift;
						}
					}
				}

			else if(*worstOuter > *bestOuter){ //case 3
				r=diffMidWorst*diffMidBestL;
				q=diffMidBest*diffMidWorstL;
				nextTry=(FLOAT_TYPE)((*mid)-(diffMidBest*q-diffMidWorst*r)/(2.0*SIGN(FMAX(fabs(q-r),TINY), q-r)));
				if(*bestOuter==effectiveMin){
					if(nextTry < effectiveMin || *mid < sweetspot || possibleZeroMLE) return 1; //case 3a2
					else {//case 3a1
						//just go with the parabolic
						}
					}
				else {
					if(possibleZeroMLE==true) nextTry=(FLOAT_TYPE)effectiveMin;
					}
				}
			}
		assert(nextTry >= effectiveMin);
		nextTryL=(*func)(thisnode, thistree, nextTry, true);
		if(nextTryL < *bestOuterL){
			if((*mid-nextTry) * (nextTry-*bestOuter)>ZERO_POINT_ZERO){//if the proposed point is between mid and bestOuter
				SHFT(dum, *worstOuter, *mid, nextTry);
				SHFT(dum, *worstOuterL, *midL, nextTryL);
				}
			else{
				SHFT(*worstOuter, *mid, *bestOuter, nextTry);
				SHFT(*worstOuterL, *midL, *bestOuterL, nextTryL);
				}
			}
		else{
			if((*mid-nextTry)*(nextTry-*bestOuter)>ZERO_POINT_ZERO){//if the proposed point is between mid and bestOuter
				assert(nextTryL < *midL);//if this isn't the case there are multiple optima
				SHFT(dum, *worstOuter, *mid, nextTry);
				SHFT(dum, *worstOuterL, *midL, nextTryL);
				}
			else{
				if(nextTryL < *midL){
					SHFT(*worstOuter, *mid, *bestOuter, nextTry);
					SHFT(*worstOuterL, *midL, *bestOuterL, nextTryL);
					}
				else{
					SHFT(dum, *mid, *bestOuter, dum);
					SHFT(dum, *midL, *bestOuterL, dum);
					*worstOuter=nextTry;
					*worstOuterL=nextTryL;
					}
				}
			}
		}while(1);
	assert(0);
	return 0;
}

//I'm reworking this a bit to better use the information that has already been generated in the bracketing function
//since we already have those function evaluations, we might as well pass them in and use them
FLOAT_TYPE DZbrent(FLOAT_TYPE ax, FLOAT_TYPE bx, FLOAT_TYPE cx, FLOAT_TYPE fa, FLOAT_TYPE fx, FLOAT_TYPE fc, FLOAT_TYPE (*f)(TreeNode *, Tree*, FLOAT_TYPE, bool), FLOAT_TYPE tol, FLOAT_TYPE *xmin, TreeNode *thisnode, Tree *thistree){
	 int iter;
 	 FLOAT_TYPE a, b, d, etemp, fu, fv, fw/*, fx*/, p, q, r, tol1, tol2, u, v, w, x, xm;
	 FLOAT_TYPE e=ZERO_POINT_ZERO;

	 if((fx<fa && fx<fc)==false){
	 	//if bx isn't the current minimum, make is so
	 	if(fa<fx){
	 		FLOAT_TYPE dummy=fa;
	 		fa=fx;
	 		fx=dummy;
	 		dummy=ax;
	 		ax=bx;
	 		bx=dummy;
	 		}
	 	else if(fc<fx){
	 		FLOAT_TYPE dummy=fc;
	 		fc=fx;
	 		fx=dummy;
	 		dummy=cx;
	 		cx=bx;
	 		bx=dummy;
	 		}
	 	}
	 assert(fx<fa && fx<fc);

	 FLOAT_TYPE paraMinlnL, paraErr, paraErrCrit;
	 paraErrCrit=(tol<.5 ? tol*10 : 5);
	 bool paraOK=false, para=false;

//	ofstream brak("brakdebug.log", ios::app);
//	brak << "node " << thisnode->nodeNum << "\n";

	 a=(ax < cx ? ax : cx); //make a the smallest of the three bracket points
	 b=(ax > cx ? ax : cx); //and b the largest
	if(ax>cx){
		FLOAT_TYPE dummy=fa;
		fa=fc;
		fc=dummy;
		}

	 x=bx;				//make x the current minimum

 	if(fa<fc){
 		w=a;
 		fw=fa;
 		v=b;
 		fv=fc;
 		}
 	else{
 		v=a;
 		fv=fa;
 		w=b;
 		fw=fc;
 		}

	xm=(FLOAT_TYPE)ZERO_POINT_FIVE*(a+b);       //xm is the midpoint of the bracket (of a and b)
	e=(x>=xm?a-x:b-x);	//set e to the larger of the two bracket intervals
	d=(FLOAT_TYPE)CGOLD*e;

//	assert(a<=x && x<=b);

//	 fw=fv=fx=(*f)(thisnode, thistree, x);

	 for(iter=1;iter<=ITMAX;iter++){

	 	xm=(FLOAT_TYPE)ZERO_POINT_FIVE*(a+b);		//xm is the midpoint of the bracket (of a and b)

	 	tol2=(FLOAT_TYPE)(2.0*(tol1=(FLOAT_TYPE)(tol*fabs(x)+ZEPS)));

/*	 	if (fabs(x-xm) <= (tol2-ZERO_POINT_FIVE*(b-a))){ //termination condition
	 		*xmin=x;
	 		return fx;
	 		}
*/// 	if (fabs(e) > tol1){	//construct a trial parabolic fit
	 		r=(x-w)*(fx-fv);
	 		q=(x-v)*(fx-fw);
	 		p=(x-v)*q-(x-w)*r;
	 		q=(FLOAT_TYPE)2.0*(q-r);
	 		if(q>ZERO_POINT_ZERO) p=-p;
	 		q=fabs(q);
	 		etemp=e;
	 		e=d;
	 		if(fabs(p) >= fabs(ZERO_POINT_FIVE*q*etemp)||p<=q*(a-x) || p>=q*(b-x)){ //determine if the parabolic fit is good
	 			d=(FLOAT_TYPE)(CGOLD*(e=(x>=xm?a-x:b-x)));  //if not
	 			u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
	 			}

	 		else{				//if so, take the parabolic step
	 			d=p/q;
	 			u=x+d;

	 			FLOAT_TYPE alph=w-u;
				FLOAT_TYPE beta=x-u;
				paraMinlnL=(((fx) * alph*alph) - ((fw) * beta*beta)) / (alph*alph - beta*beta);

	 			if(paraOK==true){
	 				//the estimation error in the parabolic step always seems to at least half each iteration,
	 				//hence the division by 2.0
	 				FLOAT_TYPE estlnL=(FLOAT_TYPE)(paraMinlnL - paraErr*ZERO_POINT_FIVE);
	 				if((fx - estlnL) < tol){
	 					*xmin=x;
	 					return fx;
	 					}
	 				}
	 			para=true;

	 			if(u-a < tol2||b-u<tol2)
	 				d=SIGN(tol1,xm-x);
	 			}
//	 		}
/*	 	else{
	 		d=CGOLD*(e=(x>=xm?a-x:b-x)); //e is the distance moved in the step before last
	 									 //d is golden section of that (.38.... times)
	 		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
	 		}
*/	 	//assert(a<=u && u<=b);
		//for some reason this occasionally proposes a new value
		//that is not within the bracket.  If that happens force
		//the new point to be within a and b
		if(!(a<=u && u<=b)){
			ofstream S("brakmiss.log", ios::app);
			S.precision(12);
			S << fa << "\t" << fw << "\t" << fv << "\t" << fu << endl;
			S << a << "\t" << b << "\t" << u << endl;
			S.close();
			u=a/b;
			}
	 	fu=(*f)(thisnode, thistree, u, false);
	 	if(para==true){
	 		paraErr=fu - paraMinlnL;
	 		if(fabs(paraErr) < paraErrCrit) paraOK=true;
	 		para=false;
	 		}

	 	if(!(fu>fx)){						//if our new try at u is better than the previous min at x

	 		if(u>=x) a=x; else b=x;		//if u is > x, x becomes a, otherwise it becomes b
	 		SHFT(v,w,x,u);				//w becomes v, x becomes w and u becomes x
	 		SHFT(fv,fw,fx,fu);
	 		}


		//DJZ 1/21/04 Rewrote this loop.  I think that it was buggy, and made no sense
		//to me previously.  Was updating variables such that w and v were not either a or b
		//(The bracket) but were older evaluations.  This resulted in terrible cases where
		//the interval {w, x, v} didn't even contain a minimum at all, but {a, x, b} did.
	 	else{						//if our new try at u is worse than the old min at x
	 		if(u<x){					//if the new try is to the left of the min
	 			a=u;
		 		if(w<x){					//if the left bracket was w
		 			w=u;						//u becomes the new w
					fw=fu;
		 			}
		 		else{						//if the left braket was v
		 			if(fu<fw){					//if fu is better than the old fw
		 				v=w;						//w becomes the new v
		 				fv=fw;
		 				w=u;						//u becomes the new w
		 				fw=fu;
		 				}
		 			else{						//if the old fw was better than fu
		 				v=u;						//u becomes the new v
		 				fv=fu;
		 				}
		 			}
	 			}
	 		else{
	 			b=u;				//if the new try is to the right of the min
	 			if(w<x){				//if the left bracket was w
	 				if(fu<fw){				//if fu is better than the old fw
	 					v=w;					//w becomes the new v
	 					fv=fw;
	 					w=u;					//u becomes the new w
	 					fw=fu;
	 					}
	 				else{					//if the old fw was better than fu
	 					v=u;						//u becomes the new v
	 					fv=fu;
	 					}
	 				}
	 			else {					//if the left bracket was v
	 				w=u;					//u becomes the new w
	 				fw=fu;
	 				}
	 			}
	 		}
/*
	 		if(u<x) a=u; else b=u;		//if u is < x, u becomes a, otherwise it becomes b

	 		if(fu<=fw||w==x){			//if the score at u is < the score of the last attempt or x was the last attempt
	 			v=w;					//w is the second best point
	 			w=u;					//v is the third best point
	 			fv=fw;
	 			fw=fu;
	 			}
	 		else if(fu<=fv||v==x||v==w){
	 			v=u;
	 			fv=fu;
	 			}
	 		}
*/	 	}
	 *xmin=x;
	 return fx;
	 }

void InferStatesFromCla(char *states, FLOAT_TYPE *cla, int nchar){
	//this function takes a cla that contains the contribution of the whole tree
	//and calculates the most probable state at each site.  The resulting array of
	//states is placed into *states, which should already be allocated.
	assert(0);//need to generalize this for n rates is it ever needs to be used again
	FLOAT_TYPE stateProbs[4];

	for(int c=0;c<nchar;c++){
		stateProbs[0]=stateProbs[1]=stateProbs[2]=stateProbs[3]=ZERO_POINT_ZERO;
		for(int i=0;i<4;i++){
			for(int j=0;j<4;j++){
				stateProbs[j] += *cla++;
				}
			}
		int max=0, next=1;

		while(next < 4){
			if(stateProbs[next] > stateProbs[max]){
				max=next;
				}
			next++;
			}

		states[c]=max;
		}
	}

vector<InternalState *> *InferStatesFromCla(FLOAT_TYPE *cla, int nchar, int nrates){
	FLOAT_TYPE stateProbs[4];

	vector<InternalState *> *stateVec = new vector<InternalState *>;

	for(int c=0;c<nchar;c++){
		stateProbs[0]=stateProbs[1]=stateProbs[2]=stateProbs[3]=ZERO_POINT_ZERO;
		for(int i=0;i<nrates;i++){
			for(int j=0;j<4;j++){
				stateProbs[j] += *cla++;
				}
			}
		InternalState *site=new InternalState(stateProbs);
		stateVec->push_back(site);


	//		out << bases[(stateProbs[max1] > stateProbs[max2] ? max1 : max2)] << "\t";
	//		out << stateProbs[0]/tot << "\t" << stateProbs[1]/tot << "\t" << stateProbs[2]/tot << "\t" << stateProbs[3]/tot << "\n";

		}
	return stateVec;
	}

FLOAT_TYPE CalculateHammingDistance(const char *str1, const char *str2, const int *counts, int nchar, int nstates){
	FLOAT_TYPE diff=0.0;
	int pos1=0, pos2=0;
	int effectiveChar=0;
	for(int i=0;i<nchar;i++){
		bool unambig1 = true;
		bool unambig2 = true;
		if(nstates == 4){
			if(str1[pos1] < 0){
				unambig1 = false;
				if(str1[pos1] == -4) pos1++;
				else{
    				int s=-str2[pos1++];
	    			for(int i=0;i<s;i++) pos1++;
					}
				}
			if(str2[pos2] < 0){
				unambig2 = false;
				if(str2[pos2] == -4) pos2++;
				else{
    				int s=-str2[pos2++];
	    			for(int i=0;i<s;i++) pos2++;
					}
				}
			}
		else{
			if(str1[pos1] == nstates){
				unambig1 = false;
				pos1++;
				}
			if(str1[pos2] == nstates){
				unambig2 = false;
				pos2++;
				}
			}
		if(unambig1 && unambig2){
			effectiveChar += counts[i];
			if(str1[pos1++] != str2[pos2++]) diff += (FLOAT_TYPE) counts[i];
			}
		}

	return diff/(FLOAT_TYPE)effectiveChar;
}

void SampleBranchLengthCurve(FLOAT_TYPE (*func)(TreeNode*, Tree*, FLOAT_TYPE, bool), TreeNode *thisnode, Tree *thistree){
	for(FLOAT_TYPE len=(FLOAT_TYPE)effectiveMin;len<(FLOAT_TYPE)effectiveMax;len*=2.0)
		(*func)(thisnode, thistree, len, true);
	}

