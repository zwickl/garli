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

extern OutputManager outman;

#undef ROOT_OPT
#define FOURTH_ROOT

#define LOG_MIN_BRLEN log(min_brlen)

//a variety of functions that don't belong to any class

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
int ReadData(const char* filename, HKYData* data)	{
	if (!FileExists(filename))	{
		throw ErrorException("data file not found: %s!", filename);
		return -1;
	}

	outman.UserMessage("Reading data file: %s...", filename);

	data->Read( filename );
	if((data->NChar() > 0) == false) throw ErrorException("problem reading data!");

	// report summary statistics about data
	data->Summarize();

	outman.UserMessage(" %d constant characters.", data->NConstant());
	outman.UserMessage(" %d parsimony-informative characters.", data->NInformative());
	outman.UserMessage(" %d autapomorphic characters.", data->NAutapomorphic());
	int total = data->NConstant() + data->NInformative() + data->NAutapomorphic();
	outman.UserMessage(" %d total characters.", total);
	outman.flush();

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
	data->DetermineConstantSites();
	if(!data->Dense()) data->Save(filename, "new");
	return 0;
	}

int ReadData(GeneralGamlConfig *conf, HKYData* data)	{

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
	if( !FileExists( conf->datafname.c_str() ) ) throw ErrorException("data file does not exist: %s", conf->datafname.c_str());

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
	return 0;

}

int RandomInt(int lb, int ub)	{
	return lb + rand() % (ub-lb+1);
}

FLOAT_TYPE RandomFrac()	{
	return (FLOAT_TYPE)rand() / RAND_MAX;
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
#define SIGN(a,b) ((b)>0.0 ? fabs(a) : -fabs(a))
#define FMAX(a,b) ((a)>(b) ? (a):(b))

//This version takes a node pointer and optimizes blens
int mnbrak(FLOAT_TYPE *ax, FLOAT_TYPE *bx, FLOAT_TYPE *cx, FLOAT_TYPE *fa, FLOAT_TYPE *fb, FLOAT_TYPE *fc, FLOAT_TYPE (*func)(TreeNode*, Tree*, FLOAT_TYPE), TreeNode *thisnode, Tree *thistree){
	FLOAT_TYPE ulim, u, r, q, fu;
	bool limited=false;
	
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
		if((*bx-u)*(u-*cx)>0.0){
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
		else if((*cx-u)*(u-ulim)>0.0){
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
		else if((u-ulim)*(ulim-*cx) >0.0){
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
	 FLOAT_TYPE e=0.0;
	 
	 a=(ax < cx ? ax : cx); //make a the smallest of the three bracket points 
	 b=(ax > cx ? ax : cx); //and b the largest
	 x=w=v=bx;				//make x the current minimum, as well as w and v
	 
	 fw=fv=fx=(*f)(thisnode, thistree, x);
	 
	 for(iter=1;iter<=ITMAX;iter++){
	 	xm=(FLOAT_TYPE)0.5*(a+b);		//xm is the midpoint of the bracket (of a and b)
	 	
	 	tol2=(FLOAT_TYPE)(2.0*(tol1=(FLOAT_TYPE)(tol*fabs(x)+ZEPS)));	
	 	
	 	if (fabs(x-xm) <= (tol2-0.5*(b-a))){ //termination condition
	 		*xmin=x;							//if the distance between x and bracket mean is < 
	 		return fx;
	 		}
	 	if (fabs(e) > tol1){	//construct a trial parabolic fit
	 		r=(x-w)*(fx-fv);
	 		q=(x-v)*(fx-fw);
	 		p=(x-v)*q-(x-w)*r;
	 		q=(FLOAT_TYPE)(2.0*(q-r));
	 		if(q>0.0) p=-p;
	 		q=fabs(q);
	 		etemp=e;
	 		e=d;
	 		if(fabs(p) >= fabs(0.5*q*etemp)||p<=q*(a-x) || p>=q*(b-x)) //determine if the parabolic fit is good
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
	bool tryPara=true, possibleZeroMLE=false;
	int numAttemptsWithBestAtMin=0;

	*worstOuterL=(*func)(thisnode, thistree, *worstOuter, true);
	if(*midL<0)
		*midL=(*func)(thisnode, thistree, *mid, true);

	if(*midL>*worstOuterL){//the min must be to the left (or maybe between) of our evals, so don't bother 
		//evaluating the bestOuter we passed in, which we know is to the right. Either evaluate the min,
		//or if worstOuter already is the the min evaluate a point between the current evals
		if(*worstOuter==effectiveMin){
			SHFT(dum, *bestOuter, *worstOuter, *mid)
			SHFT(dum, *bestOuterL, *worstOuterL, *midL)
			*mid=(FLOAT_TYPE)((*worstOuter+*bestOuter)*0.5);
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
			if((*mid-nextTry) * (nextTry-*bestOuter)>0.0){//if the proposed point is between mid and bestOuter
				SHFT(dum, *worstOuter, *mid, nextTry);
				SHFT(dum, *worstOuterL, *midL, nextTryL);
				}
			else{
				SHFT(*worstOuter, *mid, *bestOuter, nextTry);
				SHFT(*worstOuterL, *midL, *bestOuterL, nextTryL);
				}
			}
		else{
			if((*mid-nextTry)*(nextTry-*bestOuter)>0.0){//if the proposed point is between mid and bestOuter
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
	 FLOAT_TYPE e=0.0;
	 
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

	xm=(FLOAT_TYPE)0.5*(a+b);       //xm is the midpoint of the bracket (of a and b)
	e=(x>=xm?a-x:b-x);	//set e to the larger of the two bracket intervals
	d=(FLOAT_TYPE)CGOLD*e;

//	assert(a<=x && x<=b);
	 
//	 fw=fv=fx=(*f)(thisnode, thistree, x);
	 
	 for(iter=1;iter<=ITMAX;iter++){

	 	xm=(FLOAT_TYPE)0.5*(a+b);		//xm is the midpoint of the bracket (of a and b)
	 	
	 	tol2=(FLOAT_TYPE)(2.0*(tol1=(FLOAT_TYPE)(tol*fabs(x)+ZEPS)));	
	 	
/*	 	if (fabs(x-xm) <= (tol2-0.5*(b-a))){ //termination condition
	 		*xmin=x;						
	 		return fx;
	 		}
*/// 	if (fabs(e) > tol1){	//construct a trial parabolic fit
	 		r=(x-w)*(fx-fv);
	 		q=(x-v)*(fx-fw);
	 		p=(x-v)*q-(x-w)*r;
	 		q=(FLOAT_TYPE)2.0*(q-r);
	 		if(q>0.0) p=-p;
	 		q=fabs(q);
	 		etemp=e;
	 		e=d;
	 		if(fabs(p) >= fabs(0.5*q*etemp)||p<=q*(a-x) || p>=q*(b-x)){ //determine if the parabolic fit is good
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
	 				FLOAT_TYPE estlnL=(FLOAT_TYPE)(paraMinlnL - paraErr/2.0);
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
		stateProbs[0]=stateProbs[1]=stateProbs[2]=stateProbs[3]=0.0;
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
		stateProbs[0]=stateProbs[1]=stateProbs[2]=stateProbs[3]=0.0;
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

FLOAT_TYPE CalculatePDistance(const char *str1, const char *str2, int nchar){
	FLOAT_TYPE count=0.0;
	int offset1=0, offset2=0;
	int effectiveChar=nchar;
	bool skipChar=false;
	for(int i=0;i<nchar;i++){
		if(str2[i+offset2]<0){
			if(str2[i+offset2]==-4) skipChar=true;
			else{
				int s=-str2[i+offset2];
				for(int i=0;i<s;i++) offset2++;
				}
			}
		if(str1[i+offset1]<0){
			if(str1[i+offset1]==-4) skipChar=true;
			else{
				int s=-str1[i+offset1];
				for(int i=0;i<s;i++) offset1++;
				}
			}
		if(skipChar==false){
			if(str1[i+offset1]!=str2[i+offset2]){
				count += 1.0;
				}
			}
		else effectiveChar--;
		skipChar=false;
		}
	return count/(FLOAT_TYPE)effectiveChar;
	}

#ifndef GANESH
FLOAT_TYPE CalculateHammingDistance(const char *str1, const char *str2, int nchar){
	FLOAT_TYPE count=0.0;
	int offset=0;
	int effectiveChar=nchar;
	for(int i=0;i<nchar;i++){
		if(str2[i+offset]<0){
			if(str2[i+offset]==-4) effectiveChar--;
			else{
				int s=-str2[i+offset];
				for(int i=0;i<s;i++) offset++;
				}
			}
		else if(str1[i]!=str2[i+offset]) count += 1.0;
		}
	return count/(FLOAT_TYPE)effectiveChar;
	}
#else
FLOAT_TYPE CalculateHammingDistance(const char *str1, const char *str2, 
                                const int *col_count, int nchar){
	FLOAT_TYPE count=0.0;
	int offset=0;
	int effectiveChar=0;
	for(int i=0;i<nchar;i++){
        effectiveChar += col_count[i];
		if(str2[i+offset]<0){
			if(str2[i+offset]==-4) {
                effectiveChar = effectiveChar-col_count[i];
            }    
			else{
    				int s=-str2[i+offset];
	    			for(int i=0;i<s;i++) offset++;
				}
			}
		else if(str1[i]!=str2[i+offset]) count += col_count[i]*1.0;
		}
	return count/(FLOAT_TYPE)effectiveChar;
}
#endif

void SampleBranchLengthCurve(FLOAT_TYPE (*func)(TreeNode*, Tree*, FLOAT_TYPE, bool), TreeNode *thisnode, Tree *thistree){
	for(FLOAT_TYPE len=(FLOAT_TYPE)effectiveMin;len<(FLOAT_TYPE)effectiveMax;len*=2.0)
		(*func)(thisnode, thistree, len, true);
	}

void CalcFullCLAInternalInternal(CondLikeArray *destCLA, const CondLikeArray *LCLA, const CondLikeArray *RCLA, const FLOAT_TYPE *Lpr, const FLOAT_TYPE *Rpr, const int nchar, const int nRateCats){
	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.
	FLOAT_TYPE *dest=destCLA->arr;
	const FLOAT_TYPE *LCL=LCLA->arr;
	const FLOAT_TYPE *RCL=RCLA->arr;
	FLOAT_TYPE L1, L2, L3, L4, R1, R2, R3, R4;
	
#ifdef UNIX
	madvise(dest, nchar*4*nRateCats*sizeof(FLOAT_TYPE), MADV_SEQUENTIAL);
	madvise((void *)LCL, nchar*4*nRateCats*sizeof(FLOAT_TYPE), MADV_SEQUENTIAL);
	madvise((void *)RCL, nchar*4*nRateCats*sizeof(FLOAT_TYPE), MADV_SEQUENTIAL);
#endif
	
	if(nRateCats == 4){//the unrolled 4 rate version
#ifdef OMP_INTINTCLA
		#pragma omp parallel for private(dest, LCL, RCL, L1, L2, L3, L4, R1, R2, R3, R4)
		for(int i=0;i<nchar;i++){
			int index=4*4*i;
			dest = &(destCLA->arr[index]);
			LCL = &(LCLA->arr[index]); 
			RCL= &(RCLA->arr[index]);
#else
			for(int i=0;i<nchar;i++){
#endif

			L1=( Lpr[0]*LCL[0]+Lpr[1]*LCL[1]+Lpr[2]*LCL[2]+Lpr[3]*LCL[3]);
			L2=( Lpr[4]*LCL[0]+Lpr[5]*LCL[1]+Lpr[6]*LCL[2]+Lpr[7]*LCL[3]);
			L3=( Lpr[8]*LCL[0]+Lpr[9]*LCL[1]+Lpr[10]*LCL[2]+Lpr[11]*LCL[3]);
			L4=( Lpr[12]*LCL[0]+Lpr[13]*LCL[1]+Lpr[14]*LCL[2]+Lpr[15]*LCL[3]);

			R1=(Rpr[0]*RCL[0]+Rpr[1]*RCL[1]+Rpr[2]*RCL[2]+Rpr[3]*RCL[3]);
			R2=(Rpr[4]*RCL[0]+Rpr[5]*RCL[1]+Rpr[6]*RCL[2]+Rpr[7]*RCL[3]);
			R3=(Rpr[8]*RCL[0]+Rpr[9]*RCL[1]+Rpr[10]*RCL[2]+Rpr[11]*RCL[3]);			
			R4=(Rpr[12]*RCL[0]+Rpr[13]*RCL[1]+Rpr[14]*RCL[2]+Rpr[15]*RCL[3]);
			
			dest[0] = L1 * R1;
			dest[1] = L2 * R2;
			dest[2] = L3 * R3;
			dest[3] = L4 * R4;

			dest+=4;
			LCL+=4;
			RCL+=4;
			
			L1=( Lpr[16+0]*LCL[0]+Lpr[16+1]*LCL[1]+Lpr[16+2]*LCL[2]+Lpr[16+3]*LCL[3]);
			L2=( Lpr[16+4]*LCL[0]+Lpr[16+5]*LCL[1]+Lpr[16+6]*LCL[2]+Lpr[16+7]*LCL[3]);
			L3=( Lpr[16+8]*LCL[0]+Lpr[16+9]*LCL[1]+Lpr[16+10]*LCL[2]+Lpr[16+11]*LCL[3]);
			L4=( Lpr[16+12]*LCL[0]+Lpr[16+13]*LCL[1]+Lpr[16+14]*LCL[2]+Lpr[16+15]*LCL[3]);

			R1=(Rpr[16+0]*RCL[0]+Rpr[16+1]*RCL[1]+Rpr[16+2]*RCL[2]+Rpr[16+3]*RCL[3]);
			R2=(Rpr[16+4]*RCL[0]+Rpr[16+5]*RCL[1]+Rpr[16+6]*RCL[2]+Rpr[16+7]*RCL[3]);
			R3=(Rpr[16+8]*RCL[0]+Rpr[16+9]*RCL[1]+Rpr[16+10]*RCL[2]+Rpr[16+11]*RCL[3]);			
			R4=(Rpr[16+12]*RCL[0]+Rpr[16+13]*RCL[1]+Rpr[16+14]*RCL[2]+Rpr[16+15]*RCL[3]);
			
			dest[0] = L1 * R1;
			dest[1] = L2 * R2;
			dest[2] = L3 * R3;
			dest[3] = L4 * R4;

			dest+=4;
			LCL+=4;
			RCL+=4;		

			L1=( Lpr[32+0]*LCL[0]+Lpr[32+1]*LCL[1]+Lpr[32+2]*LCL[2]+Lpr[32+3]*LCL[3]);
			L2=( Lpr[32+4]*LCL[0]+Lpr[32+5]*LCL[1]+Lpr[32+6]*LCL[2]+Lpr[32+7]*LCL[3]);
			L3=( Lpr[32+8]*LCL[0]+Lpr[32+9]*LCL[1]+Lpr[32+10]*LCL[2]+Lpr[32+11]*LCL[3]);
			L4=( Lpr[32+12]*LCL[0]+Lpr[32+13]*LCL[1]+Lpr[32+14]*LCL[2]+Lpr[32+15]*LCL[3]);

			R1=(Rpr[32+0]*RCL[0]+Rpr[32+1]*RCL[1]+Rpr[32+2]*RCL[2]+Rpr[32+3]*RCL[3]);
			R2=(Rpr[32+4]*RCL[0]+Rpr[32+5]*RCL[1]+Rpr[32+6]*RCL[2]+Rpr[32+7]*RCL[3]);
			R3=(Rpr[32+8]*RCL[0]+Rpr[32+9]*RCL[1]+Rpr[32+10]*RCL[2]+Rpr[32+11]*RCL[3]);			
			R4=(Rpr[32+12]*RCL[0]+Rpr[32+13]*RCL[1]+Rpr[32+14]*RCL[2]+Rpr[32+15]*RCL[3]);
			
			dest[0] = L1 * R1;
			dest[1] = L2 * R2;
			dest[2] = L3 * R3;
			dest[3] = L4 * R4;

			dest+=4;
			LCL+=4;
			RCL+=4;

			L1=( Lpr[48+0]*LCL[0]+Lpr[48+1]*LCL[1]+Lpr[48+2]*LCL[2]+Lpr[48+3]*LCL[3]);
			L2=( Lpr[48+4]*LCL[0]+Lpr[48+5]*LCL[1]+Lpr[48+6]*LCL[2]+Lpr[48+7]*LCL[3]);
			L3=( Lpr[48+8]*LCL[0]+Lpr[48+9]*LCL[1]+Lpr[48+10]*LCL[2]+Lpr[48+11]*LCL[3]);
			L4=( Lpr[48+12]*LCL[0]+Lpr[48+13]*LCL[1]+Lpr[48+14]*LCL[2]+Lpr[48+15]*LCL[3]);

			R1=(Rpr[48+0]*RCL[0]+Rpr[48+1]*RCL[1]+Rpr[48+2]*RCL[2]+Rpr[48+3]*RCL[3]);
			R2=(Rpr[48+4]*RCL[0]+Rpr[48+5]*RCL[1]+Rpr[48+6]*RCL[2]+Rpr[48+7]*RCL[3]);
			R3=(Rpr[48+8]*RCL[0]+Rpr[48+9]*RCL[1]+Rpr[48+10]*RCL[2]+Rpr[48+11]*RCL[3]);			
			R4=(Rpr[48+12]*RCL[0]+Rpr[48+13]*RCL[1]+Rpr[48+14]*RCL[2]+Rpr[48+15]*RCL[3]);
			
			dest[0] = L1 * R1;
			dest[1] = L2 * R2;
			dest[2] = L3 * R3;
			dest[3] = L4 * R4;

			dest+=4;
			LCL+=4;
			RCL+=4;
			}
		}
	
	else{//the general N rate version
		int r;
#ifdef OMP_INTINTCLA
		int index;
		#pragma omp parallel for private(r, index, dest, LCL, RCL, L1, L2, L3, L4, R1, R2, R3, R4)
		for(int i=0;i<nchar;i++) {
			index=4*nRateCats*i;
			dest = &(destCLA->arr[index]);
			LCL = &(LCLA->arr[index]); 
			RCL= &(RCLA->arr[index]);
#else			
		for(int i=0;i<nchar;i++) {
#endif
			for(r=0;r<nRateCats;r++){
				L1=( Lpr[16*r+0]*LCL[0]+Lpr[16*r+1]*LCL[1]+Lpr[16*r+2]*LCL[2]+Lpr[16*r+3]*LCL[3]);
				L2=( Lpr[16*r+4]*LCL[0]+Lpr[16*r+5]*LCL[1]+Lpr[16*r+6]*LCL[2]+Lpr[16*r+7]*LCL[3]);
				L3=( Lpr[16*r+8]*LCL[0]+Lpr[16*r+9]*LCL[1]+Lpr[16*r+10]*LCL[2]+Lpr[16*r+11]*LCL[3]);
				L4=( Lpr[16*r+12]*LCL[0]+Lpr[16*r+13]*LCL[1]+Lpr[16*r+14]*LCL[2]+Lpr[16*r+15]*LCL[3]);

				R1=(Rpr[16*r+0]*RCL[0]+Rpr[16*r+1]*RCL[1]+Rpr[16*r+2]*RCL[2]+Rpr[16*r+3]*RCL[3]);
				R2=(Rpr[16*r+4]*RCL[0]+Rpr[16*r+5]*RCL[1]+Rpr[16*r+6]*RCL[2]+Rpr[16*r+7]*RCL[3]);
				R3=(Rpr[16*r+8]*RCL[0]+Rpr[16*r+9]*RCL[1]+Rpr[16*r+10]*RCL[2]+Rpr[16*r+11]*RCL[3]);			
				R4=(Rpr[16*r+12]*RCL[0]+Rpr[16*r+13]*RCL[1]+Rpr[16*r+14]*RCL[2]+Rpr[16*r+15]*RCL[3]);
					
				dest[0] = L1 * R1;
				dest[1] = L2 * R2;
				dest[2] = L3 * R3;
				dest[3] = L4 * R4;
#ifndef OMP_INTINTCLA
				dest+=4;
				LCL+=4;
				RCL+=4;
#endif
				}
			}
		}

	const int *left_mult=LCLA->underflow_mult;
	const int *right_mult=RCLA->underflow_mult;
	int *undermult=destCLA->underflow_mult;
	
	for(int i=0;i<nchar;i++){
		undermult[i] = left_mult[i] + right_mult[i];
		}
	destCLA->rescaleRank = 2 + LCLA->rescaleRank + RCLA->rescaleRank;
	}

void CalcFullCLATerminalTerminal(CondLikeArray *destCLA, const FLOAT_TYPE *Lpr, const FLOAT_TYPE *Rpr, const char *Ldata, const char *Rdata, const int nchar, const int nRateCats){
	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.
	FLOAT_TYPE *dest=destCLA->arr;
	
#ifdef UNIX
	madvise(dest, nchar*4*nRateCats*sizeof(FLOAT_TYPE), MADV_SEQUENTIAL);
#endif

	for(int i=0;i<nchar;i++){
		if(*Ldata > -1 && *Rdata > -1){
			for(int r=0;r<nRateCats;r++){
				*(dest++) = Lpr[(*Ldata)+16*r] * Rpr[(*Rdata)+16*r];
				*(dest++) = Lpr[(*Ldata+4)+16*r] * Rpr[(*Rdata+4)+16*r];
				*(dest++) = Lpr[(*Ldata+8)+16*r] * Rpr[(*Rdata+8)+16*r];
				*(dest++) = Lpr[(*Ldata+12)+16*r] * Rpr[(*Rdata+12)+16*r];
				}
			Ldata++;
			Rdata++;
			}
			
		else if((*Ldata == -4 && *Rdata == -4) || (*Ldata == -4 && *Rdata > -1) || (*Rdata == -4 && *Ldata > -1)){//total ambiguity of left, right or both
			
			if(*Ldata == -4 && *Rdata == -4) //total ambiguity of both
				for(int i=0;i< (4*nRateCats);i++) *(dest++) = 1.0;
			
			else if(*Ldata == -4){//total ambiguity of left
				for(int i=0;i<nRateCats;i++){
					*(dest++) = Rpr[(*Rdata)+16*i];
					*(dest++) = Rpr[(*Rdata+4)+16*i];
					*(dest++) = Rpr[(*Rdata+8)+16*i];
					*(dest++) = Rpr[(*Rdata+12)+16*i];
					assert(*(dest-4)>=0.0);
					}
				}
			else{//total ambiguity of right
				for(int i=0;i<nRateCats;i++){
					*(dest++) = Lpr[(*Ldata)+16*i];
					*(dest++) = Lpr[(*Ldata+4)+16*i];
					*(dest++) = Lpr[(*Ldata+8)+16*i];
					*(dest++) = Lpr[(*Ldata+12)+16*i];
					assert(*(dest-4)>=0.0);
					}
				}
			Ldata++;
			Rdata++;
			}	
		else {//partial ambiguity of left, right or both
			if(*Ldata>-1){//unambiguous left
				for(int i=0;i<nRateCats;i++){
					*(dest+(i*4)) = Lpr[(*Ldata)+16*i];
					*(dest+(i*4)+1) = Lpr[(*Ldata+4)+16*i];
					*(dest+(i*4)+2) = Lpr[(*Ldata+8)+16*i];
					*(dest+(i*4)+3) = Lpr[(*Ldata+12)+16*i];
					assert(*(dest)>=0.0);
					}
				Ldata++;
				}
			else{
				if(*Ldata==-4){//fully ambiguous left
					for(int i=0;i< (4*nRateCats);i++){
						*(dest+i)=1.0;
						}
					Ldata++;
					}
			 
				else{//partially ambiguous left
					int nstates=-*(Ldata++);
					for(int q=0;q< (4*nRateCats);q++) dest[q]=0;
					for(int i=0;i<nstates;i++){
						for(int r=0;r<nRateCats;r++){
							*(dest+(r*4)) += Lpr[(*Ldata)+16*r];
							*(dest+(r*4)+1) += Lpr[(*Ldata+4)+16*r];
							*(dest+(r*4)+2) += Lpr[(*Ldata+8)+16*r];
							*(dest+(r*4)+3) += Lpr[(*Ldata+12)+16*r];
//							assert(*(dest-1)>0.0);
							}
						Ldata++;
						}
					}
				}
			if(*Rdata>-1){//unambiguous right
				for(int i=0;i<nRateCats;i++){
					*(dest++) *= Rpr[(*Rdata)+16*i];
					*(dest++) *= Rpr[(*Rdata+4)+16*i];
					*(dest++) *= Rpr[(*Rdata+8)+16*i];
					*(dest++) *= Rpr[(*Rdata+12)+16*i];
//					assert(*(dest-1)>0.0);
					}
				Rdata++;		
				}
			else if(*Rdata != -4){//partially ambiguous right
				char nstates=-1 * *(Rdata++);
				//create a temporary cla to hold the results from the ambiguity of the right, 
				//which need to be +'s 
				//FLOAT_TYPE *tempcla=new FLOAT_TYPE[4*nRateCats];
				vector<FLOAT_TYPE> tempcla(4*nRateCats);
				for(int i=0;i<nstates;i++){
					for(int r=0;r<nRateCats;r++){
						tempcla[(r*4)]   += Rpr[(*Rdata)+16*r];
						tempcla[(r*4)+1] += Rpr[(*Rdata+4)+16*r];
						tempcla[(r*4)+2] += Rpr[(*Rdata+8)+16*r];
						tempcla[(r*4)+3] += Rpr[(*Rdata+12)+16*r];
//						assert(*(dest-1)>0.0);
						}
					Rdata++;
					}
				//Now multiply the temporary results against the already calced left
				for(int i=0;i<nRateCats;i++){
					*(dest++) *= tempcla[(i*4)];
					*(dest++) *= tempcla[(i*4)+1];
					*(dest++) *= tempcla[(i*4)+2];
					*(dest++) *= tempcla[(i*4)+3];
//					assert(*(dest-1)>0.0);
					}
				}
			else{//fully ambiguous right
				dest+=(4*nRateCats);
				Rdata++;
				}
			}
		}
		
		for(int site=0;site<nchar;site++){
			destCLA->underflow_mult[site]=0;
			}
		destCLA->rescaleRank=2;
	}

void CalcFullCLAInternalTerminal(CondLikeArray *destCLA, const CondLikeArray *LCLA, const FLOAT_TYPE *pr1, const FLOAT_TYPE *pr2, char *dat2, const int nchar, const int nRateCats, const unsigned *ambigMap /*=NULL*/){
	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.
	FLOAT_TYPE *des=destCLA->arr;
	FLOAT_TYPE *dest=des;
	const FLOAT_TYPE *CL=LCLA->arr;
	const FLOAT_TYPE *CL1=CL;
	const char *data2=dat2;

#ifdef UNIX	
	madvise(dest, nchar*4*nRateCats*sizeof(FLOAT_TYPE), MADV_SEQUENTIAL);
	madvise((void*)CL1, nchar*4*nRateCats*sizeof(FLOAT_TYPE), MADV_SEQUENTIAL);	
#endif

	if(nRateCats==4){//unrolled 4 rate version
#ifdef OMP_INTTERMCLA
		#pragma omp parallel for private(dest, CL1, data2)
		for(int i=0;i<nchar;i++){
			dest=&des[4*4*i];
			CL1=&CL[4*4*i];
			data2=&dat2[ambigMap[i]];
#else
		for(int i=0;i<nchar;i++){
#endif
			if(*data2 > -1){ //no ambiguity
				dest[0] = ( pr1[0]*CL1[0]+pr1[1]*CL1[1]+pr1[2]*CL1[2]+pr1[3]*CL1[3]) * pr2[(*data2)];
				dest[1] = ( pr1[4]*CL1[0]+pr1[5]*CL1[1]+pr1[6]*CL1[2]+pr1[7]*CL1[3]) * pr2[(*data2+4)];
				dest[2] = ( pr1[8]*CL1[0]+pr1[9]*CL1[1]+pr1[10]*CL1[2]+pr1[11]*CL1[3]) * pr2[(*data2+8)];
				dest[3] = ( pr1[12]*CL1[0]+pr1[13]*CL1[1]+pr1[14]*CL1[2]+pr1[15]*CL1[3]) * pr2[(*data2+12)];
				
				dest[4] = ( pr1[16]*CL1[4]+pr1[17]*CL1[5]+pr1[18]*CL1[6]+pr1[19]*CL1[7]) * pr2[(*data2)+16];
				dest[5] = ( pr1[20]*CL1[4]+pr1[21]*CL1[5]+pr1[22]*CL1[6]+pr1[23]*CL1[7]) * pr2[(*data2+4)+16];
				dest[6] = ( pr1[24]*CL1[4]+pr1[25]*CL1[5]+pr1[26]*CL1[6]+pr1[27]*CL1[7]) * pr2[(*data2+8)+16];
				dest[7] = ( pr1[28]*CL1[4]+pr1[29]*CL1[5]+pr1[30]*CL1[6]+pr1[31]*CL1[7]) * pr2[(*data2+12)+16];
			
				dest[8] = ( pr1[32]*CL1[8]+pr1[33]*CL1[9]+pr1[34]*CL1[10]+pr1[35]*CL1[11]) * pr2[(*data2)+32];
				dest[9] = ( pr1[36]*CL1[8]+pr1[37]*CL1[9]+pr1[38]*CL1[10]+pr1[39]*CL1[11]) * pr2[(*data2+4)+32];
				dest[10] = ( pr1[40]*CL1[8]+pr1[41]*CL1[9]+pr1[42]*CL1[10]+pr1[43]*CL1[11]) * pr2[(*data2+8)+32];
				dest[11] = ( pr1[44]*CL1[8]+pr1[45]*CL1[9]+pr1[46]*CL1[10]+pr1[47]*CL1[11]) * pr2[(*data2+12)+32];
			
				dest[12] = ( pr1[48]*CL1[12]+pr1[49]*CL1[13]+pr1[50]*CL1[14]+pr1[51]*CL1[15]) * pr2[(*data2)+48];
				dest[13] = ( pr1[52]*CL1[12]+pr1[53]*CL1[13]+pr1[54]*CL1[14]+pr1[55]*CL1[15]) * pr2[(*data2+4)+48];
				dest[14] = ( pr1[56]*CL1[12]+pr1[57]*CL1[13]+pr1[58]*CL1[14]+pr1[59]*CL1[15]) * pr2[(*data2+8)+48];
				dest[15] = ( pr1[60]*CL1[12]+pr1[61]*CL1[13]+pr1[62]*CL1[14]+pr1[63]*CL1[15]) * pr2[(*data2+12)+48];
#ifndef OMP_INTTERMCLA
				dest+=16;
				data2++;
#endif
				}
			else if(*data2 == -4){//total ambiguity
				dest[0] = ( pr1[0]*CL1[0]+pr1[1]*CL1[1]+pr1[2]*CL1[2]+pr1[3]*CL1[3]);
				dest[1] = ( pr1[4]*CL1[0]+pr1[5]*CL1[1]+pr1[6]*CL1[2]+pr1[7]*CL1[3]);
				dest[2] = ( pr1[8]*CL1[0]+pr1[9]*CL1[1]+pr1[10]*CL1[2]+pr1[11]*CL1[3]);
				dest[3] = ( pr1[12]*CL1[0]+pr1[13]*CL1[1]+pr1[14]*CL1[2]+pr1[15]*CL1[3]);
				
				dest[4] = ( pr1[16]*CL1[4]+pr1[17]*CL1[5]+pr1[18]*CL1[6]+pr1[19]*CL1[7]);
				dest[5] = ( pr1[20]*CL1[4]+pr1[21]*CL1[5]+pr1[22]*CL1[6]+pr1[23]*CL1[7]);
				dest[6] = ( pr1[24]*CL1[4]+pr1[25]*CL1[5]+pr1[26]*CL1[6]+pr1[27]*CL1[7]);
				dest[7] = ( pr1[28]*CL1[4]+pr1[29]*CL1[5]+pr1[30]*CL1[6]+pr1[31]*CL1[7]);
			
				dest[8] = ( pr1[32]*CL1[8]+pr1[33]*CL1[9]+pr1[34]*CL1[10]+pr1[35]*CL1[11]);
				dest[9] = ( pr1[36]*CL1[8]+pr1[37]*CL1[9]+pr1[38]*CL1[10]+pr1[39]*CL1[11]);
				dest[10] = ( pr1[40]*CL1[8]+pr1[41]*CL1[9]+pr1[42]*CL1[10]+pr1[43]*CL1[11]);
				dest[11] = ( pr1[44]*CL1[8]+pr1[45]*CL1[9]+pr1[46]*CL1[10]+pr1[47]*CL1[11]);
			
				dest[12] = ( pr1[48]*CL1[12]+pr1[49]*CL1[13]+pr1[50]*CL1[14]+pr1[51]*CL1[15]);
				dest[13] = ( pr1[52]*CL1[12]+pr1[53]*CL1[13]+pr1[54]*CL1[14]+pr1[55]*CL1[15]);
				dest[14] = ( pr1[56]*CL1[12]+pr1[57]*CL1[13]+pr1[58]*CL1[14]+pr1[59]*CL1[15]);
				dest[15] = ( pr1[60]*CL1[12]+pr1[61]*CL1[13]+pr1[62]*CL1[14]+pr1[63]*CL1[15]);
#ifndef OMP_INTTERMCLA
				dest+=16;
				data2++;
#endif
				}
			else {//partial ambiguity
				//first figure in the ambiguous terminal
				int nstates=-1 * *(data2++);
				for(int j=0;j<16;j++) dest[j]=0.0;
				for(int s=0;s<nstates;s++){
					for(int r=0;r<4;r++){
						*(dest+(r*4)) += pr2[(*data2)+16*r];
						*(dest+(r*4)+1) += pr2[(*data2+4)+16*r];
						*(dest+(r*4)+2) += pr2[(*data2+8)+16*r];
						*(dest+(r*4)+3) += pr2[(*data2+12)+16*r];
						}
					data2++;
					}
				
				//now add the internal child
				*(dest++) *= ( pr1[0]*CL1[0]+pr1[1]*CL1[1]+pr1[2]*CL1[2]+pr1[3]*CL1[3]);
				*(dest++) *= ( pr1[4]*CL1[0]+pr1[5]*CL1[1]+pr1[6]*CL1[2]+pr1[7]*CL1[3]);
				*(dest++) *= ( pr1[8]*CL1[0]+pr1[9]*CL1[1]+pr1[10]*CL1[2]+pr1[11]*CL1[3]);
				*(dest++) *= ( pr1[12]*CL1[0]+pr1[13]*CL1[1]+pr1[14]*CL1[2]+pr1[15]*CL1[3]);
				
				*(dest++) *= ( pr1[16]*CL1[4]+pr1[17]*CL1[5]+pr1[18]*CL1[6]+pr1[19]*CL1[7]);
				*(dest++) *= ( pr1[20]*CL1[4]+pr1[21]*CL1[5]+pr1[22]*CL1[6]+pr1[23]*CL1[7]);
				*(dest++) *= ( pr1[24]*CL1[4]+pr1[25]*CL1[5]+pr1[26]*CL1[6]+pr1[27]*CL1[7]);
				*(dest++) *= ( pr1[28]*CL1[4]+pr1[29]*CL1[5]+pr1[30]*CL1[6]+pr1[31]*CL1[7]);
			
				*(dest++) *= ( pr1[32]*CL1[8]+pr1[33]*CL1[9]+pr1[34]*CL1[10]+pr1[35]*CL1[11]);
				*(dest++) *= ( pr1[36]*CL1[8]+pr1[37]*CL1[9]+pr1[38]*CL1[10]+pr1[39]*CL1[11]);
				*(dest++) *= ( pr1[40]*CL1[8]+pr1[41]*CL1[9]+pr1[42]*CL1[10]+pr1[43]*CL1[11]);
				*(dest++) *= ( pr1[44]*CL1[8]+pr1[45]*CL1[9]+pr1[46]*CL1[10]+pr1[47]*CL1[11]);
			
				*(dest++) *= ( pr1[48]*CL1[12]+pr1[49]*CL1[13]+pr1[50]*CL1[14]+pr1[51]*CL1[15]);
				*(dest++) *= ( pr1[52]*CL1[12]+pr1[53]*CL1[13]+pr1[54]*CL1[14]+pr1[55]*CL1[15]);
				*(dest++) *= ( pr1[56]*CL1[12]+pr1[57]*CL1[13]+pr1[58]*CL1[14]+pr1[59]*CL1[15]);
				*(dest++) *= ( pr1[60]*CL1[12]+pr1[61]*CL1[13]+pr1[62]*CL1[14]+pr1[63]*CL1[15]);
				}
#ifndef OMP_INTTERMCLA		
			CL1+=16;
#endif
			}
		}
	else{//general N rate version
#ifdef OMP_INTTERMCLA
		#pragma omp parallel for private(dest, CL1, data2)
		for(int i=0;i<nchar;i++){
			dest=&des[4*nRateCats*i];
			CL1=&CL[4*nRateCats*i];
			data2=&dat2[ambigMap[i]];
#else
		for(int i=0;i<nchar;i++){
#endif
			if(*data2 > -1){ //no ambiguity
				for(int r=0;r<nRateCats;r++){
					dest[0] = ( pr1[16*r+0]*CL1[4*r+0]+pr1[16*r+1]*CL1[4*r+1]+pr1[16*r+2]*CL1[4*r+2]+pr1[16*r+3]*CL1[4*r+3]) * pr2[(*data2)+16*r];
					dest[1] = ( pr1[16*r+4]*CL1[4*r+0]+pr1[16*r+5]*CL1[4*r+1]+pr1[16*r+6]*CL1[4*r+2]+pr1[16*r+7]*CL1[4*r+3]) * pr2[(*data2+4)+16*r];
					dest[2] = ( pr1[16*r+8]*CL1[4*r+0]+pr1[16*r+9]*CL1[4*r+1]+pr1[16*r+10]*CL1[4*r+2]+pr1[16*r+11]*CL1[4*r+3]) * pr2[(*data2+8)+16*r];
					dest[3] = ( pr1[16*r+12]*CL1[4*r+0]+pr1[16*r+13]*CL1[4*r+1]+pr1[16*r+14]*CL1[4*r+2]+pr1[16*r+15]*CL1[4*r+3]) * pr2[(*data2+12)+16*r];
					dest+=4;
					}
				data2++;
				}
			else if(*data2 == -4){//total ambiguity
				for(int r=0;r<nRateCats;r++){
					dest[0] = ( pr1[16*r+0]*CL1[4*r+0]+pr1[16*r+1]*CL1[4*r+1]+pr1[16*r+2]*CL1[4*r+2]+pr1[16*r+3]*CL1[4*r+3]);
					dest[1] = ( pr1[16*r+4]*CL1[4*r+0]+pr1[16*r+5]*CL1[4*r+1]+pr1[16*r+6]*CL1[4*r+2]+pr1[16*r+7]*CL1[4*r+3]);
					dest[2] = ( pr1[16*r+8]*CL1[4*r+0]+pr1[16*r+9]*CL1[4*r+1]+pr1[16*r+10]*CL1[4*r+2]+pr1[16*r+11]*CL1[4*r+3]);
					dest[3] = ( pr1[16*r+12]*CL1[4*r+0]+pr1[16*r+13]*CL1[4*r+1]+pr1[16*r+14]*CL1[4*r+2]+pr1[16*r+15]*CL1[4*r+3]);
					dest+=4;
					}
				data2++;
				}
			else {//partial ambiguity
				//first figure in the ambiguous terminal
				int nstates=-1 * *(data2++);
				for(int q=0;q<4*nRateCats;q++) dest[q]=0;
				for(int s=0;s<nstates;s++){
					for(int r=0;r<nRateCats;r++){
						*(dest+(r*4)) += pr2[(*data2)+16*r];
						*(dest+(r*4)+1) += pr2[(*data2+4)+16*r];
						*(dest+(r*4)+2) += pr2[(*data2+8)+16*r];
						*(dest+(r*4)+3) += pr2[(*data2+12)+16*r];
						}
					data2++;
					}
				
				//now add the internal child
				for(int r=0;r<nRateCats;r++){
					*(dest++) *= ( pr1[16*r+0]*CL1[4*r+0]+pr1[16*r+1]*CL1[4*r+1]+pr1[16*r+2]*CL1[4*r+2]+pr1[16*r+3]*CL1[4*r+3]);
					*(dest++) *= ( pr1[16*r+4]*CL1[4*r+0]+pr1[16*r+5]*CL1[4*r+1]+pr1[16*r+6]*CL1[4*r+2]+pr1[16*r+7]*CL1[4*r+3]);
					*(dest++) *= ( pr1[16*r+8]*CL1[4*r+0]+pr1[16*r+9]*CL1[4*r+1]+pr1[16*r+10]*CL1[4*r+2]+pr1[16*r+11]*CL1[4*r+3]);
					*(dest++) *= ( pr1[16*r+12]*CL1[4*r+0]+pr1[16*r+13]*CL1[4*r+1]+pr1[16*r+14]*CL1[4*r+2]+pr1[16*r+15]*CL1[4*r+3]);
					}
				}
#ifndef OMP_INTTERMCLA	
			CL1 += 4*nRateCats;
#endif
			}
		}
		
	for(int i=0;i<nchar;i++)
		destCLA->underflow_mult[i]=LCLA->underflow_mult[i];
	
	destCLA->rescaleRank=LCLA->rescaleRank+2;
	} 

void CalcFullCLAPartialInternalRateHet(CondLikeArray *destCLA, const CondLikeArray *LCLA, const FLOAT_TYPE *pr1, CondLikeArray *partialCLA, int nchar, int nRateCats /*=4*/){
	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.
	FLOAT_TYPE *dest=destCLA->arr;
	FLOAT_TYPE *CL1=LCLA->arr;
	FLOAT_TYPE *partial=partialCLA->arr;
	
#ifdef UNIX
	madvise(dest, nchar*4*nRateCats*sizeof(FLOAT_TYPE), MADV_SEQUENTIAL);
	madvise((void*)CL1, nchar*4*nRateCats*sizeof(FLOAT_TYPE), MADV_SEQUENTIAL);
	madvise(partial, nchar*4*nRateCats*sizeof(FLOAT_TYPE), MADV_SEQUENTIAL);
#endif

	if(nRateCats==4){
		for(int i=0;i<nchar;i++){
			*(dest++) = ( pr1[0]*CL1[0]+pr1[1]*CL1[1]+pr1[2]*CL1[2]+pr1[3]*CL1[3]) * *(partial++);
			*(dest++) = ( pr1[4]*CL1[0]+pr1[5]*CL1[1]+pr1[6]*CL1[2]+pr1[7]*CL1[3]) * *(partial++);
			*(dest++) = ( pr1[8]*CL1[0]+pr1[9]*CL1[1]+pr1[10]*CL1[2]+pr1[11]*CL1[3]) * *(partial++);
			*(dest++) = ( pr1[12]*CL1[0]+pr1[13]*CL1[1]+pr1[14]*CL1[2]+pr1[15]*CL1[3]) * *(partial++);
			
			*(dest++) = ( pr1[16]*CL1[4]+pr1[17]*CL1[5]+pr1[18]*CL1[6]+pr1[19]*CL1[7]) * *(partial++);
			*(dest++) = ( pr1[20]*CL1[4]+pr1[21]*CL1[5]+pr1[22]*CL1[6]+pr1[23]*CL1[7]) * *(partial++);
			*(dest++) = ( pr1[24]*CL1[4]+pr1[25]*CL1[5]+pr1[26]*CL1[6]+pr1[27]*CL1[7]) * *(partial++);
			*(dest++) = ( pr1[28]*CL1[4]+pr1[29]*CL1[5]+pr1[30]*CL1[6]+pr1[31]*CL1[7]) * *(partial++);
		
			*(dest++) = ( pr1[32]*CL1[8]+pr1[33]*CL1[9]+pr1[34]*CL1[10]+pr1[35]*CL1[11]) * *(partial++);
			*(dest++) = ( pr1[36]*CL1[8]+pr1[37]*CL1[9]+pr1[38]*CL1[10]+pr1[39]*CL1[11]) * *(partial++);
			*(dest++) = ( pr1[40]*CL1[8]+pr1[41]*CL1[9]+pr1[42]*CL1[10]+pr1[43]*CL1[11]) * *(partial++);
			*(dest++) = ( pr1[44]*CL1[8]+pr1[45]*CL1[9]+pr1[46]*CL1[10]+pr1[47]*CL1[11]) * *(partial++);
		
			*(dest++) = ( pr1[48]*CL1[12]+pr1[49]*CL1[13]+pr1[50]*CL1[14]+pr1[51]*CL1[15]) * *(partial++);
			*(dest++) = ( pr1[52]*CL1[12]+pr1[53]*CL1[13]+pr1[54]*CL1[14]+pr1[55]*CL1[15]) * *(partial++);
			*(dest++) = ( pr1[56]*CL1[12]+pr1[57]*CL1[13]+pr1[58]*CL1[14]+pr1[59]*CL1[15]) * *(partial++);
			*(dest++) = ( pr1[60]*CL1[12]+pr1[61]*CL1[13]+pr1[62]*CL1[14]+pr1[63]*CL1[15]) * *(partial++);
			CL1+=16;
			assert(*(dest-1)>0.0);
			}
		}
	else{
		for(int i=0;i<nchar;i++){
			for(int r=0;r<nRateCats;r++){
				*(dest++) = ( pr1[16*r+0]*CL1[4*r+0]+pr1[16*r+1]*CL1[4*r+1]+pr1[16*r+2]*CL1[4*r+2]+pr1[16*r+3]*CL1[4*r+3]) * *(partial++);
				*(dest++) = ( pr1[16*r+4]*CL1[4*r+0]+pr1[16*r+5]*CL1[4*r+1]+pr1[16*r+6]*CL1[4*r+2]+pr1[16*r+7]*CL1[4*r+3]) * *(partial++);
				*(dest++) = ( pr1[16*r+8]*CL1[4*r+0]+pr1[16*r+9]*CL1[4*r+1]+pr1[16*r+10]*CL1[4*r+2]+pr1[16*r+11]*CL1[4*r+3]) * *(partial++);
				*(dest++) = ( pr1[16*r+12]*CL1[4*r+0]+pr1[16*r+13]*CL1[4*r+1]+pr1[16*r+14]*CL1[4*r+2]+pr1[16*r+15]*CL1[4*r+3]) * *(partial++);
				CL1+=4;
				assert(*(dest-1)>0.0);
				}
			}
		}
		
	for(int site=0;site<nchar;site++){
		destCLA->underflow_mult[site]=partialCLA->underflow_mult[site] + LCLA->underflow_mult[site];
		}
	}

void CalcFullCLAPartialTerminalRateHet(CondLikeArray *destCLA, const CondLikeArray *partialCLA, const FLOAT_TYPE *Lpr, char *Ldata, int nchar, int nRateCats /*=4*/){
	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.
	FLOAT_TYPE *dest=destCLA->arr;
	FLOAT_TYPE *partial=partialCLA->arr;
#ifdef UNIX
	madvise(dest, nchar*4*nRateCats*sizeof(FLOAT_TYPE), MADV_SEQUENTIAL);
	madvise((void*)partial, nchar*4*nRateCats*sizeof(FLOAT_TYPE), MADV_SEQUENTIAL);
#endif

	for(int i=0;i<nchar;i++){
		if(*Ldata > -1){ //no ambiguity
			for(int i=0;i<nRateCats;i++){
				*(dest++) = Lpr[(*Ldata)+16*i] * *(partial++);
				*(dest++) = Lpr[(*Ldata+4)+16*i] * *(partial++);
				*(dest++) = Lpr[(*Ldata+8)+16*i] * *(partial++);
				*(dest++) = Lpr[(*Ldata+12)+16*i] * *(partial++);
//				assert(*(dest-1)>0.0);
				}
			Ldata++;
			}
			
		else if(*Ldata == -4){ //total ambiguity
			for(int i=0;i<4*nRateCats;i++) *(dest++) = *(partial++);
			Ldata++;
			}
		else{ //partial ambiguity
			//first figure in the ambiguous terminal
			char nstates=-1 * *(Ldata++);
			for(int q=0;q<4*nRateCats;q++) dest[q]=0;
			for(int i=0;i<nstates;i++){
				for(int i=0;i<nRateCats;i++){
					*(dest+(i*4)) += Lpr[(*Ldata)+16*i];
					*(dest+(i*4)+1) += Lpr[(*Ldata+4)+16*i];
					*(dest+(i*4)+2) += Lpr[(*Ldata+8)+16*i];
					*(dest+(i*4)+3) += Lpr[(*Ldata+12)+16*i];
//					assert(*(dest-1)>0.0);
					}
				Ldata++;
				}
			
			//now add the partial
			for(int r=0;r<nRateCats;r++){
				*(dest++) *= *(partial++);
				*(dest++) *= *(partial++);
				*(dest++) *= *(partial++);
				*(dest++) *= *(partial++);
				}
			}
		}
	for(int i=0;i<nchar;i++)
		destCLA->underflow_mult[i]=partialCLA->underflow_mult[i];
}
