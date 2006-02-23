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

//	NOTE: Portions of this source adapted from:
//	Press, W. H., B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling.  1992. 
//	Numerical Recipes in C : The Art of Scientific Computing.  Cambridge University Press, Cambridge.

#include <unistd.h>

#include "funcs.h"
#include "population.h"
#include "parameters.h"
#include "tree.h"
#include "defs.h"

#undef ROOT_OPT
#define FOURTH_ROOT

#define LOG_MIN_BRLEN log(DEF_MIN_BRLEN)

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
	if( !FileExists( params.statefname ) )	{
	cout << "Error opening state file: " << params.statefname << endl;
		exit(0);
	}

	ifstream sf( params.statefname );
	sf >> params.prev_generations >> params.prev_time >> params.randomSeed;
	sf.close();
	
	rnd.set_seed( params.randomSeed );
}

int ReadData(const char* filename, HKYData* data)	{
	if (!FileExists(filename))	{
		cout << "ERROR - data file not found: " << filename << endl;
		return -1;
	}
#ifndef UD_VERSION
	cout  << "Reading data file: " << filename << "..." << endl;
#endif
	data->Read( filename );
#ifndef UD_VERSION
	if( data->AmbiguousStates() ) {
		cout << "  At least one IUPAC/IUB ambiguity code found." << endl;
		cout << "  Ambiguity codes are treated as missing data." << endl;
	}
#endif
	// report summary statistics about data
	data->Summarize();
#ifndef UD_VERSION
	cout << "  " << data->NConstant() << " constant characters." << endl;
	cout << "  " << data->NInformative() << " parsimony-informative characters." << endl;
	cout << "  " << data->NAutapomorphic() << " autapomorphic characters." << endl;
	int total = data->NConstant() + data->NInformative() + data->NAutapomorphic();
	cout << "  " << total << " total characters." << endl;
#endif
	// try to compress
	if (!data->Dense())	{
	#ifndef UD_VERSION
		cout << "Compressing data file..." << endl;
	#endif
		data->Collapse();
	#ifndef UD_VERSION
		cout << "  " << data->NChar() << " columns in data matrix after compression." << endl;
	#endif
	}
	else {
		#ifndef UD_VERSION
		cout << endl << "Datafile already compressed."  << endl;
		cout << data->NChar() << " columns in compressed data matrix." << endl;
		#endif
	}
	data->DetermineConstantSites();
	if(!data->Dense()) data->Save(filename, "new");
	return 0;
	}

int ReadData(const Parameters& params, HKYData* data)	{

	// regurgitate params specified
	if( params.restart ) {
		cout << "  Restarting using state file \"" << params.statefname << "\"" << endl;
		GetRestartParams( const_cast<Parameters&>(params) );
		cout << "    random number seed set to " << params.randomSeed << endl;
		cout << "    last generation from previous run was " << params.prev_generations << endl;
		cout << "    starting with previous elapsed time, which was " << params.prev_time << " seconds" << endl;
	}

	const_cast<Parameters&>(params).BriefReport( cout );
	cout << endl;

	// Check to be sure data file exists
	//
	if( !FileExists( params.datafname ) )	{
		cout << "ERROR - data file does not exist: " << params.datafname << endl;
		return -1;
	}

	// Read in the data matrix
	//
	cout.flush();
	cout << endl << "Reading data file " << params.datafname << "..." << endl;
	data->Read( params.datafname );
	if( data->AmbiguousStates() ) {
		cout << "  At least one IUPAC/IUB ambiguity code found" << endl;
		cout << "    Ambiguity codes are treated as missing data." << endl;
	}

	// report summary statistics about data
	data->Summarize();
	cout << "  " << data->NConstant() << " constant characters" << endl;
	cout << "  " << data->NInformative() << " parsimony-informative characters" << endl;
	cout << "  " << data->NAutapomorphic() << " autapomorphic characters" << endl;
	int total = data->NConstant() + data->NInformative() + data->NAutapomorphic();
	cout << "  " << total << " total characters" << endl;
	cout.flush();

	//DZ Only compress and write data to file if dense=0 (data is not already compressed)
	if(!(data->Dense())){
		cout << endl << "Compressing data file..." << endl;
		data->Collapse();
		data->Save("compdata.nex", "new");
		cout << "  " << data->NChar() << " columns in data matrix after compression" << endl;
	}
	else {
		cout << endl << "Datafile already compressed."  << endl;
		cout << data->NChar() << " columns in compressed data matrix" << endl;
		}

	data->DetermineConstantSites();
	cout << endl;
	cout.flush();
	return 0;

}

int RandomInt(int lb, int ub)	{
	return lb + rand() % (ub-lb+1);
}

double RandomFrac()	{
	return (double)rand() / RAND_MAX;
}

double RandomDouble(double lb, double ub)	{
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
int mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double (*func)(TreeNode*, Tree*, double), TreeNode *thisnode, Tree *thistree){
	double ulim, u, r, q, fu;
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
	*cx = (*cx > DEF_MIN_BRLEN ? (*cx < DEF_MAX_BRLEN ? *cx : DEF_MAX_BRLEN) : DEF_MIN_BRLEN);
	*fc=(*func)(thisnode, thistree, *cx);
*/
	while(*fb>*fc){
			
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);			
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/(2.0*SIGN(FMAX(fabs(q-r),TINY), q-r));
		u = (u > DEF_MIN_BRLEN ? (u < DEF_MAX_BRLEN ? u : DEF_MAX_BRLEN) : DEF_MIN_BRLEN);
		ulim=(*bx)+GLIMIT*(*cx-*bx);
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
			u=(*cx)+GOLD*(*cx-*bx);
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
				SHFT(*bx, *cx, u, *cx+GOLD*(*cx-*bx));
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
			u=(*cx)+GOLD*(*cx-*bx);
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
double brent(double ax, double bx, double cx, double (*f)(TreeNode *, Tree*, double), double tol, double *xmin, TreeNode *thisnode, Tree *thistree){
	 int iter;
	 double a, b, d, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
	 double e=0.0;
	 
	 a=(ax < cx ? ax : cx); //make a the smallest of the three bracket points 
	 b=(ax > cx ? ax : cx); //and b the largest
	 x=w=v=bx;				//make x the current minimum, as well as w and v
	 
	 fw=fv=fx=(*f)(thisnode, thistree, x);
	 
	 for(iter=1;iter<=ITMAX;iter++){
	 	xm=0.5*(a+b);		//xm is the midpoint of the bracket (of a and b)
	 	
	 	tol2=2.0*(tol1=tol*fabs(x)+ZEPS);	
	 	
	 	if (fabs(x-xm) <= (tol2-0.5*(b-a))){ //termination condition
	 		*xmin=x;							//if the distance between x and bracket mean is < 
	 		return fx;
	 		}
	 	if (fabs(e) > tol1){	//construct a trial parabolic fit
	 		r=(x-w)*(fx-fv);
	 		q=(x-v)*(fx-fw);
	 		p=(x-v)*q-(x-w)*r;
	 		q=2.0*(q-r);
	 		if(q>0.0) p=-p;
	 		q=fabs(q);
	 		etemp=e;
	 		e=d;
	 		if(fabs(p) >= fabs(0.5*q*etemp)||p<=q*(a-x) || p>=q*(b-x)) //determine if the parabolic fit is good
	 			d=CGOLD*(e=(x>=xm?a-x:b-x));  //if not
	 			
	 		else{				//if so, take the parabolic step
	 			d=p/q;
	 			u=x+d;
	 			if(u-a < tol2||b-u<tol2)
	 				d=SIGN(tol1,xm-x);
	 			}
	 		}
	 	else{
	 		d=CGOLD*(e=(x>=xm?a-x:b-x)); //e is the distance moved in the step before last
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
#define effectiveMin DEF_MIN_BRLEN
#define effectiveMax DEF_MAX_BRLEN
#define sweetspot 0.00000256
#define smallShift 0.00000016
#endif


//My version of the bracketing fuction that can abort and stop evaluating under certain conditions
//it assumes that the raw branch lengths are being passed in (not logs) so values are bounded
//by the minimum branch length
int DZbrak(double *worstOuter, double *mid, double *bestOuter, double *worstOuterL, double *midL, double *bestOuterL, double (*func)(TreeNode*, Tree*, double, bool), TreeNode *thisnode, Tree *thistree){
	//points are always passed in such that worstOuter < mid < bestOuter
	double nextTry, r, q, nextTryL, dum;
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
			*mid=(*worstOuter+*bestOuter)*0.5;
			*midL=(*func)(thisnode, thistree, *mid, true);
			if(*bestOuterL < *midL && !(*mid > sweetspot)) return 1;
			}
		else if(!(*worstOuter > sweetspot)){
			SHFT(dum, *mid, *worstOuter, dum)
			SHFT(dum, *midL, *worstOuterL, dum)
			*bestOuter=effectiveMin;
			*bestOuterL=(*func)(thisnode, thistree, *bestOuter, true);
			if(*bestOuterL < *midL && !(*mid > sweetspot)) return 1;
			}
		else{
			*bestOuter=sweetspot-.02;
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
				nextTry=0.16;
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
			double diffMidBestL=(*midL-*bestOuterL);
			double diffMidBest=(*mid-*bestOuter);
			double diffMidWorstL=(*worstOuterL-*midL);
			double diffMidWorst=(*worstOuter-*mid);
			if(*worstOuter < *bestOuter){ //case 2
				//check the curvature
		
				double slopeRatio=(diffMidBestL/diffMidBest) / (diffMidWorstL/diffMidWorst);
				if(slopeRatio > 0.9){
					nextTry=((*bestOuter)+2.0*(*bestOuter-*mid));//case 2a
					}
				else{ //case 2b and 2c
					r=diffMidWorst*diffMidBestL;
					q=diffMidBest*diffMidWorstL;
					nextTry=(*mid)-(diffMidBest*q-(*mid-*worstOuter)*r)/(2.0*SIGN(FMAX(fabs(q-r),TINY), q-r));
					if(/*nextTry > *bestOuter && */fabs(nextTry-*bestOuter) < smallShift){
						//if the parabolic estimate is very near our current best it tends to take
						//a while to get the bracket, so just push it a little further to the right
						nextTry += smallShift;
						}
					}
				}
			
			else if(*worstOuter > *bestOuter){ //case 3
				r=diffMidWorst*diffMidBestL;
				q=diffMidBest*diffMidWorstL;
				nextTry=(*mid)-(diffMidBest*q-diffMidWorst*r)/(2.0*SIGN(FMAX(fabs(q-r),TINY), q-r));
				if(*bestOuter==effectiveMin){
					if(nextTry < effectiveMin || *mid < sweetspot || possibleZeroMLE) return 1; //case 3a2
					else {//case 3a1
						//just go with the parabolic
						}				
					}
				else {
					if(possibleZeroMLE==true) nextTry=effectiveMin;
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
	/*
	while(*midL>*bestOuterL){
		
		if(curveDif > 0 && slopeCB < 0){
			//if this is the case, the surface is either convex or not concave enough
			//across our evals, and a parabolic estimate will suck (or fail outright), so do a GOLD
			nextTry=((*bestOuter)+2.0*(*bestOuter-*mid));
			tryPara=false;
			}
		else if(possibleZeroMLE==true){
			nextTry=effectiveMin;
			}
		else{
			r=(*mid-*worstOuter)*(*midL-*bestOuterL);
			q=(*mid-*bestOuter)*(*midL-*worstOuterL);			
			nextTry=(*mid)-((*mid-*bestOuter)*q-(*mid-*worstOuter)*r)/(2.0*SIGN(FMAX(fabs(q-r),TINY), q-r));
			}
		
		if(nextTry<effectiveMin){
			if(*bestOuter==effectiveMin){
				return 1;
				}
			nextTry=effectiveMin;
			nextTryL=(*func)(thisnode, thistree, nextTry, true);
			if(nextTryL < *bestOuterL){
				return 1;
				}
			else if(nextTryL > *midL){
				*bestOuter=nextTry;
				*bestOuterL=nextTryL;
				return 0;
				}
			}
		else{
			ulim=(*mid)+GLIMIT*(*bestOuter-*mid);
			if(((*mid-nextTry)*(nextTry-*bestOuter)>0.0)){	//if the proposed point is between mid and bestOuter
				nextTryL=(*func)(thisnode, thistree, nextTry, true);
				if(nextTryL < *midL){			//if the new score is better than our second best
					*worstOuter=*mid;			//b becomes a, nextTry becomes b
					*mid=nextTry;
					*worstOuterL=*midL;
					*midL=nextTryL;
					if(nextTryL < *bestOuterL){		//if the new score is also better that the previous best
						return 0;		//we've found our bracket.  Return.
						}
					}
				else {//this shouldn't really happen.  Implies multiple optima
					cout.precision(11);
					cout << "multiple optima????\n";
					cout << "a=" << *worstOuter << " worstOuterL=" << *worstOuterL << endl;
					cout << "b=" << *mid << " midL=" << *midL << endl;
					cout << "c=" << *bestOuter << " bestOuterL=" << *bestOuterL << endl;
					cout << "nextTry=" << nextTry << " nextTryL=" << nextTryL << endl;
					*bestOuter=nextTry;
					*bestOuterL=nextTryL;
					return 0;
					}
		/*		if(*bestOuter==effectiveMin) return 1;//if we just took a parabolic step, but our best remains at the min,
											   //assume that we have a zero MLE
				if(*bestOuter > *mid) nextTry=(*bestOuter)+2.0*(*bestOuter-*mid);
				else nextTry=(*bestOuter)+GOLD*(*bestOuter-*mid);  //the parabolic step didn't do so well, so do a gold step
				nextTry = (nextTry > effectiveMin ? (nextTry < effectiveMworstOuter ? nextTry : effectiveMworstOuter) : effectiveMin);
				nextTryL=(*func)(thisnode, thistree, nextTry, true);
		*/	/*	}
			else if((nextTry < *bestOuter && *bestOuter < *mid)){ //if the proposed parabolic estimate lies to the left of our current evaluations
				nextTryL=(*func)(thisnode, thistree, nextTry, true);
				if(nextTryL > *bestOuterL){ //we found our bracket
					SHFT(*worstOuter, *mid, *bestOuter, nextTry)
					SHFT(*worstOuterL, *midL, *bestOuterL, nextTryL)
*/	/*				*worstOuter=*mid;
					*mid=nextTry;
					*worstOuterL=*midL;
					*midL=nextTryL;*/
/*					return 0;
					}
				SHFT(*worstOuter, *mid, *bestOuter, nextTry)
				SHFT(*worstOuterL, *midL, *bestOuterL, nextTryL)
				}
			else{//the proposed parabolic estimate is to the right of our current evaluations, where it tends 
				//to undershoot.  GOLD tends to be too conservative as well.  Try a GOLD style step, but with
				//the multiplier equal to 2*GOLD
				//nextTry=((*bestOuter)+GOLD*(*bestOuter-*mid));
			//	nextTry=nextTry+GOLD*(nextTry-*mid);
			//	nextTry = (nextTry > effectiveMin ? (nextTry < effectiveMworstOuter ? nextTry : effectiveMworstOuter) : effectiveMin);
				//try the parabolic estimate
				if(nextTry>.84) tryPara=false;
				if(tryPara==false){
					nextTry=((*bestOuter)+2.0*(*bestOuter-*mid));
					}
				tryPara = !tryPara;
				nextTryL=(*func)(thisnode, thistree, nextTry, true);
				//I think these shifts need to be within this loop, not outside of the below closing bracket
				SHFT(*worstOuter, *mid, *bestOuter, nextTry)
				SHFT(*worstOuterL, *midL, *bestOuterL, nextTryL)
				}
			if(*bestOuter==effectiveMin && *bestOuterL < *midL){
				numAttemptsWithBestAtMin++;
				}
			if(numAttemptsWithBestAtMin>=1 && *mid<.08){
				return 1;
				}
			}
	/*	if(((*worstOuter < -18.42) && (*mid < -18.42)) || ((*worstOuter<-10) && (*mid<-10) && (*bestOuter>1))){
			//DZ 12-18-03 if our three best points are all < ln(1e-8), just give up and take that as a blen
			//the MLE is probably essentially 0.  Note that sometimes when worstOuter and mid are very small this 
			//func tries very large values for bestOuter, which I think is a bug.  This hack also avoids that
			return 1;
			}
		}
	return 0;
*/	}
	
//I'm reworking this a bit to better use the information that has already been generated in the bracketing function
//since we already have those function evaluations, we might as well pass them in and use them
double DZbrent(double ax, double bx, double cx, double fa, double fx, double fc, double (*f)(TreeNode *, Tree*, double, bool), double tol, double *xmin, TreeNode *thisnode, Tree *thistree){
	 int iter;
 	 double a, b, d, etemp, fu, fv, fw/*, fx*/, p, q, r, tol1, tol2, u, v, w, x, xm;
	 double e=0.0;
	 
	 if((fx<fa && fx<fc)==false){
	 	//if bx isn't the current minimum, make is so
	 	if(fa<fx){
	 		double dummy=fa;
	 		fa=fx;
	 		fx=dummy;
	 		dummy=ax;
	 		ax=bx;
	 		bx=dummy;
	 		}
	 	else if(fc<fx){
	 		double dummy=fc;
	 		fc=fx;
	 		fx=dummy;
	 		dummy=cx;
	 		cx=bx;
	 		bx=dummy;	 		
	 		}
	 	}
	 assert(fx<fa && fx<fc);
	 
	 double paraMinlnL, paraErr, paraErrCrit;
	 paraErrCrit=(tol<.5 ? tol*10 : 5);
	 bool paraOK=false, para=false;

//	ofstream brak("brakdebug.log", ios::app);
//	brak << "node " << thisnode->nodeNum << "\n";

	 a=(ax < cx ? ax : cx); //make a the smallest of the three bracket points 
	 b=(ax > cx ? ax : cx); //and b the largest
	if(ax>cx){
		double dummy=fa;
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

	xm=0.5*(a+b);       //xm is the midpoint of the bracket (of a and b)
	e=(x>=xm?a-x:b-x);	//set e to the larger of the two bracket intervals
	d=CGOLD*e;

//	assert(a<=x && x<=b);
	 
//	 fw=fv=fx=(*f)(thisnode, thistree, x);
	 
	 for(iter=1;iter<=ITMAX;iter++){

	 	xm=0.5*(a+b);		//xm is the midpoint of the bracket (of a and b)
	 	
	 	tol2=2.0*(tol1=tol*fabs(x)+ZEPS);	
	 	
/*	 	if (fabs(x-xm) <= (tol2-0.5*(b-a))){ //termination condition
	 		*xmin=x;						
	 		return fx;
	 		}
*/// 	if (fabs(e) > tol1){	//construct a trial parabolic fit
	 		r=(x-w)*(fx-fv);
	 		q=(x-v)*(fx-fw);
	 		p=(x-v)*q-(x-w)*r;
	 		q=2.0*(q-r);
	 		if(q>0.0) p=-p;
	 		q=fabs(q);
	 		etemp=e;
	 		e=d;
	 		if(fabs(p) >= fabs(0.5*q*etemp)||p<=q*(a-x) || p>=q*(b-x)){ //determine if the parabolic fit is good
	 			d=CGOLD*(e=(x>=xm?a-x:b-x));  //if not
	 			u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
	 			}
	 			
	 		else{				//if so, take the parabolic step
	 			d=p/q;
	 			u=x+d;
	 			
	 			double alph=w-u;
				double beta=x-u; 
				paraMinlnL=(((fx) * alph*alph) - ((fw) * beta*beta)) / (alph*alph - beta*beta);

	 			if(paraOK==true){
	 				//the estimation error in the parabolic step always seems to at least half each iteration,
	 				//hence the division by 2.0
	 				double estlnL=paraMinlnL - paraErr/2.0;
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
	 
void InferStatesFromCla(char *states, double *cla, int nchar){
	//this function takes a cla that contains the contribution of the whole tree
	//and calculates the most probable state at each site.  The resulting array of
	//states is placed into *states, which should already be allocated.
	double stateProbs[4];
	int nrates=4;
	
	for(int c=0;c<nchar;c++){
		stateProbs[0]=stateProbs[1]=stateProbs[2]=stateProbs[3]=0.0;
		for(int i=0;i<nrates;i++){
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

double CalculatePDistance(const char *str1, const char *str2, int nchar){
	double count=0.0;
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
	return count/(double)effectiveChar;
	}

#ifndef GANESH
double CalculateHammingDistance(const char *str1, const char *str2, int nchar){
	double count=0.0;
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
	return count/(double)effectiveChar;
	}
#else
double CalculateHammingDistance(const char *str1, const char *str2, 
                                const int *col_count, int nchar){
	double count=0.0;
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
	return count/(double)effectiveChar;
}
#endif

void SampleBranchLengthCurve(double (*func)(TreeNode*, Tree*, double, bool), TreeNode *thisnode, Tree *thistree){
	for(double len=effectiveMin;len<effectiveMax;len*=2.0)
		(*func)(thisnode, thistree, len, true);
	}

void CalcDerivCLASPartialInternalRateHet(CondLikeArray *destCLA, CondLikeArray *destD1CLA, CondLikeArray *destD2CLA, const CondLikeArray *partialCLA, const CondLikeArray *childCLA,
	const double *prmat, const double *d1mat, const double *d2mat, int nchar){
	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.
	double *dest=destCLA->arr;
	double *destD1=destD1CLA->arr;
	double *destD2=destD2CLA->arr;
	double *CL1=childCLA->arr;
	double *partial=partialCLA->arr;
	
#ifdef UNIX
	madvise(dest, nchar*16*sizeof(double), MADV_SEQUENTIAL);
	madvise((void*)CL1, nchar*16*sizeof(double), MADV_SEQUENTIAL);
	madvise(partial, nchar*16*sizeof(double), MADV_SEQUENTIAL);
#endif

//	register double c0, c1, c2, c3;
/*	for(int i=0;i<nchar;i++){
		c0=CL1[0];
		c1=CL1[1];
		c2=CL1[2];
		c3=CL1[3];
		dest[0] = ( prmat[0]*c0+prmat[1]*c1+prmat[2]*c2+prmat[3]*c3) * partial[0];
		dest[1] = ( prmat[4]*c0+prmat[5]*c1+prmat[6]*c2+prmat[7]*c3) * partial[1];
		dest[2] = ( prmat[8]*c0+prmat[9]*c1+prmat[10]*c2+prmat[11]*c3) * partial[2];
		dest[3] = ( prmat[12]*c0+prmat[13]*c1+prmat[14]*c2+prmat[15]*c3) * partial[3];

		destD1[0] = ( d1mat[0]*c0+d1mat[1]*c1+d1mat[2]*c2+d1mat[3]*c3) * partial[0];
		destD1[1] = ( d1mat[4]*c0+d1mat[5]*c1+d1mat[6]*c2+d1mat[7]*c3) * partial[1];
		destD1[2] = ( d1mat[8]*c0+d1mat[9]*c1+d1mat[10]*c2+d1mat[11]*c3) * partial[2];
		destD1[3] = ( d1mat[12]*c0+d1mat[13]*c1+d1mat[14]*c2+d1mat[15]*c3) * partial[3];

		destD2[0] = ( d2mat[0]*c0+d2mat[1]*c1+d2mat[2]*c2+d2mat[3]*c3) * partial[0];
		destD2[1] = ( d2mat[4]*c0+d2mat[5]*c1+d2mat[6]*c2+d2mat[7]*c3) * partial[1];
		destD2[2] = ( d2mat[8]*c0+d2mat[9]*c1+d2mat[10]*c2+d2mat[11]*c3) * partial[2];
		destD2[3] = ( d2mat[12]*c0+d2mat[13]*c1+d2mat[14]*c2+d2mat[15]*c3) * partial[3];

		c0=CL1[4];
		c1=CL1[5];
		c2=CL1[6];
		c3=CL1[7];
		dest[4] = ( prmat[16]*c0+prmat[17]*c1+prmat[18]*c2+prmat[19]*c3) * partial[4];
		dest[5] = ( prmat[20]*c0+prmat[21]*c1+prmat[22]*c2+prmat[23]*c3) * partial[5];
		dest[6] = ( prmat[24]*c0+prmat[25]*c1+prmat[26]*c2+prmat[27]*c3) * partial[6];
		dest[7] = ( prmat[28]*c0+prmat[29]*c1+prmat[30]*c2+prmat[31]*c3) * partial[7];

		destD1[4] = ( d1mat[16]*c0+d1mat[17]*c1+d1mat[18]*c2+d1mat[19]*c3) * partial[4];
		destD1[5] = ( d1mat[20]*c0+d1mat[21]*c1+d1mat[22]*c2+d1mat[23]*c3) * partial[5];
		destD1[6] = ( d1mat[24]*c0+d1mat[25]*c1+d1mat[26]*c2+d1mat[27]*c3) * partial[6];
		destD1[7] = ( d1mat[28]*c0+d1mat[29]*c1+d1mat[30]*c2+d1mat[31]*c3) * partial[7];

		destD2[4] = ( d2mat[16]*c0+d2mat[17]*c1+d2mat[18]*c2+d2mat[19]*c3) * partial[4];
		destD2[5] = ( d2mat[20]*c0+d2mat[21]*c1+d2mat[22]*c2+d2mat[23]*c3) * partial[5];
		destD2[6] = ( d2mat[24]*c0+d2mat[25]*c1+d2mat[26]*c2+d2mat[27]*c3) * partial[6];
		destD2[7] = ( d2mat[28]*c0+d2mat[29]*c1+d2mat[30]*c2+d2mat[31]*c3) * partial[7];

		c0=CL1[8];
		c1=CL1[9];
		c2=CL1[10];
		c3=CL1[11];
		dest[8] = ( prmat[32]*c0+prmat[33]*c1+prmat[34]*c2+prmat[35]*c3) * partial[8];
		dest[9] = ( prmat[36]*c0+prmat[37]*c1+prmat[38]*c2+prmat[39]*c3) * partial[9];
		dest[10] = ( prmat[40]*c0+prmat[41]*c1+prmat[42]*c2+prmat[43]*c3) * partial[10];
		dest[11] = ( prmat[44]*c0+prmat[45]*c1+prmat[46]*c2+prmat[47]*c3) * partial[11];

		destD1[8] = ( d1mat[32]*c0+d1mat[33]*c1+d1mat[34]*c2+d1mat[35]*c3) * partial[8];
		destD1[9] = ( d1mat[36]*c0+d1mat[37]*c1+d1mat[38]*c2+d1mat[39]*c3) * partial[9];
		destD1[10] = ( d1mat[40]*c0+d1mat[41]*c1+d1mat[42]*c2+d1mat[43]*c3) * partial[10];
		destD1[11] = ( d1mat[44]*c0+d1mat[45]*c1+d1mat[46]*c2+d1mat[47]*c3) * partial[11];

		destD2[8] = ( d2mat[32]*c0+d2mat[33]*c1+d2mat[34]*c2+d2mat[35]*c3) * partial[8];
		destD2[9] = ( d2mat[36]*c0+d2mat[37]*c1+d2mat[38]*c2+d2mat[39]*c3) * partial[9];
		destD2[10] = ( d2mat[40]*c0+d2mat[41]*c1+d2mat[42]*c2+d2mat[43]*c3) * partial[10];
		destD2[11] = ( d2mat[44]*c0+d2mat[45]*c1+d2mat[46]*c2+d2mat[47]*c3) * partial[11];

		c0=CL1[12];
		c1=CL1[13];
		c2=CL1[14];
		c3=CL1[15];
		dest[12] = ( prmat[48]*c0+prmat[49]*c1+prmat[50]*c2+prmat[51]*c3) * partial[12];
		dest[13] = ( prmat[52]*c0+prmat[53]*c1+prmat[54]*c2+prmat[55]*c3) * partial[13];
		dest[14] = ( prmat[56]*c0+prmat[57]*c1+prmat[58]*c2+prmat[59]*c3) * partial[14];
		dest[15] = ( prmat[60]*c0+prmat[61]*c1+prmat[62]*c2+prmat[63]*c3) * partial[15];

		destD1[12] = ( d1mat[48]*c0+d1mat[49]*c1+d1mat[50]*c2+d1mat[51]*c3) * partial[12];
		destD1[13] = ( d1mat[52]*c0+d1mat[53]*c1+d1mat[54]*c2+d1mat[55]*c3) * partial[13];
		destD1[14] = ( d1mat[56]*c0+d1mat[57]*c1+d1mat[58]*c2+d1mat[59]*c3) * partial[14];
		destD1[15] = ( d1mat[60]*c0+d1mat[61]*c1+d1mat[62]*c2+d1mat[63]*c3) * partial[15];
		
		destD2[12] = ( d2mat[48]*c0+d2mat[49]*c1+d2mat[50]*c2+d2mat[51]*c3) * partial[12];
		destD2[13] = ( d2mat[52]*c0+d2mat[53]*c1+d2mat[54]*c2+d2mat[55]*c3) * partial[13];
		destD2[14] = ( d2mat[56]*c0+d2mat[57]*c1+d2mat[58]*c2+d2mat[59]*c3) * partial[14];
		destD2[15] = ( d2mat[60]*c0+d2mat[61]*c1+d2mat[62]*c2+d2mat[63]*c3) * partial[15];
		partial+=16;
		CL1+=16;
		destD1+=16;
		destD2+=16;
		dest+=16;
		}
*/


	for(int i=0;i<nchar;i++){
		dest[0] = ( prmat[0]*CL1[0]+prmat[1]*CL1[1]+prmat[2]*CL1[2]+prmat[3]*CL1[3]) * partial[0];
		dest[1] = ( prmat[4]*CL1[0]+prmat[5]*CL1[1]+prmat[6]*CL1[2]+prmat[7]*CL1[3]) * partial[1];
		dest[2] = ( prmat[8]*CL1[0]+prmat[9]*CL1[1]+prmat[10]*CL1[2]+prmat[11]*CL1[3]) * partial[2];
		dest[3] = ( prmat[12]*CL1[0]+prmat[13]*CL1[1]+prmat[14]*CL1[2]+prmat[15]*CL1[3]) * partial[3];
		
		dest[4] = ( prmat[16]*CL1[4]+prmat[17]*CL1[5]+prmat[18]*CL1[6]+prmat[19]*CL1[7]) * partial[4];
		dest[5] = ( prmat[20]*CL1[4]+prmat[21]*CL1[5]+prmat[22]*CL1[6]+prmat[23]*CL1[7]) * partial[5];
		dest[6] = ( prmat[24]*CL1[4]+prmat[25]*CL1[5]+prmat[26]*CL1[6]+prmat[27]*CL1[7]) * partial[6];
		dest[7] = ( prmat[28]*CL1[4]+prmat[29]*CL1[5]+prmat[30]*CL1[6]+prmat[31]*CL1[7]) * partial[7];

		dest[8] = ( prmat[32]*CL1[8]+prmat[33]*CL1[9]+prmat[34]*CL1[10]+prmat[35]*CL1[11]) * partial[8];
		dest[9] = ( prmat[36]*CL1[8]+prmat[37]*CL1[9]+prmat[38]*CL1[10]+prmat[39]*CL1[11]) * partial[9];
		dest[10] = ( prmat[40]*CL1[8]+prmat[41]*CL1[9]+prmat[42]*CL1[10]+prmat[43]*CL1[11]) * partial[10];
		dest[11] = ( prmat[44]*CL1[8]+prmat[45]*CL1[9]+prmat[46]*CL1[10]+prmat[47]*CL1[11]) * partial[11];

		dest[12] = ( prmat[48]*CL1[12]+prmat[49]*CL1[13]+prmat[50]*CL1[14]+prmat[51]*CL1[15]) * partial[12];
		dest[13] = ( prmat[52]*CL1[12]+prmat[53]*CL1[13]+prmat[54]*CL1[14]+prmat[55]*CL1[15]) * partial[13];
		dest[14] = ( prmat[56]*CL1[12]+prmat[57]*CL1[13]+prmat[58]*CL1[14]+prmat[59]*CL1[15]) * partial[14];
		dest[15] = ( prmat[60]*CL1[12]+prmat[61]*CL1[13]+prmat[62]*CL1[14]+prmat[63]*CL1[15]) * partial[15];
			
		destD1[0] = ( d1mat[0]*CL1[0]+d1mat[1]*CL1[1]+d1mat[2]*CL1[2]+d1mat[3]*CL1[3]) * partial[0];
		destD1[1] = ( d1mat[4]*CL1[0]+d1mat[5]*CL1[1]+d1mat[6]*CL1[2]+d1mat[7]*CL1[3]) * partial[1];
		destD1[2] = ( d1mat[8]*CL1[0]+d1mat[9]*CL1[1]+d1mat[10]*CL1[2]+d1mat[11]*CL1[3]) * partial[2];
		destD1[3] = ( d1mat[12]*CL1[0]+d1mat[13]*CL1[1]+d1mat[14]*CL1[2]+d1mat[15]*CL1[3]) * partial[3];

		destD1[4] = ( d1mat[16]*CL1[4]+d1mat[17]*CL1[5]+d1mat[18]*CL1[6]+d1mat[19]*CL1[7]) * partial[4];
		destD1[5] = ( d1mat[20]*CL1[4]+d1mat[21]*CL1[5]+d1mat[22]*CL1[6]+d1mat[23]*CL1[7]) * partial[5];
		destD1[6] = ( d1mat[24]*CL1[4]+d1mat[25]*CL1[5]+d1mat[26]*CL1[6]+d1mat[27]*CL1[7]) * partial[6];
		destD1[7] = ( d1mat[28]*CL1[4]+d1mat[29]*CL1[5]+d1mat[30]*CL1[6]+d1mat[31]*CL1[7]) * partial[7];

		destD1[8] = ( d1mat[32]*CL1[8]+d1mat[33]*CL1[9]+d1mat[34]*CL1[10]+d1mat[35]*CL1[11]) * partial[8];
		destD1[9] = ( d1mat[36]*CL1[8]+d1mat[37]*CL1[9]+d1mat[38]*CL1[10]+d1mat[39]*CL1[11]) * partial[9];
		destD1[10] = ( d1mat[40]*CL1[8]+d1mat[41]*CL1[9]+d1mat[42]*CL1[10]+d1mat[43]*CL1[11]) * partial[10];
		destD1[11] = ( d1mat[44]*CL1[8]+d1mat[45]*CL1[9]+d1mat[46]*CL1[10]+d1mat[47]*CL1[11]) * partial[11];

		destD1[12] = ( d1mat[48]*CL1[12]+d1mat[49]*CL1[13]+d1mat[50]*CL1[14]+d1mat[51]*CL1[15]) * partial[12];
		destD1[13] = ( d1mat[52]*CL1[12]+d1mat[53]*CL1[13]+d1mat[54]*CL1[14]+d1mat[55]*CL1[15]) * partial[13];
		destD1[14] = ( d1mat[56]*CL1[12]+d1mat[57]*CL1[13]+d1mat[58]*CL1[14]+d1mat[59]*CL1[15]) * partial[14];
		destD1[15] = ( d1mat[60]*CL1[12]+d1mat[61]*CL1[13]+d1mat[62]*CL1[14]+d1mat[63]*CL1[15]) * partial[15];

		destD2[0] = ( d2mat[0]*CL1[0]+d2mat[1]*CL1[1]+d2mat[2]*CL1[2]+d2mat[3]*CL1[3]) * partial[0];
		destD2[1] = ( d2mat[4]*CL1[0]+d2mat[5]*CL1[1]+d2mat[6]*CL1[2]+d2mat[7]*CL1[3]) * partial[1];
		destD2[2] = ( d2mat[8]*CL1[0]+d2mat[9]*CL1[1]+d2mat[10]*CL1[2]+d2mat[11]*CL1[3]) * partial[2];
		destD2[3] = ( d2mat[12]*CL1[0]+d2mat[13]*CL1[1]+d2mat[14]*CL1[2]+d2mat[15]*CL1[3]) * partial[3];
		
		destD2[4] = ( d2mat[16]*CL1[4]+d2mat[17]*CL1[5]+d2mat[18]*CL1[6]+d2mat[19]*CL1[7]) * partial[4];
		destD2[5] = ( d2mat[20]*CL1[4]+d2mat[21]*CL1[5]+d2mat[22]*CL1[6]+d2mat[23]*CL1[7]) * partial[5];
		destD2[6] = ( d2mat[24]*CL1[4]+d2mat[25]*CL1[5]+d2mat[26]*CL1[6]+d2mat[27]*CL1[7]) * partial[6];
		destD2[7] = ( d2mat[28]*CL1[4]+d2mat[29]*CL1[5]+d2mat[30]*CL1[6]+d2mat[31]*CL1[7]) * partial[7];
		
		destD2[8] = ( d2mat[32]*CL1[8]+d2mat[33]*CL1[9]+d2mat[34]*CL1[10]+d2mat[35]*CL1[11]) * partial[8];
		destD2[9] = ( d2mat[36]*CL1[8]+d2mat[37]*CL1[9]+d2mat[38]*CL1[10]+d2mat[39]*CL1[11]) * partial[9];
		destD2[10] = ( d2mat[40]*CL1[8]+d2mat[41]*CL1[9]+d2mat[42]*CL1[10]+d2mat[43]*CL1[11]) * partial[10];
		destD2[11] = ( d2mat[44]*CL1[8]+d2mat[45]*CL1[9]+d2mat[46]*CL1[10]+d2mat[47]*CL1[11]) * partial[11];
		
		destD2[12] = ( d2mat[48]*CL1[12]+d2mat[49]*CL1[13]+d2mat[50]*CL1[14]+d2mat[51]*CL1[15]) * partial[12];
		destD2[13] = ( d2mat[52]*CL1[12]+d2mat[53]*CL1[13]+d2mat[54]*CL1[14]+d2mat[55]*CL1[15]) * partial[13];
		destD2[14] = ( d2mat[56]*CL1[12]+d2mat[57]*CL1[13]+d2mat[58]*CL1[14]+d2mat[59]*CL1[15]) * partial[14];
		destD2[15] = ( d2mat[60]*CL1[12]+d2mat[61]*CL1[13]+d2mat[62]*CL1[14]+d2mat[63]*CL1[15]) * partial[15];
		partial+=16;
		CL1+=16;
		destD1+=16;
		destD2+=16;
		dest+=16;
		}
	}


void CalcDerivCLASPartialTerminalRateHet(CondLikeArray *destCLA, CondLikeArray *destD1CLA, CondLikeArray *destD2CLA, const CondLikeArray *partialCLA,
	const double *prmat, const double *d1mat, const double *d2mat, const char *Ldata, int nchar){

	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.
	double *dest=destCLA->arr;
	double *destD1=destD1CLA->arr;
	double *destD2=destD2CLA->arr;
	double *partial=partialCLA->arr;
#ifdef UNIX
	madvise(dest, nchar*16*sizeof(double), MADV_SEQUENTIAL);
	madvise((void*)partial, nchar*16*sizeof(double), MADV_SEQUENTIAL);
#endif

	for(int i=0;i<nchar;i++){
		if(*Ldata > -1){ //no ambiguity
			for(int i=0;i<4;i++){
				*(dest++) = prmat[(*Ldata)+16*i] * *(partial);
				*(destD1++) = d1mat[(*Ldata)+16*i] * *(partial);
				*(destD2++) = d2mat[(*Ldata)+16*i] * *(partial);
				partial++;
				*(dest++) = prmat[(*Ldata+4)+16*i] * *(partial);
				*(destD1++) = d1mat[(*Ldata+4)+16*i] * *(partial);
				*(destD2++) = d2mat[(*Ldata+4)+16*i] * *(partial);
				partial++;
				*(dest++) = prmat[(*Ldata+8)+16*i] * *(partial);
				*(destD1++) = d1mat[(*Ldata+8)+16*i] * *(partial);
				*(destD2++) = d2mat[(*Ldata+8)+16*i] * *(partial);
				partial++;
				*(dest++) = prmat[(*Ldata+12)+16*i] * *(partial);
				*(destD1++) = d1mat[(*Ldata+12)+16*i] * *(partial);
				*(destD2++) = d2mat[(*Ldata+12)+16*i] * *(partial);
				partial++;
				}
			Ldata++;
			}
			
		else if(*Ldata == -4){ //total ambiguity
			for(int i=0;i<16;i++){
				*(dest++) = *(partial++);
				*(destD1++) = 0.0;
				*(destD2++) = 0.0;
				}
			Ldata++;
			}
		else{ //partial ambiguity
			//first figure in the ambiguous terminal
			char nstates=-1 * *(Ldata++);
			for(int q=0;q<16;q++){
				dest[q]=0;
				destD1[q]=0;
				destD2[q]=0;
				}
			for(int i=0;i<nstates;i++){
				for(int i=0;i<4;i++){
					*(dest+(i*4)) += prmat[(*Ldata)+16*i];
					*(dest+(i*4)+1) += prmat[(*Ldata+4)+16*i];
					*(dest+(i*4)+2) += prmat[(*Ldata+8)+16*i];
					*(dest+(i*4)+3) += prmat[(*Ldata+12)+16*i];
					
					*(destD1+(i*4)) += d1mat[(*Ldata)+16*i];
					*(destD1+(i*4)+1) += d1mat[(*Ldata+4)+16*i];
					*(destD1+(i*4)+2) += d1mat[(*Ldata+8)+16*i];
					*(destD1+(i*4)+3) += d1mat[(*Ldata+12)+16*i];

					*(destD2+(i*4)) += d2mat[(*Ldata)+16*i];
					*(destD2+(i*4)+1) += d2mat[(*Ldata+4)+16*i];
					*(destD2+(i*4)+2) += d2mat[(*Ldata+8)+16*i];
					*(destD2+(i*4)+3) += d2mat[(*Ldata+12)+16*i];
					}
				Ldata++;
				}
			
			//now add the partial
			for(int i=0;i<4;i++){
				*(dest++) *= *(partial);
				*(destD1++) *= *(partial);
				*(destD2++) *= *(partial);
				partial++;
				*(dest++) *= *(partial);
				*(destD1++) *= *(partial);
				*(destD2++) *= *(partial);
				partial++;
				*(dest++) *= *(partial);
				*(destD1++) *= *(partial);
				*(destD2++) *= *(partial);
				partial++;
				*(dest++) *= *(partial);
				*(destD1++) *= *(partial);
				*(destD2++) *= *(partial);
				partial++;			
				}
			}
	
		}
	}

void CalcFullCLAInternalInternalRateHet(CondLikeArray *destCLA, const CondLikeArray *LCLA, const CondLikeArray *RCLA, const double *Lpr, const double *Rpr, int nchar){
	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.
	double *dest=destCLA->arr;
	double *LCL=LCLA->arr;
	double *RCL=RCLA->arr;
	
#ifdef UNIX
	madvise(dest, nchar*16*sizeof(double), MADV_SEQUENTIAL);
	madvise((void *)LCL, nchar*16*sizeof(double), MADV_SEQUENTIAL);
	madvise((void *)RCL, nchar*16*sizeof(double), MADV_SEQUENTIAL);
#endif
/*
	double L1, L2, L3, L4, R1, R2, R3, R4;
	for(int i=0;i<nchar;i++){
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

		L1=( Lpr[16+0]*LCL[4+0]+Lpr[16+1]*LCL[4+1]+Lpr[16+2]*LCL[4+2]+Lpr[16+3]*LCL[4+3]);
		L2=( Lpr[16+4]*LCL[4+0]+Lpr[16+5]*LCL[4+1]+Lpr[16+6]*LCL[4+2]+Lpr[16+7]*LCL[4+3]);
		L3=( Lpr[16+8]*LCL[4+0]+Lpr[16+9]*LCL[4+1]+Lpr[16+10]*LCL[4+2]+Lpr[16+11]*LCL[4+3]);
		L4=( Lpr[16+12]*LCL[4+0]+Lpr[16+13]*LCL[4+1]+Lpr[16+14]*LCL[4+2]+Lpr[16+15]*LCL[4+3]);

		R1=(Rpr[16+0]*RCL[4+0]+Rpr[16+1]*RCL[4+1]+Rpr[16+2]*RCL[4+2]+Rpr[16+3]*RCL[4+3]);
		R2=(Rpr[16+4]*RCL[4+0]+Rpr[16+5]*RCL[4+1]+Rpr[16+6]*RCL[4+2]+Rpr[16+7]*RCL[4+3]);
		R3=(Rpr[16+8]*RCL[4+0]+Rpr[16+9]*RCL[4+1]+Rpr[16+10]*RCL[4+2]+Rpr[16+11]*RCL[4+3]);			
		R4=(Rpr[16+12]*RCL[4+0]+Rpr[16+13]*RCL[4+1]+Rpr[16+14]*RCL[4+2]+Rpr[16+15]*RCL[4+3]);
		
		dest[4] = L1 * R1;
		dest[5] = L2 * R2;
		dest[6] = L3 * R3;
		dest[7] = L4 * R4;

		L1=( Lpr[32+0]*LCL[8+0]+Lpr[32+1]*LCL[8+1]+Lpr[32+2]*LCL[8+2]+Lpr[32+3]*LCL[8+3]);
		L2=( Lpr[32+4]*LCL[8+0]+Lpr[32+5]*LCL[8+1]+Lpr[32+6]*LCL[8+2]+Lpr[32+7]*LCL[8+3]);
		L3=( Lpr[32+8]*LCL[8+0]+Lpr[32+9]*LCL[8+1]+Lpr[32+10]*LCL[8+2]+Lpr[32+11]*LCL[8+3]);
		L4=( Lpr[32+12]*LCL[8+0]+Lpr[32+13]*LCL[8+1]+Lpr[32+14]*LCL[8+2]+Lpr[32+15]*LCL[8+3]);

		R1=(Rpr[32+0]*RCL[8+0]+Rpr[32+1]*RCL[8+1]+Rpr[32+2]*RCL[8+2]+Rpr[32+3]*RCL[8+3]);
		R2=(Rpr[32+4]*RCL[8+0]+Rpr[32+5]*RCL[8+1]+Rpr[32+6]*RCL[8+2]+Rpr[32+7]*RCL[8+3]);
		R3=(Rpr[32+8]*RCL[8+0]+Rpr[32+9]*RCL[8+1]+Rpr[32+10]*RCL[8+2]+Rpr[32+11]*RCL[8+3]);			
		R4=(Rpr[32+12]*RCL[8+0]+Rpr[32+13]*RCL[8+1]+Rpr[32+14]*RCL[8+2]+Rpr[32+15]*RCL[8+3]);
		
		dest[8] = L1 * R1;
		dest[9] = L2 * R2;
		dest[10] = L3 * R3;
		dest[11] = L4 * R4;

		L1=( Lpr[48+0]*LCL[12+0]+Lpr[48+1]*LCL[12+1]+Lpr[48+2]*LCL[12+2]+Lpr[48+3]*LCL[12+3]);
		L2=( Lpr[48+4]*LCL[12+0]+Lpr[48+5]*LCL[12+1]+Lpr[48+6]*LCL[12+2]+Lpr[48+7]*LCL[12+3]);
		L3=( Lpr[48+8]*LCL[12+0]+Lpr[48+9]*LCL[12+1]+Lpr[48+10]*LCL[12+2]+Lpr[48+11]*LCL[12+3]);
		L4=( Lpr[48+12]*LCL[12+0]+Lpr[48+13]*LCL[12+1]+Lpr[48+14]*LCL[12+2]+Lpr[48+15]*LCL[12+3]);

		R1=(Rpr[48+0]*RCL[12+0]+Rpr[48+1]*RCL[12+1]+Rpr[48+2]*RCL[12+2]+Rpr[48+3]*RCL[12+3]);
		R2=(Rpr[48+4]*RCL[12+0]+Rpr[48+5]*RCL[12+1]+Rpr[48+6]*RCL[12+2]+Rpr[48+7]*RCL[12+3]);
		R3=(Rpr[48+8]*RCL[12+0]+Rpr[48+9]*RCL[12+1]+Rpr[48+10]*RCL[12+2]+Rpr[48+11]*RCL[12+3]);			
		R4=(Rpr[48+12]*RCL[12+0]+Rpr[48+13]*RCL[12+1]+Rpr[48+14]*RCL[12+2]+Rpr[48+15]*RCL[12+3]);
		
		dest[12] = L1 * R1;
		dest[13] = L2 * R2;
		dest[14] = L3 * R3;
		dest[15] = L4 * R4;

		dest+=16;
		LCL+=16;
		RCL+=16;
		}
*/

	
	double L1, L2, L3, L4, R1, R2, R3, R4;
	for(int i=0;i<nchar;i++){
		
		for(int r=0;r<4;r++){
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

			assert(dest[0] < 1e50);
			assert(dest[1] < 1e50);
			assert(dest[2] < 1e50);
			assert(dest[3] < 1e50);

			dest+=4;
			LCL+=4;
			RCL+=4;
			}
		}

/*
	for(int i=0;i<nchar;i++){
		dest[0] = ( Lpr[0]*LCL[0]+Lpr[1]*LCL[1]+Lpr[2]*LCL[2]+Lpr[3]*LCL[3]) * (Rpr[0]*RCL[0]+Rpr[1]*RCL[1]+Rpr[2]*RCL[2]+Rpr[3]*RCL[3]);
		dest[1] = ( Lpr[4]*LCL[0]+Lpr[5]*LCL[1]+Lpr[6]*LCL[2]+Lpr[7]*LCL[3]) * (Rpr[4]*RCL[0]+Rpr[5]*RCL[1]+Rpr[6]*RCL[2]+Rpr[7]*RCL[3]);
		dest[2] = ( Lpr[8]*LCL[0]+Lpr[9]*LCL[1]+Lpr[10]*LCL[2]+Lpr[11]*LCL[3]) * (Rpr[8]*RCL[0]+Rpr[9]*RCL[1]+Rpr[10]*RCL[2]+Rpr[11]*RCL[3]);
		dest[3] = ( Lpr[12]*LCL[0]+Lpr[13]*LCL[1]+Lpr[14]*LCL[2]+Lpr[15]*LCL[3]) * (Rpr[12]*RCL[0]+Rpr[13]*RCL[1]+Rpr[14]*RCL[2]+Rpr[15]*RCL[3]);
		
		dest[4] = ( Lpr[16]*LCL[4]+Lpr[17]*LCL[5]+Lpr[18]*LCL[6]+Lpr[19]*LCL[7]) * (Rpr[16]*RCL[4]+Rpr[17]*RCL[5]+Rpr[18]*RCL[6]+Rpr[19]*RCL[7]);
		dest[5] = ( Lpr[20]*LCL[4]+Lpr[21]*LCL[5]+Lpr[22]*LCL[6]+Lpr[23]*LCL[7]) * (Rpr[20]*RCL[4]+Rpr[21]*RCL[5]+Rpr[22]*RCL[6]+Rpr[23]*RCL[7]);
		dest[6] = ( Lpr[24]*LCL[4]+Lpr[25]*LCL[5]+Lpr[26]*LCL[6]+Lpr[27]*LCL[7]) * (Rpr[24]*RCL[4]+Rpr[25]*RCL[5]+Rpr[26]*RCL[6]+Rpr[27]*RCL[7]);
		dest[7] = ( Lpr[28]*LCL[4]+Lpr[29]*LCL[5]+Lpr[30]*LCL[6]+Lpr[31]*LCL[7]) * (Rpr[28]*RCL[4]+Rpr[29]*RCL[5]+Rpr[30]*RCL[6]+Rpr[31]*RCL[7]);
	
		dest[8] = ( Lpr[32]*LCL[8]+Lpr[33]*LCL[9]+Lpr[34]*LCL[10]+Lpr[35]*LCL[11]) * (Rpr[32]*RCL[8]+Rpr[33]*RCL[9]+Rpr[34]*RCL[10]+Rpr[35]*RCL[11]);
		dest[9] = ( Lpr[36]*LCL[8]+Lpr[37]*LCL[9]+Lpr[38]*LCL[10]+Lpr[39]*LCL[11]) * (Rpr[36]*RCL[8]+Rpr[37]*RCL[9]+Rpr[38]*RCL[10]+Rpr[39]*RCL[11]);
		dest[10] = ( Lpr[40]*LCL[8]+Lpr[41]*LCL[9]+Lpr[42]*LCL[10]+Lpr[43]*LCL[11]) * (Rpr[40]*RCL[8]+Rpr[41]*RCL[9]+Rpr[42]*RCL[10]+Rpr[43]*RCL[11]);
		dest[11] = ( Lpr[44]*LCL[8]+Lpr[45]*LCL[9]+Lpr[46]*LCL[10]+Lpr[47]*LCL[11]) * (Rpr[44]*RCL[8]+Rpr[45]*RCL[9]+Rpr[46]*RCL[10]+Rpr[47]*RCL[11]);
	
		dest[12] = ( Lpr[48]*LCL[12]+Lpr[49]*LCL[13]+Lpr[50]*LCL[14]+Lpr[51]*LCL[15]) * (Rpr[48]*RCL[12]+Rpr[49]*RCL[13]+Rpr[50]*RCL[14]+Rpr[51]*RCL[15]);
		dest[13] = ( Lpr[52]*LCL[12]+Lpr[53]*LCL[13]+Lpr[54]*LCL[14]+Lpr[55]*LCL[15]) * (Rpr[52]*RCL[12]+Rpr[53]*RCL[13]+Rpr[54]*RCL[14]+Rpr[55]*RCL[15]);
		dest[14] = ( Lpr[56]*LCL[12]+Lpr[57]*LCL[13]+Lpr[58]*LCL[14]+Lpr[59]*LCL[15]) * (Rpr[56]*RCL[12]+Rpr[57]*RCL[13]+Rpr[58]*RCL[14]+Rpr[59]*RCL[15]);
		dest[15] = ( Lpr[60]*LCL[12]+Lpr[61]*LCL[13]+Lpr[62]*LCL[14]+Lpr[63]*LCL[15]) * (Rpr[60]*RCL[12]+Rpr[61]*RCL[13]+Rpr[62]*RCL[14]+Rpr[63]*RCL[15]);
		assert(*(dest)<1e200);
		dest+=16;
		LCL+=16;
		RCL+=16;		
		}
*/
	int *undermult=destCLA->underflow_mult;
	int *left_mult=LCLA->underflow_mult;
	int *right_mult=RCLA->underflow_mult;
	
	for(int i=0;i<nchar;i++){
		undermult[i] = left_mult[i] + right_mult[i];
		}
	destCLA->rescaleRank = 2 + LCLA->rescaleRank + RCLA->rescaleRank;
	}


void CalcFullCLAInternalInternal(CondLikeArray *destCLA, const CondLikeArray *LCLA, const CondLikeArray *RCLA, const double *Lpr, const double *Rpr, int nchar){
	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.
	double *dest=destCLA->arr;
	double *LCL=LCLA->arr;
	double *RCL=RCLA->arr;
	
#ifdef UNIX
	madvise(dest, nchar*4*sizeof(double), MADV_SEQUENTIAL);
	madvise((void *)LCL, nchar*4*sizeof(double), MADV_SEQUENTIAL);
	madvise((void *)RCL, nchar*4*sizeof(double), MADV_SEQUENTIAL);
#endif
	for(int i=0;i<nchar;i++){
		dest[0] = ( Lpr[0]*LCL[0]+Lpr[1]*LCL[1]+Lpr[2]*LCL[2]+Lpr[3]*LCL[3]) * (Rpr[0]*RCL[0]+Rpr[1]*RCL[1]+Rpr[2]*RCL[2]+Rpr[3]*RCL[3]);
		dest[1] = ( Lpr[4]*LCL[0]+Lpr[5]*LCL[1]+Lpr[6]*LCL[2]+Lpr[7]*LCL[3]) * (Rpr[4]*RCL[0]+Rpr[5]*RCL[1]+Rpr[6]*RCL[2]+Rpr[7]*RCL[3]);
		dest[2] = ( Lpr[8]*LCL[0]+Lpr[9]*LCL[1]+Lpr[10]*LCL[2]+Lpr[11]*LCL[3]) * (Rpr[8]*RCL[0]+Rpr[9]*RCL[1]+Rpr[10]*RCL[2]+Rpr[11]*RCL[3]);
		dest[3] = ( Lpr[12]*LCL[0]+Lpr[13]*LCL[1]+Lpr[14]*LCL[2]+Lpr[15]*LCL[3]) * (Rpr[12]*RCL[0]+Rpr[13]*RCL[1]+Rpr[14]*RCL[2]+Rpr[15]*RCL[3]);

		dest+=4;
		LCL+=4;
		RCL+=4;		
		}

	int *undermult=destCLA->underflow_mult;
	int *left_mult=LCLA->underflow_mult;
	int *right_mult=RCLA->underflow_mult;
	
	for(int i=0;i<nchar;i++){
		undermult[i] = left_mult[i] + right_mult[i];
		}
	destCLA->rescaleRank = 2 + LCLA->rescaleRank + RCLA->rescaleRank;
	}

void CalcFullCLATerminalTerminalRateHet(CondLikeArray *destCLA, const double *Lpr, const double *Rpr, char *Ldata, char *Rdata, int nchar){
	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.
	double *dest=destCLA->arr;
	
#ifdef UNIX
	madvise(dest, nchar*16*sizeof(double), MADV_SEQUENTIAL);
#endif

	for(int i=0;i<nchar;i++){
		if(*Ldata > -1 && *Rdata > -1){
			for(int i=0;i<4;i++){
				*(dest++) = Lpr[(*Ldata)+16*i] * Rpr[(*Rdata)+16*i];
				*(dest++) = Lpr[(*Ldata+4)+16*i] * Rpr[(*Rdata+4)+16*i];
				*(dest++) = Lpr[(*Ldata+8)+16*i] * Rpr[(*Rdata+8)+16*i];
				*(dest++) = Lpr[(*Ldata+12)+16*i] * Rpr[(*Rdata+12)+16*i];
//				assert(*(dest-1)>0.0);
				}
			Ldata++;
			Rdata++;
			}
			
		else if((*Ldata == -4 && *Rdata == -4) || (*Ldata == -4 && *Rdata > -1) || (*Rdata == -4 && *Ldata > -1)){//total ambiguity of left, right or both
			
			if(*Ldata == -4 && *Rdata == -4) //total ambiguity of both
				for(int i=0;i<16;i++) *(dest++) = 1.0;
			
			else if(*Ldata == -4){//total ambiguity of left
				for(int i=0;i<4;i++){
					*(dest++) = Rpr[(*Rdata)+16*i];
					*(dest++) = Rpr[(*Rdata+4)+16*i];
					*(dest++) = Rpr[(*Rdata+8)+16*i];
					*(dest++) = Rpr[(*Rdata+12)+16*i];
					assert(*(dest-4)>=0.0);
					}
				}
			else{//total ambiguity of right
				for(int i=0;i<4;i++){
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
				for(int i=0;i<4;i++){
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
					for(int i=0;i<16;i++){
						*(dest+i)=1.0;
						}
					Ldata++;
					}
			 
				else{//partially ambiguous left
					int nstates=-*(Ldata++);
					for(int q=0;q<16;q++) dest[q]=0;
					for(int i=0;i<nstates;i++){
						for(int i=0;i<4;i++){
							*(dest+(i*4)) += Lpr[(*Ldata)+16*i];
							*(dest+(i*4)+1) += Lpr[(*Ldata+4)+16*i];
							*(dest+(i*4)+2) += Lpr[(*Ldata+8)+16*i];
							*(dest+(i*4)+3) += Lpr[(*Ldata+12)+16*i];
//							assert(*(dest-1)>0.0);
							}
						Ldata++;
						}
					}
				}
			if(*Rdata>-1){//unambiguous right
				for(int i=0;i<4;i++){
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
				double tempcla[16];
				for(int q=0;q<16;q++) tempcla[q]=0;
				for(int i=0;i<nstates;i++){
					for(int i=0;i<4;i++){
						tempcla[(i*4)]   += Rpr[(*Rdata)+16*i];
						tempcla[(i*4)+1] += Rpr[(*Rdata+4)+16*i];
						tempcla[(i*4)+2] += Rpr[(*Rdata+8)+16*i];
						tempcla[(i*4)+3] += Rpr[(*Rdata+12)+16*i];
//						assert(*(dest-1)>0.0);
						}
					Rdata++;
					}
				//Now multiply the temporary results against the already calced left
				for(int i=0;i<4;i++){
					*(dest++) *= tempcla[(i*4)];
					*(dest++) *= tempcla[(i*4)+1];
					*(dest++) *= tempcla[(i*4)+2];
					*(dest++) *= tempcla[(i*4)+3];
//					assert(*(dest-1)>0.0);
					}		
				}
			else{//fully ambiguous right
				dest+=16;
				Rdata++;
				}
			}
		}
		
		for(int site=0;site<nchar;site++){
			destCLA->underflow_mult[site]=0;
			}
		destCLA->rescaleRank=2;
	}

void CalcFullCLATerminalTerminal(CondLikeArray *destCLA, const double *Lpr, const double *Rpr, char *Ldata, char *Rdata, int nchar){
	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.
	double *dest=destCLA->arr;
	
#ifdef UNIX
	madvise(dest, nchar*4*sizeof(double), MADV_SEQUENTIAL);
#endif

	for(int i=0;i<nchar;i++){
		if(*Ldata > -1 && *Rdata > -1){
			*(dest++) = Lpr[(*Ldata)] * Rpr[(*Rdata)];
			*(dest++) = Lpr[(*Ldata+4)] * Rpr[(*Rdata+4)];
			*(dest++) = Lpr[(*Ldata+8)] * Rpr[(*Rdata+8)];
			*(dest++) = Lpr[(*Ldata+12)] * Rpr[(*Rdata+12)];
			Ldata++;
			Rdata++;
			}
			
		else if((*Ldata == -4 && *Rdata == -4) || (*Ldata == -4 && *Rdata > -1) || (*Rdata == -4 && *Ldata > -1)){//total ambiguity of left, right or both
			
			if(*Ldata == -4 && *Rdata == -4) //total ambiguity of both
				for(int i=0;i<4;i++) *(dest++) = 1.0;
			
			else if(*Ldata == -4){//total ambiguity of left
				*(dest++) = Rpr[(*Rdata)];
				*(dest++) = Rpr[(*Rdata+4)];
				*(dest++) = Rpr[(*Rdata+8)];
				*(dest++) = Rpr[(*Rdata+12)];
				}
			else{//total ambiguity of right
				*(dest++) = Lpr[(*Ldata)];
				*(dest++) = Lpr[(*Ldata+4)];
				*(dest++) = Lpr[(*Ldata+8)];
				*(dest++) = Lpr[(*Ldata+12)];
				}
			Ldata++;
			Rdata++;
			}	
		else {//partial ambiguity of left, right or both
			if(*Ldata>-1){//unambiguous left
				*(dest) = Lpr[(*Ldata)];
				*(dest+1) = Lpr[(*Ldata+4)];
				*(dest+2) = Lpr[(*Ldata+8)];
				*(dest+3) = Lpr[(*Ldata+12)];
				Ldata++;
				}
			else{
				if(*Ldata==-4){//fully ambiguous left
					for(int i=0;i<4;i++){
						*(dest+i)=1.0;
						}
					Ldata++;
					}
			 
				else{//partially ambiguous left
					int nstates=-*(Ldata++);
					for(int q=0;q<4;q++) dest[q]=0;
					for(int i=0;i<nstates;i++){
						*(dest) += Lpr[(*Ldata)];
						*(dest+1) += Lpr[(*Ldata+4)];
						*(dest+2) += Lpr[(*Ldata+8)];
						*(dest+3) += Lpr[(*Ldata+12)];
						Ldata++;
						}
					}
				}
			if(*Rdata>-1){//unambiguous right
				*(dest++) *= Rpr[(*Rdata)];
				*(dest++) *= Rpr[(*Rdata+4)];
				*(dest++) *= Rpr[(*Rdata+8)];
				*(dest++) *= Rpr[(*Rdata+12)];
				Rdata++;		
				}
			else if(*Rdata != -4){//partially ambiguous right
				char nstates=-1 * *(Rdata++);
				//create a temporary cla to hold the results from the ambiguity of the right, 
				//which need to be +'s 
				double tempcla[4];
				for(int q=0;q<4;q++) tempcla[q]=0;
				for(int i=0;i<nstates;i++){
					tempcla[0]   += Rpr[(*Rdata)];
					tempcla[1] += Rpr[(*Rdata+4)];
					tempcla[2] += Rpr[(*Rdata+8)];
					tempcla[3] += Rpr[(*Rdata+12)];
					Rdata++;
					}
				//Now multiply the temporary results against the already calced left
				*(dest++) *= tempcla[0];
				*(dest++) *= tempcla[1];
				*(dest++) *= tempcla[2];
				*(dest++) *= tempcla[3];

				}
			else{//fully ambiguous right
				dest+=4;
				Rdata++;
				}
			}
			}
	
		for(int site=0;site<nchar;site++){
			destCLA->underflow_mult[site]=0;
			}
		destCLA->rescaleRank=2;
	}

void CalcFullCLAInternalTerminalRateHet(CondLikeArray *destCLA, const CondLikeArray *LCLA, const double *pr1, const double *pr2, char *data2, int nchar){
	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.
	double *dest=destCLA->arr;
	double *CL1=LCLA->arr;
	
#ifdef UNIX	
	madvise(dest, nchar*16*sizeof(double), MADV_SEQUENTIAL);
	madvise((void*)CL1, nchar*16*sizeof(double), MADV_SEQUENTIAL);	
#endif

	for(int i=0;i<nchar;i++){
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
//			assert(*(dest)>0.0);
			dest+=16;
			data2++;
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
//			assert(*(dest)>0.0);
			dest+=16;
			data2++;
			}
		else {//partial ambiguity
			//first figure in the ambiguous terminal
			char nstates=-1 * *(data2++);
			for(int q=0;q<16;q++) dest[q]=0;
			for(int i=0;i<nstates;i++){
				for(int i=0;i<4;i++){
					*(dest+(i*4)) += pr2[(*data2)+16*i];
					*(dest+(i*4)+1) += pr2[(*data2+4)+16*i];
					*(dest+(i*4)+2) += pr2[(*data2+8)+16*i];
					*(dest+(i*4)+3) += pr2[(*data2+12)+16*i];
//					assert(*(dest)>0.0);
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
//			assert(*(dest-1)>0.0);
			}
			
		CL1+=16;
		}
	for(int i=0;i<nchar;i++)
		destCLA->underflow_mult[i]=LCLA->underflow_mult[i];
	
	destCLA->rescaleRank=LCLA->rescaleRank+2; 
		
	} 

void CalcFullCLAInternalTerminal(CondLikeArray *destCLA, const CondLikeArray *LCLA, const double *pr1, const double *pr2, char *data2, int nchar){
	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.
	double *dest=destCLA->arr;
	double *CL1=LCLA->arr;
	
#ifdef UNIX	
	madvise(dest, nchar*4*sizeof(double), MADV_SEQUENTIAL);
	madvise((void*)CL1, nchar*4*sizeof(double), MADV_SEQUENTIAL);	
#endif

	for(int i=0;i<nchar;i++){
		if(*data2 > -1){ //no ambiguity
			*(dest++) = ( pr1[0]*CL1[0]+pr1[1]*CL1[1]+pr1[2]*CL1[2]+pr1[3]*CL1[3]) * pr2[(*data2)];
			*(dest++) = ( pr1[4]*CL1[0]+pr1[5]*CL1[1]+pr1[6]*CL1[2]+pr1[7]*CL1[3]) * pr2[(*data2+4)];
			*(dest++) = ( pr1[8]*CL1[0]+pr1[9]*CL1[1]+pr1[10]*CL1[2]+pr1[11]*CL1[3]) * pr2[(*data2+8)];
			*(dest++) = ( pr1[12]*CL1[0]+pr1[13]*CL1[1]+pr1[14]*CL1[2]+pr1[15]*CL1[3]) * pr2[(*data2+12)];
			data2++;
			}
		else if(*data2 == -4){//total ambiguity
			*(dest++) = ( pr1[0]*CL1[0]+pr1[1]*CL1[1]+pr1[2]*CL1[2]+pr1[3]*CL1[3]);
			*(dest++) = ( pr1[4]*CL1[0]+pr1[5]*CL1[1]+pr1[6]*CL1[2]+pr1[7]*CL1[3]);
			*(dest++) = ( pr1[8]*CL1[0]+pr1[9]*CL1[1]+pr1[10]*CL1[2]+pr1[11]*CL1[3]);
			*(dest++) = ( pr1[12]*CL1[0]+pr1[13]*CL1[1]+pr1[14]*CL1[2]+pr1[15]*CL1[3]);
			data2++;
			}
		else {//partial ambiguity
			//first figure in the ambiguous terminal
			char nstates=-1 * *(data2++);
			for(int q=0;q<4;q++) dest[q]=0;
			for(int i=0;i<nstates;i++){
				*(dest) += pr2[(*data2)];
				*(dest+1) += pr2[(*data2+4)];
				*(dest+2) += pr2[(*data2+8)];
				*(dest+3) += pr2[(*data2+12)];
				data2++;
				}
			
			//now add the internal child
			*(dest++) *= ( pr1[0]*CL1[0]+pr1[1]*CL1[1]+pr1[2]*CL1[2]+pr1[3]*CL1[3]);
			*(dest++) *= ( pr1[4]*CL1[0]+pr1[5]*CL1[1]+pr1[6]*CL1[2]+pr1[7]*CL1[3]);
			*(dest++) *= ( pr1[8]*CL1[0]+pr1[9]*CL1[1]+pr1[10]*CL1[2]+pr1[11]*CL1[3]);
			*(dest++) *= ( pr1[12]*CL1[0]+pr1[13]*CL1[1]+pr1[14]*CL1[2]+pr1[15]*CL1[3]);
			}
		CL1+=4;
		}
	for(int i=0;i<nchar;i++)
		destCLA->underflow_mult[i]=LCLA->underflow_mult[i];
	
	destCLA->rescaleRank=LCLA->rescaleRank+2; 
		
	} 

void CalcFullCLAPartialInternalRateHet(CondLikeArray *destCLA, const CondLikeArray *LCLA, const double *pr1, CondLikeArray *partialCLA, int nchar){
	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.
	double *dest=destCLA->arr;
	double *CL1=LCLA->arr;
	double *partial=partialCLA->arr;
	
#ifdef UNIX
	madvise(dest, nchar*16*sizeof(double), MADV_SEQUENTIAL);
	madvise((void*)CL1, nchar*16*sizeof(double), MADV_SEQUENTIAL);
	madvise(partial, nchar*16*sizeof(double), MADV_SEQUENTIAL);
#endif
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
//		assert(*(dest-1)>0.0);
		}
	for(int site=0;site<nchar;site++){
		destCLA->underflow_mult[site]=partialCLA->underflow_mult[site] + LCLA->underflow_mult[site];
		}
	}

void CalcFullCLAPartialInternal(CondLikeArray *destCLA, const CondLikeArray *LCLA, const double *pr1, CondLikeArray *partialCLA, int nchar){
	double *dest=destCLA->arr;
	double *CL1=LCLA->arr;
	double *partial=partialCLA->arr;
	
#ifdef UNIX
	madvise(dest, nchar*4*sizeof(double), MADV_SEQUENTIAL);
	madvise((void*)CL1, nchar*4*sizeof(double), MADV_SEQUENTIAL);
	madvise(partial, nchar*4*sizeof(double), MADV_SEQUENTIAL);
#endif
	for(int i=0;i<nchar;i++){
		*(dest++) = ( pr1[0]*CL1[0]+pr1[1]*CL1[1]+pr1[2]*CL1[2]+pr1[3]*CL1[3]) * *(partial++);
		*(dest++) = ( pr1[4]*CL1[0]+pr1[5]*CL1[1]+pr1[6]*CL1[2]+pr1[7]*CL1[3]) * *(partial++);
		*(dest++) = ( pr1[8]*CL1[0]+pr1[9]*CL1[1]+pr1[10]*CL1[2]+pr1[11]*CL1[3]) * *(partial++);
		*(dest++) = ( pr1[12]*CL1[0]+pr1[13]*CL1[1]+pr1[14]*CL1[2]+pr1[15]*CL1[3]) * *(partial++);
		CL1+=4;
		}
	for(int site=0;site<nchar;site++){
		destCLA->underflow_mult[site]=partialCLA->underflow_mult[site] + LCLA->underflow_mult[site];
		}
	}

void CalcFullCLAPartialTerminalRateHet(CondLikeArray *destCLA, const CondLikeArray *partialCLA, const double *Lpr, char *Ldata, int nchar){
	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.
	double *dest=destCLA->arr;
	double *partial=partialCLA->arr;
#ifdef UNIX
	madvise(dest, nchar*16*sizeof(double), MADV_SEQUENTIAL);
	madvise((void*)partial, nchar*16*sizeof(double), MADV_SEQUENTIAL);
#endif

	for(int i=0;i<nchar;i++){
		if(*Ldata > -1){ //no ambiguity
			for(int i=0;i<4;i++){
				*(dest++) = Lpr[(*Ldata)+16*i] * *(partial++);
				*(dest++) = Lpr[(*Ldata+4)+16*i] * *(partial++);
				*(dest++) = Lpr[(*Ldata+8)+16*i] * *(partial++);
				*(dest++) = Lpr[(*Ldata+12)+16*i] * *(partial++);
//				assert(*(dest-1)>0.0);
				}
			Ldata++;
			}
			
		else if(*Ldata == -4){ //total ambiguity
			for(int i=0;i<16;i++) *(dest++) = *(partial++);
			Ldata++;
			}
		else{ //partial ambiguity
			//first figure in the ambiguous terminal
			char nstates=-1 * *(Ldata++);
			for(int q=0;q<16;q++) dest[q]=0;
			for(int i=0;i<nstates;i++){
				for(int i=0;i<4;i++){
					*(dest+(i*4)) += Lpr[(*Ldata)+16*i];
					*(dest+(i*4)+1) += Lpr[(*Ldata+4)+16*i];
					*(dest+(i*4)+2) += Lpr[(*Ldata+8)+16*i];
					*(dest+(i*4)+3) += Lpr[(*Ldata+12)+16*i];
//					assert(*(dest-1)>0.0);
					}
				Ldata++;
				}
			
			//now add the partial
			*(dest++) *= *(partial++);
			*(dest++) *= *(partial++);
			*(dest++) *= *(partial++);
			*(dest++) *= *(partial++);
			
			*(dest++) *= *(partial++);
			*(dest++) *= *(partial++);
			*(dest++) *= *(partial++);
			*(dest++) *= *(partial++);
		
			*(dest++) *= *(partial++);
			*(dest++) *= *(partial++);
			*(dest++) *= *(partial++);
			*(dest++) *= *(partial++);
		
			*(dest++) *= *(partial++);
			*(dest++) *= *(partial++);
			*(dest++) *= *(partial++);
			*(dest++) *= *(partial++);
//			assert(*(dest-1)>0.0);
			}
	
		}
	for(int i=0;i<nchar;i++)
		destCLA->underflow_mult[i]=partialCLA->underflow_mult[i];
	}

void CalcFullCLAPartialTerminal(CondLikeArray *destCLA, const CondLikeArray *partialCLA, const double *Lpr, char *Ldata, int nchar){
	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.
	double *dest=destCLA->arr;
	double *partial=partialCLA->arr;
#ifdef UNIX
	madvise(dest, nchar*4*sizeof(double), MADV_SEQUENTIAL);
	madvise((void*)partial, nchar*4*sizeof(double), MADV_SEQUENTIAL);
#endif

	for(int i=0;i<nchar;i++){
		if(*Ldata > -1){ //no ambiguity
			*(dest++) = Lpr[(*Ldata)] * *(partial++);
			*(dest++) = Lpr[(*Ldata+4)] * *(partial++);
			*(dest++) = Lpr[(*Ldata+8)] * *(partial++);
			*(dest++) = Lpr[(*Ldata+12)] * *(partial++);
			Ldata++;
	}
			
		else if(*Ldata == -4){ //total ambiguity
			for(int i=0;i<4;i++) *(dest++) = *(partial++);
			Ldata++;
			}
		else{ //partial ambiguity
			//first figure in the ambiguous terminal
			char nstates=-1 * *(Ldata++);
			for(int q=0;q<4;q++) dest[q]=0;
			for(int i=0;i<nstates;i++){
				*(dest) += Lpr[(*Ldata)];
				*(dest+1) += Lpr[(*Ldata+4)];
				*(dest+2) += Lpr[(*Ldata+8)];
				*(dest+3) += Lpr[(*Ldata+12)];
				Ldata++;
				}
			
			//now add the partial
			*(dest++) *= *(partial++);
			*(dest++) *= *(partial++);
			*(dest++) *= *(partial++);
			*(dest++) *= *(partial++);
			
			}
	
		}
	for(int i=0;i<nchar;i++)
		destCLA->underflow_mult[i]=partialCLA->underflow_mult[i];
	}
