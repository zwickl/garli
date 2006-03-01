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

//
// Some portions of the source code file for which this is
// the header were written by others and thus no copyright
// is claimed for those portions.  Methods derived from
// other works include gamln, uniform, gamain, gauinv,
// and ppchi2.
//

#ifndef __RNG_H
#define __RNG_H

#ifndef __TIME_H
#	include <time.h>
#endif
#include <cmath>
#include <cassert>
#include <iostream>

using namespace std;

class rng
{
	protected:
		long ix0, ix;
		int ifault;

	protected:
		//double loggamma( double x );
		double gamln( double x );
		double gauinv( double p );
		double gamain( double x, double p, double g );
		double ppchi2( double p, double v );

	public:
		rng();

		long seed()		{ return ix;  }
		long init_seed()	{ return ix0; }

		void randomize( int spin = 100 );
		void set_seed(long s) { ix = ix0 = s; }
		void dememorize( int spin = 100 );

		int		random_int(int);
		long	random_long(long);
		float	random_float(float);
		double	random_double(double);
		double  uniform();
		double  exponential(double);
		double  gamma( double shape ){
			double g=-1;
			do{
				g = ( ppchi2( uniform(), 2.0*shape ) / (2.0*shape) );
				}while( g < 0.0);
			assert(g > 0.0);
			return g;
			}
		//DZ 11-3-02 addition
		int random_binomial(int n, double p);
		void DirichletRandomVariable (double *alp, double *z, int n);
};

//DJZ 11-3-02 Added by me.  Stolen from ProbabLib 1.0, by Antonio Larrosa
//DJZ 3-29-04 Altering this to have the distribution mean specified, which 
//should make picking a value for different datasets a bit easier

inline int rng::random_binomial(int n, double mean){
	double p=mean/n;
	double t=p/(1.0-p);
	double u=uniform();
	double p0=pow((1.0-p),n);
	double g=p0;
	unsigned int k=0;
	while (u>g){
	    p0*=t*(n-k)/(k+1.0);
	    g+=p0;
	    k++;
	    }
	return k;
	}


/*inline int rng::random_binomial(int n, double p){
	double t=p/(1.0-p);
	double u=uniform();
	double p0=pow((1.0-p),n);
	double g=p0;
	unsigned int k=0;
	while (u>g){
	    p0*=t*(n-k)/(k+1.0);
	    g+=p0;
	    k++;
	    }
	return k;
	}
*/

#endif


