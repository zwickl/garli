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

#include "defs.h"
using namespace std;

class rng
{
	protected:
		long ix0, ix;
		int ifault;

	protected:
		//FLOAT_TYPE loggamma( FLOAT_TYPE x );
		FLOAT_TYPE gamln( FLOAT_TYPE x );
		FLOAT_TYPE gauinv( FLOAT_TYPE p );
		FLOAT_TYPE gamain( FLOAT_TYPE x, FLOAT_TYPE p, FLOAT_TYPE g );
		FLOAT_TYPE ppchi2( FLOAT_TYPE p, FLOAT_TYPE v );

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
		FLOAT_TYPE	random_FLOAT_TYPE(FLOAT_TYPE);
		FLOAT_TYPE  uniform();
		FLOAT_TYPE  exponential(FLOAT_TYPE);
		FLOAT_TYPE  gamma( FLOAT_TYPE shape ){
			FLOAT_TYPE g=-1;
			do{
				g = (FLOAT_TYPE)( ppchi2( uniform(), (FLOAT_TYPE)2.0*shape ) / (FLOAT_TYPE)(2.0*shape) );
				}while( g < 0.0);
			assert(g > 0.0);
			return g;
			}
		FLOAT_TYPE gamma_two_param(FLOAT_TYPE alpha, FLOAT_TYPE beta){
			FLOAT_TYPE g=-1;
			do{
				g = (FLOAT_TYPE)( ppchi2( uniform(), (FLOAT_TYPE)2.0*alpha ) / (FLOAT_TYPE)(2.0*beta) );
				}while( g < 0.0);
			assert(g > 0.0);
			return g;
			}
		//DZ 11-3-02 addition
		int random_binomial(int n, FLOAT_TYPE p);
		void DirichletRandomVariable (FLOAT_TYPE *alp, FLOAT_TYPE *z, int n);
};

//DJZ 11-3-02 Added by me.  Stolen from ProbabLib 1.0, by Antonio Larrosa
//DJZ 3-29-04 Altering this to have the distribution mean specified, which 
//should make picking a value for different datasets a bit easier

inline int rng::random_binomial(int n, FLOAT_TYPE mean){
	FLOAT_TYPE p=mean/n;
	FLOAT_TYPE t=(FLOAT_TYPE) p/((FLOAT_TYPE)1.0-p);
	FLOAT_TYPE u=uniform();
	FLOAT_TYPE p0=pow((FLOAT_TYPE)((FLOAT_TYPE)1.0-p),n);
	FLOAT_TYPE g=p0;
	unsigned int k=0;
	while (u>g){
	    p0*=t*(n-k)/(FLOAT_TYPE)(k+1.0);
	    g+=p0;
	    k++;
	    }
	return k;
	}


/*inline int rng::random_binomial(int n, FLOAT_TYPE p){
	FLOAT_TYPE t=p/(1.0-p);
	FLOAT_TYPE u=uniform();
	FLOAT_TYPE p0=pow((1.0-p),n);
	FLOAT_TYPE g=p0;
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


