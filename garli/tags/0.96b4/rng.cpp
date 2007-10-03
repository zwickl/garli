// GARLI version 0.96b4 source code
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
//	NOTE: Portions of this source adapted from GAML source, written by Paul O. Lewis

#include <math.h>
#include <iostream>
using namespace std;

#include "defs.h"
#include "rng.h"
rng rnd;

rng::rng() : ix0(1L), ix(1L), ifault(0)
{
	randomize();
}

void rng::dememorize( int spin /* = 100 */ )
{
        for( int k = 0; k < spin; k++ )
                uniform();
}

void rng::randomize( int spin /* = 100 */ )
{
	time_t timer;
	ix=ix0=(long) time(&timer);
        dememorize( spin );
}

long rng::random_long(long max)
{
	if (max == 0)	return 0;

	long return_val = max;

	while( return_val == max )
		return_val = (long)( (FLOAT_TYPE)max * uniform() );

	return return_val;
}

int rng::random_int(int max)
{
	if (max == 0)	return 0;
	int return_val = max;

	while( return_val == max )
		return_val = (int)( (FLOAT_TYPE)max * uniform() );

	return return_val;
}

float rng::random_float(float max)	{
	if (max == 0.0)
		return 0.0;
	float return_val = max;
	
	while (return_val == max)
		return_val = (float)( (FLOAT_TYPE)max * uniform() );
		
	return return_val;
}

FLOAT_TYPE rng::random_FLOAT_TYPE(FLOAT_TYPE max)	{
	if (max == 0.0)
		return 0.0;
	FLOAT_TYPE return_val = max;
	
	while (return_val == max)
		return_val = max * uniform();
		
	return return_val;
}

FLOAT_TYPE rng::exponential(FLOAT_TYPE lambda)
{
	FLOAT_TYPE x = 0.0;

	while( x <= 0.0 || x > 1.0 )
		x = ONE_POINT_ZERO - uniform();
	x = -log(x) / lambda;

	return x;
}

FLOAT_TYPE rng::gamln( FLOAT_TYPE x )
{
    // ====================================================================== 
    // NIST Guide to Available Math Software. 
    // Source for module GAMLN from package CMLIB. 
    // Retrieved from TIBER on Wed Apr 29 17:30:20 1998. 
    // ====================================================================== 
    //     WRITTEN BY D. E. AMOS, SEPTEMBER, 1977. 
    //
    //     REFERENCES 
    //         SAND-77-1518 
    //
    //         COMPUTER APPROXIMATIONS BY J.F.HART, ET.AL., SIAM SERIES IN 
    //         APPLIED MATHEMATICS, WILEY, 1968, P.135-136. 
    //
    //         NBS HANDBOOK OF MATHEMATICAL FUNCTIONS, AMS 55, BY 
    //         M. ABRAMOWITZ AND I.A. STEGUN, DECEMBER. 1955, P.257. 
    //
    //     ABSTRACT 
    //         GAMLN COMPUTES THE NATURAL LOG OF THE GAMMA FUNCTION FOR 
    //         X.GT.0. A RATIONAL CHEBYSHEV APPROXIMATION IS USED ON 
    //         8.LT.X.LT.1000., THE ASYMPTOTIC EXPANSION FOR X.GE.1000. AND 
    //         A RATIONAL CHEBYSHEV APPROXIMATION ON 2.LT.X.LT.3. FOR 
    //         0.LT.X.LT.8. AND X NON-INTEGRAL, FORWARD OR BACKWARD 
    //         RECURSION FILLS IN THE INTERVALS  0.LT.X.LT.2 AND 
    //         3.LT.X.LT.8. FOR X=1.,2.,...,100., GAMLN IS SET TO 
    //         NATURAL LOGS OF FACTORIALS. 
    //
    //     DESCRIPTION OF ARGUMENTS 
    //
    //         INPUT 
    //           X      - X.GT.0 
    //
    //         OUTPUT 
    //           GAMLN  - NATURAL LOG OF THE GAMMA FUNCTION AT X 
    //
    //     ERROR CONDITIONS 
    //         IMPROPER INPUT ARGUMENT - A FATAL ERROR 

    static FLOAT_TYPE xlim1 = (FLOAT_TYPE)8.;
    static FLOAT_TYPE xlim2 = (FLOAT_TYPE)1e3;
    static FLOAT_TYPE rtwpil = (FLOAT_TYPE).918938533204673;
    static FLOAT_TYPE p[5] = { (FLOAT_TYPE)7.66345188e-4,(FLOAT_TYPE)-5.9409561052e-4,(FLOAT_TYPE)
	    7.936431104845e-4,(FLOAT_TYPE)-.00277777775657725,(FLOAT_TYPE)
	    .0833333333333169 };
    static FLOAT_TYPE q[2] = { (FLOAT_TYPE)-.00277777777777778,(FLOAT_TYPE).0833333333333333 }
	    ;
    static FLOAT_TYPE pcoe[9] = { (FLOAT_TYPE).00297378664481017,(FLOAT_TYPE)
	    .0092381945590276,(FLOAT_TYPE).109311595671044,(FLOAT_TYPE).398067131020357,
	    (FLOAT_TYPE)2.15994312846059,(FLOAT_TYPE)6.33806799938727,(FLOAT_TYPE)
	    20.7824725317921,(FLOAT_TYPE)36.0367725300248,(FLOAT_TYPE)62.0038380071273 }
	    ;
    static FLOAT_TYPE qcoe[4] = { (FLOAT_TYPE)1.,(FLOAT_TYPE)-8.90601665949746,(FLOAT_TYPE)
	    9.82252110471399,(FLOAT_TYPE)62.003838007127 };
    static FLOAT_TYPE gln[100] = { (FLOAT_TYPE)0.,(FLOAT_TYPE)0.,(FLOAT_TYPE).693147180559945,(
	    FLOAT_TYPE)1.79175946922806,(FLOAT_TYPE)3.17805383034795,(FLOAT_TYPE)
	    4.78749174278205,(FLOAT_TYPE)6.5792512120101,(FLOAT_TYPE)8.52516136106541,(
	    FLOAT_TYPE)10.6046029027453,(FLOAT_TYPE)12.8018274800815,(FLOAT_TYPE)
	    15.1044125730755,(FLOAT_TYPE)17.5023078458739,(FLOAT_TYPE)19.9872144956619,(
	    FLOAT_TYPE)22.5521638531234,(FLOAT_TYPE)25.1912211827387,(FLOAT_TYPE)
	    27.8992713838409,(FLOAT_TYPE)30.6718601060807,(FLOAT_TYPE)33.5050734501369,(
	    FLOAT_TYPE)36.3954452080331,(FLOAT_TYPE)39.3398841871995,(FLOAT_TYPE)
	    42.3356164607535,(FLOAT_TYPE)45.3801388984769,(FLOAT_TYPE)48.4711813518352,(
	    FLOAT_TYPE)51.6066755677644,(FLOAT_TYPE)54.7847293981123,(FLOAT_TYPE)
	    58.0036052229805,(FLOAT_TYPE)61.261701761002,(FLOAT_TYPE)64.5575386270063,(
	    FLOAT_TYPE)67.8897431371815,(FLOAT_TYPE)71.257038967168,(FLOAT_TYPE)
	    74.6582363488302,(FLOAT_TYPE)78.0922235533153,(FLOAT_TYPE)81.557959456115,(
	    FLOAT_TYPE)85.0544670175815,(FLOAT_TYPE)88.5808275421977,(FLOAT_TYPE)
	    92.1361756036871,(FLOAT_TYPE)95.7196945421432,(FLOAT_TYPE)99.3306124547874,(
	    FLOAT_TYPE)102.968198614514,(FLOAT_TYPE)106.631760260643,(FLOAT_TYPE)
	    110.320639714757,(FLOAT_TYPE)114.034211781462,(FLOAT_TYPE)117.771881399745,(
	    FLOAT_TYPE)121.533081515439,(FLOAT_TYPE)125.317271149357,(FLOAT_TYPE)
	    129.123933639127,(FLOAT_TYPE)132.952575035616,(FLOAT_TYPE)136.802722637326,(
	    FLOAT_TYPE)140.673923648234,(FLOAT_TYPE)144.565743946345,(FLOAT_TYPE)
	    148.477766951773,(FLOAT_TYPE)152.409592584497,(FLOAT_TYPE)156.360836303079,(
	    FLOAT_TYPE)160.331128216631,(FLOAT_TYPE)164.320112263195,(FLOAT_TYPE)
	    168.327445448428,(FLOAT_TYPE)172.352797139163,(FLOAT_TYPE)176.395848406997,(
	    FLOAT_TYPE)180.456291417544,(FLOAT_TYPE)184.533828861449,(FLOAT_TYPE)
	    188.628173423672,(FLOAT_TYPE)192.739047287845,(FLOAT_TYPE)196.86618167289,(
	    FLOAT_TYPE)201.009316399282,(FLOAT_TYPE)205.168199482641,(FLOAT_TYPE)
	    209.342586752537,(FLOAT_TYPE)213.532241494563,(FLOAT_TYPE)217.736934113954,(
	    FLOAT_TYPE)221.95644181913,(FLOAT_TYPE)226.190548323728,(FLOAT_TYPE)
	    230.439043565777,(FLOAT_TYPE)234.701723442818,(FLOAT_TYPE)238.978389561834,(
	    FLOAT_TYPE)243.268849002983,(FLOAT_TYPE)247.572914096187,(FLOAT_TYPE)
	    251.890402209723,(FLOAT_TYPE)256.22113555001,(FLOAT_TYPE)260.564940971863,(
	    FLOAT_TYPE)264.921649798553,(FLOAT_TYPE)269.29109765102,(FLOAT_TYPE)
	    273.673124285694,(FLOAT_TYPE)278.067573440366,(FLOAT_TYPE)282.47429268763,(
	    FLOAT_TYPE)286.893133295427,(FLOAT_TYPE)291.32395009427,(FLOAT_TYPE)
	    295.766601350761,(FLOAT_TYPE)300.220948647014,(FLOAT_TYPE)304.686856765669,(
	    FLOAT_TYPE)309.164193580147,(FLOAT_TYPE)313.652829949879,(FLOAT_TYPE)
	    318.152639620209,(FLOAT_TYPE)322.663499126726,(FLOAT_TYPE)327.185287703775,(
	    FLOAT_TYPE)331.717887196928,(FLOAT_TYPE)336.261181979198,(FLOAT_TYPE)
	    340.815058870799,(FLOAT_TYPE)345.379407062267,(FLOAT_TYPE)349.95411804077,(
	    FLOAT_TYPE)354.539085519441,(FLOAT_TYPE)359.134205369575 };

    /* System generated locals */
    long int i__1;
    FLOAT_TYPE ret_val=0.0;

    /* Local variables */
    static FLOAT_TYPE dgam;
    static long int i__;
    static FLOAT_TYPE t, dx, px, qx, rx, xx;
    static long int ndx, nxm;
    static FLOAT_TYPE sum, rxx;

    if ( x <= (FLOAT_TYPE)0.) {
	goto L90;
    } else {
	goto L5;
    }
L5:
    ndx = (long int)x;
    t = x - (FLOAT_TYPE) ndx;
    if (t == (FLOAT_TYPE)0.) {
	goto L51;
    }
    dx = xlim1 - x;
    if (dx < (FLOAT_TYPE)0.) {
	goto L40;
    }

/*     RATIONAL CHEBYSHEV APPROXIMATION ON 2.LT.X.LT.3 FOR GAMMA(X) */

    nxm = ndx - 2;
    px = pcoe[0];
    for (i__ = 2; i__ <= 9; ++i__) {
/* L10: */
	px = t * px + pcoe[i__ - 1];
    }
    qx = qcoe[0];
    for (i__ = 2; i__ <= 4; ++i__) {
/* L15: */
	qx = t * qx + qcoe[i__ - 1];
    }
    dgam = px / qx;
    if (nxm > 0) {
	goto L22;
    }
    if (nxm == 0) {
	goto L25;
    }

/*     BACKWARD RECURSION FOR 0.LT.X.LT.2 */

    dgam /= t + (FLOAT_TYPE)1.;
    if (nxm == -1) {
	goto L25;
    }
    dgam /= t;
    ret_val = log(dgam);
    return ret_val;

/*     FORWARD RECURSION FOR 3.LT.X.LT.8 */

L22:
    xx = t + (FLOAT_TYPE)2.;
    i__1 = nxm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dgam *= xx;
/* L24: */
	xx += (FLOAT_TYPE)1.;
    }
L25:
    ret_val = log(dgam);
    return ret_val;

/*     X.GT.XLIM1 */

L40:
    rx = (FLOAT_TYPE)1. / x;
    rxx = rx * rx;
    if (x - xlim2 < (FLOAT_TYPE)0.) {
	goto L41;
    }
    px = q[0] * rxx + q[1];
    ret_val = px * rx + (x - (FLOAT_TYPE).5) * log(x) - x + rtwpil;
    return ret_val;

/*     X.LT.XLIM2 */

L41:
    px = p[0];
    sum = (x - (FLOAT_TYPE).5) * log(x) - x;
    for (i__ = 2; i__ <= 5; ++i__) {
	px = px * rxx + p[i__ - 1];
/* L42: */
    }
    ret_val = px * rx + sum + rtwpil;
    return ret_val;

/*     TABLE LOOK UP FOR INTEGER ARGUMENTS LESS THAN OR EQUAL 100. */

L51:
    if (ndx > 100) {
	goto L40;
    }
    ret_val = gln[ndx - 1];
    return ret_val;
L90:
    cerr << "GAMLN  ARGUMENT IS LESS THAN OR EQUAL TO ZERO " << endl;
    return ret_val;
}

#if 0
//
// Algorithm 291.  Pike, M. C., and I. D. Hill. 1966. Logarithm of the
// gamma function. Commun. Ass. Comput. Mach. 9: 684.
// (translated to C++ by Paul O. Lewis, November, 1996)
//
// This procedure evaluates the natural logarithm of gamma(x) for all x > 0,
// accurate to 10 decimal places.  Stirling's formula is used for the central
// polynomial part of the procedure
//
FLOAT_TYPE rng::loggamma( FLOAT_TYPE x )
{
        FLOAT_TYPE f = 0.0;
        FLOAT_TYPE z;
        if( x < 7.0 ) {
                f = 1.0;
                for( z = x - 1.0; z < 7.0; z += 1.0 ) {
                        x = z;
                        f *= z;
                }
		x += 1.0;
                f -= log(f);
        }
        z = 1.0 / (x*x);
        FLOAT_TYPE v1 = f + (x - 0.5)*log(x) - x + 0.918938533204673;

        FLOAT_TYPE v2 = -0.000595238095238*z;
        FLOAT_TYPE v3 = ( v2 - 0.000793650793651 )*z;
        FLOAT_TYPE v4 = ( v3 - 0.0002777777777778 )*z;
        FLOAT_TYPE v5 = ( v4 + 0.083333333333333 );
        FLOAT_TYPE v6 = v5 / x;

        return (v1 + v6);
}
#endif

//
// Algorithm AS 32.  Bhattacharjee, G. P. 1970. The incomplete gamma integral.
// Appl. Statist. 19: 285-187.
// (translated to C++ by Paul O. Lewis, November, 1996)
//
// Computes incomplete gamma ratio for positive values of
// arguments x and p.  g must be supplied and should be equal to
// ln( gamma(p) ).
// ifault = 1 if p <= 0. ifault = 2 if x < 0 else ifault = 0
// Uses series expansion if p > x or x <= 1, otherwise a
// continued fraction approximation.
//
FLOAT_TYPE rng::gamain( FLOAT_TYPE x, FLOAT_TYPE p, FLOAT_TYPE g )
{
	FLOAT_TYPE pn[6];

	// define accuracy and initialize
	FLOAT_TYPE acu = (FLOAT_TYPE)1.0e-8;
	const FLOAT_TYPE oflo = (FLOAT_TYPE)1.0e30;
	FLOAT_TYPE gin = 0.0;
	FLOAT_TYPE term;
	FLOAT_TYPE rn;
	ifault = 0;

        // test for admissibility of arguments
        if( p <= 0.0 ) ifault = 1;
        if( x < 0.0 ) ifault = 2;
        if( ifault > 0 || x == 0.0 )
                return gin;
        FLOAT_TYPE factor = exp( p * log(x) - x - g );
        if( x > 1.0 && x >= p )
                goto label30;

        // calculation by series expansion
        gin = 1.0;
        term = 1.0;
        rn = p;

        for(;;) {
                rn += 1.0;
                term *= (x / rn);
                gin += term;
                if( term <= acu ) break;
        }
        gin *= ( factor / p );
        return gin;

        // calculation by continued fraction
        FLOAT_TYPE a, b, an, dif;
        int i;
        label30:        // <<<<<<<<<<<<<<< label 30
        a = ONE_POINT_ZERO - p;
        b = a + x + ONE_POINT_ZERO;
        term = 0.0;
        pn[0] = 1.0;
        pn[1] = x;
        pn[2] = x + ONE_POINT_ZERO;
        pn[3] = x*b;
        gin = pn[2] / pn[3];

        label32:        // <<<<<<<<<<<<<<< label 32
        a += 1.0;
        b += 2.0;
        term += 1.0;
        an = a * term;
        pn[4] = b*pn[2] - an*pn[0];
        pn[5] = b*pn[3] - an*pn[1];
        if( pn[5] == 0.0 )
                goto label35;
        rn = pn[4] / pn[5];
        dif = fabs( gin - rn );
        if( dif > acu )
                goto label34;
        if( dif <= acu*rn ) {
                gin = ONE_POINT_ZERO - factor*gin;
                return gin;
        }

        label34:        // <<<<<<<<<<<<<<< label 34
        gin = rn;

        label35:        // <<<<<<<<<<<<<<< label 35
        for( i = 0; i < 4; i++ )
                pn[i] = pn[i+2];
        if( fabs(pn[4]) < oflo )
                goto label32;
        for( i = 0; i < 4; i++ )
                pn[i] /= oflo;
        goto label32;
}

//
// Algorithm AS 70. Odeh, R. E., and Evans, J. O. 1974. Percentage points
// of the normal distribution. Appl. Statist. 23: 96-97.
// (translated to C++ by Paul O. Lewis, November, 1996)
//
// gauinv finds percentage points of the normal distribution.
//
FLOAT_TYPE rng::gauinv( FLOAT_TYPE p )
{
#ifdef SINGLE_PRECISION_FLOATS
	const FLOAT_TYPE zero = 0.0f;
	const FLOAT_TYPE one = 1.0f;
	const FLOAT_TYPE half = 0.5f;
	const FLOAT_TYPE alimit = 1.0e-20f;
	const FLOAT_TYPE p0 = -0.322232431088f;
	const FLOAT_TYPE p1 = -1.0f;
	const FLOAT_TYPE p2 = -0.342242088547f;
	const FLOAT_TYPE p3 = -0.0204231210245f;
	const FLOAT_TYPE p4 = -0.0000453642210148f;
	const FLOAT_TYPE q0 = 0.099348462606f;
	const FLOAT_TYPE q1 = 0.58858157495f;
	const FLOAT_TYPE q2 = 0.531103462366f;
	const FLOAT_TYPE q3 = 0.10353775285f;
	const FLOAT_TYPE q4 = 0.0038560700634f;
#else
	const FLOAT_TYPE zero = 0.0;
	const FLOAT_TYPE one = 1.0;
	const FLOAT_TYPE half = 0.5;
	const FLOAT_TYPE alimit = 1.0e-20;
	const FLOAT_TYPE p0 = -0.322232431088;
	const FLOAT_TYPE p1 = -1.0;
	const FLOAT_TYPE p2 = -0.342242088547;
	const FLOAT_TYPE p3 = -0.0204231210245;
	const FLOAT_TYPE p4 = -0.0000453642210148;
	const FLOAT_TYPE q0 = 0.099348462606;
	const FLOAT_TYPE q1 = 0.58858157495;
	const FLOAT_TYPE q2 = 0.531103462366;
	const FLOAT_TYPE q3 = 0.10353775285;
	const FLOAT_TYPE q4 = 0.0038560700634;
#endif

        ifault = 1;
        FLOAT_TYPE ps = p;
        if( ps > half )
                ps = one - ps;
        if( ps < alimit )
                return zero;

        ifault = 0;
        if( ps == half )
                return zero;
        FLOAT_TYPE yi = sqrt( log( one / (ps*ps) ) );
        FLOAT_TYPE retval = yi
                + ((((yi*p4 + p3)*yi + p2)*yi + p1)*yi + p0)
                / ((((yi*q4 + q3)*yi + q2)*yi + q1)*yi + q0);
        if( p < half )
                retval = -retval;

        return retval;
}
#define A_A 16807L
#define	B_B15 32768L
#define		B_B16 65536L
#define		P_P 2147483647L
	
FLOAT_TYPE rng::uniform()
{
	//long a, p, b15, b16, xhi, xalo, leftlo, fhi, k;
	long xhi, xalo, leftlo, fhi, k;
	//
	// Uniform pseudorandom number generator
	// Provided by J. Monahan, Statistics Dept., N.C. State University
	//   From Schrage, ACM Trans. Math. Software 5:132-138 (1979)
	// Translated to C by Paul O. Lewis, Dec. 10, 1992
	//
	xhi = ix / B_B16;
	xalo = (ix - xhi * B_B16) * A_A;
	leftlo = xalo / B_B16;
	fhi = xhi * A_A + leftlo;
	k = fhi / B_B15;
	ix = (((xalo - leftlo * B_B16) - P_P) + (fhi - k * B_B15) * B_B16) + k;
	if (ix < 0) ix += P_P;
	return ix * (FLOAT_TYPE)4.6566128575e-10;
}

//
// Algorithm AS 91.  Best, D. J., and D. E. Roberts. 1975. The percentage
// points of the chi-square distribution. Appl. Statist. 24(3): 385-388.
// (translated to C++ by Paul O. Lewis, November, 1996)
//
// To evaluate the percentage points of the chi-squared
// probability distribuiton function.
//   p must lie in the range 0.000002 to 0.999998
//   v must be positive
//   g must be supplied and should be equal to ln(gamma(v/2.0))
//
//   ifault values:
//      0: everything went fine
//      1: p was not in range 0.000002 to 0.999998
//      2: v was not positive
//      3:
//
FLOAT_TYPE rng::ppchi2( FLOAT_TYPE p, FLOAT_TYPE v )
{
#ifdef SINGLE_PRECISION_FLOATS
	const FLOAT_TYPE e = 0.5e-4f;
	const FLOAT_TYPE aa = 0.69314718f;
#else
	const FLOAT_TYPE e = 0.5e-6;
	const FLOAT_TYPE aa = 0.6931471805;
#endif

	FLOAT_TYPE ch, a, b, q, p1, p2, t, x;
    FLOAT_TYPE s1, s2, s3, s4, s5, s6;

#if 0
	FLOAT_TYPE g = gammln( v / 2.0 );
#else
	FLOAT_TYPE g = gamln( v / (FLOAT_TYPE)2.0 );
#endif

	// after defining accuracy and ln(2), test arguments and initialize
	ifault = 1;
	if( p < 0.000002 || p > 0.999998 ) return -1.0;

	ifault = 2;
	if( v <= 0.0 ) return -ONE_POINT_ZERO;

	ifault = 0;
	FLOAT_TYPE xx = (FLOAT_TYPE)0.5 * v;
	FLOAT_TYPE c = xx - ONE_POINT_ZERO;

        // starting approximation for small chi-squared
        if( v >= -1.24 * log(p) )
                goto label1;
        ch = pow( (p * xx * exp( g + xx * aa)), (ONE_POINT_ZERO / xx) );
        if( ch - e < 0.0 )
                return ch;
        else
                goto label4;

        // starting approximation for v less than or equal to 0.32
	label1:
        if( v > 0.32 ) goto label3;
        ch = (FLOAT_TYPE)0.4;
        a = log( ONE_POINT_ZERO - p );

        label2:
        q = ch;
        p1 = ONE_POINT_ZERO + ch * ( (FLOAT_TYPE)4.67 + ch );
        p2 = ch * ( (FLOAT_TYPE)6.73 + ch * ( (FLOAT_TYPE)6.66 + ch ));
        t = (FLOAT_TYPE)-0.5 + ( (FLOAT_TYPE)4.67 + (FLOAT_TYPE)2.0*ch ) / p1 -
                ( (FLOAT_TYPE)6.73 + ch*( (FLOAT_TYPE)13.32 + (FLOAT_TYPE)3.0*ch )) / p2;
		ch -= ( ONE_POINT_ZERO - exp( a + g + (FLOAT_TYPE)0.5*ch + c*aa ) * p2 / p1 ) / t;
        if( fabs( q/ch - 1.0 ) - 0.01 <= 0.0 )
                goto label4;
        else
		goto label2;

        // call to algorithm AS 70 - note that p has been tested above
        label3:
	x = gauinv( p );

	// starting approximation using Wilson and Hilferty estimate
        p1 = (FLOAT_TYPE) 0.222222 / v;
        ch = v * pow( x*sqrt(p1) + ONE_POINT_ZERO - p1, 3 );

        // starting approximation for p tending to 1
        if( ch > 2.2 * v + 6.0 )
                ch = (FLOAT_TYPE)-2.0 * ( log( ONE_POINT_ZERO - p ) - c * log((FLOAT_TYPE) 0.5 * ch ) + g );

        // call to algorithm AS 32 and calculation of seven term Taylor series
	label4:
        q = ch;
        p1 = (FLOAT_TYPE)0.5 * ch;
	p2 = p - gamain( p1, xx, g );

	if( ifault != 0 ) {
		ifault = 3;
                return -ONE_POINT_ZERO;
	}

        t = p2 * exp( xx*aa + g + p1 - c * log(ch) );
        b = t / ch;
#ifdef SINGLE_PRECISION_FLOATS
        a = 0.5f*t - b*c;
        s1 = ( 210.0f + a*( 140.0f + a*( 105.0f + a*(84.0f + a*( 70.0f + a*60.0f ))))) / 420.0f;
        s2 = ( 420.0f + a*( 735.0f + a*( 966.0f + a*( 1141.0f + a*1278.0f)))) / 2520.0f;
		s3 = ( 210.0f + a*( 462.0f + a*( 707.0f + a*932.0f))) / 2520.0f;
        s4 = ( 252.0f + a*( 672.0f + a*1182.0f) + c*( 204.0f + a*( 889.0f + a*1740.0f))) / 5040.0f;
        s5 = (  84.0f + a*264.0f + c*(175.0f + a*606.0f)) / 2520.0f;
        s6 = ( 120.0f + c*( 346.0f + c*127.0f)) / 5040.0f;
        ch += t*( 1.0f + 0.5f*t*s1 - b*c*( s1 - b*( s2 - b*( s3 - b*( s4 - b*( s5 - b*s6))))));
#else
        a = 0.5*t - b*c;
        s1 = ( 210.0 + a*( 140.0 + a*( 105.0 + a*(84.0 + a*( 70.0 + a*60.0 ))))) / 420.0;
        s2 = ( 420.0 + a*( 735.0 + a*( 966.0 + a*( 1141.0 + a*1278.0)))) / 2520.0;
		s3 = ( 210.0 + a*( 462.0 + a*( 707.0 + a*932.0))) / 2520.0;
        s4 = ( 252.0 + a*( 672.0 + a*1182.0) + c*( 204.0 + a*( 889.0 + a*1740.0))) / 5040.0;
        s5 = (  84.0 + a*264.0 + c*(175.0 + a*606.0)) / 2520.0;
        s6 = ( 120.0 + c*( 346.0 + c*127.0)) / 5040.0;
        ch += t*( 1.0 + 0.5*t*s1 - b*c*( s1 - b*( s2 - b*( s3 - b*( s4 - b*( s5 - b*s6))))));
#endif
		
//		cout << q << "\t" << ch  << "\t" << fabs( q/ch - ONE_POINT_ZERO ) << endl;
        if( fabs( q/ch - ONE_POINT_ZERO ) > e )
                goto label4;
        return ch;
}

//this is from MB
void rng::DirichletRandomVariable (FLOAT_TYPE *alp, FLOAT_TYPE *z, int n){
	int		i;
	FLOAT_TYPE	sum;

	sum = 0.0;
	for(i=0; i<n; i++)
		{
		z[i]=gamma( alp[i] );
//		z[i] = RndGamma (alp[i]) / 1.0;
		sum += z[i];
		}
	for(i=0; i<n; i++)
		z[i] /= sum;
	}


