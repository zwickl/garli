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

#include <math.h>
#include <iostream>
using namespace std;

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
		return_val = (long)( (double)max * uniform() );

	return return_val;
}

int rng::random_int(int max)
{
	if (max == 0)	return 0;
	int return_val = max;

	while( return_val == max )
		return_val = (int)( (double)max * uniform() );

	return return_val;
}

float rng::random_float(float max)	{
	if (max == 0.0)
		return 0.0;
	float return_val = max;
	
	while (return_val == max)
		return_val = (float)( (double)max * uniform() );
		
	return return_val;
}

double rng::random_double(double max)	{
	if (max == 0.0)
		return 0.0;
	double return_val = max;
	
	while (return_val == max)
		return_val = max * uniform();
		
	return return_val;
}

double rng::exponential(double lambda)
{
	double x = 0.0;

	while( x <= 0.0 || x > 1.0 )
		x = 1.0 - uniform();
	x = -log(x) / lambda;

	return x;
}

double rng::gamln( double x )
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

    static double xlim1 = (double)8.;
    static double xlim2 = (double)1e3;
    static double rtwpil = (double).918938533204673;
    static double p[5] = { (double)7.66345188e-4,(double)-5.9409561052e-4,(double)
	    7.936431104845e-4,(double)-.00277777775657725,(double)
	    .0833333333333169 };
    static double q[2] = { (double)-.00277777777777778,(double).0833333333333333 }
	    ;
    static double pcoe[9] = { (double).00297378664481017,(double)
	    .0092381945590276,(double).109311595671044,(double).398067131020357,
	    (double)2.15994312846059,(double)6.33806799938727,(double)
	    20.7824725317921,(double)36.0367725300248,(double)62.0038380071273 }
	    ;
    static double qcoe[4] = { (double)1.,(double)-8.90601665949746,(double)
	    9.82252110471399,(double)62.003838007127 };
    static double gln[100] = { (double)0.,(double)0.,(double).693147180559945,(
	    double)1.79175946922806,(double)3.17805383034795,(double)
	    4.78749174278205,(double)6.5792512120101,(double)8.52516136106541,(
	    double)10.6046029027453,(double)12.8018274800815,(double)
	    15.1044125730755,(double)17.5023078458739,(double)19.9872144956619,(
	    double)22.5521638531234,(double)25.1912211827387,(double)
	    27.8992713838409,(double)30.6718601060807,(double)33.5050734501369,(
	    double)36.3954452080331,(double)39.3398841871995,(double)
	    42.3356164607535,(double)45.3801388984769,(double)48.4711813518352,(
	    double)51.6066755677644,(double)54.7847293981123,(double)
	    58.0036052229805,(double)61.261701761002,(double)64.5575386270063,(
	    double)67.8897431371815,(double)71.257038967168,(double)
	    74.6582363488302,(double)78.0922235533153,(double)81.557959456115,(
	    double)85.0544670175815,(double)88.5808275421977,(double)
	    92.1361756036871,(double)95.7196945421432,(double)99.3306124547874,(
	    double)102.968198614514,(double)106.631760260643,(double)
	    110.320639714757,(double)114.034211781462,(double)117.771881399745,(
	    double)121.533081515439,(double)125.317271149357,(double)
	    129.123933639127,(double)132.952575035616,(double)136.802722637326,(
	    double)140.673923648234,(double)144.565743946345,(double)
	    148.477766951773,(double)152.409592584497,(double)156.360836303079,(
	    double)160.331128216631,(double)164.320112263195,(double)
	    168.327445448428,(double)172.352797139163,(double)176.395848406997,(
	    double)180.456291417544,(double)184.533828861449,(double)
	    188.628173423672,(double)192.739047287845,(double)196.86618167289,(
	    double)201.009316399282,(double)205.168199482641,(double)
	    209.342586752537,(double)213.532241494563,(double)217.736934113954,(
	    double)221.95644181913,(double)226.190548323728,(double)
	    230.439043565777,(double)234.701723442818,(double)238.978389561834,(
	    double)243.268849002983,(double)247.572914096187,(double)
	    251.890402209723,(double)256.22113555001,(double)260.564940971863,(
	    double)264.921649798553,(double)269.29109765102,(double)
	    273.673124285694,(double)278.067573440366,(double)282.47429268763,(
	    double)286.893133295427,(double)291.32395009427,(double)
	    295.766601350761,(double)300.220948647014,(double)304.686856765669,(
	    double)309.164193580147,(double)313.652829949879,(double)
	    318.152639620209,(double)322.663499126726,(double)327.185287703775,(
	    double)331.717887196928,(double)336.261181979198,(double)
	    340.815058870799,(double)345.379407062267,(double)349.95411804077,(
	    double)354.539085519441,(double)359.134205369575 };

    /* System generated locals */
    long int i__1;
    double ret_val=0.0;

    /* Local variables */
    static double dgam;
    static long int i__;
    static double t, dx, px, qx, rx, xx;
    static long int ndx, nxm;
    static double sum, rxx;

    if ( x <= (double)0.) {
	goto L90;
    } else {
	goto L5;
    }
L5:
    ndx = (long int)x;
    t = x - (double) ndx;
    if (t == (double)0.) {
	goto L51;
    }
    dx = xlim1 - x;
    if (dx < (double)0.) {
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

    dgam /= t + (double)1.;
    if (nxm == -1) {
	goto L25;
    }
    dgam /= t;
    ret_val = log(dgam);
    return ret_val;

/*     FORWARD RECURSION FOR 3.LT.X.LT.8 */

L22:
    xx = t + (double)2.;
    i__1 = nxm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dgam *= xx;
/* L24: */
	xx += (double)1.;
    }
L25:
    ret_val = log(dgam);
    return ret_val;

/*     X.GT.XLIM1 */

L40:
    rx = (double)1. / x;
    rxx = rx * rx;
    if (x - xlim2 < (double)0.) {
	goto L41;
    }
    px = q[0] * rxx + q[1];
    ret_val = px * rx + (x - (double).5) * log(x) - x + rtwpil;
    return ret_val;

/*     X.LT.XLIM2 */

L41:
    px = p[0];
    sum = (x - (double).5) * log(x) - x;
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
double rng::loggamma( double x )
{
        double f = 0.0;
        double z;
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
        double v1 = f + (x - 0.5)*log(x) - x + 0.918938533204673;

        double v2 = -0.000595238095238*z;
        double v3 = ( v2 - 0.000793650793651 )*z;
        double v4 = ( v3 - 0.0002777777777778 )*z;
        double v5 = ( v4 + 0.083333333333333 );
        double v6 = v5 / x;

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
double rng::gamain( double x, double p, double g )
{
	double pn[6];

	// define accuracy and initialize
	double acu = 1.0e-8;
	const double oflo = 1.0e30;
	double gin = 0.0;
	double term;
	double rn;
	ifault = 0;

        // test for admissibility of arguments
        if( p <= 0.0 ) ifault = 1;
        if( x < 0.0 ) ifault = 2;
        if( ifault > 0 || x == 0.0 )
                return gin;
        double factor = exp( p * log(x) - x - g );
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
        double a, b, an, dif;
        int i;
        label30:        // <<<<<<<<<<<<<<< label 30
        a = 1.0 - p;
        b = a + x + 1.0;
        term = 0.0;
        pn[0] = 1.0;
        pn[1] = x;
        pn[2] = x + 1.0;
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
                gin = 1.0 - factor*gin;
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
double rng::gauinv( double p )
{
	const double zero = 0.0;
	const double one = 1.0;
	const double half = 0.5;
        const double alimit = 1.0e-20;
        const double p0 = -0.322232431088;
        const double p1 = -1.0;
        const double p2 = -0.342242088547;
        const double p3 = -0.0204231210245;
        const double p4 = -0.0000453642210148;
        const double q0 = 0.099348462606;
        const double q1 = 0.58858157495;
        const double q2 = 0.531103462366;
        const double q3 = 0.10353775285;
        const double q4 = 0.0038560700634;

        ifault = 1;
        double ps = p;
        if( ps > half )
                ps = one - ps;
        if( ps < alimit )
                return zero;

        ifault = 0;
        if( ps == half )
                return zero;
        double yi = sqrt( log( one / (ps*ps) ) );
        double retval = yi
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
	
double rng::uniform()
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
	return ix * 4.6566128575e-10;
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
double rng::ppchi2( double p, double v )
{
	const double e = 0.5e-6;
	const double aa = 0.6931471805;
	double ch, a, b, q, p1, p2, t, x;
        double s1, s2, s3, s4, s5, s6;

#if 0
	double g = gammln( v / 2.0 );
#else
	double g = gamln( v / 2.0 );
#endif

	// after defining accuracy and ln(2), test arguments and initialize
	ifault = 1;
	if( p < 0.000002 || p > 0.999998 ) return -1.0;

	ifault = 2;
	if( v <= 0.0 ) return -1.0;

	ifault = 0;
	double xx = 0.5 * v;
	double c = xx - 1.0;

        // starting approximation for small chi-squared
        if( v >= -1.24 * log(p) )
                goto label1;
        ch = pow( (p * xx * exp( g + xx * aa)), (1.0 / xx) );
        if( ch - e < 0.0 )
                return ch;
        else
                goto label4;

        // starting approximation for v less than or equal to 0.32
	label1:
        if( v > 0.32 ) goto label3;
        ch = 0.4;
        a = log( 1.0 - p );

        label2:
        q = ch;
        p1 = 1.0 + ch * ( 4.67 + ch );
        p2 = ch * ( 6.73 + ch * ( 6.66 + ch ));
        t = -0.5 + ( 4.67 + 2.0*ch ) / p1 -
                ( 6.73 + ch*( 13.32 + 3.0*ch )) / p2;
	ch -= ( 1.0 - exp( a + g + 0.5*ch + c*aa ) * p2 / p1 ) / t;
        if( fabs( q/ch - 1.0 ) - 0.01 <= 0.0 )
                goto label4;
        else
		goto label2;

        // call to algorithm AS 70 - note that p has been tested above
        label3:
	x = gauinv( p );

	// starting approximation using Wilson and Hilferty estimate
        p1 = 0.222222 / v;
        ch = v * pow( x*sqrt(p1) + 1.0 - p1, 3.0 );

        // starting approximation for p tending to 1
        if( ch > 2.2 * v + 6.0 )
                ch = -2.0 * ( log( 1.0 - p ) - c * log( 0.5 * ch ) + g );

        // call to algorithm AS 32 and calculation of seven term Taylor series
	label4:
        q = ch;
        p1 = 0.5 * ch;
	p2 = p - gamain( p1, xx, g );

	if( ifault != 0 ) {
		ifault = 3;
                return -1.0;
	}

        t = p2 * exp( xx*aa + g + p1 - c * log(ch) );
        b = t / ch;
        a = 0.5*t - b*c;
        s1 = ( 210.0 + a*( 140.0 + a*( 105.0 + a*(84.0 + a*( 70.0 + a*60.0 ))))) / 420.0;
        s2 = ( 420.0 + a*( 735.0 + a*( 966.0 + a*( 1141.0 + a*1278.0)))) / 2520.0;
	s3 = ( 210.0 + a*( 462.0 + a*( 707.0 + a*932.0))) / 2520.0;
        s4 = ( 252.0 + a*( 672.0 + a*1182.0) + c*( 204.0 + a*( 889.0 + a*1740.0))) / 5040.0;
        s5 = (  84.0 + a*264.0 + c*(175.0 + a*606.0)) / 2520.0;
        s6 = ( 120.0 + c*( 346.0 + c*127.0)) / 5040.0;
        ch += t*( 1.0 + 0.5*t*s1 - b*c*( s1 - b*( s2 - b*( s3 - b*( s4 - b*( s5 - b*s6))))));
        if( fabs( q/ch - 1.0 ) > e )
                goto label4;
        return ch;
}

//this is from MB
void rng::DirichletRandomVariable (double *alp, double *z, int n){
	int		i;
	double	sum;

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


