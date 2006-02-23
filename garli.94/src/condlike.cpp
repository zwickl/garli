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

#include "memchk.h"
#include "condlike.h"

CondLikeArray::~CondLikeArray()
{
	//don't want to delete shared CL from nodes.  Should only
	//be called from Population level if CONDLIKE SHARED is defined
	if( arr ){
		 delete []arr;
		 arr=NULL;
		 //MEM_DELETE_ARRAY(arr); //arr is of length max
		 }
	if(underflow_mult!=NULL) delete []underflow_mult;
}

void CondLikeArray::Allocate( int nk, int ns, int nr /* = 1 */ )
{
	if( arr ){
		 delete []arr;
		 arr=NULL;
		 //MEM_DELETE_ARRAY(arr); // arr is of length max
		 }
	nrates = nr;
	nsites = ns;
	nstates = nk;
	arr=new double[nk*nr*ns];
	if(arr==NULL){
		throw(1);
		}
	//DJZ 1-5-05 - the underflow mult used to be ns*nr in length,
	//but only needs to be ns
	underflow_mult=new int[ns];
	if(underflow_mult==NULL){
		throw(1);
		}
	}


