// GARLI version 0.94 source code
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

#include <iostream>
#include <cassert>

using namespace std;

#include "subset.h"

// Added by Alan for subset SPR

subset::subset(){
  total = 0;
  for(int i =0;i<1024;i++){front[i]=0;element[i]=0;}
}
void subset::setseed(int seed, double dist/*=0.0*/){
  seednumber = seed;
  total = 1;
  element[0] = seed;
  front[0] = 0;
  pathlength[0]=dist;
}
void subset::addelement(int elenumber,int elemfront, double dist){
  for(int i = 0;i<total;i++)
    {
      if(element[i]==elenumber) return;
    }
  total ++;
  element[total-1] = elenumber;
  front[total-1]   = elemfront;
  pathlength[total-1] = dist;
}

void subset::setfront(int elenumber,int elemfront){
  for(int i = 0;i<total;i++)
    {
      if(element[i]==elenumber)front[i] = elemfront;
    }
}
int subset::getfront(int elenumber){
  for(int i = 0; i<total;i++){
    if(element[i]==elenumber) return(front[i]);
  }
  assert(0);
  return 0;
}

void subset::elementremove(int elenumber){
	 for(int i = 0; i<total;i++){
	 	if(element[i]==elenumber){
	 		element[i]=0;
	 		break;
	 		}
	 	}
	}

void subset::removennis(){
	//this just removes everything with a distance of 1
	for(int i = 0; i<total;i++){
 		if(front[i]==1){
 			element[i]=0;
 			}
		}
	}

void subset::compact(){
	int offset=0;
	for(int i=0;i<total;i++){
		if(element[i+offset]==0 && (i+offset) < total)
			while(element[i+offset]==0 && (i+offset) < total) offset++;
		element[i]=element[i+offset];
		//need to add these, dummy.
		front[i]=front[i+offset];
		pathlength[i]=pathlength[i+offset];
		}
	total -= offset;
	} 
void subset::clear(){
  total = 0;
  seednumber = -1;
  for(int i = 0;i<1024;i++)
    {  element[i]=0;
    front[i] = 0;}
}
void subset::print(){
  cout <<endl<<"  Print out for debug of subset SPRmutaion:"<<endl;
  cout <<"   Total number of the range : "<< total << endl;
  cout << "  Seed number is: "<< seednumber << endl;
	}
