// GARLI version 0.952b2 source code
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



#include "topologylist.h"
#include "individual.h"


long TopologyList::ntoposexamined(0);
Individual *TopologyList::listOfInd(0L);
int TopologyList::allocIncr(100);

TopologyList::TopologyList() : sizeOfAllocation(0) , indNums(0L) , currentpos(0L) , 
	gensAlive(0) , nInds(0), exNNItried(0)
{	
}

TopologyList::~TopologyList(){
	if(indNums!=NULL)
		delete []indNums;
}

void TopologyList::SetIndL(Individual *indL)
{	listOfInd=indL;
}

void TopologyList::NewGeneration()
{	currentpos=indNums;
	gensAlive++;
}

int TopologyList::GetNumberOfUnselectedInd()
{	while(*currentpos>-1)
		{if(!listOfInd[*currentpos].willreproduce && !listOfInd[*currentpos].reproduced && !listOfInd[*currentpos].willrecombine)
			return *currentpos++;
		currentpos++;
		}
	return -1;
}

int TopologyList::Allocate(int s/*=0*/)
{	if(!s)
		s=allocIncr;
	sizeOfAllocation+=s;
	assert(sizeOfAllocation<30000);
	int i=0;
	int * temp;
	if(indNums)
		{int displace=(int)(currentpos-indNums);
		int * temptw;
		temptw=temp=new int[sizeOfAllocation];
		currentpos=indNums;
		for(;i<nInds;i++)
			*temp++=*currentpos++;
		delete [] indNums;
		indNums=temptw;
		currentpos=indNums+displace;
		}
	else
		temp=currentpos=indNums=new int[sizeOfAllocation];
	for(;i<sizeOfAllocation;i++)
		*temp++=30000;
		
	return 0;
}
void TopologyList::RemoveInd(int i)
{	int *temp=indNums;
	nInds--;
	while(*temp++!=i) i ;
	do	{*(temp-1)=*temp;
		}
	while(*(temp++-1)!=-1);
}

void TopologyList::AddInd(int i)
{	if(nInds==sizeOfAllocation-1)
		Allocate();
	*(indNums+nInds++)=i;
	*(indNums+nInds)=-1;
}

void TopologyList::DecrementTopoFieldOfInds()
{	for(int i=0;i<nInds;i++)
		{assert(indNums[i]>=0);//TEMPORARY
		listOfInd[indNums[i]].topo--;
		}
}


