// GARLI version 1.00 source code
// Copyright 2005-2010 Derrick J. Zwickl
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

#include "defs.h"
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


