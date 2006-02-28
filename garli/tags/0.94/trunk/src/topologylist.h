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


#ifndef TOPOLOGY_LIST
#define TOPOLOGY_LIST

class Individual;

class TopologyList
{	static Individual *listOfInd;
	static int allocIncr;
	int sizeOfAllocation;
	int *indNums;
	int *currentpos;
	public :
	bool exNNItried;
	int nInds;
	int gensAlive;
	static long ntoposexamined;
	
	TopologyList();
	~TopologyList();
	static void SetIndL(Individual *indL);
	int Allocate(int s=0);
	void NewGeneration();
	void Clear()	{currentpos=indNums; *currentpos=-1; nInds=0;}
	int GetNumberOfUnselectedInd();
	void AddInd(int i);
	void RemoveInd(int i);
	void DecrementTopoFieldOfInds();
	int GetIndNums(int i) {return indNums[i];}
	Individual *GetListOfInd(){return listOfInd;}
};
#endif

