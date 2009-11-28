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

