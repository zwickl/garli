// GARLI version 0.95b6 source code
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


#ifndef _RECONNODE_
#define _RECONNODE_


#include <list>
#include <algorithm>
#include "rng.h"

extern rng rnd;

using namespace std;

class ReconNode;

typedef list<ReconNode>::iterator listIt;

class ReconNode{
	public:
	unsigned short nodeNum;
	unsigned short reconDist;
	float pathlength;
	double weight;
	float chooseProb;
	bool withinCutSubtree;
	
	ReconNode(unsigned short nn, unsigned short rd, float pl, bool wcs=false) : nodeNum(nn), reconDist(rd), pathlength(pl), withinCutSubtree(wcs) {}
	};

class DistEquals:public binary_function<ReconNode, int, bool>{
	public:
	const result_type operator()(first_argument_type i, second_argument_type j) const{
		return (result_type) (i.reconDist==j);
		}
	};

class DistEqualsWithinCutSubtree:public binary_function<ReconNode, int, bool>{
	public:
	const result_type operator()(first_argument_type i, second_argument_type j) const{
		return (result_type) (i.reconDist==j && i.withinCutSubtree==true);
		}
	};

class NodeEquals:public binary_function<ReconNode, int, bool>{
	public:
	const result_type operator()(first_argument_type i, second_argument_type j) const{
		return (result_type) (i.nodeNum==j);
		}
	};

class ReconList{
	list<ReconNode> l;
	unsigned num;

	public:
	
	listIt begin(){
		return l.begin();
		}

	listIt end(){
		return l.end();
		}

	listIt GetFirstNodeAtDist(int Dist){
		return find_if(l.begin(),l.end(),bind2nd(DistEquals(), Dist));
		}

	listIt GetFirstNodeAtDistWithinCutSubtree(int Dist){
		return find_if(l.begin(),l.end(),bind2nd(DistEqualsWithinCutSubtree(), Dist));
		}

	void clear() {
		l.clear();
		num=0;
		}
	unsigned size() {
		assert(num == l.size());
		return num;
		}
		
	void print(const char *fn){
		ofstream out(fn);
		for(listIt it=l.begin();it!=l.end();it++){
			out << it->nodeNum << "\t" << it->reconDist << "\t" << it->pathlength << endl;			
			}
		out.close();
		}

	void CalcProbsFromWeights(){
		//this just fills the chooseProb field by dividing the prob between the nodes in proportion to their weight
		double weightSum = 0.0, running = 0.0;
		for(listIt it=l.begin();it!=l.end();it++){
			weightSum += (*it).weight;
			}
			
		for(listIt it=l.begin();it!=l.end();it++){
			running += (*it).weight/weightSum;
			(*it).chooseProb = (float) running;
			}
		}
		
	void RemoveNodesOfDist(int dist){
		//remove_if(l.begin(), l.end(), bind2nd(DistEquals(), dist));
		for(listIt it=l.begin();it!=l.end();){
			if((*it).reconDist==dist){
				it=l.erase(it);
				num--;
				}
			else it++;
			}
		}
		
	listIt NthElement(int index){
		listIt ret=l.begin();
		int i=0;
		while(i++<index){
			ret++;
			assert(ret != l.end());
			}
		return ret;
		}
	
	listIt RemoveNthElement(int index){
		listIt del=l.begin();
		int i=0;
		while(i++<index){
			del++;
			assert(del != l.end());
			}
		num--;
		return l.erase(del);
		}
	
	listIt RemoveElement(listIt del){
		num--;
		return l.erase(del);
		}
	
	int RandomNodeNum(){
		return (NthElement(rnd.random_int((int)l.size())))->nodeNum;
		}

	ReconNode *RandomReconNode(){
		return &(*(NthElement(rnd.random_int((int)l.size()))));
		}
		
	ReconNode *ChooseNodeByWeight(){
		double prob = rnd.uniform();
		 listIt it=l.begin();
		 for(;it!=l.end();it++){
			if(prob < (*it).chooseProb) return &(*it);
			}
		//we should only get here due to a little rounding error
		it--;
		return &(*it);
		}
	void AddNode(int nn, int rd, float pl, bool withinCutSubtree=false){	
		//first verify that we don't already have this node in the list
		if(find_if(l.begin(),l.end(),bind2nd(NodeEquals(), nn)) != l.end()) return;	
		l.push_back(ReconNode(nn, rd, pl, withinCutSubtree));
		num++;
		}
	};

class Swap;
bool SwapLessThan(const Swap &lhs, const Swap &rhs);

class Swap{
	Bipartition b;
	unsigned short count;
	unsigned short cutnum;
	unsigned short brokenum;
	unsigned short reconDist;

/*	int count;
	int cutnum;
	int brokenum;
	int reconDist;
*/
public:
	Swap(Bipartition &swap, int cut, int broke, int dist){
		b=&swap;
		count=1;
		cutnum=cut;
		brokenum=broke;
		reconDist=dist;
		}
	void Increment(){
		count++;	
		}

	int Count(){
		return count;	
		}

	void SetCount(int c){
		count = c;
		}

	void Output(ofstream &out){
		out << b.Output() << "\t" << count << "\t" << cutnum << "\t" << brokenum << "\t" << reconDist << "\n";
		}

	unsigned BipartitionBlock(int block) const{
		return b.rep[block];	
		}

	bool operator<(const Swap &rhs){
		//note that this is "less than" for sorting purposes, not in a subset sense
		//returns false is in the case of equality of bipatitions AND nodenums, but
		//otherwise sorts based on the cut node
		int i;
		for(i=0;i<Bipartition::nBlocks-1;i++){
			if(BipartitionBlock(i) > rhs.BipartitionBlock(i)) return false;
			else if(BipartitionBlock(i) < rhs.BipartitionBlock(i)) return true;
			}

		if(((BipartitionBlock(i)) & Bipartition::partialBlockMask) < ((rhs.BipartitionBlock(i)) & Bipartition::partialBlockMask)) return true;
		else{
			if(((BipartitionBlock(i)) & Bipartition::partialBlockMask) == ((rhs.BipartitionBlock(i)) & Bipartition::partialBlockMask) && cutnum < rhs.cutnum) return true;
			}
		return false;
		}

	bool operator==(const Swap &rhs){
		bool bipEqual = b.EqualsEquals(&rhs.b);
		if(reconDist == 1){//NNI's with different cuts and brokens can give the same topo, so just look at the bip
			if(bipEqual == true) return true;
			}
		else{
			if((bipEqual == true) && (cutnum == rhs.cutnum) && (brokenum == rhs.brokenum)) return true;
			}
//		if((bipEqual == true) && (cutnum == rhs.cutnum) && (brokenum == rhs.brokenum)) return true;
		return false;
		}
	};

class AttemptedSwapList{
	list<Swap> swaps;
	list<list<Swap>::iterator> indeces;
	unsigned unique;
	unsigned total;
	
public:

	AttemptedSwapList(){
		unique=total=0;	
		}

	int GetUnique() {return unique;}
	int GetTotal() {return total;}

	void ClearAttemptedSwaps(){
		swaps.clear();
		indeces.clear();
		unique=total=0;
		}

	list<Swap>::iterator end(){
		return swaps.end();
		}

	void WriteSwapCheckpoint(ofstream &out){
		for(list<Swap>::iterator it=swaps.begin();it != swaps.end(); it++){
			(*it).Output(out);
			}
		}

	void ReadSwapCheckpoint(ifstream &in, int ntax){
		assert(in.good());
		Bipartition *b;
		char *str=new char[ntax+2];
		int count, cut, broke, dist;
		in >> str;
		while(in.good() && !in.eof()){
			b=new Bipartition(str);
			in >> count;
			in >> cut;
			in >> broke;
			in >> dist;
			Swap swap(*b, cut, broke, dist);
			swap.SetCount(count);
			unique++;
			total+=count;
			swaps.push_back(swap);
			in >> str;
			}
		IndexSwaps();
		delete []str;
		}

	void IndexSwaps(){
		indeces.clear();
		int increment=(int) sqrt((double)unique);
		int count=0;
		for(list<Swap>::iterator it=swaps.begin();it != swaps.end(); it++){
			if(count % increment == 0) indeces.push_back(it);
			count++;
			}
		}

	void AddSwap(Bipartition &bip, int cut, int broke, int dist){
		//see if the bipartition already exists in the list
		//if so, increment the count, otherwise add it
		assert(bip.ContainsTaxon(1));

		Swap *swap = new Swap(bip, cut, broke, dist);

		bool found;
		list<Swap>::iterator it = FindSwap(*swap, found);

		if(found == false){
			bool reindex=false;
			//if we're adding this before the first index, be sure to reindex
			if(it == swaps.begin() && indeces.empty()==false) reindex=true;
			swaps.insert(it, *swap);
			unique++;
			total++;
			if(unique==100 || (unique % 1000)==0 || reindex==true) IndexSwaps(); 
			}
		else{
			(*it).Increment();
			total++;
			}
		assert(swaps.size() == unique);
		delete swap;
		}

	list<Swap>::iterator FindSwap(Swap &swap, bool &found){
		//this function returns the matching swap if found in the list
		//or the swap that would come immediately after it if not
		
		list<Swap>::iterator start;

		if(indeces.size() == 0) start=swaps.begin();
		else{
			for(list<list<Swap>::iterator>::iterator indexit=indeces.begin();;indexit++){
				if(indexit == indeces.end()){
					start = *(--indexit);
					break;
					}
				else if(swap < (*(*indexit))){
					if(indexit != indeces.begin()) start = *(--indexit);
					else start = *(indeces.begin());
					break;
					}
				}
			}
		
		for(list<Swap>::iterator it = start;it != swaps.end();it++){
			if(swap == (*it)){
				found=true;
				return it;
				}
			if(swap < (*it)){
				found=false;
				return it;
				}
			}


/*
		//a complete search from the start
		for(list<Swap>::iterator it = swaps.begin();it != swaps.end();it++){
			if(swap == (*it)){
				found=true;
				return it;
				}
			if(swap < (*it)){
				found=false;
				return it;
				//return swaps.end();
				}
			}
*/		found=false;
		return swaps.end();
		}
	};

#endif
