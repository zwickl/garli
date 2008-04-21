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


#ifndef _RECONNODE_
#define _RECONNODE_


#include <list>
#include <algorithm>
#include "rng.h"
#ifdef UNIX
#include "unistd.h"
#endif

extern rng rnd;

using namespace std;

class ReconNode;

typedef list<ReconNode>::iterator listIt;

class ReconNode{
	public:
	unsigned short nodeNum;
	unsigned short reconDist;
	FLOAT_TYPE pathlength;
	FLOAT_TYPE weight;
	FLOAT_TYPE chooseProb;
	bool withinCutSubtree;
	
	ReconNode(unsigned short nn, unsigned short rd, float pl, bool wcs=false) : nodeNum(nn), reconDist(rd), pathlength(pl), withinCutSubtree(wcs) {}
	void Report(ofstream &deb){
		deb << nodeNum << "\t" << reconDist << "\t" << pathlength << "\t" << weight << "\t" << chooseProb << "\t" << withinCutSubtree << "\n";
		}

	bool operator<(const ReconNode &rhs){
		return reconDist < rhs.reconDist;
		}
	};

class DistEquals:public binary_function<ReconNode, int, bool>{
	public:
	result_type operator()(first_argument_type i, second_argument_type j) const{
		return (result_type) (i.reconDist==j);
		}
	};

class DistEqualsWithinCutSubtree:public binary_function<ReconNode, int, bool>{
	public:
	result_type operator()(first_argument_type i, second_argument_type j) const{
		return (result_type) (i.reconDist==j && i.withinCutSubtree==true);
		}
	};

class NodeEquals:public binary_function<ReconNode, int, bool>{
	public:
	result_type operator()(first_argument_type i, second_argument_type j) const{
		return (result_type) (i.nodeNum==j);
		}
	};

class ReconList{
	unsigned num;
	list<ReconNode> l;

	public:

	ReconList(){
		num = 0;
		}

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
		FLOAT_TYPE weightSum = 0.0, running = 0.0;
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

	void Reverse(){
		reverse(l.begin(), l.end());
		}

	ReconNode *RandomReconNode(){
		return &(*(NthElement(rnd.random_int((int)l.size()))));
		}
		
	ReconNode *ChooseNodeByWeight(){
		FLOAT_TYPE prob = rnd.uniform();
		 listIt it=l.begin();
		 for(;it!=l.end();it++){
			if(prob < (*it).chooseProb) return &(*it);
			}
		//we should only get here due to a little rounding error
		it--;
		return &(*it);
		}
	void AddNode(ReconNode &nd){
		//this just duplicates the added ReconNode via the default copy constructor, assumably added from another list
		l.push_back(nd);
		num++;
		}

	void AddNode(int nn, int rd, float pl, bool withinCutSubtree=false){	
		//first verify that we don't already have this node in the list
		if(find_if(l.begin(),l.end(),bind2nd(NodeEquals(), nn)) != l.end()) return;	
		l.push_back(ReconNode(nn, rd, pl, withinCutSubtree));
		num++;
		}	


	void SortByDist(){
		l.sort();
		}

	void DebugReport(){
		ofstream deb("recons.log");
		for(listIt it=l.begin();it!=l.end();it++){
			(*it).Report(deb);
			}
		deb.close();
		}

	};

class Swap;
bool SwapLessThan(const Swap &lhs, const Swap &rhs);
bool SwapLessThanDist(const Swap &lhs, const Swap &rhs);

class Swap{
	Bipartition b;
	unsigned short count;
	unsigned short cutnum;
	unsigned short brokenum;
	unsigned short reconDist;

public:
	Swap(Bipartition &swap, int cut, int broke, int dist){
		b=&swap;
		count=1;
		cutnum=cut;
		brokenum=broke;
		reconDist=dist;
		}
	Swap(FILE* &in){
		b.BinaryInput(in);
		intptr_t scalarSize = (intptr_t) &(reconDist) - (intptr_t) &(count) + sizeof(reconDist);
		fread(&count, scalarSize, 1, in);
		}

	void Increment(){
		count++;	
		}

	int Count()const {
		return count;	
		}

	int ReconDist() const{
		return reconDist;
		}

	void SetCount(int c){
		count = c;
		}

	void Output(ofstream &out){
		out << b.Output() << "\t" << count << "\t" << cutnum << "\t" << brokenum << "\t" << reconDist << endl;
		}
/*
	void BinaryOutput(ofstream &out){
		b.BinaryOutput(out);
		intptr_t scalarSize = (intptr_t) &reconDist - (intptr_t) &count + sizeof(reconDist);
		out.write((char*)&count, (streamsize) scalarSize);
		}
*/

	void BinaryOutput(OUTPUT_CLASS &out){
		b.BinaryOutput(out);
		intptr_t scalarSize = (intptr_t) &reconDist - (intptr_t) &count + sizeof(reconDist);
		out.WRITE_TO_FILE(&count, (streamsize) scalarSize, 1);
		}

	unsigned BipartitionBlock(int block) const{
		return b.rep[block];	
		}

	bool operator<(const Swap &rhs){
		//note that this is "less than" for sorting purposes, not in a subset sense
		//it is a strict weak ordering, so it returns false in the case of possible equality
		//ordering is based first on bip, then on reconDist
		int i;
		for(i=0;i<Bipartition::nBlocks-1;i++){
			if(BipartitionBlock(i) > rhs.BipartitionBlock(i)) return false;
			else if(BipartitionBlock(i) < rhs.BipartitionBlock(i)) return true;
			}

		if(((BipartitionBlock(i)) & Bipartition::partialBlockMask) < ((rhs.BipartitionBlock(i)) & Bipartition::partialBlockMask)) return true;
		
		else if(((BipartitionBlock(i)) & Bipartition::partialBlockMask) == ((rhs.BipartitionBlock(i)) & Bipartition::partialBlockMask)){
			//bipartitions are equal
			if(reconDist < rhs.reconDist) return true; //dists are not
			else if(reconDist == rhs.reconDist)
				if(cutnum < rhs.cutnum) return true;//cutnum is not
			}
		return false;
		}

	bool operator==(const Swap &rhs){
		assert(rhs.b.ContainsTaxon(1));
		bool bipEqual = b.EqualsEquals(rhs.b);
		if(bipEqual == false) return false;
		//if the bips are equal but the distances are different, the pre-swap topos must be different
		//so we want to consider this a different swap
		if(reconDist != rhs.reconDist) return false;
		
		if(reconDist == 1){//NNI's with different cuts and brokens can give the same topo
			return true;
			}
		else if((cutnum == rhs.cutnum) && (brokenum == rhs.brokenum)){
			return true;
			}
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
/*
	void WriteSwapCheckpoint(ofstream &out){
		intptr_t scalarSize = (intptr_t) &total - (intptr_t) &unique + sizeof(total);
		out.write((char*) &unique, (streamsize) scalarSize);
		for(list<Swap>::iterator it=swaps.begin();it != swaps.end(); it++){
			(*it).BinaryOutput(out);
			}
		}
*/

	void WriteSwapCheckpoint(OUTPUT_CLASS &out){
		intptr_t scalarSize = (intptr_t) &total - (intptr_t) &unique + sizeof(total);
		out.WRITE_TO_FILE(&unique, scalarSize, 1);
		for(list<Swap>::iterator it=swaps.begin();it != swaps.end(); it++){
			(*it).BinaryOutput(out);
			}
		}

	void ReadBinarySwapCheckpoint(FILE* &in){
		assert(ferror(in) == false);
		intptr_t scalarSize = (intptr_t) &total - (intptr_t) &unique + sizeof(total);
		fread(&unique, scalarSize, 1, in);

		for(unsigned i=0;i<unique;i++){
			Swap s(in);
			swaps.push_back(s);
			}
		IndexSwaps();

		assert(swaps.size() == unique);
		int tot=0;
		for(list<Swap>::iterator it=swaps.begin();it != swaps.end(); it++) tot += (*it).Count();
		if(tot != total) throw ErrorException("problem reading swap checkpoint!");
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
		int increment=(int) sqrt((FLOAT_TYPE)unique);
		int count=0;
		for(list<Swap>::iterator it=swaps.begin();it != swaps.end(); it++){
			if(count % increment == 0) indeces.push_back(it);
			count++;
			}
		}

	bool AddSwap(Bipartition &bip, int cut, int broke, int dist){
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
		return (found == false);//return value is true if the swap is _unique_
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
	
	void SwapReport(ofstream &swapLog){
		unsigned int distTotCounts[200];
		unsigned int distUniqueCounts[200];
		for(int i=0;i<200;i++){
			distTotCounts[i]=distUniqueCounts[i]=0;
			}
		for(list<Swap>::iterator it = swaps.begin();it != swaps.end();it++){
			distUniqueCounts[(*it).ReconDist() - 1]++;
			distTotCounts[(*it).ReconDist() - 1] += (*it).Count();
			}

		swapLog << "\t" << GetUnique() << "\t" << GetTotal() << "\t" ;
		
		for(int i=0;i<200;i++){
			if(i > 5 && distUniqueCounts[i] == 0) break;
			swapLog << distUniqueCounts[i] << "\t" << distTotCounts[i] << "\t";
			}
		swapLog << endl;
		}

	void AttemptedSwapDump(ofstream &deb){
		deb << "\t" << GetUnique() << "\t" << GetTotal() << "\n" ;
		for(list<Swap>::iterator it = swaps.begin();it != swaps.end();it++){
			(*it).Output(deb);
			}
		}

};

#endif
