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
	int nodeNum;
	int reconDist;
	float pathlength;
	
	ReconNode(int nn, int rd, float pl) : nodeNum(nn), reconDist(rd), pathlength(pl) {}
	};

class DistEquals:public binary_function<ReconNode, int, bool>{
	public:
	const result_type operator()(first_argument_type i, second_argument_type j) const{
		return (result_type) (i.reconDist==j);
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
		
	void clear() {l.clear();}
	int size() {return (int) l.size();}
		
	void AddNode(int nn, int rd, float pl);
	
	void print(char *fn){
		ofstream out(fn);
		for(listIt it=l.begin();it!=l.end();it++){
			out << it->nodeNum << "\t" << it->reconDist << "\t" << it->pathlength << endl;			
			}
		out.close();
		}
		
	void RemoveNodesOfDist(int dist){
		//remove_if(l.begin(), l.end(), bind2nd(DistEquals(), dist));
		for(listIt it=l.begin();it!=l.end();){
			if((*it).reconDist==dist) it=l.erase(it);
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
		return l.erase(del);
		}
	
	listIt RemoveElement(listIt del){
		return l.erase(del);
		}
	
	int RandomNodeNum(){
		return (NthElement(rnd.random_int((int)l.size())))->nodeNum;
		}
	
	};

inline void ReconList::AddNode(int nn, int rd, float pl){
	//first verify that we don't already have this node in the list
	if(find_if(l.begin(),l.end(),bind2nd(NodeEquals(), nn)) != l.end()) return;
	l.push_back(ReconNode(nn, rd, pl));
	}

#endif