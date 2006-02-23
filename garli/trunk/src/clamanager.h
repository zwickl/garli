#ifndef CLA_MANAGER
#define CLA_MANAGER

#include <vector>
#include <algorithm>
#include <cassert>
#include "memchk.h"
#include "condlike.h"
using namespace std;

extern int memLevel;

#undef CLA_DEBUG

class ClaManager{
	int numNodes;//the number of nodes in each tree
	int numRates;
	int numClas;
	int numHolders;
	int maxUsed;
	CondLikeArray **allClas; //these are the actual arrays to be used in calculations, but will assigned to 
							 //nodes via a CondLikeArrayHolder.  There may be a limited number
							 
	CondLikeArrayHolder *holders; //there will be enough of these such that every node and direction could
								  //have a unique one, although many will generally be shared
								  
//	CondLikeArray **tempClas;
	vector<CondLikeArray *> claStack;
	vector<int> holderStack;
	
	public:
	ClaManager(int nnod, int nClas, int nHolders, int nchar, int nrates) : numNodes(nnod), numClas(nClas), numHolders(nHolders), numRates(nrates){
		maxUsed=0;
		allClas=new CondLikeArray*[numClas];
		for(int i=numClas-1;i>=0;i--){
			allClas[i]=new CondLikeArray;
			allClas[i]->Allocate(4, nchar, numRates);
			claStack.push_back(allClas[i]);
			}
		holders = new CondLikeArrayHolder[numHolders];
		for(int i=numHolders-1;i>=0;i--)
			holderStack.push_back(i);
		}

	~ClaManager(){
		if(allClas!=NULL){
			for(int i=0;i<numClas;i++){
				delete allClas[i];
				}
			delete []allClas;
			}
		delete []holders;
		}
	
	int NumClas(){
		return numClas;
		}	
	
	int NumFreeClas(){
		return claStack.size();
		}
	
	int NumFreeHolders(){
		return holderStack.size();
		}
	
	
	int AssignClaHolder(){
		int index=holderStack[holderStack.size()-1];
		IncrementCla(index);
		holderStack.pop_back();
		return index;
		}
	
	void FillHolder(int index, int dir){
		holders[index].theArray = AssignFreeCla();
		holders[index].reclaimLevel=dir;
		}
	
	CondLikeArray* AssignFreeCla(){
		#ifdef CLA_DEBUG
		ofstream deb("cladebug.log", ios::app);
		#endif

		if(claStack.empty() == true) RecycleClas();
		
		CondLikeArray *arr=claStack[claStack.size()-1];
		assert(arr != NULL);
		claStack.pop_back();
		if(numClas - claStack.size() > maxUsed) maxUsed=numClas - claStack.size();
		
		return arr;
		}
/*
	void MarkReclaimable(int index, int val, bool observeCounts=true){
		assert(0);
		#ifdef CLA_DEBUG
		ofstream deb("cladebug.log", ios::app);
		deb << index << "\tmarked reclaimed\n";
		#endif
//		if(holders[index].theArray != NULL && (observeCounts==false || holders[index].nodes.size()==1)){
//		if(holders[index].theArray != NULL && (observeCounts==false || holders[index].numAssigned==1)){
//		if(holders[index].theArray != NULL && holders[index].reserved==false){
//			holders[index].SetReclaimLevel(val);

//		if(holders[index].theArray != NULL){
//			if(holders[index].reserved==false) holders[index].SetReclaimLevel(val);
			holders[index].tempReserved = false;
//			}
		}
*/	
	int GetReclaimLevel(int index){
		if(holders[index].theArray == NULL) return -1;
		return holders[index].GetReclaimLevel();
		}

	int GetNumAssigned(int index){
		return holders[index].numAssigned;
		}

	void ReserveCla(int index, bool temp=true){
		if(temp==true) holders[index].tempReserved=true;
		else holders[index].reserved=true;
		}

	void ClearTempReservation(int index){
		holders[index].tempReserved=false;
		}

	void UnreserveCla(int index){
//		holders[index].tempReserved=false;
		holders[index].reserved=false;
		if(memLevel>1)
			holders[index].SetReclaimLevel(1);
		}

	bool IsClaReserved(int index){
		return holders[index].reserved;
		}

	bool IsClaTempReserved(int index){
		return holders[index].tempReserved;
		}

	void ReclaimSingleCla(int index){
		//this simply removes the cla from a holder.  It is equivalent to just
		//dirtying it if only a single tree shares the holder
		if(holders[index].theArray==NULL) return;
		claStack.push_back(holders[index].theArray);
		holders[index].SetReclaimLevel(0);
		holders[index].theArray=NULL;				
		}

	void CountClaTotals(int &clean, int &tempres, int &res){
		for(int i=0;i<numHolders;i++){
			if(holders[i].theArray != NULL) clean++;
			if(holders[i].tempReserved ==true) tempres++;
			if(holders[i].reserved==true) res++;
			}		
		}


	void RecycleClas(){
		int numReclaimed=0;
		for(int i=0;i<numHolders;i++){
			if(holders[i].theArray != NULL){
				if(holders[i].GetReclaimLevel() == 2 && holders[i].tempReserved == false && holders[i].reserved == false){
					assert(holders[i].theArray->NStates()==4);
					claStack.push_back(holders[i].theArray);
					holders[i].SetReclaimLevel(0);
					holders[i].theArray=NULL;
					numReclaimed++;
					}
				}
			if(memLevel<2) if(numReclaimed>50) return;
			}
		if(numReclaimed>10) return;
		for(int i=0;i<numHolders;i++){
			if(holders[i].theArray != NULL){
				if((holders[i].GetReclaimLevel() == 1 && holders[i].tempReserved == false && holders[i].reserved == false)){
					assert(holders[i].theArray->NStates()==4);
					claStack.push_back(holders[i].theArray);
					holders[i].SetReclaimLevel(0);
					holders[i].theArray=NULL;
					numReclaimed++;
					}
				}
			if(numReclaimed==20) return;
			}
		if(numReclaimed==0){
			throw(2);
			}
		assert(numReclaimed > 0);
		}
	
	int GetClaNumber(int index){
		//this is ugly, but should only be called for debugging
		if(holders[index].theArray == NULL) return -1;
		for(int i=0;i<numClas;i++)
			if(holders[index].theArray == allClas[i]) return i;
		assert(0);
		return -1;
		}

	int CountClasInUse(int recLevel){
		int num=0;
		for(int i=0;i<numHolders;i++){
			if(holders[i].theArray != NULL)
				if(holders[i].GetReclaimLevel() == recLevel) num++;
			}
		return num;
		}

	CondLikeArray *GetCla(int index){
		assert(holders[index].theArray != NULL);
		return holders[index].theArray;
		}
	
	const CondLikeArrayHolder *GetHolder(int index){
		return &holders[index];
		}
	
	bool IsDirty(int index){
		//dirtyness is now synonymous with a null cla pointer in the holder
		assert(index > -1);
		return (holders[index].theArray == NULL);	
		}

	int SetDirty(int index){
		//there are only two options here:
		//1. Cla is being made dirty, and only node node points to it 
		//	->null the holder's cla pointer and return the same index
		//2. Cla is being made dirty, and multiple nodes point to it
		//	->remove this node from the holder (decrement) and assign a new one	
	
		if(holders[index].numAssigned==1){
			if(holders[index].theArray != NULL){
				holders[index].SetReclaimLevel(0);
				claStack.push_back(holders[index].theArray);
				holders[index].theArray=NULL;
				}
			}
		else{
			DecrementCla(index);
			index=holderStack[holderStack.size()-1];
			holderStack.pop_back();
			IncrementCla(index);
			}
		return index;
		}

	void IncrementCla(int index){
		holders[index].numAssigned++;
		}

	void DecrementCla(int index){
		assert(index != -1);
		if(holders[index].numAssigned==1){
			holderStack.push_back(index);
			if(holders[index].theArray != NULL){
				assert(find(claStack.begin(), claStack.end(), holders[index].theArray) == claStack.end());
				assert(holders[index].theArray->NStates()==4);
				claStack.push_back(holders[index].theArray);
				}
			holders[index].Reset();
			}
		else{
			holders[index].numAssigned--; 
			//this is important!
			holders[index].tempReserved=false;
			}
		}

/*		
	void CheckClaManager(int checktot){
		int tot=0;
		//verify that the total number of clas assigned for any particular node is less than the number of trees
		for(int i=1;i<size;i++){
			tot+=assignedClaArray[i];
			}
		}
*/
/*	
	void CheckAssignedNumber(int chk, int node, int index){
		assert(chk==assignedClaArray[node*numCopies+index]);
		}
*.		
//	int NumCopies(){ return numCopies;}
	int NumNodes(){return numNodes;}
/*	
	void OutputAssignedClaArray(){
		static int count = 0;
		//ofstream out("claindeces.txt", ios::app);
		ofstream out("claindeces.txt");
		out << "\ncount " << count++ << endl;
		for(int i=0;i<numNodes;i++){
			for(int j=0;j<numCopies;j++){
				out << i << "\t" << j << "\t" << assignedClaArray[i*numCopies+j] << endl;
				}
			}
		out << endl;
		}
	void ClearAllClas(){
		for(int i=0;i<numNodes*numCopies;i++){
			allClas[i]->SetDirty(true);
			assignedClaArray[i]=0;
			}
		}
*/	
/*	void OutputClaReport(){
		ofstream cla("clareport.log");
		for(int i=1;i<numClas;i++){
			cla << i << "\t" << assignedClaArray[i] << "\t" << allClas[i]->nodeNum << "\t" << allClas[i]->IsDirty() << "\n";
			}
		}
*/
/*
	int NumAssigned(int index){
		return assignedClaArray[index];
		}	
*/
/*
	void OutputClaInfo(ostream &str, int index, int nd){
		str << nd << "\t" << index << "\t" << allClas[index]->nodeNum << "\t" << allClas[index]->IsDirty(nd) << "\t" << assignedClaArray[index] << "\t" << allClas[index]->GetReclaimLevel() << "\n";
		}
*/

	
	void CheckClaHolders(){
		int used=0;
		int reclaim2=0;
		for(int i=0;i<numHolders;i++){
			if(holders[i].theArray != NULL){
				used++;
				if(holders[i].GetReclaimLevel() == 2) reclaim2++;
				}
			}
		assert(used == numClas - claStack.size());
		}
	
	void MakeAllHoldersDirty(){
		
		for(int i=0;i<numHolders;i++){
			if(holders[i].theArray != NULL){
				claStack.push_back(holders[i].theArray);
				holders[i].theArray=NULL;
				}
			
			
			}
		}
	};		
#endif
