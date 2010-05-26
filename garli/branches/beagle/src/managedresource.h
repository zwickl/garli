#ifndef MANAGED_RESOURCE
#define MANAGED_RESOURCE

extern int memLevel;
extern ModelSpecification modSpec;

/*
class object{
	int index;
public:
	object(){index =0;}
	};

class derivedObject : public object{
	
	};

class base{
protected:
	object *basePointer;
	object **baseDoublePointer;
public:
	base(){
		basePointer = NULL;
		baseDoublePointer = NULL;
		}
	};

class derived : public base{
public:
	derived(){
		basePointer = new derivedObject[1];
		baseDoublePointer = (object **) new derivedObject*[1];
		}
	};
*/

class Resource{
protected:
	int index;
public:
	Resource(){
		index = -1;
		}
	int Index(){return index;}
	};

class CondLikeArray2 : public Resource{
protected:
	unsigned nsites, nrates, nstates;

	FLOAT_TYPE *arr;
	int* underflow_mult;
	unsigned rescaleRank;

public:
	CondLikeArray2() : Resource(){
		nsites = nrates = nrates = 0;
		arr = NULL;
		underflow_mult = NULL;
		rescaleRank = 1;
		}
	~CondLikeArray2();
	int NStates() {return nstates;}
	int NSites() {return nsites;}
	int Index() {return index;}

	void Allocate(int ind, int nk, int ns, int nr);
	};	

class ResourceHolder{
	friend class ResourceManager;
protected:
	int numAssigned;
	bool reserved;
	bool tempReserved;
	int depLevel;
	int reclaimLevel;
	Resource *theResource;
public:
	ResourceHolder(){
		numAssigned = 0;
		reserved = tempReserved = false;
		depLevel = -1;
		reclaimLevel = -1;
		theResource = NULL;
		}
	virtual void Reset() = 0;
	int ResourceIndex(){return theResource->Index();}
	bool IsDirty() const {return theResource == NULL;}
	int GetReclaimLevel() {return reclaimLevel;}
	};

class CondLikeArrayHolder2 : public ResourceHolder{
	friend class ClaManager2;

	int holderDep1;
	int holderDep2;
	int transMatDep1;	
	int transMatDep2;

	void SetHolderDependencies(int depIndex1, int pDepIndex1, int depIndex2, int pDepIndex2);
public:
	CondLikeArrayHolder2() : ResourceHolder(){
		holderDep1 = holderDep2 = transMatDep1 = transMatDep2 = -1;
		theResource = NULL;
		}
	void Reset(){
		numAssigned = 0;
		reserved = tempReserved = false;
		depLevel = reclaimLevel = -1;

		holderDep1 = holderDep2 = transMatDep1 = transMatDep2 = -1;
		theResource = NULL;
		}
	}; 

class ResourceManager{
protected:
	int numHolders;
	int numResources;

	int curReclaimLevel;
	int memLevel;

	vector<int> resourceIndexStack;
	vector<int> holderIndexStack;

	//probably shouldn't use this, but carry over from old clamanager
	vector<Resource*> resourceStack;

	Resource **resources;
	ResourceHolder *holders;

	bool debug_resources;
	int maxResourcesUsed;
public:
	ResourceManager(){
		numHolders = -1;
		numResources = -1;
		curReclaimLevel = -1;
	
		resources = NULL;
		holders = NULL;

		debug_resources = false;
		maxResourcesUsed = 0;
		}
	~ResourceManager(){
		if(resources != NULL){
			for(int i = 0;i < numResources;i++){
				delete resources[i];
				}
			delete []resources;
			}
		delete []holders;
		resourceIndexStack.clear();
		holderIndexStack.clear();
		resourceStack.clear();
		}
protected:
	int _GetFreeHIndex(){
		assert(holderIndexStack.size() > 0);
		int newIndex = holderIndexStack[holderIndexStack.size()-1];		
		holderIndexStack.pop_back();
		IncrementHolder(newIndex);
		return newIndex;
		}

	Resource *_GetFreeR(){
		if(resourceStack.empty() == true)	
			RecycleResources();

		Resource *arr=resourceStack[resourceStack.size()-1];
		resourceStack.pop_back();

		if(numResources - (int)resourceStack.size() > maxResourcesUsed) 
			maxResourcesUsed = numResources - (int)resourceStack.size();
		if(debug_resources){
			ReportHolderTotals("assignedC\0", arr->Index());
			}
		return arr;		
		}

public:
	const ResourceHolder *GetHolder(int index){
		assert(index > -1 && index < numHolders);
		return &holders[index];
		}

	ResourceHolder *GetMutableHolder(int index){
		assert(index > -1 && index < numHolders);
		return &holders[index];
		}

	int GetReclaimLevel(int index){
		assert(index > -1 && index < numHolders);
		//TODO
//		if(holders[index].theArray == NULL) 
//			return -1;
		return holders[index].GetReclaimLevel();
		}

	//basic setting and getting
	void SetDepLevel(int index, int lvl) {holders[index].depLevel = lvl;}
	void RemoveNormalReservation(int index);
	void RemoveTempReservation(int index);
	void ReserveResource(int index, bool temp=true);
	void SetCurReclaimLevel(int lvl){curReclaimLevel = lvl;}

	bool IsHolderDirty(int index) {return holders[index].IsDirty();}
	bool IsResourceReserved(int index) {return holders[index].reserved;}
	bool IsResourceTempReserved(int index) {return holders[index].tempReserved;};

	//basic operations
	void IncrementHolder(int index);
	void DecrementHolder(int index);
	void RecycleResources();
	int AssignFreeHolder();
	int TradeInHolder(int oldIndex);
	void FillHolder(int index, int dir);
	void EmptyHolder(int index);
	Resource *GetHolderFillIfNecessary(int index);
	void ClaimHolderFillIfNecessary(int index, int depLevel);
	int SetHolderDirty(int index);

	//for debugging
	void ReportHolderTotals(const char *mess, int index) const;
	void CheckHolders() const;
	void OutputResourceReport() const;
	};

class ClaManager2 : public ResourceManager{
protected:	
	int numNodes;//the number of nodes in each tree
	int numRates;
	int numStates;
	int numChar;

public:
	ClaManager2() : ResourceManager() {
		numNodes = numRates = numStates = -1;		
		}

	ClaManager2(int nnod, int nClas, int nHolders, int nchar, int nrates, int nstates){
		curReclaimLevel = -1;
		numNodes = nnod;
		numResources = nClas;
		numHolders = nHolders;
		numRates = nrates;
		numStates = nstates;
		numChar = nchar;

		holders = new CondLikeArrayHolder2[numHolders];
		resources = (Resource **) new CondLikeArray2*[numResources];
		
		for(int i=numResources-1;i>=0;i--){
			resources[i] = new CondLikeArray2;
			((CondLikeArray2*) resources[i])->Allocate(i, numStates, numChar, numRates);
			resourceIndexStack.push_back(i);
			resourceStack.push_back(resources[i]);
			}
		
		for(int i=numHolders-1;i>=0;i--)
			holderIndexStack.push_back(i);
#ifdef CLA_DEBUG
		debug_resources = true;
#else
		debug_resources = false;
#endif
		}
	
	void SetHolderDependencies(int index, int depIndex1, int pDepIndex1, int depIndex2, int pDepIndex2);

	int GetClaIndexForBeagle(int index) const{
		assert(index > -1 && index < numHolders);
		assert(holders[index].IsDirty());
		return holders[index].ResourceIndex();
		}

	};

#endif
