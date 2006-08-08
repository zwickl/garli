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


#include "tree.h"
#include "model.h"
#include "funcs.h"

//a bunch of functions from the Tree class, relating to optimization

extern double globalBest;

extern int optCalcs;

#undef OPT_DEBUG

#define FOURTH_ROOT


#ifdef OPT_DEBUG
#include "optimizationinfo.h"
OptimizationInfo optInfo;
ofstream opt("optimization.log");
//ofstream der("derivs.log");
ofstream optsum("optsummary.log");
ofstream curves("curves.log");
#endif

#ifdef FOURTH_ROOT
#define effectiveMin 0.01
#define effectiveMax 1.77827941
#elif ROOT_OPT
#define effectiveMin 0.0001
#define effectiveMax 3.16227766
#else
#define effectiveMin=DEF_MIN_BRLEN
#define effectiveMax=DEF_MAX_BRLEN
#endif

typedef double GAMLfloat;

inline GAMLfloat CallBranchLike(TreeNode *thisnode, Tree *thistree, GAMLfloat blen, bool brak /*=false*/){
	brak;

#ifdef FOURTH_ROOT
	thisnode->dlen=blen*blen*blen*blen;
#elif ROOT_OPT
	thisnode->dlen=blen*blen;
#else
	thisnode->dlen=blen;
#endif
	double like=thistree->BranchLike(thisnode)*-1;
	
	optCalcs++;

#ifdef OPT_DEBUG
	if(brak) optInfo.BrakAdd(blen, like);
	else optInfo.BrentAdd(blen, like);
#endif

	return like;
	}

void Tree::OptimizeBranchesInArray(int *nodes, int numNodes, double optPrecision){
	//this takes an array of nodeNums (branches) to be optimized and does so
	for(int i=0;i<numNodes;i++){
		BrentOptimizeBranchLength(optPrecision, allNodes[nodes[i]], true);
		}	
	}

double Tree::OptimizeAllBranches(double optPrecision){
	double improve=0.0;
	improve = RecursivelyOptimizeBranches(root->left, optPrecision, 0, numNodesTotal, true, improve, true);
	improve = RecursivelyOptimizeBranches(root->left->next, optPrecision, 0, numNodesTotal, true, improve, true);
	improve = RecursivelyOptimizeBranches(root->right, optPrecision, 0, numNodesTotal, true, improve, true);

	return improve;
	}

double Tree::OptimizeTreeScale(){
	if(lnL==-1.0) Score();
	Score();
	double start=lnL;
	double prev=lnL;
	double cur;
	double scale;

/*
	ofstream deb("debug.log");
	deb.precision(20);
	for(int s=0;s<30;s++){
		double scale=.85 + s*.01;
		ScaleWholeTree(scale);
		Score();
		deb << scale << "\t" << lnL << endl;
		ScaleWholeTree(1.0/scale);	
		}
	deb.close();
*/

	while(1){
		double incr=0.0001;
		scale=1.0 + incr;
		ScaleWholeTree(scale);
		Score();
		cur=lnL;
		double d11=(cur-prev)/incr;

		//return the tree to its original scale		
		ScaleWholeTree(1.0/scale);
		scale=1.0-incr;

		ScaleWholeTree(scale);
		Score();
		cur=lnL;
		double d12=(cur-prev)/-incr;
		
		double d1=(d11+d12)*.5;
		double d2=(d11-d12)/incr;
		
		double est=-d1/d2;
		
		//return the tree to its original scale		
		ScaleWholeTree(1.0/scale);
		if(abs(est) < 0.001 && d2 < 0.0) return prev-start;
		
		double t;
		if(d2 < 0.0){
			t=1.0 + est;
			if(t<.95) t=.95;
			}
		else{
			if(d1 > 0.0) t=1.01;
			else t=.99;
			}
		
		scale=t;
		ScaleWholeTree(scale);
		Score();
		cur=lnL;
		prev=cur;
		}
	return -1;
	}

double Tree::OptimizeAlpha(){

/*
	double initVal=mod->Alpha();

	ofstream deb("debug.log");
	deb.precision(20);
	for(int s=0;s<30;s++){
		double scale=.45 + s*.01;
		mod->SetAlpha(scale);
		MakeAllNodesDirty();
		Score();
		deb << scale << "\t" << lnL << endl;
		}
	deb.close();
	mod->SetAlpha(initVal);
	MakeAllNodesDirty();
	Score();	
*/	

	if(lnL==-1) Score();
	double start=lnL;
	double prev=lnL;
	double cur;
	double prevVal=mod->Alpha();;
	
	while(1){
		double incr=0.001;
		mod->SetAlpha(prevVal+incr);
		
		MakeAllNodesDirty();
		Score();
		cur=lnL;
		double d11=(cur-prev)/incr;
		
		mod->SetAlpha(prevVal-incr);
		MakeAllNodesDirty();
		Score();
		cur=lnL;
		double d12=(cur-prev)/-incr;
		
		double d1=(d11+d12)*.5;
		double d2=(d11-d12)/incr;
		
		double est=-d1/d2;
		
		if((abs(est) < 0.001 && d2 < 0.0) || (abs(d1) < 50.0)){
			mod->SetAlpha(prevVal);
			MakeAllNodesDirty();			
			return prev-start;
			}
		
		double t;
		if(d2 < 0.0){
			//don't move too much in any one step, or it tends to way overshoot
			if(abs(est) > 0.1) est *= 0.5;
			t=prevVal+est;
			if(t < 0.05) t=prevVal*.5;
			}
		else{
			if(d1 > 0.0) t=prevVal*1.1;
			else t=prevVal*0.91;
			}
		
		mod->SetAlpha(t);
		assert((prevVal==0.05 && mod->Alpha()==0.05)==false);
		MakeAllNodesDirty();			
		Score();
		prev=lnL;
		prevVal=t;
		}
	return -1;
	}

double Tree::OptimizePinv(){


/*	ofstream deb("debug.log");
	deb.precision(20);
	for(int s=0;s<30;s++){
		double scale=.15 + s*.01;
		mod->SetPinv(scale);
		MakeAllNodesDirty();
		Score();
		deb << scale << "\t" << lnL << endl;
		}
	deb.close();
*/
	if(lnL==-1) Score();
	double start=lnL;
	double prev=lnL;
	double cur;
	double prevVal=mod->ProportionInvariant();
	
	while(1){
		double incr=0.001;
		mod->SetPinv(prevVal+incr);
		
		MakeAllNodesDirty();
		Score();
		cur=lnL;
		double d11=(cur-prev)/incr;
		
		mod->SetPinv(prevVal-incr);
		MakeAllNodesDirty();
		Score();
		cur=lnL;
		double d12=(cur-prev)/-incr;
		
		double d1=(d11+d12)*.5;
		double d2=(d11-d12)/incr;
		
		double est=-d1/d2;
		
		if((abs(est) < 0.001 && d2 < 0.0) || (d1>0.0 && prevVal==mod->MaxPinv())){
			mod->SetPinv(prevVal);
			MakeAllNodesDirty();			
			return prev-start;
			}
		
		double t;
		if(d2 < 0.0){
			t=prevVal+est;
			if(t < 0.0) t=prevVal*.5;
			//if(t > mod->MaxPinv()) t=(prevVal+mod->MaxPinv())*.5;
			if(t > mod->MaxPinv()) t=mod->MaxPinv();
			}
		else{
			if(d1 > 0.0) t=prevVal*1.1;
			else t=prevVal*0.9;
	//		if(t > mod->MaxPinv()) t=(prevVal+mod->MaxPinv())*.5;
			if(t > mod->MaxPinv()) t=mod->MaxPinv();
			}
			
		mod->SetPinv(t);
		MakeAllNodesDirty();			
		Score();
		prev=lnL;
		prevVal=t;	
		}
	return -1;
	}

double Tree::OptimizeBranchLength(double optPrecision, TreeNode *nd, bool goodGuess){
	double improve;
	
#ifdef OPT_DEBUG
	optsum << nd->nodeNum << "\t" << nd->dlen << "\t";
#endif
		
#ifdef BRENT
	improve = BrentOptimizeBranchLength(optPrecision, nd, goodGuess);
#else
	improve = NewtonRaphsonOptimizeBranchLength(optPrecision, nd, goodGuess);
#endif

#ifdef OPT_DEBUG
	optsum << nd->dlen << "\t" << improve << endl;
#endif

	return improve;
	}

void Tree::OptimizeBranchesWithinRadius(TreeNode *nd, double optPrecision, int subtreeNode, TreeNode *prune/*=NULL*/){
	nodeOptVector.clear();
	
	
	double totalIncrease=0.0, prunePointIncrease=0.0, thisIncr;
	totalIncrease += OptimizeBranchLength(optPrecision, nd->left, false);
	totalIncrease += OptimizeBranchLength(optPrecision, nd, false);
	totalIncrease += OptimizeBranchLength(optPrecision, nd->right, false);	

	Score(nd->nodeNum);
	if(lnL < globalBest - treeRejectionThreshold){
		#ifdef OPT_DEBUG
		optsum << "3 branch total\t" << totalIncrease << endl;
		optsum << "bailing early\t" << globalBest-lnL << endl;
		#endif
		return;
		}
	
	nodeOptVector.push_back(nd->left);
	nodeOptVector.push_back(nd);
	nodeOptVector.push_back(nd->right);

	if(prune!=NULL){
		thisIncr = OptimizeBranchLength(optPrecision, prune->left, true);
		//if(thisIncr > optPrecision) nodeOptVector.push_back(prune->left);
		nodeOptVector.push_back(prune->left);
		prunePointIncrease += thisIncr;
		if(prune->anc==NULL){
			thisIncr = OptimizeBranchLength(optPrecision, prune->left->next, true);
			//if(thisIncr > optPrecision) nodeOptVector.push_back(prune);
			nodeOptVector.push_back(prune->left->next);
			prunePointIncrease += thisIncr;
			}
		else{
			thisIncr = OptimizeBranchLength(optPrecision, prune, true);
			//if(thisIncr > optPrecision) nodeOptVector.push_back(prune);
			nodeOptVector.push_back(prune);
			prunePointIncrease += thisIncr;
			}	
		thisIncr = OptimizeBranchLength(optPrecision, prune->right, true);
		//if(thisIncr > optPrecision) nodeOptVector.push_back(prune->right);
		nodeOptVector.push_back(prune->right);
		prunePointIncrease += thisIncr;
		totalIncrease+=prunePointIncrease;
		}

	//now spread out from there
	int rad=10;

	if(rad>0){
		if(nd->right->left!=NULL) totalIncrease += RecursivelyOptimizeBranches(nd->right->left, optPrecision, subtreeNode, rad, false, 0.0);
		if(nd->left->left!=NULL) totalIncrease += RecursivelyOptimizeBranches(nd->left->left, optPrecision, subtreeNode, rad, false, 0.0);
		if(nd->anc!=root) totalIncrease += RecursivelyOptimizeBranchesDown(nd->anc, nd, optPrecision, subtreeNode, rad, 0.0);
		else{
			if(root->left != nd) totalIncrease +=  RecursivelyOptimizeBranches(root->left, optPrecision, subtreeNode, rad, true, 0.0);
			if(root->left->next != nd) totalIncrease +=  RecursivelyOptimizeBranches(root->left->next, optPrecision, subtreeNode, rad, true, 0.0);
			if(root->right != nd) totalIncrease +=  RecursivelyOptimizeBranches(root->right, optPrecision, subtreeNode, rad, true, 0.0);
			}
		
		if(prunePointIncrease > 0.0 ){
			if(prune->right->left!=NULL) totalIncrease += RecursivelyOptimizeBranches(prune->right->left, optPrecision, subtreeNode, rad, false, 0.0);
			if(prune->left->left!=NULL) totalIncrease += RecursivelyOptimizeBranches(prune->left->left, optPrecision, subtreeNode, rad, false, 0.0);
			if(prune==root){
				if(prune->left->next->left!=NULL) totalIncrease += RecursivelyOptimizeBranches(prune->left->next->left, optPrecision, subtreeNode, rad, false, 0.0);
				}
			else{
				if(prune->anc!=root) totalIncrease += RecursivelyOptimizeBranchesDown(prune->anc, prune, optPrecision, subtreeNode, rad, 0.0);
				else{
					if(prune->nodeNum!=subtreeNode){
						if(root->left != prune) totalIncrease +=  RecursivelyOptimizeBranches(root->left, optPrecision, subtreeNode, rad, true, 0.0);
						if(root->left->next != prune) totalIncrease +=  RecursivelyOptimizeBranches(root->left->next, optPrecision, subtreeNode, rad, true, 0.0);
						if(root->right != prune) totalIncrease +=  RecursivelyOptimizeBranches(root->right, optPrecision, subtreeNode, rad, true, 0.0);
						}
					}
				}
			}
		}
#ifdef OPT_DEBUG
optsum << "radiusopt intial pass total\t" << totalIncrease << endl;
#endif

	double postIncr=0.0;
	TreeNode *finalNode;
	
		if(nodeOptVector.empty()==false && prune != NULL){//remove any duplicate entries caused by an overlaping radius around connector and prune
			for(list<TreeNode*>::iterator fir=nodeOptVector.begin();fir!=nodeOptVector.end();fir++){
				list<TreeNode*>::iterator sec=fir;
				sec++;
				if(sec==nodeOptVector.end()) break;
				for(;sec!=nodeOptVector.end();){
					if(*fir == *sec){
						list<TreeNode*>::iterator del=sec;
						if(sec != nodeOptVector.end()) sec++;
						nodeOptVector.erase(del);
						}
					else 
						if(sec != nodeOptVector.end()) sec++;
					}
				}
			}
	while(nodeOptVector.empty() == false){
		list<TreeNode*>::iterator it=nodeOptVector.begin();
		while(it!=nodeOptVector.end()){
			if(nodeOptVector.size() == 1) finalNode=*it;
			thisIncr=OptimizeBranchLength(optPrecision, *(it), true);
			postIncr+= thisIncr;
			if(!(thisIncr > optPrecision)){
				list<TreeNode*>::iterator del=it;
				it++;
				nodeOptVector.erase(del);
				}
			else it++;
			}
		}
	if(finalNode != root)
		Score(finalNode->anc->nodeNum);
	else Score(0);


#ifdef OPT_DEBUG
optsum << "postopt total\t" << postIncr << endl;
#endif
	}

double Tree::RecursivelyOptimizeBranches(TreeNode *nd, double optPrecision, int subtreeNode, int radius, bool dontGoNext, double scoreIncrease, bool ignoreDelta/*=false*/){
	double delta = OptimizeBranchLength(optPrecision, nd, true);
	scoreIncrease += delta;
	
	if(!(delta < optPrecision))
		nodeOptVector.push_back(nd);
//	if(radius==0) cout << "hit max radius!" <<endl;
	if(nd->left!=NULL && radius>1 && (!(delta < optPrecision) || ignoreDelta == true)){
		scoreIncrease += RecursivelyOptimizeBranches(nd->left, optPrecision, subtreeNode, radius-1, false, 0, ignoreDelta);
		}
	if(nd->next!=NULL && dontGoNext==false){
		scoreIncrease += RecursivelyOptimizeBranches(nd->next, optPrecision, subtreeNode, radius, false, 0, ignoreDelta);
		}
	if(memLevel > 1) RemoveTempClaReservations();
	return scoreIncrease;
	}

double Tree::RecursivelyOptimizeBranchesDown(TreeNode *nd, TreeNode *calledFrom, double optPrecision, int subtreeNode, int radius, double scoreIncrease){
		
	double delta = OptimizeBranchLength(optPrecision, nd, true);
	scoreIncrease += delta;

	if(nd->nodeNum == subtreeNode){
		return scoreIncrease;
		}

	if(!(delta < optPrecision))
		nodeOptVector.push_back(nd);
//	if(radius==0) cout << "hit max radius!" <<endl;
	if(nd->left!=NULL && nd->left!=calledFrom && radius>1) scoreIncrease += RecursivelyOptimizeBranches(nd->left, optPrecision, subtreeNode, radius, true, 0);
	else if(radius>1) scoreIncrease += RecursivelyOptimizeBranches(nd->left->next, optPrecision, subtreeNode, radius, false, 0);
	if(nd->anc!=root && radius>1 && !(delta < optPrecision)){
		scoreIncrease += RecursivelyOptimizeBranchesDown(nd->anc, nd, optPrecision, subtreeNode, radius-1, 0);
		}
	if(nd->anc==root){
		if(radius>1 && !(delta < optPrecision)){
			if(nd->next!=NULL) scoreIncrease += RecursivelyOptimizeBranches(nd->next, optPrecision, subtreeNode, radius-1, true, 0);
			else scoreIncrease += RecursivelyOptimizeBranches(nd->prev->prev, optPrecision, subtreeNode, radius-1, true, 0);
			if(nd->prev!=NULL) scoreIncrease += RecursivelyOptimizeBranches(nd->prev, optPrecision, subtreeNode, radius-1, true, 0);
			else scoreIncrease += RecursivelyOptimizeBranches(nd->next->next, optPrecision, subtreeNode, radius-1, true, 0);
			}
		}
	if(memLevel > 1) RemoveTempClaReservations();
	return scoreIncrease;
	}

pair<double, double> Tree::CalcDerivativesRateHet(TreeNode *nd1, TreeNode *nd2){
	//nd1 and nd2 are the nodes on either side of the branch of interest
	//nd1 will always be the "lower" one, and will always be internal, while
	//nd2 can be internal or terminal
	int nsites=data->NChar();
	
//	CondLikeArray *destCLA=claMan->GetTempCla(0);
//	CondLikeArray *destDeriv1=claMan->GetTempCla(1);
//	CondLikeArray *destDeriv2=claMan->GetTempCla(2);
	
	const CondLikeArray *claOne;
	if(nd1->left == nd2)
		claOne=GetClaUpLeft(nd1, true);
	else if(nd1->right == nd2)
		claOne=GetClaUpRight(nd1, true);
	else //nd1 must be the root, and nd2 it's middle des
		claOne=GetClaDown(nd1, true);
	
	//this must happen BEFORE the derivs are calced, or the prmat won't be current for this branch!
	CondLikeArray *claTwo=NULL;
	if(nd2->left != NULL)
		claTwo=GetClaDown(nd2, true);
			
	double ***deriv1, ***deriv2, ***prmat;
	
	mod->CalcDerivatives(nd2->dlen, prmat, deriv1, deriv2);

	double d1=0.0, d2=0.0;

	if(nd2->left == NULL){
		const char *childData=nd2->tipData;
		GetDerivsPartialTerminalFlexRates(claOne, **prmat, **deriv1, **deriv2, childData, d1, d2);
		}
	else {
		GetDerivsPartialInternalFlexRates(claOne, claTwo, **prmat, **deriv1, **deriv2, d1, d2);
		}

	return pair<double, double>(d1, d2);
}

double Tree::BranchLike(TreeNode *optNode){

	bool scoreOK=true;
	do{
		try{
			if(optNode->anc->left==optNode){
				optNode->anc->claIndexDown = claMan->SetDirty(optNode->anc->claIndexDown);
				optNode->anc->claIndexUR = claMan->SetDirty(optNode->anc->claIndexUR);		
				GetClaUpLeft(optNode->anc);
				}
			else if(optNode->anc->right==optNode){
				optNode->anc->claIndexDown = claMan->SetDirty(optNode->anc->claIndexDown);
				optNode->anc->claIndexUL = claMan->SetDirty(optNode->anc->claIndexUL);
				GetClaUpRight(optNode->anc);
				}
			else {
				optNode->anc->claIndexUL = claMan->SetDirty(optNode->anc->claIndexUL);
				optNode->anc->claIndexUR = claMan->SetDirty(optNode->anc->claIndexUR);
				GetClaDown(optNode->anc);
				}
			
			//now sum as if this were the root
			if(mod->NRateCats()>1) ConditionalLikelihoodRateHet(ROOT, optNode->anc);
			else ConditionalLikelihood(ROOT, optNode->anc);
			return lnL;
			}
		catch(int){
			scoreOK=false;
			MakeAllNodesDirty();
			rescaleEvery -= 2;
			ofstream resc("rescale.log", ios::app);
			resc << "rescale reduced to " << rescaleEvery << endl;
			resc.close();
			if(rescaleEvery<2) rescaleEvery=2;
			}			
		}while(scoreOK==false);
	return 0;
	}


void Tree::SampleBlenCurve(TreeNode *nd, ofstream &out){
	
	double initialLen=nd->dlen;
	Score();
	
	out << nd->dlen << "\t" << lnL << "\n";
	
	SetBranchLength(nd, 1e-4);
	for(int i=0;i<15;i++){
		Score();
		out << nd->dlen << "\t" << lnL << "\n";
		SetBranchLength(nd, nd->dlen * 2.0);	
		}
	SetBranchLength(nd, initialLen);	
	} 

 double Tree::NewtonRaphsonOptimizeBranchLength(double precision1, TreeNode *nd, bool goodGuess){

	if(goodGuess==false && (nd->dlen < 0.0001 || nd->dlen > .1)){
		SetBranchLength(nd, .001);
		}

/*if(nd->dlen==DEF_MIN_BRLEN){
//if(nd->nodeNum == 20){
	ofstream scr("NRcurve.log");
	scr.precision(20);
	assert(scr.good());
	scr.precision(15);
	double initDlen = nd->dlen;
	for(double d=1e-8;d<.1;d*=1.33){
		nd->dlen = d;
		SweepDirtynessOverTree(nd);
		Score();
		scr << d << "\t" << lnL << endl;
		}
	nd->dlen=initDlen;
	SweepDirtynessOverTree(nd);
	scr.close();
	}	
*/
//	nd->dlen=.3254;
//	SweepDirtynessOverTree(nd);

	//DEBUG
/*	
	ofstream deb("curves.log");
	SampleBlenCurve(nd, deb);
	deb.close();
*/	

//	MakeAllNodesDirty();
/*
	const double startDLen = nd->dlen;
	double incr;
	if(nd->nodeNum==1){
		int poo=1;
		}
	if(nd->dlen > 1e-4)
		incr=.00001;
	else incr=.0000001;
	
	incr=nd->dlen/1000.0;
	
	Score();
	double start=lnL;

	SetBranchLength(nd, startDLen + incr);
	Score();

	double empD11= (lnL - start)/incr;

//	SetBranchLength(nd, prevDLen);
//	Score();

	SetBranchLength(nd, startDLen - incr);
	Score();

	double empD12 = (lnL - start)/-incr;

	double empD1=(empD11+empD12)*.5;
	double empD2=(empD11-empD12)/incr;

	SetBranchLength(nd, startDLen);
*/
//	MakeAllNodesDirty();

#ifdef OPT_DEBUG
//	ofstream log("optimization.log", ios::app);
//	log.precision(10);

	opt << nd->nodeNum << "\t" << nd->dlen << "\n";
	
//	ofstream scr("impVSd1.log", ios::app);
	if(lnL > -2){
		Score(nd->anc->nodeNum);
		}
	double delta;
	
#endif
	double totalEstImprove=0.0;
	bool continueOpt=false;
	int iter=0;
	double abs_d1_prev=1e200;;
	const double v_onEntry=nd->dlen;
	double v=nd->dlen;
	double v_prev = nd->dlen;				/* in case we don't like the new value (see below) */
	bool moveOn = false;
	double prevScore=lnL;
	double curScore=lnL;
	int negProposalNum=0;

	do{
		bool scoreOK;
		int sweeps=0;
		pair<double, double> derivs;
		do{		//this part just catches the exception that could be thrown by the rescaling 
				//function if it decides that the current rescaleEvery is too large
			try{
				scoreOK=true;
				derivs = CalcDerivativesRateHet(nd->anc, nd);
				}
			catch(int err){
				scoreOK=false;
				if(err==1){
					MakeAllNodesDirty();
					rescaleEvery -= 2;
					ofstream resc("rescale.log", ios::app);
					resc << "rescale reduced to " << rescaleEvery << endl;
					resc.close();
					if(rescaleEvery<2) rescaleEvery=2;
					}
				else if(err==2){
					//this is necessary because rarely it is possible that attempted optimization at nodes
					//across the tree causes more than a single set of clas to be in use, which can cause 
					//clas to run out if we are in certain memory situations
					assert(sweeps==0);
					SweepDirtynessOverTree(nd);
					sweeps++;
					}
				}
			}while(scoreOK==false);	
		
		
		double d1=derivs.first; //* 4.0 * pow(nd->dlen, 0.75);
		double d2=derivs.second;// * 3.0 * pow(nd->dlen, -0.25);

		//DEBUG
/*		if(iter==0 && nd->dlen > 1e-7 && (abs(d1-empD1) > 1.0)){// || (abs(d2-empD2) > 50.0))){

			ofstream bad("badD1.log", ios::app);
			bad << nd->nodeNum << "\t" << nd->dlen << "\t" << d1 << "\t" << d2 << "\t" << empD1 << "\t" << empD2 << "\t";
			mod->OutputGamlFormattedModel(bad);
			bad << "\n";
			bad.close();
			}
*/
		double estDeltaNR=-d1/d2;
		
		//this was my original ad hoc estimated score change, which was always an overestimate
//		double estScoreDelta = .6666 * d1*estDeltaNR;
		
		//estimated change in score by a Taylor series (always an underestimate)
		double estScoreDelta = d1*estDeltaNR + (d2 * estDeltaNR * estDeltaNR * 0.5);
			
		if(iter==0 && estScoreDelta > precision1) continueOpt=true;
		
														#ifdef OPT_DEBUG			
														opt << d1 << "\t" << d2 << "\t" << estScoreDelta << "\t";		
														#endif
		double abs_d1 = fabs(d1);
		if (d2 >= 0.0){
														#ifdef OPT_DEBUG			
														opt << "d2 > 0, try minbrlen?\t";				
														#endif
			//the only time this seems to happen is when the peak is at the min and the
			//current point is far from that
			assert(d1 < 0.0);
//			if(abs(nd->dlen - DEF_MIN_BRLEN) < DEF_MIN_BRLEN){				
//			if(nd->dlen==DEF_MIN_BRLEN){
			if(((DEF_MIN_BRLEN - nd->dlen)*d1) < precision1){
				#ifdef OPT_DEBUG
				opt << "no, return\n";
				#endif
				return totalEstImprove;
				}
			v=DEF_MIN_BRLEN;
			totalEstImprove += precision1;
			}
		else{
			if(estScoreDelta < precision1){
														#ifdef OPT_DEBUG			
														opt << "delta < prec, return\n";
														if(curScore==-1.0){
															Score(nd->anc->nodeNum);
															}		
														#endif
				return totalEstImprove;
				}
			else{/* Take the Newton-Raphson step */
					v += estDeltaNR;
														#ifdef OPT_DEBUG			
														opt << v << "\t";			
														#endif
				}
			if ((iter != 0) && (abs_d1 > abs_d1_prev)){
				//not doing anything special here.  This generally means that we overshot the peak, but
				//should get it from the other side
														#ifdef OPT_DEBUG			
														opt << "d1 increased!\t";	
														#endif
				}
			if (v < DEF_MIN_BRLEN){
				negProposalNum++;
				double deltaToMin=DEF_MIN_BRLEN-nd->dlen;
				double scoreDeltaToMin = (deltaToMin * d1 + (deltaToMin*deltaToMin*d2*.5));
				if(scoreDeltaToMin < precision1){

//				if(nd->dlen == DEF_MIN_BRLEN){
													#ifdef OPT_DEBUG
													assert(curScore != -1.0);		
													opt << nd->dlen << "\t" << curScore  <<"\n";			
													#endif
					//if the blen was already the min and the proposed blen is negative,
					//we don't need to continue sweeping over nodes
					return totalEstImprove;
					}
				else if(v < -1.0){
					//if the proposed blen is this negative, just try the min brlen
					v=DEF_MIN_BRLEN;
//					totalEstImprove += precision1;
					totalEstImprove += scoreDeltaToMin;
													#ifdef OPT_DEBUG
													Score(nd->anc->nodeNum);		
													opt << nd->dlen << "\t" << lnL << "\n";			
													#endif
					}
				else if(negProposalNum==1 && nd->dlen > 1e-4 && v_prev != 1e-4){
					//try a somewhat smaller length before going all the way to the min
					v = nd->dlen * .1;
					double delta=v - nd->dlen;
					totalEstImprove += (delta * d1 + (delta*delta*d2*.5));
//					totalEstImprove += precision1;
					}
				else{
					v = DEF_MIN_BRLEN;
					totalEstImprove += scoreDeltaToMin;
//					totalEstImprove += precision1;
					}
				}
			else if (v > DEF_MAX_BRLEN){
				double deltaToMax=DEF_MAX_BRLEN - nd->dlen;
				double scoreDeltaToMax = (deltaToMax * d1 + (deltaToMax*deltaToMax*d2*.5));
				if(scoreDeltaToMax < precision1) return totalEstImprove;
				else{
					v = DEF_MAX_BRLEN;
					totalEstImprove += scoreDeltaToMax;
	//				totalEstImprove += precision1;
					}
				}
			else totalEstImprove += estScoreDelta;

			abs_d1_prev = abs_d1;
			}
		assert(v >= DEF_MIN_BRLEN);
		SetBranchLength(nd, v);
#ifdef OPT_DEBUG
		Score(nd->anc->nodeNum);
		
		if(curScore != -1.0){
			if(lnL < curScore){
				cout << lnL << "\t" << curScore << endl;
				if(curScore - lnL < .005){
					//don't want to have different logic when OPT_DEBUG is on
/*					SetBranchLength(nd, v_prev);
					Score(nd->anc->nodeNum);
					return lnL;
*/					}
				else {//assert(0);
					double poo=lnL;
					SetBranchLength(nd, v_prev);
					MakeAllNodesDirty();
					Score(nd->anc->nodeNum);
//					assert(fabs(prevScore - lnL) < .0001);
					
					SetBranchLength(nd, v);
					MakeAllNodesDirty();
					Score(nd->anc->nodeNum);
//					assert(fabs(poo - lnL) < .0001);
					}
				}
			}
				
		curScore=lnL;	
		delta=prevScore - lnL;
opt.precision(7);
opt << v << "\t" << lnL << "\n";
#endif
		prevScore=lnL;
		v_prev=v;
		
		iter++;
		if(iter>50){
			ofstream deb("optdeb.log");
			deb << "initial length " << v_onEntry << endl;
			deb << "current length " << nd->dlen << endl;
			deb << "prev length " << v_prev << endl;
			deb << "d1 " << d1 << " d2 " << d2 << endl;
			deb << "neg proposal num " << negProposalNum << endl;			
			deb.close();
			assert(iter<=50);
			}		
		}while(moveOn==false);
#ifdef OPT_DEBUG
	opt << "final\t" << nd->dlen << "\t" << lnL << endl;
#endif
	assert(nd->dlen > 0.0);
//	assert(curScore != -1.0);
	return continueOpt;
	}
/*
void Tree::RecursivelyOptimizeBranches(TreeNode *nd, double optPrecision, int subtreeNode, int radius, int centerNode, bool dontGoNext){
	double prevScore=lnL;
#ifdef BRENT	
	BrentOptimizeBranchLength(optPrecision, nd, false);
	double delta=lnL - prevScore;
	bool continueOpt=(delta*2.0 > optPrecision ? true : false);
#else
	bool continueOpt = NewtonRaphsonOptimizeBranchLength(optPrecision, nd);
//	continueOpt=true;
#endif
	
	if(nd->left!=NULL && radius>1 && continueOpt) RecursivelyOptimizeBranches(nd->left, optPrecision, subtreeNode, radius-1, centerNode, false);
	if(nd->next!=NULL && dontGoNext==false){
		RecursivelyOptimizeBranches(nd->next, optPrecision, subtreeNode, radius, centerNode, false);
		}
	}

void Tree::RecursivelyOptimizeBranchesDown(TreeNode *nd, TreeNode *calledFrom, double optPrecision, int subtreeNode, int radius, int ){
	double prevScore=lnL;
#ifdef BRENT	
	BrentOptimizeBranchLength(optPrecision, nd, false);
	double delta=lnL - prevScore;
	bool continueOpt=(delta*2.0 > optPrecision ? true : false); 
#else
	bool continueOpt = NewtonRaphsonOptimizeBranchLength(optPrecision, nd);
//	continueOpt=true;
#endif
	
	if(nd->left!=NULL && nd->left!=calledFrom && radius>1) RecursivelyOptimizeBranches(nd->left, optPrecision, subtreeNode, radius, 0, true);
	else if(radius>1) RecursivelyOptimizeBranches(nd->left->next, optPrecision, subtreeNode, radius, 0, false);
	if(nd->anc!=root && radius>1 && continueOpt){
		RecursivelyOptimizeBranchesDown(nd->anc, nd, optPrecision, subtreeNode, radius-1, 0);
		}
	if(nd->anc==root){
		if(radius>1 && continueOpt){
			if(nd->next!=NULL) RecursivelyOptimizeBranches(nd->next, optPrecision, subtreeNode, radius-1, 0, true);
			else RecursivelyOptimizeBranches(nd->prev->prev, optPrecision, subtreeNode, radius-1, 0, true);
			if(nd->prev!=NULL) RecursivelyOptimizeBranches(nd->prev, optPrecision, subtreeNode, radius-1, 0, true);
			else RecursivelyOptimizeBranches(nd->next->next, optPrecision, subtreeNode, radius-1, 0, true);
			}
		}
	}
*/
/*
void Tree::OptimizeBranchesAroundNode(TreeNode *nd, double optPrecision, int subtreeNode){
	//this function will optimize the three branches (2 descendents and one anc) connecting
	//to it.  It assumes that everything that is dirty has been marked so.
	//by default there is only a single optimization pass over the three nodes
	double precision1, precision2;

	if(subtreeNode==0) SetAllTempClasDirty();
	
	precision1=optPrecision;// * 0.5;
	if(optPrecision > .2) precision2=0.0;
	else precision2=precision1 * 0.5;
	
	if(nd != root){
		BrentOptimizeBranchLength(precision1, nd, false);
		BrentOptimizeBranchLength(precision1, nd->left, false);
		BrentOptimizeBranchLength(precision1, nd->right, false);
		}
	else{
		BrentOptimizeBranchLength(precision1, nd->left, false);
		BrentOptimizeBranchLength(precision1, nd->left->next, false);
		BrentOptimizeBranchLength(precision1, nd->right, false);	
		}

	
	if(precision2 > 0){
		//if were're doing multiple optimization passes, only this stuff needs to be set dirty
		claMan->SetDirty(nd->nodeNum, nd->claIndex, true);
		claMan->SetTempDirty(nd->nodeNum, true);
		if(nd != root) claMan->SetTempDirty(nd->anc->nodeNum, true);

		if(nd != root){
			BrentOptimizeBranchLength(precision2, nd, false);
			BrentOptimizeBranchLength(precision2, nd->left, false);
			BrentOptimizeBranchLength(precision2, nd->right, false);
			}
		else {
			BrentOptimizeBranchLength(precision2, nd->left, false);
			BrentOptimizeBranchLength(precision2, nd->left->next, false);
			BrentOptimizeBranchLength(precision2, nd->right, false);			
			}
		}
		
	//these must be called after all optimization passes are done around this node
	TraceDirtynessToRoot(nd);
	if(subtreeNode==0)
		SetAllTempClasDirty();
	else SetTempClasDirtyWithinSubtree(subtreeNode);
	}
*/

/*
inline double CallBranchLike(TreeNode *thisnode, Tree *thistree, double blen){
	thisnode->dlen=exp(blen);
	return thistree->BranchLike(thisnode)*-1;
	}
	
inline double CallBranchLikeRateHet(TreeNode *thisnode, Tree *thistree, double blen){

	thisnode->dlen=blen;
	double like=thistree->BranchLikeRateHet(thisnode)*-1;

#ifdef OPT_DEBUG
	ofstream opt("optimization.log" ,ios::app);
	opt.precision(11);	
	opt << thisnode->dlen << "\t" << like << "\n\t";
	opt.close();

	ofstream opttrees("opttrees.tre", ios::app);
		char treeString[20000];
		thistree->root->MakeNewick(treeString, false);
		opttrees <<  "utree tree1=" << treeString << ";" << endl;
		opttrees.close();
	//if(thisnode->left!=NULL) thistree->TraceDirtynessToRoot(thisnode);
	
	ofstream scr("optscores.log", ios::app);
	scr.precision(10);
	scr << like << "\t" << blen << endl;
	scr.close(); 
#endif

	thistree->RerootHere(thisnode->nodeNum);
	thistree->MakeAllNodesDirty();
	thistree->Score(thistree->data);
	
	return like;
	}

double Tree::BrentOptimizeBranchLength(double accuracy_cutoff, TreeNode *here, bool firstPass){
	//we pass the node whose branch length whose blen we want to optimize, but note that the
	//calculations occur at the node below that
	//if firstPass is true, we have no idea what a reasonable value for the blen is, so
	//use a wide bracket.  If it is false, try a fairly tight bracket around the current val
	double a, b, c, fa, fb, fc, minimum, minScore=0.0;
	double blen=here->dlen;

	assert(blen>=DEF_MIN_BRLEN);

	if(here->anc){
		if(firstPass){
			if(blen<1e-6){
				a=DEF_MIN_BRLEN;
				if(blen!=DEF_MIN_BRLEN){
					b=blen;
					}
				else{
					b=DEF_MIN_BRLEN*100;
					lnL=-1;
					}
				c=DEF_MIN_BRLEN*10000.0;
				}				
			else{
				if(blen<.0001){
					a=.000001;
					b=blen;
					c=.01;
					}
				else if(blen<.1){
					a=.0001;
					b=blen;
					c=.1;
					}
				else {
					a=.1;
					b=blen;
					c=.75;					
					}

		}
			}
		else{
			//tighter 
			if(blen > 1e-6){
				a=blen*.66;
				b=blen;
				c=blen*1.5;
				}
			else{
				a=DEF_MIN_BRLEN;
				if(blen!=DEF_MIN_BRLEN){
					b=blen;
					}
				else{
					b=DEF_MIN_BRLEN*100;
					lnL=-1;
					}
				c=DEF_MIN_BRLEN*10000.0;				
				}			
			}

#ifdef OPT_DEBUG
		ofstream opt("optimization.log" ,ios::app);
		opt << "node " << here->nodeNum << "\t" << here->dlen  << "\n";
//		opt << "\t" << a << "\t" << b << "\t" << c << "\n";
#endif
			
		if(mod->NRateCats()==1){
			mnbrak(&a, &b, &c, &fa, &fb, &fc, CallBranchLike, here, this);
	//		opt << a << "\t" << b << "\t" << c << "\t";
			brent(a, b, c, CallBranchLike, accuracy_cutoff, &minimum, here, this);
			}
		else{
#ifdef OPT_DEBUG
			opt << "brak\t";
			opt.close();
#endif
			fb=lnL;
			int zeroMLE = DZbrak(&a, &b, &c, &fa, &fb, &fc, CallBranchLikeRateHet, here, this);
			
			bool flatSurface=false;
			if(fa-fb + fc-fb < .000001) flatSurface=true;
			//braka=fa;
			//brakb=fb;
#ifdef OPT_DEBUG
			ofstream opt("optimization.log" ,ios::app);
//			opt << "bracket\t" << a << "\t" << fa << "\n\t" << b << "\t" << fb << "\n\t" << c << "\t" << fc << endl;
			opt << "brent\t";
			opt.close();
#endif
			if(zeroMLE==0 && flatSurface==false) //if the bracket suggests that the MLE is very near 0, don't bother calling brent
				minScore=DZbrent(a, b, c, fa, fb, fc, CallBranchLikeRateHet, accuracy_cutoff, &minimum, here, this);
			else if(zeroMLE==1){
				minimum=(min_brlen);
				if(a==min_brlen) minScore=fa;
				else if(c==min_brlen) minScore=fc;
				else minScore=-1;
				}
			else{
				minimum=b;
				minScore=fb;
				}
			}
		double min_len=minimum;		
//		double min_len=exp(minimum);
		here->dlen = (min_len > min_brlen ? (min_len < max_brlen ? min_len : max_brlen) : min_brlen);
		
#ifdef OPT_DEBUG
		opt.open("optimization.log" ,ios::app);
		opt.precision(9);
		opt << "final " << "\t" << minScore << "\t" << here->dlen << "\n";
		opt.close();
#endif
	//	claMan->SetTempDirty(-1, true);

/*		MakeAllNodesDirty();
		SetAllTempClasDirty();
	if(minScore!=0.0){
//		TraceDirtynessToRoot(here);
		Score(Tree::data);
		assert(abs(lnL+minScore) <.001);
		}
	}

	SweepDirtynessOverTree(here);
	lnL=minScore;
	return minScore;
	}
*/	

double Tree::BrentOptimizeBranchLength(double accuracy_cutoff, TreeNode *here, bool goodGuess){
	//we pass the node whose branch length whose blen we want to optimize, but note that the
	//calculations occur at the node below that
	//if firstPass is true, we have no idea what a reasonable value for the blen is, so
	//use a wide bracket.  If it is false, try a fairly tight bracket around the current val
	GAMLfloat a, b, c, fa, fb, fc, minimum, minScore=0.0;
	double blen=here->dlen;
	
	double min_len;
	
	assert(blen>=DEF_MIN_BRLEN);

	double initialScore;
	fb=initialScore=CallBranchLike(here, this, sqrt(sqrt(here->dlen)), true);

if(here->anc){
#ifndef FOURTH_ROOT
//	if(here->anc){
		if(firstPass){
			if(!(blen>1e-6)){
				a=DEF_MIN_BRLEN;
				if(blen!=DEF_MIN_BRLEN){
					b=blen;
					}
				else{
					b=DEF_MIN_BRLEN*100;
					lnL=-1;
					}
				c=DEF_MIN_BRLEN*10000.0;
				}				
			else{
				if((blen>0.0001)){
					a=.000001;
					b=blen;
					c=.01;
					}
				else if(!(blen>0.1)){
					a=.0001;
					b=blen;
		//			c=blen*2.0;
					c=blen*16.0;
					}
				else {
					a=.01;
					b=blen;
					c=blen*2.0;
					//c=.75;					
					}
				}
			}
		else{
			//tighter 
			if(blen >= 1e-6){
				a=blen*.66;
				b=blen;
				c=blen*1.5;
				}
			else{
				a=DEF_MIN_BRLEN;
				if(blen!=DEF_MIN_BRLEN){
					b=blen;
					}
				else{
					b=DEF_MIN_BRLEN*100;
					lnL=-1;
					}
				c=DEF_MIN_BRLEN*10000.0;				
				}			
			}
#endif
#ifdef FOURTH_ROOT
	if(blen < DEF_MIN_BRLEN*10){
		a=.01;
		b=a+.05;
		fb=-1;
		c=b+.05;
		}
	else{
		b=sqrt(sqrt(blen));
		if(goodGuess==false){
			a=(b <= 0.06 ? .01 : b-0.05);
			c=b+0.05;
			}
		else{
			a=(b <= 0.026 ? .01 : b-0.025);
			c=b+0.025;		
			}
		}

#elif ROOT_OPT
	a=sqrt(a);
	b=sqrt(b);
	c=sqrt(c);
#endif

#ifdef OPT_DEBUG
	optInfo.Setup(here->nodeNum, blen, accuracy_cutoff, goodGuess, a, b, c);
//	SampleBranchLengthCurve(CallBranchLike, here, this);
	optInfo.Report(curves);
	bool trueMin=optInfo.IsMinAtMinAllowableLength();
	curves.flush();
	optInfo.Setup(here->nodeNum, blen, accuracy_cutoff, goodGuess, a, b, c);

#endif		
		int zeroMLE = DZbrak(&a, &b, &c, &fa, &fb, &fc, CallBranchLike, here, this);
		
#ifdef OPT_DEBUG
/*		if(trueMin != zeroMLE){
			assert(0);
			} */
#endif
		bool flatSurface=false;
		if(fa-fb + fc-fb < .000001){
			flatSurface=true;
			}

		if(zeroMLE==0 && flatSurface==false) //if the bracket suggests that the MLE is very near 0, don't bother calling brent
			minScore=DZbrent(a, b, c, fa, fb, fc, CallBranchLike, accuracy_cutoff, &minimum, here, this);
		else if(zeroMLE==1){
#ifdef FOURTH_ROOT
			assert(c==effectiveMin);
			minimum=c;
			minScore=fc;
#elif ROOT_OPT
			double sqrtmin=sqrt(min_brlen);
			minimum=sqrtmin;
			if(a==sqrtmin) minScore=fa;
			else if(c==sqrtmin) minScore=fc;
			else minScore=-1;
#else
			minimum=(min_brlen);
			if(a==min_brlen) minScore=fa;
			else if(c==min_brlen) minScore=fc;
			else minScore=-1;
#endif
			}
		else{
			minimum=b;
			minScore=fb;
			}

#ifdef FOURTH_ROOT
		min_len=minimum*minimum*minimum*minimum;
#elif ROOT_OPT
		if(zeroMLE)
			min_len=minimum;
		else 
			min_len=minimum*minimum;
#else
		min_len=minimum;
#endif	
		}

//	if(here->dlen != min_len){
		here->dlen = (min_len > min_brlen ? (min_len < max_brlen ? min_len : max_brlen) : min_brlen);
		SweepDirtynessOverTree(here);
//		}
	assert(minScore!=-1);
/*	if(minScore == -1){
		minScore=CallBranchLike(here, this, here->dlen, false);
		}
*/	lnL=-minScore;
	
	#ifdef OPT_DEBUG
	optInfo.Report(opt);
	opt << "final\t" << minimum << "\t" << minScore << endl;
	
//	optsum << here->nodeNum << "\t" << blen << "\t" << min_len << "\t" << initialScore - minScore << endl;
	
	#endif
	
	return initialScore - minScore;
	}

void Tree::GetDerivsPartialTerminalFlexRates(const CondLikeArray *partialCLA, const double *prmat, const double *d1mat, const double *d2mat, const char *Ldata, double &d1Tot, double &d2Tot){

	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.
	double *partial=partialCLA->arr;

	int nchar=data->NChar();

#ifdef UNIX
	madvise(partial, nchar*16*sizeof(double), MADV_SEQUENTIAL);
#endif

	double siteL, tempD1, D1tot, D2tot;
	D1tot=D2tot=0.0;
	double La, Lc, Lg, Lt;
	double D1a, D1c, D1g, D1t;
	double D2a, D2c, D2g, D2t;
	
	const int *countit=data->GetCounts();
	const double *rateProb=mod->GetRateProbs();
	const int nRateCats=mod->NRateCats();
	
	int lastConst=data->LastConstant();
	const int *conBases=data->GetConstBases();
	double prI=mod->ProportionInvariant();

	for(int i=0;i<nchar;i++){
		if((*countit) != 0){
			La=Lc=Lg=Lt=D1a=D1c=D1g=D1t=D2a=D2c=D2g=D2t=0.0;
			if(*Ldata > -1){ //no ambiguity
				for(int i=0;i<nRateCats;i++){
					La  += prmat[(*Ldata)+16*i] * partial[0] * rateProb[i];
					D1a += d1mat[(*Ldata)+16*i] * partial[0] * rateProb[i];
					D2a += d2mat[(*Ldata)+16*i] * partial[0] * rateProb[i];
					Lc  += prmat[(*Ldata+4)+16*i] * partial[1]* rateProb[i];
					D1c += d1mat[(*Ldata+4)+16*i] * partial[1]* rateProb[i];
					D2c += d2mat[(*Ldata+4)+16*i] * partial[1]* rateProb[i];
					Lg  += prmat[(*Ldata+8)+16*i] * partial[2]* rateProb[i];
					D1g += d1mat[(*Ldata+8)+16*i] * partial[2]* rateProb[i];
					D2g += d2mat[(*Ldata+8)+16*i] * partial[2]* rateProb[i];
					Lt  += prmat[(*Ldata+12)+16*i] * partial[3]* rateProb[i];
					D1t += d1mat[(*Ldata+12)+16*i] * partial[3]* rateProb[i];
					D2t += d2mat[(*Ldata+12)+16*i] * partial[3]* rateProb[i];
					partial += 4;
					}
				Ldata++;
				}
				
			else if(*Ldata == -4){ //total ambiguity
				for(int i=0;i<nRateCats;i++){
					La += partial[0]* rateProb[i];
					Lc += partial[1]* rateProb[i];
					Lg += partial[2]* rateProb[i];
					Lt += partial[3]* rateProb[i];
					partial += 4;
					}
				Ldata++;
				}
			else{ //partial ambiguity
				char nstates=-1 * *(Ldata++);
				for(int i=0;i<nstates;i++){
					for(int i=0;i<nRateCats;i++){
						La += prmat[(*Ldata)+16*i]  * partial[4*i]* rateProb[i];
						D1a += d1mat[(*Ldata)+16*i] * partial[4*i]* rateProb[i];		
						D2a += d2mat[(*Ldata)+16*i] * partial[4*i]* rateProb[i];
											
						Lc += prmat[(*Ldata+4)+16*i] * partial[1+4*i]* rateProb[i];
						D1c += d1mat[(*Ldata+4)+16*i]* partial[1+4*i]* rateProb[i];
						D2c += d2mat[(*Ldata+4)+16*i]* partial[1+4*i]* rateProb[i];
											
						Lg += prmat[(*Ldata+8)+16*i]* partial[2+4*i]* rateProb[i];
						D1g += d1mat[(*Ldata+8)+16*i]* partial[2+4*i]* rateProb[i];
						D2g += d2mat[(*Ldata+8)+16*i]* partial[2+4*i]* rateProb[i];
						
						Lt += prmat[(*Ldata+12)+16*i]* partial[3+4*i]* rateProb[i];
						D1t += d1mat[(*Ldata+12)+16*i]* partial[3+4*i]* rateProb[i];
						D2t += d2mat[(*Ldata+12)+16*i]* partial[3+4*i]* rateProb[i];
						}
					Ldata++;
					}
				partial+=4*nRateCats;
				}
			if((mod->NoPinvInModel() == false) && (i<=lastConst)){
				double btot=0.0;
				if(conBases[i]&1) btot+=mod->Pi(0);
				if(conBases[i]&2) btot+=mod->Pi(1);
				if(conBases[i]&4) btot+=mod->Pi(2);
				if(conBases[i]&8) btot+=mod->Pi(3);
				//6-27-05 fixed this to calc derivs correctly if constant site has been rescaled
				siteL  = ((La*mod->Pi(0)+Lc*mod->Pi(1)+Lg*mod->Pi(2)+Lt*mod->Pi(3)) + (prI*btot)*exp((double)partialCLA->underflow_mult[i]));
				}
			else
				siteL  = ((La*mod->Pi(0)+Lc*mod->Pi(1)+Lg*mod->Pi(2)+Lt*mod->Pi(3)));

			tempD1 = (((D1a*mod->Pi(0)+D1c*mod->Pi(1)+D1g*mod->Pi(2)+D1t*mod->Pi(3))) / siteL);
			d1Tot += *countit * tempD1;
			assert(d1Tot == d1Tot);
			double siteD2=((D2a*mod->Pi(0)+D2c*mod->Pi(1)+D2g*mod->Pi(2)+D2t*mod->Pi(3)));
			d2Tot += *countit * ((siteD2 / siteL) - tempD1*tempD1);
			assert(d2Tot == d2Tot);
			}
		else{
			partial+=4*nRateCats;
			if(!(*Ldata < 0)) Ldata++;
			else if(*Ldata == -4) Ldata++;
			else{
				char nstates=-1 * *(Ldata++);
				for(int i=0;i<nstates;i++) Ldata++;
				}
			}
		countit++;
		}
	}

void Tree::GetDerivsPartialTerminalRateHet(const CondLikeArray *partialCLA, const double *prmat, const double *d1mat, const double *d2mat, const char *Ldata, double &d1Tot, double &d2Tot){

	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.
	double *partial=partialCLA->arr;

	int nchar=data->NChar();

#ifdef UNIX
	madvise(partial, nchar*4*mod->NRateCats()*sizeof(double), MADV_SEQUENTIAL);
#endif

	double siteL, tempD1, D1tot, D2tot;
	D1tot=D2tot=0.0;
	double La, Lc, Lg, Lt;
	double D1a, D1c, D1g, D1t;
	double D2a, D2c, D2g, D2t;
	
	//gamma and invariants
	const int *countit=data->GetCounts();
	int lastConst=data->LastConstant();
	const int *conBases=data->GetConstBases();
	double prI=mod->ProportionInvariant();
	int nRateCats=mod->NRateCats();
	double scaledGammaProp=(1.0-prI) / nRateCats;
	
	for(int i=0;i<nchar;i++){
		if((*countit) != 0){
			La=Lc=Lg=Lt=D1a=D1c=D1g=D1t=D2a=D2c=D2g=D2t=0.0;
			if(*Ldata > -1){ //no ambiguity
				for(int i=0;i<nRateCats;i++){
					La  += prmat[(*Ldata)+16*i] * partial[0];
					D1a += d1mat[(*Ldata)+16*i] * partial[0];
					D2a += d2mat[(*Ldata)+16*i] * partial[0];
					Lc  += prmat[(*Ldata+4)+16*i] * partial[1];
					D1c += d1mat[(*Ldata+4)+16*i] * partial[1];
					D2c += d2mat[(*Ldata+4)+16*i] * partial[1];
					Lg  += prmat[(*Ldata+8)+16*i] * partial[2];
					D1g += d1mat[(*Ldata+8)+16*i] * partial[2];
					D2g += d2mat[(*Ldata+8)+16*i] * partial[2];
					Lt  += prmat[(*Ldata+12)+16*i] * partial[3];
					D1t += d1mat[(*Ldata+12)+16*i] * partial[3];
					D2t += d2mat[(*Ldata+12)+16*i] * partial[3];
					partial += 4;
					}
				Ldata++;
				}
				
			else if(*Ldata == -4){ //total ambiguity
				for(int i=0;i<nRateCats;i++){
					La += partial[0];
					Lc += partial[1];
					Lg += partial[2];
					Lt += partial[3];
					partial += 4;
					}
				Ldata++;
				}
			else{ //partial ambiguity
				char nstates=-1 * *(Ldata++);
				for(int i=0;i<nstates;i++){
					for(int i=0;i<nRateCats;i++){
						La += prmat[(*Ldata)+16*i]  * partial[4*i];
						D1a += d1mat[(*Ldata)+16*i] * partial[4*i];		
						D2a += d2mat[(*Ldata)+16*i] * partial[4*i];
											
						Lc += prmat[(*Ldata+4)+16*i] * partial[1+4*i];
						D1c += d1mat[(*Ldata+4)+16*i]* partial[1+4*i];
						D2c += d2mat[(*Ldata+4)+16*i]* partial[1+4*i];
											
						Lg += prmat[(*Ldata+8)+16*i]* partial[2+4*i];
						D1g += d1mat[(*Ldata+8)+16*i]* partial[2+4*i];
						D2g += d2mat[(*Ldata+8)+16*i]* partial[2+4*i];
						
						Lt += prmat[(*Ldata+12)+16*i]* partial[3+4*i];
						D1t += d1mat[(*Ldata+12)+16*i]* partial[3+4*i];
						D2t += d2mat[(*Ldata+12)+16*i]* partial[3+4*i];
						}
					Ldata++;
					}
				partial+=16;
				}
			if((mod->NoPinvInModel() == false) && (i<=lastConst)){
				double btot=0.0;
				if(conBases[i]&1) btot+=mod->Pi(0);
				if(conBases[i]&2) btot+=mod->Pi(1);
				if(conBases[i]&4) btot+=mod->Pi(2);
				if(conBases[i]&8) btot+=mod->Pi(3);
				//6-27-05 fixed this to calc derivs correctly if constant site has been rescaled
				siteL  = ((La*mod->Pi(0)+Lc*mod->Pi(1)+Lg*mod->Pi(2)+Lt*mod->Pi(3)) * scaledGammaProp + (prI*btot)*exp((double)partialCLA->underflow_mult[i]));
				}
			else
				siteL  = ((La*mod->Pi(0)+Lc*mod->Pi(1)+Lg*mod->Pi(2)+Lt*mod->Pi(3)) * scaledGammaProp);

			tempD1 = (((D1a*mod->Pi(0)+D1c*mod->Pi(1)+D1g*mod->Pi(2)+D1t*mod->Pi(3)) * scaledGammaProp) / siteL);
			d1Tot += *countit * tempD1;
			assert(d1Tot == d1Tot);
			double siteD2=((D2a*mod->Pi(0)+D2c*mod->Pi(1)+D2g*mod->Pi(2)+D2t*mod->Pi(3)) * scaledGammaProp);
			d2Tot += *countit * ((siteD2 / siteL) - tempD1*tempD1);
			assert(d2Tot == d2Tot);
			}
		else{
			partial+=16;
			if(!(*Ldata < 0)) Ldata++;
			else if(*Ldata == -4) Ldata++;
			else{
				char nstates=-1 * *(Ldata++);
				for(int i=0;i<nstates;i++) Ldata++;
				}
			}
		countit++;
		}
	}
	
void Tree::GetDerivsPartialInternalFlexRates(const CondLikeArray *partialCLA, const CondLikeArray *childCLA, const double *prmat, const double *d1mat, const double *d2mat, double &d1Tot, double &d2Tot){
	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.
	double *CL1=childCLA->arr;
	double *partial=partialCLA->arr;

	int nchar=data->NChar();

#ifdef UNIX
	madvise(partial, nchar*16*sizeof(double), MADV_SEQUENTIAL);
	madvise(CL1, nchar*16*sizeof(double), MADV_SEQUENTIAL);
#endif

	double siteL, tempD1;
	double La, Lc, Lg, Lt;
	double D1a, D1c, D1g, D1t;
	double D2a, D2c, D2g, D2t;
	
	const int *countit=data->GetCounts();
	const double *rateProb=mod->GetRateProbs();
	const int nRateCats=mod->NRateCats();
	
	int lastConst=data->LastConstant();
	const int *conBases=data->GetConstBases();
	double prI=mod->ProportionInvariant();	

	for(int i=0;i<nchar;i++){
		if((*countit) != 0){
			La=Lc=Lg=Lt=D1a=D1c=D1g=D1t=D2a=D2c=D2g=D2t=0.0;
			for(int r=0;r<nRateCats;r++){
				int rOff=r*16;
				La += ( prmat[rOff ]*CL1[0]+prmat[rOff + 1]*CL1[1]+prmat[rOff + 2]*CL1[2]+prmat[rOff + 3]*CL1[3]) * partial[0] * rateProb[r];
				Lc += ( prmat[rOff + 4]*CL1[0]+prmat[rOff + 5]*CL1[1]+prmat[rOff + 6]*CL1[2]+prmat[rOff + 7]*CL1[3]) * partial[1] * rateProb[r];
				Lg += ( prmat[rOff + 8]*CL1[0]+prmat[rOff + 9]*CL1[1]+prmat[rOff + 10]*CL1[2]+prmat[rOff + 11]*CL1[3]) * partial[2] * rateProb[r];
				Lt += ( prmat[rOff + 12]*CL1[0]+prmat[rOff + 13]*CL1[1]+prmat[rOff + 14]*CL1[2]+prmat[rOff + 15]*CL1[3]) * partial[3] * rateProb[r];
				
				D1a += ( d1mat[rOff ]*CL1[0]+d1mat[rOff + 1]*CL1[1]+d1mat[rOff + 2]*CL1[2]+d1mat[rOff + 3]*CL1[3]) * partial[0] * rateProb[r];
				D1c += ( d1mat[rOff + 4]*CL1[0]+d1mat[rOff + 5]*CL1[1]+d1mat[rOff + 6]*CL1[2]+d1mat[rOff + 7]*CL1[3]) * partial[1] * rateProb[r];
				D1g += ( d1mat[rOff + 8]*CL1[0]+d1mat[rOff + 9]*CL1[1]+d1mat[rOff + 10]*CL1[2]+d1mat[rOff + 11]*CL1[3]) * partial[2] * rateProb[r];
				D1t += ( d1mat[rOff + 12]*CL1[0]+d1mat[rOff + 13]*CL1[1]+d1mat[rOff + 14]*CL1[2]+d1mat[rOff + 15]*CL1[3]) * partial[3] * rateProb[r];		

				D2a += ( d2mat[rOff ]*CL1[0]+d2mat[rOff + 1]*CL1[1]+d2mat[rOff + 2]*CL1[2]+d2mat[rOff + 3]*CL1[3]) * partial[0] * rateProb[r];
				D2c += ( d2mat[rOff + 4]*CL1[0]+d2mat[rOff + 5]*CL1[1]+d2mat[rOff + 6]*CL1[2]+d2mat[rOff + 7]*CL1[3]) * partial[1] * rateProb[r];
				D2g += ( d2mat[rOff + 8]*CL1[0]+d2mat[rOff + 9]*CL1[1]+d2mat[rOff + 10]*CL1[2]+d2mat[rOff + 11]*CL1[3]) * partial[2] * rateProb[r];
				D2t += ( d2mat[rOff + 12]*CL1[0]+d2mat[rOff + 13]*CL1[1]+d2mat[rOff + 14]*CL1[2]+d2mat[rOff + 15]*CL1[3]) * partial[3] * rateProb[r];
				
				partial+=4;
				CL1+=4;
				}
			if((mod->NoPinvInModel() == false) && (i<=lastConst)){
				double btot=0.0;
				if(conBases[i]&1) btot+=mod->Pi(0);
				if(conBases[i]&2) btot+=mod->Pi(1);
				if(conBases[i]&4) btot+=mod->Pi(2);
				if(conBases[i]&8) btot+=mod->Pi(3);
				//6-27-05 fixed this to calc derivs correctly if constant site has been rescaled
				double underTot=childCLA->underflow_mult[i]+partialCLA->underflow_mult[i];
				siteL  = ((La*mod->Pi(0)+Lc*mod->Pi(1)+Lg*mod->Pi(2)+Lt*mod->Pi(3)) + (prI*btot)*exp(underTot));
				}
			else
				siteL  = ((La*mod->Pi(0)+Lc*mod->Pi(1)+Lg*mod->Pi(2)+Lt*mod->Pi(3)));
			tempD1 = (((D1a*mod->Pi(0)+D1c*mod->Pi(1)+D1g*mod->Pi(2)+D1t*mod->Pi(3))) / siteL);
			d1Tot += *countit * tempD1;
			assert(d1Tot == d1Tot);
			double siteD2=((D2a*mod->Pi(0)+D2c*mod->Pi(1)+D2g*mod->Pi(2)+D2t*mod->Pi(3)));
			d2Tot += *countit * ((siteD2 / siteL) - tempD1*tempD1);
			assert(d2Tot == d2Tot);
			}
		else{
			partial+=4*nRateCats;
			CL1+=4*nRateCats;
			}
		countit++;
		}
	}

void Tree::GetDerivsPartialInternalRateHet(const CondLikeArray *partialCLA, const CondLikeArray *childCLA, const double *prmat, const double *d1mat, const double *d2mat, double &d1Tot, double &d2Tot){
	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.
	double *CL1=childCLA->arr;
	double *partial=partialCLA->arr;

	int nchar=data->NChar();

#ifdef UNIX
	madvise(partial, nchar*mod->NRateCats()*4*sizeof(double), MADV_SEQUENTIAL);
	madvise(CL1, nchar*mod->NRateCats()*4*sizeof(double), MADV_SEQUENTIAL);
#endif

	double siteL, tempD1;
	double La, Lc, Lg, Lt;
	double D1a, D1c, D1g, D1t;
	double D2a, D2c, D2g, D2t;
	
	//gamma and invariants
	const int *countit=data->GetCounts();
	int lastConst=data->LastConstant();
	const int *conBases=data->GetConstBases();
	double prI=mod->ProportionInvariant();
	int nRateCats=mod->NRateCats();
	double scaledGammaProp=(1.0-prI) / nRateCats;

	for(int i=0;i<nchar;i++){
		if((*countit) != 0){
			La=Lc=Lg=Lt=D1a=D1c=D1g=D1t=D2a=D2c=D2g=D2t=0.0;
			for(int r=0;r<nRateCats;r++){
				int rOff=r*16;
				La += ( prmat[rOff ]*CL1[0]+prmat[rOff + 1]*CL1[1]+prmat[rOff + 2]*CL1[2]+prmat[rOff + 3]*CL1[3]) * partial[0];
				Lc += ( prmat[rOff + 4]*CL1[0]+prmat[rOff + 5]*CL1[1]+prmat[rOff + 6]*CL1[2]+prmat[rOff + 7]*CL1[3]) * partial[1];
				Lg += ( prmat[rOff + 8]*CL1[0]+prmat[rOff + 9]*CL1[1]+prmat[rOff + 10]*CL1[2]+prmat[rOff + 11]*CL1[3]) * partial[2];
				Lt += ( prmat[rOff + 12]*CL1[0]+prmat[rOff + 13]*CL1[1]+prmat[rOff + 14]*CL1[2]+prmat[rOff + 15]*CL1[3]) * partial[3];
				
				D1a += ( d1mat[rOff ]*CL1[0]+d1mat[rOff + 1]*CL1[1]+d1mat[rOff + 2]*CL1[2]+d1mat[rOff + 3]*CL1[3]) * partial[0];
				D1c += ( d1mat[rOff + 4]*CL1[0]+d1mat[rOff + 5]*CL1[1]+d1mat[rOff + 6]*CL1[2]+d1mat[rOff + 7]*CL1[3]) * partial[1];
				D1g += ( d1mat[rOff + 8]*CL1[0]+d1mat[rOff + 9]*CL1[1]+d1mat[rOff + 10]*CL1[2]+d1mat[rOff + 11]*CL1[3]) * partial[2];
				D1t += ( d1mat[rOff + 12]*CL1[0]+d1mat[rOff + 13]*CL1[1]+d1mat[rOff + 14]*CL1[2]+d1mat[rOff + 15]*CL1[3]) * partial[3];		

				D2a += ( d2mat[rOff ]*CL1[0]+d2mat[rOff + 1]*CL1[1]+d2mat[rOff + 2]*CL1[2]+d2mat[rOff + 3]*CL1[3]) * partial[0];
				D2c += ( d2mat[rOff + 4]*CL1[0]+d2mat[rOff + 5]*CL1[1]+d2mat[rOff + 6]*CL1[2]+d2mat[rOff + 7]*CL1[3]) * partial[1];
				D2g += ( d2mat[rOff + 8]*CL1[0]+d2mat[rOff + 9]*CL1[1]+d2mat[rOff + 10]*CL1[2]+d2mat[rOff + 11]*CL1[3]) * partial[2];
				D2t += ( d2mat[rOff + 12]*CL1[0]+d2mat[rOff + 13]*CL1[1]+d2mat[rOff + 14]*CL1[2]+d2mat[rOff + 15]*CL1[3]) * partial[3];
				
				partial+=4;
				CL1+=4;
				}
			if((mod->NoPinvInModel() == false) && (i<=lastConst)){
				double btot=0.0;
				if(conBases[i]&1) btot+=mod->Pi(0);
				if(conBases[i]&2) btot+=mod->Pi(1);
				if(conBases[i]&4) btot+=mod->Pi(2);
				if(conBases[i]&8) btot+=mod->Pi(3);
				//6-27-05 fixed this to calc derivs correctly if constant site has been rescaled
				double underTot=childCLA->underflow_mult[i]+partialCLA->underflow_mult[i];
				siteL  = ((La*mod->Pi(0)+Lc*mod->Pi(1)+Lg*mod->Pi(2)+Lt*mod->Pi(3)) * scaledGammaProp + (prI*btot)*exp(underTot));
				}
			else
				siteL  = ((La*mod->Pi(0)+Lc*mod->Pi(1)+Lg*mod->Pi(2)+Lt*mod->Pi(3)) * scaledGammaProp);
			tempD1 = (((D1a*mod->Pi(0)+D1c*mod->Pi(1)+D1g*mod->Pi(2)+D1t*mod->Pi(3)) * scaledGammaProp) / siteL);
			d1Tot += *countit * tempD1;
			assert(d1Tot == d1Tot);
			double siteD2=((D2a*mod->Pi(0)+D2c*mod->Pi(1)+D2g*mod->Pi(2)+D2t*mod->Pi(3)) * scaledGammaProp);
			d2Tot += *countit * ((siteD2 / siteL) - tempD1*tempD1);
			assert(d2Tot == d2Tot);
			}
		else{
			partial+=16;
			CL1+=16;
			}
		countit++;
		}
	}

/*DEPRECATED
void Tree::FillSiteLikes(const CondLikeArray *fullCla, double *dest, bool addPinv){
	
	double *cla=fullCla->arr;
	
	int nSites=data->NChar();
	const int *countit=data->GetCounts();
	double Lk;
	//gamma and invariants
	int lastConst=data->LastConstant();
	const int *conBases=data->GetConstBases();
	double prI=mod->ProportionInvariant();
	double scaledGammaProp=(1.0-prI) / mod->NRateCats();
	
	//if we're doing these calcs for the first or second derivatives, the contribution
	//of pinv should not be included.  This is indicated by addPinv being false, so we 
	//just jump to the second loop for all sites
	int lastPinvIndex = (addPinv ? lastConst : -1);
		
	for( int k = 0; k <= lastPinvIndex; k++ ){
		Lk = mod->Pi(0)*(cla[0] + cla[4] + cla[8] + cla[12]) + mod->Pi(1)*(cla[1] + cla[5] + cla[9] + cla[13]) + mod->Pi(2)*(cla[2] + cla[6] + cla[10] + cla[14]) + mod->Pi(3)*(cla[3] + cla[7] + cla[11] + cla[15]);	
		cla+=16;
		double btot=0.0;
		if(conBases[k]&1) btot+=mod->Pi(0);
		if(conBases[k]&2) btot+=mod->Pi(1);
		if(conBases[k]&4) btot+=mod->Pi(2);
		if(conBases[k]&8) btot+=mod->Pi(3);
		dest[k] = Lk * scaledGammaProp + (prI * btot);
		}			
	
	for( int k = lastPinvIndex+1; k < nSites; k++ ){
		Lk = mod->Pi(0)*(cla[0] + cla[4] + cla[8] + cla[12]) + mod->Pi(1)*(cla[1] + cla[5] + cla[9] + cla[13]) + mod->Pi(2)*(cla[2] + cla[6] + cla[10] + cla[14]) + mod->Pi(3)*(cla[3] + cla[7] + cla[11] + cla[15]);	
		cla+=16;
		dest[k] = Lk * scaledGammaProp;
		}
	}
*/	

	
	
	
