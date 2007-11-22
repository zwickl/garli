// GARLI version 0.951 source code
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

#include "defs.h"
#include "tree.h"
#include "model.h"
#include "funcs.h"
#include "outputman.h"

//a bunch of functions from the Tree class, relating to optimization

#include "utility.h"
Profiler ProfIntDeriv ("IntDeriv      ");
Profiler ProfTermDeriv("TermDeriv     ");
Profiler ProfModDeriv ("ModDeriv      ");
Profiler ProfNewton   ("Newton-Raphson");
extern Profiler ProfEQVectors;
				   

extern FLOAT_TYPE globalBest;

extern int optCalcs;

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
#define effectiveMin=min_brlen
#define effectiveMax=max_brlen
#endif

inline FLOAT_TYPE CallBranchLike(TreeNode *thisnode, Tree *thistree, FLOAT_TYPE blen, bool brak /*=false*/){
	brak;

#ifdef FOURTH_ROOT
	thisnode->dlen=blen*blen*blen*blen;
#elif ROOT_OPT
	thisnode->dlen=blen*blen;
#else
	thisnode->dlen=blen;
#endif
	FLOAT_TYPE like=thistree->BranchLike(thisnode)*-1;
	
	optCalcs++;

#ifdef OPT_DEBUG
	if(brak) optInfo.BrakAdd(blen, like);
	else optInfo.BrentAdd(blen, like);

#endif

	return like;
	}

void Tree::OptimizeBranchesInArray(int *nodes, int numNodes, FLOAT_TYPE optPrecision){
	//this takes an array of nodeNums (branches) to be optimized and does so
	for(int i=0;i<numNodes;i++){
		BrentOptimizeBranchLength(optPrecision, allNodes[nodes[i]], true);
		}	
	}

FLOAT_TYPE Tree::OptimizeAllBranches(FLOAT_TYPE optPrecision){
	FLOAT_TYPE improve=ZERO_POINT_ZERO;
	SetNodesUnoptimized();
	improve = RecursivelyOptimizeBranches(root->left, optPrecision, 0, numNodesTotal, true, improve, true);
	improve = RecursivelyOptimizeBranches(root->left->next, optPrecision, 0, numNodesTotal, true, improve, true);
	improve = RecursivelyOptimizeBranches(root->right, optPrecision, 0, numNodesTotal, true, improve, true);

	return improve;
	}

FLOAT_TYPE Tree::OptimizeTreeScale(FLOAT_TYPE optPrecision){
	if(lnL==-ONE_POINT_ZERO) Score();
	Score();
	FLOAT_TYPE start=lnL;
	FLOAT_TYPE prev=lnL;
	FLOAT_TYPE cur;
	FLOAT_TYPE scale;
	FLOAT_TYPE t;
	FLOAT_TYPE lastChange=(FLOAT_TYPE)9999.9;
	FLOAT_TYPE effectiveScale = ONE_POINT_ZERO; //this measures the change in scale relative to what it began at.
	FLOAT_TYPE upperBracket = FLT_MAX;   //the smallest value we know of with a negative d1 (relative to inital scale of 1.0!)
	FLOAT_TYPE lowerBracket = FLT_MIN;   //the largest value we know of with a positive d1 (relative to inital scale of 1.0!)
	FLOAT_TYPE incr;

#ifdef DEBUG_SCALE_OPT
	ofstream deb("scaleTrace.log");
	deb.precision(20);
	for(int s=0;s<50;s++){
		FLOAT_TYPE scale=0.5 + s*.025;
		ScaleWholeTree(scale);
		Score();
		deb << scale << "\t" << lnL << endl;
		ScaleWholeTree(ONE_POINT_ZERO/scale);	
		}
	deb.close();
#endif

	while(1){
		//reversed this now so the reduction in scale is done first when getting the 
		//derivs.  This works better if some blens are at DEF_MAX_BLEN because the 
		//scaling up causes them to hit the max and the relative blens to change

#ifdef SINGLE_PRECISION_FLOATS
		incr=0.005f;
#else
		incr=0.0001;
#endif

		scale=ONE_POINT_ZERO-incr;

		ScaleWholeTree(scale);
		Score();
		cur=lnL;
		ScaleWholeTree(ONE_POINT_ZERO/scale);//return the tree to its original scale	
		FLOAT_TYPE d12=(cur-prev)/-incr;

		scale=ONE_POINT_ZERO + incr;
		ScaleWholeTree(scale);
		Score();
		cur=lnL;
		ScaleWholeTree(ONE_POINT_ZERO/scale);//return the tree to its original scale
		FLOAT_TYPE d11=(cur-prev)/incr;

		FLOAT_TYPE d1=(d11+d12)*ZERO_POINT_FIVE;
		FLOAT_TYPE d2=(d11-d12)/incr;
		
		FLOAT_TYPE est = -d1/d2;
		FLOAT_TYPE estImprove = d1*est + d2*(est*est*ZERO_POINT_FIVE);

		//return conditions
		if(estImprove < optPrecision && d2 < ZERO_POINT_ZERO){
			return prev-start;
			}
		
		if(d2 < ZERO_POINT_ZERO){
			est = max(min((FLOAT_TYPE)0.1, est), (FLOAT_TYPE)-0.1);
			t=ONE_POINT_ZERO + est;
			}
		else{//if we have lots of data, move
			//very slowly here
			//DEBUG
			//if(data->NInformative() > 500){
			if(0){
				if(d1 > ZERO_POINT_ZERO) t=(FLOAT_TYPE)1.01;
				else t=(FLOAT_TYPE)0.99;
				}
			else{
				if(d1 > ZERO_POINT_ZERO) t=(FLOAT_TYPE)1.05;
				else t=(FLOAT_TYPE)0.95;
				}
			}
		
		//update the brackets
		if(d1 <= ZERO_POINT_ZERO && effectiveScale < upperBracket)
			upperBracket = effectiveScale;
		else if(d1 > ZERO_POINT_ZERO && effectiveScale > lowerBracket)
			lowerBracket = effectiveScale;

		//if the surface is wacky and we are going to shoot past one of our brackets
		//take evasive action by going halfway to the bracket
		if((effectiveScale * t) <= lowerBracket){
			t = (lowerBracket + effectiveScale) * ZERO_POINT_FIVE / effectiveScale;
			}
		else if((effectiveScale * t) >= upperBracket){
			t = (upperBracket + effectiveScale) * ZERO_POINT_FIVE / effectiveScale;
			}

		scale=t;
		effectiveScale *= scale;
		ScaleWholeTree(scale);
		Score();
		cur=lnL;
		lastChange = cur - prev;
		prev=cur;
		}
	return -1;
	}

FLOAT_TYPE Tree::OptimizeAlpha(FLOAT_TYPE optPrecision){

#ifdef DEBUG_ALPHA_OPT
	FLOAT_TYPE initVal=mod->Alpha();

	ofstream deb("alphaTrace.log");
	deb.precision(20);
	for(int s=0;s<50;s++){
		FLOAT_TYPE a=.3 + s*.025;
		mod->SetAlpha(a, false);
		MakeAllNodesDirty();
		Score();
		deb << a << "\t" << lnL << endl;
		}
	deb.close();
	mod->SetAlpha(initVal, false);
	MakeAllNodesDirty();
	Score();	
#endif

	if(lnL==-1) Score();
	FLOAT_TYPE start, prev, cur;
	prev = start = cur = lnL;
	FLOAT_TYPE prevVal=mod->Alpha();
	FLOAT_TYPE lastChange=(FLOAT_TYPE)9999.9;
	FLOAT_TYPE upperBracket = FLT_MAX;   //the smallest value we know of with a negative d1 
	FLOAT_TYPE lowerBracket = FLT_MIN;   //the largest value we know of with a positive d1 
	FLOAT_TYPE incr;

	while(1){
#ifdef SINGLE_PRECISION_FLOATS
		incr=0.005f;
#else
		incr=0.001;
#endif
		mod->SetAlpha(prevVal+incr, false);
		
		MakeAllNodesDirty();
		Score();
		cur=lnL;
		FLOAT_TYPE d11=(cur-prev)/incr;
		
		mod->SetAlpha(prevVal-incr, false);
		MakeAllNodesDirty();
		Score();
		cur=lnL;
		FLOAT_TYPE d12=(cur-prev)/-incr;
		
		FLOAT_TYPE d1=(d11+d12)*ZERO_POINT_FIVE;
		FLOAT_TYPE d2=(d11-d12)/incr;
		
		FLOAT_TYPE est=-d1/d2;
		FLOAT_TYPE estImprove = d1*est + d2*(est*est*ZERO_POINT_FIVE);
		
		if(estImprove < optPrecision && d2 <= ZERO_POINT_ZERO){
			mod->SetAlpha(prevVal, false);
			MakeAllNodesDirty();			
			return prev-start;
			}
		
		FLOAT_TYPE t;
		if(d2 < ZERO_POINT_ZERO){
			t=prevVal+est;
			if((t < 0.05) && (d1 < 0) && (prevVal*ZERO_POINT_FIVE > t)){
				t=prevVal*ZERO_POINT_FIVE;
				}
			}
		else{
			if(d1 > ZERO_POINT_ZERO) t=prevVal*(FLOAT_TYPE)1.25;
			else t=prevVal*(FLOAT_TYPE)0.8;
			}
		
		//update the brackets
		if(d1 <= ZERO_POINT_ZERO && prevVal < upperBracket)
			upperBracket = prevVal;
		else if(d1 > ZERO_POINT_ZERO && prevVal > lowerBracket)
			lowerBracket = prevVal;

		//if the surface is wacky and we are going to shoot past one of our brackets
		//take evasive action by going halfway to the bracket
		if(t <= lowerBracket){
			t = (lowerBracket + prevVal) * ZERO_POINT_FIVE;
			}
		else if(t >= upperBracket){
			t = (upperBracket + prevVal) * ZERO_POINT_FIVE;
			}

		mod->SetAlpha(t, false);
		assert((prevVal==0.05 && mod->Alpha()==0.05)==false);
		MakeAllNodesDirty();			
		Score();
		lastChange = lnL - prev;
		prev=lnL;
		prevVal=t;
		}
	return -1;
	}

FLOAT_TYPE Tree::OptimizeBoundedParameter(FLOAT_TYPE optPrecision, FLOAT_TYPE prevVal, int which, FLOAT_TYPE lowBound, FLOAT_TYPE highBound, void (Model::*SetParam)(int, FLOAT_TYPE)){

#ifdef DEBUG_ALPHA_OPT
	FLOAT_TYPE initVal=mod->Alpha();

	ofstream deb("alphaTrace.log");
	deb.precision(20);
	for(int s=0;s<50;s++){
		FLOAT_TYPE a=.3 + s*.025;
		mod->SetAlpha(a, false);
		MakeAllNodesDirty();
		Score();
		deb << a << "\t" << lnL << endl;
		}
	deb.close();
	mod->SetAlpha(initVal, false);
	MakeAllNodesDirty();
	Score();	
#endif

	if(lnL==-1) Score();
	FLOAT_TYPE start, prev, cur;
	prev = start = cur = lnL;
	FLOAT_TYPE lastChange=(FLOAT_TYPE)9999.9;
	FLOAT_TYPE upperBracket = highBound;   //the smallest value we know of with a negative d1, or the minimum allowed value
	FLOAT_TYPE lowerBracket = lowBound;   //the largest value we know of with a positive d1 , or the maximum allowed value
	FLOAT_TYPE incr;
	FLOAT_TYPE epsilon = 1e-5;
	int lowBoundOvershoot = 0;
	int upperBoundOvershoot = 0;

	while(1){
#ifdef SINGLE_PRECISION_FLOATS
		incr=0.005f;
#else
		incr=min(0.001, min(prevVal - lowerBracket, upperBracket - prevVal));
#endif
		//mod->SetAlpha(prevVal+incr, false);
		CALL_SET_PARAM_FUNCTION(*mod, SetParam)(which, prevVal+incr);

		MakeAllNodesDirty();
		Score();
		cur=lnL;
		FLOAT_TYPE d11=(cur-prev)/incr;
		
		//mod->SetAlpha(prevVal-incr, false);
		CALL_SET_PARAM_FUNCTION(*mod, SetParam)(which, prevVal-incr);
		MakeAllNodesDirty();
		Score();
		cur=lnL;
		FLOAT_TYPE d12=(cur-prev)/-incr;
		
		FLOAT_TYPE d1=(d11+d12)*ZERO_POINT_FIVE;
		if((d11 - d12) == ZERO_POINT_ZERO){
			CALL_SET_PARAM_FUNCTION(*mod, SetParam)(which, prevVal);;
			MakeAllNodesDirty();			
			return prev-start;
			}
		
		FLOAT_TYPE d2=(d11-d12)/incr;
		
		FLOAT_TYPE est=-d1/d2;
		
		FLOAT_TYPE proposed = prevVal + est;

		if(d2 > ZERO_POINT_ZERO){
			if(d1 > ZERO_POINT_ZERO) proposed=prevVal*(FLOAT_TYPE)1.25;
			else proposed=prevVal*(FLOAT_TYPE)0.8;
			}
		if(proposed < (lowerBracket + epsilon)){
			if(prevVal - lowerBracket - epsilon < epsilon * ZERO_POINT_FIVE){
				CALL_SET_PARAM_FUNCTION(*mod, SetParam)(which, prevVal);;
				MakeAllNodesDirty();			
				return prev-start;
				}
			lowBoundOvershoot++;
			if(lowBoundOvershoot > 1)
				proposed = lowerBracket + epsilon;
			else
				proposed = (prevVal + lowerBracket) * ZERO_POINT_FIVE;
			}
		else if(proposed > upperBracket - epsilon){
			if(upperBracket - epsilon - prevVal < epsilon * ZERO_POINT_FIVE){
				CALL_SET_PARAM_FUNCTION(*mod, SetParam)(which, prevVal);;
				MakeAllNodesDirty();			
				return prev-start;
				}
			upperBoundOvershoot++;
			if(upperBoundOvershoot > 1)
				proposed = upperBracket - epsilon;
			else
				proposed = (prevVal + upperBracket) * ZERO_POINT_FIVE;
			}

		FLOAT_TYPE estImprove;
		if(d2 < ZERO_POINT_ZERO) estImprove = d1*(proposed - prevVal) + (d2 * (proposed - prevVal) * (proposed - prevVal)) * ZERO_POINT_FIVE;
		else estImprove = 9999.9;

		if(estImprove < optPrecision){
			//mod->SetAlpha(prevVal, false);
			CALL_SET_PARAM_FUNCTION(*mod, SetParam)(which, prevVal);;
			MakeAllNodesDirty();			
			return prev-start;
			}
		
		//update the brackets
		if(d1 <= ZERO_POINT_ZERO && prevVal < upperBracket)
			upperBracket = prevVal;
		else if(d1 > ZERO_POINT_ZERO && prevVal > lowerBracket)
			lowerBracket = prevVal;

		//mod->SetAlpha(t, false);
		CALL_SET_PARAM_FUNCTION(*mod, SetParam)(which, proposed);;
		MakeAllNodesDirty();			
		Score();
		lastChange = lnL - prev;
		prev=lnL;
		prevVal=proposed;
		}
	return -1;
	}


FLOAT_TYPE Tree::OptimizePinv(){


/*	ofstream deb("debug.log");
	deb.precision(20);
	for(int s=0;s<30;s++){
		FLOAT_TYPE scale=.15 + s*.01;
		mod->SetPinv(scale);
		MakeAllNodesDirty();
		Score();
		deb << scale << "\t" << lnL << endl;
		}
	deb.close();
*/
	if(lnL==-1) Score();
	FLOAT_TYPE start=lnL;
	FLOAT_TYPE prev=lnL;
	FLOAT_TYPE cur;
	FLOAT_TYPE prevVal=mod->PropInvar();
	
	while(1){
		FLOAT_TYPE incr=(FLOAT_TYPE)0.001;
		mod->SetPinv(prevVal+incr, false);
		
		MakeAllNodesDirty();
		Score();
		cur=lnL;
		FLOAT_TYPE d11=(cur-prev)/incr;
		
		mod->SetPinv(prevVal-incr, false);
		MakeAllNodesDirty();
		Score();
		cur=lnL;
		FLOAT_TYPE d12=(cur-prev)/-incr;
		
		FLOAT_TYPE d1=(d11+d12)*ZERO_POINT_FIVE;
		FLOAT_TYPE d2=(d11-d12)/incr;
		
		FLOAT_TYPE est=-d1/d2;
		
		if((abs(est) < 0.001 && d2 < ZERO_POINT_ZERO) || (d1>ZERO_POINT_ZERO && prevVal==mod->MaxPinv())){
			mod->SetPinv(prevVal, false);
			MakeAllNodesDirty();			
			return prev-start;
			}
		
		FLOAT_TYPE t;
		if(d2 < ZERO_POINT_ZERO){
			t=prevVal+est;
			if(t < ZERO_POINT_ZERO) t=prevVal*ZERO_POINT_FIVE;
			//if(t > mod->MaxPinv()) t=(prevVal+mod->MaxPinv())*.5;
			if(t > mod->MaxPinv()) t=mod->MaxPinv();
			}
		else{
			if(d1 > ZERO_POINT_ZERO) t=prevVal*(FLOAT_TYPE)1.1;
			else t=prevVal*(FLOAT_TYPE)0.9;
	//		if(t > mod->MaxPinv()) t=(prevVal+mod->MaxPinv())*ZERO_POINT_FIVE;
			if(t > mod->MaxPinv()) t=mod->MaxPinv();
			}
			
		mod->SetPinv(t, false);
		MakeAllNodesDirty();			
		Score();
		prev=lnL;
		prevVal=t;	
		}
	return -1;
	}

int num=1;

FLOAT_TYPE Tree::OptimizeBranchLength(FLOAT_TYPE optPrecision, TreeNode *nd, bool goodGuess){
	nd->alreadyOptimized=true;
	FLOAT_TYPE improve;
	
#ifdef OPT_DEBUG
	optsum << nd->nodeNum << "\t" << nd->dlen << "\t";
#endif
		
#ifdef BRENT
	improve = BrentOptimizeBranchLength(optPrecision, nd, goodGuess);
#else
	//improve = NewtonRaphsonOptimizeBranchLength(optPrecision, nd, goodGuess);
	//abandoning use of goodGuess.  Doesn't seem to be reducing opt passes, which
	//was the point.
	//optCalcs++;
	ProfNewton.Start();
	improve = NewtonRaphsonOptimizeBranchLength(optPrecision, nd, true);
	ProfNewton.Stop();
#endif

#ifdef OPT_DEBUG
	optsum << nd->dlen << "\t" << improve << endl;
	
/*	ofstream opttrees;
	if(num == 1) opttrees.open("everyTree.tre");
	else opttrees.open("everyTree.tre", ios::app);
	char treeString[20000];
	root->MakeNewick(treeString, false, true);
	opttrees <<  "utree tree" << num++ << "_" << nd->nodeNum << "=" << treeString << ";" << endl;
	opttrees.close();
*/
#endif

	return improve;
	}

void Tree::SetNodesUnoptimized(){
	root->left->SetUnoptimized();
	root->left->next->SetUnoptimized();
	root->right->SetUnoptimized();
	}

void Tree::OptimizeBranchesWithinRadius(TreeNode *nd, FLOAT_TYPE optPrecision, int subtreeNode, TreeNode *prune){
	nodeOptVector.clear();
	SetNodesUnoptimized();

#ifdef EQUIV_CALCS
	if(dirtyEQ){
		ProfEQVectors.Start();
		root->SetEquivalentConditionalVectors(data);
		ProfEQVectors.Stop();
		dirtyEQ=false;
		}
#endif

	FLOAT_TYPE totalIncrease=ZERO_POINT_ZERO, prunePointIncrease=ZERO_POINT_ZERO, thisIncr, pruneRadIncrease=ZERO_POINT_ZERO;
	
#ifdef CHECK_LNL_BEFORE_RAD
	FLOAT_TYPE leftIncrease=ZERO_POINT_ZERO, rightIncrease=ZERO_POINT_ZERO, ancIncrease=ZERO_POINT_ZERO;
	leftIncrease = OptimizeBranchLength(optPrecision, nd->left, false);
	ancIncrease = OptimizeBranchLength(optPrecision, nd, false);
	rightIncrease = OptimizeBranchLength(optPrecision, nd->right, false);	

	if(leftIncrease > ZERO_POINT_ZERO) nodeOptVector.push_back(nd->left);
	if(ancIncrease > ZERO_POINT_ZERO) nodeOptVector.push_back(nd);
	if(rightIncrease > ZERO_POINT_ZERO) nodeOptVector.push_back(nd->right);
	totalIncrease = leftIncrease + rightIncrease + ancIncrease;
#else
	totalIncrease += OptimizeBranchLength(optPrecision, nd->left, false);
	totalIncrease+= OptimizeBranchLength(optPrecision, nd, false);
	totalIncrease += OptimizeBranchLength(optPrecision, nd->right, false);	

	nodeOptVector.push_back(nd->left);
	nodeOptVector.push_back(nd);
	nodeOptVector.push_back(nd->right);

#endif

	if(prune!=NULL){
		prunePointIncrease = OptimizeBranchLength(optPrecision, prune, true);
		if(prunePointIncrease > ZERO_POINT_ZERO) nodeOptVector.push_back(prune);
		totalIncrease+=prunePointIncrease;
		Score(prune->anc->nodeNum);

		#ifdef OPT_DEBUG
		optsum << "prune total\t" << prunePointIncrease << endl;
		#endif
		}
	else Score(nd->nodeNum);

	#ifdef OPT_DEBUG
	optsum << "3/4 branch total\t" << totalIncrease << endl;
	if(lnL < globalBest - treeRejectionThreshold)
		optsum << "bailing early\t";
	optsum << "Scores:" << globalBest << "\t" << lnL << "\t" << globalBest-lnL << endl;
	#endif

#ifdef CHECK_LNL_BEFORE_RAD
	bool fullOpt = false;
	if(lnL > globalBest){
		fullOpt = true;
		}
#endif

	if(lnL < globalBest - treeRejectionThreshold){
		return;
		}

	//now spread out
	int rad=10;

	if(rad>0){
#ifdef CHECK_LNL_BEFORE_RAD
		if(((rightIncrease > ZERO_POINT_ZERO) || fullOpt) && nd->right->left!=NULL && nd->right->left->alreadyOptimized==false) totalIncrease += RecursivelyOptimizeBranches(nd->right->left, optPrecision, subtreeNode, rad, false, ZERO_POINT_ZERO);
		if(((leftIncrease > ZERO_POINT_ZERO) || fullOpt) && nd->left->left!=NULL && nd->left->left->alreadyOptimized==false) totalIncrease += RecursivelyOptimizeBranches(nd->left->left, optPrecision, subtreeNode, rad, false, ZERO_POINT_ZERO);
		if(((ancIncrease > ZERO_POINT_ZERO)) || fullOpt){
#else
		if(nd->right->left!=NULL && nd->right->left->alreadyOptimized==false) totalIncrease += RecursivelyOptimizeBranches(nd->right->left, optPrecision, subtreeNode, rad, false, ZERO_POINT_ZERO);
		if(nd->left->left!=NULL && nd->left->left->alreadyOptimized==false) totalIncrease += RecursivelyOptimizeBranches(nd->left->left, optPrecision, subtreeNode, rad, false, ZERO_POINT_ZERO);
		if(1){
#endif
			if(nd->anc!=root && nd->anc->alreadyOptimized==false) totalIncrease += RecursivelyOptimizeBranchesDown(nd->anc, nd, optPrecision, subtreeNode, rad, ZERO_POINT_ZERO);
			else{
				if(root->left != nd && root->left->alreadyOptimized==false) totalIncrease +=  RecursivelyOptimizeBranches(root->left, optPrecision, subtreeNode, rad, true, ZERO_POINT_ZERO);
				if(root->left->next != nd && root->left->next->alreadyOptimized==false) totalIncrease +=  RecursivelyOptimizeBranches(root->left->next, optPrecision, subtreeNode, rad, true, ZERO_POINT_ZERO);
				if(root->right != nd && root->right->alreadyOptimized==false) totalIncrease +=  RecursivelyOptimizeBranches(root->right, optPrecision, subtreeNode, rad, true, ZERO_POINT_ZERO);
				}
			}
		
		if(prunePointIncrease > ZERO_POINT_ZERO){//now doing a radius opt at the prune point starting from the 4 branches attached to that branch
			//in other words, this is no longer centered on a node, but on a branch
			if(prune->left != NULL){
				if(prune->right->alreadyOptimized==false) pruneRadIncrease += RecursivelyOptimizeBranches(prune->right, optPrecision, subtreeNode, rad, true, ZERO_POINT_ZERO);
				if(prune->left->alreadyOptimized==false) pruneRadIncrease += RecursivelyOptimizeBranches(prune->left, optPrecision, subtreeNode, rad, true, ZERO_POINT_ZERO);
				}
			if(prune == root){
				assert(0);
				//pruneRadIncrease += RecursivelyOptimizeBranches(prune->left->next, optPrecision, subtreeNode, rad, false, ZERO_POINT_ZERO);
				}
			else{//this RecursivelyOptimizeBranchesDown will implicitly also optimize prune's next or prev
				if(prune->anc!=root){
					pruneRadIncrease += RecursivelyOptimizeBranchesDown(prune->anc, prune, optPrecision, subtreeNode, rad, ZERO_POINT_ZERO);
					}
				else{
					if(root->left != prune && root->left->alreadyOptimized==false) totalIncrease +=  RecursivelyOptimizeBranches(root->left, optPrecision, subtreeNode, rad, true, ZERO_POINT_ZERO);
					if(root->left->next != prune && root->left->next->alreadyOptimized==false) totalIncrease +=  RecursivelyOptimizeBranches(root->left->next, optPrecision, subtreeNode, rad, true, ZERO_POINT_ZERO);
					if(root->right != prune && root->right->alreadyOptimized==false) totalIncrease +=  RecursivelyOptimizeBranches(root->right, optPrecision, subtreeNode, rad, true, ZERO_POINT_ZERO);
					}
				}
//				if(prune->next != NULL) pruneRadIncrease += RecursivelyOptimizeBranches(prune->next, optPrecision, subtreeNode, rad, true, ZERO_POINT_ZERO);
//				if(prune->prev != NULL) pruneRadIncrease += RecursivelyOptimizeBranches(prune->prev, optPrecision, subtreeNode, rad, true, ZERO_POINT_ZERO);
//				}
			/*
			if(prune->right->left!=NULL) totalIncrease += RecursivelyOptimizeBranches(prune->right->left, optPrecision, subtreeNode, rad, false, ZERO_POINT_ZERO);
			if(prune->left->left!=NULL) totalIncrease += RecursivelyOptimizeBranches(prune->left->left, optPrecision, subtreeNode, rad, false, ZERO_POINT_ZERO);
			if(prune==root){
				if(prune->left->next->left!=NULL) totalIncrease += RecursivelyOptimizeBranches(prune->left->next->left, optPrecision, subtreeNode, rad, false, ZERO_POINT_ZERO);
				}
			else{
				if(prune->anc!=root) totalIncrease += RecursivelyOptimizeBranchesDown(prune->anc, prune, optPrecision, subtreeNode, rad, ZERO_POINT_ZERO);
				else{
					if(prune->nodeNum!=subtreeNode){
						if(root->left != prune) totalIncrease +=  RecursivelyOptimizeBranches(root->left, optPrecision, subtreeNode, rad, true, ZERO_POINT_ZERO);
						if(root->left->next != prune) totalIncrease +=  RecursivelyOptimizeBranches(root->left->next, optPrecision, subtreeNode, rad, true, ZERO_POINT_ZERO);
						if(root->right != prune) totalIncrease +=  RecursivelyOptimizeBranches(root->right, optPrecision, subtreeNode, rad, true, ZERO_POINT_ZERO);
						}
					}
				}
			*/
			totalIncrease += pruneRadIncrease;
			#ifdef OPT_DEBUG
			optsum << "pruneRadOpt total\t" << pruneRadIncrease << endl;
			#endif
			}
		}
#ifdef OPT_DEBUG
optsum << "radiusopt intial pass total\t" << totalIncrease << endl;
#endif

	FLOAT_TYPE postIncr=ZERO_POINT_ZERO;
	TreeNode *finalNode;

	//no longer need to do this because checking if nodes have already been optimized before doing it again
/*		if(nodeOptVector.empty()==false && prune != NULL){//remove any duplicate entries caused by an overlaping radius around connector and prune
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
*/	
	if(nodeOptVector.empty()) finalNode = nd->right;
	else{
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
		}
	if(finalNode != root)
		Score(finalNode->anc->nodeNum);
	else Score(0);
	totalIncrease += postIncr;

//	if(fourBranchTot > treeRejectionThreshold ) cout << "r\t" << (totalIncrease - fourBranchTot) << endl;
//	else if(treeRejectionThreshold < (totalIncrease - fourBranchTot)) cout << (totalIncrease - fourBranchTot) << endl;

#ifdef OPT_DEBUG
optsum << "postopt total\t" << postIncr << endl;
optsum << "total\t" << totalIncrease << endl;
#endif
	}

FLOAT_TYPE Tree::RecursivelyOptimizeBranches(TreeNode *nd, FLOAT_TYPE optPrecision, int subtreeNode, int radius, bool dontGoNext, FLOAT_TYPE scoreIncrease, bool ignoreDelta/*=false*/){
	FLOAT_TYPE delta = ZERO_POINT_ZERO;
	
	if(nd->alreadyOptimized == false) delta = OptimizeBranchLength(optPrecision, nd, true);
	scoreIncrease += delta;
	
	if(!(delta < optPrecision))
		nodeOptVector.push_back(nd);
//	if(radius==0) cout << "hit max radius!" <<endl;
	if(nd->left!=NULL && radius>1 && (!(delta < optPrecision) || ignoreDelta == true)){
		/*if(nd->left->alreadyOptimized == false)*/ scoreIncrease += RecursivelyOptimizeBranches(nd->left, optPrecision, subtreeNode, radius-1, false, 0, ignoreDelta);
		}
	if(nd->next!=NULL && dontGoNext==false){
		if(nd->next->alreadyOptimized == false) scoreIncrease += RecursivelyOptimizeBranches(nd->next, optPrecision, subtreeNode, radius, false, 0, ignoreDelta);
		}
	if(memLevel > 1) RemoveTempClaReservations();
	return scoreIncrease;
	}

FLOAT_TYPE Tree::RecursivelyOptimizeBranchesDown(TreeNode *nd, TreeNode *calledFrom, FLOAT_TYPE optPrecision, int subtreeNode, int radius, FLOAT_TYPE scoreIncrease){
		
	FLOAT_TYPE delta = ZERO_POINT_ZERO;
	if(nd->alreadyOptimized==false)//because the next or prev of calledFrom could be unoptimized
		//even if nd has been, this check needs to be done here, rather than before calling this func
		delta = OptimizeBranchLength(optPrecision, nd, true);
	scoreIncrease += delta;

	if(nd->nodeNum == subtreeNode){
		return scoreIncrease;
		}

	if(!(delta < optPrecision))
		nodeOptVector.push_back(nd);
//	if(radius==0) cout << "hit max radius!" <<endl;
	if(nd->left!=NULL && nd->left!=calledFrom && radius>1){
		if(nd->left->alreadyOptimized == false) scoreIncrease += RecursivelyOptimizeBranches(nd->left, optPrecision, subtreeNode, radius, true, 0);
		}
	else if(radius>1 && nd->left->next->alreadyOptimized == false) scoreIncrease += RecursivelyOptimizeBranches(nd->left->next, optPrecision, subtreeNode, radius, false, 0);
	if(nd->anc!=root && radius>1 && !(delta < optPrecision)){
		if(nd->anc->alreadyOptimized == false) scoreIncrease += RecursivelyOptimizeBranchesDown(nd->anc, nd, optPrecision, subtreeNode, radius-1, 0);
		}
	if(nd->anc==root){
		if(radius>1 && !(delta < optPrecision)){
			if(nd->next!=NULL){
				if(nd->next->alreadyOptimized == false) scoreIncrease += RecursivelyOptimizeBranches(nd->next, optPrecision, subtreeNode, radius-1, true, 0);
				}
			else if(nd->prev->prev->alreadyOptimized == false) scoreIncrease += RecursivelyOptimizeBranches(nd->prev->prev, optPrecision, subtreeNode, radius-1, true, 0);
			if(nd->prev!=NULL){
				if(nd->prev->alreadyOptimized == false) scoreIncrease += RecursivelyOptimizeBranches(nd->prev, optPrecision, subtreeNode, radius-1, true, 0);
				}
			else if(nd->next->next->alreadyOptimized == false) scoreIncrease += RecursivelyOptimizeBranches(nd->next->next, optPrecision, subtreeNode, radius-1, true, 0);
			}
		}
	if(memLevel > 1) RemoveTempClaReservations();
	return scoreIncrease;
	}

pair<FLOAT_TYPE, FLOAT_TYPE> Tree::CalcDerivativesRateHet(TreeNode *nd1, TreeNode *nd2, FLOAT_TYPE &like){
	//nd1 and nd2 are the nodes on either side of the branch of interest
	//nd1 will always be the "lower" one, and will always be internal, while
	//nd2 can be internal or terminal
	int nsites=data->NChar();
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

	FLOAT_TYPE ***deriv1, ***deriv2, ***prmat;

	ProfModDeriv.Start();
	mod->CalcDerivatives(nd2->dlen, prmat, deriv1, deriv2);
	ProfModDeriv.Stop();

	FLOAT_TYPE d1=ZERO_POINT_ZERO, d2=ZERO_POINT_ZERO;

	if(nd2->left == NULL){
		char *childData=nd2->tipData;
		ProfTermDeriv.Start();
#ifdef OPEN_MP
		if(modSpec.IsNucleotide() == false)
			GetDerivsPartialTerminalNState(claOne, **prmat, **deriv1, **deriv2, childData, d1, d2, like);
		else
			GetDerivsPartialTerminal(claOne, **prmat, **deriv1, **deriv2, childData, d1, d2, nd2->ambigMap);
#else
		if(modSpec.IsNucleotide() == false)
			GetDerivsPartialTerminalNState(claOne, **prmat, **deriv1, **deriv2, childData, d1, d2, like);
		else
			GetDerivsPartialTerminal(claOne, **prmat, **deriv1, **deriv2, childData, d1, d2);
#endif
		ProfTermDeriv.Stop();
		}
	else {
		ProfIntDeriv.Start();
#ifdef EQUIV_CALCS
		GetDerivsPartialInternalEQUIV(claOne, claTwo, **prmat, **deriv1, **deriv2, d1, d2, nd2->tipData);
#else
		if(modSpec.IsNucleotide() == false)
			GetDerivsPartialInternalNState(claOne, claTwo, **prmat, **deriv1, **deriv2, d1, d2, like);
		else
			GetDerivsPartialInternal(claOne, claTwo, **prmat, **deriv1, **deriv2, d1, d2);
#endif
		ProfIntDeriv.Stop();
		}

	return pair<FLOAT_TYPE, FLOAT_TYPE>(d1, d2);
}

FLOAT_TYPE Tree::BranchLike(TreeNode *optNode){

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
			ConditionalLikelihoodRateHet(ROOT, optNode->anc);
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
	
	FLOAT_TYPE initialLen=nd->dlen;
	Score();
	
	out << nd->dlen << "\t" << lnL << "\n";
	
	SetBranchLength(nd, (FLOAT_TYPE)1e-4);
	for(int i=0;i<15;i++){
		Score();
		out << nd->dlen << "\t" << lnL << "\n";
		SetBranchLength(nd, nd->dlen * (FLOAT_TYPE)2.0);	
		}
	SetBranchLength(nd, initialLen);	
	} 

void Tree::CalcEmpiricalDerivatives(TreeNode *nd, FLOAT_TYPE &D1, FLOAT_TYPE &D2){
	FLOAT_TYPE start_blen = nd->dlen;

	FLOAT_TYPE incr;
	FLOAT_TYPE blen_used = start_blen;

	if(blen_used > 1e-6)
		incr= blen_used * 0.001;
	else if(blen_used > min_brlen + 1.0e-8) incr = blen_used * 0.01; 
	else{
		incr =1.0e-8;
		blen_used = max(start_blen, min_brlen + incr);
		SetBranchLength(nd, blen_used);
		}
	
	MakeAllNodesDirty();
	Score();
	FLOAT_TYPE start=lnL;

	SetBranchLength(nd, blen_used + incr);
//	MakeAllNodesDirty();
	Score();

	FLOAT_TYPE empD11= (lnL - start)/incr;

//	SetBranchLength(nd, prevDLen);
//	Score();

	SetBranchLength(nd, blen_used - incr);
//	MakeAllNodesDirty();
	Score();

	FLOAT_TYPE empD12 = (lnL - start)/-incr;

	D1=(empD11+empD12)*.5;
	D2=(empD11-empD12)/incr;

	SetBranchLength(nd, start_blen);
//	MakeAllNodesDirty();
	}

#ifdef SPOOF_NEWTON_RAPHSON
//this allows the ability to play with optimization without actually disrupting program flow
FLOAT_TYPE Tree::NewtonRaphsonOptimizeBranchLength(FLOAT_TYPE precision1, TreeNode *nd, bool goodGuess){
	FLOAT_TYPE origLen =  nd->dlen;
	Score();
	FLOAT_TYPE origScore = lnL;

	FLOAT_TYPE estNRImprove = NewtonRaphsonSpoof(precision1, nd, goodGuess);

	FLOAT_TYPE nrLen = nd->dlen;
	Score();
	FLOAT_TYPE nrScore = lnL;


	SetBranchLength(nd, origLen);
	FLOAT_TYPE estBestImprove = NewtonRaphsonSpoof(0.0001, nd, goodGuess);

	FLOAT_TYPE bestLen = nd->dlen;
	Score();
	FLOAT_TYPE bestScore = lnL;
	
	FLOAT_TYPE trueNRImprove = nrScore - origScore;
	FLOAT_TYPE trueBestImprove = bestScore - origScore;

	SetBranchLength(nd, nrLen);

	ofstream spoof("optspoof.log", ios::app);
	spoof << nd->nodeNum << "\t" << origLen << "\t" << origScore << "\t" << goodGuess << "\n";
	spoof << "\t" << nrLen << "\t" << nrScore << "\t" << estNRImprove << "\t" << trueNRImprove << "\n";
	spoof << "\t" << bestLen << "\t" << bestScore << "\t" << estBestImprove << "\t" << trueBestImprove << "\n";
	spoof.close();

	return estNRImprove;
}

FLOAT_TYPE Tree::NewtonRaphsonSpoof(FLOAT_TYPE precision1, TreeNode *nd, bool goodGuess){
#else
 FLOAT_TYPE Tree::NewtonRaphsonOptimizeBranchLength(FLOAT_TYPE precision1, TreeNode *nd, bool goodGuess){
#endif
/*	if(goodGuess==false && (nd->dlen < 0.0001 || nd->dlen > .1)){
		SetBranchLength(nd, (FLOAT_TYPE).001);
		}
*/
//if(nd->dlen==min_brlen){
#ifdef OPT_DEBUG
/*
if(nd->nodeNum == 8){
	ofstream scr("NRcurve.log");
	scr.precision(20);
	assert(scr.good());
	scr.precision(15);
	FLOAT_TYPE initDlen = nd->dlen;
	for(FLOAT_TYPE d=1e-8;d<.5;d*=1.33){
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
#endif
//	nd->dlen=.3254;
//	SweepDirtynessOverTree(nd);
/*
	if(nd->nodeNum==4){	
		ofstream deb("curves.log");
		SampleBlenCurve(nd, deb);
		deb.close();
		 }
*/	

//	MakeAllNodesDirty();
/*	FLOAT_TYPE start, empD11, empD12, empD1, empD2;
	 if(nd->nodeNum == 4 && nd->dlen < 1e-7 && nd->anc->nodeNum==12 && nd->next != NULL && nd->next->nodeNum==10){
	 //if(0){
		SetBranchLength(nd, 1.0e-5);
		MakeAllNodesDirty();
		Score();
		 }
*/

#ifdef OPT_DEBUG
//	ofstream log("optimization.log", ios::app);
//	log.precision(10);

	opt << nd->nodeNum << "\t" << nd->dlen << "\t" << lnL <<endl;
	
//	ofstream scr("impVSd1.log", ios::app);
	if(lnL > -2){
		Score(nd->anc->nodeNum);
		}

	FLOAT_TYPE delta;
	
#endif
	FLOAT_TYPE totalEstImprove=ZERO_POINT_ZERO;
	int iter=0;
	FLOAT_TYPE abs_d1_prev=FLT_MAX;
	const FLOAT_TYPE v_onEntry=nd->dlen;
	FLOAT_TYPE v=nd->dlen;
	FLOAT_TYPE v_prev = nd->dlen;				/* in case we don't like the new value (see below) */
	bool moveOn = false;
	FLOAT_TYPE prevScore=lnL;
	FLOAT_TYPE curScore=lnL;
	int negProposalNum=0;
	FLOAT_TYPE knownMin=min_brlen, knownMax=max_brlen;
	FLOAT_TYPE d1, d2, estScoreDelta, estDeltaNR;

	FLOAT_TYPE L, initialL;
	do{
		bool scoreOK;
		int sweeps=0;
		pair<FLOAT_TYPE, FLOAT_TYPE> derivs;
		do{		//this part just catches the exception that could be thrown by the rescaling 
				//function if it decides that the current rescaleEvery is too large
			try{
				scoreOK=true;
				derivs = CalcDerivativesRateHet(nd->anc, nd, L);
				if(iter == 0) initialL = L;
				}catch(int err){
				scoreOK=false;
				if(err==1){
					MakeAllNodesDirty();
					rescaleEvery -= 2;
					ofstream resc("rescale.log", ios::app);
					resc << "rescale reduced to " << rescaleEvery << endl;
					resc.close();
					if(rescaleEvery<2) throw(ErrorException("Problem with rescaling in branchlength optimization.  Please report this error to zwickl@nescent.org.\nDetails: nd=%d init=%f cur=%f prev=%d iter=%d neg=%d", nd->nodeNum, v_onEntry, v_prev, nd->dlen, iter, negProposalNum));
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
		
		d1=derivs.first;
		d2=derivs.second;
		//DEBUG
		optCalcs++;

#ifdef OPT_DEBUG
		FLOAT_TYPE empD1, empD2;
//		if(nd->nodeNum == 4){// && nd->anc->nodeNum==12 && nd->next != NULL && nd->next->nodeNum==10){
			CalcEmpiricalDerivatives(nd, empD1, empD2);
			opt << empD1 << "\t" << empD2 << "\t" << nd->dlen + (-empD1/empD2) << "\t";
//			d1 = empD1;
//			d2 = empD2;
//			}
#endif

		estDeltaNR=-d1/d2;
		//estimated change in score by a Taylor series
		estScoreDelta = d1*estDeltaNR + (d2 * estDeltaNR * estDeltaNR * ZERO_POINT_FIVE);

		if(d1 <= ZERO_POINT_ZERO && nd->dlen < knownMax) knownMax = nd->dlen;
		else if(d1 > ZERO_POINT_ZERO && nd->dlen > knownMin) knownMin = nd->dlen;

		//DEBUG
/*		ofstream poo("poo.log", ios::app);
		poo << nd->nodeNum << "\t" << iter << "\t" << d1 << "\t" << d2 << "\t" << estScoreDelta << "\t" << initialL << "\t" << L << endl;
		poo.close();
*/
														#ifdef OPT_DEBUG			
														opt << d1 << "\t" << d2 << "\t" << estScoreDelta << "\t";		
														#endif
		FLOAT_TYPE abs_d1 = fabs(d1);
		if (d2 >= ZERO_POINT_ZERO){//curvature is wrong for NR use 
			//this does NOT only happen when the peak is at the min, as I used to think
														#ifdef OPT_DEBUG			
														opt << "d2 > 0\t";				
														#endif
			//Not allowing this escape anymore
			if(fabs(d1) < ONE_POINT_ZERO){//don't bother doing anything if the surface is this flat
														#ifdef OPT_DEBUG			
														opt << "very small d1.  Return.\n";				
														#endif				
//				return totalEstImprove;
				}

			if(d1 <= ZERO_POINT_ZERO){//if d1 is negative, try shortening arbitrarily, or go halfway to the knownMin
				FLOAT_TYPE proposed;
#ifdef SINGLE_PRECISION_FLOATS
				if(knownMin == min_brlen){
					if(nd->dlen <= 1.0e-4f) proposed = min_brlen;
					else if(nd->dlen <= 0.005f) proposed = 1.0e-4f;
					else if(nd->dlen <= 0.05f) proposed = nd->dlen * 0.1f;
					else proposed = nd->dlen * 0.25f;
					}
#else
				if(knownMin == min_brlen){
					if(nd->dlen <= 1.0e-4) proposed = min_brlen;
					else if(nd->dlen <= 0.005) proposed = 1.0e-4;
					else if(nd->dlen <= 0.05) proposed = nd->dlen * 0.1;
					else proposed = nd->dlen * 0.25;
					}
#endif
				else proposed = (knownMin + nd->dlen) * ZERO_POINT_FIVE;

				if(iter > 0 || proposed == min_brlen){//don't let this bail out on the first iteration based on the estimated
					//change if we are jumping to an arbitrary point, because we are just trying to get to a point
					//where we can actually trust the derivs
					FLOAT_TYPE estImp = d1*(proposed - nd->dlen) + (d2 * (proposed - nd->dlen) * (proposed - nd->dlen) * ZERO_POINT_FIVE);
					if(estImp < precision1){
						#ifdef OPT_DEBUG
						opt << "imp to proposed " << proposed << " < prec, return\n";
						#endif
						return totalEstImprove;
						}
					}
				v=proposed;
				totalEstImprove += precision1;
				}

			else{//d1 > 0.0, try increasing the blen by 10 or 2 if knownMax==max_brlen, otherwise try a step half-way to the knownMax
				FLOAT_TYPE proposed;
				if(knownMax == max_brlen){
#ifdef SINGLE_PRECISION_FLOATS
					if(nd->dlen < 0.1f) proposed = nd->dlen * 10.0f;
					else proposed = min(nd->dlen * 2.0f, max_brlen);
#else
					if(nd->dlen < 0.1) proposed = nd->dlen * 10.0;
					else proposed = min(nd->dlen * 2.0, max_brlen);
#endif
				}
				else proposed = (knownMax + nd->dlen) * ZERO_POINT_FIVE;

				if(iter > 0){//don't let this bail out on the first iteration based on the estimated
					//change if we are jumping to an arbitrary point, because we are just trying to get to a point
					//where we can actually trust the derivs
					FLOAT_TYPE estImp = d1*(proposed - nd->dlen) + (d2 * (proposed - nd->dlen) * (proposed - nd->dlen) * ZERO_POINT_FIVE);
					if(estImp < precision1){
						#ifdef OPT_DEBUG
						opt << "imp to prop < prec, return\n";
						#endif
						return totalEstImprove;
						}
					}
				v=proposed;
				totalEstImprove += precision1;				
				}
			}
		else{//trying NR is feasible
			//DEBUG - keep going if the abs_d1 increased, regardless of what the estScoreDelta is
#define CHECK_LIKE

#ifdef IGNORE_INCREASED_D1
			if(estScoreDelta < precision1){
			
#elif defined(IGNORE_SMALL_INCREASED_D1)
			if(estScoreDelta < precision1 && (iter == 0 || abs_d1 < abs_d1_prev * 2.0) && L >= initialL){

#elif defined(CHECK_LIKE)
			//if(estScoreDelta < precision1 && L < initialL) outman.UserMessage("Worsened!: %f %f", initialL, L);
			
			if(estScoreDelta < precision1 && (mod->NStates() > 4 ? L >= initialL : 1)){
	//		if(estScoreDelta < precision1 && L >= initialL){
#else
			if(estScoreDelta < precision1 && abs_d1 < abs_d1_prev){
#endif
														#ifdef OPT_DEBUG			
														opt << "delta < prec, return\n";
														if(curScore==-ONE_POINT_ZERO){
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
			if (v <= knownMin){
				negProposalNum++;
				if(knownMin == min_brlen){
					FLOAT_TYPE deltaToMin=min_brlen-nd->dlen;
					FLOAT_TYPE scoreDeltaToMin = (deltaToMin * d1 + (deltaToMin*deltaToMin*d2*ZERO_POINT_FIVE));
					if(scoreDeltaToMin < precision1){
													#ifdef OPT_DEBUG
													opt << "imp to MIN < prec, return\n";			
													#endif
						return totalEstImprove;
						}
#ifdef SINGLE_PRECISION_FLOATS
					else if(negProposalNum==1 && nd->dlen > 1e-4f && v_prev != 1e-4f){
						//try a somewhat smaller length before going all the way to the min
						if(nd->dlen < .005f ) v = 1e-4f;
						else if(nd->dlen < 0.05f) v = nd->dlen * 0.1f;
						else v = nd->dlen * .25f;
#else
					else if(negProposalNum==1 && nd->dlen > 1e-4 && v_prev != 1e-4){
						//try a somewhat smaller length before going all the way to the min
						if(nd->dlen < .005 ) v = 1e-4;
						else if(nd->dlen < 0.05) v = nd->dlen * 0.1;
						else v = nd->dlen * 0.25;
#endif
						FLOAT_TYPE delta=v - nd->dlen;
						totalEstImprove += (delta * d1 + (delta*delta*d2*ZERO_POINT_FIVE));
						}
					else{
						v = min_brlen;
						totalEstImprove += scoreDeltaToMin;
						}
					}
				else{//knownMin is > absolute min, so we must already have a better guess
					//go half way to that guess
					FLOAT_TYPE proposed = (knownMin + nd->dlen) * ZERO_POINT_FIVE;
					FLOAT_TYPE deltaToMin=proposed-nd->dlen;
					FLOAT_TYPE scoreDeltaToMin = (deltaToMin * d1 + (deltaToMin*deltaToMin*d2*ZERO_POINT_FIVE));
					if(scoreDeltaToMin < precision1){
						#ifdef OPT_DEBUG
						opt << "imp to knownMIN < prec, return\n";
						#endif
						return totalEstImprove;
						}
					v=proposed;
					totalEstImprove += scoreDeltaToMin;				
					}
				}
			
			else if (v >= knownMax){
				if(knownMax == max_brlen){
					FLOAT_TYPE deltaToMax=max_brlen - nd->dlen;
					FLOAT_TYPE scoreDeltaToMax = (deltaToMax * d1 + (deltaToMax*deltaToMax*d2*ZERO_POINT_FIVE));
					if(scoreDeltaToMax < precision1){
						#ifdef OPT_DEBUG
						opt << "imp to MAX < prec, return\n";
						#endif
						return totalEstImprove;
						}
					else{
						v = max_brlen;
						totalEstImprove += scoreDeltaToMax;
						}
					}
				else{//knownMax is < absolute max, so we must already have a better guess
					//go half way to that guess
					FLOAT_TYPE proposed = (knownMax + nd->dlen) * ZERO_POINT_FIVE;
					FLOAT_TYPE deltaToMax=proposed-nd->dlen;
					FLOAT_TYPE scoreDeltaToMax = (deltaToMax * d1 + (deltaToMax*deltaToMax*d2*ZERO_POINT_FIVE));
					if(scoreDeltaToMax < precision1){
						#ifdef OPT_DEBUG
						opt << "imp to knownMAX < prec, return\n";
						#endif
						return totalEstImprove;
						}
					v=proposed;
					totalEstImprove += scoreDeltaToMax;
					}
				}
			else totalEstImprove += estScoreDelta;

			abs_d1_prev = abs_d1;
			}
		assert(v >= min_brlen);
		assert(v >= knownMin);
		assert(v <= knownMax);

		SetBranchLength(nd, v);

#ifdef OPT_DEBUG
		Score(nd->anc->nodeNum);

		if(curScore != -ONE_POINT_ZERO){
			if(lnL < curScore){
				cout << lnL << "\t" << curScore << endl;
				if(curScore - lnL < .005){
					//don't want to have different logic when OPT_DEBUG is on
//					SetBranchLength(nd, v_prev);
//					Score(nd->anc->nodeNum);
//					return lnL;
					}
				else {//assert(0);
					FLOAT_TYPE poo=lnL;
					SetBranchLength(nd, v_prev);
					MakeAllNodesDirty();
					Score(nd->anc->nodeNum);
					assert(fabs(prevScore - lnL) < .01);
					
					SetBranchLength(nd, v);
					MakeAllNodesDirty();
					Score(nd->anc->nodeNum);
					assert(fabs(poo - lnL) < .01);
					}
				}
			}
				
		curScore=lnL;	
		delta=prevScore - lnL;

opt << v << "\t" << lnL << "\n";
opt.flush();
#endif

		prevScore=lnL;
		v_prev=v;
		
		iter++;
		if(iter>50){
/*			ofstream deb("optdeb.log");
			deb << "initial length " << v_onEntry << endl;
			deb << "current length " << nd->dlen << endl;
			deb << "prev length " << v_prev << endl;
			deb << "d1 " << d1 << " d2 " << d2 << endl;
			deb << "neg proposal num " << negProposalNum << endl;			
			deb.close();
*/
			throw(ErrorException("Problem with branchlength optimization.  Please report this error to zwickl@nescent.org.\nDetails: nd=%d init=%f cur=%f prev=%d d1=%f d2=%f neg=%d", nd->nodeNum, v_onEntry, v_prev, nd->dlen, d1, d2, negProposalNum));
/*
				ofstream scr("NRcurve.log");
				scr.precision(20);
				assert(scr.good());
				scr.precision(15);
				FLOAT_TYPE initDlen = nd->dlen;
				for(FLOAT_TYPE d=1e-8;d<.5;d*=1.33){
					nd->dlen = d;
					SweepDirtynessOverTree(nd);
					Score();
					scr << d << "\t" << lnL << endl;
					}
				nd->dlen=initDlen;
				SweepDirtynessOverTree(nd);
				scr.close();

			bool poo=true;
			outman.UserMessage("long opt: %d", iter);
			ofstream deb("longopt.log", ios::app);
			deb << iter << "\t" << precision1 << "\t" << nd->nodeNum << "\t" << v_onEntry << "\t" << nd->dlen << "\t" << d1 << "\t" << d2 << "\t" << estDeltaNR << "\t" << estScoreDelta << "\t" << negProposalNum << "\n";
			deb.close();
*/
//			while(poo){
//				opt.close();
//				}
			//assert(iter<=50);
			}		
		}while(moveOn==false);
#ifdef OPT_DEBUG
	opt << "final\t" << nd->dlen << "\t" << lnL << endl;
#endif
	assert(0);//shouldn't be exiting this way
	return totalEstImprove;
	}
/*
void Tree::RecursivelyOptimizeBranches(TreeNode *nd, FLOAT_TYPE optPrecision, int subtreeNode, int radius, int centerNode, bool dontGoNext){
	FLOAT_TYPE prevScore=lnL;
#ifdef BRENT	
	BrentOptimizeBranchLength(optPrecision, nd, false);
	FLOAT_TYPE delta=lnL - prevScore;
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

void Tree::RecursivelyOptimizeBranchesDown(TreeNode *nd, TreeNode *calledFrom, FLOAT_TYPE optPrecision, int subtreeNode, int radius, int ){
	FLOAT_TYPE prevScore=lnL;
#ifdef BRENT	
	BrentOptimizeBranchLength(optPrecision, nd, false);
	FLOAT_TYPE delta=lnL - prevScore;
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
void Tree::OptimizeBranchesAroundNode(TreeNode *nd, FLOAT_TYPE optPrecision, int subtreeNode){
	//this function will optimize the three branches (2 descendents and one anc) connecting
	//to it.  It assumes that everything that is dirty has been marked so.
	//by default there is only a single optimization pass over the three nodes
	FLOAT_TYPE precision1, precision2;

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
inline FLOAT_TYPE CallBranchLike(TreeNode *thisnode, Tree *thistree, FLOAT_TYPE blen){
	thisnode->dlen=exp(blen);
	return thistree->BranchLike(thisnode)*-1;
	}
	
inline FLOAT_TYPE CallBranchLikeRateHet(TreeNode *thisnode, Tree *thistree, FLOAT_TYPE blen){

	thisnode->dlen=blen;
	FLOAT_TYPE like=thistree->BranchLikeRateHet(thisnode)*-1;

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

FLOAT_TYPE Tree::BrentOptimizeBranchLength(FLOAT_TYPE accuracy_cutoff, TreeNode *here, bool firstPass){
	//we pass the node whose branch length whose blen we want to optimize, but note that the
	//calculations occur at the node below that
	//if firstPass is true, we have no idea what a reasonable value for the blen is, so
	//use a wide bracket.  If it is false, try a fairly tight bracket around the current val
	FLOAT_TYPE a, b, c, fa, fb, fc, minimum, minScore=0.0;
	FLOAT_TYPE blen=here->dlen;

	assert(blen>=min_brlen);

	if(here->anc){
		if(firstPass){
			if(blen<1e-6){
				a=min_brlen;
				if(blen!=min_brlen){
					b=blen;
					}
				else{
					b=min_brlen*100;
					lnL=-1;
					}
				c=min_brlen*10000.0;
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
				a=min_brlen;
				if(blen!=min_brlen){
					b=blen;
					}
				else{
					b=min_brlen*100;
					lnL=-1;
					}
				c=min_brlen*10000.0;				
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
		FLOAT_TYPE min_len=minimum;		
//		FLOAT_TYPE min_len=exp(minimum);
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

FLOAT_TYPE Tree::BrentOptimizeBranchLength(FLOAT_TYPE accuracy_cutoff, TreeNode *here, bool goodGuess){
	//we pass the node whose branch length whose blen we want to optimize, but note that the
	//calculations occur at the node below that
	//if firstPass is true, we have no idea what a reasonable value for the blen is, so
	//use a wide bracket.  If it is false, try a fairly tight bracket around the current val
	FLOAT_TYPE a, b, c, fa, fb, fc, minimum, minScore=ZERO_POINT_ZERO;
	FLOAT_TYPE blen=here->dlen;
	
	FLOAT_TYPE min_len;
	
	assert(blen>=min_brlen);

	FLOAT_TYPE initialScore;
	fb=initialScore=CallBranchLike(here, this, sqrt(sqrt(here->dlen)), true);

if(here->anc){
#ifndef FOURTH_ROOT
//	if(here->anc){
		if(firstPass){
			if(!(blen>1e-6)){
				a=min_brlen;
				if(blen!=min_brlen){
					b=blen;
					}
				else{
					b=min_brlen*100;
					lnL=-1;
					}
				c=min_brlen*10000.0;
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
				a=min_brlen;
				if(blen!=min_brlen){
					b=blen;
					}
				else{
					b=min_brlen*100;
					lnL=-1;
					}
				c=min_brlen*10000.0;				
				}			
			}
#endif
#ifdef FOURTH_ROOT
/*
	if(blen < min_brlen*10){
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
*/
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
			FLOAT_TYPE sqrtmin=sqrt(min_brlen);
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

void Tree::GetDerivsPartialTerminal(const CondLikeArray *partialCLA, const FLOAT_TYPE *prmat, const FLOAT_TYPE *d1mat, const FLOAT_TYPE *d2mat, const char *Ldat, FLOAT_TYPE &d1Tot, FLOAT_TYPE &d2Tot, const unsigned *ambigMap /*=NULL*/){

	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.
	FLOAT_TYPE *partial=partialCLA->arr;
	const int nchar=data->NChar();
	const int nRateCats=mod->NRateCats();

#ifdef UNIX
	madvise(partial, nchar*4*nRateCats*sizeof(FLOAT_TYPE), MADV_SEQUENTIAL);
#endif

	FLOAT_TYPE siteL;
	FLOAT_TYPE La, Lc, Lg, Lt;
	FLOAT_TYPE D1a, D1c, D1g, D1t;
	FLOAT_TYPE D2a, D2c, D2g, D2t;
	FLOAT_TYPE tot1=ZERO_POINT_ZERO, tot2=ZERO_POINT_ZERO;//can't use d1Tot and d2Tot in OMP reduction because they are references

	const char *Ldata=Ldat;
	const int *countit=data->GetCounts();
	const FLOAT_TYPE *rateProb=mod->GetRateProbs();
	const int lastConst=data->LastConstant();
	const int *conBases=data->GetConstStates();
	const FLOAT_TYPE prI=mod->PropInvar();

	FLOAT_TYPE freqs[4];
	for(int i=0;i<4;i++) freqs[i]=mod->StateFreq(i);

#ifdef OMP_TERMDERIV
	#pragma omp parallel for private(La, Lc, Lg, Lt, D1a, D1c, D1g, D1t, D2a, D2c, D2g, D2t, partial, Ldata, siteL) reduction(+ : tot1, tot2)
	for(int i=0;i<nchar;i++){
		Ldata = &Ldat[ambigMap[i]];
		partial = &partialCLA->arr[i*4*nRateCats];
#else
	for(int i=0;i<nchar;i++){
#endif
#ifdef USE_COUNTS_IN_BOOT
		if(countit[i] > 0){
#else
		if(1){
#endif
			La=Lc=Lg=Lt=D1a=D1c=D1g=D1t=D2a=D2c=D2g=D2t=ZERO_POINT_ZERO;
			if(*Ldata > -1){ //no ambiguity
				for(int r=0;r<nRateCats;r++){
					La  += prmat[(*Ldata)+16*r] * partial[0] * rateProb[r];
					D1a += d1mat[(*Ldata)+16*r] * partial[0] * rateProb[r];
					D2a += d2mat[(*Ldata)+16*r] * partial[0] * rateProb[r];
					Lc  += prmat[(*Ldata+4)+16*r] * partial[1]* rateProb[r];
					D1c += d1mat[(*Ldata+4)+16*r] * partial[1]* rateProb[r];
					D2c += d2mat[(*Ldata+4)+16*r] * partial[1]* rateProb[r];
					Lg  += prmat[(*Ldata+8)+16*r] * partial[2]* rateProb[r];
					D1g += d1mat[(*Ldata+8)+16*r] * partial[2]* rateProb[r];
					D2g += d2mat[(*Ldata+8)+16*r] * partial[2]* rateProb[r];
					Lt  += prmat[(*Ldata+12)+16*r] * partial[3]* rateProb[r];
					D1t += d1mat[(*Ldata+12)+16*r] * partial[3]* rateProb[r];
					D2t += d2mat[(*Ldata+12)+16*r] * partial[3]* rateProb[r];
					partial += 4;
					}
				Ldata++;
				}
				
			else if(*Ldata == -4){ //total ambiguity
				for(int r=0;r<nRateCats;r++){
					La += partial[0]* rateProb[r];
					Lc += partial[1]* rateProb[r];
					Lg += partial[2]* rateProb[r];
					Lt += partial[3]* rateProb[r];
					partial += 4;
					}
				Ldata++;
				}
			else{ //partial ambiguity
				int nstates=-1 * *(Ldata++);
				for(int s=0;s<nstates;s++){
					for(int r=0;r<nRateCats;r++){
						La += prmat[(*Ldata)+16*r]  * partial[4*r]* rateProb[r];
						D1a += d1mat[(*Ldata)+16*r] * partial[4*r]* rateProb[r];		
						D2a += d2mat[(*Ldata)+16*r] * partial[4*r]* rateProb[r];
											
						Lc += prmat[(*Ldata+4)+16*r] * partial[1+4*r]* rateProb[r];
						D1c += d1mat[(*Ldata+4)+16*r]* partial[1+4*r]* rateProb[r];
						D2c += d2mat[(*Ldata+4)+16*r]* partial[1+4*r]* rateProb[r];
											
						Lg += prmat[(*Ldata+8)+16*r]* partial[2+4*r]* rateProb[r];
						D1g += d1mat[(*Ldata+8)+16*r]* partial[2+4*r]* rateProb[r];
						D2g += d2mat[(*Ldata+8)+16*r]* partial[2+4*r]* rateProb[r];
						
						Lt += prmat[(*Ldata+12)+16*r]* partial[3+4*r]* rateProb[r];
						D1t += d1mat[(*Ldata+12)+16*r]* partial[3+4*r]* rateProb[r];
						D2t += d2mat[(*Ldata+12)+16*r]* partial[3+4*r]* rateProb[r];
						}
					Ldata++;
					}
				partial+=4*nRateCats;
				}
			if((mod->NoPinvInModel() == false) && (i<=lastConst)){
				FLOAT_TYPE btot=ZERO_POINT_ZERO;
				if(conBases[i]&1) btot+=freqs[0];
				if(conBases[i]&2) btot+=freqs[1];
				if(conBases[i]&4) btot+=freqs[2];
				if(conBases[i]&8) btot+=freqs[3];
				//6-27-05 fixed this to calc derivs correctly if constant site has been rescaled
				siteL  = ((La*freqs[0]+Lc*freqs[1]+Lg*freqs[2]+Lt*freqs[3]) + (prI*btot)*exp((FLOAT_TYPE)partialCLA->underflow_mult[i]));
				}
			else
				siteL  = ((La*freqs[0]+Lc*freqs[1]+Lg*freqs[2]+Lt*freqs[3]));

			assert(La > 0.0f || Lc > 0.0f || Lg > 0.0f || Lt > 0.0f);
			assert(La < 1.0e30 && Lc < 1.0e30 && Lg < 1.0e30 && Lt < 1.0e30);

#ifdef SINGLE_PRECISION_FLOATS
			FLOAT_TYPE tempD1 = (((D1a*freqs[0]+D1c*freqs[1]+D1g*freqs[2]+D1t*freqs[3])) / siteL);
			if(fabs(tempD1) < 1.0e8f){
				tot1+= countit[i] * tempD1;
				FLOAT_TYPE siteD2=((D2a*freqs[0]+D2c*freqs[1]+D2g*freqs[2]+D2t*freqs[3]));
				tot2 += countit[i] * ((siteD2 / siteL) - tempD1*tempD1);
				}
#else
			FLOAT_TYPE tempD1 = (((D1a*freqs[0]+D1c*freqs[1]+D1g*freqs[2]+D1t*freqs[3])) / siteL);
			tot1+= countit[i] * tempD1;
			FLOAT_TYPE siteD2=((D2a*freqs[0]+D2c*freqs[1]+D2g*freqs[2]+D2t*freqs[3]));
			tot2 += countit[i] * ((siteD2 / siteL) - tempD1*tempD1);
#endif
//			assert(tot1 < 1.0e10 && tot2 < 1.0e10);
			}
#ifndef OMP_TERMDERIV
		else{
			//partial+=4*nRateCats;
			if(!(*Ldata < 0)) Ldata++;
			else if(*Ldata == -4) Ldata++;
			else{
				int nstates=-1 * *(Ldata++);
				for(int s=0;s<nstates;s++) Ldata++;
				}
			}
#endif
		}

	d1Tot = tot1;
	d2Tot = tot2;
	}
	
void Tree::GetDerivsPartialTerminalNState(const CondLikeArray *partialCLA, const FLOAT_TYPE *prmat, const FLOAT_TYPE *d1mat, const FLOAT_TYPE *d2mat, const char *Ldat, FLOAT_TYPE &d1Tot, FLOAT_TYPE &d2Tot, FLOAT_TYPE &like){
	//this function assumes that the pmat is arranged with nstates^2 entries for the
	//first rate, followed by nstates^2 for the second, etc.
	FLOAT_TYPE *partial=partialCLA->arr;
	const int nRateCats=mod->NRateCats();

	int nchar = data->NChar();
	const int *countit = data->GetCounts();
	int nstates = mod->NStates();
	const char *Ldata = Ldat;

	const FLOAT_TYPE *rateProb=mod->GetRateProbs();
	const int lastConst=data->LastConstant();
	const int *conStates=data->GetConstStates();
	const FLOAT_TYPE prI=mod->PropInvar();

#ifdef UNIX
	madvise(partial, nchar*nstates*nRateCats*sizeof(FLOAT_TYPE), MADV_SEQUENTIAL);
#endif

	FLOAT_TYPE tot1=ZERO_POINT_ZERO, tot2=ZERO_POINT_ZERO, totL=ZERO_POINT_ZERO;//can't use d1Tot and d2Tot in OMP reduction because they are references

	FLOAT_TYPE *freqs = new FLOAT_TYPE[nstates];
	for(int i=0;i<nstates;i++) freqs[i]=mod->StateFreq(i);

	FLOAT_TYPE siteL, siteD1, siteD2;

#undef OUTPUT_DERIVS

#ifdef OUTPUT_DERIVS
	ofstream out("derivs.log");

	ofstream pout("pmat.log");
	ofstream d1out("d1mat.log");
	ofstream d2out("d2mat.log");

	for(int a = 0;a < nstates;a++){
		for(int b = 0;b < nstates;b++){
			pout << prmat[b + nstates*a] << "\t";
			d1out << d1mat[b + nstates*a] << "\t";
			d2out << d2mat[b + nstates*a] << "\t";
			}
		pout << endl;
		d1out << endl;
		d2out << endl;
		}
	pout.close();
	d1out.close();
	d2out.close();

#endif

	if(nRateCats == 1){

	#ifdef OMP_TERMDERIV_NSTATE
		#pragma omp parallel for private(partial, Ldata, siteL, siteD1, siteD2) reduction(+ : tot1, tot2, totL)
		for(int i=0;i<nchar;i++){
			Ldata = &Ldat[i];
			partial = &partialCLA->arr[i*nstates*nRateCats];
	#else
		for(int i=0;i<nchar;i++){
	#endif
			siteL = ZERO_POINT_ZERO;
			siteD1 = ZERO_POINT_ZERO;
			siteD2 = ZERO_POINT_ZERO;
#ifdef USE_COUNTS_IN_BOOT
			if((countit[i]) > 0){//this check speeds us up in the case of bootstrapping
#else
			if(1){
#endif
				if(*Ldata != nstates){ //no ambiguity
					for(int from=0;from<nstates;from++){
						siteL += prmat[(*Ldata)+nstates*from] * partial[from] * freqs[from];
						siteD1 += d1mat[(*Ldata)+nstates*from] * partial[from] * freqs[from];
						siteD2 += d2mat[(*Ldata)+nstates*from] * partial[from] * freqs[from];
						}
					}
					
				else if(*Ldata == nstates){ //total ambiguity
					for(int from=0;from<nstates;from++){
						siteL += partial[from] * freqs[from];
						}
					}
				else{ //partial ambiguity
					assert(0);
					}
				siteL *= rateProb[0]; //multiply by (1-pinv)
				
				if((mod->NoPinvInModel() == false) && (i<=lastConst)){
					siteL += (prI*freqs[conStates[i]] * exp((FLOAT_TYPE)partialCLA->underflow_mult[i]));
					}
				
				totL += log(siteL) * countit[i];
				
				siteD1 /= siteL;
				tot1 += countit[i] * siteD1;
				tot2 += countit[i] * ((siteD2 / siteL) - siteD1*siteD1);
				assert(tot1 == tot1);
				assert(tot2 == tot2);

#ifdef OUTPUT_DERIVS
				out << "L\t" << siteL << "\tD1\t" << siteD1 << "\tD2\t" << siteD2 << "\tcount\t" << countit[i] << "\tunderM\t" << partialCLA->underflow_mult[i] << endl;
#endif

				Ldata++;
				partial+=nstates*nRateCats;
				}
			else{
#ifdef OPEN_MP
				//this is a little strange, but the arrays needs to be advanced in the case of OMP (if this function is not OMP enabled)
				//because sections of the CLAs corresponding to sites with count=0 are skipped
				//over in OMP instead of being eliminated
				partial+=nstates*nRateCats;
#endif
				Ldata++;
				}
			}
		}
	else{

#ifdef OMP_TERMDERIV_NSTATE
		#pragma omp parallel for private(partial, Ldata, siteL, siteD1, siteD2) reduction(+ : tot1, tot2, totL)
		for(int i=0;i<nchar;i++){
			Ldata = &Ldat[i];
			partial = &partialCLA->arr[i*nstates*nRateCats];
#else
		for(int i=0;i<nchar;i++){
#endif
#ifdef USE_COUNTS_IN_BOOT
			if((countit[i]) != 0){//this check speeds us up in the case of bootstrapping
#else
			if(1){
#endif
				siteL = ZERO_POINT_ZERO;
				siteD1 = ZERO_POINT_ZERO;
				siteD2 = ZERO_POINT_ZERO;
				vector<FLOAT_TYPE> stateContribL(nstates);
				vector<FLOAT_TYPE> stateContribD1(nstates);
				vector<FLOAT_TYPE> stateContribD2(nstates);
				for(int from=0;from<nstates;from++){
					stateContribL[from] = 0.0;
					stateContribD1[from] = 0.0;
					stateContribD2[from] = 0.0;
					}

				if(*Ldata < nstates){ //no ambiguity
					for(int rate=0;rate<nRateCats;rate++){
						for(int from=0;from<nstates;from++){
							stateContribL[from] += prmat[rate*nstates*nstates + (*Ldata)+nstates*from] * partial[from] * rateProb[rate];
							stateContribD1[from] += d1mat[rate*nstates*nstates + (*Ldata)+nstates*from] * partial[from] * rateProb[rate];
							stateContribD2[from] += d2mat[rate*nstates*nstates + (*Ldata)+nstates*from] * partial[from] * rateProb[rate];
							}
						partial += nstates;
						}
					Ldata++;
					}
					
				else if(*Ldata == nstates){ //total ambiguity
					for(int rate=0;rate<nRateCats;rate++){
						for(int from=0;from<nstates;from++){
							stateContribL[from] += partial[from];
							}
						partial += nstates;
						}
					Ldata++;
					}
				else{ //partial ambiguity
					assert(0);
					}

				for(int from=0;from<nstates;from++){
					siteL += stateContribL[from] * freqs[from];
					siteD1 += stateContribD1[from] * freqs[from];
					siteD2 += stateContribD2[from] * freqs[from];
					}

				if((mod->NoPinvInModel() == false) && (i<=lastConst)){
					siteL += (prI*freqs[conStates[i]] * exp((FLOAT_TYPE)partialCLA->underflow_mult[i]));
					}

				totL += log(siteL) * countit[i];
				siteD1 /= siteL;
				tot1 += countit[i] * siteD1;
				tot2 += countit[i] * ((siteD2 / siteL) - siteD1*siteD1);
				assert(siteL == siteL);
				assert(totL == totL);
				assert(tot1 == tot1);
				assert(tot2 == tot2);
				}
			else{
#ifdef OPEN_MP	//this needs to be advanced in the case of openmp, regardless of whether
				//this function actually has OMP enabled or not.
				partial += nstates * nRateCats;
#endif
				Ldata++;
				}
			}
		}

#ifdef OUTPUT_DERIVS
	out.close();
#endif

	d1Tot = tot1;
	d2Tot = tot2;
	like = totL;
	delete []freqs;
	}

void Tree::GetDerivsPartialInternal(const CondLikeArray *partialCLA, const CondLikeArray *childCLA, const FLOAT_TYPE *prmat, const FLOAT_TYPE *d1mat, const FLOAT_TYPE *d2mat, FLOAT_TYPE &d1Tot, FLOAT_TYPE &d2Tot){
	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.
	FLOAT_TYPE *CL1=childCLA->arr;
	FLOAT_TYPE *partial=partialCLA->arr;

	const int nchar=data->NChar();
	const int nRateCats=mod->NRateCats();

#ifdef UNIX
	madvise(partial, nchar*4*nRateCats*sizeof(FLOAT_TYPE), MADV_SEQUENTIAL);
	madvise(CL1, nchar*4*nRateCats*sizeof(FLOAT_TYPE), MADV_SEQUENTIAL);
#endif

	FLOAT_TYPE siteL;
	FLOAT_TYPE La, Lc, Lg, Lt;
	FLOAT_TYPE D1a, D1c, D1g, D1t;
	FLOAT_TYPE D2a, D2c, D2g, D2t;
	FLOAT_TYPE tot1=ZERO_POINT_ZERO, tot2=ZERO_POINT_ZERO;//can't use d1Tot and d2Tot in OMP reduction because they are references
	
	const int *countit=data->GetCounts();
	const FLOAT_TYPE *rateProb=mod->GetRateProbs();
		
	const int lastConst=data->LastConstant();
	const int *conBases=data->GetConstStates();
	const FLOAT_TYPE prI=mod->PropInvar();	

	FLOAT_TYPE freqs[4];
	for(int i=0;i<4;i++) freqs[i]=mod->StateFreq(i);

#ifdef OMP_INTDERIV
#pragma omp parallel for private(La, Lc, Lg, Lt, D1a, D1c, D1g, D1t, D2a, D2c, D2g, D2t, partial, CL1, siteL) reduction(+ : tot1, tot2)
	for(int i=0;i<nchar;i++){
		partial = &(partialCLA->arr[4*i*nRateCats]);
		CL1		= &(childCLA->arr[4*i*nRateCats]);
#else
	for(int i=0;i<nchar;i++){
#endif
#ifdef USE_COUNTS_IN_BOOT
		if(countit[i] > 0){
#else
		if(1){
#endif
			La=Lc=Lg=Lt=D1a=D1c=D1g=D1t=D2a=D2c=D2g=D2t=ZERO_POINT_ZERO;
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
				FLOAT_TYPE	btot=ZERO_POINT_ZERO;
				if(conBases[i]&1) btot+=freqs[0];
				if(conBases[i]&2) btot+=freqs[1];
				if(conBases[i]&4) btot+=freqs[2];
				if(conBases[i]&8) btot+=freqs[3];
				//6-27-05 fixed this to calc derivs correctly if constant site has been rescaled
				siteL  = ((La*freqs[0]+Lc*freqs[1]+Lg*freqs[2]+Lt*freqs[3]) + (prI*btot)*exp((FLOAT_TYPE)childCLA->underflow_mult[i]+partialCLA->underflow_mult[i]));
				}
			else
				siteL  = ((La*freqs[0]+Lc*freqs[1]+Lg*freqs[2]+Lt*freqs[3]));
			FLOAT_TYPE tempD1 = (((D1a*freqs[0]+D1c*freqs[1]+D1g*freqs[2]+D1t*freqs[3])) / siteL);
#ifdef SINGLE_PRECISION_FLOATS
			if(fabs(tempD1) < 1.0e8f){
				assert(d1Tot == d1Tot);
				FLOAT_TYPE siteD2=((D2a*freqs[0]+D2c*freqs[1]+D2g*freqs[2]+D2t*freqs[3]));
				tot1 += countit[i] * tempD1;
				tot2 += countit[i] * ((siteD2 / siteL) - tempD1*tempD1);
				}
#else
			assert(d1Tot == d1Tot);
			FLOAT_TYPE siteD2=((D2a*freqs[0]+D2c*freqs[1]+D2g*freqs[2]+D2t*freqs[3]));
			tot1 += countit[i] * tempD1;
			tot2 += countit[i] * ((siteD2 / siteL) - tempD1*tempD1);
#endif
			assert(d2Tot == d2Tot);
//			assert(tot1 < 1.0e10 && tot2 < 1.0e10);
			}
#ifndef OMP_INTDERIV
		else{
	//		partial+=4*nRateCats;
	//		CL1+=4*nRateCats;			
			}
#endif
		}
	d1Tot = tot1;
	d2Tot = tot2;
	}

void Tree::GetDerivsPartialInternalNState(const CondLikeArray *partialCLA, const CondLikeArray *childCLA, const FLOAT_TYPE *prmat, const FLOAT_TYPE *d1mat, const FLOAT_TYPE *d2mat, FLOAT_TYPE &d1Tot, FLOAT_TYPE &d2Tot, FLOAT_TYPE &like){
	//this function assumes that the pmat is arranged with the nstates^2 entries for the
	//first rate, followed by nstates^2 for the second, etc.
	FLOAT_TYPE *CL1=childCLA->arr;
	FLOAT_TYPE *partial=partialCLA->arr;

	int nchar = data->NChar();
	const int *countit = data->GetCounts();
	int nstates = mod->NStates();
	const int nRateCats = mod->NRateCats();

	const FLOAT_TYPE *rateProb=mod->GetRateProbs();

	const int lastConst=data->LastConstant();
	const int *conStates=data->GetConstStates();
	const FLOAT_TYPE prI=mod->PropInvar();	

#ifdef UNIX
	madvise(partial, nchar*nstates*nRateCats*sizeof(FLOAT_TYPE), MADV_SEQUENTIAL);
	madvise(CL1, nchar*nstates*nRateCats*sizeof(FLOAT_TYPE), MADV_SEQUENTIAL);
#endif

	FLOAT_TYPE tot1=ZERO_POINT_ZERO, tot2=ZERO_POINT_ZERO, totL = ZERO_POINT_ZERO;//can't use d1Tot and d2Tot in OMP reduction because they are references

	FLOAT_TYPE *freqs = new FLOAT_TYPE[nstates];
	for(int i=0;i<nstates;i++) freqs[i]=mod->StateFreq(i);

	FLOAT_TYPE siteL = 0.0;
	FLOAT_TYPE siteD1 = 0.0;
	FLOAT_TYPE siteD2 = 0.0;
	FLOAT_TYPE tempL = 0.0;
	FLOAT_TYPE tempD1 = 0.0;
	FLOAT_TYPE tempD2 = 0.0;

	if(nRateCats == 1){
#ifdef OMP_INTDERIV_NSTATE
	#pragma omp parallel for private(siteL, siteD1, siteD2, tempL, tempD1, tempD2, partial, CL1) reduction(+ : tot1, tot2, totL)
		for(int i=0;i<nchar;i++){
			partial = &(partialCLA->arr[nstates*i]);
			CL1		= &(childCLA->arr[nstates*i]);
#else
		for(int i=0;i<nchar;i++){
#endif
#ifdef USE_COUNTS_IN_BOOT
			if(countit[i] > 0){//this check speeds us up in the case of bootstrapping
#else
			if(1){
#endif
				siteL = 0.0;
				siteD1 = 0.0;
				siteD2 = 0.0;
				for(int from=0;from<nstates;from++){
					tempL = prmat[from*nstates]*CL1[0];
					tempD1 = d1mat[from*nstates]*CL1[0];
					tempD2 = d2mat[from*nstates]*CL1[0];
					for(int to=1;to<nstates;to++){
						tempL += prmat[from*nstates + to]*CL1[to];
						tempD1 += d1mat[from*nstates + to]*CL1[to];
						tempD2 += d2mat[from*nstates + to]*CL1[to];
						}
					siteL += tempL * partial[from] * freqs[from];
					siteD1 += tempD1 * partial[from] * freqs[from];
					siteD2 += tempD2 * partial[from] * freqs[from];
					}
				siteL *= rateProb[0]; //multiply by (1-pinv)

				if((mod->NoPinvInModel() == false) && (i<=lastConst)){
					siteL += (prI*freqs[conStates[i]] * exp((FLOAT_TYPE)partialCLA->underflow_mult[i]));
					}

				totL += log(siteL) * countit[i];
				siteD1 /= siteL;
				tot1 += countit[i] * siteD1;
				tot2 += countit[i] * ((siteD2 / siteL) - siteD1*siteD1);
				assert(tot1 == tot1);
				assert(tot2 == tot2);

				partial += nstates;
				CL1 += nstates;
				}
			}
		}
	else{

#ifdef OMP_INTDERIV_NSTATE
		#pragma omp parallel for private(siteL, siteD1, siteD2, tempL, tempD1, tempD2, partial, CL1) reduction(+ : tot1, tot2, totL)
		for(int i=0;i<nchar;i++){
			partial = &(partialCLA->arr[nRateCats * nstates * i]);
			CL1		= &(childCLA->arr[nRateCats * nstates * i]);
#else
		for(int i=0;i<nchar;i++){
#endif
#ifdef USE_COUNTS_IN_BOOT
			if(countit[i] > 0){//this check speeds us up in the case of bootstrapping
#else
			if(1){
#endif
				siteL = 0.0;
				siteD1 = 0.0;
				siteD2 = 0.0;

				vector<FLOAT_TYPE> stateContribL(nstates);
				vector<FLOAT_TYPE> stateContribD1(nstates);
				vector<FLOAT_TYPE> stateContribD2(nstates);

				for(int from=0;from<nstates;from++){
					stateContribL[from] = 0.0;
					stateContribD1[from] = 0.0;
					stateContribD2[from] = 0.0;
					}
				for(int rate=0;rate<nRateCats;rate++){
					for(int from=0;from<nstates;from++){
						tempL = prmat[rate*nstates*nstates + from*nstates]*CL1[0];
						tempD1 = d1mat[rate*nstates*nstates + from*nstates]*CL1[0];
						tempD2 = d2mat[rate*nstates*nstates + from*nstates]*CL1[0];
						for(int to=1;to<nstates;to++){
							tempL += prmat[rate*nstates*nstates + from*nstates + to]*CL1[to];
							tempD1 += d1mat[rate*nstates*nstates + from*nstates + to]*CL1[to];
							tempD2 += d2mat[rate*nstates*nstates + from*nstates + to]*CL1[to];
							}
						stateContribL[from] += tempL * partial[from] * rateProb[rate];
						stateContribD1[from] += tempD1 * partial[from] * rateProb[rate];
						stateContribD2[from] += tempD2 * partial[from] * rateProb[rate];
						}
					partial += nstates;
					CL1 += nstates;
					}
				for(int from=0;from<nstates;from++){
					siteL += stateContribL[from] * freqs[from];
					siteD1 += stateContribD1[from] * freqs[from];
					siteD2 += stateContribD2[from] * freqs[from];
					}
				if((mod->NoPinvInModel() == false) && (i<=lastConst)){
					siteL += (prI*freqs[conStates[i]] * exp((FLOAT_TYPE)partialCLA->underflow_mult[i]));
					}
				totL += log(siteL) * countit[i];
				siteD1 /= siteL;
				tot1 += countit[i] * siteD1;
				tot2 += countit[i] * ((siteD2 / siteL) - siteD1*siteD1);
				assert(tot1 == tot1);
				assert(tot2 == tot2);
				}
			}
		}

	d1Tot = tot1;
	d2Tot = tot2;
	like = totL;
	delete []freqs;
	}

void Tree::GetDerivsPartialInternalEQUIV(const CondLikeArray *partialCLA, const CondLikeArray *childCLA, const FLOAT_TYPE *prmat, const FLOAT_TYPE *d1mat, const FLOAT_TYPE *d2mat, FLOAT_TYPE &d1Tot, FLOAT_TYPE &d2Tot, char *equiv){
	//this function assumes that the pmat is arranged with the 16 entries for the
	//first rate, followed by 16 for the second, etc.
	FLOAT_TYPE *CL1=childCLA->arr;
	FLOAT_TYPE *partial=partialCLA->arr;

	const int nchar=data->NChar();
	const int nRateCats=mod->NRateCats();

#ifdef UNIX
	madvise(partial, nchar*4*nRateCats*sizeof(FLOAT_TYPE), MADV_SEQUENTIAL);
	madvise(CL1, nchar*4*nRateCats*sizeof(FLOAT_TYPE), MADV_SEQUENTIAL);
#endif

	FLOAT_TYPE siteL;
	FLOAT_TYPE La, Lc, Lg, Lt;
	FLOAT_TYPE D1a, D1c, D1g, D1t;
	FLOAT_TYPE D2a, D2c, D2g, D2t;

	FLOAT_TYPE eLa, eLc, eLg, eLt;
	FLOAT_TYPE eD1a, eD1c, eD1g, eD1t;
	FLOAT_TYPE eD2a, eD2c, eD2g, eD2t;

	FLOAT_TYPE tot1=ZERO_POINT_ZERO, tot2=ZERO_POINT_ZERO;//can't use d1Tot and d2Tot in OMP reduction because they are references
	
	const int *countit=data->GetCounts();
	const FLOAT_TYPE *rateProb=mod->GetRateProbs();
	assert(nRateCats  == 1);

	const int lastConst=data->LastConstant();
	const int *conBases=data->GetConstStates();
	const FLOAT_TYPE prI=mod->PropInvar();	

	FLOAT_TYPE freqs[4];
	for(int i=0;i<4;i++) freqs[i]=mod->StateFreq(i);

	int rOff =0;
	for(int i=0;i<nchar;i++){
		if((countit[i]) != 0){
			if(equiv[i] == false){
				eLa = ( prmat[rOff ]*CL1[0]+prmat[rOff + 1]*CL1[1]+prmat[rOff + 2]*CL1[2]+prmat[rOff + 3]*CL1[3])  * rateProb[0];
				eLc = ( prmat[rOff + 4]*CL1[0]+prmat[rOff + 5]*CL1[1]+prmat[rOff + 6]*CL1[2]+prmat[rOff + 7]*CL1[3])  * rateProb[0];
				eLg = ( prmat[rOff + 8]*CL1[0]+prmat[rOff + 9]*CL1[1]+prmat[rOff + 10]*CL1[2]+prmat[rOff + 11]*CL1[3])  * rateProb[0];
				eLt = ( prmat[rOff + 12]*CL1[0]+prmat[rOff + 13]*CL1[1]+prmat[rOff + 14]*CL1[2]+prmat[rOff + 15]*CL1[3])  * rateProb[0];
				
				eD1a = ( d1mat[rOff ]*CL1[0]+d1mat[rOff + 1]*CL1[1]+d1mat[rOff + 2]*CL1[2]+d1mat[rOff + 3]*CL1[3])  * rateProb[0];
				eD1c = ( d1mat[rOff + 4]*CL1[0]+d1mat[rOff + 5]*CL1[1]+d1mat[rOff + 6]*CL1[2]+d1mat[rOff + 7]*CL1[3])  * rateProb[0];
				eD1g = ( d1mat[rOff + 8]*CL1[0]+d1mat[rOff + 9]*CL1[1]+d1mat[rOff + 10]*CL1[2]+d1mat[rOff + 11]*CL1[3])  * rateProb[0];
				eD1t = ( d1mat[rOff + 12]*CL1[0]+d1mat[rOff + 13]*CL1[1]+d1mat[rOff + 14]*CL1[2]+d1mat[rOff + 15]*CL1[3])  * rateProb[0];		

				eD2a = ( d2mat[rOff ]*CL1[0]+d2mat[rOff + 1]*CL1[1]+d2mat[rOff + 2]*CL1[2]+d2mat[rOff + 3]*CL1[3])  * rateProb[0];
				eD2c = ( d2mat[rOff + 4]*CL1[0]+d2mat[rOff + 5]*CL1[1]+d2mat[rOff + 6]*CL1[2]+d2mat[rOff + 7]*CL1[3])  * rateProb[0];
				eD2g = ( d2mat[rOff + 8]*CL1[0]+d2mat[rOff + 9]*CL1[1]+d2mat[rOff + 10]*CL1[2]+d2mat[rOff + 11]*CL1[3])  * rateProb[0];
				eD2t = ( d2mat[rOff + 12]*CL1[0]+d2mat[rOff + 13]*CL1[1]+d2mat[rOff + 14]*CL1[2]+d2mat[rOff + 15]*CL1[3])  * rateProb[0];
				}
			La = eLa * partial[0];
			Lc = eLc * partial[1];
			Lg = eLg * partial[2];
			Lt = eLt * partial[3];

			D1a = eD1a * partial[0];
			D1c = eD1c * partial[1];
			D1g = eD1g * partial[2];
			D1t = eD1t * partial[3];

			D2a = eD2a * partial[0];
			D2c = eD2c * partial[1];
			D2g = eD2g * partial[2];
			D2t = eD2t * partial[3];

			partial+=4;
			CL1+=4;
			
			if((mod->NoPinvInModel() == false) && (i<=lastConst)){
				FLOAT_TYPE	btot=ZERO_POINT_ZERO;
				if(conBases[i]&1) btot+=freqs[0];
				if(conBases[i]&2) btot+=freqs[1];
				if(conBases[i]&4) btot+=freqs[2];
				if(conBases[i]&8) btot+=freqs[3];
				//6-27-05 fixed this to calc derivs correctly if constant site has been rescaled
				siteL  = ((La*freqs[0]+Lc*freqs[1]+Lg*freqs[2]+Lt*freqs[3]) + (prI*btot)*exp((FLOAT_TYPE)childCLA->underflow_mult[i]+partialCLA->underflow_mult[i]));
				}
			else
				siteL  = ((La*freqs[0]+Lc*freqs[1]+Lg*freqs[2]+Lt*freqs[3]));
			FLOAT_TYPE tempD1 = (((D1a*freqs[0]+D1c*freqs[1]+D1g*freqs[2]+D1t*freqs[3])) / siteL);
#ifdef SINGLE_PRECISION_FLOATS
			if(fabs(tempD1) < 1.0e8f){
				assert(d1Tot == d1Tot);
				FLOAT_TYPE siteD2=((D2a*freqs[0]+D2c*freqs[1]+D2g*freqs[2]+D2t*freqs[3]));
				tot1 += countit[i] * tempD1;
				tot2 += countit[i] * ((siteD2 / siteL) - tempD1*tempD1);
				}
#else
			assert(d1Tot == d1Tot);
			FLOAT_TYPE siteD2=((D2a*freqs[0]+D2c*freqs[1]+D2g*freqs[2]+D2t*freqs[3]));
			tot1 += countit[i] * tempD1;
			tot2 += countit[i] * ((siteD2 / siteL) - tempD1*tempD1);
#endif
			assert(d2Tot == d2Tot);
//			assert(tot1 < 1.0e10 && tot2 < 1.0e10);
			}
		else{
			partial+=4*nRateCats;
			CL1+=4*nRateCats;			
			}
		}
	d1Tot = tot1;
	d2Tot = tot2;
	}
