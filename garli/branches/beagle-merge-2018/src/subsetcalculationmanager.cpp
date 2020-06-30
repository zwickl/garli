#include "defs.h"
#include "treenode.h"
#include "utility.h"
#include "calculationmanager.h"

ClaManager* SubsetCalculationManager::claMan = NULL;
PmatManager* SubsetCalculationManager::pmatMan = NULL;

#define BATCHED_CALLS

void SubsetCalculationManager::CheckBeagleReturnValue(int err, const char *funcName) const{

	//BPART TODO, need to find a better place for this 
	bool termOnBeagleError = true;

	string mess;
	if (err == BEAGLE_SUCCESS || err >= 0)
		return;

	else if (err == BEAGLE_ERROR_GENERAL)
		mess = "General Beagle error in ";

	else if (err == BEAGLE_ERROR_OUT_OF_MEMORY)
		mess = "Beagle out-of-memory error in ";

	else if (err == BEAGLE_ERROR_UNIDENTIFIED_EXCEPTION)
		mess = "Beagle unidentified exception error in ";

	else if (err == BEAGLE_ERROR_UNINITIALIZED_INSTANCE)
		mess = "Beagle uninitialized instance error in ";

	else if (err == BEAGLE_ERROR_OUT_OF_RANGE)
		mess = "Beagle out-of-range error in ";

	else if (err == BEAGLE_ERROR_NO_RESOURCE)
		mess = "Beagle was unable to find a suitable resource in ";

	else if (err == BEAGLE_ERROR_NO_IMPLEMENTATION)
		mess = "Beagle no-implementation error in ";

	else if (err == BEAGLE_ERROR_FLOATING_POINT)
		mess = "Beagle floating point error in ";

	else
		mess = "Unknown Beagle error in ";

	mess += funcName;
	if (termOnBeagleError)
		throw ErrorException("%s", mess.c_str());

	else
		outman.UserMessage("%s", mess.c_str());

	return;
}


//this could be improved by calculating some of these arguments in this function 
BeagleInstanceDetails SubsetCalculationManager::InitializeSubset(int nClas, int nHolders, int pref_flags, int req_flags, SequenceData *data, ModelSpecification *modSpec, int modelInd) {
	//set these class members
	modelIndex = modelInd;
	subsetModSpec = modSpec;
	subsetData = data;
	ntax = data->NTax();
	nstates = subsetModSpec->nstates;
	nchar = data->NChar();
	nrates = (subsetModSpec->numRateCats + (subsetModSpec->includeInvariantSites ? 1 : 0));

	//to allow ambiguity we need to create partials for the tips with ambiguity, and normal tip states for the others
	ambigTips = subsetData->NumTaxaWithPartialAmbig();
	compactCount = ntax - ambigTips;
	partialsCount = nClas + ambigTips;

#ifndef DONT_SEND_TRANSMATS
	//num freqs = numEigens in beagle, so need at least one
	int eigCount = 1;
#else
	//BMERGE - don't know why this was nclas, at most should be # categories (nrates)
	//int eigCount = nClas;
	int eigCount = nrates;
#endif
	//this is a gross overestimate of how many matrices are needed, but haven't worked out recycling quite yet for mats
	int matrixCount = (nHolders * 3) * 2;//x3 for the pmats, d1mats and d2mats, x2 for both internals and terms
	//int matrixCount = (nClas * 3) * 2;//x3 for the pmats, d1mats and d2mats, x2 for both internals and terms

	//try one scaler per cla, as in normal garli.  These are doubles rather than ints though, so larger.
	//scaler for a given cla will share same index, and will be cumulative, as mine are now
	//add one for a destinationScaleWrite which will always be the last and will be used in all calls for scratch
	//Note that although there is a 1to1 correspondence between the cla and scaler indeces, this is NOT true for what is
	//passed to beagle, so beagle partial index = my partial index + #tips, beagle scaler = my scaler

	//DEBUG - it is hard to tell whether we will later be forced in a SP instance, in which case we'd like to turn on rescaling
	//so, just allocate enough scalers to do it whether we will or not
	//int scalerCount = (IsRescaling() ? nClas + 1 : 0);
	int scalerCount = nClas + 1;
	int resourceList[1] = { 0 };
	int resourceListCount = 0;

	/*TODO_BEAGLE_PART
	if (beagleDeviceNum > -1) {
		resourceList[0] = beagleDeviceNum;
		resourceListCount = 1;
	}
	*/

	BeagleInstanceDetails det;

	//this returns either the instance number or a beagle error, which is somewhat annoying
	int beagleInstNum = beagleCreateInstance(ntax,
		partialsCount,
		compactCount,
		nstates,
		nchar,
		eigCount,
		matrixCount,
		nrates,
		scalerCount,
		resourceList,
		resourceListCount,
		pref_flags,
		req_flags,
		&det);

	beagleInst = beagleInstNum;

	CheckBeagleReturnValue(
		beagleInstNum,
		"beagleCreateInstance");

	ParseInstanceDetails(&det);

	/*
	int fpSize = (IsSinglePrecision() ? 4 : 8);

	outman.DebugMessage("BEAGLE ALLOCATIONS:");
	outman.DebugMessage("\tstates: %d char: %d rates: %d", nstates, nchar, nrates);
	outman.DebugMessage("\ttips: %d", tipCount);
	outman.DebugMessage("\tcompact tips arrays (no partial ambiguity): %d (%.1f KB)", compactCount, compactCount * nchar * 4 / 1024.0);
	outman.DebugMessage("\tpartial tips arrays (some partial ambiguity): %d (%.1f KB)", ambigTips, ambigTips * nstates * nchar * fpSize / 1024.0);
	outman.DebugMessage("\ttotal partial arrays: %d (%.1f KB)", partialsCount, partialsCount * nstates * nchar * nrates * fpSize / 1024.0);
	outman.DebugMessage("\teigen solutions: %d (%.1f KB (?))", eigCount, eigCount * nstates * 2 * fpSize / 1024.0);
	outman.DebugMessage("\ttransition matrices (pmats and derivs): %d (%.1f KB)", matrixCount, matrixCount * nstates * nstates * nrates * fpSize / 1024.0);
	outman.DebugMessage("\trescaling arrays %d (%.1f KB)", scalerCount, scalerCount * nchar * fpSize / 1024.0);

	double totMem = (compactCount * nchar * 4 / 1024.0) + (ambigTips * nstates * nchar * fpSize / 1024.0) + (partialsCount * nstates * nchar * nrates * fpSize / 1024.0) +
		(eigCount * nstates * 2 * fpSize / 1024.0) + (eigCount, eigCount * nstates * 2 * fpSize / 1024.0) + (matrixCount * nstates * nstates * nrates * fpSize / 1024.0) +
		(scalerCount * nchar * fpSize / 1024.0);

	outman.DebugMessage("Total beagle memory allocations approx %.1f MB)", totMem / 1024.0);

	vector<double> counts;
	for (int pat = 0; pat < data->NChar(); pat++)
		counts.push_back((double)data->Count(pat));
#ifdef OUTPUT_OTHER_BEAGLE
	outman.DebugMessageNoCR("counts ");
	for (vector<double>::iterator it = counts.begin(); it != counts.end(); it++) {
		outman.DebugMessageNoCR("%d ", (int)*it);
	}
	outman.DebugMessage("");
#endif
#ifdef OUTPUT_BEAGLE_SITELIKES
	string name = ofprefix + ".Bsitelikes.log";
	ofstream file(name.c_str());
	file.close();
#endif
#ifdef OUTPUT_BEAGLE_PARTIALS
	string name2 = ofprefix + ".partials.log";
	ofstream part(name2.c_str());
	part.close();
#endif
	*/
	vector<double> counts;
	for (int pat = 0; pat < subsetData->NChar(); pat++)
		counts.push_back((double)subsetData->Count(pat));
#ifdef OUTPUT_OTHER_BEAGLE
	outman.DebugMessageNoCR("counts ");
	for (vector<double>::iterator it = counts.begin(); it != counts.end(); it++) {
		outman.DebugMessageNoCR("%d ", (int)* it);
	}
	outman.DebugMessage("");
#endif
	beagleSetPatternWeights(beagleInstNum, &(counts[0]));
	SendTipDataToBeagle();

	outman.UserMessage("#######################################################");
	
	return det;
}


void SubsetCalculationManager::SendTipDataToBeagle() {
	outman.DebugMessage("SENDING DATA");
	char convert[16] = { -1, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };

	for (int t = 0; t < ntax; t++) {
		bool partialAmbig = subsetData->TaxonHasPartialAmbig(t);
		outman.DebugMessageNoCR("t %d ", t);

		const unsigned char* dataString = subsetData->GetRow(t);
		if (partialAmbig) {
			//ambiguity if currently for nuc only (all versions, not just beagle)
			assert(nstates == 4);
			outman.DebugMessage("(some ambig)");
			vector<double> tipPartial;
			for (int c = 0; c < subsetData->NChar(); c++) {
				for (int s = 0; s < nstates; s++) {
					((dataString[c] & (1 << s)) ? tipPartial.push_back(1.0) : tipPartial.push_back(0.0));
				}
			}
			CheckBeagleReturnValue(
				beagleSetTipPartials(beagleInst, t, &(tipPartial[0])),
				"beagleSetTipStates");
		}
		else {
			outman.DebugMessage("(no ambig)");
			vector<int> dat;

			for (int c = 0; c < subsetData->NChar(); c++) {
				if (nstates == 4) {
					//nucleotide data needs to be converted from the bitwise format to 0, 1, 2, 3, 4
					assert(dataString[c] < 16);
					assert(convert[dataString[c]] < 5);
					dat.push_back(convert[dataString[c]]);
				}
				else {
					//for non-nuc the data is already in the correct indexing scheme, but must be converted from
					//chars to ints
					dat.push_back(dataString[c]);
				}
			}
			CheckBeagleReturnValue(
				beagleSetTipStates(beagleInst, t, &(dat[0])),
				"beagleSetTipStates");
		}
	}
	outman.DebugMessage("DATA SENT");
}

void SubsetCalculationManager::PerformTransMatOperationBatch(const list<TransMatOperation> &theOps) {
	//for now assuming that all transmat ops in a set have same eigen solution and rate multipliers

	//currently just using a single eigen index and sending it every time
	int eigenIndex = 0;

	//this call will calculate the eigen solution first if necessary
	//currently it will always be resent (sol.changed will always be true) because
	//it is non-trivial to figure out when it has changed since it was last sent to beagle
	ModelEigenSolution sol;

	pmatMan->GetEigenSolution(theOps.begin()->destTransMatIndex, sol, modelIndex);

	//omega cats
	if (nstates > 60 && nrates > 1) {
		if (sol.changed) {
			for (int cat = 0; cat < nrates; cat++) {
				CheckBeagleReturnValue(
					beagleSetEigenDecomposition(
						beagleInst,
						cat,
						&sol.eigenVecs[cat * nstates * nstates],
						&sol.invEigenVecs[cat * nstates * nstates],
						&sol.eigenVals[cat * nstates]),
					"beagleSetEigenDecomposition");
			}
		}
	}
	else {
		CheckBeagleReturnValue(
			beagleSetEigenDecomposition(
				beagleInst,
				0,
				sol.eigenVecs,
				sol.invEigenVecs,
				sol.eigenVals),
			"beagleSetEigenDecomposition");
	}

	vector<FLOAT_TYPE> categRates;
	pmatMan->GetCategoryRatesForBeagle(theOps.begin()->destTransMatIndex, categRates, modelIndex);

	//as above, always resend currently
	if (categRates.size() > 0) {
		CheckBeagleReturnValue(
			beagleSetCategoryRates(beagleInst,
				&(categRates[0])),
			"beagleSetCategoryRates");
	}

	vector<int> pmatInd;
	vector<int> d1MatInd;
	vector<int> d2MatInd;
	vector<double> edgeLens;

	for (list<TransMatOperation>::const_iterator it = theOps.begin(); it != theOps.end(); it++) {
		//DEBUG
		//outman.DebugMessage("transmat %d mod %d modadr %d edge %.4f", (*it).destTransMatIndex, (*it).modelIndex, pmatMan->GetCorrespondingModel((*it).destTransMatIndex), (*it).edgeLength);
		//pmatMan->GetCorrespondingModel(theOps.begin()->destTransMatIndex)->NumRateCatsForBeagle();
		pmatMan->ClaimTransmatsetFillIfNecessary((*it).destTransMatIndex);
		pmatInd.push_back(PmatIndexForBeagle((*it).destTransMatIndex));
		d1MatInd.push_back(D1MatIndexForBeagle((*it).destTransMatIndex));
		d2MatInd.push_back(D2MatIndexForBeagle((*it).destTransMatIndex));
		edgeLens.push_back(pmatMan->GetScaledEdgelen((*it).destTransMatIndex, modelIndex));
	}

	int count = pmatInd.size();
	bool calcDerivs = theOps.begin()->calcDerivs;

#ifdef OUTPUT_OTHER_BEAGLE
	outman.DebugMessageNoCR("UPDATING TRANS MAT ");
	outman.DebugMessage("%d (%d), eigen %d, blen %f", PmatIndexForBeagle(theOps.front().destTransMatIndex), theOps.front().destTransMatIndex, eigenIndex, edgeLens[0]);
#endif

	if (nstates > 60 && nrates > 1) {
		vector<int> eigenIndeces;
		for (int r = 0; r < nrates; r++) eigenIndeces.push_back(r);
		//eigenIndeces = { 0, 1 };
		CheckBeagleReturnValue(
			beagleUpdateTransitionMatricesWithModelCategories(
				beagleInst,
				&eigenIndeces[0],
				&pmatInd[0],
				(calcDerivs ? &d1MatInd[0] : NULL),
				(calcDerivs ? &d2MatInd[0] : NULL),
				&edgeLens[0],
				count),
			"beagleUpdateTransitionMatricesWithModelCategories");
	}
	else {
		CheckBeagleReturnValue(
			beagleUpdateTransitionMatrices(
				beagleInst,
				eigenIndex,
				&pmatInd[0],
				(calcDerivs ? &d1MatInd[0] : NULL),
				(calcDerivs ? &d2MatInd[0] : NULL),
				&edgeLens[0],
				count),
			"beagleUpdateTransitionMatrices");
	}


#ifdef OUTPUT_PMATS

	ofstream* deb = outman.GetDebugStream();
	//int nstates = data->NStates();

	const Model* localDebugModel = pmatMan->GetModelForTransMatHolder(theOps.begin()->destTransMatIndex, modelIndex);
	//int nrates = pmatMan->GetModelForTransMatHolder(theOps.begin()->destTransMatIndex)->NumRateCatsForBeagle();
//	if(calcDerivs){

	/*wrong	*/
	//double inmat[16] = {0.9070614,	0.02927732,	0.02870608,	0.0349794, 0.03213807,	0.90418679,	0.02868981,	0.03497981, 0.03213807,	0.02926432,	0.90361226, 0.03497981, 0.03213831,	0.02927236,	0.02868767,	0.90990174};
	//double inmat[16] = {0.9070613980, 0.0292773210, 0.0287060775, 0.0349793993, 0.0321380682, 0.9041867852, 0.0286898147, 0.0349798091, 0.0321380682, 0.0292643178, 0.9036122561, 0.0349798091, 0.0321383066, 0.0292723589, 0.0286876746, 0.9099017382}; 
	//wrong, normalized rows sum to zero
	//double inmat[16] = {0.9070394515, 0.0292766126, 0.0287053830, 0.0349785530, 0.0321382457, 0.9041917789, 0.0286899731, 0.0349800023, 0.0321382465, 0.0292644802, 0.9036172701, 0.0349800032, 0.0321383041, 0.0292723566, 0.0286876724, 0.9099016670};
	//wrong, rescaled such that column totals are same, three lower elements in each column match
	//double inmat[16] = {0.9070613980, 0.0292713326, 0.0286945223, 0.0349796725, 0.0321381477, 0.9041867852, 0.0286945223, 0.0349796725, 0.0321381477, 0.0292713326, 0.9036122561, 0.0349796725, 0.0321381477, 0.0292713326, 0.0286945223, 0.9099017382};
	/*right*/
	//double inmat[16] = {0.9070611, 0.0292719, 0.02868779, 0.03497917, 0.03213868, 0.90419436, 0.02868777, 0.03497917, 0.03213868, 0.0292719, 0.90361023, 0.03497917, 0.03213868, 0.02927192, 0.02868777, 0.90990162};
	//double inmat[16] = {0.9070611396, 0.0292718949, 0.0286877828, 0.0349791827, 0.0321386719, 0.9041943625, 0.0286877828, 0.0349791827, 0.0321386719, 0.0292718949, 0.9036102504, 0.0349791827, 0.0321386719, 0.0292718949, 0.0286877828, 0.9099016503};
	//pmat from DP beagle
	double inmat[16] = { 0.907060771,  0.02927234,  0.028687678,  0.034979238,  0.032138303,  0.904194759,  0.028687707,  0.034979216,  0.032138303,  0.029272277,  0.903610189,  0.034979216,  0.032138297,  0.029272304,  0.028687689,  0.909901709 };

	vector<double> outMat(nstates * nstates * nrates);
	for (list<TransMatOperation>::const_iterator pit = theOps.begin(); pit != theOps.end(); pit++) {
		//		outman.UserMessage("MANUALLY SETTING PMATS");
		//		beagleSetTransitionMatrix(beagleInst, PmatIndexForBeagle((*pit).destTransMatIndex), inmat);

		outMat.clear();
		for (int p = 0; p < nstates * nstates * nrates; p++)
			outMat.push_back(p);
		beagleGetTransitionMatrix(beagleInst, PmatIndexForBeagle((*pit).destTransMatIndex), &(outMat[0]));
		//beagleGetTransitionMatrix(beagleInst, D1MatIndexForBeagle((*pit).destTransMatIndex), &(outMat[0]));
		vector<double>::iterator it = outMat.begin();
		outman.DebugMessage("Pmat: transmatIndex %d", (*pit).destTransMatIndex);

		localDebugModel->OutputPmat(*deb, outMat);
		
		for (int r = 0; r < nrates; r++) {
			outman.DebugMessage("rate %d", r);
			for (int i = 0; i < nstates; i++) {
				for (int j = 0; j < nstates; j++) {
					//outman.DebugMessageNoCR("%.1f\t", *it++);
					outman.DebugMessageNoCR("%.6f\t", *it++);
				}
				outman.DebugMessage("");
			}
		}
		if (theOps.front().calcDerivs) {
			beagleGetTransitionMatrix(beagleInst, D1MatIndexForBeagle((*pit).destTransMatIndex), &(outMat[0]));
			it = outMat.begin();
			outman.DebugMessage("D1: transmatIndex %d", (*pit).destTransMatIndex);
			localDebugModel->OutputPmat(*deb, outMat);

			for (int r = 0; r < nrates; r++) {
				outman.DebugMessage("rate %d", r);
				for (int i = 0; i < nstates; i++) {
					for (int j = 0; j < nstates; j++) {
						//outman.DebugMessageNoCR("%.1f\t", *it++);
						outman.DebugMessageNoCR("%.6f\t", *it++);
					}
					outman.DebugMessage("");
				}
			}


			beagleGetTransitionMatrix(beagleInst, D2MatIndexForBeagle((*pit).destTransMatIndex), &(outMat[0]));
			it = outMat.begin();
			outman.DebugMessage("D2: transmatIndex %d", (*pit).destTransMatIndex);
			localDebugModel->OutputPmat(*deb, outMat);

			for (int r = 0; r < nrates; r++) {
				outman.DebugMessage("rate %d", r);
				for (int i = 0; i < nstates; i++) {
					for (int j = 0; j < nstates; j++) {
						//outman.DebugMessageNoCR("%.1f\t", *it++);
						outman.DebugMessageNoCR("%.6f\t", *it++);
					}
					outman.DebugMessage("");
				}
			}
		}
		
	}
	//		}
#endif
}
void SubsetCalculationManager::OutputBeagleSiteValues(ofstream &out, bool derivs) const {

	//there is only one set of sitelikes stored by an instance, the most recent ones		
	vector<double> siteLikesOut(subsetData->NChar());
	vector<double> siteD1Out(subsetData->NChar());
	vector<double> siteD2Out(subsetData->NChar());

	const int *counts = subsetData->GetCounts();

	beagleGetSiteLogLikelihoods(beagleInst, &(siteLikesOut[0]));
	if (derivs)
		beagleGetSiteDerivatives(beagleInst, &(siteD1Out[0]), &(siteD2Out[0]));

	//ScoreSet mySummation = SumSiteValues(&siteLikesOut[0], (derivs ? &siteD1Out[0] : NULL), (derivs ? &siteD2Out[0] : NULL));
	//out << "mine\t" << mySummation.lnL << "\t" << mySummation.d1 << "\t" << mySummation.d2 << "\n";
	//outman.DebugMessage("myL\t%f\tD1\t%f\tD2\t%f", mySummation.lnL, mySummation.d1, mySummation.d2);
	//outman.UserMessage("mine = %.4f beag = %.4f", mySummation.lnL, results.lnL);

	for (int s = 0; s < subsetData->GapsIncludedNChar(); s++) {
		int packed = subsetData->Number(s);
		if (packed >= 0) {
			out << s << "\t" << packed << "\t" << counts[packed] << "\t" << siteLikesOut[packed];
			if (derivs) {
				out << "\t" << siteD1Out[packed] << "\t" << siteD2Out[packed];
			}
		}
		else
			out << s << "\t" << packed << "\t" << counts[packed] << "\tNA";
		out << "\n";
	}
}

//For scale arrays my indexing scheme and Beagle's happen to be the same
//further accumulate the cumulative rescalings from the two children.  In the batched case
//the dest array will already have the rescaling done at this node included, so don't reset it
void SubsetCalculationManager::AccumulateRescalers(int destIndex, int childIndex1, int childIndex2) {

#ifdef OUTPUT_OTHER_BEAGLE
	//DEBUG
	//outman.DebugMessage("accum %d + %d = %d", childIndex1, childIndex2, destIndex);
#endif
#ifndef BATCHED_CALLS
	//clear out the destination scaler first
	CheckBeagleReturnValue(
		beagleResetScaleFactors(beagleInst,
			ScalerIndexForBeagle(destIndex)),
		"beagleResetScaleFactors");
#endif

	vector<int> childScalers;
	if (!(childIndex1 < 0))
		childScalers.push_back(ScalerIndexForBeagle(childIndex1));
	if (!(childIndex2 < 0))
		childScalers.push_back(ScalerIndexForBeagle(childIndex2));

	if (childScalers.size() > 0) {
		CheckBeagleReturnValue(
			beagleAccumulateScaleFactors(
				beagleInst,
				&(childScalers[0]),
				(int)childScalers.size(),
				ScalerIndexForBeagle(destIndex)),
			"beagleAccumulateScaleFactors");
	}
}

//this just performs all of the transmat and cla ops that were queued
//will generally be followed by a call to PerformScoring to get lnL or derivs
void SubsetCalculationManager::UpdateAllConditionals(list<BlockingOperationsSet> operationSetQueue) {
	//DEBUG
#ifdef DEBUG_OPS
	OutputOperationsSummary();
#endif
	
	//this just yoinks the pmat ops from each BlockingOpSet into a new list
	//note that PerformTransMatOperationBatch assumes that all the ops either
	//need derivs or not.  So, avoid adding ops that need derivs here.
	list<TransMatOperation> tops;
	for(list<BlockingOperationsSet>::iterator it = operationSetQueue.begin();it != operationSetQueue.end();it++){
		if(! (*it).pmatOps.begin()->calcDerivs)
			tops.splice(tops.end(), (*it).pmatOps);
		}

	if(tops.size() > 0){
#ifndef DONT_SEND_TRANSMATS
		SendTransMatsToBeagle(tops);
#else
		PerformTransMatOperationBatch(tops);
#endif
		}

	//Any transmats that need derivs will still be included in a blockingOpSet, so will get done here.  Otherwise
	//no transmat calls will be done here because the ops were pulled out above to be done in one big batch
	
	/*Calculation of all necessary partials will go like this:
		0. Before getting here all cla ops will have their cla and pmat dependencies set, and
		   the dependencies of the first destination clas will have been reserved if clean clas
		1. (calc any transmats necessary - currently this will actually happen above)
		2. Claim the first set of destination clas (this could in theory cause a call to Recycle, but
		   shouldn't be necessary unless something is wrong, I think)
		3. Calc the first level
		4. Unreserve deps of first level (which would be clas clean to begin with), don't add to freeable
		5. Claim second (could call recyc, but wouldn't touch anything used in this entire scoring/deriv call)
		6. Calc second
		7. Mark dependencies of second (i.e., the first level) as freeable and push to claman freeable queue
		6. Claim third (could call recyc, AND would first reclaim the first level from claman freeable queue)
		7. Calc third
		6. Mark dependencies of third (i.e., the second level) as freeable and push to claman freeable queue
		7. Claim fourth (could call recyc, AND would first reclaim the first and second level from claman freeable queue)
		6. Calc fourth
		8. Mark deps of fourth as freeable and push to claman freeable queue
		etc...
		9. Unreserve anything remaining, clear freeable queue and op queue

	*/

	claMan->SetCurReclaimableLevel(0);
	list<BlockingOperationsSet>::iterator toFree = operationSetQueue.begin();

	for (list<BlockingOperationsSet>::iterator it = operationSetQueue.begin(); it != operationSetQueue.end();) {
		if ((*it).pmatOps.size() > 0) {
#ifndef DONT_SEND_TRANSMATS
			SendTransMatsToBeagle((*it).pmatOps);
#else
			PerformTransMatOperationBatch((*it).pmatOps);
#endif
		}
		if ((*it).claOps.size() > 0) {
			//this will reserve the DESTINATION clas. The child cla dependencies should already have
			//been reserved if they are clas (not tips).  If they were clean in the first place they
			//should have been reserved in AccumulateOps, otherwise when they were reserved when previous
			//destinations
			ReserveDestinationClas(*it);
			PerformClaOperationBatch((*it).claOps);

			//No need to do all this unreservation if the number of clas is large enough to simultaneously hold
			//all clas of the tree. We'll still unreserve them all at the end
			//DEBUG
			//if(memLevel > 0){
			if (1) {
				//remove the temp reservations from the level 1 DEPENDENCIES (if they are clean clas and not tips)
				//but DO NOT add them to the freeable queue (since they must have been clean in the first place)
				if ((*it).opSetDepLevel == 1) {
					UnreserveDependencyClas(*it);
					claMan->SetCurReclaimableLevel((*it).opSetDepLevel - 1);
					it = operationSetQueue.erase(it);
				}
				//when we've finished level 2 or greater, we can now tell the claMan that the dependencies are ripe for the taking
				else if ((*it).opSetDepLevel > 1) {
					QueueDependencyClasForFreeing(*it);
					claMan->SetCurReclaimableLevel((*it).opSetDepLevel - 1);
					it = operationSetQueue.erase(it);
				}
			}
			else it++;
		}
		else
			it++;
	}

//TODO - what should happen with reclaim level after all the calcs?  
//trying to allow freeing of first level dependencies in general
claMan->SetCurReclaimableLevel(1);
}

	void SubsetCalculationManager::PerformClaOperationBatch(const list<ClaOperation> &theOps) {
		//not sure if this is right - will always use a single scale array for destWrite (essentially
		//scratch space, I think) and then pass a cumulative scaler to actually keep track of the scaling
		//int destinationScaleWrite = (rescaleBeagle ? claMan->NumClas() : BEAGLE_OP_NONE);

		//scale read is if you are telling it to use precomputed scaling factors
		int	destinationScaleRead = BEAGLE_OP_NONE;

		vector<int> tupleList;
		for (list<ClaOperation>::const_iterator it = theOps.begin(); it != theOps.end(); it++) {
			//This is happening a level higher now
			//this will fill the destination holder with a cla index if necessary
			//claMan->ClaimClaFillIfNecessary((*it).destClaIndex, (*it).depLevel);

			tupleList.push_back(PartialIndexForBeagle((*it).destClaIndex));
			//if rescaling, for each partial op pass in an array to store the amount of rescaling
			//at this node, in my case stored as logs
			if (IsRescaling())
				tupleList.push_back(ScalerIndexForBeagle((*it).destClaIndex));
			else
				tupleList.push_back(BEAGLE_OP_NONE);
			tupleList.push_back(destinationScaleRead);
			tupleList.push_back(PartialIndexForBeagle((*it).childClaIndex1));
			tupleList.push_back(PmatIndexForBeagle((*it).transMatIndex1));
			tupleList.push_back(PartialIndexForBeagle((*it).childClaIndex2));
			tupleList.push_back(PmatIndexForBeagle((*it).transMatIndex2));
	}

		int instanceCount = 1;
		int operationCount = (int)theOps.size();

		//accumulate rescaling factors - For scale arrays my indexing scheme and Beagle's happen to be the same, 
		//and my negative (tip) corresponds to a NULL in beagle
		int cumulativeScaleIndex = BEAGLE_OP_NONE;

		CheckBeagleReturnValue(
			beagleUpdatePartials(
				BeagleInstance(),
				(BeagleOperation*)&tupleList[0],
				operationCount,
				cumulativeScaleIndex),
			"beagleUpdatePartials");

#ifdef OUTPUT_PARTIALS
		string name2 = ofprefix + ".partials.log";
		ofstream part(name2.c_str());
		part.precision(5);

		vector<double> outPartials(nstates * subsetData->NChar() * nrates);
		vector<double> outScalers(subsetData->NChar());
		for (list<ClaOperation>::const_iterator it = theOps.begin(); it != theOps.end(); it++) {
			//could pass a scaler index here to have the values unscaled, I think
			beagleGetPartials(beagleInst, PartialIndexForBeagle((*it).destClaIndex), BEAGLE_OP_NONE, &outPartials[0]);
			beagleGetScaleFactors(beagleInst, (*it).destClaIndex, &outScalers[0]);
			for (int site = 0; site < subsetData->NChar(); site++) {
				part << site << "\t";
				for (int rate = 0; rate < nrates; rate++) {
					for (int state = 0; state < 4; state++)
						part << outPartials[site * 4 + state] << "\t";
				}
				part << outScalers[site];
				part << "\n";
			}
		}
		part.close();
#endif

		//this is a bit sneaky. each partial now has a rescaling array associated with it, holding
		//the amount done at that node alone.  Take the rescaling amounts from the two children,
		//and add then to the rescaling already done at this node.  In the next dependency level pass
		//the dest here will become the child there.

		if (IsRescaling()) {
			for (list<ClaOperation>::const_iterator it = theOps.begin(); it != theOps.end(); it++) {
				AccumulateRescalers(
					(*it).destClaIndex,
					(*it).childClaIndex1,
					(*it).childClaIndex2);
			}
		}
		/*
			//to indicate that the partial is now clean, fill the holder that represents it, despite the fact that the memory there is never being used
			//DEBUG - need to figure out second argument here (direction) which I think determined how likely a holder is to be recycled.
			for(list<ClaOperation>::const_iterator it = theOps.begin();it != theOps.end();it++){
				claMan->FillHolder((*it).destClaIndex, 0);
				}
		*/
	}

/*beagle rescaling indeces are in the same order as mine, and they match my partial indeces*/
int SubsetCalculationManager::ScalerIndexForBeagle(int ind) const {
#ifdef PASS_HOLDER_INDECES		
	//this is the code used for the cumulative rescaling index which is reused again and again
//there aren't actually that many clas, so don't pass this in G functions
	if (ind == GARLI_FINAL_SCALER_INDEX)
		return claMan->NumClas();
	else
		return ind;
#else
	//this is the code used for the cumulative rescaling index which is reused again and again
	//there aren't actually that many clas, so don't pass this in G functions
	if (ind == GARLI_FINAL_SCALER_INDEX)
		return claMan->NumClas();
	else if (ind < 0)
		return ind;
	else
		return claMan->GetClaIndexForBeagle(ind);
#endif
}

//passed in a holder dependency index
//return beagle partial index (might actually represent tip)
int SubsetCalculationManager::PartialIndexForBeagle(int ind) const {

	//BEAGLEMERGE DEBUG
	//#define PASS_HOLDER_INDECES

#ifdef PASS_HOLDER_INDECES
	return (ind < 0 ? ((-ind) - 1) : ind + data->NTax());
#else
	return (ind < 0 ? ((-ind) - 1) : claMan->GetClaIndexForBeagle(ind) + ntax);
#endif
}

//As a hack, duplicating these in SubsetCalcMan and CalcMan to get all ported functions working
void SubsetCalculationManager::ReserveDestinationClas(BlockingOperationsSet &set) {
	//outman.DebugMessage("reserving dests of level %d", set.opSetDepLevel);
	for (list<ClaOperation>::const_iterator op = set.claOps.begin(); op != set.claOps.end(); op++)
		claMan->ClaimClaFillIfNecessary((*op).destClaIndex, (*op).opDepLevel);
}


void SubsetCalculationManager::UnreserveDependencyClas(const BlockingOperationsSet &set) {
	//outman.DebugMessage("unreserving deps of level %d", set.opSetDepLevel);
	for (list<ClaOperation>::const_iterator op = set.claOps.begin(); op != set.claOps.end(); op++) {
		//This ONLY removes the temp reservations from the dependencies
		//It does NOT add them to the freeableQueue
		if ((*op).childClaIndex1 > -1)
			claMan->RemoveTempReservation((*op).childClaIndex1);
		if ((*op).childClaIndex2 > -1)
			claMan->RemoveTempReservation((*op).childClaIndex2);
	}
}

void SubsetCalculationManager::QueueDependencyClasForFreeing(const BlockingOperationsSet &set) {
	//outman.DebugMessage("queuing deps level %d for freeing", set.opSetDepLevel);
	for (list<ClaOperation>::const_iterator op = set.claOps.begin(); op != set.claOps.end(); op++) {
		//Adds DEPENDENCY clas to the freeable queue AND removes their temp reservations
		claMan->AddToFreeableQueue((*op).childClaIndex1);
		claMan->AddToFreeableQueue((*op).childClaIndex2);
	}
}
/*
void SubsetCalculationManager::QueueDestinationClasForFreeing(const BlockingOperationsSet &set) {
	//outman.DebugMessage("queuing dests level %d for freeing", set.opSetDepLevel);
	int dontFree = scoreOps.begin()->childClaIndex1;
	int dontFree2 = scoreOps.begin()->childClaIndex2;
	for (list<ClaOperation>::const_iterator op = set.claOps.begin(); op != set.claOps.end(); op++) {
		int thisIndex = (*op).destClaIndex;
		//Adds destination clas to the freeable queue AND removes their temp reservations
		//this is annoying, but can't do with dests needed for final scoring.  They will
		//eventually get freed in ResetDepLevelsAndRes
		if (thisIndex != dontFree && thisIndex != dontFree2)
			claMan->AddToFreeableQueue(thisIndex);
	}
}
*/

//print actual beagle resources used by specific instance
void SubsetCalculationManager::ParseInstanceDetails(const BeagleInstanceDetails *det) {
	outman.UserMessageNoCR("Instance details:\n");
	outman.UserMessageNoCR("\t\tnumber: %d\n", det->resourceNumber);
	outman.UserMessageNoCR("\t\tresource name: %s\n", det->resourceName);
	outman.UserMessageNoCR("\t\timplementation name: %s\n", det->implName);
	outman.UserMessageNoCR("\t\tFlags: ");

	beagle_instance_flags = det->flags;
	//BPART todo - figure out what to do with these at SubMan level
	//string flag_string;
	//InterpretBeagleResourceFlags(beagle_instance_flags, flag_string);
	//outman.UserMessage("%s", flag_string.c_str());
}

//this does the last operation, beagleCalculateEdgeLogLikelihoods, and expects that the 
//required partials have been calculated as well as the transmats necessary for this operation
ScoreSet SubsetCalculationManager::PerformScoringOperation(const ScoringOperation *theOp) {
	ScoreSet results = { 0.0, 0.0, 0.0 };
#ifdef OUTPUT_OTHER_BEAGLE
	//	outman.DebugMessage("ENTER PERF SCR");
#endif
#ifdef USE_BEAGLE
	//state freqs - these are always being sent and stored in the same slot
	vector<FLOAT_TYPE> freqs(nstates);
	//BEAGLE MERGE
	pmatMan->GetModelForTransMatHolder(theOp->transMatIndex, modelIndex)->GetStateFreqs(&(freqs[0]));
	//pmatMan->GetModelForTransMatHolder(theOp->transMatIndex)->stateFreqs;

#ifdef OUTPUT_OTHER_BEAGLE
	outman.DebugMessageNoCR("freqs sent\t");
	for(int i = 0; i < 4;i++)
		outman.DebugMessageNoCR("%.6f\t", freqs[i]);
	outman.DebugMessage("");
	
#endif
	int freqIndex = 0;
	CheckBeagleReturnValue(
		beagleSetStateFrequencies(
			beagleInst,
			freqIndex,
			&(freqs[0])),
		"beagleSetStateFrequencies");

	//category weights (or probs) - these are also always being sent.  annoying to figure out when they need to be updated
	vector<FLOAT_TYPE> inWeights;
	pmatMan->GetCategoryWeightsForBeagle(theOp->transMatIndex, inWeights, modelIndex);
#ifdef OUTPUT_OTHER_BEAGLE
	if (inWeights.size() > 1) {
		outman.DebugMessageNoCR("rate props sent\t");
		for (int i = 0; i < inWeights.size(); i++)
			outman.DebugMessageNoCR("%.6f\t", inWeights[i]);
		outman.DebugMessage("");
	}
#endif
	int weightIndex = 0;
	if (inWeights.size() > 0) {
		CheckBeagleReturnValue(
			beagleSetCategoryWeights(
				beagleInst,
				weightIndex,
				&(inWeights[0])),
			"beagleSetCategoryWeights");
	}

	//the two children to be included- they should have been calculated and reserved
	//earlier in UpdateClas
	int buffer1[1] = { PartialIndexForBeagle(theOp->childClaIndex1) };
	int buffer2[1] = { PartialIndexForBeagle(theOp->childClaIndex2) };
	int numInput = 2;

	//arrays to hold site lnLs and derivs returned from Beagle
	//note that later beagle versions don't return the array anymore, just the sum
	vector<double> siteLikesOut(nchar);
	vector<double> siteD1Out(nchar);
	vector<double> siteD2Out(nchar);

	int pmatIndeces[1] = { PmatIndexForBeagle(theOp->transMatIndex) };
	int	d1MatIndeces[1] = { D1MatIndexForBeagle(theOp->transMatIndex) };
	int	d2MatIndeces[1] = { D2MatIndexForBeagle(theOp->transMatIndex) };

#ifdef OUTPUT_OTHER_BEAGLE
	outman.DebugMessageNoCR("Scr:\t");
	outman.DebugMessageNoCR("\tC\t%d (%d)\t%d (%d)", buffer2[0], theOp->childClaIndex2, buffer1[0], theOp->childClaIndex1);
	outman.DebugMessage("\tP\t%d (%d)", pmatIndeces[0], theOp->transMatIndex);
#endif

	int count = 1;

	//accumulate rescaling factors of the two clas that are being combined (one might be a tip)
	//For scale arrays my indexing scheme and Beagle's happen to be the same, and my negative (tip) corresponds to a NULL in beagle
	//use this scratch index to hold the final accumulated scalers 
	int cumulativeBeagleScaleIndex = BEAGLE_OP_NONE;
	if (IsRescaling()) {
		cumulativeBeagleScaleIndex = ScalerIndexForBeagle(GARLI_FINAL_SCALER_INDEX);
		beagleResetScaleFactors(beagleInst, cumulativeBeagleScaleIndex);
		//add all of the rescaling up to the two involved children.  This will be passed to
		//calcEdgeLike, which will give the final scores
		AccumulateRescalers(GARLI_FINAL_SCALER_INDEX, theOp->childClaIndex1, theOp->childClaIndex2);
	}

	CheckBeagleReturnValue(
		beagleCalculateEdgeLogLikelihoods(
			beagleInst,
			buffer2,
			buffer1,
			pmatIndeces,
			(theOp->derivatives ? d1MatIndeces : NULL),
			(theOp->derivatives ? d2MatIndeces : NULL),
			&weightIndex,
			&freqIndex,
			&cumulativeBeagleScaleIndex,
			count,
			&siteLikesOut[0],
			(theOp->derivatives ? &siteD1Out[0] : NULL),
			(theOp->derivatives ? &siteD2Out[0] : NULL)),
		"beagleCalculateEdgeLogLikelihoods");

	bool beagleReturnsSums = true;

	if (beagleReturnsSums) {

		//The relative rate of this subset is needed to rescale the derivatives
		double subsetRate = pmatMan->GetSubsetRate(theOp->transMatIndex, modelIndex);

		results.lnL = siteLikesOut[0];
		results.d1 = siteD1Out[0] * subsetRate;
		results.d2 = siteD2Out[0] * subsetRate * subsetRate;

#ifdef OUTPUT_BEAGLE_SITELIKES
		string name = ofprefix + ".Bsitelikes.log";
		ofstream file(name.c_str(), ios::app);
		file.precision(10);
		OutputBeagleSiteValues(file, theOp->derivatives);
		file.close();
		//mySummation = SumSiteValues(&siteLikesOut[0], (theOp->derivatives ? &siteD1Out[0] : NULL), (theOp->derivatives ? &siteD2Out[0] : NULL));
		//outman.UserMessage("mine = %.4f beag = %.4f", mySummation.lnL, results.lnL);
#endif
	}
	//else
	//	results = SumSiteValues(&siteLikesOut[0], (theOp->derivatives ? &siteD1Out[0] : NULL), (theOp->derivatives ? &siteD2Out[0] : NULL));
	assert(results.lnL < 0.0 && results.lnL > -10.0e10);
	assert(results.d1 < 10.0e15 && results.d1 > -10.0e15);
	assert(results.d2 < 10.0e25 && results.d2 > -10.0e25);
#ifdef OUTPUT_OTHER_BEAGLE
	outman.DebugMessage("L\t%f\tD1\t%f\tD2\t%f", results.lnL, results.d1, results.d2);
#endif
#ifdef TERMINATE_AFTER_SCORE
	throw ErrorException("TERMINATING AFTER ONE SCORE");
#endif

#endif
	return results;
}

void SubsetCalculationManager::GetBeagleSiteLikelihoods(double *likes) {
	beagleGetSiteLogLikelihoods(beagleInst, &likes[0]);
}
