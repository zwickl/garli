// GARLI version 0.96b8 source code
// Copyright 2005-2008 Derrick J. Zwickl
// email zwickl@nescent.org
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
//
//	NOTE: Portions of this source adapted from GAML source, written by Paul O. Lewis
#include "adaptation.h"
#include "population.h"

#include "garlireader.h"
extern FLOAT_TYPE globalBest; // defined in population.cpp
extern int optCalcs ; // defined in population.cpp

bool ShouldWriteResults(bool prematureTermination, Population::output_details od, unsigned currRep, unsigned nReps);
unsigned RefillTreeBuffer(GarliReader &reader, unsigned treeNum);

unsigned RefillTreeBuffer(GarliReader &, unsigned treeNum) {
	return treeNum;
}
	
bool ShouldWriteResults(bool prematureTermination, Population::output_details od, unsigned currRep, unsigned nReps) {
	if (prematureTermination) 
		return (od & Population::WRITE_PREMATURE);
	if (od & Population::WRITE_REP_TERM)
		return true;
	return (currRep == nReps) && (od & Population::WRITE_REPSET_TERM);
	}

#if defined(SUBROUTINE_GARLI)

void Population::AddTaxonRunMode() {
	throw ErrorException("MPI support for the AddTaxonRunMode has not been added");
}

#elif defined (SWAP_BASED_TERMINATION)

void Population::AddTaxonRunMode() {
	throw ErrorException("SWAP_BASED_TERMINATION support for the AddTaxonRunMode has not been added");
}

#elif defined (BOINC) 

void Population::AddTaxonRunMode() {
	throw ErrorException("BOINC support for the AddTaxonRunMode has not been added");
}

#else
void Population::AddTaxonRunMode() {
	if (conf->restart)
		throw ErrorException("Checkpointing is not supported in AddTaxonRunMode");
	if (conf->searchReps > 1)
		throw ErrorException("The multiple searchReps option is not supported in AddTaxonRunMode");
	const GeneralGamlConfig::StartingTree startMode = conf->GetStartMode();
	if (startMode != GeneralGamlConfig::INCOMPLETE_FROM_FILE_START)
		throw ErrorException("streefname must be set to incomplete when using AddTaxonRunMode");
	if (!startingTreeInNCL)
		throw ErrorException("Incomplete starting trees are only supported when using NEXUS starting trees");
	currentSearchRep = 1;
	string s;
	unsigned attachmentsPerTaxonVar = conf->attachmentsPerTaxon;
	FLOAT_TYPE branchOptPrecisionVar = adap->branchOptPrecision;

	//ensure that the user can ctrl-c kill the program during creation of each stepwise addition tree
	//the signal handling will be returned to the custom message below
	signal( SIGINT, SIG_DFL );

	GetRepNums(s);
	if(s.length() > 0)
		outman.UserMessage("\n>>>%s<<<", s.c_str());
	
	
	GarliReader & reader = GarliReader::GetInstance();
	Individual scratchIndividual;// used in INCOMPLETE_FROM_FILE_START mode
	
	
	for(unsigned i=0;i<total_size;i++){
		if(indiv[i].treeStruct != NULL)
			indiv[i].treeStruct->RemoveTreeFromAllClas();
		if(newindiv[i].treeStruct != NULL)
			newindiv[i].treeStruct->RemoveTreeFromAllClas();
		}

	InitializeOutputStreams();

	//create the first indiv, and then copy the tree and clas
	indiv[0].mod->SetDefaultModelParameters(data);

	unsigned treeNum = 0;
	const NxsFullTreeDescription * treeDescription = 0L;
	unsigned numTrees = reader.GetNumTrees();
	if(reader.FoundModelString())
		startingModelInNCL = true;
	string modelString;

	unsigned totalNumTrees = numTrees;
	for (;!prematureTermination;) {
		if (treeNum >= numTrees) {
			treeNum = RefillTreeBuffer(reader, treeNum);
			std::string newModelString = reader.GetModelString();
			if (newModelString != modelString) {
				modelString = newModelString;
				indiv[0].mod->ReadGarliFormattedModelString(modelString);
				outman.UserMessage("Model parameter values set.");
				outman.UserMessage("MODEL REPORT - Parameters are at their INITIAL values (not yet optimized)");
				indiv[0].mod->OutputHumanReadableModelReportWithParams();
				}
				
			numTrees = reader.GetNumTrees();
			totalNumTrees += numTrees;
		}
		
		treeDescription = reader.GetNxsFullTreeDescription(treeNum);
		if (treeDescription == 0L) {
			outman.UserMessage("No more trees to evaluate. Exiting...");
			break;
			}
		outman.UserMessage("Obtained incomplete starting tree %d from Nexus", treeNum+1);
		
		this->NextAddTaxonRound(*treeDescription, attachmentsPerTaxonVar, branchOptPrecisionVar, scratchIndividual, numTrees, totalNumTrees);
		//this needs to be set here so that the population is reset at the top of this loop before the next rep
		conf->restart = false;
		++treeNum;
		}
	//finalize anything that needs it at rep end
	FinalizeOutputStreams(0);
	//finalize anything that needs it at the end of the repset
	if(currentSearchRep == conf->searchReps)
		FinalizeOutputStreams(1);
	ClearStoredTrees();
}
#endif //defined(SUBROUTINE_GARLI)


// should only be called by AddTaxonRunMode
void Population::NextAddTaxonRound(const NxsFullTreeDescription & treeDescription, 
								   unsigned attachmentsPerTaxonVar,
								   FLOAT_TYPE branchOptPrecisionVar,
								   Individual & scratchIndividual,
								   unsigned repN,
								   unsigned nReps)
	{
	stopwatch.Restart();
	ResetTerminationVariables();
	string s;
	const unsigned nTax = data->NTax();
	scratchIndividual.ReadNxsFullTreeDescription(treeDescription, false);
	scratchIndividual.SetDirty();
	indiv[0].SetDirty();

	outman.UserMessage("Starting with seed=%d\n", rnd.seed());
	if(Tree::constraints.empty())
		outman.UserMessage("using stepwise addition to complete the starting tree...");
	else 
		outman.UserMessage("using stepwise addition to complete the starting tree (compatible with constraints)...");

	globalBest = ZERO_POINT_ZERO;
	indiv[0].FinishIncompleteTreeByStepwiseAddition(nTax, attachmentsPerTaxonVar, branchOptPrecisionVar, scratchIndividual);
	indiv[0].SetTopo(0);

	const char * msg = 0L;
	if(modSpec.IsNucleotide() && modSpec.IsUserSpecifiedStateFrequencies() && !modSpec.gotStateFreqsFromFile)
		msg  = "state frequencies specified as fixed, but no\n\tparameter values found in %s or %s!";
	else if(modSpec.fixAlpha && !modSpec.gotAlphaFromFile)
		msg = "alpha parameter specified as fixed, but no\n\tparameter values found in %s or %s!";
	else if(modSpec.fixInvariantSites && !modSpec.gotPinvFromFile) 
		msg = "proportion of invariant sites specified as fixed, but no\n\tparameter values found in %s or %s!";
	else if(modSpec.IsUserSpecifiedRateMatrix() && !modSpec.gotRmatFromFile) 
		msg = "relative rate matrix specified as fixed, but no\n\tparameter values found in %s or %s!";
	if (msg)
		throw ErrorException(msg, conf->GetTreeFilename(), conf->datafname.c_str());

	assert(indiv[0].treeStruct != NULL);
	bool foundPolytomies = indiv[0].treeStruct->ArbitrarilyBifurcate();
	if(foundPolytomies) outman.UserMessage("WARNING: Polytomies found in start tree.  These were arbitrarily resolved.");

	indiv[0].treeStruct->root->CheckTreeFormation();
	indiv[0].treeStruct->root->CheckforPolytomies();

	indiv[0].treeStruct->CheckBalance();
	indiv[0].treeStruct->SetModel(indiv[0].mod);
	indiv[0].CalcFitness(0);

	FLOAT_TYPE ePrec = 0.0;
	if (repN == 0)
		ePrec = Population::CheckPrecision();
	outman.UserMessage("expected likelihood precision = %.4e", ePrec);

	//if there are not mutable params in the model, remove any weight assigned to the model
	if(indiv[0].mod->NumMutatableParams() == 0) {
		if((conf->bootstrapReps == 0 && currentSearchRep == 1) || (currentBootstrapRep == 1 && currentSearchRep == 1))
			outman.UserMessage("NOTE: Model contains no mutable parameters!\nSetting model mutation weight to zero.\n");
		adap->modelMutateProb=ZERO_POINT_ZERO;
		adap->UpdateProbs();
		}

	outman.precision(10);
	outman.UserMessage("Initial ln Likelihood: %.4f", indiv[0].Fitness());

#	ifdef SCORE_INITIAL_ONLY
		exit(0);
#	endif

#	ifdef MAC_FRONTEND
		NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
		[[MFEInterfaceClient sharedClient] didBeginInitializingSearch];
		[pool release];
#	endif

	if(conf->refineStart==true){
		//12/26/07 now only passing the first argument here ("optModel") as false if no model muts are used
		//if single parameters are fixed that will be checked in the Refine function itself
		indiv[0].RefineStartingConditions(adap->modWeight != ZERO_POINT_ZERO, adap->branchOptPrecision);
		indiv[0].CalcFitness(0);
		outman.UserMessage("lnL after optimization: %.4f", indiv[0].Fitness());
		}

	globalBest = bestFitness = prevBestFitness = indiv[0].Fitness();

#	ifndef INPUT_RECOMBINATION
		for(unsigned i=1;i<total_size;i++){
			if(indiv[i].treeStruct==NULL)
				indiv[i].treeStruct=new Tree();
			indiv[i].CopySecByRearrangingNodesOfFirst(indiv[i].treeStruct, &indiv[0]);
			indiv[i].treeStruct->SetModel(indiv[i].mod);
			}
#	else
		for(unsigned i=1;i<conf->nindivs;i++){
			if(indiv[i].treeStruct==NULL)
				indiv[i].treeStruct=new Tree();
			indiv[i].CopySecByRearrangingNodesOfFirst(indiv[i].treeStruct, &indiv[0]);
			indiv[i].treeStruct->SetModel(indiv[i].mod);
			}
		for(unsigned i=conf->nindivs;i<total_size;i++){
			indiv[i].GetStartingConditionsFromFile(conf->streefname.c_str(), i-conf->nindivs, nTax);
			indiv[i].treeStruct->SetModel(indiv[i].mod);
			indiv[i].SetDirty();
			indiv[i].CalcFitness(0);
			}
#	endif

	UpdateTopologyList(indiv);
	CalcAverageFitness();

	//3/24/08 moving this after SeedPop, since it disallows normal ctrl-c killing of runs during stepwise
	CatchInterrupt();
	RunImplForAddTaxonRunMode();

	//this rep is over
	if (prematureTermination == false) {
		if (s.length() > 0)
			outman.UserMessage(">>>Completed %s<<<\n", s.c_str());
		//not sure where this should best go
		outman.UserMessage("MODEL REPORT - Parameter values are FINAL");
		indiv[bestIndiv].mod->OutputHumanReadableModelReportWithParams();
		if(Tree::outgroup != NULL)
			OutgroupRoot(&indiv[bestIndiv], bestIndiv);
		Individual *repResult = new Individual(&indiv[bestIndiv]);
		if (conf->collapseBranches) {
			int numCollapsed = 0;
			repResult->treeStruct->root->CollapseMinLengthBranches(numCollapsed);
			outman.UserMessage("\nNOTE: Collapsing of minimum length branches was requested (collapsebranches = 1)");\
			if(numCollapsed == 0)
				outman.UserMessage("    No branches were short enough to be collapsed.");
			else
				outman.UserMessage("    %d branches were collapsed.", numCollapsed);
			}
		storedTrees.push_back(repResult);
		}
	else{
		if(s.length() > 0)
			outman.UserMessage(">>>Terminated %s<<<\n", s.c_str());
		outman.UserMessage("NOTE: ***Run was terminated before termination condition was reached!\nLikelihood scores, topologies and model estimates obtained may not\nbe fully optimal!***");
		}

	int best=0;
	if(storedTrees.size() > 1) {
		best=EvaluateStoredTrees(true);
		}

	if (ShouldWriteResults(prematureTermination, all_best_output, currentSearchRep, nReps)) {
		if(storedTrees.size() > 0)
			WriteStoredTrees(besttreefile.c_str());
		}
	if (ShouldWriteResults(prematureTermination, best_output, currentSearchRep, nReps))
		WriteTreeFile(besttreefile.c_str());

	if(conf->bootstrapReps > 0){
		if (ShouldWriteResults(prematureTermination, bootlog_output, currentSearchRep, nReps)) {
			Individual * bt;
			if(conf->searchReps > 1 && storedTrees.size() > 0){
				outman.UserMessage("Saving best search rep (#%d) to bootstrap file", best+1);
				bt = storedTrees[best];
				}
			//this was a bug - when collapse was on and bootstrapping was being done with one search
			//rep, the best tree was being written to the boot file.  The tree in the storedTrees is
			//the one that was actually collapsed though
			//else FinishBootstrapRep(&indiv[bestIndiv], currentBootstrapRep);
			else
				bt = (storedTrees.size() == 1 ? storedTrees[0] : &indiv[bestIndiv]);
			FinishBootstrapRep(bt, currentBootstrapRep);
			}
		else if (prematureTermination && !(bootlog_output & WRITE_PREMATURE))
			outman.UserMessage("Not saving search rep to bootstrap file due to early termination");
		}

	if(conf->inferInternalStateProbs == true){
		if(prematureTermination == false && currentSearchRep == conf->searchReps){
			if(storedTrees.size() > 0){//careful here, the trees in the storedTrees array don't have clas assigned
				outman.UserMessage("Inferring internal state probabilities on best tree....");
				storedTrees[best]->treeStruct->InferAllInternalStateProbs(conf->ofprefix.c_str());
				}
			}
		else if(prematureTermination){
			outman.UserMessage(">>>Internal state probabilities not inferred due to premature termination<<<");
			}
		}
	if(prematureTermination == true)
		return;

	//write a checkpoint that will indicate that the rep is done and results have been written to file
	//the gen will be UINT_MAX, as it is after a rep has terminated, which will tell the function that reads
	//the checkpoint to set finishedrep = true.  This automatically happens in the BOINC case
#	ifndef BOINC
		if(conf->checkpoint)
			WriteStateFiles();
#	else
		WriteStateFiles();
#	endif
	}


void Population::RunImplForAddTaxonRunMode() {
	optCalcs=0;

#	ifdef MAC_FRONTEND
		NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
		[[MFEInterfaceClient sharedClient] didBeginRun];
		[pool release];
#	endif

	outman.precision(6);
	outman.UserMessage("%-10s%-15s%-10s%-15s", "gen", "current_lnL", "precision", "last_tree_imp");
	outman.UserMessage("%-10d%-15.4f%-10.3f\t%-15d", gen, BestFitness(), adap->branchOptPrecision, lastTopoImprove);
	OutputLog();
	if(conf->outputMostlyUselessFiles) OutputFate();

	CatchInterrupt();

	gen++;
	for (; gen < conf->stopgen+1; ++gen){

		NextGeneration();

		keepTrack();
		if(conf->outputMostlyUselessFiles)
			OutputFate();
		if(conf->logevery > 0 && !(gen % conf->logevery))
			OutputLog();
		if(conf->saveevery > 0 && !(gen % conf->saveevery)) 
			OutputSave();
			

		prematureTermination = CheckForUserSignal();
		if(prematureTermination)
			break;

#		ifdef PERIODIC_SCORE_DEBUG
			if(gen % 500 == 0 ||gen==1)
				OutputFilesForScoreDebugging(&indiv[bestIndiv], tempGlobal++);
#		endif

#		ifdef NNI_SPECTRUM
			if(gen % 1000 == 0 || gen==1)
				NNISpectrum(bestIndiv);
#		endif

		if(!(gen%adap->intervalLength)) {
			const ContinuationCode x = this->AdaptPrec();
			if (x == FINISH_RUN)
				break;
			if (x == ABORT_IMMEDIATELY)
				return;
		}
		if(conf->checkpoint==true && ((gen % conf->saveevery) == 0))
			WriteStateFiles();

		if(stopwatch.SplitTime() > (int)conf->stoptime){
			outman.UserMessage("NOTE: ****Specified time limit (%d seconds) reached...", conf->stoptime);
			prematureTermination = true;
			break;
			}
		if(gen == conf->stopgen)
			outman.UserMessage("NOTE: ****Specified generation limit (%d) reached...", conf->stopgen);
#		ifdef INCLUDE_PERTURBATION
			if(pertMan->pertAbandoned==true && pertMan->restartAfterAbandon==true && (gen - pertMan->lastPertGeneration > pertMan->gensBeforeRestart)){
				params->starting_tree="";
				pertMan->lastPertGeneration=gen;
				pertMan->pertAbandoned=false;
				pertMan->numPertsNoImprove=0;
				bestSinceRestart.SetFitness(-1e100);
				if(BestFitness() > allTimeBest.Fitness()) StoreAllTimeBest();
				SeedPopulationWithStartingTree();
				outman.UserMessage("restarting ....");
				}
#		endif
		}

	FinalOptimization();
	gen = UINT_MAX;
	OutputLog();

	if(conf->bootstrapReps==0)
		outman.UserMessage("finished");
}
