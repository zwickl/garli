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
#include <csignal>
#include <sstream>
#include <ctime>
#include "adaptation.h"
#include "population.h"

#include "garlireader.h"
extern FLOAT_TYPE globalBest; // defined in population.cpp
extern int optCalcs ; // defined in population.cpp


bool ShouldWriteResults(bool prematureTermination, Population::output_details od, unsigned currRep, unsigned nReps);
void treeToWrapper(char * treeString, Individual & ind);
void writeGarliIndividualDescription(char * treeString, Individual & ind);


enum InteractiveModeActionEnum {FINISH_STEPADD_ACTION, SEARCH_ACTION, SWAP_ACTION};

////////////////////////////////////////////////////////////////////////////////
// globals hacked in for the AddTaxonRunMode version:
InteractiveModeActionEnum gInteractive_mode_action =  FINISH_STEPADD_ACTION;
MultiFormatReader * gConstraintReader = 0L;
bool gModelIsFixed;
unsigned gCurrIGarliResultIndex = 0;
bool gGetTreesFromReader = true;
std::string gModelString;
typedef std::pair<double, std::string> ScoreStringPair;
typedef std::vector< ScoreStringPair > ScoreStringPairList;
ScoreStringPairList gSuboptimalTreeList;
char * gTreeBufferString = 0L;
bool gOptimizePlausibleDuringStepwise = true;
bool gDoOutputSiteLike = false;
bool gIsFinalScoring = false;
unsigned gMaxSuboptimalTreesToStore = 0;
std::vector<unsigned> gPrevCharWt;
// end globals hacked in for the AddTaxonRunMode version:
////////////////////////////////////////////////////////////////////////////////

MultiFormatReader * createConstraintFileReader() {
	MultiFormatReader * mfr = new MultiFormatReader();
	return mfr;
}

void suboptimalTreesToWrapper();
void ResetSuboptimalTreesList() {
	gSuboptimalTreeList.clear();
	gSuboptimalTreeList.reserve(gMaxSuboptimalTreesToStore + 1);
}

bool Population::GetOutputSiteLikes() {
	return gDoOutputSiteLike && gIsFinalScoring;
}

void Population::RecordTreeFoundDuringSearch(const Tree &tree) {
	assert(gTreeBufferString);
	if (gMaxSuboptimalTreesToStore == 0)
		return;
	ScoreStringPair ssp(-tree.lnL, std::string() );
	const ScoreStringPairList::iterator gstlEnd =  gSuboptimalTreeList.end();
	ScoreStringPairList::iterator lb = std::lower_bound(gSuboptimalTreeList.begin(), gstlEnd, ssp);
	if (lb == gstlEnd) {
		if (gSuboptimalTreeList.size() < gMaxSuboptimalTreesToStore) {
			tree.root->MakeNewick(gTreeBufferString, false, true);
			ssp.second.assign(gTreeBufferString);
			gSuboptimalTreeList.push_back(ssp);
		}
	}
	else {
		tree.root->MakeNewick(gTreeBufferString, false, true);
		ssp.second.assign(gTreeBufferString);
		gSuboptimalTreeList.insert(lb, ssp);
		if (gSuboptimalTreeList.size() > gMaxSuboptimalTreesToStore)
			gSuboptimalTreeList.pop_back();
	}
}

NxsFullTreeDescription readConstraintTreeDesc(const char * nexusContent, NxsTaxaBlock *tb) {
	if (gConstraintReader == 0L)
		gConstraintReader = createConstraintFileReader();
	else {
		assert(gConstraintReader->GetNumTaxaBlocks() == 1);
		assert(gConstraintReader->GetTaxaBlock(N) == 1);
	}
	gConstraintReader->ClearContent();
	if (tb)
		gConstraintReader->AddReadTaxaBlock(tb);

	gConstraintReader->ReadStringAsNexusContent(nexusContent);
	unsigned nt = gConstraintReader->GetNumTreesBlocks(tb);
	if (nt == 0)
		throw ErrorException("No trees block found!");
	if (nt > 1)
		throw ErrorException("More than one trees block found");
	NxsTreesBlock * treesB = gConstraintReader->GetTreesBlock(tb, 0);
	const unsigned ntrees = treesB->GetNumTrees();
	if (ntrees > 1)
		throw ErrorException("More than one trees block found");
	NxsFullTreeDescription t = treesB->GetFullTreeDescription(0);
	delete treesB;
	gConstraintReader->ClearContent();
	return t;
}

void Population::PatternCountsToWrapper()
{
	if (!data) {
		std::cerr << "No data in memory\n";
		return;
	}
	const unsigned nc = data->NChar();
	const int * countit = data->GetCounts();
	std::cerr << "[igarlipatterncounts "; 
	for (unsigned i = 0; i < nc; ++i)
		std::cerr << countit[i] << ' ';
	std::cerr << "]\n";
}
std::pair<unsigned, unsigned> Population::RefillTreeBuffer(GarliReader &reader, unsigned treeNum) {
	unsigned endTreeNum = UINT_MAX;
	int nSleeps = 0;
	std::string nextLine;
	for (;;) {
		GeneralGamlConfig c(*(this->conf));
		std::cerr << "iGarli[0 - " << gCurrIGarliResultIndex  << "]>" << std::endl;
		nextLine.clear();
		std::getline(std::cin, nextLine);
		if (nextLine.empty()) {
			std::cerr << "getline failed\n";
			if (nSleeps++ > 100) {
				sleep(1);
				if (nSleeps == 1) {
					//finalize anything that needs it at rep end
					FinalizeOutputStreams(0);
					//finalize anything that needs it at the end of the repset
					if(currentSearchRep == conf->searchReps)
						FinalizeOutputStreams(1);
				}
				continue;
			}
		}
		else {
			if (nSleeps == 1) {
				InitializeOutputStreams();
			}
			nSleeps = 0;
		}
		try {
			if (NxsString::case_insensitive_equals(nextLine.c_str(), "run") || NxsString::case_insensitive_equals(nextLine.c_str(), "r")) {
				if (gGetTreesFromReader) {
					if (reader.GetNxsFullTreeDescription(treeNum) != 0L)
						break;
				}
				else {
					if (treeNum < storedTrees.size())
						break;
				}
				throw ErrorException("no operation queued -- the run command is inappropriate at this point.\n Use \"quit\" to exit");
			}
			else if (NxsString::case_insensitive_equals(nextLine.c_str(), "quit")
					|| NxsString::case_insensitive_equals(nextLine.c_str(), "q")
					|| NxsString::case_insensitive_equals(nextLine.c_str(), "exit")) {
				if (gGetTreesFromReader)
					return std::pair<unsigned, unsigned>(reader.GetNumTrees(), 0);
				return std::pair<unsigned, unsigned>(this->storedTrees.size(), 0);
			}
			const AttemptedParseResult parseResult = c.ParseLineIntoConfigObject(nextLine);
			if (parseResult.first) {
				*(this->conf) = c;
			}
			else {
				const char * key = parseResult.second.first.c_str();
				const std::string & valueStr = parseResult.second.second;
				const char * value = valueStr.c_str();
				if (NxsString::case_insensitive_equals(key, "clear")) {
					NxsTaxaBlock * currTaxa = reader.GetTaxaBlock(0);
					assert (currTaxa);
					if (currTaxa) {
						reader.RemoveBlockFromUsedBlockList(currTaxa);
						reader.DeleteBlocksFromFactories();
						reader.AddReadTaxaBlock(currTaxa);
						assert(reader.GetTaxaBlock(0));
					}
					treeNum = 0;
				}
				else if (NxsString::case_insensitive_equals(nextLine.c_str(), "run") || NxsString::case_insensitive_equals(nextLine.c_str(), "r")) {
					if (gGetTreesFromReader) {
						if (reader.GetNxsFullTreeDescription(treeNum) != 0L)
							break;
					}
					else {
						if (treeNum < storedTrees.size())
							break;
					}
					throw ErrorException("no operation queued -- the run command is inappropriate at this point.\n Use \"quit\" to exit");
				}
				else if (NxsString::case_insensitive_equals(nextLine.c_str(), "quit")
						|| NxsString::case_insensitive_equals(nextLine.c_str(), "q")
						|| NxsString::case_insensitive_equals(nextLine.c_str(), "exit")) {
					if (gGetTreesFromReader)
						return std::pair<unsigned, unsigned>(reader.GetNumTrees(), 0);
					return std::pair<unsigned, unsigned>(this->storedTrees.size(), 0);
				}
				else if (NxsString::case_insensitive_equals(key, "model")) {
					gModelString.assign(value);
					std::cerr << "model is now = " << gModelString << std::endl;
				}
				else if (NxsString::case_insensitive_equals(key, "tree")) {
					const NxsTaxaBlock * queryTaxa = reader.GetTaxaBlock(0);
					if (queryTaxa) {
						const std::string & taxaTitle = queryTaxa->GetTitle();
						std::ostringstream outF;
						outF << "#NEXUS\nBegin Trees;\n Link taxa = " << NxsString::GetEscaped(taxaTitle);
						outF << ";\n Tree stdintree = " << value << " ;\nend;\n";
						reader.ReadStringAsNexusContent(outF.str());
					}
					else {
						assert (queryTaxa);
						throw ErrorException("Tree command is not available when a Taxa block has not been read.");
					}
				}
				else if (NxsString::case_insensitive_equals(key, "stepadd"))
					gInteractive_mode_action = FINISH_STEPADD_ACTION;
				else if (NxsString::case_insensitive_equals(key, "treenum") || NxsString::case_insensitive_equals(key, "endtreenum") ) {
					long tn = 0;
					if (!NxsString::to_long(value, &tn)) {
						throw ErrorException("Expecting an index to follow treenum");
					}
					if (tn < 0 || tn >= reader.GetNumTrees()) {
						throw ErrorException("treenum index is out of range");
					}
					if (NxsString::case_insensitive_equals(key, "treenum"))
						treeNum = (unsigned) tn;
					else
						endTreeNum = (unsigned)tn;
				}
				else if (NxsString::case_insensitive_equals(key, "nbest")) {
					long tn = -1;
					if (!NxsString::to_long(value, &tn) || tn < 0) {
						throw ErrorException("Expecting a non-negative integer to follow \"nbest\"");
					}
					gMaxSuboptimalTreesToStore  = (unsigned) tn;
				}
				else if (NxsString::case_insensitive_equals(key, "optplausible")) {
					long tn = -1;
					if (!NxsString::to_long(value, &tn) || tn < 0) {
						throw ErrorException("Expecting 1 or 0 after \"optplausible\"");
					}
					gOptimizePlausibleDuringStepwise  = bool(tn != 0);
				}
				else if (NxsString::case_insensitive_equals(key, "sitelikes")) {
					long tn = -1;
					if (!NxsString::to_long(value, &tn) || tn < 0) {
						throw ErrorException("Expecting 1 or 0 after \"sitelikes\"");
					}
					gDoOutputSiteLike  = bool(tn != 0);
				}
				else if (NxsString::case_insensitive_equals(key, "patterncounts")) {
					long tn = -1;
					if (!NxsString::to_long(value, &tn) || tn < 0) {
						throw ErrorException("Expecting 1 or 0 after \"patterncounts\"");
					}
					if (bool(tn != 0) )
						PatternCountsToWrapper();
				}
				else if (NxsString::case_insensitive_equals(key, "setpatterncounts")) {
					std::list<std::string> cList;
					NxsString::split(valueStr, &cList);
					std::vector<unsigned> countsRead;
					countsRead.reserve(cList.size());
					for (std::list<std::string>::const_iterator cIt = cList.begin(); cIt != cList.end(); ++cIt) {
						long tn = -1;
						if (!NxsString::to_long(cIt->c_str(), &tn) || tn < 0)
							throw ErrorException("Expecting non-negative integers after setpatterncounts");
						countsRead.push_back((unsigned) tn);
					}
					if (!data)
						throw ErrorException("No data in memory");
					if (data->NChar() != countsRead.size())
						throw ErrorException("Incorrect number of pattern counts specified");
					unsigned patIndex = 0;
					std::vector<unsigned>::const_iterator ucIt = countsRead.begin();
					for (; ucIt != countsRead.end(); ++ucIt, ++patIndex)
						data->SetCount(patIndex, *ucIt);
				}
				else if (NxsString::case_insensitive_equals(key, "setcharintwts")){
					std::list<std::string> cList;
					NxsString::split(valueStr, &cList);
					std::vector<unsigned> countsRead;
					countsRead.reserve(cList.size());
					for (std::list<std::string>::const_iterator cIt = cList.begin(); cIt != cList.end(); ++cIt) {
						long tn = -1;
						if (!NxsString::to_long(cIt->c_str(), &tn) || tn < 0)
							throw ErrorException("Expecting non-negative integers after setcharintwts");
						countsRead.push_back((unsigned) tn);
					}
					if (!data)
						throw ErrorException("No data in memory");
					if (data->GapsIncludedNChar() != countsRead.size())
						throw ErrorException("Incorrect number of character weights specified");
					unsigned charIndex = 0;
					std::vector<unsigned>::const_iterator ucIt = countsRead.begin();
					if (gPrevCharWt.empty())
						gPrevCharWt.assign(data->GapsIncludedNChar(), 1);
					for (; ucIt != countsRead.end(); ++ucIt, ++charIndex) {
						const int newWt = (int)*ucIt;
						const int oldWt = (int)gPrevCharWt[charIndex];
						if (newWt != oldWt) {
							int number = data->Number(charIndex); // great name!
							if (number >= 0) {
								const int oldCount = data->Count(number);
								const int diff = newWt - oldWt ;
								const int newCount = oldCount + diff;
								data->SetCount(number, newCount);
							}
						}						
						gPrevCharWt[charIndex] = newWt;
					}
				}
				else if (NxsString::case_insensitive_equals(key, "search"))
					gInteractive_mode_action = SEARCH_ACTION;
				else if (NxsString::case_insensitive_equals(key, "swap"))
					gInteractive_mode_action = SWAP_ACTION;
				else if (NxsString::case_insensitive_equals(key, "clearconstraints")) {
					Tree::ClearConstraints();
					std::cerr << "constraints cleared\n";
				}
				else if (NxsString::case_insensitive_equals(key, "treesource")) {
					if (NxsString::case_insensitive_equals(value, "results"))
						gGetTreesFromReader = false;
					else if (NxsString::case_insensitive_equals(value, "input"))
						gGetTreesFromReader = true;
					else
						throw ErrorException("The \"TreeSource\" value must be \"Results\" or \"Input\"");
				}
				else if (NxsString::case_insensitive_equals(key, "posconstraint") || NxsString::case_insensitive_equals(key, "negconstraint")) {
					const bool negConstraint = NxsString::case_insensitive_equals(key, "negconstraint");
					NxsTaxaBlock * queryTaxa = reader.GetTaxaBlock(0);
					if (queryTaxa) {
						const std::string & taxaTitle = queryTaxa->GetTitle();
						std::ostringstream outF;
						outF << "#NEXUS\nBegin Trees;\n Link taxa = " << NxsString::GetEscaped(taxaTitle);
						outF << ";\n Tree stdintree = " << value << " ;\nend;\n";
						std::string s = outF.str();
						NxsFullTreeDescription newConTreeDescription = readConstraintTreeDesc(s.c_str(), queryTaxa);
						std::string constraintNewick = newConTreeDescription.GetNewick();
						Tree::ReadNewickConstraint(constraintNewick.c_str(), true, !negConstraint);
					}
					else {
						assert (queryTaxa);
						throw ErrorException("Constraint command is not available when a Taxa block has not been read.");
					}
				}
				else {
					throw ErrorException("Unrecognized command.");
				}
			}
		}
		catch (ErrorException & x) {
			x.Print(std::cerr);
		}
	}
	return std::pair<unsigned, unsigned>(treeNum, endTreeNum);
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

	const FLOAT_TYPE taxsize = log10((FLOAT_TYPE) ((FLOAT_TYPE)data->NTax())*data->NTax()*2);
	const unsigned stringSize=(int)((data->NTax()*2)*(10+DEF_PRECISION)+taxsize);
	if (gTreeBufferString)
		delete gTreeBufferString;
	gTreeBufferString = new char[stringSize];

	GarliReader & reader = GarliReader::GetInstance();
	Individual scratchIndividual;// used in INCOMPLETE_FROM_FILE_START mode


	for(unsigned i = 0; i < total_size; i++){
		if(indiv[i].treeStruct != NULL)
			indiv[i].treeStruct->RemoveTreeFromAllClas();
		if(newindiv[i].treeStruct != NULL)
			newindiv[i].treeStruct->RemoveTreeFromAllClas();
		}

	InitializeOutputStreams();

	//create the first indiv, and then copy the tree and clas
	indiv[0].mod->SetDefaultModelParameters(data);

	unsigned numTrees = reader.GetNumTrees();
	unsigned treeNum = numTrees;
	unsigned endTreeNum = numTrees;
	if(reader.FoundModelString())
		startingModelInNCL = true;
	string modelString;

	unsigned totalNumTrees = numTrees;
	bool filled = false;
	for (;!prematureTermination;) {
		std::cerr << storedTrees.size() << " trees in storedTrees\n";
		std::cerr << numTrees << " trees in the input tree reader\n";

		if (treeNum >= numTrees || treeNum >= endTreeNum) {
			if (!filled)
				treeNum = 0;
			filled = true;
			std::pair<unsigned, unsigned> treeRange = RefillTreeBuffer(reader, treeNum);
			treeNum = treeRange.first;
			endTreeNum = treeRange.second;
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

		stopwatch.Restart();
		ResetTerminationVariables();
		if (gGetTreesFromReader) {
			const NxsFullTreeDescription * treeDescription = reader.GetNxsFullTreeDescription(treeNum);
			if (treeDescription == 0L) {
				outman.UserMessage("No more trees to evaluate. Exiting...");
				break;
			}
			scratchIndividual.ReadNxsFullTreeDescription(*treeDescription, false);
		}
		else {
			if (treeNum >= storedTrees.size() || storedTrees[treeNum] == 0L) {
				outman.UserMessage("No more trees to evaluate. Exiting...");
				break;
			}
			else {
				scratchIndividual = *storedTrees[treeNum];
			}
		}

		this->Reset(); // this triggers the creation of a new adap object -- among other things

		if (!gModelString.empty()) {
			indiv[0].mod->ReadGarliFormattedModelString(gModelString);
			gModelString.clear();
		}
		scratchIndividual.SetDirty();

		ResetSuboptimalTreesList();

		outman.UserMessage("Obtained incomplete starting tree %d from Nexus", treeNum+1);
		if (gInteractive_mode_action == FINISH_STEPADD_ACTION)
			this->NextAddTaxonRound(scratchIndividual, attachmentsPerTaxonVar, branchOptPrecisionVar, numTrees, totalNumTrees);
		else if (gInteractive_mode_action == SEARCH_ACTION) {
			throw ErrorException("Not implemented yet");
		}
		else if (gInteractive_mode_action == SWAP_ACTION) {
			assert(0);
			throw ErrorException("Not implemented yet");
			this->AddTaxonSwap(scratchIndividual, attachmentsPerTaxonVar, branchOptPrecisionVar, numTrees, totalNumTrees);
		}

		AfterRunHook(totalNumTrees);

		if(s.length() > 0)
			outman.UserMessage("\n>>>Completed %s<<<", s.c_str());


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

void Population::AfterRunHook(unsigned nReps) {
	globalBest = bestFitness = prevBestFitness = indiv[0].Fitness();
	//this rep is over
	if (prematureTermination == false) {
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
		outman.UserMessage("NOTE: ***Run was terminated before termination condition was reached!\nLikelihood scores, topologies and model estimates obtained may not\nbe fully optimal!***");
		}

	//if(storedTrees.size() > 1) {
	//	EvaluateStoredTrees(true);
	//	}

	WriteTreeFile(besttreefile.c_str());
	
	if (gDoOutputSiteLike) {
		gIsFinalScoring = true;
		std::cerr << "in gDoOutputSiteLike conditional branch\n";
		assert(this->indiv[bestIndiv].treeStruct);
		Tree & bestTree = *(this->indiv[bestIndiv].treeStruct);
		bestTree.ConditionalLikelihoodRateHet(Tree::ROOT, bestTree.root, false);
//		bestTree.AssignCLAsFromMaster();
//		bestTree.RecursivelyCalculateInternalStateProbs(bestTree.root, 0L); // this triggers the site likelihood output

		gIsFinalScoring = false;
		
	}
	treeToWrapper(this->treeString, this->indiv[bestIndiv]);
	suboptimalTreesToWrapper();
	
	if(prematureTermination == true)
		return;
	if (conf->inferInternalStateProbs == true)
		throw ErrorException("inferInternalStateProbs is not supported in the addtaxonrunmode");
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

#endif //defined(SUBROUTINE_GARLI)


void Population::ReconfigureAdaptationParams() {
	//if there are not mutable params in the model, remove any weight assigned to the model
	if (gModelIsFixed || indiv[0].mod->NumMutatableParams() == 0) {
		if(currentSearchRep == 1 && (conf->bootstrapReps == 0 || currentBootstrapRep == 1)) {
			const char * msg;
			if (gModelIsFixed)
				msg = "NOTE: Model is fixed.\nSetting model mutation weight to zero.\n";
			else
				msg = "NOTE: Model contains no mutable parameters!\nSetting model mutation weight to zero.\n";
			outman.UserMessage(msg);
		}
		adap->modelMutateProb = ZERO_POINT_ZERO;
		adap->UpdateProbs();
	}
	else {
		adap->NormalizeMutateProbs();
	}
}
void Population::AddTaxonSwap(Individual & scratchIndividual,
								   unsigned attachmentsPerTaxonVar,
								   FLOAT_TYPE branchOptPrecisionVar,
								   unsigned repN,
								   unsigned nReps) {
	indiv[0] = scratchIndividual;
	outman.UserMessage("Starting swapping seed=%d\n", rnd.seed());
	globalBest = ZERO_POINT_ZERO;
	bool foundPolytomies = indiv[0].treeStruct->ArbitrarilyBifurcate();
	if(foundPolytomies)
		outman.UserMessage("WARNING: Polytomies found in start tree.  These were arbitrarily resolved.");

	indiv[0].treeStruct->root->CheckTreeFormation();
	indiv[0].treeStruct->CheckBalance();
	indiv[0].treeStruct->SetModel(indiv[0].mod);
	indiv[0].CalcFitness(0);

	FLOAT_TYPE ePrec = 0.0;
	if (repN == 0)
		ePrec = Population::CheckPrecision();
	outman.UserMessage("expected likelihood precision = %.4e", ePrec);

	outman.precision(10);
	outman.UserMessage("Initial ln Likelihood: %.4f", indiv[0].Fitness());


#	ifdef MAC_FRONTEND
		NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
		[[MFEInterfaceClient sharedClient] didBeginInitializingSearch];
		[pool release];
#	endif

	if(conf->refineStart == true){
		//12/26/07 now only passing the first argument here ("optModel") as false if no model muts are used
		//if single parameters are fixed that will be checked in the Refine function itself
		indiv[0].RefineStartingConditions(adap->modWeight != ZERO_POINT_ZERO, adap->branchOptPrecision);
		indiv[0].CalcFitness(0);
		outman.UserMessage("lnL after optimization: %.4f", indiv[0].Fitness());
		}

	ReconfigureAdaptationParams();

	this->SwapToCompletion(SWAPPER_BY_DIST_NOT_FURTHEST_RUN_MODE, branchOptPrecisionVar);
}

void treeToWrapper(char * treeString, Individual & ind) {
	std::cerr << "[iGarli "<< gCurrIGarliResultIndex++ << " ] ";
	writeGarliIndividualDescription(treeString, ind);
}

void suboptimalTreesToWrapper() {
	std::cerr.setf( ios::floatfield, ios::fixed );
	std::cerr.setf( ios::showpoint );
	ScoreStringPairList::const_iterator it = gSuboptimalTreeList.begin();
	unsigned suboptN = 1;
	for (; it != gSuboptimalTreeList.end() ; ++it, ++suboptN) {
		std::cerr << "[iGarli "<< gCurrIGarliResultIndex++;
		std::cerr << " ] tree subopt" << suboptN << " = [&U][!GarliScore " << -(it->first) << "] " << it->second << " ;\n";
	}
	if (gSuboptimalTreeList.empty())
		std::cerr << "iGarli -- no suboptimal trees recorded"<<  std::endl;
}

void writeGarliIndividualDescription(char * treeString, Individual & ind) {
	std::cerr << "tree best = [&U][!GarliScore ";
	if (ind.IsDirty())
		std::cerr << "-0.0" ;
	else
		std::cerr << ind.Fitness() ;
	std::string modstr;
	ind.mod->FillGarliFormattedModelString(modstr);
	std::cerr << "][!GarliModel " <<  modstr <<  "] ";
	std::cerr.setf( ios::floatfield, ios::fixed );
	std::cerr.setf( ios::showpoint );
	ind.treeStruct->root->MakeNewick(treeString, false, true);
	std::cerr << treeString << " ;\n";
}

// should only be called by AddTaxonRunMode
void Population::NextAddTaxonRound(Individual & scratchIndividual,
								   unsigned attachmentsPerTaxonVar,
								   FLOAT_TYPE branchOptPrecisionVar,
								   unsigned repN,
								   unsigned nReps)
	{

	outman.UserMessage("Starting with seed=%d\n", rnd.seed());
	if(Tree::IsUsingConstraints())
		outman.UserMessage("using stepwise addition to complete the starting tree...");
	else
		outman.UserMessage("using stepwise addition to complete the starting tree (compatible with constraints)...");

	globalBest = ZERO_POINT_ZERO;

	std::cerr << "\nStarting individual:\n";
	writeGarliIndividualDescription(this->treeString, scratchIndividual);
	const unsigned nTax = data->NTax();
	
	
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
		adap->modelMutateProb = ZERO_POINT_ZERO;
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

	if(conf->refineStart == true){
		//12/26/07 now only passing the first argument here ("optModel") as false if no model muts are used
		//if single parameters are fixed that will be checked in the Refine function itself
		indiv[0].RefineStartingConditions(adap->modWeight != ZERO_POINT_ZERO, adap->branchOptPrecision);
		indiv[0].CalcFitness(0);
		outman.UserMessage("lnL after optimization: %.4f", indiv[0].Fitness());
		}


#	ifndef INPUT_RECOMBINATION
		for(unsigned i = 1; i < total_size; i++){
			if(indiv[i].treeStruct == NULL)
				indiv[i].treeStruct = new Tree();
			indiv[i].CopySecByRearrangingNodesOfFirst(indiv[i].treeStruct, &indiv[0]);
			indiv[i].treeStruct->SetModel(indiv[i].mod);
			}
#	else
		for(unsigned i = 1; i < conf->nindivs; i++){
			if(indiv[i].treeStruct == NULL)
				indiv[i].treeStruct = new Tree();
			indiv[i].CopySecByRearrangingNodesOfFirst(indiv[i].treeStruct, &indiv[0]);
			indiv[i].treeStruct->SetModel(indiv[i].mod);
			}
		for(unsigned i = conf->nindivs; i < total_size; i++){
			indiv[i].GetStartingConditionsFromFile(conf->streefname.c_str(), i-conf->nindivs, nTax);
			indiv[i].treeStruct->SetModel(indiv[i].mod);
			indiv[i].SetDirty();
			indiv[i].CalcFitness(0);
			}
#	endif

	UpdateTopologyList(indiv);
	CalcAverageFitness();
	ReconfigureAdaptationParams();

	//3/24/08 moving this after SeedPop, since it disallows normal ctrl-c killing of runs during stepwise
	CatchInterrupt();
	RunImplForAddTaxonRunMode();
	}


void Population::RunImplForAddTaxonRunMode() {
	optCalcs = 0;

#	ifdef MAC_FRONTEND
		NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
		[[MFEInterfaceClient sharedClient] didBeginRun];
		[pool release];
#	endif

	outman.precision(6);
	outman.UserMessage("%-10s%-15s%-10s%-15s", "gen", "current_lnL", "precision", "last_tree_imp");
	outman.UserMessage("%-10d%-15.4f%-10.3f\t%-15d", gen, BestFitness(), adap->branchOptPrecision, GetLastTopoImprove());
	OutputLog();
	if(conf->outputMostlyUselessFiles)
		OutputFate();

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
			if(gen % 500 == 0 ||gen == 1)
				OutputFilesForScoreDebugging(&indiv[bestIndiv], tempGlobal++);
#		endif

#		ifdef NNI_SPECTRUM
			if(gen % 1000 == 0 || gen == 1)
				NNISpectrum(bestIndiv);
#		endif

		if(!(gen%adap->intervalLength)) {
			const ContinuationCode x = this->AdaptPrec();
			if (x == FINISH_RUN)
				break;
			if (x == ABORT_IMMEDIATELY)
				return;
		}
		if(conf->checkpoint == true && ((gen % conf->saveevery) == 0))
			WriteStateFiles();

		if(stopwatch.SplitTime() > (int)conf->stoptime){
			outman.UserMessage("NOTE: ****Specified time limit (%d seconds) reached...", conf->stoptime);
			prematureTermination = true;
			break;
			}
		if(gen == conf->stopgen)
			outman.UserMessage("NOTE: ****Specified generation limit (%d) reached...", conf->stopgen);
#		ifdef INCLUDE_PERTURBATION
			if(pertMan->pertAbandoned == true && pertMan->restartAfterAbandon == true && (gen - pertMan->lastPertGeneration > pertMan->gensBeforeRestart)){
				params->starting_tree = "";
				pertMan->lastPertGeneration = gen;
				pertMan->pertAbandoned = false;
				pertMan->numPertsNoImprove = 0;
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

	if(conf->bootstrapReps == 0)
		outman.UserMessage("finished");
}
