// GARLI version 2.1 source code
// Copyright 2005-2014 Derrick J. Zwickl
// email: garli.support@gmail.com
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

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>

using namespace std;

#include "defs.h"
#include "datamatr.h"
#include "rng.h"
#include "nxsstring.h"
#include "errorexception.h"
#include "outputman.h"
#include "model.h"
#include "garlireader.h"
#include "stopwatch.h"

//extern ModelSpecification modSpec;

#define MAX_TAXON_LABEL		100

extern rng rnd;
extern OutputManager outman;

bool my_pair_compare(pair<int, int> fir, pair<int,int> sec) {return fir.second < sec.second;}

int SitePattern::numTax;
int SitePattern::maxNumStates;

int		numCompares;

bool SitePattern::operator<(const SitePattern &rhs) const{
	//zero state sites (all missing) will now be shuffled to the start (previously the end) and removed later
	//potentially constant sites always need to come just after that
	//sorting will first be by number of states (fast), then by the state vectors (slow)
	assert(numStates > -1 && rhs.numStates > -1);
	numCompares++;
	
	if(numStates < rhs.numStates)
		return true;
	if(numStates > rhs.numStates)
		return false;

	assert(stateVec.empty() == false);
	assert(stateVec.size() == rhs.stateVec.size());

	//this lexigraphically compares the vector contents
	if(stateVec < rhs.stateVec)
		return true;

	return false;
	}

bool SitePattern::operator==(const SitePattern &rhs) const{
	return (stateVec == rhs.stateVec);
	}


//CalcPatternTypeAndNumStates determines whether pattern a is constant, informative, or missing
//The passed in vector is used as scratch, and is assumed to already be of size maxNumStates
//This ALSO has the side effect of filling in the SitePattern::numStates field, which is necessary for sorting.
int SitePattern::CalcPatternTypeAndNumStates( vector<unsigned int> &stateCounts ){
	bool ambig = false;	//any total or partial ambiguity
	int nStates = 0;
	bool informative = false;
	bool constant = false;
	bool missing = false;

	//fill the scratch array with zeros
	std::fill(stateCounts.begin(), stateCounts.end(), 0);

	//count the number of times each state occurs, and whether there are any partially
	//ambiguous characters (currently only allowed for nuc data)
	unsigned char full_ambig = (maxNumStates == 4 ? 15 : maxNumStates);
	if(maxNumStates == 4){
		for(StateVector::iterator sit = stateVec.begin();sit != stateVec.end();sit++){
			unsigned char c = *sit;
			if(c != full_ambig && (c & (c - 1))){
				ambig = true;
				break;
				}
			else if(c != full_ambig){
				stateCounts[(c > 1) + (c > 2) + (c > 4)]++;
				}
			}
		}
	else {
		for(StateVector::iterator sit = stateVec.begin();sit != stateVec.end();sit++){
			unsigned char c = *sit;
			if(c != full_ambig){
				stateCounts[c]++;
				}
			}
		}

	if(!ambig){
		//no partial ambiguity (all AA and codon will come this way) 
		//without ambiguity, having 2+ states with 2+ counts means informativeness
		int numDoubles = 0;
		for(int s = 0; s < maxNumStates; s++ ){
			if(stateCounts[s] > 0){
				nStates++;
				if(stateCounts[s] > 1){
					numDoubles++;
					}
				}
			}

		if(nStates == 0){
			missing = true;
			assert(numDoubles == 0);
			}
		else if(nStates == 1){
			constant = true;
			assert(numDoubles < 2);
			}
		else{
			if(numDoubles > 1){
				informative = true;
				}
			}
		}
	else{
		assert(maxNumStates == 4);
		//this very convoluted scheme (worked out by Mark) must be used to determine informativeness
		//if partial ambiguity is allowed (only for nuc data currently)
		multiset<unsigned char> pat;
		unsigned char conStates = 15;
		for(StateVector::iterator sit = stateVec.begin();sit != stateVec.end();sit++){
			unsigned char c = *sit;
			pat.insert(c);
			conStates &= c;
			}

		//constant sites are possible with partial ambiguity if some resolution gives a single state
		if(conStates){
			if(conStates == 15) 
				missing = true;
			else 
				constant = true;
			}
		else{
			vector< pair<int, int> > stateScores;
			for(unsigned state=0;state < 4;state++){
				int sc = 0;
				for(multiset<unsigned char>::iterator it=pat.begin();it != pat.end();it++){
					if(!((*it) & (1 << state))){
						sc++;
						}
					}
				stateScores.push_back(pair<int, int>(state, sc));
				}
			sort(stateScores.begin(), stateScores.end(), my_pair_compare);
			int minStar = stateScores[0].second;
			if(minStar > 1){
				set<unsigned char> uPat;
				for(multiset<unsigned char>::iterator it=pat.begin();it != pat.end();it++)
					uPat.insert(*it);
				int minScore = MinScore(uPat, minStar);
				if(minScore < minStar){
					informative = true;
					}
				}
			}
		}

	if(missing){
		type = MISSING;
		nStates = 0;
		}
	else if(constant){
		type = CONSTANT;
		nStates = 1;
		}
	else if(informative){
		type = INFORMATIVE;
		nStates = max(2, nStates);
		}
	else{
		type = UNINFORM_VARIABLE;
		nStates = max(2, nStates);
		}

	//Note that numStates here may not be the true number of states in the
	//case of ambiguity, but it really only matters that it is accurate in 
	//discriminating 0/1/1+ states because code elsewhere depends on it.
	numStates = nStates;
	return type;
	}

//this is used for determining informative sites when there is partial ambiguity
int SitePattern::MinScore(set<unsigned char> patt, int bound, unsigned char bits/*=15*/, int prevSc/*=0*/) const{
	if(patt.size() == 0) return 0;
	int min_sc_this_lvl = 9999;
	int curr_sc_this_lvl = 9999;
	for(unsigned s2 = 0;s2 < 4;s2++){
		unsigned char thisBit = (1 << s2);
		if(bits & thisBit){
			set<unsigned char> remaining;
			for(set<unsigned char>::iterator it=patt.begin();it != patt.end();it++){
				if(!(*it & thisBit)) remaining.insert(*it);
				}
			if(remaining.size() > 0){
				if(prevSc + 1 < bound)
					curr_sc_this_lvl = 1 + MinScore(remaining, bound, bits & ~thisBit, prevSc+1);
				else 
					curr_sc_this_lvl = bound - prevSc;
				}
			else return 0;
			
			if(curr_sc_this_lvl < min_sc_this_lvl)
				min_sc_this_lvl = curr_sc_this_lvl;
			if(min_sc_this_lvl == 0 || min_sc_this_lvl + prevSc < bound)
				return min_sc_this_lvl;
			}
		}
	return min_sc_this_lvl;
	}

//Collapse merges like patterns, transfering over the counts and site numbers represented by each sucessive identical column.
//Patterns that are assigned zero counts here will be removed in Pack(), but will still contribute to the numNonMissingRealSitesInOrigMatrix, except
//for those with zero states (= missing)
void PatternManager::NewCollapse(){
	list<SitePattern>::iterator first;
	list<SitePattern>::iterator second  = patterns.begin();

	while(second != patterns.end()){
		first = second++;
		while(second != patterns.end() && (*first == *second)){
			(*first).count += (*second).count;
			(*first).siteNumbers.insert((*first).siteNumbers.end(), (*second).siteNumbers.begin(), (*second).siteNumbers.end());
			//if a wtset was used, this definitely doesn't need to be the case
			//assert((*first).count == (*first).siteNumbers.size());
			(*second).count = 0;
			second++;
			}
		}
	//this will zero the count of all missing pats, which will make them not get put into the uniquePatterns
	//list in NewPack
	for(list<SitePattern>::iterator pit = patterns.begin();pit != patterns.end();pit++){
		if((*pit).numStates == 0){
			(*pit).count = 0;
			}
		}
	}

void PatternManager::NewSort(){
	//this is the stl list sort function, using SitePattern::operator<
	patterns.sort();
	}

// This version of pack copies unique patterns from the patterns list into the uniquePatterns list
void PatternManager::NewPack(){
	for(list<SitePattern>::iterator pit = patterns.begin();pit != patterns.end();pit++){
		if(pit->numStates > 0){
			if(pit->count > 0){ 
				uniquePatterns.push_back(*pit);
				}
			}
		}
	pman_numPatterns = uniquePatterns.size();
	compressed = true;
	}

//This does all necessary processing in the patman (assuming that it has already been filled with data)
//up to the point when the compressed matrix can be copied back into 
//this will only be used for nuc/AA/codon data
void PatternManager::ProcessPatterns(){
	CalcPatternTypesAndNumStates();
	NewSort();
	NewCollapse();
	NewPack();
	NewDetermineConstantSites();
	}

//it would really make more sense to do this after packing, but the number of states
//is needed in pattern comparison in sorting.  This also does what Summarize used to, 
//filling the counts of various types of patterns
//THIS DOES NOT CURRENTLY SUPPORT CONDITIONING PATTERNS!
void PatternManager::CalcPatternTypesAndNumStates(){
	//this is just a scratch array to be used repeatedly in PatternType	
	vector<unsigned int> s(maxNumStates);

	pman_numMissingChars = pman_numConstantChars = pman_numInformativeChars = pman_numUninformVariableChars = pman_numNonMissingRealCountsInOrigMatrix = 0;

	pman_numRealSitesInOrigMatrix = patterns.size();
	for(list<SitePattern>::iterator pit = patterns.begin();pit != patterns.end();pit++){
		int t = pit->CalcPatternTypeAndNumStates(s);
		//Fixed 2 bugs - It is important to calculate numNonMissingRealCountsInOrigMatrix here from counts of the  
		//the generally unpacked data, because it could effectively be partially packed due to the use of a wtset, 
		//but also need to keep separate track of the number of columns in the orig matrix with numRealSitesInOrigMatrix
		if( t != SitePattern::MISSING )
			pman_numNonMissingRealCountsInOrigMatrix += pit->count;

		if( t == SitePattern::MISSING )
			pman_numMissingChars += pit->count;
		else if( t == SitePattern::CONSTANT )
			pman_numConstantChars += pit->count;
		else if( t == SitePattern::INFORMATIVE )
			pman_numInformativeChars += pit->count;
		else{
			assert(t == SitePattern::UNINFORM_VARIABLE);
			pman_numUninformVariableChars += pit->count;
			}
		}
	pman_numNonMissingChars = pman_numRealSitesInOrigMatrix - pman_numMissingChars;
	if( pman_numNonMissingChars == 0 ){
		throw ErrorException("Matrix is made up entirely of missing characters (?, -, or N)!");
		}
	}

//note where all of the constant sites are, and what state they are.
//this is kind of ugly, but will never be rate limiting
void PatternManager::NewDetermineConstantSites(){
	assert(compressed);
	lastConstant=-1;
	list<SitePattern>::iterator pat = uniquePatterns.begin();
	assert(pat->numStates > 0);
	while(pat != uniquePatterns.end() && pat->numStates == 1){
		lastConstant++;
		pat++;
		}
	
	int t = 0;
	int thisCon = 0;
	if(maxNumStates == 4){
		for(pat = uniquePatterns.begin();thisCon++ <= lastConstant;pat++){
			t = 0;
			char c=15;
			while(t < numTax){
				char ch = pat->stateVec[t];
				c = c & ch;
				t++;
				}
			assert(c != 0);
			pat->constStates = c;
			}
		}
	else{//not allowing ambiguity for codon/AA's, so this is a bit easier
		for(pat = uniquePatterns.begin();thisCon++ <= lastConstant;pat++){
			t = 0;
			char c = maxNumStates;
			do{
				c = pat->stateVec[t];
				t++;
				}while(c == maxNumStates && t < numTax);
			assert(t <= numTax);
			pat->constStates = c;
			}
		}
	}

//The following are for copying the results of the pattern processing back into the old fields of DataMatrix

//This takes the unique pattern types and uses their siteNumbers vector to map back to the original
//ordering of sites, as used to tbe stored in the number array.
void PatternManager::FillNumberVector(vector<int> &nums) const{
	if(nums.size() != patterns.size()){
		nums.clear();
		nums.resize(patterns.size());
		}
	
	//this is necessary so that all missing patterns, which should already have been removed from
	//uniquePatterns, will properly show up as -1 in the number array, since they will not be overwritten
	//with other values below
	for(vector<int>::iterator nit = nums.begin();nit != nums.end();nit++)
		(*nit) = -1;

	int p = 0;
	for(list<SitePattern>::const_iterator pit = uniquePatterns.begin();pit != uniquePatterns.end();pit++){
		for(vector<int>::const_iterator nit = (*pit).siteNumbers.begin(); nit != (*pit).siteNumbers.end();nit++)
			nums[*nit] = p;
		p++;
		}
	}

void PatternManager::FillCountVector(vector<int> &counts) const{
	counts.clear();
	for(list<SitePattern>::const_iterator pit = uniquePatterns.begin();pit != uniquePatterns.end();pit++){
		counts.push_back((*pit).count);
		}
	}

void PatternManager::FillNumStatesVector(vector<int> &ns) const{
	ns.clear();
	for(list<SitePattern>::const_iterator pit = uniquePatterns.begin();pit != uniquePatterns.end();pit++){
		ns.push_back((*pit).numStates);
		}
	}

void PatternManager::FillConstStatesVector(vector<int> &cs) const{
	int c = 0;
	for(list<SitePattern>::const_iterator pit = uniquePatterns.begin();pit != uniquePatterns.end();pit++){
		cs.push_back((*pit).constStates);
		c++;
		}
	}

//Takes the data out of the SitePattern list and copies into the DataMatrix 2d matrix
void PatternManager::FillTaxaXCharMatrix(unsigned char **mat) const{
	for(int t = 0;t < numTax;t++){
		int c = 0;
		for(list<SitePattern>::const_iterator cit = uniquePatterns.begin();cit != uniquePatterns.end();cit++){
			mat[t][c++] = (*cit).stateVec[t];
			}
		}
	}

void PatternManager::FillIntegerValues(int &_nMissing, int &_nConstant, int &_nVarUninform, int &_nInformative, int &_lastConst, int &_numRealSitesInOrigMatrix, int &_numNonMissingRealCountsInOrigMatrix, int &_numNonMissingRealSitesInOrigMatrix, int &_numPatterns) const {
	_numRealSitesInOrigMatrix = pman_numRealSitesInOrigMatrix;
	_numNonMissingRealCountsInOrigMatrix = pman_numNonMissingRealCountsInOrigMatrix;
	_numNonMissingRealSitesInOrigMatrix = pman_numNonMissingChars;
	_numPatterns = pman_numPatterns;

	_nMissing = pman_numMissingChars;
	_nConstant = pman_numConstantChars;
	_nVarUninform = pman_numUninformVariableChars;
	_nInformative = pman_numInformativeChars;
	_lastConst = lastConstant;
	}

vector<IdenticalColumnRange> PatternManager::FindIdenticalAlignmentColumns(const PatternManager &other) const{
	vector<SitePattern> basePos1(patterns.size());
	vector<SitePattern> basePos2(other.patterns.size());

	//assuming nucleotide for now
	unsigned char ambigState = 15;
	//make dummy site patterns consisting of the base numbers allong each sequence, begining at 1.
	//gaps are the negative of the next base number i.e.
	//-ACG-T = -1123-34
	//GAGG-- = 1234-5-5
	//identity of site patterns for different alignments means exactly identical alignment columns, NOT just 
	//in the observed states, including gaps, i.e. for the following two alignments NO columns would be identical
	//AAA		AAA
	//--A		A--
	//AAA		AAA
	int baseNum, colNum;
	for(int tax = 0;tax < numTax;tax++){
		baseNum = 1;
		colNum = 0;
		for(list<SitePattern>::const_iterator sit = patterns.begin();sit != patterns.end();sit++)
			basePos1[colNum++].AddChar( (*sit).stateVec[tax] == ambigState ? -baseNum : baseNum++); 
			
		baseNum = 1;
		colNum = 0;
		for(list<SitePattern>::const_iterator sit = other.patterns.begin();sit != other.patterns.end();sit++)
			basePos2[colNum++].AddChar( (*sit).stateVec[tax] == ambigState ? -baseNum : baseNum++); 
		}

	//Things get confusing below.  IdenticalColumnPair consists of one column index in one alignment and one in another. 
	vector<IdenticalColumnPair> columnMatches;
	int lastMatch, index1, index2;
	lastMatch = index1 = index2 = 0;
	for(vector<SitePattern>::iterator pat1 = basePos1.begin();pat1 != basePos1.end();pat1++){
		index2 = lastMatch;
		for(vector<SitePattern>::iterator pat2 = basePos2.begin() + lastMatch;pat2 != basePos2.end();pat2++){
			if(*pat2 == *pat1){
				columnMatches.push_back(IdenticalColumnPair(index1, index2));
				lastMatch = index2;
				break;
				}
			else
				index2++;
			}
		index1++;
		}

	//IdenticalColumnRange is a pair of ColumnRanges, with each pair being WITHIN one of the alignments, i.e.
	//((firstStart, firstEnd), (secondStart, secondEnd))
	vector<IdenticalColumnRange> identicalRanges;
	for(vector<IdenticalColumnPair>::iterator startPair = columnMatches.begin();startPair != columnMatches.end();){
		vector<IdenticalColumnPair>::iterator currentPair = startPair;
		vector<IdenticalColumnPair>::iterator nextPair = startPair + 1;
		while((*nextPair).first == (*currentPair).first + 1 && (*nextPair).second == (*currentPair).second + 1){
			currentPair++;
			nextPair++;
			if(nextPair == columnMatches.end())
				break;
			}
		//identicalRanges.push_back(IdenticalColumnRange(*startPair, *currentPair));
		identicalRanges.push_back(IdenticalColumnRange(ColumnRange((*startPair).first, (*currentPair).first), ColumnRange((*startPair).second, (*currentPair).second)));
		startPair = nextPair;
		}

	return identicalRanges;
	}

void DataMatrix::CreateMatrixFromOtherMatrix(const DataMatrix &other, int startIndex, int endIndex){
	NewMatrix(other.patman.numTax, other.patman.patterns.size());
	patman.Initialize(other.patman.numTax, other.patman.maxNumStates);

	list<SitePattern>::const_iterator fromStart = other.patman.patterns.begin();
	list<SitePattern>::const_iterator fromEnd = fromStart;
	advance(fromStart, startIndex);
	advance(fromEnd, endIndex + 1);
	for(list<SitePattern>::const_iterator copypat = fromStart;copypat != fromEnd;copypat++){
		SitePattern thisPat(*copypat);
		for(vector<int>::iterator sit = thisPat.siteNumbers.begin();sit != thisPat.siteNumbers.end();sit++)
			(*sit) -= startIndex;
		patman.AddPattern(thisPat);
		}

	for(int tax = 0;tax < other.nTax;tax++)
		this->SetTaxonLabel(tax, other.TaxonLabel(tax));


	}

void DataMatrix::OutputDataSummary() const{
	//outman.UserMessage("\n#######################################################");
	outman.UserMessage("\tSummary of data:");
	outman.UserMessage("\t  %d sequences.", NTax());
	outman.UserMessage("\t  %d constant characters.", NConstant() - numConditioningPatterns);
	outman.UserMessage("\t  %d parsimony-informative characters.", NInformative());
	outman.UserMessage("\t  %d uninformative variable characters.", NVarUninform());
	int total = NConstant() + NInformative() + NVarUninform() - numConditioningPatterns;
	if(NMissing() > 0){
		outman.UserMessage("\t  %d characters were completely missing or ambiguous (removed).", NMissing());
		//outman.UserMessage("\t  %d total characters (%d before removing empty columns).", total, GapsIncludedNChar() - numConditioningPatterns);
		outman.UserMessage("\t  %d total characters (%d before removing empty columns).", total, numNonMissingRealCountsInOrigMatrix + NMissing());
		}
	else outman.UserMessage("\t  %d total characters.", total);

	assert(total == numNonMissingRealCountsInOrigMatrix);

	outman.UserMessage("\t  %d unique patterns in compressed data matrix.", NChar() - numConditioningPatterns);
	outman.flush();
	}

void DataMatrix::ProcessPatterns() {
	Stopwatch stoppy;
	stoppy.Start();
	if(usePatternManager){
		patman.ProcessPatterns();
		GetDataFromPatternManager();
		patman.Reset();
		}
	else{
		Summarize();
		Collapse();
		DetermineConstantSites();
		}
	CalcEmpiricalFreqs();
	ReserveOriginalCounts();
	OutputDataSummary();
	int t = stoppy.SplitTime();
	/*There isn't much point in outputting all of this clutter
	if(t == 0)
		outman.UserMessage("\tPattern processing required < 1 second");
	else
		outman.UserMessage("\tPattern processing required %d second(s)", stoppy.SplitTime());
	if(numCompares > 0)
		outman.UserMessage("\t%d pattern comparisons were needed", numCompares);
	outman.UserMessage("");
	*/
	}

//this pulls all of the processed data back out of the patman into the old fields of DataMatrix
void DataMatrix::GetDataFromPatternManager(){
	ResizeCharacterNumberDependentVariables(patman.NChar()) ;
	patman.FillNumberVector(newNumber);
	patman.FillCountVector(newCount);
	patman.FillNumStatesVector(newNumStates);
	patman.FillConstStatesVector(newConstStates);
	patman.FillIntegerValues(numMissingChars, numConstantChars, numVariableUninformChars, numInformativeChars, lastConstant, numRealSitesInOrigMatrix, numNonMissingRealCountsInOrigMatrix, numNonMissingRealSitesInOrigMatrix, numPatterns);
	patman.FillTaxaXCharMatrix(matrix);
	if(patman.compressed)
		dense = 1;
	}

DataMatrix::~DataMatrix(){
	if( count ) MEM_DELETE_ARRAY(count); // count is of length numPatterns
	if( numStates ) MEM_DELETE_ARRAY(numStates); // numStates is of length numPatterns
	if( number ) MEM_DELETE_ARRAY(number); // number is of length numPatterns
	if( origDataNumber ) MEM_DELETE_ARRAY(origDataNumber); // origDataNumber is of length numPatterns
	if( taxonLabel ) {
		for( int j = 0; j < nTaxAllocated; j++ )
			MEM_DELETE_ARRAY( taxonLabel[j] ); // taxonLabel[j] is of length strlen(taxonLabel[j])+1
	    MEM_DELETE_ARRAY(taxonLabel); // taxonLabel is of length nTax
	}
	if( matrix ) {
		for( int j = 0; j < nTaxAllocated; j++ )
			MEM_DELETE_ARRAY(matrix[j]); // matrix[j] is of length numPatterns
		MEM_DELETE_ARRAY(matrix); // matrix is of length nTax
	}
	if(constStates!=NULL) delete []constStates;
	if(origCounts!=NULL) delete []origCounts;
	}

void DataMatrix::SetTaxonLabel(int i, const char* s)
{
	if( taxonLabel && (i < nTax) )
		ReplaceTaxonLabel(i, s);
}

void DataMatrix::ReplaceTaxonLabel( int i, const char* s )
{
	assert( taxonLabel );
	if( taxonLabel[i] ) {
		MEM_DELETE_ARRAY(taxonLabel[i]); // taxonLabel[i] is of length strlen(taxonLabel[i])+1
	}
	int newLength = (strlen(s)+1);
	if(newLength > MAX_TAXON_LABEL) throw ErrorException("Sorry, taxon name %s for taxon #%d is too long (max length=%d)", s, i+1, MAX_TAXON_LABEL);
	MEM_NEW_ARRAY(taxonLabel[i],char,newLength);
	strcpy(taxonLabel[i], s);
}

//
// PositionOf returns position (starting from 0) of taxon whose name
// matches the string s in the taxonLabel list
//
int DataMatrix::PositionOf( char* s ) const
{
	int i;

	for( i = 0; i < nTax; i++ ) {
		if( strcmp( taxonLabel[i], s ) == 0 ) break;
	}

	assert( i < nTax );

	return i;
}

int DataMatrix::MinScore(set<unsigned char> patt, int bound, unsigned char bits/*=15*/, int prevSc/*=0*/) const{
	if(patt.size() == 0) return 0;
	int min_sc_this_lvl = 9999;
	int curr_sc_this_lvl = 9999;
	for(unsigned s2 = 0;s2 < 4;s2++){
		unsigned char thisBit = (1 << s2);
		if(bits & thisBit){
			set<unsigned char> remaining;
			for(set<unsigned char>::iterator it=patt.begin();it != patt.end();it++){
				if(!(*it & thisBit)) remaining.insert(*it);
				}
			if(remaining.size() > 0){
				if(prevSc + 1 < bound)
					curr_sc_this_lvl = 1 + MinScore(remaining, bound, bits & ~thisBit, prevSc+1);
				else 
					curr_sc_this_lvl = bound - prevSc;
				}
			else return 0;
			
			if(curr_sc_this_lvl < min_sc_this_lvl)
				min_sc_this_lvl = curr_sc_this_lvl;
			if(min_sc_this_lvl == 0 || min_sc_this_lvl + prevSc < bound)
				return min_sc_this_lvl;
			}
		}
	return min_sc_this_lvl;
	}

//
// PatternType determines whether pattern k is constant, informative, or missing
//it used to try to determine autapomorphies, although not correctly
//
int DataMatrix::PatternType( int k , unsigned int *stateCounts) const{
	assert(k < numPatterns);
	if( k >= numPatterns )
		return 0;
	int retval;

	bool ambig = false;	//any total or partial ambiguity
	int nStates = 0; 
	bool informative = false;
	bool constant = false;
	bool missing = false;

	//fill the scratch array with zeros
	memset(stateCounts, 0x00, maxNumStates * sizeof(*stateCounts));

	//count the number of times each state occurs, and whether there are any partially
	//ambiguous characters (currently only allowed for nuc data)
	unsigned char full_ambig = (maxNumStates == 4 ? 15 : maxNumStates);
	if(maxNumStates == 4){
		for(int t = 0; t < nTax; t++ ){
			unsigned char c = Matrix( t, k );
			if(c != full_ambig && (c & (c-1))){
				ambig = true;
				break;
				}
			else if(c != full_ambig)
				stateCounts[(c > 1) + (c > 2) + (c > 4)]++;
			}
		}
	else {
		for(int t = 0; t < nTax; t++ ){
			unsigned char c = Matrix( t, k );
			if(c != full_ambig)
				stateCounts[c]++;
			}
		}

	if(!ambig){
		//no partial ambiguity (all AA and codon will come this way) 
		//without ambiguity, having 2+ states with 2+ counts means informativeness
		int numDoubles = 0;
		for(int s = 0; s < maxNumStates; s++ ){
			if(stateCounts[s] > 0){
				nStates++;
				if(stateCounts[s] > 1) numDoubles++;
				}
			}

		if(nStates == 0){
			missing = true;
			assert(numDoubles == 0);
			}
		else if(nStates == 1){
			constant = true;
			assert(numDoubles < 2);
			}
		else{
			if(numDoubles > 1) informative = true;
			}
		}
	else{
		//this very convoluted scheme must be used to determine informativeness
		//if ambiguity is allowed (only for nuc data currently)
		multiset<unsigned char> patt;
		unsigned char conStates = 15;
		for(int t = 0;t < nTax;t++){
			unsigned char c = Matrix( t, k );
			patt.insert(c);
			conStates &= c;
			}

		//constant sites are possible with ambiguity of some resolution gives a single state
		if(conStates){
			if(conStates == 15) missing = true;
			else constant = true;
			}
		else{
			vector< pair<int, int> > stateScores;
			for(unsigned state=0;state < 4;state++){
				int sc = 0;
				for(multiset<unsigned char>::iterator it=patt.begin();it != patt.end();it++){
					if(!((*it) & (1 << state))){
						sc++;
						}
					}
				stateScores.push_back(pair<int, int>(state, sc));
				}
			sort(stateScores.begin(), stateScores.end(), my_pair_compare);
			int minStar = stateScores[0].second;
			if(minStar > 1){
				set<unsigned char> uPatt;
				for(multiset<unsigned char>::iterator it=patt.begin();it != patt.end();it++)
					uPatt.insert(*it);
				int minScore = MinScore(uPatt, minStar);
				if(minScore < minStar){
					informative = true;
					}
				}
			}
		}

	if(missing){
		retval = PT_MISSING;
		nStates = 0;
		}
	else if(constant){
		retval = PT_CONSTANT;
		nStates = 1;
		}
	else if(informative){
		retval = PT_INFORMATIVE | PT_VARIABLE;
		nStates = max(2, nStates);
		}
	else{
		retval = PT_VARIABLE;
		nStates = max(2, nStates);
		}

/*	ofstream deb;
	if(k==0) deb.open("pat.log");
	else deb.open("pat.log", ios::app);
	deb << k << "\t" << constant << "\t" << informative << "\t" << nStates << "\n";
	deb.close();
*/
	//Note that numStates here may not be the true number of states in the
	//case of ambiguity, but it really only matters that it is accurate in 
	//discriminating 0/1/1+ states because code elsewhere depends on it.
	numStates[k] = nStates;
	return retval;
	}

//
// Summarize tallies number of constant, informative, and autapomorphic characters
//
void DataMatrix::Summarize(){
	assert( numPatterns > 0 );

	numMissingChars = numConstantChars = numInformativeChars = numVariableUninformChars = numNonMissingRealCountsInOrigMatrix = 0;

    //this is just a scratch array to be used repeatedly in PatternType
    vector<unsigned int> s(maxNumStates);
	
	numRealSitesInOrigMatrix = numPatterns - numConditioningPatterns;
	for(int k = 0; k < numPatterns; k++ ) {
		int ptFlags = PatternType(k, &s[0]);
		//Fixed 2 bugs - It is important to calculate numNonMissingRealCountsInOrigMatrix here from counts of the  
		//the generally unpacked data, because it could effectively be partially packed due to the use of a wtset, 
		//but also need to keep separate track of the number of columns in the orig matrix with numRealSitesInOrigMatrix
		if( ptFlags != PT_MISSING && k >= numConditioningPatterns)
			numNonMissingRealCountsInOrigMatrix += count[k];
		if( ptFlags == PT_MISSING )
			numMissingChars += count[k];
		else if( ptFlags & PT_CONSTANT )
			numConstantChars += count[k];
		else if( ptFlags & PT_INFORMATIVE )
			numInformativeChars += count[k];
		else{
			assert(ptFlags & PT_VARIABLE);
			numVariableUninformChars += count[k];
			}
		}
	numNonMissingRealSitesInOrigMatrix = numRealSitesInOrigMatrix - numMissingChars;
	if( numConstantChars + numInformativeChars + numVariableUninformChars == 0 ){
		throw ErrorException("Matrix is made up entirely of missing characters (?, -, or N)!");
		}
	}

//
// NewMatrix deletes old matrix, taxonLabel, count, and number
// arrays and creates new ones
//
void DataMatrix::NewMatrix( int taxa, int sites )
{
	//allocate an extra taxon unless there previously wasn't one
	int extraTax = 1;
	if(nTaxAllocated > 0){
		extraTax = nTaxAllocated - nTax;
		}

	// delete taxon labels
	if( taxonLabel ) {
		int i;
		for( i = 0; i < nTaxAllocated; i++ )
			MEM_DELETE_ARRAY(taxonLabel[i]); // taxonLabel[i] is of length strlen(taxonLabel[i])+1
		MEM_DELETE_ARRAY(taxonLabel); // taxonLabel is of length nTax
		}

	// create new array of taxon label pointers
	if( taxa > 0 ) {
		MEM_NEW_ARRAY(taxonLabel,char*,taxa + extraTax);
		for( int i = 0; i < taxa + extraTax; i++ )
			taxonLabel[i] = NULL;
		}

	// delete data matrix and count and number arrays
	if( matrix ) {
		int j;
		for( j = 0; j < taxa + extraTax; j++ )
			MEM_DELETE_ARRAY(matrix[j]); // matrix[j] has length numPatterns
		MEM_DELETE_ARRAY(matrix); // matrix has length nTax
		}

	if(usePatternManager == false){
		if( count ) {
			MEM_DELETE_ARRAY(count); //count has length numPatterns
			}
		if( numStates ) {
			MEM_DELETE_ARRAY(numStates); // numStates has length numPatterns
			}
		}
	//yarg - this is used even when usePatMan == true, when converting from nuc to AA matrix
	if( number ) {
		MEM_DELETE_ARRAY(number); // number has length numPatterns
		}
	if( origDataNumber ) {
		MEM_DELETE_ARRAY(origDataNumber); // origDataNumber has length numPatterns
		}
	// create new data matrix, and new count and number arrays
	// all counts are initially 1, and characters are numbered
	// sequentially from 0 to numPatterns-1
	if( taxa > 0 && sites > 0 ) {
		MEM_NEW_ARRAY(matrix,unsigned char*,taxa + extraTax);
		MEM_NEW_ARRAY(number,int,sites);
		MEM_NEW_ARRAY(origDataNumber,int,sites);
		for( int j = 0; j < sites; j++ ) {
			number[j] = j;
			//number[j] = ( j < numConditioningPatterns ? -1 : j - numConditioningPatterns);
			//in the case of conditioning patterns or partitioning this will be updated later anyway
			origDataNumber[j] = j;
			}
		if(usePatternManager == false){
			MEM_NEW_ARRAY(count,int,sites);
			MEM_NEW_ARRAY(numStates,int,sites);

			for( int j = 0; j < sites; j++ ) {
				count[j] = 1;
				numStates[j] = 1;
				}
			}
		for( int i = 0; i < taxa + extraTax; i++ ) {
			matrix[i]=new unsigned char[sites];
			//MEM_NEW_ARRAY(matrix[i],unsigned char,sites);
			//memset( matrix[i], 0xff, taxa*sizeof(unsigned char) );
			memset( matrix[i], 0xff, sites*sizeof(unsigned char) );
			}
		}

	// set dimension variables to new values
	nTax = taxa;
	nTaxAllocated = nTax + extraTax;
	//these will likely be updated later
	numRealSitesInOrigMatrix = numNonMissingRealSitesInOrigMatrix = sites - numConditioningPatterns;
	nonZeroCharCount = numPatterns = sites;
	}

void DataMatrix::ResizeCharacterNumberDependentVariables(int nCh) {
	//ONLY CALL THIS BEFORE GETTING DATA FROM PATMAN, OTHERWISE SOME VARIABLES
	//WILL BE WRONG!!!
	numPatterns = nCh;

	// delete data matrix and count and number arrays
	if( matrix ) {
		int j;
		for( j = 0; j < nTaxAllocated; j++ )
			MEM_DELETE_ARRAY(matrix[j]); // matrix[j] has length numPatterns
		MEM_DELETE_ARRAY(matrix); // matrix has length nTax
	}

	if( count ) {
		MEM_DELETE_ARRAY(count); //count has length numPatterns
	}
	if( numStates ) {
		MEM_DELETE_ARRAY(numStates); // numStates has length numPatterns
	}

	// create new data matrix, and new count and number arrays
	// all counts are initially 1, and characters are numbered
	// sequentially from 0 to numPatterns-1
	if( numPatterns > 0 ) {
		MEM_NEW_ARRAY(matrix,unsigned char*,nTaxAllocated);
		MEM_NEW_ARRAY(count,int,numPatterns);
		MEM_NEW_ARRAY(numStates,int,numPatterns);

		for( int j = 0; j < numPatterns; j++ ) {
			count[j] = 1;
			numStates[j] = 1;
		}
		for( int i = 0; i < nTaxAllocated; i++ ) {
			matrix[i]=new unsigned char[numPatterns];
			memset( matrix[i], 0xff, numPatterns*sizeof(unsigned char) );
			}
		}

	//set dimension variables to new values, which actually MUST be updated elsewhere to be correct
	//see note at top of func
	nonZeroCharCount = numRealSitesInOrigMatrix = numNonMissingRealSitesInOrigMatrix = numPatterns;
	}

//deprecated
DataMatrix& DataMatrix::operator =(const DataMatrix& d){
	assert(0);
	NewMatrix( d.NTax(), d.NChar() );

	int i, j;
	for( i = 0; i < nTax; i++ ) {
		SetTaxonLabel(i, d.TaxonLabel(i) );
	}

	for( j = 0; j < numPatterns; j++ ) {
		SetCount(j, d.Count(j) );
		origCounts[j] = d.origCounts[j];
		number[j] = d.Number(j);
		numStates[j] = d.NumStates(j);
	}

	for( i = 0; i < nTax; i++ ) {
		for( j = 0; j < numPatterns; j++ )
			SetMatrix(i, j, d.Matrix(i, j));
	}

	return *this;
	}

//
// Pack simply deletes sites having a count of zero
//
void DataMatrix::Pack(){

	int i, j, newNChar = 0;

	// determine dimensions of new matrix
	for( j = 0; j < numPatterns; j++ ) {
		if( count[j] )
			newNChar++;
	}

	//DEBUG - something was going wrong and causing crashes in some cases (only AA's?) when a new matrix
	//was created with the same dimensions as the original.  Haven't figured out why yet, 
	//but this avoids the crash at least.
	if(newNChar == numPatterns){
		dense = true;
		return;
		}

	// create new matrix and count arrays and fill
	unsigned char** newMatrix;
        MEM_NEW_ARRAY(newMatrix,unsigned char*,nTaxAllocated);
	int* newCount;
        MEM_NEW_ARRAY(newCount,int,newNChar);
	int* newNumStates;
        MEM_NEW_ARRAY(newNumStates,int,newNChar);

	for( i = 0; i < nTaxAllocated; i++ )
		 MEM_NEW_ARRAY(newMatrix[i],unsigned char,newNChar);


	i = 0;
	for( j = 0; j < numPatterns; j++ ) {
		if( count[j] ) {
			for( int k = 0; k < nTax; k++ )
				newMatrix[k][i] = matrix[k][j];
			newCount[i] = count[j];
			newNumStates[i] = numStates[j];
			i++;
			}
		else{//as we remove columns, shift all the greater numbers over
			for(int c=0;c < numPatterns;c++){
				if(number[c] >= i) number[c]--;
				}
			}
		}

	// delete old matrix and count arrays
	if( count ) MEM_DELETE_ARRAY(count); // count has length numPatterns
	if( numStates ) MEM_DELETE_ARRAY(numStates); // numStates has length numPatterns
	if( matrix ) {
		for( i = 0; i < nTaxAllocated; i++ )
			MEM_DELETE_ARRAY(matrix[i]); // matrix[i] has length numPatterns
		MEM_DELETE_ARRAY(matrix); // matrix has length nTax
        }

	// set count, matrix and numStates to their new counterparts
	count = newCount;
	numStates = newNumStates;
	matrix = newMatrix;
	numPatterns = newNChar;
	nonZeroCharCount = numPatterns;
	}


void DataMatrix::DetermineConstantSites(){
	//note where all of the constant sites are, and what they are
	//this is kind of ugly, but will never be rate limiting
	lastConstant=-1;
	assert(numStates[0] > 0);
	while(numStates[lastConstant+1]==1) lastConstant++;
	
	constStates=new int[lastConstant+1];
	int t;
	if(maxNumStates == 4){
		for(int i=0;i<lastConstant+1;i++){
			t=0;
			char c=15;
			while(t<nTax){
				char ch=Matrix(t, i);
				c = c & ch;
				t++;
				}
			assert(c!=0);
			constStates[i]=c;
			}
		}
	else{//not allowing ambiguity for codon/AA's, so this is a bit easier
		for(int i=0;i<lastConstant+1;i++){
			t=0;
			char c = maxNumStates;
			do{
				c = Matrix(t, i);
				t++;
				}while(c == maxNumStates && t < nTax);
			assert(t <= nTax);
			constStates[i]=c;
			}
		}
	}

//
//	SwapCharacters swaps matrix column i with column j
//
void DataMatrix::SwapCharacters( int i, int j ){
	unsigned char tmp;
	for( int k = 0; k < nTax; k++ ) {
		tmp = Matrix( k, i );
		SetMatrix( k, i, Matrix( k, j ) );
		SetMatrix( k, j, tmp );
		}

	//swap pattern counts
	int c = count[i];
	count[i] = count[j];
	count[j] = c;

	//and the nStates
	int s=numStates[i];
	numStates[i]=numStates[j];
	numStates[j]=s;

	//and the number
	for(int c=0;c<numRealSitesInOrigMatrix + numConditioningPatterns;c++){
		if(number[c] == i) 
			number[c] = j;
		else if(number[c] == j) 
			number[c] = i;
		}
	}

void DataMatrix::BeginNexusTreesBlock(string &trans) const{
	//this outputs everything up through the translate table
	trans = "#NEXUS\n\nbegin trees;\ntranslate\n";

	char temp[500];
	for(int k=0;k<nTax;k++){
		//DO NOT call GetEscaped() or BlanksToUnderscores() on the names here - they are stored 
		//exactly as they should be output upon initial reading
		NxsString tnstr = TaxonLabel(k);
		if( k == nTax-1 )
			sprintf(temp, " %d %s;\n", (k+1), tnstr.c_str());
		else
			sprintf(temp, " %d %s,\n", (k+1), tnstr.c_str());
		trans += temp;
		}
	}

void DataMatrix::BeginNexusTreesBlock(ofstream &treeout) const{
	//this outputs everything up through the translate table
	treeout << "#NEXUS\n\nbegin trees;\ntranslate\n";
 	for(int k=0;k<nTax;k++){
		//DO NOT call GetEscaped() or BlanksToUnderscores() on the names here - they are stored 
		//exactly as they should be output upon initial reading
		treeout << "  " << (k+1);
		NxsString tnstr = TaxonLabel(k);
		treeout << "  " << tnstr.c_str();
		if( k == nTax-1 )
			treeout << ";\n";
		else
			treeout << ",\n";
		}
	treeout.flush();
	}

int DataMatrix::TaxonNameToNumber(const NxsString &name) const{\
	string nameStr = NxsToken::Tokenize(name)[0].GetToken();
	for(int i=0;i<nTax;i++){
		if(nameStr == (NxsToken::Tokenize(TaxonLabel(i)))[0].GetToken())
			return i+1;//indeces run 0->ntax-1, taxon numbers 1->ntax
			}
	return -1;
	}

//
// ComparePatterns returns:
//	 0		complete identity
//	-1		if i less than j
//	 1		if i greater than j
//
int DataMatrix::ComparePatterns( const int i, const int j ) const{		
	//DJZ 10/28/03 altering this to always put constant patterns at the start, which will
	//make implementing invariant sites much easier.  

	int cmp = 0;
		
	if(numStates[i]==1){
		if(numStates[j]==1){
			if(Matrix(0,i) < Matrix(0,j)) return -1;
			//else return 1;
			else{
				for( int k = 0; k < nTax; k++ ) {
					int same = ( Matrix( k, i ) == Matrix( k, j ) );
					if( !same )	{
						FLOAT_TYPE diff = ( (FLOAT_TYPE)Matrix( k, i ) - (FLOAT_TYPE)Matrix( k, j ) );
						cmp = ( diff < 0.0 ? -1 : 1 );
						break;
					}
				}
				return cmp;
				}
			}
		else return -1;
		}
	else if(numStates[j]==1){
		return 1;
		}
	
	
	for( int k = 0; k < nTax; k++ ) {
		int same = ( Matrix( k, i ) == Matrix( k, j ) );
		if( !same )	{
			FLOAT_TYPE diff = ( (FLOAT_TYPE)Matrix( k, i ) - (FLOAT_TYPE)Matrix( k, j ) );
			cmp = ( diff < 0.0 ? -1 : 1 );
			break;
		}
	}
	return cmp;
}

//
// Collapse merges like patterns
//
void DataMatrix::Collapse(){
	int i = 0, j = 1;
	assert(nonZeroCharCount == numPatterns);
	Sort();

	while( i < numPatterns ) {
		while( j < numPatterns && ComparePatterns( i, j ) == 0 ) {
			// pattern j same as pattern i
			count[i] += count[j];
			count[j] = 0;
			j++;
			}
		i = j++;
		}
		
	//DJZ 10/28/03 get rid of all missing patterns	
	int q=numPatterns-1;
	while(numStates[q] == 0){
		//This sets the number to -1 for an all missing site, indicating that none of the packed
		//matrix columns corresponds to it
		for(i = 0;i < numPatterns;i++){ 
			if(number[i]==q) 
				number[i]=-1;
			}
		count[q--]=0;
		//NO, this is now done in Summarize()!!
		//when all missing columns are deleted, remove them from the total number of characters
		//numNonMissingRealSitesInOrigMatrix--;
		}
	
	Pack();
	assert(nonZeroCharCount == numPatterns);
	}

//
// EliminateAdjacentIdenticalColumns sets the count of successive identical patterns
// in the original alignment to zero (usually applied to gaps)
// i.e., adjacent identical patterns count as only one observation
// 5/17/12 Ooops, changed to only do this for non-constant columns
void DataMatrix::EliminateAdjacentIdenticalColumns(){
	//this needs to happen here to know the number of state counts, but will be redone later
	Summarize();
	
	int i = 0, j = 1;
	assert(nonZeroCharCount == numPatterns);

	int numCombined = 0;
	while( i < NChar() ) {
		//need to avoid subtracting zero state chars here (blank cols) since the will be removed already
		//oops, and one state characters, since we don't want to collapse all present columns
		while( numStates[i] > 1 && j < NChar() && ComparePatterns( i, j ) == 0 ) {
			// pattern j same as pattern i
			count[j] = 0;
			//when columns are eliminated, remove them from the total number of characters
			numNonMissingRealSitesInOrigMatrix--;
			numCombined++;
			j++;
			}
		i = j++;
		}
	outman.UserMessage("	***%d IDENTICAL ADJACENT CHARACTERS ELIMINATED***", numCombined);
	}

//
//  BSort implements a simple bubblesort
//
void DataMatrix::BSort( int byCounts /* = 0 */ ){
	int swap, k;
	for( int i = 0; i < numPatterns-1; i++ ) {
		for( int j = i+1; j < numPatterns; j++ ) {
			if( byCounts )
				swap = ( count[i] < count[j] ? 1 : 0 );
			else
				swap = ( ComparePatterns( i, j ) > 0 ? 1 : 0 );
			if( swap ) {
				SwapCharacters( i, j );

				k = count[i];
				count[i] = count[j];
				count[j] = k;

				k = numStates[i];
				numStates[i] = numStates[j];
				numStates[j] = k;

				k = number[i];
				number[i] = number[j];
				number[j] = k;
			}
		}
	}
}

void DataMatrix::DebugSaveQSortState( int top, int bottom, int ii, int jj, int xx, const char* title )
{
	ofstream qsf( "qsstate.txt", ios::out | ios::app );
	qsf << endl << title << endl;

	int i, j;
	for( j = 0; j < numPatterns; j++ )
	{
		qsf << setw(6) << j << "  ";
		for( i = 0; i < nTax; i++ )
			qsf << DatumToChar( Matrix( i, j ) );
		if( j == top )
			qsf << " <-- top   ";
		if( j == ii )
			qsf << " <-- i     ";
		if( j == bottom )
			qsf << " <-- bottom";
		if( j == jj )
			qsf << " <-- j     ";
		if( j == xx )
			qsf << " <-- x     ";
		qsf << endl;
	}

	qsf.close();
}

//
//  QSort implements the quicksort algorithm
//
void DataMatrix::QSort( int top, int bottom )
{
	int i = top;
	int j = bottom;
	int x = ( top + bottom ) / 2;
	//DebugSaveQSortState( top, bottom, i, j, x, "Entering QSort" );
	do {
		while( ComparePatterns( i, x ) < 0  &&  i < bottom ) i++ ;
		while( ComparePatterns( x, j ) < 0  &&  j > top    ) j-- ;

		if( i <= j ) {
			//DebugSaveQSortState( top, bottom, i, j, x, "Just about to swap i and j" );
			SwapCharacters( i, j );

			if( x == i )		// keep track of the reference pattern!
				x = j;
			else if( x == j )	
				x = i;
			i++;
			if(j) j--;
			//DebugSaveQSortState( top, bottom, i, j, x, "Just after swapping" );
		}

	} while( i <= j );

	if( top <  j      ) QSort( top, j      );
	if( i   <  bottom ) QSort(   i, bottom );
}

int DataMatrix::GetToken( istream& in, char* tokenbuf, int maxlen, bool acceptComments /*=true*/ )
{
	int ok = 1;

	int i;
	char ch = ' ';

	// skip leading whitespace
	while( in && ( isspace(ch) || ch == '[' ) ){
		in.get(ch);
		if(ch == '[' && acceptComments==false) return -1;
		}
	if( !in ) return 0;

	tokenbuf[0] = ch;
	tokenbuf[1] = '\0';
	tokenbuf[maxlen-1] = '\0';
		
	for( i = 1; i < maxlen-1; i++ ) {
		in.get(ch);
		if( isspace(ch) || ch == ']' )
			break;
		tokenbuf[i] = ch;
		tokenbuf[i+1] = '\0';
	}

	if( i >= maxlen-1 )
		ok = 0;

	return ok;
}

int DataMatrix::GetToken( FILE *in, char* tokenbuf, int maxlen){
	int ok = 1;

	int i;
	char ch = ' ';

	// skip leading whitespace
	while( !ferror(in) && ( isspace(ch) || ch == '[' ) ){
		ch = getc(in);
		}
	if( ferror(in) ) return 0;

	tokenbuf[0] = ch;
	tokenbuf[1] = '\0';
	tokenbuf[maxlen-1] = '\0';
		
	for( i = 1; i < maxlen-1; i++ ) {
		ch = getc(in);
		if( isspace(ch) || ch == ']' )
			break;
		tokenbuf[i] = ch;
		tokenbuf[i+1] = '\0';
	}

	if( i >= maxlen-1 )
		ok = 0;

	return ok;
}
//
// Read reads in data from a file

int DataMatrix::ReadPhylip( const char* infname){

	//PARTITION
	ModelSpecification *modSpec = modSpecSet.GetModSpec(0);

	char ch;
	bool isNexus=false;

	FILE *inf;
#ifdef BOINC
	char input_path[512];
	boinc_resolve_filename(infname, input_path, sizeof(input_path));
    inf = boinc_fopen(input_path, "r");
#else
	inf = fopen(infname, "r");
#endif
	if(ferror(inf)) throw ErrorException("problem opening datafile %s for reading", infname);
	
	// get comments (note: comments only allowed at the beginning of the file)
	int end_of_comments = 0;
	while( !end_of_comments ){
		ch = getc(inf);
		if( ch != '/' ) {
			ungetc(ch, inf);
			end_of_comments = 1;
			}
		else {
			// ch is a slash, ignore rest of this line
			while( ch != '\n' && ch != '\r' && ch != EOF) {
				ch = getc(inf);
				}
			}
		}

	// get the dimensions of the data file
	int num_taxa=0, num_chars=0;
	fscanf(inf, "%d  %d", &num_taxa, &num_chars);

	NewMatrix( num_taxa, num_chars );

	// read in the data, including taxon names
	int blockStartNum = 0, charNum, i;
	bool firstPass = true;
	bool allDataRead = false;
	bool interleaved = false;

	while(allDataRead == false){
		//loop over the taxa, doing so multiple times for interleaved data
		for( i = 0; i < num_taxa; i++ ) {
			if(firstPass){
				// get name for taxon i
				char taxon_name[ MAX_TAXON_LABEL ];
				int ok = GetToken(inf, taxon_name, MAX_TAXON_LABEL);
				if( !ok ) {
					throw ErrorException("problem reading data: name for taxon #%d too long", i+1);
					}
				SetTaxonLabel( i, taxon_name );
				}

			// get data for taxon i
			unsigned char datum;
			for( charNum = blockStartNum; charNum < num_chars; charNum++ ) {
				if(firstPass == false && charNum == blockStartNum){
					do{
						ch = getc(inf);
						}while(isspace(ch) && ch != EOF);
					}
				else{
					do{
						ch = getc(inf);
						}while(ch == ' ' || ch == '\t');
					}
				if(ch == '['){//if there is a comment here, which is how the "color" used to be represented
					while (ch != ']' && ch != EOF) ch = getc(inf);
					ch = getc(inf);
					}
				if( ch == '.' ){
	    			datum = Matrix( 0, charNum );
					}
				else if(ch == '\n' || ch == '\r'){
					//file must be interleaved (or broken)
					if(!interleaved && i != 0)
						throw ErrorException("Unexpected line break found while reading data for taxon %s", TaxonLabel(i));
					else{
						interleaved = true;
						break;
						}
					}
	 			else{ 
					if(modSpec->IsAminoAcid() && modSpec->IsCodonAminoAcid() == false)
						datum = CharToDatum(ch);
					else 
						datum = CharToBitwiseRepresentation(ch);
					}
				SetMatrix( i, charNum, datum );
				}
			}
		if(charNum == num_chars && i == num_taxa) allDataRead = true;
		else{
			firstPass = false;
			blockStartNum = charNum;
			}
		}

	// read in the line containing the counts
	do{
		ch = getc(inf);
		}while(ch != EOF && isspace(ch));
	if( !feof(inf) ) {
		if(isdigit(ch) == false) throw ErrorException("Found extraneous information at end of phylip formatted datafile");
		ungetc(ch, inf);
		int i;
		char buf[10];
		for( i = 0; i < num_chars; i++ ) {
			int ok = GetToken( inf, buf, 10);
			if(feof(inf)) break;
			int cnt = atoi(buf);			
			SetCount( i, cnt );
			}
		if(i != num_chars) throw ErrorException("problem reading pattern counts");
		else dense = 1;
		//DJZ 9-13-06
		//It is very important to properly set the numNonMissingRealSitesInOrigMatrix variable now
		//to be the sum of the counts, otherwise bootstrapping after reading
		//a .cond file will give wrong resampling!!!!!
		numNonMissingRealSitesInOrigMatrix=0;
		for(int i=0;i<num_chars;i++){
			numNonMissingRealSitesInOrigMatrix += count[i];
			}
		}

	// read in the line containing the number of states for each character
	if( ferror(inf) == false ) {
		int i;
		char buf[10];
		for( i = 0; i < num_chars; i++ ) {
			int nstates;
			GetToken(inf, buf, 10);
			if( !inf ) break;
			nstates = atoi(buf);
			SetNumStates( i, nstates );
		}
	}

	fclose(inf);
	return 1;
}

int DataMatrix::ReadFasta( const char* infname){
	char ch;
	bool isNexus=false;

	//PARTITION
	ModelSpecification *modSpec = modSpecSet.GetModSpec(0);

	FILE *inf;
#ifdef BOINC
	char input_path[512];
	boinc_resolve_filename(infname, input_path, sizeof(input_path));
    inf = boinc_fopen(input_path, "r");
#else
	inf = fopen(infname, "r");
#endif
	if(ferror(inf)) throw ErrorException("problem opening datafile %s for reading", infname);
	
	//we don't know in advance what the number of characters or taxa is with fasta files
	//So, look through the file to get that info, then read it for real
	
	int num_taxa=0, num_chars=0, cur_char;
	char taxon_name[ MAX_TAXON_LABEL ];
	while( !feof(inf) ){
		GetToken(inf, taxon_name, MAX_TAXON_LABEL);
		num_taxa++;
		cur_char = 0;
		do{
			ch = getc(inf);
			if( !isspace(ch) && ch != '>' && ch != EOF) cur_char++;
			}while(ch != '>' && ch != EOF);
		if(num_taxa == 1) num_chars = cur_char;
		else if(cur_char != num_chars) throw ErrorException("# of characters for taxon %s (%d) not equal\n\tto the # of characters for first taxon (%d)", taxon_name, cur_char, num_chars);
		}
	rewind(inf);

	NewMatrix( num_taxa, num_chars );

	for(int i = 0; i < num_taxa; i++ ) {
		// get name for taxon i
		char taxon_name[ MAX_TAXON_LABEL ];
		int ok = GetToken(inf, taxon_name, MAX_TAXON_LABEL);
		if( !ok ) {
			throw ErrorException("problem reading data: name for taxon #%d too long", i+1);
			}
		SetTaxonLabel( i, taxon_name+1 );
			

		// get data for taxon i
		unsigned char datum;
		for(int charNum = 0; charNum < num_chars; charNum++ ) {
			do{
				ch = getc(inf);
				}while(isspace(ch));
			if(modSpec->IsAminoAcid() && modSpec->IsCodonAminoAcid() == false)
				datum = CharToDatum(ch);
			else 
				datum = CharToBitwiseRepresentation(ch);

			SetMatrix( i, charNum, datum );
			}
		}

	fclose(inf);
	return 1;
}

void DataMatrix::DumpCounts( const char* s )
{
	ofstream tmpf( "tmpfile.txt", ios::out | ios::app );
   tmpf << endl << endl;
   if(s) { tmpf << s << endl; }
   for( int j = 0; j < numPatterns; j++ ) {
      tmpf << j << "  " << Count(j) << endl;
   }
   tmpf << endl;
}

//
// saves data under the name s but with extension changed to '.mlt'
// if third argument supplied, a NEXUS file ending in '.nex' is saved also
//
int DataMatrix::Save( const char* path, char* newfname /* = 0 */, char*
#if defined( AUTOSAVE_NEXUS )
	nxsfname /* = 0 */
#endif
   )
{
	int i, j;//, nchar_total;
	char newpath[ MAXPATH ];

	strcpy( newpath, path );

#if defined( AUTOSAVE_NEXUS )
	//  ________________________________________
	// |                                        |
	// | save uncompressed data to file nxspath |
	// |________________________________________|
	//
   int k;
	char nxspath[ MAXPATH ];
	strcat( nxspath, ".nex" );

	cerr << endl << "Opening file '" << nxspath << "' for saving..." << endl;

	ofstream nxsf( nxspath );
	if( !nxsf ) {
		cerr << endl << "Error: could not open file '" << nxspath << "' for saving" << endl;
		return 0;
	}

	nchar_total = 0;
	for( j = 0; j < numPatterns; j++ )
		nchar_total += Count(j);

	nxsf << "#nexus" << endl << endl;
	nxsf << "begin data;" << endl;
	nxsf << "  dimensions ntax=" << nTax << "  nchar=" << nchar_total << ";" << endl;
	nxsf << "  format missing=? datatype=standard;" << endl;
	nxsf << "  matrix" << endl;

	for( i = 0; i < nTax; i++ ) {
		nxsf << TaxonLabel(i) << "  ";
		nxsf << " [" << TaxonColor(i) << "]  ";

		for( j = 0; j < numPatterns; j++ ) {
			for( k = 0; k < Count(j); k++ ) {
				nxsf << DatumToChar( Matrix( i, j ) );
			}
		}
		nxsf << endl;
	}

	nxsf << ";" << endl;
	nxsf << "end;" << endl << endl;

	if( !nxsf ) {
		cerr << endl << "Error saving data to file '" << nxspath << "':  disk full?" << endl;
		return 0;
	}

	nxsf.close();

	if( nxsfname ) {
		strcpy( nxsfname, nxspath );
	}
#endif

	//  _______________________________________
	// |                                       |
	// | save compressed data to file newpath  |
	// |_______________________________________|
	//

	//strcat( newpath, ".comp" );
	outman.UserMessage("Opening file \"%s\" for saving...", newpath);

	ofstream outf( newpath );
	if( !outf ) throw ErrorException("Error: could not open file \"%s\"", newpath);

/*	nchar_total = 0;
	for( j = 0; j < numPatterns; j++ ) {
		int k = PatternType(j);
		if( (k & PT_CONSTANT) && !InvarCharsExpected() ) continue;
		nchar_total++;
	}
*/

	outf << "#NEXUS\nbegin data;\ndimensions ntax=" <<nTax << " nchar=" << numPatterns << ";\n";
	outf << "format datatype=dna missing=? gap=-;\n";
	outf << "matrix" << endl;

	//outf << nTax << "  " << numPatterns << endl;
	for( i = 0; i < nTax; i++ ) {
		outf << TaxonLabel(i) << "  ";
		for( j = 0; j < numPatterns; j++ ) {
//			int k = PatternType(j);
//			if( (k & PT_CONSTANT) && !InvarCharsExpected() ) continue;
			outf << DatumToChar( Matrix( i, j ) );
		}
		outf << endl;
	}

	outf << ";" << endl;
	outf << "end;";
	outf << "begin assumptions;\n";
	string str;
	this->MakeWeightSetString(str, "packed");
	outf << str.c_str() << "\n;end;\n";
	outf.close();
	return 1;

	// save a line containing the counts for each character
	for( j = 0; j < numPatterns; j++ ) {
//		int k = PatternType(j);
//		if( (k & PT_CONSTANT) && !InvarCharsExpected() ) continue;
		outf << Count(j) << ' ';
	}
	outf << endl;

	// save a line containing the number of states for each character
	for( j = 0; j < numPatterns; j++ ) {
//		int k = PatternType(j);
//		if( (k & PT_CONSTANT) && !InvarCharsExpected() ) continue;
		outf << NumStates(j) << ' ';
	}
	outf << endl;

	if( !outf ) {
		cerr << endl << "Error saving data to file '" << newpath << "':  disk full?" << endl;
		return 0;
	}

	outf.close();

	/* cjb
	if( newfname ) {
		strcpy( newfname, newpath );
	}
	*/

	return 1;
}

void DataMatrix::WriteCollapsedData(){
		
//write the data matrix
	for(int i=0;i<nTax;i++){
		}

	}

void DataMatrix::Reweight(FLOAT_TYPE prob){
	for(int i=0;i<numPatterns;i++){
		FLOAT_TYPE r=rnd.uniform();
		if(r * 2.0 < prob) count[i]++;
		else if(r < prob) count[i]--;
		}
	}

//4-15-08 adding the option to specify a resample proportion, for jackknifing and what Cecile Ane called
//the "multidimentional bootstrap"
int DataMatrix::BootstrapReweight(int seedToUse, FLOAT_TYPE resampleProportion){
	//a seed is passed in and used for the reweighting - Either for restarting or not
	//Either way we'll return the seed at the end of the reweighting, to be stored as the Population::nextBootstrapSeed
	//which allows exactly the same bootstraped datasets to be used in multiple runs, but with different
	//settings for the actual search
	if(resampleProportion >= 5.0) outman.UserMessage("WARNING: The resampleproportion setting is the proportion to resample,\nNOT the percentage (1.0 = 100%%).\nThe value you specified (%.2f) is a very large proportion.", resampleProportion);

	int originalSeed = rnd.seed();
	rnd.set_seed(seedToUse);

	//This is a little dumb, but since there are parallel counts and origCounts variables depending on whether the new PatternManager
	//is being used, need to alias them so that the remainder of this function works unchanged
	const int *origCountsAlias;
	if(newOrigCounts.size() > 0){
		origCountsAlias = &newOrigCounts[0];
		}
	else
		origCountsAlias = origCounts;

	int *countsAlias;
	if(newCount.size() > 0){
		countsAlias = &newCount[0];
		}
	else
		countsAlias = count;

	FLOAT_TYPE *cumProbs = new FLOAT_TYPE[numPatterns];
	
	FLOAT_TYPE p=0.0;
	cumProbs[0]=(FLOAT_TYPE) origCountsAlias[0] / ((FLOAT_TYPE) numNonMissingRealCountsInOrigMatrix);

	countsAlias[0] = 0;
	for(int i = 1;i < numPatterns;i++){
		cumProbs[i] = cumProbs[i-1] + (FLOAT_TYPE) origCountsAlias[i] / ((FLOAT_TYPE) numNonMissingRealCountsInOrigMatrix);
		countsAlias[i] = 0;
		}
	cumProbs[numPatterns - 1] = 1.0;

	//ofstream deb("counts.log", ios::app);
	//ofstream deb("counts.log");

	//round to nearest int
	int numToSample = (int) (((FLOAT_TYPE)numNonMissingRealCountsInOrigMatrix * resampleProportion) + 0.5);
	if(numToSample != numNonMissingRealCountsInOrigMatrix) outman.UserMessage("Resampling %d characters (%.2f%%).\n", numToSample, resampleProportion*100);

	for(int c=0;c<numToSample;c++){
		FLOAT_TYPE p=rnd.uniform();
		int pat=0; 
		while(p > cumProbs[pat]) 
			pat++;
		countsAlias[pat]++;
		}
/*
	for(int i = 0;i < numPatterns;i++)
		deb << i << "\t" << origCountsAlias[i] << "\t" << countsAlias[i] <<  endl;
*/
	//take a count of the number of chars that were actually resampled
	nonZeroCharCount = 0;
	int numZero = 0;
	int totCounts = 0;
	for(int d=0;d<numPatterns;d++){
		if(countsAlias[d] > 0) {
			nonZeroCharCount++;
			totCounts += countsAlias[d];
			}
		else 
			numZero++;
		}
	delete []cumProbs;
	assert(totCounts == numNonMissingRealCountsInOrigMatrix);
	assert(nonZeroCharCount + numZero == numPatterns);

	int nextSeed = rnd.seed();
	rnd.set_seed(originalSeed);
	return nextSeed;
	}

void DataMatrix::CheckForIdenticalTaxonNames(){
	const char *name1, *name2;
	vector< pair<int, int> > identicals;

	for(int t1=0;t1<nTax-1;t1++){
		for(int t2=t1+1;t2<nTax;t2++){
			name1 = TaxonLabel(t1);
			name2 = TaxonLabel(t2);
			if(_stricmp(name1, name2) == 0) identicals.push_back(make_pair(t1, t2));
			}
		}
	
	if(identicals.size() > 0){
		outman.UserMessage("Error! Multiple sequences with same name encountered!:");
		for(vector< pair<int, int> >::iterator it=identicals.begin() ; it != identicals.end() ; it++){
			outman.UserMessage("\t%s : numbers %d and %d", TaxonLabel((*it).first), (*it).first+1, (*it).second+1);
			}
		throw(ErrorException("Terminating.  Please make all sequence names unique!"));
		}
	}

void DataMatrix::GetStringOfOrigDataColumns(string &str) const{
	//note that GetSetAsNexusString takes zero offset indeces and converts them to
	//char nums, ie adds 1 to each
	NxsUnsignedSet chars;
	for(int c = numConditioningPatterns;c < numRealSitesInOrigMatrix + numConditioningPatterns;c++)
		chars.insert(origDataNumber[c]);
	str = NxsSetReader::GetSetAsNexusString(chars);
	}

void DataMatrix::CountMissingCharsByColumn(vector<int> &vec){
	for(int c = 0;c < numPatterns;c++){
		int missing = 0;
		for(int t = 0;t < nTax;t++){
			if(Matrix(t, c) == fullyAmbigChar) 
				missing++;
			}
		vec.push_back(missing);
		}
	}

void DataMatrix::MakeWeightSetString(NxsCharactersBlock &charblock, std::string &wtstring, string name){
	NxsTransformationManager &transformer = charblock.GetNxsTransformationManagerRef();
	//this is a list of IntWeightToIndexSet objects
	NxsTransformationManager::ListOfIntWeights intWeights;
	
	NxsUnsignedSet dummy;
	//the charset was empty, implying that all characters in this block will go into a single matrix
	for(int i = 0;i < charblock.GetNumChar();i++)
		dummy.insert(i);

	for(int countNum = 0;dummy.size() > 0;countNum++){
		//this is a pair<int, std::set<unsigned> >
		NxsTransformationManager::IntWeightToIndexSet weightToIndex;
		weightToIndex.first = countNum;
		for(NxsUnsignedSet::iterator it = dummy.begin();it != dummy.end();){
			int thisCount = Count(*it);
			if(thisCount == countNum){
				weightToIndex.second.insert(*it);
				int err = *it++;
				dummy.erase(err);
				}
			else it++;
			}
		if(weightToIndex.second.size() > 0)
			intWeights.push_back(weightToIndex);
		}
	
	transformer.AddIntWeightSet("bootstrapped", intWeights, true);
	ostringstream out;
	transformer.WriteWtSet(out);
	wtstring = out.str();
	}

void DataMatrix::MakeWeightSetString(std::string &wtstring, string name){
	NxsTransformationManager transformer;// = charblock.GetNxsTransformationManagerRef();
	//this is a list of IntWeightToIndexSet objects
	NxsTransformationManager::ListOfIntWeights intWeights;
	
	NxsUnsignedSet dummy;
	//the charset was empty, implying that all characters in this block will go into a single matrix
	for(int i = 0;i < numPatterns;i++)
		dummy.insert(i);

	for(int countNum = 0;dummy.size() > 0;countNum++){
		//this is a pair<int, std::set<unsigned> >
		NxsTransformationManager::IntWeightToIndexSet weightToIndex;
		weightToIndex.first = countNum;
		for(NxsUnsignedSet::iterator it = dummy.begin();it != dummy.end();){
			int thisCount = Count(*it);
			if(thisCount == countNum){
				weightToIndex.second.insert(*it);
				int err = *it++;
				dummy.erase(err);
				}
			else it++;
			}
		if(weightToIndex.second.size() > 0)
			intWeights.push_back(weightToIndex);
		}
	
	transformer.AddIntWeightSet(name.c_str(), intWeights, true);
	ostringstream out;
	transformer.WriteWtSet(out);
	wtstring = out.str();
	}

