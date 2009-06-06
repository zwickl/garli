// GARLI version 0.96b8 source code
// Copyright 2005-2008 Derrick J. Zwickl
// email: zwickl@nescent.org
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

#include "defs.h"
#include "bipartition.h"
#include "reconnode.h"
#include "tree.h"

int Bipartition::nBlocks;
int Bipartition::blockBits;
int Bipartition::ntax;
unsigned int Bipartition::largestBlockDigit;
unsigned int Bipartition::allBitsOn;
unsigned int Bipartition::partialBlockMask;

Bipartition::Bipartition(const char *c) {
	//construct from a ***.... string
	Clear();
	size_t len= (size_t) strlen(c);
	assert(len == (size_t)ntax);
	Bipartition tmp;
	for(unsigned i=0;i<len;i++){
		if(c[i] == '*') *this += tmp.TerminalBipart(i+1);
		}
	Standardize();
}

std::string Bipartition::Output() const  {
	std::string s;
	s.reserve(nBlocks*blockBits);
	for (int i = 0; i < nBlocks; i++) {
		unsigned int t = rep[i];
		for(int j=0;j<blockBits;j++) {
			if(i*blockBits+j >= ntax)
				break;
			const char c = (t&largestBlockDigit ? '*' : '-');
			s.append(1,c);
			t = t << 1;
		}
	}
	return s;
}

void Bipartition::NumericalOutput(string &out, const Bipartition *mask) const {
		string left = "(";
		string right = "(";
		char temp[100];
		if(mask != NULL){
			for(int i=0;i<nBlocks;i++){
				unsigned int t=rep[i];
				unsigned int m=mask->rep[i];
				unsigned int bit = largestBlockDigit;
				for(int j=0;j<blockBits;j++){
					if(i*blockBits+j >= ntax) break;
					if(bit & m){
						if(bit & t){
							if(strcmp(left.c_str(), "(")) left += ',';
							sprintf(temp, "%d", (i*blockBits + j + 1));
							left += temp;
							}
						else{
							if(strcmp(right.c_str(), "(")) right += ',';
							sprintf(temp, "%d", (i*blockBits + j + 1));
							right += temp;
							}
						}
					bit = bit >> 1;	
					}
				}
			}
		else{
			for(int i=0;i<nBlocks;i++){
				unsigned int t=rep[i];
				unsigned int bit = largestBlockDigit;
				for(int j=0;j<blockBits;j++){
					if(i*blockBits+j >= ntax) break;
					if(bit & t){
						if(strcmp(left.c_str(), "(")) left += ',';
						sprintf(temp, "%d", (i*blockBits + j + 1));
						left += temp;
						}
					else{
						if(strcmp(right.c_str(), "(")) right += ',';
						sprintf(temp, "%d", (i*blockBits + j + 1));
						right += temp;
						}
					bit = bit >> 1;	 
					}
				}
			}
		left += ')';
		right += ')';
		out = left + " | " + right;
		}


//note that this is "less than" for sorting purposes, not in a subset sense
bool BipartitionLessThan(const Bipartition &lhs, const Bipartition &rhs){
	int i;
	for(i=0;i<Bipartition::nBlocks-1;i++){
		if(lhs.rep[i] > rhs.rep[i]) return false;
		else if(lhs.rep[i] < rhs.rep[i]) return true;
		}
		
	if(((lhs.rep[i]) & lhs.partialBlockMask) > ((rhs.rep[i]) & lhs.partialBlockMask)) return false;
	else return true;
	}

void Bipartition::SetBipartitionStatics(int nt){
	Bipartition::blockBits=sizeof(int)*8;
	Bipartition::ntax=nt;
	Bipartition::nBlocks=(int)ceil((FLOAT_TYPE)nt/(FLOAT_TYPE)Bipartition::blockBits);
	Bipartition::largestBlockDigit=1<<(Bipartition::blockBits-1);
	Bipartition::allBitsOn=(unsigned int)(pow(2.0, Bipartition::blockBits)-1);
	Bipartition::SetPartialBlockMask();	
	}

void Bipartition::SetPartialBlockMask(){
		partialBlockMask=0;
		unsigned int bit=largestBlockDigit;
		for(int b=0;b<ntax%blockBits;b++){
			partialBlockMask += bit;
			bit = bit >> 1;
			}
		if(ntax%blockBits == 0) partialBlockMask=allBitsOn;
		}

//this function does 2 things
//1. Fills this bipartition with the bitwise intersection of a backbone mask and a mask
//representing a subset of taxa in a growing tree.  Note that it is safe to call this when the
//constraint is not a backbone and/or when the partialMask is NULL.  In that case it will fill
//the bipartition with one or the other, or with all bits on if their if neither 
//2. Checks if there is a meaningful intersection between the created joint mask and 
//the constraint.  That means at least 2 bits are "on" on each site of the constrained bipartition
bool Bipartition::MakeJointMask(const Constraint &constr, const Bipartition *partialMask){
	if(constr.IsBackbone()){
		//this just uses Bipartition::Operator=()
		*this = constr.GetBackboneMask();
		if(partialMask != NULL)
			this->AndEquals(*partialMask);
		}
	else if(partialMask != NULL){
		*this = *(partialMask);
		}
	else FillAllBits();

	Bipartition temp;
	temp = constr.GetBipartition();
	temp.AndEquals(*this);
	if(temp.CountOnBits() < 2) return false;
	temp = constr.GetBipartition();
	temp.Complement();
	temp.AndEquals(*this);
	if(temp.CountOnBits() < 2) return false;
	return true;
	}

//check that any defined constraints are present in the tree.  If, so return false
//	and set violated to a description of the violated split. 
bool Tree::ObeysConstraints(std::string *violated) const {	
	vector<Constraint>::const_iterator conIt = Tree::constraintsVec.begin();
 	for(; conIt != Tree::constraintsVec.end(); conIt++) {
 		if (! conIt->CheckTree(*this)) {
 			if (violated)
 				*violated = conIt->GetBipartition().Output();
			return false;
 		}
 	}
 	return true;
}
 
bool Constraint::CheckTree(const Tree &tree) const {
	const TreeNode *check = NULL;
	if (this->IsBackbone())
		check = tree.ContainsMaskedBipartitionOrComplement(this->GetBipartition(), this->GetBackboneMask());
	else
		check = tree.ContainsBipartitionOrComplement(this->GetBipartition());
	if (this->IsPositive())
		return check != 0L;
	return check == 0L;
}

void Tree::AddConstraint(const Constraint & constraint) {
	constraintsVec.push_back(constraint);
}


void Tree::LoadConstraints(ifstream &con, int nTaxa) {
	string temp;//=new char[numTipsTotal + 100];
	Constraint constr;
	do {
		temp.clear();
		char c;
		con.get(c);
		do {
			temp += c;
			con.get(c);
		}
		while(c != '\n' && c!= '\r' && con.eof() == false);
		
		while((con.peek() == '\n' || con.peek() == '\r') && con.eof() == false) {
			con.get(c);
		}

		//getline works strangely on some compilers.  temp should end with ; or \0 , but
		//might end with \r or \n
		size_t len=temp.length();
		char last=temp.c_str()[len-1];
		while(last == '\r' || last == '\n' || last == ' ') {
			temp.erase(len-1, 1);
			len--;
			last=temp.c_str()[len-1];
		}
		if(temp[0] != '\0') {
			if(temp[0] != '+' && temp[0] != '-')
				throw ErrorException("constraint string must start with \'+\' (positive constraint) or \'-\' (negative constraint)");
			if(temp[1] == '.' || temp[1] == '*') {//if individual biparts are specified in *. format
				//while(temp[temp.length()-1] == ' ') temp.erase(temp.length()-1);//eat any spaces at the end
				if(len != nTaxa+1) {
					const char * msg = "constraint # %d does not have the correct number of characters!\n(has %d) constraint strings must start with \n\'+\' (positive constraint) or \'-\' (negative constraint)\nfollowed by either a ...*** type specification\nor a constraint in newick format.  \nNote that backbone constraints cannot be specified in ...*** format.";
					throw ErrorException(msg, Tree::GetConstraints().size(), len);
				}
				constr.ReadDotStarConstraint(temp.c_str());
				Tree::AddConstraint(constr);
			}
			else if(temp[1] == '(') {//if a constraint tree in parenthetical notation is used
				bool numericalTaxa=true;
				for(unsigned i=0;i<len;i++) {//see if we are dealing with a treestring with taxa as # or names
					if(isalpha(temp[i])) {
						numericalTaxa=false;
						break;
					}
				}
				const bool isPositive = (temp[0] == '+');
				ReadNewickConstraint(temp.c_str() + 1, numericalTaxa, isPositive);
			}
			else {
				const char * msg = "problem with constraint # %d\nconstraint strings must start with \n\'+\' (positive constraint) or \'-\' (negative constraint)\nfollowed by either a ...*** type specification\nor a constraint in newick format";
				throw ErrorException(msg, GetConstraints().size(), len);
			}
		}
	}
	while(!con.eof());


	//make sure the constraints are compatible with each other!
	if(Tree::IsUsingConstraints()) {
		for(vector<Constraint>::const_iterator first=GetConstraints().begin(); first != GetConstraints().end();first++) {
			for(vector<Constraint>::const_iterator sec=first+1;sec!=GetConstraints().end();sec++) {
				if((*first).IsPositive() != (*sec).IsPositive())
					throw ErrorException("cannot mix positive and negative constraints!");
				if(((*first).IsPositive()==false) && ((*sec).IsPositive()==false))
					throw ErrorException("Sorry, GARLI can currently only handle a single negatively (conversely) constrainted branch :-(");
				if((*first).ConstraintIsCompatibleWithConstraint((*sec)) == false)
					throw ErrorException("constraints are not compatible with one another!");
			}
		}
	}
	//summarize the constraint info to the screen
	string str;
	int num=1;
	if(GetConstraints()[0].IsPositive()) {
		outman.UserMessage("Found %d positively constrained bipartition(s)", GetConstraints().size());
		for(vector<Constraint>::const_iterator first=GetConstraints().begin();first!=GetConstraints().end();first++) {
			(*first).NumericalOutput(str);
			if((*first).IsBackbone())
				outman.UserMessage("     Bipartition %d (backbone): %s", num, str.c_str());
			else
				outman.UserMessage("     Bipartition %d: %s", num, str.c_str());
			num++;
		}
	}
	else {
		outman.UserMessage("Found 1 negatively (conversely) constrained bipartition");
		Tree::GetConstraints()[0].NumericalOutput(str);
		if(Tree::GetConstraints()[0].IsBackbone())
			outman.UserMessage("     Bipartition %d (backbone): %s", num, str.c_str());
		else
			outman.UserMessage("     Bipartition %d: %s", num, str.c_str());
	}
}
