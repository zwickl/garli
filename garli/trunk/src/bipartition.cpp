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

#ifdef BETTER_BIPART
unsigned Bipartition::nBlocks;
unsigned Bipartition:: blockBits;
unsigned Bipartition::ntax;
#else
int Bipartition::nBlocks;
int Bipartition:: blockBits;
int Bipartition::ntax;
#endif

//class statics
unsigned int Bipartition::largestBlockDigit;
unsigned int Bipartition::allBitsOn;
char * Bipartition::str;
unsigned int Bipartition::partialBlockMask;

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
	Bipartition::allBitsOn=(unsigned int)(pow(2.0, (int) Bipartition::blockBits)-1);
	Bipartition::str=new char[nt+1];
	Bipartition::str[nt] = '\0';
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
		*this = *(constr.GetBackboneMask());
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


void Bipartition::ClearBipartition(){
#ifdef BETTER_BIPART
	memset(rep, 0L, sizeof(int) * nBlocks);
#else
	for(int i=0;i<nBlocks;i++)
	rep[i]=0;
#endif
	}

void Bipartition::operator+=(const Bipartition *rhs){
	for(int i=0;i<nBlocks;i++){
		rep[i]|=rhs->rep[i];
		}
	}

void Bipartition::operator-=(const Bipartition *rhs){
	//note that this assumes that rhs is a subset of this or vice versa!!!
	if(rhs->IsASubsetOf(*this)){
		for(int i=0;i<nBlocks;i++){
			rep[i] ^= rhs->rep[i];
			}
		}
	else{
		for(int i=0;i<nBlocks;i++){
			rep[i] |= ~rhs->rep[i];
			}
		}
	}

void Bipartition::operator=(const Bipartition *rhs){
	memcpy(rep, rhs->rep, nBlocks*sizeof(int));
	}

void Bipartition::operator=(const Bipartition &rhs){
	memcpy(rep, rhs.rep, nBlocks*sizeof(int));
	}

int Bipartition::CountOnBits() const{
	int num=0;
	for(int i=1;i<=ntax;i++)
		if(ContainsTaxon(i)) num++;
	return num;
	}

void Bipartition::FlipBits(const Bipartition &flip){
	//this just flips the bits in the passed in bipartition
	EqualsXORComplement(flip);
	}

void Bipartition::FillAllBits(){
#ifdef BETTER_BIPART
	memset(rep, allBitsOn, sizeof(int) * nBlocks);
#else
	for(int i=0;i<nBlocks;i++)
		rep[i] = allBitsOn;
#endif
	}

bool Bipartition::EqualsEquals(const Bipartition &rhs) const{
	//assert(this->ContainsTaxon(1));
	//assert(rhs.ContainsTaxon(1));
	int i;
	for(i=0;i<nBlocks-1;i++){
		if(rep[i]!=rhs.rep[i]) return false;
		}
	if((rep[i]&partialBlockMask)!=((rhs.rep[i])&partialBlockMask)) return false;
	return true;
	}

bool Bipartition::MaskedEqualsEquals(const Bipartition &rhs, const Bipartition &mask) const{
	//assert(this->ContainsTaxon(1));
	//assert(rhs.ContainsTaxon(1));
	int i;
	for(i=0;i<nBlocks-1;i++){
		if((rep[i] & mask.rep[i]) != (rhs.rep[i] & mask.rep[i])) return false;
		}
	if((rep[i] & partialBlockMask & mask.rep[i]) != ((rhs.rep[i]) & partialBlockMask & mask.rep[i])) return false;
	return true;
	}

void Bipartition::Complement(){
	for(int i=0;i<nBlocks;i++){
		rep[i] = ~rep[i];
		}
	}

bool Bipartition::ComplementEqualsEquals(const Bipartition &rhs) const{
	int i;
	for(i=0;i<nBlocks-1;i++){
		if(~rep[i]!=rhs.rep[i]) return false;
		}
		
	if(((~rep[i])&partialBlockMask)!=(rhs.rep[i]&partialBlockMask)) return false;
	
	return true;
	}
	
bool Bipartition::MaskedComplementEqualsEquals(const Bipartition &rhs, const Bipartition &mask) const{
	//assert(this->ContainsTaxon(1));
	//assert(rhs.ContainsTaxon(1));
	int i;
	for(i=0;i<nBlocks-1;i++){
		if((~rep[i] & mask.rep[i]) != (rhs.rep[i] & mask.rep[i])) return false;
		}
	if((~rep[i] & partialBlockMask & mask.rep[i]) != ((rhs.rep[i]) & partialBlockMask & mask.rep[i])) return false;
	return true;
	}

void Bipartition::Standardize(){
#ifdef BETTER_BIPART
	if(ContainsFirstTaxon()==true) return;
#else
	if(ContainsTaxon(1)==true) return;
#endif
	else{
		for(int i=0;i<nBlocks;i++){
			rep[i]=~rep[i];
			}
		}
	}

int Bipartition::FirstPresentTaxon() const{
	int blk=0;
	unsigned int tmp=rep[blk];
	while(tmp==0) tmp=rep[++blk];
	
	int t=blockBits*blk+1;
	while( ! (tmp & largestBlockDigit)){
		tmp = tmp << 1;
		t++;
		}
	return t;
	}

int Bipartition::FirstNonPresentTaxon() const{
	int blk=0;
	unsigned int tmp=rep[blk];
	while(tmp==allBitsOn) tmp=rep[++blk];
	
	int t=blockBits*blk+1;
	while((tmp & largestBlockDigit)){
		tmp = tmp << 1;
		t++;
		}
	return t;
	}

#ifdef BETTER_BIPART
//if you know the taxon is in the first block, it is uneccesary and very slow to do the divisions
//in ContainsTaxon
bool Bipartition::ContainsFirstTaxon() const{
	if(rep[0] & largestBlockDigit) 
		return true;
	return false;
	}

bool Bipartition::ContainsTaxon(int t) const{
	int index = 0;
	while(((index+1) * blockBits) < t)
		index++;
	assert(index < nBlocks);
	//unsigned int tmp=rep[(t-1)/blockBits];
	const unsigned int tmp=rep[index];
	if(tmp & largestBlockDigit >> ((t - 1) - (index * blockBits))) return true;
	//if(tmp & largestBlockDigit>>((t-1)%blockBits)) return true;
	return false;
	}
#else
bool Bipartition::ContainsTaxon(int t) const{
	unsigned int tmp=rep[(t-1)/blockBits];
	if(tmp & largestBlockDigit>>((t-1)%blockBits)) return true;
	return false;
	}
#endif

void Bipartition::FillWithXORComplement(const Bipartition &lhs, const Bipartition &rhs){
	for(int i=0;i<nBlocks;i++){
		rep[i] = lhs.rep[i] ^ (~rhs.rep[i]);
		}
	}

void Bipartition::EqualsXORComplement(const Bipartition &rhs){
	for(int i=0;i<nBlocks;i++){
		rep[i] = rep[i] ^ (~rhs.rep[i]);
		}
	}

void Bipartition::AndEquals(const Bipartition &rhs){
	//sets the calling bipart to be the bitwise AND of it and the argument
	for(int i=0;i<nBlocks;i++){
		rep[i] = rep[i] & rhs.rep[i];
		}
	}

bool Bipartition::IsCompatibleWithBipartition(const Bipartition &constr) const{
	//using buneman's 4 point condition.  At least one of the four intersections must be empty
	bool compat=true;
	int i;
	//A & B
	for(i=0;i<nBlocks-1;i++){
		if((rep[i] & constr.rep[i]) == 0) continue;
		else{
			compat=false;
			break;
			}
		}
	if(compat==true && (((rep[i] & constr.rep[i]) & partialBlockMask) == 0)) return true;
	//A & ~B
	compat=true;
	for(i=0;i<nBlocks-1;i++){
		if((rep[i] & ~(constr.rep[i])) == 0) continue;
		else{
			compat=false;
			break;
			}
		}
	if(compat==true && (((rep[i] & ~(constr.rep[i])) & partialBlockMask) == 0)) return true;
	//~A & B
	compat=true;
	for(i=0;i<nBlocks-1;i++){
		if((~rep[i] & constr.rep[i]) == 0) continue;
		else{
			compat=false;
			break;
			}
		}
	if(compat==true && (((~rep[i] & constr.rep[i]) & partialBlockMask) == 0)) return true;
	//~A & ~B
	compat=true;
	for(i=0;i<nBlocks-1;i++){
		if((~rep[i] & ~(constr.rep[i])) == 0) continue;
		else{
			compat=false;
			break;
			}
		}
	if(compat==true && (((~rep[i] & ~(constr.rep[i])) & partialBlockMask) == 0)) return true;
	return false;
	}

bool Bipartition::IsCompatibleWithBipartitionWithMask(const Bipartition &bip, const Bipartition &mask) const{
	//using buneman's 4 point condition.  At least one of the four intersections must be empty
	//typically this will be called with 'this' as a POSITIVE constraint
	bool compat=true;
	int i;

	//A & B
	for(i=0;i<nBlocks-1;i++){
		if(((rep[i] & bip.rep[i]) & mask.rep[i]) == 0) continue;
		else{
			compat=false;
			break;
			}
		}
	if(compat==true && ((((rep[i] & bip.rep[i]) & mask.rep[i]) & partialBlockMask) == 0)) return true;
	//A & ~B
	compat=true;
	for(i=0;i<nBlocks-1;i++){
		if(((rep[i] & ~bip.rep[i]) & mask.rep[i]) == 0) continue;
		else{
			compat=false;
			break;
			}
		}
	if(compat==true && ((((rep[i] & ~bip.rep[i]) & mask.rep[i]) & partialBlockMask) == 0)) return true;
	//~A & B
	compat=true;
	for(i=0;i<nBlocks-1;i++){
		if(((~rep[i] & bip.rep[i]) & mask.rep[i]) == 0) continue;
		else{
			compat=false;
			break;
			}
		}
	if(compat==true && ((((~rep[i] & bip.rep[i]) & mask.rep[i]) & partialBlockMask) == 0)) return true;
	//~A & ~B
	compat=true;
	for(i=0;i<nBlocks-1;i++){
		if(((~rep[i] & ~bip.rep[i]) & mask.rep[i]) == 0) continue;
		else{
			compat=false;
			break;
			}
		}
	if(compat==true && ((((~rep[i] & ~bip.rep[i]) & mask.rep[i]) & partialBlockMask) == 0)) return true;
	return false;
	}

bool Bipartition::IsCompatibleWithNegativeBipartition(const Bipartition &bip) const{
	//To be consistent with a negative bipartition (constraint) neither the bipartition
	//or its complement can BE the constraint bipartition
	if(this->EqualsEquals(bip)==false && this->ComplementEqualsEquals(bip)==false) return true;
	else return false;
	}

bool Bipartition::IsCompatibleWithNegativeBipartitionWithMask(const Bipartition &bip, const Bipartition &mask) const{
	//To be consistent with a negative bipartition (constraint) neither the masked bipartition
	//or its masked complement can BE the constraint bipartition
	if(this->MaskedEqualsEquals(bip, mask)==false && this->MaskedComplementEqualsEquals(bip, mask)==false) return true;
	else return false;
	}

bool Bipartition::IsIncompatibleWithBipartitionWithMask(const Bipartition &bip, const Bipartition &mask) const{
	//using buneman's 4 point condition.  At none of the four intersections must be empty
	//typically this will be called with 'this' as a NEGATIVE constraint
	//2/28/08 - I think that I had the true/falses backwards here.  As soon as an
	//empty intersection is found the bipartition IS compatable with the constraint
	//so it should return FALSE.  If there are 4 intersections, it is INcompatable, so return TRUE

	bool compat=true;
	int i;

	//first check if there is no intersection between the mask and the constraint
	//or the complement of the constraint. If not, the question of incompatability
	//is moot (but it is consistent with incompatability, so return true)
	bool intersect=false;
	for(i=0;i<nBlocks-1;i++){
		if((rep[i] & mask.rep[i]) == 0) continue;
		else{
			intersect=true;
			break;
			}
		}
	if(intersect==false && ((((rep[i] & mask.rep[i]) & partialBlockMask) == 0))) return true;

	intersect=false;
	for(i=0;i<nBlocks-1;i++){
		if((~rep[i] & mask.rep[i]) == 0) continue;
		else{
			intersect=true;
			break;
			}
		}
	if(intersect==false && ((((~rep[i] & mask.rep[i]) & partialBlockMask) == 0))) return true;

	//A & B
	for(i=0;i<nBlocks-1;i++){
		if(((rep[i] & mask.rep[i]) & (bip.rep[i] & mask.rep[i])) == 0) continue;
		else{
			compat=false;
			break;
			}
		}
	if(compat==true && ((((rep[i] & mask.rep[i]) & (bip.rep[i] & mask.rep[i])) & partialBlockMask) == 0)) return false;
	//A & ~B
	compat=true;
	for(i=0;i<nBlocks-1;i++){
		if(((rep[i] & mask.rep[i]) & (~bip.rep[i] & mask.rep[i])) == 0) continue;
		else{
			compat=false;
			break;
			}
		}
	if(compat==true && ((((rep[i] & mask.rep[i]) & (~bip.rep[i]& mask.rep[i])) & partialBlockMask) == 0)) return false;
	//~A & B
	compat=true;
	for(i=0;i<nBlocks-1;i++){
		if(((~rep[i] & mask.rep[i]) & (bip.rep[i] & mask.rep[i])) == 0) continue;
		else{
			compat=false;
			break;
			}
		}
	if(compat==true && ((((~rep[i] & mask.rep[i]) & (bip.rep[i] & mask.rep[i])) & partialBlockMask) == 0)) return false;
	//~A & ~B
	compat=true;
	for(i=0;i<nBlocks-1;i++){
		if(((~rep[i] & mask.rep[i]) & (~bip.rep[i] & mask.rep[i])) == 0) continue;
		else{
			compat=false;
			break;
			}
		}
	if(compat==true && ((((~rep[i] & mask.rep[i]) & (~bip.rep[i] & mask.rep[i])) & partialBlockMask) == 0)) return false;
	return true;
	}

bool Bipartition::IsASubsetOf(const Bipartition &target) const{
	//returns true if this is a subset of target
	int i;
	for(i=0;i<nBlocks-1;i++){
		if((target.rep[i] | rep[i]) != target.rep[i]) return false;
		}
	if(((target.rep[i] | rep[i]) & partialBlockMask) != (target.rep[i] & partialBlockMask)) return false;
	return true;
	}

bool Bipartition::MaskedIsASubsetOf(const Bipartition &target, const Bipartition &mask) const{
	//returns true if this is a subset of target
	int i;
	for(i=0;i<nBlocks-1;i++){
		if(((target.rep[i] | rep[i]) & mask.rep[i]) != (target.rep[i] & mask.rep[i])) return false;
		}
	if(((target.rep[i] | rep[i]) & partialBlockMask & mask.rep[i]) != (target.rep[i] & partialBlockMask & mask.rep[i])) return false;
	return true;
	}

bool Bipartition::ComplementIsASubsetOf(const Bipartition &target) const{
	//returns true if this is NOT a subset of target's complement
	int i;
	for(i=0;i<nBlocks-1;i++){
		if((target.rep[i] | ~rep[i]) != ~rep[i]) return true;
		}
	if(((target.rep[i] | ~rep[i]) & partialBlockMask) != (~rep[i] & partialBlockMask)) return false;
	return true;
	}

Bipartition * Bipartition::TerminalBipart(int taxNum){
	assert(taxNum > 0 && taxNum <= ntax);
	for(int i=0;i<nBlocks;i++) rep[i]=0;
	//8-19-05 fixed this bit of stupidity
	//rep[(taxNum)/blockBits]|=(largestBlockDigit>>((taxNum-1)%blockBits));
	rep[(taxNum-1)/blockBits]|=(largestBlockDigit>>((taxNum-1)%blockBits));
	return this;
	}
char * Bipartition::Output(){
	for(int i=0;i<nBlocks;i++){
		unsigned int t=rep[i];
		for(int j=0;j<blockBits;j++){
			if(i*blockBits+j >= ntax) break;
			if(t&largestBlockDigit) str[i*blockBits+j]='*';
			else str[i*blockBits+j]='-';
			t=t<<1;
			}
		}
	return str;
	}

void Bipartition::NumericalOutput(string &out, const Bipartition *mask){
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
/*
void BinaryOutput(ofstream &out){
	int size = nBlocks * sizeof(unsigned int);
	out.write((char*) rep, size);
	}
*/
void Bipartition::BinaryOutput(OUTPUT_CLASS &out){
	out.WRITE_TO_FILE(rep, sizeof(unsigned int), nBlocks);
	}

void Bipartition::BinaryInput(FILE* &in){
	fread((char*) rep, sizeof(unsigned int), nBlocks, in);
	}

vector<int> Bipartition::NodenumsFromBipart(){
	vector<int> nodes;
	for(int i=1;i<=ntax;i++) if(ContainsTaxon(i)) nodes.push_back(i);
	return nodes; 
	}

void Bipartition::BipartFromNodenums(const vector<int> & nodes){
	ClearBipartition();
	Bipartition temp;
	for(vector<int>::const_iterator it = nodes.begin();it != nodes.end();it++)
		*this += temp.TerminalBipart(*it);
	}

void Bipartition::BipartFromNodenums(const std::set<unsigned> & nodes){
	ClearBipartition();
	Bipartition temp;
	for(set<unsigned>::const_iterator it = nodes.begin();it != nodes.end();it++)
		*this += temp.TerminalBipart(*it);
	}