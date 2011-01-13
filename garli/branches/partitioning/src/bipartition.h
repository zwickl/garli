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


#ifndef BIPARTITION
#define BIPARTITION

#include <cmath>
#include <string>
#include <cstring>
#include <vector>
#include <cassert>
#include <set>

using namespace std;

#include "errorexception.h"

class Constraint;

class Bipartition{
	public:	
	unsigned int *rep;
	static int nBlocks;
	static int blockBits;
	static int ntax;
	static unsigned int largestBlockDigit;
	static unsigned int allBitsOn;
	static char* str;
	static unsigned int partialBlockMask;//this can be used to mask out the bits that
										// aren't used in the last block.  This becomes
										//important if we start doing complements. Bits
										//that represent actual taxa are ON

	Bipartition(){
		rep=new unsigned int[nBlocks];
		ClearBipartition();
		}
	Bipartition(const Bipartition &b){//copy constructor
		rep=new unsigned int[nBlocks];
		memcpy(rep, b.rep, nBlocks*sizeof(int));
		//for(int i=0;i<nBlocks;i++) rep[i] = b.rep[i];
		}

	Bipartition(const char *c){//construct from a ***.... string
		rep=new unsigned int[nBlocks];
		ClearBipartition();
		size_t len=strlen(c);
		assert(len == ntax);
		Bipartition tmp;
		for(unsigned i=0;i<len;i++){
			if(c[i] == '*') *this += tmp.TerminalBipart(i+1);
			}
		Standardize();
		}

	~Bipartition(){
		if(rep!=NULL) delete []rep;
		rep=NULL;
		}	

	bool IsCompatibleWithBipartition(const Bipartition &constr) const;
	bool OldIsCompatibleWithBipartition(const Bipartition &constr) const;
	
	bool IsCompatibleWithBipartitionWithMask(const Bipartition &constr, const Bipartition &mask) const;
	bool OldIsCompatibleWithBipartitionWithMask(const Bipartition &bip, const Bipartition &mask) const;

	//this is not being used
	bool IsIncompatibleWithBipartitionWithMask(const Bipartition &bip, const Bipartition &mask) const;

	bool MakeJointMask(const Constraint &constr, const Bipartition *partialMask);
	void NumericalOutput(string &out, const Bipartition *mask);
	static void SetBipartitionStatics(int);
	static void SetPartialBlockMask();
	
	void ClearBipartition(){
		memset(rep, 0L, sizeof(int) * nBlocks);
		}
	
	void operator+=(const Bipartition *rhs){
		for(int i=0;i<nBlocks;i++){
			rep[i] |= rhs->rep[i];
			}
		}

	void operator-=(const Bipartition *rhs){
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

	void operator=(const Bipartition *rhs){
		memcpy(rep, rhs->rep, nBlocks*sizeof(int));
		}

	void operator=(const Bipartition &rhs){
		memcpy(rep, rhs.rep, nBlocks*sizeof(int));
		}

	int CountOnBits() const{
		int num=0;
		for(int i=1;i<=ntax;i++)
			if(ContainsTaxon(i)) num++;
		return num;
		}

	bool MoreThanOneBitSet() const{
		int num = 0;
		int i = 0;
		for(i=0;i<nBlocks-1;i++){
			//check if there is more than one bit set in this block
			if(rep[i] & (rep[i] - 1))
				return true;
			//keep track of how many bits we've seen on
			if(rep[i])
				num++;
			if(num > 1)
				return true;
			}
		if((rep[i] & partialBlockMask) & ((rep[i] & partialBlockMask) - 1))
			return true;
		if(rep[i] & partialBlockMask)
			num++;
		if(num > 1)
			return true;
		return false;
		}

	void FlipBits(const Bipartition &flip){
		//this just flips the bits in the passed in bipartition
		EqualsXORComplement(flip);
		}

	void FillAllBits(){
		//the argument here is what to fill each _byte_ of the ints with
		memset(rep, 0xFF, sizeof(int) * nBlocks);
		}

	bool EqualsEquals(const Bipartition &rhs) const{
		//assert(this->ContainsTaxon(1));
		//assert(rhs.ContainsTaxon(1));
		int i;
		for(i=0;i<nBlocks-1;i++){
			if(rep[i] != rhs.rep[i]) 
				return false;
			}
		if((rep[i]&partialBlockMask)!=((rhs.rep[i])&partialBlockMask)) 
			return false;
		return true;
		}

	bool MaskedEqualsEquals(const Bipartition &rhs, const Bipartition &mask) const{
		//assert(this->ContainsTaxon(1));
		//assert(rhs.ContainsTaxon(1));
		int i;
		for(i=0;i<nBlocks-1;i++){
			if((rep[i] & mask.rep[i]) != (rhs.rep[i] & mask.rep[i])) 
				return false;
			}
		if((rep[i] & partialBlockMask & mask.rep[i]) != ((rhs.rep[i]) & partialBlockMask & mask.rep[i])) 
			return false;
		return true;
		}

	//this flips bits in the unused portion of the last block, but since it is
	//always masked it doesn't matter
	void Complement(){
		for(int i=0;i<nBlocks;i++){
			rep[i] = ~rep[i];
			}
		}

	bool ComplementEqualsEquals(const Bipartition &rhs) const{
		int i;
		for(i=0;i<nBlocks-1;i++){
			if(~rep[i]!=rhs.rep[i]) 
				return false;
			}
			
		if(((~rep[i])&partialBlockMask)!=(rhs.rep[i]&partialBlockMask)) 
			return false;
		
		return true;
		}
	
	bool MaskedComplementEqualsEquals(const Bipartition &rhs, const Bipartition &mask) const{
		//assert(this->ContainsTaxon(1));
		//assert(rhs.ContainsTaxon(1));
		int i;
		for(i=0;i<nBlocks-1;i++){
			if((~rep[i] & mask.rep[i]) != (rhs.rep[i] & mask.rep[i])) 
				return false;
			}
		if((~rep[i] & partialBlockMask & mask.rep[i]) != ((rhs.rep[i]) & partialBlockMask & mask.rep[i])) 
			return false;
		return true;
		}

	void Standardize(){
		if(ContainsFirstTaxon())
			return;
		else
			Complement();
		}
	
	int FirstPresentTaxon() const{
		int blk=0;
		unsigned int tmp=rep[blk];
		while(tmp == 0) 
			tmp=rep[++blk];
		
		int t=blockBits*blk+1;
		while( ! (tmp & largestBlockDigit)){
			tmp = tmp << 1;
			t++;
			}
		return t;
		}

	int FirstNonPresentTaxon() const{
		int blk=0;
		unsigned int tmp=rep[blk];
		while(tmp==allBitsOn) 
			tmp=rep[++blk];
		
		int t=blockBits*blk+1;
		while((tmp & largestBlockDigit)){
			tmp = tmp << 1;
			t++;
			}
		return t;
		}

	bool ContainsTaxon(int t) const{
		unsigned int tmp=rep[(t-1)/blockBits];
		if(tmp & largestBlockDigit>>((t-1)%blockBits)) 
			return true;
		return false;
		}

	//if you know the taxon is in the first block, it is unnecessary and very slow to do 
	//the divisions in ContainsTaxon
	inline bool ContainsFirstTaxon() const{
		if(rep[0] & largestBlockDigit) 
			return true;
		return false;
		}

	void FillWithXORComplement(const Bipartition &lhs, const Bipartition &rhs){
		for(int i=0;i<nBlocks;i++){
			rep[i] = lhs.rep[i] ^ (~rhs.rep[i]);
			}
		}

	void EqualsXORComplement(const Bipartition &rhs){
		for(int i=0;i<nBlocks;i++){
			rep[i] = rep[i] ^ (~rhs.rep[i]);
			}
		}

	void AndEquals(const Bipartition &rhs){
		//sets the calling bipart to be the bitwise AND of it and the argument
		for(int i=0;i<nBlocks;i++){
			rep[i] = rep[i] & rhs.rep[i];
			}
		}

//NOTE that the _BETTER_BIPART sections here are usually (much) slower unless intersections are rare.
//They don't have the branch within the loop, but they do have an extra bitwise operation and don't
//allow breaking out early.  With branch prediction the branch cost is only paid upon early return
//(I think).  So, when the number of taxa is small and the loop executes few times it might be more
//helpful.  For nine blocks (329 taxa) with a mask it is much slower.  The story might be different
//without the mask, when the per loop cost is lower but the misprediction cost is the same.
bool HasIntersection(const Bipartition &rhs, const Bipartition *mask) const
	{
		int i;
		if(!mask){
		for(i=0;i<nBlocks-1;i++){
				if(rep[i] & rhs.rep[i]){
					return true;
					}
				}
			if(((rep[i] & rhs.rep[i]) & partialBlockMask))
				return true;
			}
#ifdef _BETTER_BIPART
			else{
			unsigned sum = 0;
			for(i=0;i<nBlocks-1;i++){
				sum |= ((rep[i] & rhs.rep[i]) & mask->rep[i]);
				}
			if(sum)
				return true;
			if((rep[i] & rhs.rep[i]) & (mask->rep[i] & partialBlockMask))
				return true;

			}
#else
		else{
		for(i=0;i<nBlocks-1;i++){
				if((rep[i] & rhs.rep[i]) & mask->rep[i]){
					return true;
				}
			}
			if(((rep[i] & rhs.rep[i]) & (mask->rep[i] & partialBlockMask)))
				return true;
			}
#endif
		return false;
		}

	bool HasIntersectionWithComplement(const Bipartition &rhs, const Bipartition *mask) const{
		int i;
		if(!mask){
		for(i=0;i<nBlocks-1;i++){
				if(rep[i] & ~rhs.rep[i]){
					return true;
					}
				}
			if(((rep[i] & ~rhs.rep[i]) & partialBlockMask))
				return true;
			}
#ifdef _BETTER_BIPART
			else{
			unsigned sum = 0;
			for(i=0;i<nBlocks-1;i++){
				sum |= ((rep[i] & ~rhs.rep[i]) & mask->rep[i]); 
				}
			if(sum)
				return true;
			if((rep[i] & ~rhs.rep[i]) & (mask->rep[i] & partialBlockMask))
				return true;
			}
#else
		else{
		for(i=0;i<nBlocks-1;i++){
				if((rep[i] & ~rhs.rep[i]) & mask->rep[i]){
					return true;
					}
				}
			if(((rep[i] & ~rhs.rep[i]) & (mask->rep[i] & partialBlockMask)))
				return true;
			}
#endif
		return false;
		}

	bool ComplementHasIntersection(const Bipartition &rhs, const Bipartition *mask) const{
		int i;
		if(!mask){
		for(i=0;i<nBlocks-1;i++){
				if(~rep[i] & rhs.rep[i]){
					return true;
					}
				}
			if(((~rep[i] & rhs.rep[i]) & partialBlockMask))
				return true;
			}
#ifdef _BETTER_BIPART
			else{
			unsigned sum = 0;
			for(i=0;i<nBlocks-1;i++){
				sum |= ((~rep[i] & rhs.rep[i]) & mask->rep[i]); 
				}
			if(sum)
				return true;
			if((~rep[i] & rhs.rep[i]) & (mask->rep[i] & partialBlockMask))
				return true;
			}
#else
		else{
		for(i=0;i<nBlocks-1;i++){
				if((~rep[i] & rhs.rep[i]) & mask->rep[i]){
					return true;
					}
				}
			if(((~rep[i] & rhs.rep[i]) & (mask->rep[i] & partialBlockMask)))
				return true;
				}
#endif
		return false;
			}

	bool ComplementHasIntersectionWithComplement(const Bipartition &rhs, const Bipartition *mask) const{
		int i;
		if(!mask){
		for(i=0;i<nBlocks-1;i++){
				if(~rep[i] & ~rhs.rep[i]){
					return true;
					}
				}
			if(((~rep[i] & ~rhs.rep[i]) & partialBlockMask))
				return true;
			}
#ifdef _BETTER_BIPART
			else{
			unsigned sum = 0;
			for(i=0;i<nBlocks-1;i++){
				sum |= ((~rep[i] & ~rhs.rep[i]) & mask->rep[i]); 
				}
			if(sum)
				return true;
			if((~rep[i] & ~rhs.rep[i]) & (mask->rep[i] & partialBlockMask))
				return true;
			}
#else
		else{
		for(i=0;i<nBlocks-1;i++){
				if((~rep[i] & ~rhs.rep[i]) & mask->rep[i]){
					return true;
				}
			}
			if(((~rep[i] & ~rhs.rep[i]) & (mask->rep[i] & partialBlockMask)))
				return true;
			}
#endif
		return false;
		}

	bool IsCompatibleWithNegativeBipartition(const Bipartition &bip) const{
		//To be consistent with a negative bipartition (constraint) neither the bipartition
		//or its complement can BE the constraint bipartition
		if(this->EqualsEquals(bip)==false && this->ComplementEqualsEquals(bip)==false) return true;
 		else return false;
		}

	bool IsCompatibleWithNegativeBipartitionWithMask(const Bipartition &bip, const Bipartition &mask) const{
		//To be consistent with a negative bipartition (constraint) neither the masked bipartition
		//or its masked complement can BE the constraint bipartition
		if(this->MaskedEqualsEquals(bip, mask)==false && this->MaskedComplementEqualsEquals(bip, mask)==false) return true;
		else return false;
		}

	bool IsASubsetOf(const Bipartition &target) const{
		//returns true if this is a subset of target
		int i;
		for(i=0;i<nBlocks-1;i++){
			if((target.rep[i] | rep[i]) != target.rep[i]) return false;
			}
		if(((target.rep[i] | rep[i]) & partialBlockMask) != (target.rep[i] & partialBlockMask)) return false;
		return true;
		}

	bool MaskedIsASubsetOf(const Bipartition &target, const Bipartition &mask) const{
		//returns true if this is a subset of target
		int i;
		for(i=0;i<nBlocks-1;i++){
			if(((target.rep[i] | rep[i]) & mask.rep[i]) != (target.rep[i] & mask.rep[i])) return false;
			}
		if(((target.rep[i] | rep[i]) & partialBlockMask & mask.rep[i]) != (target.rep[i] & partialBlockMask & mask.rep[i])) return false;
		return true;
		}

	bool ComplementIsASubsetOf(const Bipartition &target) const{
		//returns true if this is NOT a subset of target's complement
		int i;
		for(i=0;i<nBlocks-1;i++){
			if((target.rep[i] | ~rep[i]) != ~rep[i]) return true;
			}
		if(((target.rep[i] | ~rep[i]) & partialBlockMask) != (~rep[i] & partialBlockMask)) return false;
		return true;
		}

	Bipartition *TerminalBipart(int taxNum){
		assert(taxNum > 0 && taxNum <= ntax);
		ClearBipartition();
		//8-19-05 fixed this bit of stupidity
		//rep[(taxNum)/blockBits]|=(largestBlockDigit>>((taxNum-1)%blockBits));
		rep[(taxNum-1)/blockBits]|=(largestBlockDigit>>((taxNum-1)%blockBits));
		return this;
		}
	char * Output(){
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
	
/*
	void BinaryOutput(ofstream &out){
		int size = nBlocks * sizeof(unsigned int);
		out.write((char*) rep, size);
		}
*/
	void BinaryOutput(OUTPUT_CLASS &out){
		out.WRITE_TO_FILE(rep, sizeof(unsigned int), nBlocks);
		}

	void BinaryInput(FILE* &in){
		fread((char*) rep, sizeof(unsigned int), nBlocks, in);
		}

	vector<int> NodenumsFromBipart(){
		vector<int> nodes;
		for(int i=1;i<=ntax;i++) if(ContainsTaxon(i)) nodes.push_back(i);
		return nodes; 
		}
	
	void BipartFromNodenums(const vector<int> & nodes){
		ClearBipartition();
		Bipartition temp;
		for(vector<int>::const_iterator it = nodes.begin();it != nodes.end();it++)
			*this += temp.TerminalBipart(*it);
		}

	void BipartFromNodenums(const std::set<unsigned> & nodes){
		ClearBipartition();
		Bipartition temp;
		for(set<unsigned>::const_iterator it = nodes.begin();it != nodes.end();it++)
			*this += temp.TerminalBipart(*it);
		}
	};

bool BipartitionLessThan(const Bipartition &lhs, const Bipartition &rhs);

class Constraint{
	//a very simple class that just contains a bipartition and whether
	//it is a negative or positive constraint
	Bipartition con;
	bool positive;

	//now adding the potential for a backbone constraint
	//BACKBONE
	bool backbone;
	Bipartition backboneMask;

	//an outgroup is a special type of positive constraint
	bool outgroup;
public:
	//if these are true (and they almost always will be), it makes checking constraint validity much easier
	static bool allBackbone;
	static bool anyBackbone;
	static bool sharedMask;

	Constraint() {
		backbone=false;
		};
	Constraint(const Bipartition *b, bool pos){
		positive=pos;
		con=b;
		backbone=false;
		}

	Constraint(const Bipartition *b, const Bipartition *m, bool pos){
		positive=pos;
		con=b;
		backbone=true;
		backboneMask = m;
		}

	static void SetConstraintStatics(bool allBack, bool anyBack, bool oneMask);

	bool IsPositive() const {return positive;}

	bool IsBackbone() const {return backbone;}

	bool IsOutgroup() const {return outgroup;}

	const Bipartition *GetBipartition() const {return &con;}

	const Bipartition *GetBackboneMask() const {return &backboneMask;}

	void SetAsOutgroup() {outgroup = true;}

	void Standardize() {con.Standardize();}

	void ReadDotStarConstraint(const char *c){
		Bipartition b;
		size_t len=strlen(c);
		con.ClearBipartition();
		if(c[0] == '+') positive=true;
		else if(c[0] == '-') positive=false;
		else throw ErrorException("constraint string must start with \'+\' (positive constraint) or \'-\' (negative constraint)");
		for(unsigned i=1;i<len;i++){
			if(c[i] == '*') con += b.TerminalBipart(i);
			}
		con.Standardize();
		}
	bool ConstraintIsCompatibleWithConstraint(const Constraint &other) const{
		//this version takes a constraint and checks if the two are compatible with one another
		if(positive==true && other.IsPositive() == true){
			return con.IsCompatibleWithBipartition(*(other.GetBipartition()));		
			}
		return true;
		//need to figure out how to verify compatibility when at least one of the constraints is negative
		assert(0);
		return false;
		}

	//Checks whether a single bipartition is compatible with the constraint.  The mask passed in should already
	//contain the and'ing of a backbone mask and/or a mask representing taxa present in a growing tree
	bool BipartitionIsCompatibleWithConstraint(const Bipartition &other, const Bipartition *mask) const{
		if(positive){//positive constraint - determine bipartition compatibility with 4-point condition
			//DEBUG
			if(mask != NULL)
				return con.IsCompatibleWithBipartitionWithMask(other, *mask);
			else 
				return con.IsCompatibleWithBipartition(other);

/*			if(backbone){
				if(mask != NULL){
					Bipartition jointMask = *mask;
					jointMask.AndEquals(backboneMask);
					if(jointMask.CountOnBits() < 4) return true;
					return con.IsCompatibleWithBipartitionWithMask(other, jointMask);
					}
				else 
					return con.IsCompatibleWithBipartitionWithMask(other, backboneMask);
				}
			else{
				if(mask != NULL)
					return con.IsCompatibleWithBipartitionWithMask(other, *mask);
				else 
					return con.IsCompatibleWithBipartition(other);
				}
*/			}
		else{//negative constraints - for the bipartition to be incompatible with the negative
			//constraint it needs to acutally BE the constrained bipartition (potentially with a mask)
			if(mask != NULL) 
				return con.IsCompatibleWithNegativeBipartitionWithMask(other, *mask);
			else 
				return con.IsCompatibleWithNegativeBipartition(other);

/*			if(backbone){
				if(mask != NULL){
					Bipartition jointMask = *mask;
					jointMask.AndEquals(backboneMask);
					if(jointMask.CountOnBits() < 4) return true;
					return con.IsCompatibleWithNegativeBipartitionWithMask(other, jointMask);
					}
				else
					return con.IsCompatibleWithNegativeBipartitionWithMask(other, backboneMask);
				}
			else{
				if(mask != NULL)
					return con.IsCompatibleWithNegativeBipartitionWithMask(other, *mask);
				else
					return con.IsCompatibleWithNegativeBipartition(other);
				}
*/			}
		}

	void NumericalOutput(string &out){
		if(IsBackbone()) 
			con.NumericalOutput(out, &backboneMask);
		else 
			con.NumericalOutput(out, NULL);
		}
	};

#endif

