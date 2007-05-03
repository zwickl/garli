// GARLI version 0.952b2 source code
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


#ifndef BIPARTITION
#define BIPARTITION

#include <cmath>
#include <string>
#include <vector>
#include <cassert>

using namespace std;

#include "errorexception.h"

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

	~Bipartition();

	void ClearBipartition(){
		for(int i=0;i<nBlocks;i++)
			rep[i]=0;
		}

	static void SetBipartitionStatics(int);
	static void SetPartialBlockMask();
	
	void operator+=(const Bipartition *rhs){
		for(int i=0;i<nBlocks;i++){
			rep[i]|=rhs->rep[i];
			}
		}

	void operator-=(const Bipartition *rhs){
		//note that this assumes that rhs is a subset of this or vice versa!!!
		if(rhs->IsASubsetOf(this)){
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

	void FlipBits(Bipartition *flip){
		//this just flips the bits in the passed in bipartition
		EqualsXORComplement(flip);
		}

	bool EqualsEquals(const Bipartition *rhs){
		//assert(this->ContainsTaxon(1));
		//assert(rhs->ContainsTaxon(1));
		int i;
		for(i=0;i<nBlocks-1;i++){
			if(rep[i]!=rhs->rep[i]) return false;
			}
		if((rep[i]&partialBlockMask)!=((rhs->rep[i])&partialBlockMask)) return false;
		return true;
		}

	bool MaskedEqualsEquals(const Bipartition *rhs, const Bipartition *mask){
		//assert(this->ContainsTaxon(1));
		//assert(rhs->ContainsTaxon(1));
		int i;
		for(i=0;i<nBlocks-1;i++){
			if((rep[i] & mask->rep[i]) != (rhs->rep[i] & mask->rep[i])) return false;
			}
		if((rep[i] & partialBlockMask & mask->rep[i]) != ((rhs->rep[i]) & partialBlockMask & mask->rep[i])) return false;
		return true;
		}

	void Complement(){
		for(int i=0;i<nBlocks;i++){
			rep[i] = ~rep[i];
			}
		}

	bool ComplementEqualsEquals(const Bipartition *rhs){
		int i;
		for(i=0;i<nBlocks-1;i++){
			if(~rep[i]!=rhs->rep[i]) return false;
			}
			
		if(((~rep[i])&partialBlockMask)!=(rhs->rep[i]&partialBlockMask)) return false;
		
		return true;
		}
	
	void Standardize(){
		if(ContainsTaxon(1)==true) return;
		else{
			for(int i=0;i<nBlocks;i++){
				rep[i]=~rep[i];
				}
			}
		}
	
	int FirstPresentTaxon() const{
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

	int FirstNonPresentTaxon() const{
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

	bool ContainsTaxon(int t) const{
		unsigned int tmp=rep[(t-1)/blockBits];
		if(tmp & largestBlockDigit>>((t-1)%blockBits)) return true;
		return false;
		}

	void FillWithXORComplement(const Bipartition *lhs, const Bipartition *rhs){
		for(int i=0;i<nBlocks;i++){
			rep[i] = lhs->rep[i] ^ (~rhs->rep[i]);
			}
		}

	void EqualsXORComplement(const Bipartition *rhs){
		for(int i=0;i<nBlocks;i++){
			rep[i] = rep[i] ^ (~rhs->rep[i]);
			}
		}

	void AndEquals(const Bipartition *rhs){
		//sets the calling bipart to be the bitwise AND of it and the argument
		for(int i=0;i<nBlocks;i++){
			rep[i] = rep[i] & rhs->rep[i];
			}
		}

	bool IsCompatibleWithBipartition(const Bipartition *constr){
		//using buneman's 4 point condition.  At least one of the four intersections must be empty
		bool compat=true;
		int i;
		//A & B
		for(i=0;i<nBlocks-1;i++){
			if((rep[i] & constr->rep[i]) == 0) continue;
			else{
				compat=false;
				break;
				}
			}
		if(compat==true && (((rep[i] & constr->rep[i]) & partialBlockMask) == 0)) return true;
		//A & ~B
		compat=true;
		for(i=0;i<nBlocks-1;i++){
			if((rep[i] & ~(constr->rep[i])) == 0) continue;
			else{
				compat=false;
				break;
				}
			}
		if(compat==true && (((rep[i] & ~(constr->rep[i])) & partialBlockMask) == 0)) return true;
		//~A & B
		compat=true;
		for(i=0;i<nBlocks-1;i++){
			if((~rep[i] & constr->rep[i]) == 0) continue;
			else{
				compat=false;
				break;
				}
			}
		if(compat==true && (((~rep[i] & constr->rep[i]) & partialBlockMask) == 0)) return true;
		//~A & ~B
		compat=true;
		for(i=0;i<nBlocks-1;i++){
			if((~rep[i] & ~(constr->rep[i])) == 0) continue;
			else{
				compat=false;
				break;
				}
			}
		if(compat==true && (((~rep[i] & ~(constr->rep[i])) & partialBlockMask) == 0)) return true;
		return false;
		}

	bool IsCompatibleWithNegativeBipartition(const Bipartition *bip){
		//To be consistent with a negative bipartition (constraint) neither the bipartition
		//or its complement can equal the constraint bipartition
		if(this->EqualsEquals(bip)==false && this->ComplementEqualsEquals(bip)==false) return true;
		else return false;
		}

	bool IsCompatibleWithBipartitionWithMask(const Bipartition *bip, const Bipartition *mask) const{
		//using buneman's 4 point condition.  At least one of the four intersections must be empty
		//typically this will be called with 'this' as a POSITIVE constraint
		bool compat=true;
		int i;

		//A & B
		for(i=0;i<nBlocks-1;i++){
			if(((rep[i] & bip->rep[i]) & mask->rep[i]) == 0) continue;
			else{
				compat=false;
				break;
				}
			}
		if(compat==true && ((((rep[i] & bip->rep[i]) & mask->rep[i]) & partialBlockMask) == 0)) return true;
		//A & ~B
		compat=true;
		for(i=0;i<nBlocks-1;i++){
			if(((rep[i] & ~bip->rep[i]) & mask->rep[i]) == 0) continue;
			else{
				compat=false;
				break;
				}
			}
		if(compat==true && ((((rep[i] & ~bip->rep[i]) & mask->rep[i]) & partialBlockMask) == 0)) return true;
		//~A & B
		compat=true;
		for(i=0;i<nBlocks-1;i++){
			if(((~rep[i] & bip->rep[i]) & mask->rep[i]) == 0) continue;
			else{
				compat=false;
				break;
				}
			}
		if(compat==true && ((((~rep[i] & bip->rep[i]) & mask->rep[i]) & partialBlockMask) == 0)) return true;
		//~A & ~B
		compat=true;
		for(i=0;i<nBlocks-1;i++){
			if(((~rep[i] & ~bip->rep[i]) & mask->rep[i]) == 0) continue;
			else{
				compat=false;
				break;
				}
			}
		if(compat==true && ((((~rep[i] & ~bip->rep[i]) & mask->rep[i]) & partialBlockMask) == 0)) return true;
		return false;
/*
		//A & B
		for(i=0;i<nBlocks-1;i++){
			if(((rep[i] & mask->rep[i]) & (bip->rep[i] & mask->rep[i])) == 0) continue;
			else{
				compat=false;
				break;
				}
			}
		if(compat==true && ((((rep[i] & mask->rep[i]) & (bip->rep[i] & mask->rep[i])) & partialBlockMask) == 0)) return true;
		//A & ~B
		compat=true;
		for(i=0;i<nBlocks-1;i++){
			if(((rep[i] & mask->rep[i]) & (~bip->rep[i] & mask->rep[i])) == 0) continue;
			else{
				compat=false;
				break;
				}
			}
		if(compat==true && ((((rep[i] & mask->rep[i]) & (~bip->rep[i]& mask->rep[i])) & partialBlockMask) == 0)) return true;
		//~A & B
		compat=true;
		for(i=0;i<nBlocks-1;i++){
			if(((~rep[i] & mask->rep[i]) & (bip->rep[i] & mask->rep[i])) == 0) continue;
			else{
				compat=false;
				break;
				}
			}
		if(compat==true && ((((~rep[i] & mask->rep[i]) & (bip->rep[i] & mask->rep[i])) & partialBlockMask) == 0)) return true;
		//~A & ~B
		compat=true;
		for(i=0;i<nBlocks-1;i++){
			if(((~rep[i] & mask->rep[i]) & (~bip->rep[i] & mask->rep[i])) == 0) continue;
			else{
				compat=false;
				break;
				}
			}
		if(compat==true && ((((~rep[i] & mask->rep[i]) & (~bip->rep[i] & mask->rep[i])) & partialBlockMask) == 0)) return true;
		return false;
*/
		}

	bool IsIncompatibleWithBipartitionWithMask(const Bipartition *bip, const Bipartition *mask) const{
		//using buneman's 4 point condition.  At none of the four intersections must be empty
		//typically this will be called with 'this' as a NEGATIVE constraint
		bool compat=true;
		int i;

		//first check if there is no intersection between the mask and the constraint
		//or the complement of the constraint. If not, the question of incompatability
		//is moot (but it is consistent with incompatability, so return true)
		bool intersect=false;
		for(i=0;i<nBlocks-1;i++){
			if((rep[i] & mask->rep[i]) == 0) continue;
			else{
				intersect=true;
				break;
				}
			}
		if(intersect==false && ((((rep[i] & mask->rep[i]) & partialBlockMask) == 0))) return true;

		intersect=false;
		for(i=0;i<nBlocks-1;i++){
			if((~rep[i] & mask->rep[i]) == 0) continue;
			else{
				intersect=true;
				break;
				}
			}
		if(intersect==false && ((((~rep[i] & mask->rep[i]) & partialBlockMask) == 0))) return true;

		//A & B
		for(i=0;i<nBlocks-1;i++){
			if(((rep[i] & mask->rep[i]) & (bip->rep[i] & mask->rep[i])) == 0) continue;
			else{
				compat=false;
				break;
				}
			}
		if(compat==true && ((((rep[i] & mask->rep[i]) & (bip->rep[i] & mask->rep[i])) & partialBlockMask) == 0)) return true;
		//A & ~B
		compat=true;
		for(i=0;i<nBlocks-1;i++){
			if(((rep[i] & mask->rep[i]) & (~bip->rep[i] & mask->rep[i])) == 0) continue;
			else{
				compat=false;
				break;
				}
			}
		if(compat==true && ((((rep[i] & mask->rep[i]) & (~bip->rep[i]& mask->rep[i])) & partialBlockMask) == 0)) return true;
		//~A & B
		compat=true;
		for(i=0;i<nBlocks-1;i++){
			if(((~rep[i] & mask->rep[i]) & (bip->rep[i] & mask->rep[i])) == 0) continue;
			else{
				compat=false;
				break;
				}
			}
		if(compat==true && ((((~rep[i] & mask->rep[i]) & (bip->rep[i] & mask->rep[i])) & partialBlockMask) == 0)) return true;
		//~A & ~B
		compat=true;
		for(i=0;i<nBlocks-1;i++){
			if(((~rep[i] & mask->rep[i]) & (~bip->rep[i] & mask->rep[i])) == 0) continue;
			else{
				compat=false;
				break;
				}
			}
		if(compat==true && ((((~rep[i] & mask->rep[i]) & (~bip->rep[i] & mask->rep[i])) & partialBlockMask) == 0)) return true;
		return false;
		}

	bool IsASubsetOf(const Bipartition *target) const{
		//returns true if this is a subset of target
		int i;
		for(i=0;i<nBlocks-1;i++){
			if((target->rep[i] | rep[i]) != target->rep[i]) return false;
			}
		if(((target->rep[i] | rep[i]) & partialBlockMask) != (target->rep[i] & partialBlockMask)) return false;
		return true;
		}

	bool MaskedIsASubsetOf(const Bipartition *target, const Bipartition *mask) const{
		//returns true if this is a subset of target
		int i;
		for(i=0;i<nBlocks-1;i++){
			if(((target->rep[i] | rep[i]) & mask->rep[i]) != (target->rep[i] & mask->rep[i])) return false;
			}
		if(((target->rep[i] | rep[i]) & partialBlockMask & mask->rep[i]) != (target->rep[i] & partialBlockMask & mask->rep[i])) return false;
		return true;
		}

	bool ComplementIsASubsetOf(const Bipartition *target) const{
		//returns true if this is NOT a subset of target's complement
		int i;
		for(i=0;i<nBlocks-1;i++){
			if((target->rep[i] | ~rep[i]) != ~rep[i]) return true;
			}
		if(((target->rep[i] | ~rep[i]) & partialBlockMask) != (~rep[i] & partialBlockMask)) return false;
		return true;
		}

	Bipartition *TerminalBipart(int taxNum){
		assert(taxNum > 0 && taxNum <= ntax);
		for(int i=0;i<nBlocks;i++) rep[i]=0;
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
		
	vector<int> NodenumsFromBipart(){
		vector<int> nodes;
		for(int i=1;i<=ntax;i++) if(ContainsTaxon(i)) nodes.push_back(i);
		return nodes; 
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

public:
	Constraint() {};
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

	const bool IsPositive() {return positive;}

	const bool IsBackbone() {return backbone;}

	const Bipartition *GetBipartition(){return &con;}

	const Bipartition *GetBackboneMask(){return &backboneMask;}

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
	bool IsCompatibleWithConstraint(Constraint *other){
		//this version takes a constraint and checks if the two are compatible with one another
		if(positive==true && other->IsPositive() == true){
			return con.IsCompatibleWithBipartition(other->GetBipartition());		
			}
		return true;
		//need to figure out how to verify compatibility when at least one of the constraints is negative
		assert(0);
		return false;
		}
	bool IsCompatibleWithConstraint(const Bipartition *other){
		if(positive == true){
			//BACKBONE
			if(backbone == true) return con.IsCompatibleWithBipartitionWithMask(other, &backboneMask);		
			else return con.IsCompatibleWithBipartition(other);		
			}
		else return con.IsCompatibleWithNegativeBipartition(other);
		return false;
		}
	bool IsCompatibleWithConstraintWithMask(const Bipartition *other, const Bipartition *mask){
		if(positive == true)
			return con.IsCompatibleWithBipartitionWithMask(other, mask);		
		else return con.IsIncompatibleWithBipartitionWithMask(other, mask);
		}
	};

#endif

