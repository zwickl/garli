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
class Tree;
class Bipartition{

	static int nBlocks;
	static int blockBits;
	static int ntax;
	static unsigned int largestBlockDigit;
	static unsigned int allBitsOn;
	static unsigned int partialBlockMask;//this can be used to mask out the bits that
										// aren't used in the last block.  This becomes
										//important if we start doing complements. Bits
										//that represent actual taxa are ON

public:	
	static void SetBipartitionStatics(int);
	static void SetPartialBlockMask();

	std::vector<unsigned int> rep;

	Bipartition() {
		Clear();
	}

	Bipartition(const Bipartition &b) {
		rep = b.rep;
	}

	Bipartition(const char *c);

	void Clear() {
		rep.assign(nBlocks, 0);
	}

	
	void operator+=(const Bipartition &rhs) {
		for(int i = 0; i < nBlocks; i++) {
			rep[i] |= rhs.rep[i];
		}
	}

	void operator-=(const Bipartition &rhs) {
		//note that this assumes that rhs is a subset of this or vice versa!!!
		if(rhs.IsASubsetOf(*this)){
			for(int i = 0; i < nBlocks; i++) {
				rep[i] ^= rhs.rep[i];
			}
		}
		else {
			for(int i = 0; i < nBlocks; i++) {
				rep[i] |= ~(rhs.rep[i]);
			}
		}
	}

	void operator=(const Bipartition &rhs) {
		rep = rhs.rep;
	}

	int CountOnBits() const {
		int num=0;
		for(int i = 1;i <= ntax; i++) {
			if(ContainsTaxon(i))
				num++;
		}
		return num;
	}

	void FlipBits(const Bipartition & flip){
		//this just flips the bits in the passed in bipartition
		EqualsXORComplement(flip);
	}


	void FillAllBits(){
		for(int i=0;i<nBlocks;i++)
			rep[i] = allBitsOn;
	}

	bool EqualsEquals(const Bipartition &rhs) const{
		//assert(this->ContainsTaxon(1));
		//assert(rhs.ContainsTaxon(1));
		int i;
		for(i=0;i<nBlocks-1;i++){
			if(rep[i]!=rhs.rep[i])
				return false;
		}
		if ((rep[i]&partialBlockMask)!=((rhs.rep[i])&partialBlockMask))
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
		for(i=0;i<nBlocks-1;i++) {
			if((~rep[i] & mask.rep[i]) != (rhs.rep[i] & mask.rep[i]))
				return false;
		}
		if((~rep[i] & partialBlockMask & mask.rep[i]) != ((rhs.rep[i]) & partialBlockMask & mask.rep[i]))
			return false;
		return true;
	}

	void Standardize() {
		if (!ContainsTaxon(1))
			Complement();
	}
	
	int FirstPresentTaxon() const {
		int blk=0;
		unsigned int tmp = rep[blk];
		while (tmp==0)
			tmp = rep[++blk];
		
		int t = blockBits*blk + 1;
		while( ! (tmp & largestBlockDigit)) {
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

	bool ContainsTaxon(int t) const {
		unsigned int tmp = rep[(t-1)/blockBits];
		return (tmp & largestBlockDigit>>((t-1)%blockBits));
	}

	void FillWithXORComplement(const Bipartition &lhs, const Bipartition &rhs) {
		for(int i=0;i<nBlocks;i++) {
			rep[i] = lhs.rep[i] ^ (~rhs.rep[i]);
		}
	}

	void EqualsXORComplement(const Bipartition &rhs) {
		for(int i=0;i<nBlocks;i++){
			rep[i] = rep[i] ^ (~rhs.rep[i]);
		}
	}

	void AndEquals(const Bipartition &rhs) {
		//sets the calling bipart to be the bitwise AND of it and the argument
		for(int i=0;i<nBlocks;i++){
			rep[i] = rep[i] & rhs.rep[i];
		}
	}

	bool IsCompatibleWithBipartition(const Bipartition &constr) const {
		//using buneman's 4 point condition.  At least one of the four intersections must be empty
		bool compat=true;
		int i;
		//A & B
		for (i=0;i<nBlocks-1;i++) {
			if ((rep[i] & constr.rep[i]) != 0) {
				compat=false;
				break;
			}
		}
		if (compat && (((rep[i] & constr.rep[i]) & partialBlockMask) == 0))
			return true;
		//A & ~B
		compat=true;
		for (i=0;i<nBlocks-1;i++){
			if((rep[i] & ~(constr.rep[i])) != 0) {
				compat=false;
				break;
			}
		}
		if (compat==true && (((rep[i] & ~(constr.rep[i])) & partialBlockMask) == 0))
			return true;
		//~A & B
		compat=true;
		for(i=0;i<nBlocks-1;i++) {
			if((~rep[i] & constr.rep[i]) != 0) {
				compat=false;
				break;
			}
		}
		if(compat==true && (((~rep[i] & constr.rep[i]) & partialBlockMask) == 0))
			return true;
		//~A & ~B
		for(i=0;i<nBlocks-1;i++){
			if((~rep[i] & ~(constr.rep[i])) != 0)
				return false;
		}
		return (((~rep[i] & ~(constr.rep[i])) & partialBlockMask) == 0);
		}

	bool IsCompatibleWithBipartitionWithMask(const Bipartition &bip, const Bipartition &mask) const{
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

	bool IsCompatibleWithNegativeBipartition(const Bipartition &bip) const{
		//To be consistent with a negative bipartition (constraint) neither the bipartition
		//or its complement can BE the constraint bipartition
		return (this->EqualsEquals(bip)==false && this->ComplementEqualsEquals(bip)==false);
	}

	bool IsCompatibleWithNegativeBipartitionWithMask(const Bipartition &bip, const Bipartition &mask) const{
		//To be consistent with a negative bipartition (constraint) neither the masked bipartition
		//or its masked complement can BE the constraint bipartition
		return (this->MaskedEqualsEquals(bip, mask)==false && this->MaskedComplementEqualsEquals(bip, mask)==false);
	}

	bool IsIncompatibleWithBipartitionWithMask(const Bipartition &bip, const Bipartition &mask) const{
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

	Bipartition &TerminalBipart(int taxNum) {
		assert(taxNum > 0 && taxNum <= ntax);
		Clear();
		rep[(taxNum-1)/blockBits] |= (largestBlockDigit>>((taxNum-1)%blockBits));
		return *this;
	}

	std::string Output() const ;
	
	void NumericalOutput(string &out, const Bipartition *mask) const;

	void BinaryOutput(OUTPUT_CLASS &out) const {
		out.WRITE_TO_FILE(&rep[0], sizeof(unsigned int), nBlocks);
		}

	void BinaryInput(FILE* &in){
		rep.resize(nBlocks);
		fread((char*) &rep[0], sizeof(unsigned int), nBlocks, in);
		}

	vector<int> NodenumsFromBipart() const {
		vector<int> nodes;
		for(int i=1;i<=ntax;i++) {
			if(ContainsTaxon(i))
				nodes.push_back(i);
		}
		return nodes; 
	}
	
	void BipartFromNodenums(const vector<int> & nodes) {
		Clear();
		Bipartition temp;
		for(vector<int>::const_iterator it = nodes.begin();it != nodes.end();it++)
			*this += temp.TerminalBipart(*it);
		}


	void BipartFromNodenums(const std::set<unsigned> & nodes) {
		Clear();
		Bipartition temp;
		for(set<unsigned>::const_iterator it = nodes.begin();it != nodes.end();it++)
			*this += temp.TerminalBipart(*it);
		}

	bool MakeJointMask(const Constraint &constr, const Bipartition *partialMask);
	friend class Swap;
	friend bool BipartitionLessThan(const Bipartition &lhs, const Bipartition &rhs);
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
	Constraint() {
		backbone=false;
		};
	Constraint(const Bipartition &b, bool pos){
		positive=pos;
		con=b;
		backbone=false;
		}

	Constraint(const Bipartition &b, const Bipartition &m, bool pos){
		positive=pos;
		con=b;
		backbone=true;
		backboneMask = m;
		}

	bool IsPositive() const {return positive;}

	bool IsBackbone() const {return backbone;}

	bool IsOutgroup() const {return outgroup;}

	const Bipartition &GetBipartition() const {return con;}

	const Bipartition &GetBackboneMask() const {return backboneMask;}

	void SetAsOutgroup() {outgroup = true;}

	void Standardize() {con.Standardize();}

	void ReadDotStarConstraint(const char *c) {
		Bipartition b;
		size_t len=strlen(c);
		con.Clear();
		if(c[0] == '+')
			positive=true;
		else if(c[0] == '-')
			positive=false;
		else
			throw ErrorException("constraint string must start with \'+\' (positive constraint) or \'-\' (negative constraint)");
		for(unsigned i=1;i<len;i++){
			if(c[i] == '*')
				con += b.TerminalBipart(i);
		}
		con.Standardize();
	}

	bool ConstraintIsCompatibleWithConstraint(const Constraint &other) const{
		//this version takes a constraint and checks if the two are compatible with one another
		if (positive==true && other.IsPositive() == true) {
			return con.IsCompatibleWithBipartition(other.GetBipartition());		
		}
		//need to figure out how to verify compatibility when at least one of the constraints is negative
		assert(0);
		return true;
		return false;
	}

	//Checks whether a single bipartition is compatible with the constraint.  The mask passed in should already
	//contain the and'ing of a backbone mask and/or a mask representing taxa present in a growing tree
	bool BipartitionIsCompatibleWithConstraint(const Bipartition &other, const Bipartition *mask) const{
		if (positive) {
			if(mask != NULL)
				return con.IsCompatibleWithBipartitionWithMask(other, *mask);
			return con.IsCompatibleWithBipartition(other);
		}
		if(mask != NULL)
			return con.IsCompatibleWithNegativeBipartitionWithMask(other, *mask);
		return con.IsCompatibleWithNegativeBipartition(other);
	}

	void NumericalOutput(string &out) const{
		if(IsBackbone())
			con.NumericalOutput(out, &backboneMask);
		con.NumericalOutput(out, NULL);
	}

	bool CheckTree(const Tree &) const;

	};

#endif

