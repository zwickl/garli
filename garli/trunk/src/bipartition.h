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
#ifdef BETTER_BIPART
	static unsigned nBlocks;
	static unsigned blockBits;
	static unsigned ntax;
#else
	static int nBlocks;
	static int blockBits;
	static int ntax;
#endif
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

	Bipartition::~Bipartition(){
		if(rep!=NULL) delete []rep;
		rep=NULL;
		}	

	void BipartFromNodenums(const std::set<unsigned> & nodes);
	void BipartFromNodenums(const vector<int> & nodes);
	vector<int> NodenumsFromBipart();

	void BinaryOutput(OUTPUT_CLASS &out);
	void BinaryInput(FILE* &in);
	void NumericalOutput(string &out, const Bipartition *mask);
	char * Output();

	//manipulations
	void Complement();
	Bipartition *TerminalBipart(int taxNum);
	bool MakeJointMask(const Constraint &constr, const Bipartition *partialMask);
	void AndEquals(const Bipartition &rhs);
	void EqualsXORComplement(const Bipartition &rhs);
	void FillWithXORComplement(const Bipartition &lhs, const Bipartition &rhs);
	void Standardize();
	void FillAllBits();
	void FlipBits(const Bipartition &flip);
	void operator=(const Bipartition *rhs);
	void operator=(const Bipartition &rhs);
	static void SetBipartitionStatics(int);
	static void SetPartialBlockMask();
	void operator+=(const Bipartition *rhs);
	void operator-=(const Bipartition *rhs);
	void ClearBipartition();

	//tests
	bool MaskedIsASubsetOf(const Bipartition &target, const Bipartition &mask) const;
	bool ComplementIsASubsetOf(const Bipartition &target) const;
	bool IsASubsetOf(const Bipartition &target) const;
	bool IsIncompatibleWithBipartitionWithMask(const Bipartition &bip, const Bipartition &mask) const;
	bool IsCompatibleWithNegativeBipartitionWithMask(const Bipartition &bip, const Bipartition &mask) const;
	bool IsCompatibleWithNegativeBipartition(const Bipartition &bip) const;
	bool IsCompatibleWithBipartitionWithMask(const Bipartition &bip, const Bipartition &mask) const;
	bool IsCompatibleWithBipartition(const Bipartition &constr) const;
	bool ContainsTaxon(int t) const;
	bool ContainsFirstTaxon() const;
	bool MaskedComplementEqualsEquals(const Bipartition &rhs, const Bipartition &mask) const;
	bool ComplementEqualsEquals(const Bipartition &rhs) const;
	bool MaskedEqualsEquals(const Bipartition &rhs, const Bipartition &mask) const;
	bool EqualsEquals(const Bipartition &rhs) const;
	int FirstPresentTaxon() const;
	int FirstNonPresentTaxon() const;
	int CountOnBits() const;
	};

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
			if(mask != NULL) return con.IsCompatibleWithBipartitionWithMask(other, *mask);
			else return con.IsCompatibleWithBipartition(other);

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
			if(mask != NULL) return con.IsCompatibleWithNegativeBipartitionWithMask(other, *mask);
			else return con.IsCompatibleWithNegativeBipartition(other);

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
		if(IsBackbone()) con.NumericalOutput(out, &backboneMask);
		else con.NumericalOutput(out, NULL);
		}
	};

#endif

