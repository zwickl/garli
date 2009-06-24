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

int Bipartition::nBlocks;
int Bipartition:: blockBits;
int Bipartition::ntax;
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
	Bipartition::allBitsOn=(unsigned int)(pow(2.0, Bipartition::blockBits)-1);
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
		if(partialMask != NULL)//in this case we'll need to test for meaningful intersection below
			this->AndEquals(*partialMask);
		else//here we don't need to check, since a backbone constraint and its own mask must be meaningful
			return true;
		}
	else if(partialMask != NULL){
		//in this case we'll need to test for meaningful intersection below
		*this = *(partialMask);
		}
	else{
		FillAllBits();
		return true;
		}

	Bipartition temp;
	temp = constr.GetBipartition();
	temp.AndEquals(*this);
	if(temp.MoreThanOneBitSet() == false)
		return false;

	temp = constr.GetBipartition();
	temp.Complement();
	temp.AndEquals(*this);
	if(temp.MoreThanOneBitSet() == false)
		return false;

	return true;
	}
	
bool Bipartition::IsCompatibleWithBipartition(const Bipartition &constr) const{
	//using buneman's 4 point condition.  At least one of the four intersections must be empty
	//A & B
	if(HasIntersection(constr, NULL) == false)
		return true;
	//A & ~B
	if(HasIntersectionWithComplement(constr, NULL) == false)
		return true;
	//~A & B
	if(ComplementHasIntersection(constr, NULL) == false)
		return true;
	//~A & ~B
	if(ComplementHasIntersectionWithComplement(constr, NULL) == false)
		return true;
	return false;
	}


bool Bipartition::OldIsCompatibleWithBipartition(const Bipartition &constr) const{
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

bool Bipartition::IsCompatibleWithBipartitionWithMask(const Bipartition &constr, const Bipartition &mask) const{
	//using buneman's 4 point condition.  At least one of the four intersections must be empty
	//A & B
	if(HasIntersection(constr, &mask) == false)
		return true;
	//A & ~B
	if(HasIntersectionWithComplement(constr, &mask) == false)
		return true;
	//~A & B
	if(ComplementHasIntersection(constr, &mask) == false)
		return true;
	//~A & ~B
	if(ComplementHasIntersectionWithComplement(constr, &mask) == false)
		return true;
	return false;
	}

bool Bipartition::OldIsCompatibleWithBipartitionWithMask(const Bipartition &bip, const Bipartition &mask) const{
	//using buneman's 4 point condition.  At least one of the four intersections must be empty
	//typically this will be called with 'this' as a POSITIVE constraint
	bool compat=true;
	int i;

	//A & B
	for(i=0;i<nBlocks-1;i++){
		if(((rep[i] & bip.rep[i]) & mask.rep[i]) == 0) 
			continue;
		else{
			compat=false;
			break;
			}
		}
	if(compat==true && ((((rep[i] & bip.rep[i]) & mask.rep[i]) & partialBlockMask) == 0)) 
		return true;
	//A & ~B
	compat=true;
	for(i=0;i<nBlocks-1;i++){
		if(((rep[i] & ~bip.rep[i]) & mask.rep[i]) == 0) 
			continue;
		else{
			compat=false;
			break;
			}
		}
	if(compat==true && ((((rep[i] & ~bip.rep[i]) & mask.rep[i]) & partialBlockMask) == 0)) 
		return true;
	//~A & B
	compat=true;
	for(i=0;i<nBlocks-1;i++){
		if(((~rep[i] & bip.rep[i]) & mask.rep[i]) == 0) 
			continue;
		else{
			compat=false;
			break;
			}
		}
	if(compat==true && ((((~rep[i] & bip.rep[i]) & mask.rep[i]) & partialBlockMask) == 0)) 
		return true;
	//~A & ~B
	compat=true;
	for(i=0;i<nBlocks-1;i++){
		if(((~rep[i] & ~bip.rep[i]) & mask.rep[i]) == 0) 
			continue;
		else{
			compat=false;
			break;
			}
		}
	if(compat==true && ((((~rep[i] & ~bip.rep[i]) & mask.rep[i]) & partialBlockMask) == 0)) 
		return true;
	return false;
	}

//This isn't used currently and may not be working
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