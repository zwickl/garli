// GARLI version 1.00 source code
// Copyright 2005-2010 Derrick J. Zwickl
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
#include <iostream>
#include <sstream>

using namespace std;

#include "defs.h"
#include "utility.h"
#include "linalg.h"
#include "model.h"
#include "individual.h"
#include "sequencedata.h"
#include "rng.h"

#undef ALIGN_MODEL

Profiler ProfCalcPmat("CalcPmat      ");
Profiler ProfCalcEigen("CalcEigen     ");
					 
extern rng rnd;
FLOAT_TYPE Model::mutationShape;

FLOAT_TYPE Model::maxPropInvar;

FLOAT_TYPE PointNormal (FLOAT_TYPE prob);
FLOAT_TYPE IncompleteGamma (FLOAT_TYPE x, FLOAT_TYPE alpha, FLOAT_TYPE LnGamma_alpha);
FLOAT_TYPE PointChi2 (FLOAT_TYPE prob, FLOAT_TYPE v);

GeneticCode *Model::code = NULL;
int Model::qmatLookup[62*62];

Model::~Model(){
	if(stateFreqs.empty() == false){
		for(int i=0;i<(int)stateFreqs.size();i++)
			delete stateFreqs[i];
		}

	if(relNucRates.empty() == false){
		if(nst==6  || nst == -1){
			//3/25/08 this needed to change a bit for arbitrary matrices
			//since some of the elements might be aliased
			for(int i=0;i<(int)relNucRates.size();i++){
				if(relNucRates[i] != NULL){
					for(int j=i+1;j<(int)relNucRates.size();j++){
						if(relNucRates[j] == relNucRates[i]) relNucRates[j] = NULL;
						}
					delete relNucRates[i];
					relNucRates[i] = NULL;
					}
				}
			}
		else if(nst==2){
			delete relNucRates[0];
			delete relNucRates[1];
			}
		else if(nst==1) delete relNucRates[0];
		}

	if(modSpec.IsCodon()){
		for(int r=0;r<NRateCats();r++){
			delete omegas[r];
			delete omegaProbs[r];
			}
		}

	if(propInvar != NULL) delete propInvar;

	if(alpha != NULL) delete alpha;

	for(vector<BaseParameter*>::iterator delit=paramsToMutate.begin();delit!=paramsToMutate.end();delit++)
		delete *(delit);

	Delete2DArray(eigvals);
	delete []eigvalsimag;
	delete []iwork;
	delete []work;
	delete []col;
	delete []indx;
	if(c_ijk != NULL)
		Delete2DArray(c_ijk);
	delete []EigValexp;
	delete []EigValderiv;
	delete []EigValderiv2;
	delete []blen_multiplier;

#ifndef ALIGN_MODEL
	Delete3DArray(eigvecs);
	Delete2DArray(teigvecs);
	Delete3DArray(inveigvecs);
	//Delete3DArray(pmat);
	Delete3DArray(pmat1);
	Delete3DArray(pmat2);
	Delete3DArray(qmat);
	Delete3DArray(tempqmat);
	Delete3DArray(deriv1);
	Delete3DArray(deriv2);

	#ifdef SINGLE_PRECISION_FLOATS
	Delete3DArray(fpmat1);
	Delete3DArray(fpmat2);	
	Delete3DArray(fderiv1);
	Delete3DArray(fderiv2);
	#endif
#else
	Delete2DAlignedArray(eigvecs);
	Delete2DAlignedArray(teigvecs);
	Delete2DAlignedArray(inveigvecs);
	Delete3DAlignedArray(pmat);
	Delete2DAlignedArray(qmat);
	Delete2DAlignedArray(tempqmat);
	Delete3DAlignedArray(deriv1);
	Delete3DAlignedArray(deriv2);
#endif
	}

void Model::AllocateEigenVariables(){
#ifndef ALIGN_MODEL
	//a bunch of allocation here for all of the qmatrix->eigenvector->pmatrix related variables
	eigvalsimag=new MODEL_FLOAT[nstates];
	iwork=new int[nstates];
	work=new MODEL_FLOAT[nstates];
	col=new MODEL_FLOAT[nstates];
	indx=new int[nstates];
	EigValexp=new MODEL_FLOAT[nstates*NRateCats()];
	EigValderiv=new MODEL_FLOAT[nstates*NRateCats()];
	EigValderiv2=new MODEL_FLOAT[nstates*NRateCats()];

	//create the matrix for the eigenvectors
	eigvecs=New3DArray<MODEL_FLOAT>(NRateCats(), nstates, nstates);

	//create a temporary matrix to hold the eigenvectors that will be destroyed during the invertization
	teigvecs=New2DArray<MODEL_FLOAT>(nstates,nstates);

	//create the matrix for the inverse eigenvectors
	inveigvecs=New3DArray<MODEL_FLOAT>(NRateCats(), nstates, nstates);	

	//allocate the pmats
	pmat1=New3DArray<MODEL_FLOAT>(NRateCats(), nstates, nstates);
	pmat2=New3DArray<MODEL_FLOAT>(NRateCats(), nstates, nstates);
	
#ifdef SINGLE_PRECISION_FLOATS
	//allocate single precision versions of the matrices
	fpmat1=New3DArray<FLOAT_TYPE>(NRateCats(), nstates, nstates);
	fpmat2=New3DArray<FLOAT_TYPE>(NRateCats(), nstates, nstates);
	fderiv1=New3DArray<FLOAT_TYPE>(NRateCats(), nstates, nstates);
	fderiv2=New3DArray<FLOAT_TYPE>(NRateCats(), nstates, nstates);

#endif
	
	//it is actually less efficient to precalc the c_ijk for codon models due to the immense
	//size of the matrix.  So don't allocate it at all.
	if(modSpec.IsCodon() == false){
		c_ijk=New2DArray<MODEL_FLOAT>(1,nstates*nstates*nstates);
		}
	else c_ijk = NULL;

	//allocate qmat and tempqmat
	//if this is a model with multiple qmats (like multi-omega models or mixtures)
	//it needs to be bigger
	if(modSpec.IsNonsynonymousRateHet() == false){
		qmat=New3DArray<MODEL_FLOAT>(1, nstates,nstates);
		tempqmat=New3DArray<MODEL_FLOAT>(1, nstates,nstates);
		blen_multiplier = new FLOAT_TYPE[1];
		eigvals=New2DArray<MODEL_FLOAT>(1, nstates);//eigenvalues
		}
	else{
		qmat=New3DArray<MODEL_FLOAT>(NRateCats(), nstates, nstates);
		tempqmat=New3DArray<MODEL_FLOAT>(NRateCats(), nstates, nstates);
		blen_multiplier = new FLOAT_TYPE[NRateCats()];
		eigvals=New2DArray<MODEL_FLOAT>(NRateCats(), nstates);//eigenvalues
		}

	deriv1=New3DArray<MODEL_FLOAT>(NRateCats(), nstates, nstates);
	deriv2=New3DArray<MODEL_FLOAT>(NRateCats(), nstates, nstates);
#else

	//a bunch of allocation here for all of the qmatrix->eigenvector->pmatrix related variables
	eigvals=new MODEL_FLOAT[nstates];//eigenvalues
	eigvalsimag=new MODEL_FLOAT[nstates];
	iwork=new int[nstates];
	work=new MODEL_FLOAT[nstates];
	col=new MODEL_FLOAT[nstates];
	indx=new int[nstates];
	c_ijk=new MODEL_FLOAT[nstates*nstates*nstates];	
	EigValexp=new MODEL_FLOAT[nstates*NRateCats()];	

	//create the matrix for the eigenvectors
	eigvecs=New2DAlignedArray<MODEL_FLOAT>(nstates,nstates);

	//create a temporary matrix to hold the eigenvectors that will be destroyed during the invertization
	teigvecs=New2DAlignedArray<MODEL_FLOAT>(nstates,nstates);

	//create the matrix for the inverse eigenvectors
	inveigvecs=New2DAlignedArray<MODEL_FLOAT>(nstates,nstates);	

	//allocate the pmat
	pmat=New3DAlignedArray<MODEL_FLOAT>(NRateCats(), nstates, nstates);

	//allocate qmat and tempqmat
	qmat=New2DAlignedArray<MODEL_FLOAT>(nstates,nstates);
	tempqmat=New2DAlignedArray<MODEL_FLOAT>(nstates,nstates);

	deriv1=New3DAlignedArray<MODEL_FLOAT>(NRateCats(), nstates, nstates);
	deriv2=New3DAlignedArray<MODEL_FLOAT>(NRateCats(), nstates, nstates);
#endif
	}

void Model::UpdateQMat(){
	//recalculate the qmat from the statefreqs and rates

	if(modSpec.IsCodon()){
		UpdateQMatCodon();
		return;
		}
	else if(modSpec.IsAminoAcid()){
		UpdateQMatAminoAcid();
		return;
		}
	
	if(nstates==4){
		qmat[0][0][1]=*relNucRates[0] * *stateFreqs[1];  //a * piC
		qmat[0][0][2]=*relNucRates[1] * *stateFreqs[2];  //b * piG
		qmat[0][0][3]=*relNucRates[2] * *stateFreqs[3];  //c * piT
		qmat[0][1][2]=*relNucRates[3] * *stateFreqs[2];  //d * piG
		qmat[0][1][3]=*relNucRates[4] * *stateFreqs[3];  //e * piT
		qmat[0][2][3]=*stateFreqs[3];  			//f(=1) * piT 
		qmat[0][1][0]=*relNucRates[0] * *stateFreqs[0];  //a * piA
		qmat[0][2][0]=*relNucRates[1] * *stateFreqs[0];  //b * piA
		qmat[0][2][1]=*relNucRates[3] * *stateFreqs[1];  //d * piC
		qmat[0][3][0]=*relNucRates[2] * *stateFreqs[0];  //c * piA
		qmat[0][3][1]=*relNucRates[4] * *stateFreqs[1];  //e * piC
		qmat[0][3][2]=*stateFreqs[2]; 			//f(=1) * piG
		}
	else {//this isn't being used - see UpdateQmatCodon and UpdateQmatAminoAcid
		//general nstate x nstate method	
		int rnum=0;
		for(int i=0;i<nstates;i++){
			for(int j=i+1;j<nstates;j++){
				qmat[0][i][j]=*relNucRates[rnum] * *stateFreqs[j];
				qmat[0][j][i]=*relNucRates[rnum] * *stateFreqs[i];
				rnum++;
				}
			}
		}
		
	//set diags to sum rows to 0
	MODEL_FLOAT sum;
	for(int x=0;x<nstates;x++){
		sum=ZERO_POINT_ZERO;
		for(int y=0;y<nstates;y++){
			if(x!=y) sum+=qmat[0][x][y];
			}
		qmat[0][x][x]=-sum;
		}

	//calculate the branch length rescaling factor
	blen_multiplier[0]=(ZERO_POINT_FIVE/((qmat[0][0][1]**stateFreqs[0])+(qmat[0][0][2]**stateFreqs[0])+(qmat[0][0][3]**stateFreqs[0])+(qmat[0][1][2]**stateFreqs[1])+(qmat[0][1][3]**stateFreqs[1])+(qmat[0][2][3]**stateFreqs[2])));
	}

void Model::FillQMatLookup(){
	//the code here is:
	//bit 1 = viable 1 nuc change path
	//bit 2 = transition
	//bit 4 = Nonsynonymous
	//bit 8 = A->C or C->A nuc change
	//bit 16 = A->G or G->A nuc change
	//bit 32 = A->T or T->A nuc change
	//bit 64 = C->G or G->C nuc change
	//bit 128 = C->T or T->C nuc change
	//bit 256 = G->T or T->G nuc change
	
	//although it seems a little wacky, I'm going to fill the 64x64 matrix, and then eliminate the
	//rows and columns that are stop codons, since they differ for different codes and the following 
	//stuff would be hell without regularity.  The static qmatLookup is only calculated once anyway.
	
	int tempqmatLookup[64*64];
	for(int q=0;q<64*64;q++) tempqmatLookup[q]=0;
	//its easier to do this in 4 x 4 blocks
	for(int i=0;i<16;i++){
		for(int j=0;j<16;j++){
			for(int ii=0;ii<4;ii++){
				for(int jj=0;jj<4;jj++){
					if(i==j){//on diagonal 4x4
						if(ii!=jj){
							//all the cells in this subsection are 1 nuc change away
							tempqmatLookup[64*(i*4+ii) + (j*4+jj)] = 1;
							if((ii+jj)%2 == 0)
								tempqmatLookup[64*(i*4+ii) + (j*4+jj)] |= 2;
							if((ii==0 && jj==1) || (ii==1 && jj==0)) tempqmatLookup[64*(i*4+ii) + (j*4+jj)] |= 8;
							else if((ii==0 && jj==2) || (ii==2 && jj==0)) tempqmatLookup[64*(i*4+ii) + (j*4+jj)] |= 16;
							else if((ii==0 && jj==3) || (ii==3 && jj==0)) tempqmatLookup[64*(i*4+ii) + (j*4+jj)] |= 32;
							else if((ii==1 && jj==2) || (ii==2 && jj==1)) tempqmatLookup[64*(i*4+ii) + (j*4+jj)] |= 64;
							else if((ii==1 && jj==3) || (ii==3 && jj==1)) tempqmatLookup[64*(i*4+ii) + (j*4+jj)] |= 128;
							else if((ii==2 && jj==3) || (ii==3 && jj==2)) tempqmatLookup[64*(i*4+ii) + (j*4+jj)] |= 256;
							}
						}
					else if(floor(i/4.0)==floor(j/4.0)){//near diagonal 4x4, some cells differ at the 2nd pos
						if(ii==jj){
							//the diagonal cells in this subsection are 1 nuc change away
							tempqmatLookup[64*(i*4+ii) + (j*4+jj)] = 1;
							if(abs(i-j) == 2)
								tempqmatLookup[64*(i*4+ii) + (j*4+jj)] |= 2;
							if((i%4==0 && j%4==1) || (i%4==1 && j%4==0)) tempqmatLookup[64*(i*4+ii) + (j*4+jj)] |= 8;
							else if((i%4==0 && j%4==2) || (i%4==2 && j%4==0)) tempqmatLookup[64*(i*4+ii) + (j*4+jj)] |= 16;
							else if((i%4==0 && j%4==3) || (i%4==3 && j%4==0)) tempqmatLookup[64*(i*4+ii) + (j*4+jj)] |= 32;
							else if((i%4==1 && j%4==2) || (i%4==2 && j%4==1)) tempqmatLookup[64*(i*4+ii) + (j*4+jj)] |= 64;
							else if((i%4==1 && j%4==3) || (i%4==3 && j%4==1)) tempqmatLookup[64*(i*4+ii) + (j*4+jj)] |= 128;
							else if((i%4==2 && j%4==3) || (i%4==3 && j%4==2)) tempqmatLookup[64*(i*4+ii) + (j*4+jj)] |= 256;							
							}
						}
					else{//far from diagonal 4x4, some cells differ at the 1nd pos
						if(i%4 == j%4){
							if(ii==jj){
								//the diagonal cells in this subsection are 1 nuc change away
								tempqmatLookup[64*(i*4+ii) + (j*4+jj)] = 1;
								if(abs(i-j) ==8)
									tempqmatLookup[64*(i*4+ii) + (j*4+jj)] |= 2;
								if((i/4==0 && j/4==1) || (i/4==1 && j/4==0)) tempqmatLookup[64*(i*4+ii) + (j*4+jj)] |= 8;
								else if((i/4==0 && j/4==2) || (i/4==2 && j/4==0)) tempqmatLookup[64*(i*4+ii) + (j*4+jj)] |= 16;
								else if((i/4==0 && j/4==3) || (i/4==3 && j/4==0)) tempqmatLookup[64*(i*4+ii) + (j*4+jj)] |= 32;
								else if((i/4==1 && j/4==2) || (i/4==2 && j/4==1)) tempqmatLookup[64*(i*4+ii) + (j*4+jj)] |= 64;
								else if((i/4==1 && j/4==3) || (i/4==3 && j/4==1)) tempqmatLookup[64*(i*4+ii) + (j*4+jj)] |= 128;
								else if((i/4==2 && j/4==3) || (i/4==3 && j/4==2)) tempqmatLookup[64*(i*4+ii) + (j*4+jj)] |= 256;	
								}
							}
						}
					}
				}
			}
		}

	//mark the nonsynonymous changes with |4 and the stops with -1
	for(int from=0;from<64;from++){
		for(int to=0;to<64;to++){
			int fromAA = code->CodonLookup(from);
			int toAA = code->CodonLookup(to);
			//if one of the codons is a stop
			if(fromAA == 20 || toAA == 20)
				tempqmatLookup[64*from + to] = -1;
			//if this is a viable 1 nucleotide change
			else if(tempqmatLookup[64*from + to] & 1){
				//if this is a nonsynonymous change
				if(fromAA != toAA) 
					tempqmatLookup[64*from + to] |= 4;
				}
			}
		}
#ifdef CODON_QMAT_HACK
//WHEN PLAYING HERE, REMEMBER THAT THESE INDECES ARE WITH RESPECT TO THE
//WHOLE 64X64 MATRIX
	//Put in whatever ad hoc alterations to the codon matrix
	//Hack in single changes for serine -> serine double hits
	//giving them a rate of omega X the first pos change
	//AGC->TCC
/*	tempqmatLookup[629] = 1 | 4 | 32;
	//AGT->TCT
	tempqmatLookup[759] = 1 | 4 | 32;
	//TCT->AGT
	tempqmatLookup[3531] = 1 | 4 | 32;
	//TCC->AGC
	tempqmatLookup[3401] = 1 | 4 | 32;
*/
//all transitions between
//the two sets of serine codons
	//AGC->TCN
	tempqmatLookup[628] = 1 | 4 | 32;
	tempqmatLookup[629] = 1 | 4 | 32;
	tempqmatLookup[630] = 1 | 4 | 32;
	tempqmatLookup[631] = 1 | 4 | 32;
	//AGT->TCN
	tempqmatLookup[756] = 1 | 4 | 32;
	tempqmatLookup[757] = 1 | 4 | 32;
	tempqmatLookup[758] = 1 | 4 | 32;
	tempqmatLookup[759] = 1 | 4 | 32;
	//TCN->AGC
	tempqmatLookup[3337] = 1 | 4 | 32;
	tempqmatLookup[3401] = 1 | 4 | 32;
	tempqmatLookup[3465] = 1 | 4 | 32;
	tempqmatLookup[3529] = 1 | 4 | 32;
	//TCN->AGT
	tempqmatLookup[3339] = 1 | 4 | 32;
	tempqmatLookup[3403] = 1 | 4 | 32;
	tempqmatLookup[3467] = 1 | 4 | 32;
	tempqmatLookup[3531] = 1 | 4 | 32;
	
#endif

	//remove the columns and rows representing stops
	int reducedCell = 0;
	for(int fullCell=0;fullCell<64*64;fullCell++){
		if(tempqmatLookup[fullCell] != -1)
			qmatLookup[reducedCell++] = tempqmatLookup[fullCell];
		}
	//assert(reducedCell == nstates*nstates);

/*	ofstream deb("debugQmat.log");
	for(int from=0;from<nstates;from++){
		for(int to=0;to<nstates;to++){
			deb << qmatLookup[from*nstates + to] << "\t";
			}
		deb << endl;
		}
*/	}


void Model::UpdateQMatCodon(){
//	ofstream dout("AAdist.log");

//	double composition[61]={0.12, 0.484, 0.12, 0.484, 0.0727, 0.0727, 0.0727, 0.0727, 0.236, 0.516, 0.236, 0.516, 0, 0, 0, 0, 0.324, 0.211, 0.324, 0.211, 0.142, 0.142, 0.142, 0.142, 0.236, 0.236, 0.236, 0.236, 0, 0, 0, 0, 0.335, 0.502, 0.335, 0.502, 0, 0, 0, 0, 0.269, 0.269, 0.269, 0.269, 0, 0, 0, 0, 0.0727, 0.0727, 0.516, 0.516, 0.516, 0.516, 1, 0.0473, 1, 0, 0, 0, 0};
//	double polarity[61]={0.79, 0.827, 0.79, 0.827, 0.457, 0.457, 0.457, 0.457, 0.691, 0.531, 0.691, 0.531, 0.037, 0.037, 0.0988, 0.037, 0.691, 0.679, 0.691, 0.679, 0.383, 0.383, 0.383, 0.383, 0.691, 0.691, 0.691, 0.691, 0, 0, 0, 0, 0.914, 1, 0.914, 1, 0.395, 0.395, 0.395, 0.395, 0.506, 0.506, 0.506, 0.506, 0.123, 0.123, 0.123, 0.123, 0.16, 0.16, 0.531, 0.531, 0.531, 0.531, 0.0741, 0.0617, 0.0741, 0, 0, 0, 0};
//	double molvol[61]={0.695, 0.317, 0.695, 0.317, 0.347, 0.347, 0.347, 0.347, 0.725, 0.647, 0.725, 0.647, 0.647, 0.647, 0.611, 0.647, 0.491, 0.557, 0.491, 0.557, 0.177, 0.177, 0.177, 0.177, 0.725, 0.725, 0.725, 0.725, 0.647, 0.647, 0.647, 0.647, 0.479, 0.305, 0.479, 0.305, 0.168, 0.168, 0.168, 0.168, 0, 0, 0, 0, 0.485, 0.485, 0.485, 0.485, 0.796, 0.796, 0.647, 0.647, 0.647, 0.647, 0.311, 1, 0.311, 0.647, 0.772, 0.647, 0.772};
	
	for(int w=0;w<NRateCats();w++){
		for(int i=0;i<nstates;i++){
			for(int j=0;j<nstates;j++){
				if(qmatLookup[i*nstates+j]==0){//if two codons are >1 appart
					qmat[w][i][j]=0.0;
					}
				else{
					qmat[w][i][j] = *stateFreqs[j];
					if(nst == 2){
						if(qmatLookup[i*nstates+j] & 2){//if the difference is a transition
							qmat[w][i][j] *= *relNucRates[1];
							}
						}
					else if(nst == 6){
						if(qmatLookup[i*nstates+j] & 8) qmat[w][i][j] *= *relNucRates[0];
						else if(qmatLookup[i*nstates+j] & 16) qmat[w][i][j] *= *relNucRates[1];
						else if(qmatLookup[i*nstates+j] & 32) qmat[w][i][j] *= *relNucRates[2];
						else if(qmatLookup[i*nstates+j] & 64) qmat[w][i][j] *= *relNucRates[3];
						else if(qmatLookup[i*nstates+j] & 128) qmat[w][i][j] *= *relNucRates[4];
						else if(qmatLookup[i*nstates+j] & 256) qmat[w][i][j] *= *relNucRates[5];
						}
					if(qmatLookup[i*nstates+j]&4){
						//this is where omega or AA property stuff will go
						//double phi=pow(abs((composition[i]-scalerC->val)/(composition[j]-scalerC->val)),Vc->val);
						//double phi=pow(abs((polarity[i]-scalerC->val)/(polarity[j]-scalerC->val)),Vc->val);
						//double phi=pow(exp(abs(polarity[i]-polarity[j])),Vc->val);
						//double phi=pow(exp(abs(composition[i]-composition[j])),Vc->val);
						//double phi=pow(1-abs(polarity[i]-polarity[j]),Vc->val);

	/*
						double compdist=abs(composition[i]-composition[j]);
						double poldist=abs(polarity[i]-polarity[j]);
						double voldist=abs(molvol[i]-molvol[j]);
						
						double comp, pol, vol;
	*/
						/*
						//this is essentially Yang's geometric relationship for each distance separately
						//there is a separate 'b' variable for each property, and a single 'a'
						double comp=exp(-Vc->val*abs(composition[i]-composition[j]));
						double pol=exp(-Vc->val*abs(composition[i]-composition[j]));
						double vol=exp(-Vc->val*abs(composition[i]-composition[j]));
						qmat[i][j] *= comp;
						qmat[i][j] *= scalerC->val; //overall omega type thing 
						qmat[i][j] *= pol;
						qmat[i][j] *= vol;	
	*/
						int nParams = 0;

						if(nParams==5){
							//this is essentially Yang's linear relationship for each distance separately
							//there is a separate 'b' variable for each property, and a single 'a'
	/*						comp=(1-scalerC->val*compdist);
							pol=(1-scalerP->val*poldist);
							vol=(1-scalerM->val*voldist);
							qmat[i][j] *= comp;

							qmat[i][j] *= omega->val; //overall omega type thing 
							qmat[i][j] *= pol;
							qmat[i][j] *= vol;
	*/						}

						if(nParams==0){
							qmat[w][i][j] *= *omegas[w]; //overall omega
							}

						if(nParams==1){
							//raw distances with overall omega
	/*						comp=1.0-compdist;
							pol=1.0-poldist;
							vol=1.0-voldist;
							qmat[i][j] *= comp;
							qmat[i][j] *= omega->val; //overall omega
							qmat[i][j] *= pol;
							qmat[i][j] *= vol;
	*/						}

						else if(nParams==4){
							//powered distances with overall omega
	/*						comp=pow((1.001-compdist),powerC->val-1.0);
							pol=pow((1.001-poldist),powerP->val-1.0);
							vol=pow((1.001-voldist),powerM->val-1.0);
							qmat[i][j] *= comp;
							qmat[i][j] *= omega->val; //overall omega type thing 
							qmat[i][j] *= pol;
							qmat[i][j] *= vol;
	*/						}

						else if(nParams==7){
							//powered distances with scalers and overall omega
	/*						comp=pow(1-scalerC->val * abs(composition[i]-composition[j]),powerC->val);
							pol=pow(1-scalerP->val * abs(polarity[i]-polarity[j]),powerP->val);
							vol=pow(1-scalerM->val * abs(molvol[i]-molvol[j]),powerM->val);
							qmat[i][j] *= comp;
							qmat[i][j] *= omega->val; //overall omega type thing 
							qmat[i][j] *= pol;
							qmat[i][j] *= vol;
	*/						}
					
	/*					//raw distances with overall omega
						double comp=1-abs(composition[i]-composition[j]);
						double pol=1-abs(polarity[i]-polarity[j]);
						double vol=1-abs(molvol[i]-molvol[j]);
						qmat[i][j] *= comp;
						qmat[i][j] *= scalerC->val; //overall omega type thing 
						qmat[i][j] *= pol;
						qmat[i][j] *= vol;	
	*/
	/*					//raw distances converted to single euclidian and not rescaled to new max, with overall omega
						double eucdist=1 - sqrt(comp*comp + pol*pol + vol*vol);
						if(eucdist<0){
							eucdist=0;
							}
						qmat[i][j] *= eucdist;
						qmat[i][j] *= scalerC->val; //overall omega type thing 					
	*/
	/*
						//powered distances converted to a euclidian with overall omega, not rescaled
						double comp=pow(1-abs(composition[i]-composition[j]),scalerC->val);
						double pol=pow(1-abs(polarity[i]-polarity[j]),powerP->val);
						double vol=pow(1-abs(molvol[i]-molvol[j]),powerM->val);
						double eucdist=1 - sqrt((1-comp)*(1-comp) + (1-pol)*(1-pol) + (1-vol)*(1-vol));
						qmat[i][j] *= eucdist;
						qmat[i][j] *= Vc->val; //overall omega type thing 
	*/					
	/*					//powered distances converted to a euclidian with overall omega, rescaled to max
						double comp=pow(1-abs(composition[i]-composition[j]),scalerC->val);
						double pol=pow(1-abs(polarity[i]-polarity[j]),powerP->val);
						double vol=pow(1-abs(molvol[i]-molvol[j]),powerM->val);
						double eucdist=1 - (sqrt((1-comp)*(1-comp) + (1-pol)*(1-pol) + (1-vol)*(1-vol)))/sqrt(3.0);
						qmat[i][j] *= eucdist;
						qmat[i][j] *= Vc->val; //overall omega type thing 
	*/
	/*
						//raw distances converted to single euclidian rescaled to euc max, with overall omega
						double eucdist=1.0 - sqrt(comp*comp + pol*pol + vol*vol) / sqrt(3)); //make sure to rescale the max here
						qmat[i][j] *= eucdist;
						qmat[i][j] *= scalerC->val; //overall omega type thing
	*/					
						//raw distances with estimated weights, converted to single euclidian, with overall omega
						//the weights are really only relative, with the composition weight fixed to 1
	/*					//distance rescaled to new max value
						double eucdist=sqrt(comp*comp + powerP->val*pol*pol + powerM->val*vol*vol);
						double eucmax=sqrt(1 + powerP->val*powerP->val + powerM->val*powerM->val);
						qmat[i][j] *= 1.0 - eucdist/eucmax;
						qmat[i][j] *= scalerC->val; //overall omega type thing
	*/
	/*					//raw distances converted to single euclidian rescaled to euc max, with overall omega
						//also with a parameter like the 'b' in Yang's linear model
						double eucdist=1.0 - Vc->val*sqrt(comp*comp + pol*pol + vol*vol) / sqrt(3); //make sure to rescale the max here
						if(eucdist<0) eucdist=0.0;
						qmat[i][j] *= eucdist;
						qmat[i][j] *= scalerC->val; //overall omega type thing
	*/
						//raw distances converted to single euclidian rescaled to euc max, with overall omega
	/*					//raised to an estimated power
						double eucdist=pow(1.0 - sqrt(comp*comp + pol*pol + vol*vol) / sqrt(3), powerP->val); //make sure to rescale the max here
						if(eucdist<0) eucdist=0.0;
						qmat[i][j] *= eucdist;
						qmat[i][j] *= scalerC->val; //overall omega type thing
	*/
	/*
						//alternative forumulation allowing the curve to go above 1.
						//infers "intercept" and power
						//double eucdist=sqrt(comp*comp + pol*pol + vol*vol) / sqrt(3.0); //make sure to rescale the max here
						double eucdist=comp; //make sure to rescale the max here
						eucdist=pow((1.0-(eucdist)*(eucdist-powerP->val)), scalerC->val);
						qmat[i][j] *= eucdist;
						qmat[i][j] *= Vc->val; //overall omega type thing
	*/
						//raw distances converted to single euclidian rescaled to euc max, with overall omega
	/*					//raised to an estimated power, with a scaler (ie 'b' on the euc) <0 set to 0
						double eucdist=1.0 - powerP->val * sqrt(comp*comp + pol*pol + vol*vol) / sqrt(3); //make sure to rescale the max here
						if(eucdist<0){
							eucdist=0.0;
							}
						eucdist=pow(eucdist, Vc->val);
						qmat[i][j] *= eucdist;
						qmat[i][j] *= scalerC->val; //overall omega type thing
	*/
						//raw distances with estimated weights, converted to single euclidian, with overall omega
						//the weights are really only relative, with the composition weight fixed to 1
						//distance rescaled to new max value, also with the b of Yang
	/*					double eucdist=scalerC->val * sqrt(comp*comp + powerP->val*pol*pol + powerM->val*vol*vol);
						double eucmax=sqrt(1 + powerP->val*powerP->val + powerM->val*powerM->val);
						eucdist = 1.0 - eucdist/eucmax;
						if(eucdist<0) eucdist=0.0;
						qmat[i][j] *= eucdist;
						qmat[i][j] *= Vc->val; //overall omega type thing
	*/
						//raw distances with estimated weights, converted to single euclidian, with overall omega
						//the weights are really only relative, with the composition weight fixed to 1
	/*					//distance rescaled to new max value, also with the b of Yang
						double eucdist=scalerC->val * sqrt(comp*comp + powerP->val*pol*pol + powerM->val*vol*vol);
						double eucmax=sqrt(1 + powerP->val*powerP->val + powerM->val*powerM->val);
						eucdist=1.0 - eucdist/eucmax;
						if(eucdist<0) eucdist=0.0;
						eucdist = pow(eucdist, scalerP->val);
						qmat[i][j] *= eucdist;
						qmat[i][j] *= Vc->val; //overall omega type thing
	*/
						/*
						//this is called "3propPowerDenom"
						double comp=pow(1-(abs(composition[i]-composition[j])/(1+scalerC->val)),Vc->val);
						double pol=pow(1-(abs(polarity[i]-polarity[j])/(1+scalerP->val)),powerP->val);
						double vol=pow(1-(abs(molvol[i]-molvol[j])/(1+scalerM->val)),powerM->val);
						qmat[i][j] *= comp;
						qmat[i][j] *= pol;
						qmat[i][j] *= vol;
						*/
						}
					}
				}
			}
		}

	//set diags to sum rows to 0 and calculate the branch length rescaling factor
	//note the there is really only a single rescaler, but it is stored separately 
	//for each omega model
	double sum, weightedDiagSum;
	blen_multiplier[0] = 0.0;

	for(int w=0;w<NRateCats();w++){
		weightedDiagSum = 0.0;
		for(int x=0;x<nstates;x++){
			sum = 0.0;
			for(int y=0;y<nstates;y++){
				if(x!=y) sum+=qmat[w][x][y];
				}
			qmat[w][x][x]=-sum;
			weightedDiagSum += sum * *stateFreqs[x];
			}
		blen_multiplier[0] += weightedDiagSum * *omegaProbs[w];
		}
	//note that although there si one blen multiplier per matrix, they are all set to the same value
	blen_multiplier[0] = ONE_POINT_ZERO / blen_multiplier[0];
	for(int i=1;i<NRateCats();i++)
		blen_multiplier[i] = blen_multiplier[0];
	}

//This just duplicates what happens at the end of UpdateQmatCodon, where the total rate is summed
//across the matrix to calc the blens scaler.  Here it sums the rates for S and NS cells separately
//and returns a vector with (S rate sum) / ((S rate sum) + (NS rate sum)) for each w set and then
//over all categories
void Model::CalcSynonymousBranchlengthProportions(vector<FLOAT_TYPE> &results){
	results.clear();
	UpdateQMatCodon();

//calc the S and NS blens separately
	vector<double> sumS, sumNS;
	sumS.resize(NRateCats());
	sumNS.resize(NRateCats());	
	
	double weightedSumS, weightedSumNS;
	double tempSumS, tempSumNS;
	for(int w=0;w<NRateCats();w++){
		weightedSumS = weightedSumNS = 0.0;
		for(int x=0;x<nstates;x++){
			tempSumS = tempSumNS = 0.0;
			for(int y=0;y<nstates;y++){
				if(x!=y){
					if(qmatLookup[x*nstates+y]&4)
						tempSumNS += qmat[w][x][y];
					else
						tempSumS += qmat[w][x][y];
					}
				}
			//qmat[w][x][x]=-sum;
			weightedSumS += tempSumS * *stateFreqs[x];
			weightedSumNS += tempSumNS * *stateFreqs[x];
			}
		sumS[w] = weightedSumS * *omegaProbs[w];
		sumNS[w] = weightedSumNS * *omegaProbs[w];
		results.push_back((sumS[w] / (sumS[w] + sumNS[w])));
		}
		
	double totSumS = 0.0, totSumNS = 0.0;
	for(int w=0;w<NRateCats();w++){
		totSumS += sumS[w];
		totSumNS += sumNS[w];
		}
	//verify that this all makes sense given the already calc'ed blen mults 
	assert(FloatingPointEquals(blen_multiplier[0], (ONE_POINT_ZERO / (totSumS + totSumNS)), 1e-3));	
	//outman.UserMessage("w = %f S = %f NS = %f, propS = %f", *omegas[0], totSumS, totSumNS, (totSumS / (totSumS + totSumNS)));
	results.push_back(totSumS / (totSumS + totSumNS));
	}

void Model::UpdateQMatAminoAcid(){

	for(int from=0;from<20;from++)
		for(int to=0;to<20;to++)
			qmat[0][from][to] = *stateFreqs[to];

	if(modSpec.IsJonesAAMatrix()) MultiplyByJonesAAMatrix();
	else if(modSpec.IsDayhoffAAMatrix()) MultiplyByDayhoffAAMatrix();
	else if(modSpec.IsWAGAAMatrix()) MultiplyByWAGAAMatrix();
	else if(modSpec.IsMtMamAAMatrix()) MultiplyByMtMamAAMatrix();
	else if(modSpec.IsMtRevAAMatrix()) MultiplyByMtRevAAMatrix();
	else if(modSpec.IsEstimateAAMatrix() || modSpec.IsUserSpecifiedRateMatrix()){
		vector<FLOAT_TYPE *>::iterator r = relNucRates.begin();
		for(int from=0;from<19;from++){
			for(int to=from+1;to<20;to++){
				qmat[0][from][to] *= **r;
				r++;
				}
			}
		assert(r == relNucRates.end());
		r = relNucRates.begin();
		for(int to=0;to<19;to++){
			for(int from=to+1;from<20;from++){
				qmat[0][from][to] *= **r;
				r++;
				}
			}
		assert(r == relNucRates.end());
		}
	
	//set diags to sum rows to 0 and calculate the branch length rescaling factor
	double sum, weightedDiagSum = 0.0;
	blen_multiplier[0] = 0.0;

	for(int from=0;from<20;from++){
		//qmat[0][from][from] = 0.0;
		sum = 0.0;
		for(int to=0;to<20;to++){
			if(from != to) sum += qmat[0][from][to];
			}
		qmat[0][from][from] = -sum;
		weightedDiagSum += sum * *stateFreqs[from];
		}
	blen_multiplier[0] = ONE_POINT_ZERO / weightedDiagSum;
	}

void Model::CalcEigenStuff(){
	ProfCalcEigen.Start();
	//if rate params or statefreqs have been altered, requiring the recalculation of the eigenvectors and c_ijk
	//NOTE that the calculation of the blen_multiplier (rate matrix scaler) now occurs in UpdateQMat()
	UpdateQMat();
	
	int effectiveModels = modSpec.IsNonsynonymousRateHet() ? NRateCats() : 1;
	memcpy(**tempqmat, **qmat, effectiveModels*nstates*nstates*sizeof(MODEL_FLOAT));
	for(int m=0;m<effectiveModels;m++){
		EigenRealGeneral(nstates, tempqmat[m], &eigvals[m][0], eigvalsimag, eigvecs[m], iwork, work);

		memcpy(*teigvecs, *eigvecs[m], nstates*nstates*sizeof(MODEL_FLOAT));
		InvertMatrix(teigvecs, nstates, col, indx, inveigvecs[m]);
		
		//DEBUG - rescaling the eigenvals by blen_multiplier here _should_ give same result as using it later
		//in transition matrix calculating function.  This will be easier for passing to beagle prescaled as well
		for(int ev = 0;ev < nstates;ev++)
			eigvals[m][ev] *= blen_multiplier[m];
		blen_multiplier[m] = 1.0;

		//For codon models using this precalculation actually makes things things slower in CalcPmat (cache thrashing,
		//I think) so don't bother doing it here.  In fact, don't even allocate it in the model
		if(modSpec.IsCodon() == false)
			CalcCijk(&c_ijk[m][0], nstates, (const MODEL_FLOAT**) eigvecs[m], (const MODEL_FLOAT**) inveigvecs[m]);
		}

	eigenDirty=false;
	ProfCalcEigen.Stop();
	}

//this just copies elements from a double precision matrix into a single precision one
void ChangeMatrixPrecision(int elements, double ***pmat, float ***fpmat){
	for(int e=0;e<elements;e++)
		fpmat[0][0][e] = (float) pmat[0][0][e];
	}

//usually this will be called with 2 branch lengths to caluclate two pmats, but if only one 
//is needed the other blen with be -1
void Model::CalcPmats(FLOAT_TYPE blen1, FLOAT_TYPE blen2, FLOAT_TYPE *&mat1, FLOAT_TYPE *&mat2){
	ProfCalcPmat.Start();
//	if(NStates() > 4){
		if(!(blen1 < ZERO_POINT_ZERO)){
			AltCalcPmat(blen1, pmat1);
#ifdef SINGLE_PRECISION_FLOATS
			ChangeMatrixPrecision(modSpec.nstates * modSpec.nstates * modSpec.numRateCats, pmat1, fpmat1);
			mat1 = **fpmat1;
#else
			mat1 = **pmat1;
#endif
			}
		if(!(blen2 < ZERO_POINT_ZERO)){
			AltCalcPmat(blen2, pmat2);
#ifdef SINGLE_PRECISION_FLOATS
			ChangeMatrixPrecision(modSpec.nstates * modSpec.nstates * modSpec.numRateCats, pmat2, fpmat2);
			mat2 = **fpmat2;
#else
			mat2 = **pmat2;
#endif
			}
//		}

/*		for(int i=0;i<nstates;i++)
		for(int j=0;j<nstates;j++)
			assert(FloatingPointEquals(metaPmat[i*nstates+j], pmat[0][i][j], 1e-5));
*/
	ProfCalcPmat.Stop();
	return;
	
	}

//usually this will be called with 2 branch lengths to caluclate two pmats, but if only one 
//is needed the other blen with be -1
//NOTE: unlike the previous version of this, the arguments mat1 and mat2 are the actual destination matrices that
//should be used to hold the resulting pmats, NOT addresses that are to be returned aliased to the two pmats
//that are data members of Model
void Model::CalcPmatsInProvidedMatrices(FLOAT_TYPE blen1, FLOAT_TYPE blen2, FLOAT_TYPE ***mat1, FLOAT_TYPE ***mat2){
	ProfCalcPmat.Start();
//	if(NStates() > 4){
		if(!(blen1 < ZERO_POINT_ZERO)){
			//AltCalcPmat(blen1, pmat1);
			AltCalcPmat(blen1, mat1);
/*
#ifdef SINGLE_PRECISION_FLOATS
			ChangeMatrixPrecision(modSpec.nstates * modSpec.nstates * modSpec.numRateCats, pmat1, fpmat1);
			mat1 = **fpmat1;
#else
			mat1 = **pmat1;
#endif
*/
			}
		if(!(blen2 < ZERO_POINT_ZERO)){
			//AltCalcPmat(blen2, p+mat2);
			AltCalcPmat(blen2, mat2);
/*
#ifdef SINGLE_PRECISION_FLOATS
			ChangeMatrixPrecision(modSpec.nstates * modSpec.nstates * modSpec.numRateCats, pmat2, fpmat2);
			mat2 = **fpmat2;
#else
			mat2 = **pmat2;
#endif
*/
			}
//		}

/*		for(int i=0;i<nstates;i++)
		for(int j=0;j<nstates;j++)
			assert(FloatingPointEquals(metaPmat[i*nstates+j], pmat[0][i][j], 1e-5));
*/
	ProfCalcPmat.Stop();
	return;
	
	}

void Model::CalcPmat(MODEL_FLOAT blen, MODEL_FLOAT *metaPmat, bool flip /*=false*/){
	assert(0);
	/*
	ProfCalcPmat.Start();
	assert(flip == false);

	//this is a bit of a hack to avoid requiring the fuction calling this one to know if 
	//this is a nucleotide, AA or codon model
	if(NStates() > 4){
		if(NStates() == 20){
			CalcPmatNState(blen, metaPmat);
			}
		else{
			FLOAT_TYPE ***ptr;
			AltCalcPmat(blen, ptr);
			memcpy(metaPmat, **ptr, nstates*nstates*NRateCats()*sizeof(FLOAT_TYPE));
			}

		ProfCalcPmat.Stop();
		return;
		}

	//this will be a wacky pmat calculation that combines the pmats for all of the rates
	FLOAT_TYPE tmpFreqs[4];
	for(int i=0;i<nstates;i++) tmpFreqs[i] = *stateFreqs[i];

	for(int r=0;r<NRateCats();r++){
		if(nst==6){
			if(eigenDirty==true)
				CalcEigenStuff();

			FLOAT_TYPE tempblen;
			if(NoPinvInModel()==true || modSpec.IsFlexRateHet())//if we're using flex rates, pinv should already be included
				//in the rate normalization, and doesn't need to be figured in here
				tempblen=(blen * blen_multiplier[0] * rateMults[r]);
			else
				tempblen=(blen * blen_multiplier[0] * rateMults[r]) / (ONE_POINT_ZERO-*propInvar);
			
			CalcPij(c_ijk[0], nstates, eigvals[0], 1, tempblen, pmat[0], EigValexp);
			}
		else if(nst==2 || modSpec.IsEqualStateFrequencies() == false){
			//remember that relNucRates[1] is kappa for nst=2 models
			FLOAT_TYPE PI, A, K=*relNucRates[1];
			FLOAT_TYPE R=tmpFreqs[0]+tmpFreqs[2];
			FLOAT_TYPE Y=ONE_POINT_ZERO - R;
			blen_multiplier[0]=(ZERO_POINT_FIVE/((R*Y)+K*((tmpFreqs[0])*((tmpFreqs[2]))+(tmpFreqs[1])*((tmpFreqs[3])))));
			FLOAT_TYPE tempblen ;
			if(NoPinvInModel()==true || modSpec.IsFlexRateHet())//if we're using flex rates, pinv should already be included
				//in the rate normalization, and doesn't need to be figured in here
				tempblen=(blen * blen_multiplier[0] * rateMults[r]);
			else
				tempblen=(blen * blen_multiplier[0] * rateMults[r]) / (ONE_POINT_ZERO-*propInvar);
			FLOAT_TYPE expblen=exp(-tempblen);

			for(register int f=0;f<4;f++){
				for(register int t=0;t<4;t++){	
					if(f==t){
						if(t==0||t==2) PI = R;
						else PI = Y;
						A=ONE_POINT_ZERO + PI * (K - ONE_POINT_ZERO);
						(**pmat)[f*4+t]=(tmpFreqs[t])+(tmpFreqs[t])*((ONE_POINT_ZERO/PI)-ONE_POINT_ZERO)*expblen+((PI-(tmpFreqs[t]))/PI)*exp(-A*tempblen);
						assert((**pmat)[f*4+t] > ZERO_POINT_ZERO);
						assert((**pmat)[f*4+t] < ONE_POINT_ZERO);
						}
					else if((f+t)%2){
						(**pmat)[f*4+t]=((tmpFreqs[t]))*(ONE_POINT_ZERO-expblen);//tranversion
						assert((**pmat)[f*4+t] > ZERO_POINT_ZERO);
						assert((**pmat)[f*4+t] < ONE_POINT_ZERO);
						}
					else{
						if(t==0||t==2) PI=R;
						else PI = Y;
						A=ONE_POINT_ZERO + PI * (K - ONE_POINT_ZERO);
						(**pmat)[f*4+t]=(tmpFreqs[t])+(tmpFreqs[t])*((ONE_POINT_ZERO/PI)-ONE_POINT_ZERO)*expblen-((tmpFreqs[t])/PI)*exp(-A*tempblen);//transition
						assert((**pmat)[f*4+t] > ZERO_POINT_ZERO);
						assert((**pmat)[f*4+t] < ONE_POINT_ZERO);
						}
					}
				}
			}
		else if(nst==1){
			blen_multiplier[0]=(FLOAT_TYPE)(4.0/3.0);
			//	}
			FLOAT_TYPE tempblen ;
			if(NoPinvInModel()==true || modSpec.IsFlexRateHet())//if we're using flex rates, pinv should already be included
				//in the rate normalization, and doesn't need to be figured in here
				tempblen=(blen * blen_multiplier[0] * rateMults[r]);
			else
				tempblen=(blen * blen_multiplier[0] * rateMults[r]) / (ONE_POINT_ZERO-*propInvar);
			FLOAT_TYPE expblen=exp(-tempblen);			
			for(register int f=0;f<4;f++){
				for(register int t=0;t<4;t++){
					if(f==t)
						(**pmat)[f*4+t]=expblen+(FLOAT_TYPE) 0.25*(ONE_POINT_ZERO-expblen);
					else 
						(**pmat)[f*4+t]=(FLOAT_TYPE)0.25*((ONE_POINT_ZERO)-expblen);
					}
				}
			}
	
		if(flip==true){
			//copy and flip the calculated pmat into the metaPmat
			for(int i=0;i<4;i++)
				for(int j=0;j<4;j++){
					metaPmat[i*16+j+r*4]=pmat[0][0][i+j*4];
					}
			}
		else{
			//Copy the pmats into the metaPmat in order
			for(int i=0;i<4;i++)
				for(int j=0;j<4;j++){
					metaPmat[i*4+j+r*16]=pmat[0][0][i*4+j];
					}
			}
		}	
	ProfCalcPmat.Stop();
*/
	}	

void Model::CalcPmatNState(FLOAT_TYPE blen, MODEL_FLOAT *metaPmat){
	assert(0);
/*
	if(eigenDirty==true)
		CalcEigenStuff();

	if(modSpec.IsNonsynonymousRateHet()){
		for(int w=0;w<NRateCats();w++){
			FLOAT_TYPE tempblen;
			tempblen=blen * blen_multiplier[w];

			CalcPij(&c_ijk[w][0], nstates, &eigvals[w][0], 1, tempblen, pmat[0], EigValexp);

			//Copy the pmats into the metaPmat in order
			for(int i=0;i<nstates;i++)
				for(int j=0;j<nstates;j++)
					metaPmat[w*nstates*nstates + i*nstates + j]=pmat[0][i][j];
			}
		}
	else{
		for(int r=0;r<NRateCats();r++){
			FLOAT_TYPE tempblen;
			if(NoPinvInModel()==true || modSpec.IsFlexRateHet())//if we're using flex rates, pinv should already be included
				//in the rate normalization, and doesn't need to be figured in here
				tempblen=(blen * blen_multiplier[0] * rateMults[r]);
			else
				tempblen=(blen * blen_multiplier[0] * rateMults[r]) / (ONE_POINT_ZERO-*propInvar);

			CalcPij(c_ijk[0], nstates, eigvals[0], 1, tempblen, pmat[0], EigValexp);

			//Copy the pmats into the metaPmat in order
			for(int i=0;i<nstates;i++)
				for(int j=0;j<nstates;j++)
					metaPmat[r*nstates*nstates + i*nstates + j]=pmat[0][0][i*nstates + j];
			}
		}
*/
/*	char filename[50];
	char temp[10];
	sprintf(temp, "%.2f.log", blen);	
	strcpy(filename, "pmatdebug");
	strcat(filename, temp);	
	*/
	
	//output pmat
/*	ofstream pmd(filename);
	for(int i=0;i<nstates;i++){
		for(int j=0;j<nstates;j++){
			pmd << p[i][j] << "\t";
			}
		pmd << endl;
		}
*/	}


void Model::OutputPmats(ofstream &deb){
	for(int r=0;r<NRateCats();r++){
		deb << "pmat1 rate" << r << endl;
		for(int f=0;f<nstates;f++){
			for(int t=0;t<nstates;t++){
				deb << pmat1[r][f][t] << "\t";
				}
			deb << endl;
			}
		deb << endl;
		}
	deb << endl;
	for(int r=0;r<NRateCats();r++){
		deb << "pmat2 rate" << r << endl;
		for(int f=0;f<nstates;f++){
			for(int t=0;t<nstates;t++){
				deb << pmat2[r][f][t] << "\t";
				}
			deb << endl;
			}
		deb << endl;
		}

	}

void Model::CalcDerivatives(FLOAT_TYPE dlen, FLOAT_TYPE ***&pr, FLOAT_TYPE ***&one, FLOAT_TYPE ***&two){
/*	double before = *omegas[0];
	if(dlen < 0.011 && dlen > 0.009)
		SetOmega(0, 1.0);
*/
	if(eigenDirty==true)
		CalcEigenStuff();
/*
	if(dlen < 0.011 && dlen > 0.009)
		SetOmega(0, before);
*/
	for(int rate=0;rate<NRateCats();rate++){
		const unsigned rateOffset = nstates*rate; 
		for(int k=0; k<nstates; k++){
			MODEL_FLOAT scaledEigVal;
			if(modSpec.IsNonsynonymousRateHet() == false){
				if(NoPinvInModel()==true || modSpec.IsFlexRateHet())//if we're using flex rates, pinv should already be included
					//in the rate normalization, and doesn't need to be figured in here
					scaledEigVal = eigvals[0][k]*rateMults[rate]*blen_multiplier[0];	
				else
					scaledEigVal = eigvals[0][k]*rateMults[rate]*blen_multiplier[0]/(ONE_POINT_ZERO-*propInvar);
				}
			else{
				scaledEigVal = eigvals[rate][k]*blen_multiplier[rate];
				}
			EigValexp[k+rateOffset] = exp(scaledEigVal * dlen);
			EigValderiv[k+rateOffset] = scaledEigVal*EigValexp[k+rateOffset];
			EigValderiv2[k+rateOffset] = scaledEigVal*EigValderiv[k+rateOffset];
			}
		}

	if(NStates() > 59){//using precalced eigvecs X inveigvecs (c_ijk) is less efficient for codon models, and I
					//don't want a conditional in the inner loop
		for(int rate=0;rate<NRateCats();rate++){
			int model=0;
			if(modSpec.IsNonsynonymousRateHet())
				model = rate;
			const unsigned rateOffset = nstates*rate;
			for (int i = 0; i < nstates; i++){
				for (int j = 0; j < nstates; j++){
					MODEL_FLOAT sum_p=ZERO_POINT_ZERO;
					MODEL_FLOAT sum_d1p=ZERO_POINT_ZERO;
					MODEL_FLOAT sum_d2p = ZERO_POINT_ZERO;
					for (int k = 0; k < nstates; k++){ 
						const MODEL_FLOAT x = eigvecs[model][i][k]*inveigvecs[model][k][j];
						sum_p   += x*EigValexp[k+rateOffset];
						sum_d1p += x*EigValderiv[k+rateOffset];
						sum_d2p += x*EigValderiv2[k+rateOffset];
						}
					pmat1[rate][i][j] = (sum_p > ZERO_POINT_ZERO ? sum_p : ZERO_POINT_ZERO);
					deriv1[rate][i][j] = sum_d1p;
					deriv2[rate][i][j] = sum_d2p;
					}
				}
			}
		}
	else{ // aminoacids or nucleotides
		for(int rate=0;rate<NRateCats();rate++){
			const unsigned rateOffset = nstates*rate;
			for (int i = 0; i < nstates; i++){
				for (int j = 0; j < nstates; j++){
					MODEL_FLOAT sum_p=ZERO_POINT_ZERO;
					MODEL_FLOAT sum_d1p=ZERO_POINT_ZERO;
					MODEL_FLOAT sum_d2p = ZERO_POINT_ZERO;
					for (int k = 0; k < nstates; k++){ 
						MODEL_FLOAT x = c_ijk[0][i*nstates*nstates + j*nstates +k];
						sum_p   += x*EigValexp[k+rateOffset];
						sum_d1p += x*EigValderiv[k+rateOffset];
						sum_d2p += x*EigValderiv2[k+rateOffset];
						}
					pmat1[rate][i][j] = (sum_p > ZERO_POINT_ZERO ? sum_p : ZERO_POINT_ZERO);
					deriv1[rate][i][j] = sum_d1p;
					deriv2[rate][i][j] = sum_d2p;
					}
				}
			}
		}
#ifdef SINGLE_PRECISION_FLOATS
	ChangeMatrixPrecision(nstates * nstates * NRateCats(), deriv1, fderiv1);
	ChangeMatrixPrecision(nstates * nstates * NRateCats(), deriv2, fderiv2);
	ChangeMatrixPrecision(nstates * nstates * NRateCats(), pmat1, fpmat1);

	one=fderiv1;
	two=fderiv2;
	pr=fpmat1;
#else
	one=deriv1;
	two=deriv2;
	pr=pmat1;
#endif
	}


bool DoubleAbsLessThan(double &first, double &sec){return fabs(first) <= fabs(sec);}

void Model::AltCalcPmat(FLOAT_TYPE dlen, MODEL_FLOAT ***&pmat){
/*	double before = *omegas[0];
	if(dlen < 0.011 && dlen > 0.009)
		SetOmega(0, 1.0);
*/
	if(eigenDirty==true)
		CalcEigenStuff();
/*
	if(dlen < 0.011 && dlen > 0.009)
		SetOmega(0, before);
*/
	for(int rate=0;rate<NRateCats();rate++){
		const unsigned rateOffset = nstates*rate; 
		for(int k=0; k<nstates; k++){
			MODEL_FLOAT scaledEigVal;
			if(modSpec.IsNonsynonymousRateHet() == false){
				if(NoPinvInModel()==true || modSpec.IsFlexRateHet())//if we're using flex rates, pinv should already be included
					//in the rate normalization, and doesn't need to be figured in here
					scaledEigVal = eigvals[0][k]*rateMults[rate]*blen_multiplier[0];	
				else
					scaledEigVal = eigvals[0][k]*rateMults[rate]*blen_multiplier[0]/(ONE_POINT_ZERO-*propInvar);
				}
			else{
				scaledEigVal = eigvals[rate][k]*blen_multiplier[rate];
				}
			EigValexp[k+rateOffset] = exp(scaledEigVal * dlen);
			}
		}

	if(NStates() == 20){
		for(int rate=0;rate<NRateCats();rate++){
			int model=0;
			const unsigned rateOffset = 20*rate;
			for (int i = 0; i < 20; i++){
				for (int j = 0; j < 20; j++){
					MODEL_FLOAT sum_p=ZERO_POINT_ZERO;
					for (int k = 0; k < 20; k++){ 
						const MODEL_FLOAT x = c_ijk[0][model*20*20*20 + i*20*20 + j*20 +k];
						sum_p   += x*EigValexp[k+rateOffset];
						}
					pmat[rate][i][j] = (sum_p > ZERO_POINT_ZERO ? sum_p : ZERO_POINT_ZERO);
					}
				}
			}
		}
	else if(NStates()>59){
		for(int rate=0;rate<NRateCats();rate++){
			int model=0;
			if(modSpec.IsNonsynonymousRateHet())
				model = rate;
			const unsigned rateOffset = nstates*rate;
			for (int i = 0; i < nstates; i++){
				for (int j = 0; j < nstates; j++){
					MODEL_FLOAT sum_p=ZERO_POINT_ZERO;

/*					
					FLOAT_TYPE sum_pBig=ZERO_POINT_ZERO;
					FLOAT_TYPE sum_pSmall=ZERO_POINT_ZERO;
					FLOAT_TYPE sum_pBig2=ZERO_POINT_ZERO;
					FLOAT_TYPE sum_pSmall2=ZERO_POINT_ZERO;
					for (int k = 0; k < nstates; k++){ 
						const FLOAT_TYPE x = eigvecs[model][i][k]*inveigvecs[model][k][j];
					
						if(x < ZERO_POINT_ZERO){
							if(x < -1e-4)
								sum_pSmall   += x*EigValexp[k+rateOffset];
							else 
								sum_pSmall2   += x*EigValexp[k+rateOffset];
							}
						else{
							if(x > 1e-4)
								sum_pBig   += x*EigValexp[k+rateOffset];
							else
								sum_pBig2   += x*EigValexp[k+rateOffset];
							}
						}
//					FLOAT_TYPE tot = sum_pBig2 + sum_pSmall2; 
//					tot += sum_pBig + sum_pSmall;
					FLOAT_TYPE tot = sum_pBig2 + sum_pBig; 
					tot += (sum_pSmall2 + sum_pSmall);
					sum_p = tot;
*/
					for (int k = 0; k < nstates; k++){ 
						const MODEL_FLOAT x = eigvecs[model][i][k]*inveigvecs[model][k][j];
						sum_p   += x*EigValexp[k+rateOffset];
						}
					pmat[rate][i][j] = (sum_p > ZERO_POINT_ZERO ? sum_p : ZERO_POINT_ZERO);
					}
				}
			}
		}
	else{
		for(int rate=0;rate<NRateCats();rate++){
			int model=0;
			const unsigned rateOffset = 4*rate;
			for (int i = 0; i < 4; i++){
				for (int j = 0; j < 4; j++){
					MODEL_FLOAT sum_p=ZERO_POINT_ZERO;
					for (int k = 0; k < 4; k++){ 
						const MODEL_FLOAT x = c_ijk[0][model*4*4*4 + i*4*4 + j*4 +k];
						sum_p   += x*EigValexp[k+rateOffset];
						}
					pmat[rate][i][j] = (sum_p > ZERO_POINT_ZERO ? sum_p : ZERO_POINT_ZERO);
					}
				}
			}
		}
	}

void Model::SetDefaultModelParameters(const SequenceData *data){
	//some of these depend on having read the data already
	//also note that this resets the values in the case of 
	//bootstrapping.  Any of this could be overridden by
	//values specified in a start file

	for(vector<BaseParameter*>::iterator pit=paramsToMutate.begin();pit != paramsToMutate.end();pit++){
		(*pit)->SetToDefaultValues();
		}
	if(modSpec.numRateCats > 1 && modSpec.IsNonsynonymousRateHet() == false){
		if(modSpec.IsFlexRateHet()){
			//if alpha is only being used to manipulate the flex rates, it wouldn't be reset above
			SetAlpha(0, 0.5);
			}
		DiscreteGamma(rateMults, rateProbs, *alpha);
		}

	if((modSpec.IsEqualStateFrequencies() == false && (modSpec.IsCodon() && modSpec.IsUserSpecifiedStateFrequencies()) == false && modSpec.IsDayhoffAAFreqs() == false && modSpec.IsWAGAAFreqs() == false && modSpec.IsJonesAAFreqs() == false && modSpec.IsMtMamAAFreqs() == false && modSpec.IsMtRevAAFreqs() == false)
		|| (modSpec.IsF3x4StateFrequencies() || modSpec.IsF1x4StateFrequencies())){
		//if the state freqs aren't equal, they will either start at the empirical values 
		//or be fixed at them
		//if using the F3x4 or F1x4 flavors, they should have already be calculated and stored in the data empirical frequency field
		FLOAT_TYPE *f = new FLOAT_TYPE[modSpec.nstates];
		data->GetEmpiricalFreqs(f);
		SetPis(f, false, true);
		delete []f;
		}

	if(modSpec.includeInvariantSites==false){
		SetPinv(ZERO_POINT_ZERO, false);
		SetMaxPinv(ZERO_POINT_ZERO);
		}
	else{
		//if there are no constant sites, warn user that Pinv should not be used
		//if(data->NConstant() == 0) throw(ErrorException("This dataset contains no constant characters!\nInference of the proportion of invariant sites is therefore meaningless.\nPlease set invariantsites to \"none\""));
		if(data->NConstant() == 0){
			outman.UserMessage("This dataset contains no constant characters!\nInference of the proportion of invariant sites is therefore meaningless.\nSetting invariantsites to \"none\".");
			SetPinv(ZERO_POINT_ZERO, false);
			SetMaxPinv(ZERO_POINT_ZERO);
			modSpec.includeInvariantSites = false;
			}
		else{
			SetPinv((FLOAT_TYPE)0.25 * ((FLOAT_TYPE)data->NConstant()/(data->NConstant()+data->NInformative()+data->NVarUninform())), false);
			SetMaxPinv((FLOAT_TYPE)data->NConstant()/(data->NConstant()+data->NInformative()+data->NVarUninform()));
			if(modSpec.IsFlexRateHet()) 
				NormalizeRates();
			else AdjustRateProportions();
			}
		}
	eigenDirty = true;
	}

void Model::MutateRates(){
	//paramsToMutate[1]->Mutator(Model::mutationShape);
	//assert(Rates(0) == Rates(2));

/*	int rateToChange=int(rnd.uniform()*(nst));
	
	if(rateToChange<nst-1){
		rates[rateToChange] *= rnd.gamma( Model::mutationShape );
//		rates[rateToChange] *= exp(MODEL_CHANGE_SCALER * (params->rnd.uniform()-.5));
		if(rates[rateToChange]>99.9) rates[rateToChange]=99.9;
		}

	else{//if we alter the reference rate GT (fixed to 1.0)
		//scale all of the other rates
		//FLOAT_TYPE scaler=exp(MODEL_CHANGE_SCALER * (params->rnd.uniform()-.5));
		FLOAT_TYPE scaler= rnd.gamma( Model::mutationShape );
		for(int i=0;i<nst-1;i++){
			rates[i] /= scaler;
			}
		}

	// don't let rates[0] become greater than 99.0
	// if this upper limit is changed, be sure to check consequences
	// in AllocMigrantStrings function in gamlmain.C (i.e., currently
	// only allow 3 characters plus number needed for precision)
*/	eigenDirty=true;
	}

void Model::MutatePis(){
	assert(0);
	//basetest->Mutator(Model::mutationShape);
//	paramsToMutate[0]->Mutator(Model::mutationShape);
//	dirty=true;

	//alternative:change one pi with a multiplier and rescale the rest
/*	int piToChange=int(rnd.uniform()*4.0);
	
	FLOAT_TYPE newPi=pi[piToChange] * rnd.gamma( Model::mutationShape );
	for(int b=0;b<4;b++)
		if(b!=piToChange) pi[b] *= (1.0-newPi)/(1.0-pi[piToChange]);
	pi[piToChange]=newPi;
	dirty=true;
*/	}
/*
void Model::MutateRateProbs(){
	int ProbToChange=int(rnd.uniform()*(FLOAT_TYPE) NRateCats());
	
	FLOAT_TYPE newProb=rateProbs[ProbToChange] * rnd.gamma( Model::mutationShape / 10.0 );
	for(int b=0;b<NRateCats();b++)
		if(b!=ProbToChange) rateProbs[b] *= (1.0-newProb)/(1.0-rateProbs[ProbToChange]);
	rateProbs[ProbToChange]=newProb;
	NormalizeRates();
	}

void Model::MutateRateMults(){
	int rateToChange=int(rnd.uniform()*NRateCats());
	rateMults[rateToChange] *= rnd.gamma( Model::mutationShape / 10.0);
	NormalizeRates();
	}
	
void Model::MutateAlpha(){
//	alpha *= exp(MODEL_CHANGE_SCALER * (params->rnd.uniform()-.5));
	*alpha *=rnd.gamma( Model::mutationShape );
	DiscreteGamma(rateMults, rateProbs, *alpha);
	//change the proportion of rates in each gamma cat
	}
	
void Model::MutatePropInvar(){
//	propInvar *= exp(MODEL_CHANGE_SCALER * (params->rnd.uniform()-.5));
	FLOAT_TYPE mult=rnd.gamma( Model::mutationShape );
	if(*propInvar == maxPropInvar && (mult > 1.0)) mult=1.0/mult;
	*propInvar *= mult;
	*propInvar = (*propInvar > maxPropInvar ? maxPropInvar : *propInvar);
	//change the proportion of rates in each gamma cat
	for(int i=0;i<NRateCats();i++){
		rateProbs[i]=(1.0-*propInvar)/NRateCats();
		}
	}
*/
void Model::CopyModel(const Model *from){
	if(modSpec.IsCodon()){
		for(int i=0;i<omegas.size();i++)
			*omegas[i]=*(from->omegas[i]);
		for(int i=0;i<omegaProbs.size();i++)
			*omegaProbs[i]=*(from->omegaProbs[i]);
		}

	if(modSpec.IsAminoAcid() == false || modSpec.IsEstimateAAMatrix() || (modSpec.IsAminoAcid() && modSpec.IsUserSpecifiedRateMatrix()))
		for(int i=0;i<relNucRates.size();i++)
			*relNucRates[i]=*(from->relNucRates[i]);
	
	for(int i=0;i<nstates;i++)
		*stateFreqs[i]=*(from->stateFreqs[i]);

	//memcpy(pi, from->pi, sizeof(FLOAT_TYPE)*4);

	memcpy(rateMults, from->rateMults, sizeof(FLOAT_TYPE)*NRateCats());
	memcpy(rateProbs, from->rateProbs, sizeof(FLOAT_TYPE)*NRateCats());

	if(modSpec.IsGammaRateHet())
		*alpha=*(from->alpha);
	*propInvar=*(from->propInvar);

	if(from->eigenDirty == false){
		//copy the already calculated eigen variables, which are nontrivial to 
		//calculate for non-nucleotide models
		CopyEigenVariables(from);
		eigenDirty = false;
		}
	else 
		eigenDirty=true;
	}	

void Model::CopyEigenVariables(const Model *from){
	int effectiveModels = modSpec.IsNonsynonymousRateHet() ? NRateCats() : 1;
	memcpy(**qmat, **from->qmat, effectiveModels*nstates*nstates*sizeof(MODEL_FLOAT));
	memcpy(**eigvecs, **from->eigvecs, NRateCats()*nstates*nstates*sizeof(MODEL_FLOAT));
	memcpy(**inveigvecs, **from->inveigvecs, NRateCats()*nstates*nstates*sizeof(MODEL_FLOAT));
	memcpy(*eigvals, *from->eigvals, effectiveModels * nstates * sizeof(MODEL_FLOAT));
	memcpy(blen_multiplier, from->blen_multiplier, effectiveModels * sizeof(FLOAT_TYPE));
	//c_ijk isn't allocated or used for codon models
	if(c_ijk != NULL)
		memcpy(*c_ijk, *from->c_ijk, effectiveModels*nstates*nstates*nstates*sizeof(MODEL_FLOAT));	
	}

void Model::SetModel(FLOAT_TYPE *model_string){
	int slot=0;
	for(int i=0;i<nst-1;i++)
		*relNucRates[i]=model_string[slot++];
	for(int j=0;j<4;j++)
		*stateFreqs[j]=model_string[slot++];
		
	if(NRateCats()>1) *alpha=model_string[slot++];
	DiscreteGamma(rateMults, rateProbs, *alpha);
	//using whether or not this individual had a PI of >0 in the first
	//place to decide whether we should expect one in the string.
	//Seems safe.
	if(*propInvar!=ZERO_POINT_ZERO) *propInvar=model_string[slot++];
	eigenDirty=true;
	}

FLOAT_TYPE Model::TRatio() const{
	FLOAT_TYPE numerator = *relNucRates[1] * ( *stateFreqs[0]**stateFreqs[2] + *stateFreqs[1]**stateFreqs[3] );
	FLOAT_TYPE denominator = ( *stateFreqs[0] + *stateFreqs[2] ) * ( *stateFreqs[1] + *stateFreqs[3] );
	return ( numerator / denominator );
	}

bool Model::IsModelEqual(const Model *other) const {
	assert(0);
	//this will need to be generalized if other models are introduced
	for(int i=0;i<6;i++)
		if(!FloatingPointEquals(*relNucRates[i], *(other->relNucRates[i]), 1e-15)) return false;
	
	for(int i=0;i<nstates;i++)
		if(!FloatingPointEquals(*stateFreqs[i], *(other->stateFreqs[i]), 1e-15)) return false;

	if(!modSpec.IsCodon() && NRateCats() > 1){
		for(int i=0;i<this->NRateCats();i++){
			if(!FloatingPointEquals(rateMults[i], other->rateMults[i], 1e-15)) return false;
			if(!FloatingPointEquals(rateProbs[i], other->rateProbs[i], 1e-15)) return false;
			}
		}
	else if(modSpec.IsCodon()){
		for(int i=0;i<this->NRateCats();i++){
			if(!FloatingPointEquals(Omega(i), other->Omega(i), 1e-15)) return false;
			if(!FloatingPointEquals(OmegaProb(i), other->OmegaProb(i), 1e-15)) return false;
			}
		}

/*
	if(rateMults[0] != other->rateMults[0]) return false;
	if(rateMults[1] != other->rateMults[1]) return false;
	if(rateMults[2] != other->rateMults[2]) return false;
	if(rateMults[3] != other->rateMults[3]) return false;
	
	if(rateProbs[0] != other->rateProbs[0]) return false;
	if(rateProbs[1] != other->rateProbs[1]) return false;
	if(rateProbs[2] != other->rateProbs[2]) return false;
	if(rateProbs[3] != other->rateProbs[3]) return false;
*/
	if(alpha!=other->alpha) return false;
	if(propInvar!=other->propInvar) return false;

	

	return true;
	}

//a bunch of the gamma rate het machinery
//from MrBayes

/*-------------------------------------------------------------------------------
|                                                                               |
|  Discretization of gamma distribution with equal proportions in each          |
|  category.                                                                    |
|                                                                               |
-------------------------------------------------------------------------------*/ 
#ifdef SINGLE_PRECISION_FLOATS
#define POINTGAMMA(prob,alpha,beta) 		PointChi2(prob,2.0f*(alpha))/(2.0f*(beta))
#else
#define POINTGAMMA(prob,alpha,beta) 		PointChi2(prob,2.0*(alpha))/(2.0*(beta))
#endif

/* ------------------------------------------------------------------------------
|                                                                               |
|  Returns z so That Prob{x<z} = prob where x ~ N(0,1) and                      |
|  (1e-12) < prob < 1-(1e-12).  Returns (-9999) if in error.                    |
|                                                                               |
|  Odeh, R. E. and J. O. Evans.  1974.  The percentage points of the normal     |
|     distribution.  Applied Statistics, 22:96-97 (AS70)                        |
|                                                                               |
|  Newer methods:                                                               |
|                                                                               |
|  Wichura, M. J.  1988.  Algorithm AS 241: The percentage points of the        |
|     normal distribution.  37:477-484.                                         |
|  Beasley, JD & S. G. Springer.  1977.  Algorithm AS 111: The percentage       |
|     points of the normal distribution.  26:118-121.                           |
|                                                                               |
-------------------------------------------------------------------------------*/   
FLOAT_TYPE PointNormal (FLOAT_TYPE prob){
#ifdef SINGLE_PRECISION_FLOATS
	FLOAT_TYPE 		a0 = -0.322232431088f, a1 = -1.0f, a2 = -0.342242088547f, a3 = -0.0204231210245f,
 					a4 = -0.453642210148e-4f, b0 = 0.0993484626060f, b1 = 0.588581570495f,
 					b2 = 0.531103462366f, b3 = 0.103537752850f, b4 = 0.0038560700634f,
 					y, z = 0.0f, p = prob, p1;
#else
	FLOAT_TYPE 		a0 = -0.322232431088, a1 = -1.0, a2 = -0.342242088547, a3 = -0.0204231210245,
 					a4 = -0.453642210148e-4, b0 = 0.0993484626060, b1 = 0.588581570495,
 					b2 = 0.531103462366, b3 = 0.103537752850, b4 = 0.0038560700634,
 					y, z = 0, p = prob, p1;
#endif

	p1 = (p<ZERO_POINT_FIVE ? p : 1-p);
	if (p1<1e-20) 
	   return (-9999);
	y = sqrt (log(1/(p1*p1)));   
	z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
	return (p<ZERO_POINT_FIVE ? -z : z);

}

/*-------------------------------------------------------------------------------
|                                                                               |
|  Returns the incomplete gamma ratio I(x,alpha) where x is the upper           |
|  limit of the integration and alpha is the shape parameter.  Returns (-1)     |
|  if in error.                                                                 |
|  LnGamma_alpha = ln(Gamma(alpha)), is almost redundant.                      |
|  (1) series expansion     if (alpha>x || x<=1)                                |
|  (2) continued fraction   otherwise                                           |
|                                                                               |
|  RATNEST FORTRAN by                                                           |
|  Bhattacharjee, G. P.  1970.  The incomplete gamma integral.  Applied         |
|     Statistics, 19:285-287 (AS32)                                             |
|                                                                               |
-------------------------------------------------------------------------------*/   
FLOAT_TYPE IncompleteGamma (FLOAT_TYPE x, FLOAT_TYPE alpha, FLOAT_TYPE LnGamma_alpha){
	int 			i;
#ifdef SINGLE_PRECISION_FLOATS
	FLOAT_TYPE 		p = alpha, g = LnGamma_alpha,
					accurate = GARLI_FP_EPS, overflow = 1e30f,
					factor, gin = 0.0f, rn = 0.0f, a = 0.0f, b = 0.0f, an = 0.0f, 
					dif = 0.0f, term = 0.0f, pn[6];
#else
	FLOAT_TYPE 		p = alpha, g = LnGamma_alpha,
					accurate = 1e-8, overflow = 1e30,
					factor, gin = 0.0, rn = 0.0, a = 0.0, b = 0.0, an = 0.0, 
					dif = 0.0, term = 0.0, pn[6];
#endif

	if (x == ZERO_POINT_ZERO) 
		return (ZERO_POINT_ZERO);
	if (x < 0 || p <= 0) 
		return (-ONE_POINT_ZERO);

	factor = exp(p*log(x)-x-g);   
	if (x>1 && x>=p) 
		goto l30;
	gin = ONE_POINT_ZERO;  
	term = ONE_POINT_ZERO;  
	rn = p;
	l20:
		rn++;
		term *= x/rn;   
		gin += term;
		if (term > accurate) 
			goto l20;
		gin *= factor/p;
		goto l50;
	l30:
		a = ONE_POINT_ZERO-p;   
		b = a+x+ONE_POINT_ZERO;  
		term = ZERO_POINT_ZERO;
		pn[0] = ONE_POINT_ZERO;  
		pn[1] = x;  
		pn[2] = x+1;  
		pn[3] = x*b;
		gin = pn[2]/pn[3];
	l32:
		a++;  
		b += 2.0;  
		term++;   
		an = a*term;
		for (i=0; i<2; i++) 
			pn[i+4] = b*pn[i+2]-an*pn[i];
		if (pn[5] == 0) 
			goto l35;
		rn = pn[4]/pn[5];   
		dif = fabs(gin-rn);
		if (dif>accurate) 
			goto l34;
		if (dif<=accurate*rn) 
			goto l42;
	l34:
		gin = rn;
	l35:
		for (i=0; i<4; i++) 
			pn[i] = pn[i+2];
		if (fabs(pn[4]) < overflow) 
			goto l32;
		for (i=0; i<4; i++) 
			pn[i] /= overflow;
		goto l32;
	l42:
		gin = ONE_POINT_ZERO-factor*gin;
	l50:
		return (gin);

}

inline FLOAT_TYPE LnGamma (FLOAT_TYPE alp){
/*	FLOAT_TYPE cof[6];
	cof[0]=76.18009172947146;
    cof[1]=-86.50532032941677;
    cof[2]=24.01409824083091;
    cof[3]=-1.231739572450155;
    cof[4]=0.1208650973866179e-2;
    cof[5]=-0.5395239384953e-5;	
	FLOAT_TYPE xx=alp;
	FLOAT_TYPE yy=alp;
	FLOAT_TYPE tmp=xx + 5.5 - (xx + 0.5) * log(xx + 5.5);
	FLOAT_TYPE ser = 1.000000000190015;
	for(int j=0;j<5;j++){
		ser += (cof[j] / ++yy);
		}
	return log(2.5066282746310005*ser/xx)-tmp;
	}
*/
	FLOAT_TYPE x = alp, f=ZERO_POINT_ZERO, z;
	
	if (x < 7) 
		{
		f = ONE_POINT_ZERO;  
		z = x-ONE_POINT_ZERO;
		while (++z < 7.0)  
			f *= z;
		x = z;   
		f = -log(f);
		}
	z = ONE_POINT_ZERO/(x*x);
#ifdef SINGLE_PRECISION_FLOATS
	return  (f + (x-0.5f)*log(x) - x + 0.918938533204673f + 
			(((-0.000595238095238f*z+0.000793650793651f)*z-0.002777777777778f)*z +
			0.083333333333333f)/x);
#else
	return  (f + (x-0.5)*log(x) - x + 0.918938533204673 + 
			(((-0.000595238095238*z+0.000793650793651)*z-0.002777777777778)*z +
			0.083333333333333)/x);
#endif
	}

FLOAT_TYPE PointChi2 (FLOAT_TYPE prob, FLOAT_TYPE v){
#ifdef SINGLE_PRECISION_FLOATS
	//potential error e needs to be increased here
	//because of lesser decimal precision of floats
	FLOAT_TYPE 		e = 0.5e-4f, aa = 0.6931471805f, p = prob, g,
					xx, c, ch, a = 0.0f, q = 0.0f, p1 = 0.0f, p2 = 0.0f, t = 0.0f, 
					x = 0.0f, b = 0.0f, s1, s2, s3, s4, s5, s6;
	if (p < 0.000002f || p > 0.999998f || v <= 0.0f) 
		return (-1.0f);

	g = LnGamma (v*ZERO_POINT_FIVE);
	xx = v/2.0f;   
	c = xx - ONE_POINT_ZERO;
	if (v >= -1.24f*log(p)) 
		goto l1;
#else
	FLOAT_TYPE 		e = 0.5e-6, aa = 0.6931471805, p = prob, g,
					xx, c, ch, a = 0.0, q = 0.0, p1 = 0.0, p2 = 0.0, t = 0.0, 
					x = 0.0, b = 0.0, s1, s2, s3, s4, s5, s6;
	if (p < 0.000002 || p > 0.999998 || v <= 0.0) 
		return (-ONE_POINT_ZERO);
	
	g = LnGamma (v*ZERO_POINT_FIVE);
	xx = v/2.0;   
	c = xx - ONE_POINT_ZERO;
	if (v >= -1.24*log(p)) 
		goto l1;
#endif

	ch = pow((p*xx*exp(g+xx*aa)), ONE_POINT_ZERO/xx);
	if (ch-e < ZERO_POINT_ZERO) 
		return (ch);
	goto l4;
#ifdef SINGLE_PRECISION_FLOATS
	l1:
		if (v > 0.32f) 
			goto l3;
		ch = 0.4f;
		a = log(ONE_POINT_ZERO-p);
	l2:
		q = ch;  
		p1 = ONE_POINT_ZERO+ch*(4.67f+ch);  
		p2 = ch*(6.73f+ch*(6.66f+ch));
		t = -0.5f+(4.67f+2.0f*ch)/p1 - (6.73f+ch*(13.32f+3.0f*ch))/p2;
		ch -= (ONE_POINT_ZERO-exp(a+g+0.5f*ch+c*aa)*p2/p1)/t;
		if (fabs(q/ch-ONE_POINT_ZERO)-0.01f <= ZERO_POINT_ZERO) 
			goto l4;
		else                       
			goto l2;
	l3: 
		x = PointNormal (p);
		p1 =  0.222222f/v;   
		ch = v*pow((x*sqrt(p1)+ONE_POINT_ZERO-p1), 3.0f);
		if (ch > 2.2f*v+6.0f)  
			ch =  -2.0f*(log(ONE_POINT_ZERO-p)-c*log(0.5f*ch)+g);
#else
	l1:
		if (v > 0.32) 
			goto l3;
		ch = 0.4;
		a = log(ONE_POINT_ZERO-p);
	l2:
		q = ch;  
		p1 = ONE_POINT_ZERO+ch*(4.67+ch);  
		p2 = ch*(6.73+ch*(6.66+ch));
		t = -0.5+(4.67+2.0*ch)/p1 - (6.73+ch*(13.32+3.0*ch))/p2;
		ch -= (ONE_POINT_ZERO-exp(a+g+0.5*ch+c*aa)*p2/p1)/t;
		if (fabs(q/ch-ONE_POINT_ZERO)-0.01 <= ZERO_POINT_ZERO) 
			goto l4;
		else                       
			goto l2;
	l3: 
		x = PointNormal (p);
		p1 =  0.222222/v;   
		ch = v*pow((x*sqrt(p1)+ONE_POINT_ZERO-p1), 3.0);
		if (ch > 2.2*v+6.0)  
			ch =  -2.0*(log(ONE_POINT_ZERO-p)-c*log(0.5*ch)+g);
#endif
	l4:
		q = ch;
		p1 = ZERO_POINT_FIVE*ch;
		if ((t = IncompleteGamma (p1, xx, g)) < ZERO_POINT_ZERO) 
			{
			printf ("\nerr IncompleteGamma");
			return (-ONE_POINT_ZERO);
			}
		p2 = p-t;
		t = p2*exp(xx*aa+g+p1-c*log(ch));   
		b = t/ch;
#ifdef SINGLE_PRECISION_FLOATS
		a = 0.5f*t-b*c;
		s1 = (210.0f+a*(140.0f+a*(105.0f+a*(84.0f+a*(70.0f+60.0f*a))))) / 420.0f;
		s2 = (420.0f+a*(735.0f+a*(966.0f+a*(1141.0f+1278.0f*a))))/2520.0f;
		s3 = (210.0f+a*(462.0f+a*(707.0f+932.0f*a)))/2520.0f;
		s4 = (252.0f+a*(672.0f+1182.0f*a)+c*(294.0f+a*(889.0f+1740.0f*a)))/5040.0f;
		s5 = (84.0f+264.0f*a+c*(175.0f+606.0f*a))/2520.0f;
		s6 = (120.0f+c*(346.0f+127.0f*c))/5040.0f;
		ch += t*(1+0.5f*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
#else
		a = 0.5*t-b*c;
		s1 = (210.0+a*(140.0+a*(105.0+a*(84.0+a*(70.0+60.0*a))))) / 420.0;
		s2 = (420.0+a*(735.0+a*(966.0+a*(1141.0+1278.0*a))))/2520.0;
		s3 = (210.0+a*(462.0+a*(707.0+932.0*a)))/2520.0;
		s4 = (252.0+a*(672.0+1182.0*a)+c*(294.0+a*(889.0+1740.0*a)))/5040.0;
		s5 = (84.0+264.0*a+c*(175.0+606.0*a))/2520.0;
		s6 = (120.0+c*(346.0+127.0*c))/5040.0;
		ch += t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
#endif
		if (fabs(q/ch-ONE_POINT_ZERO) > e) 
			goto l4;
		return (ch);

}

//function taken from MB and hard wired for use here	
void Model::DiscreteGamma(FLOAT_TYPE *rates, FLOAT_TYPE *props, FLOAT_TYPE shape){
	bool median=false;
	int 	i;
	FLOAT_TYPE 	gap05 = ZERO_POINT_FIVE/NRateCats(), t, factor = shape/shape*NRateCats(), lnga1;

	if (median){
		for (i=0; i<NRateCats(); i++) 
			rates[i] = POINTGAMMA((i/ZERO_POINT_FIVE+ONE_POINT_ZERO)*gap05, shape, shape);
		for (i=0,t=0; i<NRateCats(); i++) 
			t += rates[i];
		for (i=0; i<NRateCats(); i++)     
			rates[i] *= factor/t;
		}
	else {
		lnga1 = LnGamma(shape+1);
		
		//DZ HACK
		//I don't think that these lines are needed, since the frequencies are fixed at .25 anyway.
		for (i=0; i<NRateCats()-1; i++) 
			props[i] = POINTGAMMA((i+ONE_POINT_ZERO)/NRateCats(), shape, shape);
		for (i=0; i<NRateCats()-1; i++) 
			props[i] = IncompleteGamma(props[i]*shape, shape+1, lnga1);
		
		
		rates[0] = props[0]*factor;
		rates[NRateCats()-1] = (1-props[NRateCats()-2])*factor;
		for (i=1; i<NRateCats()-1; i++)  
			rates[i] = (props[i]-props[i-1])*factor;
		}
	for (i=0; i<NRateCats(); i++) 
		props[i]=(ONE_POINT_ZERO-*propInvar)/NRateCats();
	}	
	
void Model::OutputPaupBlockForModel(ofstream &outf, const char *treefname) const{
	assert(modSpec.IsNucleotide());
	outf << "begin paup;\nclear;\ngett file=" << treefname << " storebr;\nlset userbr ";
	if(modSpec.Nst() == 2) outf << "nst=2 trat= " << TRatio();
	else if(modSpec.Nst() == 1) outf << "nst=1 ";
	else{
		if(modSpec.IsArbitraryRateMatrix()) outf << "nst=6 rclass=" << modSpec.arbitraryRateMatrixString.c_str() << " rmat=(" << Rates(0) << " " << Rates(1) << " " << Rates(2) << " " << Rates(3) << " " << Rates(4) << ")";
		else outf << "nst=6 rmat=(" << Rates(0) << " " << Rates(1) << " " << Rates(2) << " " << Rates(3) << " " << Rates(4) << ")";
		}
	
	if(modSpec.IsEqualStateFrequencies() == true) outf << " base=eq ";
	else if(modSpec.IsEmpiricalStateFrequencies() == true) outf << " base=emp ";
	else outf << " base=(" << StateFreq(0) << " " << StateFreq(1) << " " << StateFreq(2) << ")";
	
	if(modSpec.IsFlexRateHet() == false){
		if(NRateCats()>1) outf << " rates=gamma shape= " << Alpha() << " ncat=" << NRateCats();
		else outf << " rates=equal";
		outf << " pinv= " << PropInvar();
		outf << ";\nend;\n";
		}
	else{
		outf << " pinv= " << PropInvar();
		outf << " [FLEX RATES:\t";
		for(int i=0;i<NRateCats();i++){
			outf << rateMults[i] << "\t";
			outf << rateProbs[i] << "\t";
			}
		outf << "];\nend;\n";
		outf << "[!THIS TREE INFERRED UNDER FLEX RATE MODEL WITH GARLI.\nNO COMPARABLE MODEL IS AVAILABLE IN PAUP!]" << endl;
		}
	}

void Model::FillPaupBlockStringForModel(string &str, const char *treefname) const{
	char temp[200];
	sprintf(temp, "begin paup;\nclear;\ngett file=%s storebr;\nlset userbr ", treefname);
	str += temp;
	if(modSpec.Nst() == 2){
		sprintf(temp, "nst=2 trat=%f ", TRatio());
		str += temp;
		}
	else if(modSpec.Nst() == 1) str += "nst=1 ";
	else{
		if(modSpec.IsArbitraryRateMatrix())
			sprintf(temp,"nst=6 rclass=%s rmat=(%f %f %f %f %f)", modSpec.arbitraryRateMatrixString.c_str(), Rates(0), Rates(1), Rates(2), Rates(3), Rates(4));
		else
			sprintf(temp,"nst=6 rmat=(%f %f %f %f %f)", Rates(0), Rates(1), Rates(2), Rates(3), Rates(4));
		str += temp;
		}
	if(modSpec.IsEqualStateFrequencies()) str +=" base=eq ";
	else if(modSpec.IsEmpiricalStateFrequencies()) str += " base=emp ";
	else{
		sprintf(temp," base=( %f %f %f)", StateFreq(0), StateFreq(1), StateFreq(2));
		str += temp;
		}

	if(modSpec.IsFlexRateHet()==false){
		if(NRateCats()>1){
			sprintf(temp, " rates=gamma shape=%f ncat=%d", Alpha(), NRateCats());
			str += temp;
			}
		else str += " rates=equal";
		sprintf(temp, " pinv=%f;\nend;\n", PropInvar());;
		str += temp;
		}
	else{
		sprintf(temp, " pinv=%f  [FLEX RATES:\t", PropInvar());
		str += temp;
		for(int i=0;i<NRateCats();i++){
			sprintf(temp, "%f\t%f\t", rateMults[i], rateProbs[i]);
			str += temp;
			}
		str += "];\nend;\n[!THIS TREE INFERRED UNDER FLEX RATE MODEL WITH GARLI.\nNO COMPARABLE MODEL IS AVAILABLE IN PAUP!]\n";
		}
	}

void Model::OutputGarliFormattedModel(ostream &outf) const{
	//no reason to have different versions of the same thing, so just use the fill string function
	string s;
	this->FillGarliFormattedModelString(s);
	outf << s.c_str();
	return;
/*
	if(modSpec.IsCodon()){
		outf << "o ";
		for(int i=0;i<omegas.size();i++){
			outf << *omegas[i] << " " << *omegaProbs[i] << " ";
			}
		}

	if(modSpec.IsAminoAcid() == false)
		//outf << " r " << Rates(0) << " " << Rates(1) << " " << Rates(2) << " " << Rates(3) << " " << Rates(4);

	outf << " e " ;
	for(int i=0;i<nstates;i++)
		outf << StateFreq(i) << " ";;
	
	if(modSpec.IsFlexRateHet()){
		outf << " f ";
		for(int i=0;i<NRateCats();i++){
			outf << " " << rateMults[i] << "\t";
			outf << rateProbs[i] << "\t";
			}
		}
	else{
		if(NRateCats()>1 && modSpec.IsNonsynonymousRateHet() == false) outf << " a " << Alpha();
		}
	if(PropInvar()!=ZERO_POINT_ZERO) outf << " p " << PropInvar();
	outf << " ";
*/	}

void Model::FillModelOrHeaderStringForTable(string &s, bool model) const{	
	s.clear();
	char cStr[500];
	if(modSpec.IsCodon()){
		for(int i=0;i<omegas.size();i++){
			if(model){
				sprintf(cStr," %5.3f %5.3f", *omegas[i], *omegaProbs[i]);
				s += cStr;
				}
			else{
				char oStr[50];
				sprintf(oStr, "w(%d)", i);
				sprintf(cStr," %5s", oStr);
				s += cStr;
				sprintf(oStr, "p(%d)", i);
				sprintf(cStr," %5s", oStr);
				s += cStr;
				}
			}
		}
	if(modSpec.IsNucleotide() || modSpec.IsCodon() || modSpec.IsEstimateAAMatrix()){
		if(model){
			//sprintf(cStr, " %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f", Rates(0), Rates(1), Rates(2), Rates(3), Rates(4), 1.0);
			for(int st = 0;st < relNucRates.size();st++){
				sprintf(cStr," %6.3f", Rates(st));
				s += cStr;
				}
			}
		else{
			string states;
			if(modSpec.IsAminoAcid())
				states="ACDEFGHIKLMNPQRSTVWY";
			else
				states="ACGT";
			char rStr[50];
			for(int from=0;from<(modSpec.IsCodon() ? 4 - 1 : NStates() - 1);from++){
				for(int to=from+1;to<(modSpec.IsCodon() ? 4 : NStates());to++){
					sprintf(rStr, "r(%c%c)", states[from], states[to]);
					sprintf(cStr," %6s", rStr);
					s += cStr;
					}
				}
/*
			char rStr[50];
			sprintf(rStr, "r(AC)");
			sprintf(cStr," %5s", rStr);
			s += cStr;
			sprintf(rStr, "r(AG)");
			sprintf(cStr," %5s", rStr);
			s += cStr;
			sprintf(rStr, "r(AT)");
			sprintf(cStr," %5s", rStr);
			s += cStr;
			sprintf(rStr, "r(CG)");
			sprintf(cStr," %5s", rStr);
			s += cStr;
			sprintf(rStr, "r(CT)");
			sprintf(cStr," %5s", rStr);
			s += cStr;
			sprintf(rStr, "r(GT)");
			sprintf(cStr," %5s", rStr);
			s += cStr;
*/			}
		}
/*	if(modSpec.IsNucleotide()){
		if(model){
			sprintf(cStr," %5.3f %5.3f %5.3f %5.3f ", StateFreq(0), StateFreq(1), StateFreq(2), StateFreq(3));
			s += cStr;
			}
		else{
			char pStr[50];
			sprintf(pStr, "pi(A)");
			sprintf(cStr,"%5s ", pStr);
			s += cStr;
			sprintf(pStr, "pi(C)");
			sprintf(cStr,"%5s ", pStr);
			s += cStr;
			sprintf(pStr, "pi(G)");
			sprintf(cStr,"%5s ", pStr);
			s += cStr;
			sprintf(pStr, "pi(T)");
			sprintf(cStr,"%5s ", pStr);
			s += cStr;
			}
		}
*/	//else if(modSpec.IsAminoAcid()){
		if(modSpec.IsNucleotide() || (modSpec.IsAminoAcid() && (modSpec.fixStateFreqs == false && modSpec.IsEqualStateFrequencies() == false && modSpec.IsEmpiricalStateFrequencies() == false))){
			if(model){
				for(int st = 0;st < stateFreqs.size();st++){
					sprintf(cStr," %5.3f", StateFreq(st));
					s += cStr;
					}
				}
			else{
				char pStr[50];
				string states;
				if(modSpec.IsAminoAcid())
					states="ACDEFGHIKLMNPQRSTVWY";
				else
					states="ACGT";
				for(int st = 0;st < stateFreqs.size();st++){
					sprintf(pStr,"pi(%c)", states[st]);
					sprintf(cStr," %5s", pStr);
					s += cStr;
					}
				}
			}

	if(modSpec.IsFlexRateHet()){
		for(int i=0;i<NRateCats();i++){
			if(model){
				sprintf(cStr, " %5.3f %5.3f", rateMults[i], rateProbs[i]);
				s += cStr;
				}
			else{
				char fStr[50];
				sprintf(fStr, "fr(%d)", i);
				sprintf(cStr," %5s", fStr);
				s += cStr;
				sprintf(fStr, "p(%d)", i);
				sprintf(cStr," %5s", fStr);
				s += cStr;
				}
			}
		}
	else{
		if(modSpec.IsGammaRateHet()){
			if(model)
				sprintf(cStr, " %5.3f", Alpha());
			else{
				sprintf(cStr, " %5s", "alpha");
				}
			s += cStr;
			}
		}
	if(PropInvar()!=ZERO_POINT_ZERO){
		if(model)
			sprintf(cStr, " %5.3f", PropInvar());
		else{
			sprintf(cStr, " %5s", "pinv");
			}
		s += cStr;
		}
	}

void Model::OutputAminoAcidRMatrixArray(ostream &out){
	//assert(el.size() == 400);
	//first make a full 20x20 matrix
	assert(modSpec.IsAminoAcid());
	vector<FLOAT_TYPE> el(400, ZERO_POINT_ZERO);
	vector<FLOAT_TYPE *>::iterator r = relNucRates.begin();
	FLOAT_TYPE tot = ZERO_POINT_ZERO;
	for(int from=0;from<20;from++){
		for(int to=from;to<20;to++){
			if(from == to)
				el[from * 20 + to] = 0.0;
			else{
				el[from * 20 + to] = **r;
				el[to * 20 + from] = **r;
				tot += **r;
				r++;
				}
			}
		}
	assert(r == relNucRates.end());
	char str[100];

	out << "Estimated AA rate matrices:" << endl;;
	out << "NOTE THAT THIS FUNCTION IS FAIRLY EXPERIMENTAL, SO CHECK YOUR OUTPUT AND LET ME KNOW OF ANY PROBLEMS\n" << endl;;
	out << "GARLI's order of AA's is alphabetically BY SINGLE LETTER CODE, i.e.:\n ACDEFGHIKLMNPQRSTVWY" << endl;
	out << "The correspondence with the 3-letter codes and full names is this:" << endl;

	out << "A\tAla\tAlanine\nC\tCys\tCysteine\nD\tAsp\tAspartic Acid\nE\tGlu\tGlutamic Acid\nF\tPhe\tPhenylalanine\nG\tGly\tGlycine\nH\tHis\tHistidine\n";
	out << "I\tIle\tIsoleucine\nK\tLys\tLysine\nL\tLeu\tLeucine\nM\tMet\tMethionine\nN\tAsn\tAsparagine\nP\tPro\tProline\nQ\tGln\tGlutamine\nR\tArg\tArginine\n";
	out << "S\tSer\tSerine\nT\tThr\tThreonine\nV\tVal\tValine\nW\tTrp\tTryptophan\nY\tTyr\tTyrosine\n" << endl;

	out << "Unfortunately, I beleive that GARLI, PAML, and MrBayes all have different orderings of the amino acids.  PAML" << endl;
	out << "is alphabetical by three-letter code, MrBayes is alphabetical by full name (same as PAML, but swap Gln and Glu), GARLI" << endl;
	out << "is alphabetical by single letter code.  Additionally, I believe that PAML takes the below diagonal matrix as input,"<< endl;
	out << "while GARLI and MrBayes take the upper." << endl;
	out << "I COULD BE WRONG ABOUT THIS, AND YOU SHOULD VERIFY THAT THE ABOVE FACTS ARE TRUE BEFORE USING THE BELOW MATRICES" << endl;
	out << "IN ANOTHER PROGRAM" << endl;
	
	out << "Following are the matrix inferred by GARLI in GARLI's order, then the same matrices ordered by the other systems" << endl;
	out << "described above.  Both the above and below diagonal versions appear for each." << endl;

	out << "The entries are scaled such that the mean rate is 100.  It can be rescaled by any constant factor without" << endl;
	out << "changing its meaning. Entries on the diagonal are all zero." << endl;
	out << "\nThe ABOVE diagonal SINGLE LETTER order is what would be fed back into GARLI as a starting condition to use this matrix" << endl;
	out << "in future analyses.  Here is a GARLI block that could be used to do this.  The values could be fixed for further analyses" << endl;
	out << "by setting \"ratematrix = fixed\" in the configuration file, or it could be used as starting values for another run estimating" << endl;
	out << "the full matrix by leaving \"ratematrix = estimate\".  The block itself could be put in the same file as a NEXUS" << endl;
	out << "data matrix, or put in a file (which must start with #NEXUS) specified on the streefname line of the configuarion file.\n" << endl;

	out << "begin garli;" << endl;
	out << "[this specifies an amino acid rate matrix, with AA's ordered alphabetically by SINGLE LETTER CODE]" << endl;
	out << "[it is the above diagonal portion of the matrix, in order across each row]" << endl;
	out << "r ";
	
	for(int from=0;from<19;from++){
		for(int to=from+1;to<20;to++){
			sprintf(str, "%.3f", (el[from * 20 + to] * (19000.0/tot)));
			out << str << " ";
			}
		}
	out << ";\nend;\n" << endl;

	int corThree[20] = {0, 14, 11, 2, 1, 13, 3, 5, 6, 7, 9, 8, 10, 4, 12, 15, 16, 18, 19, 17};
	int corFull[20] =  {0, 14, 11, 2, 1, 3, 13, 5, 6, 7, 9, 8, 10, 4, 12, 15, 16, 18, 19, 17};

	out << "This is the SINGLE LETTER order (GARLI), above diagonal matrix\n" << endl;
	out << "(this is what appears in the above GARLI block)\n";
	for(int from=0;from<19;from++){
		for(int to=0;to<20;to++){
			if(to <= from)
				out << "\t";
			else{
				sprintf(str, "%.3f", (el[from * 20 + to] * (19000.0/tot)));
				out << str;
				if(to != from -1)
					out << "\t";
				}
			}
		out << endl;
		}
	out << endl;

	out << "\nThis is the SINGLE LETTER order (GARLI) below diagonal matrix\n" << endl;

	for(int from=1;from<20;from++){
		for(int to=0;to<from;to++){
			sprintf(str, "%.3f", (el[from * 20 + to] * (19000.0/tot)));
			out << str;
			if(to != from -1)
				out << "\t";
			}
		out << endl;
		}
	out << endl;

	out << "This is the THREE LETTER code order (PAML) above diagonal matrix\n" << endl;
	//above diagonal
	for(int from=0;from<19;from++){
		for(int to=0;to<20;to++){
			if(to <= from)
				out << "\t";
			else{
				sprintf(str, "%.3f", (el[corThree[from] * 20 + corThree[to]] * (19000.0/tot)));
				out << str;
				if(to != from -1)
					out << "\t";
				}
			}
		out << endl;
		}
	out << endl;

	out << "\nThis is the THREE LETTER order (PAML), below diagonal matrix\n" << endl;

	for(int from=1;from<20;from++){
		for(int to=0;to<from;to++){
			sprintf(str, "%.3f", (el[corThree[from] * 20 + corThree[to]] * (19000.0/tot)));
			out << str;
			if(to != from -1)
				out << "\t";
			}
		out << endl;
		}
	out << endl;

	out << "This is the FULL NAME code order (MrBayes) above diagonal matrix\n" << endl;
	//above diagonal
	for(int from=0;from<19;from++){
		for(int to=0;to<20;to++){
			if(to <= from)
				out << "\t";
			else{
				sprintf(str, "%.3f", (el[corFull[from] * 20 + corFull[to]] * (19000.0/tot)));
				out << str;
				if(to != from -1)
					out << "\t";
				}
			}
		out << endl;
		}
	out << endl;

	out << "\nThis is the FULL NAME order (MrBayes), below diagonal matrix\n" << endl;

	for(int from=1;from<20;from++){
		for(int to=0;to<from;to++){
			sprintf(str, "%.3f", (el[corFull[from] * 20 + corFull[to]] * (19000.0/tot)));
			out << str;
			if(to != from -1)
				out << "\t";
			}
		out << endl;
		}
	out << endl;

	
	out << "These are AA frequencies that were used, which may have been estimated or not." << endl;
	string states="ACDEFGHIKLMNPQRSTVWY";
	out << "Single letter order" << endl;
	for(int st = 0;st < stateFreqs.size();st++){
		out << states[st] << "\t";
		}
	out << endl;
	for(int st = 0;st < stateFreqs.size();st++){
		out << StateFreq(st) << "\t";
		}
	out << endl;
	out << "Three letter order" << endl;
	for(int st = 0;st < stateFreqs.size();st++){
		out << StateFreq(corThree[st]) << "\t";
		}
	out << endl;
	out << "Full name order" << endl;
	for(int st = 0;st < stateFreqs.size();st++){
		out << StateFreq(corFull[st]) << "\t";
		}
	out << endl;
	}

void Model::OutputHumanReadableModelReportWithParams() const{
	//Report on the model setup and parameter values - like a beefed up version of Population::ModelReport
	if(modSpec.IsCodon()){
		if(modSpec.IsVertMitoCode()) outman.UserMessage("  Number of states = 60 (codon data, vertebrate mitochondrial code)");
		else if(modSpec.IsInvertMitoCode()) outman.UserMessage("  Number of states = 62 (codon data, invertebrate mitochondrial code)");
		else outman.UserMessage("  Number of states = 61 (codon data, standard code)");
		}
	else if(modSpec.IsAminoAcid())
		outman.UserMessage("  Number of states = 20 (amino acid data)");
	else 
		outman.UserMessage("  Number of states = 4 (nucleotide data)");
	
	if(modSpec.IsAminoAcid() == false){
		if(modSpec.IsCodon() && modSpec.numRateCats == 1){
			if(!modSpec.fixOmega)
				outman.UserMessageNoCR("  One estimated dN/dS ratio (aka omega) = %f\n", Omega(0));
			else
				outman.UserMessageNoCR("  One estimated dN/dS ratio (aka omega).\n    Value provided by user (fixed) = %f\n", Omega(0));
			}
		if(modSpec.IsCodon()) outman.UserMessage("  Nucleotide Relative Rate Matrix Assumed by Codon Model:     ");
		else outman.UserMessage("  Nucleotide Relative Rate Matrix: ");
		if(modSpec.Nst() == 6){
			if(modSpec.IsArbitraryRateMatrix()) outman.UserMessageNoCR("    User specified matrix type: %s ", modSpec.arbitraryRateMatrixString.c_str());
			else outman.UserMessageNoCR("    6 rates ");
			if(modSpec.fixRelativeRates == true) outman.UserMessage(" values specified by user (fixed)");
			else outman.UserMessage("");
			outman.UserMessage("    AC = %.3f, AG = %.3f, AT = %.3f, CG = %.3f, CT = %.3f, GT = %.3f", Rates(0), Rates(1), Rates(2), Rates(3), Rates(4), 1.0);
			}
		else if(modSpec.Nst() == 2){
			outman.UserMessageNoCR("    2 rates (transition and transversion) K param = %.4f", Rates(1));
			if(modSpec.IsCodon() == false) outman.UserMessage(" (ti/tv = %.4f)",  TRatio());
			else outman.UserMessage("");
			}
		else outman.UserMessage("    1 rate");
		}
	else{
		outman.UserMessageNoCR("  Amino Acid Rate Matrix: ");
		if(modSpec.IsJonesAAMatrix()) outman.UserMessage("Jones");
		else if(modSpec.IsDayhoffAAMatrix()) outman.UserMessage("Dayhoff");
		else if(modSpec.IsPoissonAAMatrix()) outman.UserMessage("Poisson");
		else if(modSpec.IsWAGAAMatrix()) outman.UserMessage("WAG");
		else if(modSpec.IsMtMamAAMatrix()) outman.UserMessage("MtMam");
		else if(modSpec.IsMtRevAAMatrix()) outman.UserMessage("MtRev");
		else if(modSpec.IsEstimateAAMatrix()) outman.UserMessage("Estimated (189 free parameters)");
		else if(modSpec.IsUserSpecifiedRateMatrix()) outman.UserMessage(" values specified by user (fixed)");
		}

	outman.UserMessageNoCR("  Equilibrium State Frequencies: ");
	if(modSpec.IsEqualStateFrequencies()){
		if(modSpec.IsCodon()){
			if(modSpec.IsVertMitoCode()) outman.UserMessage("equal (1/60 = 0.01667, fixed)");
			else if(modSpec.IsInvertMitoCode()) outman.UserMessage("equal (1/62 = 0.01613, fixed)");
			else outman.UserMessage("equal (1/61 = 0.01639, fixed)");
			}
		else if(modSpec.IsAminoAcid())
			outman.UserMessage("equal (0.05, fixed)");
		else 
			outman.UserMessage("equal (0.25, fixed)");
		}
	else if(modSpec.IsF3x4StateFrequencies()) outman.UserMessage("\n    empirical values calculated by F3x4 method (fixed)");
	else if(modSpec.IsF1x4StateFrequencies()) outman.UserMessage("\n    empirical values calculated by F1x4 method (fixed)");
	else if(modSpec.IsEmpiricalStateFrequencies()){
		if(modSpec.IsAminoAcid()) outman.UserMessage("empirical (observed) values (+F)");
		else outman.UserMessage("empirical (observed) values, fixed:");
		}
	else if(modSpec.IsJonesAAFreqs()) outman.UserMessage("Jones");
	else if(modSpec.IsWAGAAFreqs()) outman.UserMessage("WAG");
	else if(modSpec.IsMtMamAAFreqs()) outman.UserMessage("MtMam");
	else if(modSpec.IsMtRevAAFreqs()) outman.UserMessage("MtRev");
	else if(modSpec.IsDayhoffAAFreqs()) outman.UserMessage("Dayhoff");
	else if(modSpec.IsUserSpecifiedStateFrequencies()) outman.UserMessage("specified by user (fixed)");
	else outman.UserMessage("estimated");
		
	if(!modSpec.IsEqualStateFrequencies()){
		if(modSpec.IsCodon())  outman.UserMessageNoCR("    (AAA, AAC, AAG, AAT, ACA, ... etc)\n    ");
		else if(modSpec.IsAminoAcid()) outman.UserMessageNoCR("    (ACDEFGHIKLMNPQRSTVWY)\n    ");
		else outman.UserMessageNoCR("    (ACGT) ");
		for(int i=0;i<nstates;i++){
			outman.UserMessageNoCR("%.4f ", StateFreq(i));
			if(i>0 && (i+1)!= nstates && !((i+1)%5)) outman.UserMessageNoCR("\n    ");
			}
		outman.UserMessage("");
		}

	outman.UserMessage("  Rate Heterogeneity Model:");
	if(modSpec.numRateCats == 1){
		if(modSpec.includeInvariantSites == false) outman.UserMessage("    no rate heterogeneity");
		else{
			if(modSpec.fixInvariantSites == true) outman.UserMessage("    only an invariant (invariable) site category,\n    proportion specified by user (fixed)\n    %.4f", PropInvar());
			else outman.UserMessage("    only an invariant (invariable) site category, proportion estimated\n    %.4f", PropInvar());
			}
		}
	else{
		outman.UserMessageNoCR("    %d ", modSpec.numRateCats);
		if(modSpec.IsNonsynonymousRateHet()){
			if(!modSpec.fixOmega){
				outman.UserMessage("nonsynonymous rate categories, rate and proportion of each estimated\n     (this is effectively the M3 model of PAML)");
				}
			else{
				outman.UserMessage("nonsynonymous rate categories, rate and proportion of each provided by user (fixed)\n     (this is effectively the M3 model of PAML)");
				}
			outman.UserMessage("      dN/dS\tProportion");
			for(int i=0;i<modSpec.numRateCats;i++)
				outman.UserMessage("      %5.4f\t%5.4f", Omega(i), OmegaProb(i));
			}
		else if(modSpec.IsFlexRateHet() == false){
			if(modSpec.fixAlpha == true) outman.UserMessage("discrete gamma distributed rate categories,\n      alpha param specified by user (fixed)\n      %.4f", Alpha());
			else outman.UserMessage("discrete gamma distributed rate categories, alpha param estimated\n      %.4f", Alpha());
			if(modSpec.includeInvariantSites == true){
				if(modSpec.fixInvariantSites == true) outman.UserMessage("    with an invariant (invariable) site category,\n    proportion specified by user (fixed)\n      %.4f", PropInvar());				
				else outman.UserMessage("    with an invariant (invariable) site category, proportion estimated\n      %.4f", PropInvar());	
				}
			outman.UserMessage("    Substitution rate categories under this model:\n      rate\tproportion");
			if(modSpec.includeInvariantSites == true) outman.UserMessage("      %5.4f\t%5.4f", 0.0, PropInvar());
			for(int r=0;r<modSpec.numRateCats;r++)
				outman.UserMessage("      %5.4f\t%5.4f", rateMults[r], rateProbs[r]);
			}
		else{
			outman.UserMessage("FLEX rate categories, rate and proportion of each estimated");
			if(modSpec.includeInvariantSites == true){
				if(modSpec.fixInvariantSites == true) outman.UserMessage("    with an invariant (invariable) site category,\n    proportion specified by user (fixed)");				
				else outman.UserMessage("    with an invariant (invariable) site category, proportion estimated");
				}
			outman.UserMessage("      Estimated substitution rate categories:\n      rate\tproportion");
			for(int r=0;r<modSpec.numRateCats;r++)
				outman.UserMessage("      %5.4f\t%5.4f", rateMults[r], rateProbs[r]);
			}
		}
	outman.UserMessage("");

	}

#define MODEL_OUTPUT_PREC 5
void Model::FillGarliFormattedModelString(string &s) const{
	char temp[1000];
	int prec = MODEL_OUTPUT_PREC;
	if(modSpec.IsCodon()){
		s += " o";
		for(int i=0;i<omegas.size();i++){
			sprintf(temp," %.*f %.*f",  prec, *omegas[i],  prec, *omegaProbs[i]);
			s += temp;
			}
		}
	if(modSpec.IsAminoAcid() == false || modSpec.IsEstimateAAMatrix() || (modSpec.IsAminoAcid() && modSpec.IsUserSpecifiedRateMatrix())){
		//sprintf(temp," r %.*f %.*f %.*f %.*f %.*f", prec, Rates(0), prec, Rates(1), prec, Rates(2), prec, Rates(3), prec, Rates(4));
		//s += temp;
		s += " r ";
		for(int st = 0;st < relNucRates.size();st++){
			sprintf(temp," %.*f", prec, Rates(st));
			s += temp;
			}
		}
	if(modSpec.IsNucleotide()){
		sprintf(temp," e %.*f %.*f %.*f %.*f",  prec, StateFreq(0),  prec, StateFreq(1),  prec, StateFreq(2),  prec, StateFreq(3));
		s += temp;
		}
	else{
		sprintf(temp," e ");
		s += temp;
		for(int i=0;i<nstates;i++){
			sprintf(temp," %.*f ",  prec, StateFreq(i));
			s += temp;
			}
		}

	if(modSpec.IsFlexRateHet()){
		s += " f ";
		for(int i=0;i<NRateCats();i++){
			sprintf(temp, " %.*f %.*f ",  prec, rateMults[i],  prec, rateProbs[i]);
			s += temp;
			}
		}
	else{
		if(modSpec.IsGammaRateHet()){
			sprintf(temp, " a %.*f",  prec, Alpha());
			s += temp;
			}
		}
	if(PropInvar()!=ZERO_POINT_ZERO){
		sprintf(temp, " p %.*f",  prec, PropInvar());
		s += temp;
		}
	s += " ";
	}

/*	
void Model::ReadModelFromFile(NexusToken &token){
	token.GetNextToken();

	do{
		if(token.Equals("r")){//rate parameters
			token.GetNextToken();
			FLOAT_TYPE r[5];
			for(int i=0;i<5;i++){
				r[i]=atof(token.GetToken().c_str());
				token.GetNextToken();
				}
			SetRmat(r);
			if(token.IsNumericalToken()) token.GetNextToken();//this is necessary incase GT is included
			}
		else if(token.Equals("b")){
			token.GetNextToken();
			FLOAT_TYPE b[3];
			for(int i=0;i<3;i++){
				b[i]=atof(token.GetToken().c_str());
				token.GetNextToken();
				}
			SetPis(b);
			if(token.IsNumericalToken()) token.GetNextToken();//this is necessary incase T is included						
			}
		else if(token.Equals("a")){
			token.GetNextToken();
			SetAlpha(atof(token.GetToken().c_str()));
			token.GetNextToken();
			}				
		else if(token.Equals("p")){
			token.GetNextToken();
			SetPinv(atof(token.GetToken().c_str()));
			token.GetNextToken();
			}
		else if(token.Begins("(") == false){
			token.GetNextToken();
			}
		}while(token.Begins("(") == false);
	UpdateQMat();
}
*/	
	
void Model::ReadGarliFormattedModelString(string &modString){
	istringstream stf(modString, stringstream::in);

	char c;
	NxsString temp;
	c=stf.get();
	do{//read parameter values identified by single letter identifier.  Each section should
		//take care of advancing to the following letter 
		if(c == 'R' || c == 'r'){//rate parameters
			if(modSpec.IsAminoAcid() && modSpec.IsEstimateAAMatrix() == false && modSpec.IsUserSpecifiedRateMatrix() == false) 
				throw ErrorException("Amino acid rate matrix parameters cannot be specified unless \"ratematrix = fixed\" or \"ratematrix = estimate\" are used.");
			//FLOAT_TYPE r[6];
			vector<FLOAT_TYPE> r;
			//for(int i=0;i<5;i++){
			for(int i=0;i<relNucRates.size() - 1;i++){
				temp.clear();
				stf >> temp;

				if(temp.size() == 0)
					throw(ErrorException("Unexpected end of model string while reading rate matrix parameters.\nExamine file and check manual for format.\n"));
				if(temp[0] != '.' && (!isdigit(temp[0]))) 
					throw(ErrorException("Problem reading rate matrix parameters from file (maybe too few rates?).\n\tFor amino acid models 190 rates should be specified, (or 189 rates if the last rate is assumed to be 1.0).\n\tFor nucleotide models 6 should be specified (or 5 if the last rate is assumed to be 1.0).\n\tExamine file and check manual or website for format.\n"));
				//r[i]=(FLOAT_TYPE)atof(temp.c_str());
				r.push_back((FLOAT_TYPE)atof(temp.c_str()));
				}
			do{
				c=stf.get();
				}while(!stf.eof() && (c == ' ' || c == '\t'));	
			if(isdigit(c) || c=='.'){//read the GT rate, if specified
				string v;
				v = c;
				while((!isalpha(c) || c=='e' || c == 'E') && !stf.eof()){
					c=stf.get();
					if(isdigit(c) || c=='.' || c=='e' || c=='E' || c == '-')
						v += c;
					else if(c == ' ' || c == '\t'){
						c=stf.get();
						if(isdigit(c) || c=='.')
							throw ErrorException("It appears that too many relative rates was specified in the model string.\n\tFor amino acid models 190 rates should be specified, (or 189 rates if the last rate is assumed to be 1.0).\n\tFor nucleotide models 6 should be specified (or 5 if the last rate is assumed to be 1.0).");
						break;
						}
					}
				//r[5] = atof(v.c_str());
				r.push_back((FLOAT_TYPE)atof(v.c_str()));
				}
			//else r[5] = ONE_POINT_ZERO;
			else 
				r.push_back(ONE_POINT_ZERO);
			if(r.size() != relNucRates.size()){
				if(modSpec.IsAminoAcid())
					throw ErrorException("It appears that too few relative rates were specified in the model string (found %d).\n\tFor amino acid models 190 rates should be specified, (or 189 rates if the last rate is assumed to be 1.0).", r.size());
				else
					throw ErrorException("Incorrect number of relative rates specified in the model string.\t6 rates should be specified, (or 5 rates if the G-T rate is assumed to be 1.0).");
				}
			SetRmat(&r[0], true, true);
			modSpec.gotRmatFromFile=true;
			}
		else if(c == 'E' || c == 'e' || c == 'b' || c == 'B'){//base freqs
			//7/12/07 changing this to pay attention to the 4th state, if specified
			//although it should be calcuable from the other three, having exact restartability
			//sometimes requires that it is taken as is
			//FLOAT_TYPE b[4];
			int nstates = modSpec.nstates;
			vector<FLOAT_TYPE> b(nstates);
			for(int i=0;i<nstates-1;i++){
				temp.clear();
				stf >> temp;
				if(temp.size() == 0)
					throw(ErrorException("Unexpected end of model string while reading equilibrium frequency parameters.\nExamine file and check manual for format.\n"));
				if(temp[0] != '.' && (!isdigit(temp[0]))) 
					throw(ErrorException("Problem reading equilibrium state frequency parameters from file.\nExamine file and check manual for format.\n"));
				b[i]=(FLOAT_TYPE)atof(temp.c_str());
				}
			do{
				c=stf.get();
				}while(!stf.eof() && (c == ' ' || c == '\t'));	
			if(isdigit(c) || c=='.'){
				string v;
				v = c;
				while(!isalpha(c) && !stf.eof()){
					c=stf.get();
					if(isdigit(c) || c=='.')
						v += c;
					else if(c == ' ' || c == '\t'){
						c=stf.get();
						if(isdigit(c) || c=='.')
							throw ErrorException("It appears that too many equilibrium frequencies were specified in the model string.\n\tFor amino acid models 20 should be specified, (or 19 if the last is assumed to make the sum 1.0).\n\tFor nucleotide models 4 (or 3 if the last is assumed to make the sum 1.0).");
						break;
						}
					}
				b[nstates-1]=(FLOAT_TYPE)atof(v.c_str());
				}
			else{
				FLOAT_TYPE tot = ZERO_POINT_ZERO;
				for(int i=0;i<nstates-1;i++) 
					tot += b[i];
				b[nstates-1] = ONE_POINT_ZERO - tot;
				}
			//in this case we're "forcing" estimation of state frequencies but providing starting values, 
			//and because this is rather a hack we can't actually do the validation without crapping out 
			if(modSpec.IsCodon() && modSpec.fixStateFreqs == false && modSpec.IsEmpiricalStateFrequencies())
				SetPis(&b[0], false, true);
			else
				SetPis(&b[0], true, true);
			modSpec.gotStateFreqsFromFile=true;
			}
		else if(c == 'A' || c == 'a'){//alpha shape
			if(modSpec.IsFlexRateHet()) 
				throw(ErrorException("Config file specifies ratehetmodel = flex, but starting model contains alpha!\n"));
			if(modSpec.IsNonsynonymousRateHet()) 
				throw(ErrorException("Config file specifies ratehetmodel = nonsynonymous, but starting model contains alpha!\n"));
			temp.clear();
			stf >> temp;
			if(temp.size() == 0)
				throw(ErrorException("Unexpected end of model string while reading alpha parameter.\nExamine file and check manual for format.\n"));
			if(temp[0] != '.' && (!isdigit(temp[0]))) 
				throw(ErrorException("Problem reading alpha parameter from file.\nExamine file and check manual for format.\n"));
			SetAlpha((FLOAT_TYPE)atof(temp.c_str()), true);
			c=stf.get();
			modSpec.gotAlphaFromFile=true;
			}				
		else if(c == 'P' || c == 'p' || c == 'i' || c == 'I'){//proportion invariant
			temp.clear();
			stf >> temp;
			if(temp.size() == 0)
				throw(ErrorException("Unexpected end of model string while reading proportion of invariant sites parameter.\nExamine file and check manual for format.\n"));
			if(temp[0] != '.' && (!isdigit(temp[0]))) 
				throw(ErrorException("Problem reading proportion of invariant sites parameter from file.\nExamine file and check manual for format.\n"));
			FLOAT_TYPE p=(FLOAT_TYPE)atof(temp.c_str());
			SetPinv(p, true);
			c=stf.get();
			modSpec.gotPinvFromFile=true;
			}
		else if(c == 'F' || c == 'f'){//flex rates
			if(modSpec.IsFlexRateHet()==false) 
				throw(ErrorException("Flex rate parameters specified, but ratehetmodel is not flex!\n"));
			FLOAT_TYPE rates[20];
			FLOAT_TYPE probs[20];
			for(int i=0;i<NRateCats();i++){
				temp.clear();
				stf >> temp;
				if(temp.size() == 0)
					throw(ErrorException("Unexpected end of model string while reading flex rate parameters.\nExamine file and check manual for format.\n"));
				if(isalpha(temp[0])) 
					throw ErrorException("Problem with flex rates specification in starting condition file");
				rates[i]=(FLOAT_TYPE)atof(temp.c_str());
				temp.clear();
				stf >> temp;
				if(temp.size() == 0)
					throw(ErrorException("Unexpected end of model string while reading flex rate parameters.\nExamine file and check manual for format.\n"));
				if(isalpha(temp[0])) 
					throw ErrorException("Problem with flex rates specification in starting condition file");
				probs[i]=(FLOAT_TYPE)atof(temp.c_str());
				}		
			SetFlexRates(rates, probs);
			NormalizeRates();
			c=stf.get();
			modSpec.gotFlexFromFile=true;
			}
		else if(c == 'O' || c == 'o'){//omega parameters
			if(modSpec.IsCodon() == false) 
				throw ErrorException("Omega parameters specified for non-codon model?");
			FLOAT_TYPE rates[20];
			FLOAT_TYPE probs[20];
			if(NRateCats() == 1){//just a single omega value to get, maybe with a proportion of 1.0 following it
				temp.clear();
				stf >> temp;
				if(temp.size() == 0)
					throw(ErrorException("Unexpected end of model string while reading omega parameters.\nExamine file and check manual for format.\n"));
				if(isalpha(temp[0])) 
					throw ErrorException("Problem with omega parameter specification in starting condition file");
				rates[0]=(FLOAT_TYPE)atof(temp.c_str());
				do{
					c=stf.get();
					}while(!stf.eof() && (c == ' ' || c == '\t'));	
				if(isdigit(c) || c=='.'){
					string v;
					v = c;
					temp.clear();
					stf >> temp;
					if(temp.size() == 0)
						throw(ErrorException("Unexpected end of model string while reading omega parameters.\nExamine file and check manual for format.\n"));
					v += temp;
					if(FloatingPointEquals(atof(v.c_str()), ONE_POINT_ZERO, 1.0e-5) == false)
						throw ErrorException("Problem with omega parameter specification in starting condition file\n(wrong number of rate cats specified in config?)");
					do{
						c=stf.get();
						}while(!stf.eof() && (c == ' ' || c == '\t'));		
					if(isdigit(c) || c == '.') 
						throw ErrorException("Problem with omega parameter specification in starting condition file");
					}
				probs[0] = ONE_POINT_ZERO;
				SetOmegas(rates, probs);
				}
			else{
				for(int i=0;i<NRateCats();i++){
					temp.clear();
					stf >> temp;
					if(temp.size() == 0)
						throw(ErrorException("Unexpected end of model string while reading omega parameters.\nExamine file and check manual for format.\n"));
					if(isalpha(temp[0])) 
						throw ErrorException("Problem with omega parameter specification in starting condition file");
					rates[i]=(FLOAT_TYPE)atof(temp.c_str());
					temp.clear();
					stf >> temp;
					if(temp.size() == 0)
						throw(ErrorException("Unexpected end of model string while reading omega parameters.\nExamine file and check manual for format.\n"));
					if(isalpha(temp[0])) 
						throw ErrorException("Problem with omega parameter specification in starting condition file");
					probs[i]=(FLOAT_TYPE)atof(temp.c_str());
					}
				do{
					c=stf.get();
					}while(!stf.eof() && (c == ' ' || c == '\t'));			
				if(isdigit(c) || c == '.') throw ErrorException("Problem with omega parameter specification in starting condition file");
				SetOmegas(rates, probs);
				}
			modSpec.gotOmegasFromFile=true;
			}
		else if(c == 'n'){
			//the number of cats should now be set in the config file
			c=stf.get();
			assert(0);
			}
		else if(isalpha(c)) 
			throw(ErrorException("Unknown model parameter specification in file.\nExamine file and check manual for format.\n"));
		else if(c != '(') c=stf.get();
		}while(c != '(' && c != '\r' && c != '\n' && !stf.eof());
	}

void Model::CreateModelFromSpecification(int modnum){
	nstates = modSpec.nstates;
	if(modSpec.IsNucleotide() || modSpec.IsCodon())
		nst = modSpec.Nst();
	
	else nst = -1;
	
	//deal with rate het models
	propInvar = new FLOAT_TYPE;
	if(modSpec.includeInvariantSites){
		assert(modSpec.IsCodon() == false);
		*propInvar=(FLOAT_TYPE)0.2;
		if(modSpec.fixInvariantSites == false){
			ProportionInvariant *pi = new ProportionInvariant("proportion invariant", (FLOAT_TYPE **) &propInvar);
			pi->SetWeight(1);
			paramsToMutate.push_back(pi);
			}			
		}
	else *propInvar=ZERO_POINT_ZERO;

	if(NRateCats() > 1 && modSpec.IsNonsynonymousRateHet() == false){
		//assert(modSpec.IsNucleotide() || modSpec.IsAminoAcid());
		alpha = new FLOAT_TYPE;
		*alpha = ZERO_POINT_FIVE;
		
		if(modSpec.IsFlexRateHet() == false){
			DiscreteGamma(rateMults, rateProbs, *alpha);
			if(modSpec.fixAlpha == false){
				AlphaShape *a= new AlphaShape("alpha", &alpha);
				a->SetWeight(1);
				paramsToMutate.push_back(a);
				}
			}
		else{
			//start the flex rates out being equivalent to
			//a gamma with alpha=.5
			DiscreteGamma(rateMults, rateProbs, ZERO_POINT_FIVE);
			if(modSpec.includeInvariantSites == true) NormalizeRates();

			vector<FLOAT_TYPE*> dummy;
			dummy.reserve(NRateCats());
			
			for(int i=0;i<NRateCats();i++)
				dummy.push_back(&rateProbs[i]);
			RateProportions *rateP=new RateProportions(&dummy[0], NRateCats());
			rateP->SetWeight((FLOAT_TYPE)NRateCats());
			paramsToMutate.push_back(rateP);

			dummy.clear();
			for(int i=0;i<NRateCats();i++)
				dummy.push_back(&rateMults[i]);
			RateMultipliers *rateM=new RateMultipliers(&dummy[0], NRateCats());
			rateM->SetWeight((FLOAT_TYPE)NRateCats());
			paramsToMutate.push_back(rateM);
			}
		}
	else{
		rateMults[0]=ONE_POINT_ZERO;
		rateProbs[0]=ONE_POINT_ZERO;
		alpha=NULL;
		}

	//deal with the state frequencies
	for(int i=0;i<nstates;i++){
		FLOAT_TYPE *f=new FLOAT_TYPE;
		*f=(ONE_POINT_ZERO/(FLOAT_TYPE) nstates);
		stateFreqs.push_back(f);
		}
	if(modSpec.IsEqualStateFrequencies() == false && modSpec.fixStateFreqs == false){
		StateFrequencies *s=new StateFrequencies(&stateFreqs[0], nstates);
		s->SetWeight(nstates);
		paramsToMutate.push_back(s);
		}
	if(modSpec.IsAminoAcid()){
		if(modSpec.IsJonesAAFreqs()) SetJonesAAFreqs();
		if(modSpec.IsDayhoffAAFreqs()) SetDayhoffAAFreqs();
		if(modSpec.IsWAGAAFreqs()) SetWAGAAFreqs();
		if(modSpec.IsMtMamAAFreqs()) SetMtMamAAFreqs();
		if(modSpec.IsMtRevAAFreqs()) SetMtRevAAFreqs();
		}

	//deal with the relative rates

	if(modSpec.IsAminoAcid() == false){
		if(nst==6){
			if(modSpec.IsArbitraryRateMatrix()){
				//user specified rate matrix type, like rclass = (a b c d e f) in paup
				//trying to do this as generically as possible
				string matrixSpec = modSpec.GetArbitraryRateMatrixString();
				int pos = 0;
				char characters[10];
		//		int usedCharacters = 0;
				FLOAT_TYPE **params = new FLOAT_TYPE*[6]; 

				for(int r=0;r<6;r++){
					while (pos < matrixSpec.size() && !isalnum(matrixSpec[pos])) pos++;
					bool newChar = true;
					char thisChar = matrixSpec[pos];
					for(int c=0;c<r;c++){
						if(thisChar == characters[c]){
							params[r] = params[c];
							newChar = false;
							arbitraryMatrixIndeces[r] = c;
							characters[r] = thisChar;
							break;
							}
						}
					if(newChar){
						params[r] = new FLOAT_TYPE;
						*params[r] = ONE_POINT_ZERO;
						characters[r] = thisChar;
						arbitraryMatrixIndeces[r] = r;
	//					usedCharacters++;
						}
					pos++;
					if(matrixSpec[pos] != ' ' && matrixSpec[pos] != '\t' && matrixSpec[pos] != ')') throw ErrorException("Problem parsing rate matrix specification.\n\tIt should look something like this (a b a c b c)");
					}
				while (pos < matrixSpec.size() && (matrixSpec[pos] == ' ' || matrixSpec[pos] == '\t')) pos++;
				if(matrixSpec[pos] != ')') throw ErrorException("Problem parsing rate matrix specification.\n\tIt should look something like this (a b a c b c)");
				for(int i=0;i<6;i++){
					relNucRates.push_back(params[i]);
					}

				delete []params;
				//it is hard to know what values to start the arbitrary matrix at, but try to make the transitions higher
				//it doesn't really matter if some are aliased to one another
				*relNucRates[1] = 4.0;
				*relNucRates[4] = 4.0;
				*relNucRates[5] = ONE_POINT_ZERO;
				RelativeRates *r=new RelativeRates("Rate matrix", &relNucRates[0], 6, 1e-5, 999.9);
				r->SetWeight(6);
				paramsToMutate.push_back(r);
				}
			else{//normal GTR
				//make the transitions higher to begin with
				for(int i=0;i<6;i++){
					FLOAT_TYPE *d=new FLOAT_TYPE;
					relNucRates.push_back(d);
					}
				*relNucRates[0]=*relNucRates[2]=*relNucRates[3]=*relNucRates[5] = ONE_POINT_ZERO;
				*relNucRates[1]=*relNucRates[4] = 4.0;
				}
			if(modSpec.fixRelativeRates == false){
				RelativeRates *r=new RelativeRates("Rate matrix", &relNucRates[0], 6, 1e-3, 999.9);
	
				r->SetWeight(6);
				paramsToMutate.push_back(r);
				}
			}
		else if(nst==2){
			FLOAT_TYPE *a=new FLOAT_TYPE;
			FLOAT_TYPE *b=new FLOAT_TYPE;
			*a=ONE_POINT_ZERO;
			*b=4.0;
			relNucRates.push_back(a);
			relNucRates.push_back(b);
			relNucRates.push_back(a);
			relNucRates.push_back(a);
			relNucRates.push_back(b);
			relNucRates.push_back(a);
			
			if(modSpec.fixRelativeRates == false){
				RelativeRates *r=new RelativeRates("Rate matrix", &b, 1, 1e-3, 999.9);
				r->SetWeight(2);
				paramsToMutate.push_back(r);
				}
			}
		else if(nst==1){
			FLOAT_TYPE *a=new FLOAT_TYPE;
			*a=ONE_POINT_ZERO;
			for(int i=0;i<6;i++)
				relNucRates.push_back(a);
			}
		}
	else{//estimating or fixing the aminoacid rate matrix
		if(modSpec.fixRelativeRates == false || modSpec.IsUserSpecifiedRateMatrix()){
			int seed = rnd.seed();
			for(int i=0;i<190;i++){
				FLOAT_TYPE *d=new FLOAT_TYPE;
				//*d = ONE_POINT_ZERO;
				if(i == 189)
					*d = 1.0;
				else
					*d = max(rnd.gamma(1), MIN_REL_RATE);
				relNucRates.push_back(d);
				}
			rnd.set_seed(seed);
#ifdef SUM_AA_REL_RATES	
			this->NormalizeSumConstrainedRelativeRates(true, -1);
#endif
			if(! modSpec.IsUserSpecifiedRateMatrix()){
#ifdef SUM_AA_REL_RATES
				SumConstrainedRelativeRates *r = new SumConstrainedRelativeRates("Rate matrix", &relNucRates[0], 190, SUM_TO * 1.0e-6/190.0, SUM_TO * 1.0e6/190.0, SUM_TO);
#else
				RelativeRates *r=new RelativeRates("Rate matrix", &relNucRates[0], 190, 1e-3, 9999.9);
#endif
				
				r->SetWeight(190);
				paramsToMutate.push_back(r);
				}
			}
		}

	AllocateEigenVariables();//these need to be allocated regardless of
		//nst because I don't feel like simplifying the deriv calcs for simpler
		//models.  Pmat calcs for simpler models are simplified, and don't
		//require the Eigen stuff	

	if(modSpec.IsCodon() == false) 
		UpdateQMat();
	else{
		FLOAT_TYPE *d;
		for(int i=0;i<NRateCats();i++){
			d = new FLOAT_TYPE;
			*d = 0.25 * (FLOAT_TYPE) (i + 1);
			omegas.push_back(d);
			d = new FLOAT_TYPE;
			*d = 1.0 / NRateCats();
			omegaProbs.push_back(d);

			//changes for beagle, should have no side effects
			//The overall rate multipliers are all equal for different NS rate cats, even if omega varies
			rateMults[i] = 1.0;
			rateProbs[i] = *omegaProbs[i];
			}

/*		*omegas[0] = 0.0000;
		*omegas[1] = 1.24023;
		*omegas[2] = 2.99539;
*/
/*		if(NRateCats() == 1){
			*omegas[0] = 1.0000;
			*omegaProbs[0] = 1.0;
			}
		else{
			*omegas[0] = 0.00000;
			*omegas[1] = 0.79116;
			*omegas[2] = 1.96505;

			*omegaProbs[0] = 0.64547;
			*omegaProbs[1] = 0.21651;
			*omegaProbs[2] = 0.13802;
			}
*/
/*		*omegas[0] = 0.8;
		*omegas[1] = 1.0;
		*omegas[2] = 1.2;
*/

/*		*omegaProbs[0] = 0.68280;
		*omegaProbs[1] = 0.28284;
		*omegaProbs[2] = 0.03436;
*/

/*		rateProbs[0] = 0.68280;
		rateProbs[1] = 0.28284;
		rateProbs[2] = 0.03436;
*/
		//*relNucRates[1] = 2.89288;

		if(!modSpec.fixOmega){
			if(NRateCats() > 1){
				RateProportions *omegaP=new RateProportions(&omegaProbs[0], NRateCats());
				omegaP->SetWeight((FLOAT_TYPE)NRateCats());
				paramsToMutate.push_back(omegaP);
				}
				
			RateMultipliers *omegaM=new RateMultipliers(&omegas[0], NRateCats());
			omegaM->SetWeight((FLOAT_TYPE)NRateCats());
			paramsToMutate.push_back(omegaM);
			}

/*		FLOAT_TYPE *NS=new FLOAT_TYPE;
		*NS = 0.5;
		FLOAT_TYPE *S=new FLOAT_TYPE;
		*S = 1.0;
		omegas.push_back(NS);
		omegas.push_back(S);
		RelativeRates *o=new RelativeRates("Omega", &omega[0], 2);

		o->SetWeight(2);
		paramsToMutate.push_back(o);
*/
		UpdateQMatCodon();
		}

	eigenDirty=true;
	}

void Model::SetMtMamAAFreqs(){
	*stateFreqs[0] 	=	0.0692	;
	*stateFreqs[14]	=	0.0184	;
	*stateFreqs[11]	=	0.0400	;
	*stateFreqs[2]	=	0.0186	;
	*stateFreqs[1]	=	0.0065	;
	*stateFreqs[13]	=	0.0238	;
	*stateFreqs[3]	=	0.0236	;
	*stateFreqs[5]	=	0.0557	;
	*stateFreqs[6]	=	0.0277	;
	*stateFreqs[7]	=	0.0905	;
	*stateFreqs[9]	=	0.1675	;
	*stateFreqs[8]	=	0.0221	;
	*stateFreqs[10]	=	0.0561	;
	*stateFreqs[4]	=	0.0611	;
	*stateFreqs[12]	=	0.0536	;
	*stateFreqs[15]	=	0.0725	;
	*stateFreqs[16]	=	0.0870	;
	*stateFreqs[18]	=	0.0293	;
	*stateFreqs[19]	=	0.0340	;
	*stateFreqs[17]	=	0.0428	;
	}

void Model::SetMtRevAAFreqs(){
	*stateFreqs[0] 	=	0.0720	;
	*stateFreqs[14]	=	0.0190	;
	*stateFreqs[11]	=	0.0390	;
	*stateFreqs[2]	=	0.0190	;
	*stateFreqs[1]	=	0.0060	;
	*stateFreqs[13]	=	0.0250	;
	*stateFreqs[3]	=	0.0240	;
	*stateFreqs[5]	=	0.0560	;
	*stateFreqs[6]	=	0.0280	;
	*stateFreqs[7]	=	0.0880	;
	*stateFreqs[9]	=	0.1690	;
	*stateFreqs[8]	=	0.0230	;
	*stateFreqs[10]	=	0.0540	;
	*stateFreqs[4]	=	0.0610	;
	*stateFreqs[12]	=	0.0540	;
	*stateFreqs[15]	=	0.0720	;
	*stateFreqs[16]	=	0.0860	;
	*stateFreqs[18]	=	0.0290	;
	*stateFreqs[19]	=	0.0330	;
	*stateFreqs[17]	=	0.0430	;
	}

void Model::SetJonesAAFreqs(){
		*stateFreqs[0] =0.076748;
		*stateFreqs[14]=0.051691;
		*stateFreqs[11]=0.042645;
		*stateFreqs[2]=0.051544;
		*stateFreqs[1]=0.019803;
		*stateFreqs[13]=0.040752;
		*stateFreqs[3]=0.06183;
		*stateFreqs[5]=0.073152;
		*stateFreqs[6]=0.022944;
		*stateFreqs[7]=0.053761;
		*stateFreqs[9]=0.091904;
		*stateFreqs[8]=0.058676;
		*stateFreqs[10]=0.023826;
		*stateFreqs[4]=0.040126;
		*stateFreqs[12]=0.050901;
		*stateFreqs[15]=0.068765;
		*stateFreqs[16]=0.058565;
		*stateFreqs[18]=0.014261;
		*stateFreqs[19]=0.032101;
		*stateFreqs[17]=0.066005;
		}
		
void Model::SetDayhoffAAFreqs(){
	*stateFreqs[0]		=0.087127;
	*stateFreqs[14]	=0.040904;
	*stateFreqs[11]	=0.040432;
	*stateFreqs[2]		=0.046872;
	*stateFreqs[1]		=0.033474;
	*stateFreqs[13]	=0.038255;
	*stateFreqs[3]		=0.04953;
	*stateFreqs[5]		=0.088612;
	*stateFreqs[6]		=0.033618;
	*stateFreqs[7]		=0.036886;
	*stateFreqs[9]		=0.085357;
	*stateFreqs[8]		=0.080482;
	*stateFreqs[10]	=0.014753;
	*stateFreqs[4]		=0.039772;
	*stateFreqs[12]	=0.05068;
	*stateFreqs[15]	=0.069577;
	*stateFreqs[16]	=0.058542;
	*stateFreqs[18]	=0.010494;
	*stateFreqs[19]	=0.029916;
	*stateFreqs[17]	=0.064718;
	}		

void Model::SetWAGAAFreqs(){
	*stateFreqs[0]=0.0866279;
	*stateFreqs[14]=0.043972;
	*stateFreqs[11]=0.0390894;
	*stateFreqs[2]=0.0570451;
	*stateFreqs[1]=0.0193078;
	*stateFreqs[13]=0.0367281;
	*stateFreqs[3]=0.0580589;
	*stateFreqs[5]=0.0832518;
	*stateFreqs[6]=0.0244313;
	*stateFreqs[7]=0.048466;
	*stateFreqs[9]=0.086209;
	*stateFreqs[8]=0.0620286;
	*stateFreqs[10]=0.0195027;
	*stateFreqs[4]=0.0384319;
	*stateFreqs[12]=0.0457631;
	*stateFreqs[15]=0.0695179;
	*stateFreqs[16]=0.0610127;
	*stateFreqs[18]=0.0143859;
	*stateFreqs[19]=0.0352742;
	*stateFreqs[17]=0.0708956;
	}

int Model::PerformModelMutation(){
	if(paramsToMutate.empty()) return 0;
	BaseParameter *mut = SelectModelMutation();
	assert(mut != NULL);
	mut->Mutator(mutationShape);
	int retType;

	if(mut->Type() == RELATIVERATES){
		UpdateQMat();
		retType=Individual::rates;
		eigenDirty=true;
		}
	else if(mut->Type() == STATEFREQS){
		UpdateQMat();
		retType=Individual::pi;
		eigenDirty=true;
		}
	
	else if(mut->Type() == PROPORTIONINVARIANT){
		//this max checking should really be rolled into the parameter class
		*propInvar = (*propInvar > maxPropInvar ? maxPropInvar : *propInvar);
		//the non invariant rates need to be rescaled even if there is only 1
		if(modSpec.IsFlexRateHet() == false) AdjustRateProportions();
		else NormalizeRates();
		retType=Individual::pinv;
		}
	else if(mut->Type() == ALPHASHAPE){
		DiscreteGamma(rateMults, rateProbs, *alpha);
		retType=Individual::alpha;
		}
	else if(mut->Type() == RATEPROPS || mut->Type() == RATEMULTS){
		//flex rates and omega muts come through here

		//enforce an ordering of the rate multipliers, so that they can't "cross" one another
		if(NRateCats() > 1)
			CheckAndCorrectRateOrdering();

		if(modSpec.IsFlexRateHet() == true)
			NormalizeRates();
		else if(modSpec.IsCodon()){
			//this normalization could really be taken care of in the mutator, but this general purpose
			//function does a better job of enforcing minimum values
			NormalizeSumConstrainedValues(&omegaProbs[0], NRateCats(), ONE_POINT_ZERO, 1.0e-5, -1);
			//eigen stuff needs to be recalced for changes to nonsynonymous rates
			eigenDirty = true;
			}
		retType=Individual::alpha;
		}
	return retType;
	}


BaseParameter *Model::SelectModelMutation(){
	CalcMutationProbsFromWeights();
	if(paramsToMutate.empty() == true) return NULL;
	FLOAT_TYPE r=rnd.uniform();
	vector<BaseParameter*>::iterator it;
	for(it=paramsToMutate.begin();it!=paramsToMutate.end();it++){
		if((*it)->GetProb() > r) return *it;
		}
	it--;
	return *it;
	}

void Model::CalcMutationProbsFromWeights(){
	FLOAT_TYPE tot=ZERO_POINT_ZERO, running=ZERO_POINT_ZERO;
	for(vector<BaseParameter*>::iterator it=paramsToMutate.begin();it!=paramsToMutate.end();it++){
		tot += (*it)->GetWeight();
		}
	for(vector<BaseParameter*>::iterator it=paramsToMutate.begin();it!=paramsToMutate.end();it++){
		running += (*it)->GetWeight() / tot;
		(*it)->SetProb(running);
		}
	}

/*
void Model::OutputBinaryFormattedModel(OUTPUT_CLASS &out) const{
	FLOAT_TYPE *r = new FLOAT_TYPE;
	for(int i=0;i<5;i++){
		*r = Rates(i);
		out.write((char *) r, sizeof(FLOAT_TYPE));
		}
	for(int i=0;i<NStates();i++){
		*r = StateFreq(i);
		out.write((char *) r, sizeof(FLOAT_TYPE));
		}
	
	if(modSpec.flexRates==true){
		for(int i=0;i<NRateCats();i++){
			out.write((char *) &rateMults[i], sizeof(FLOAT_TYPE));
			out.write((char *) &rateProbs[i], sizeof(FLOAT_TYPE));
			}
		}
	else{
		if(NRateCats()>1){
			*r = Alpha();
			out.write((char *) r, sizeof(FLOAT_TYPE));
			}
		}
	if(PropInvar()!=ZERO_POINT_ZERO){
		*r = PropInvar();
		out.write((char *) r, sizeof(FLOAT_TYPE));
		}
	delete r;
	}
*/

void Model::OutputBinaryFormattedModel(OUTPUT_CLASS &out) const{
	FLOAT_TYPE *r = new FLOAT_TYPE;
	if(modSpec.IsAminoAcid() == false || modSpec.IsUserSpecifiedRateMatrix() || modSpec.IsEstimateAAMatrix()){
		if(modSpec.IsAminoAcid())
			assert(NumRelRates() == 190);
		else
			assert(NumRelRates() == 6);
		for(int i=0;i<NumRelRates();i++){
			*r = Rates(i);
			out.WRITE_TO_FILE(r, sizeof(FLOAT_TYPE), 1);
			}
		}
	//for codon models, output omega(s)
	if(modSpec.IsCodon()){
		for(int i=0;i<omegas.size();i++){
			*r = *omegas[i];
			out.WRITE_TO_FILE(r, sizeof(FLOAT_TYPE), 1);
			*r = *omegaProbs[i];
			out.WRITE_TO_FILE(r, sizeof(FLOAT_TYPE), 1);
			}
		}

	//these may not actually be free params, but output them anyway
	for(int i=0;i<NStates();i++){
		*r = StateFreq(i);
		out.WRITE_TO_FILE(r, sizeof(FLOAT_TYPE), 1);
		}
	
	if(modSpec.IsFlexRateHet()){
		for(int i=0;i<NRateCats();i++){
			out.WRITE_TO_FILE(&rateMults[i], sizeof(FLOAT_TYPE), 1);
			out.WRITE_TO_FILE(&rateProbs[i], sizeof(FLOAT_TYPE), 1);
			}
		}
	else if(modSpec.IsGammaRateHet()){
		*r = Alpha();
		out.WRITE_TO_FILE(r, sizeof(FLOAT_TYPE), 1);
		}
	if(PropInvar()!=ZERO_POINT_ZERO){
		*r = PropInvar();
		out.WRITE_TO_FILE(r, sizeof(FLOAT_TYPE), 1);
		}
	delete r;
	}

void Model::ReadBinaryFormattedModel(FILE *in){
	if(modSpec.IsAminoAcid() == false || modSpec.IsUserSpecifiedRateMatrix() || modSpec.IsEstimateAAMatrix()){
		if(modSpec.IsAminoAcid())
			assert(NumRelRates() == 190);
		else
			assert(NumRelRates() == 6);
		FLOAT_TYPE *r = new FLOAT_TYPE[NumRelRates()];
		for(int i=0;i<NumRelRates();i++){
			assert(ferror(in) == false);
			fread(r+i, sizeof(FLOAT_TYPE), 1, in);
			}
		SetRmat(r, false, false);
		delete []r;
		}

	if(modSpec.IsCodon()){
		FLOAT_TYPE o;
		for(int i=0;i<omegas.size();i++){
			fread(&o, sizeof(FLOAT_TYPE), 1, in);
			*omegas[i] = o;
			fread(&o, sizeof(FLOAT_TYPE), 1, in);
			*omegaProbs[i] = o;
			}
		}	

	FLOAT_TYPE *b = new FLOAT_TYPE[NStates()];
	for(int i=0;i<NStates();i++){
		fread((char*) &(b[i]), sizeof(FLOAT_TYPE), 1, in);
		}
	SetPis(b, false, false);
	delete []b;

	if(modSpec.IsFlexRateHet()){
		for(int i=0;i<NRateCats();i++){
			fread((char*) &(rateMults[i]), sizeof(FLOAT_TYPE), 1, in);
			fread((char*) &(rateProbs[i]), sizeof(FLOAT_TYPE), 1, in);
			}
		}
	else{
		if(modSpec.IsGammaRateHet()){
			FLOAT_TYPE a;
			assert(ferror(in) == false);
			fread((char*) &a, sizeof(FLOAT_TYPE), 1, in);
			SetAlpha(a, false);
			}
		}
	if(PropInvar()!=ZERO_POINT_ZERO){
		FLOAT_TYPE p;
		fread((char*) &p, sizeof(FLOAT_TYPE), 1, in);
		SetPinv(p, false);
		}
	}

void Model::MultiplyByJonesAAMatrix(){
	int modNum=0;
	MODEL_FLOAT **qmatOffset = qmat[modNum];

	qmatOffset[0][1] *= 0.056; qmatOffset[1][0] *= 0.056; qmatOffset[0][2] *= 0.081; qmatOffset[2][0] *= 0.081; qmatOffset[0][3] *= 0.105; qmatOffset[3][0] *= 0.105; 
	qmatOffset[0][4] *= 0.015; qmatOffset[4][0] *= 0.015; qmatOffset[0][5] *= 0.179; qmatOffset[5][0] *= 0.179; qmatOffset[0][6] *= 0.027; qmatOffset[6][0] *= 0.027; 
	qmatOffset[0][7] *= 0.036; qmatOffset[7][0] *= 0.036; qmatOffset[0][8] *= 0.035; qmatOffset[8][0] *= 0.035; qmatOffset[0][9] *= 0.03; qmatOffset[9][0] *= 0.03; 
	qmatOffset[0][10] *= 0.054; qmatOffset[10][0] *= 0.054; qmatOffset[0][11] *= 0.054; qmatOffset[11][0] *= 0.054; qmatOffset[0][12] *= 0.194; qmatOffset[12][0] *= 0.194; 
	qmatOffset[0][13] *= 0.057; qmatOffset[13][0] *= 0.057; qmatOffset[0][14] *= 0.058; qmatOffset[14][0] *= 0.058; qmatOffset[0][15] *= 0.378; qmatOffset[15][0] *= 0.378; 
	qmatOffset[0][16] *= 0.475; qmatOffset[16][0] *= 0.475; qmatOffset[0][17] *= 0.298; qmatOffset[17][0] *= 0.298; qmatOffset[0][18] *= 0.009; qmatOffset[18][0] *= 0.009; 
	qmatOffset[0][19] *= 0.011; qmatOffset[19][0] *= 0.011; qmatOffset[1][2] *= 0.01; qmatOffset[2][1] *= 0.01; qmatOffset[1][3] *= 0.005; qmatOffset[3][1] *= 0.005; 
	qmatOffset[1][4] *= 0.078; qmatOffset[4][1] *= 0.078; qmatOffset[1][5] *= 0.059; qmatOffset[5][1] *= 0.059; qmatOffset[1][6] *= 0.069; qmatOffset[6][1] *= 0.069; 
	qmatOffset[1][7] *= 0.017; qmatOffset[7][1] *= 0.017; qmatOffset[1][8] *= 0.007; qmatOffset[8][1] *= 0.007; qmatOffset[1][9] *= 0.023; qmatOffset[9][1] *= 0.023; 
	qmatOffset[1][10] *= 0.031; qmatOffset[10][1] *= 0.031; qmatOffset[1][11] *= 0.034; qmatOffset[11][1] *= 0.034; qmatOffset[1][12] *= 0.014; qmatOffset[12][1] *= 0.014; 
	qmatOffset[1][13] *= 0.009; qmatOffset[13][1] *= 0.009; qmatOffset[1][14] *= 0.113; qmatOffset[14][1] *= 0.113; qmatOffset[1][15] *= 0.223; qmatOffset[15][1] *= 0.223; 
	qmatOffset[1][16] *= 0.042; qmatOffset[16][1] *= 0.042; qmatOffset[1][17] *= 0.062; qmatOffset[17][1] *= 0.062; qmatOffset[1][18] *= 0.115; qmatOffset[18][1] *= 0.115; 
	qmatOffset[1][19] *= 0.209; qmatOffset[19][1] *= 0.209; qmatOffset[2][3] *= 0.767; qmatOffset[3][2] *= 0.767; qmatOffset[2][4] *= 0.004; qmatOffset[4][2] *= 0.004; 
	qmatOffset[2][5] *= 0.13; qmatOffset[5][2] *= 0.13; qmatOffset[2][6] *= 0.112; qmatOffset[6][2] *= 0.112; qmatOffset[2][7] *= 0.011; qmatOffset[7][2] *= 0.011; 
	qmatOffset[2][8] *= 0.026; qmatOffset[8][2] *= 0.026; qmatOffset[2][9] *= 0.007; qmatOffset[9][2] *= 0.007; qmatOffset[2][10] *= 0.015; qmatOffset[10][2] *= 0.015; 
	qmatOffset[2][11] *= 0.528; qmatOffset[11][2] *= 0.528; qmatOffset[2][12] *= 0.015; qmatOffset[12][2] *= 0.015; qmatOffset[2][13] *= 0.049; qmatOffset[13][2] *= 0.049; 
	qmatOffset[2][14] *= 0.016; qmatOffset[14][2] *= 0.016; qmatOffset[2][15] *= 0.059; qmatOffset[15][2] *= 0.059; qmatOffset[2][16] *= 0.038; qmatOffset[16][2] *= 0.038; 
	qmatOffset[2][17] *= 0.031; qmatOffset[17][2] *= 0.031; qmatOffset[2][18] *= 0.004; qmatOffset[18][2] *= 0.004; qmatOffset[2][19] *= 0.046; qmatOffset[19][2] *= 0.046; 
	qmatOffset[3][4] *= 0.005; qmatOffset[4][3] *= 0.005; qmatOffset[3][5] *= 0.119; qmatOffset[5][3] *= 0.119; qmatOffset[3][6] *= 0.026; qmatOffset[6][3] *= 0.026; 
	qmatOffset[3][7] *= 0.012; qmatOffset[7][3] *= 0.012; qmatOffset[3][8] *= 0.181; qmatOffset[8][3] *= 0.181; qmatOffset[3][9] *= 0.009; qmatOffset[9][3] *= 0.009; 
	qmatOffset[3][10] *= 0.018; qmatOffset[10][3] *= 0.018; qmatOffset[3][11] *= 0.058; qmatOffset[11][3] *= 0.058; qmatOffset[3][12] *= 0.018; qmatOffset[12][3] *= 0.018; 
	qmatOffset[3][13] *= 0.323; qmatOffset[13][3] *= 0.323; qmatOffset[3][14] *= 0.029; qmatOffset[14][3] *= 0.029; qmatOffset[3][15] *= 0.03; qmatOffset[15][3] *= 0.03; 
	qmatOffset[3][16] *= 0.032; qmatOffset[16][3] *= 0.032; qmatOffset[3][17] *= 0.045; qmatOffset[17][3] *= 0.045; qmatOffset[3][18] *= 0.01; qmatOffset[18][3] *= 0.01; 
	qmatOffset[3][19] *= 0.007; qmatOffset[19][3] *= 0.007; qmatOffset[4][5] *= 0.005; qmatOffset[5][4] *= 0.005; qmatOffset[4][6] *= 0.04; qmatOffset[6][4] *= 0.04; 
	qmatOffset[4][7] *= 0.089; qmatOffset[7][4] *= 0.089; qmatOffset[4][8] *= 0.004; qmatOffset[8][4] *= 0.004; qmatOffset[4][9] *= 0.248; qmatOffset[9][4] *= 0.248; 
	qmatOffset[4][10] *= 0.043; qmatOffset[10][4] *= 0.043; qmatOffset[4][11] *= 0.01; qmatOffset[11][4] *= 0.01; qmatOffset[4][12] *= 0.017; qmatOffset[12][4] *= 0.017; 
	qmatOffset[4][13] *= 0.004; qmatOffset[13][4] *= 0.004; qmatOffset[4][14] *= 0.005; qmatOffset[14][4] *= 0.005; qmatOffset[4][15] *= 0.092; qmatOffset[15][4] *= 0.092; 
	qmatOffset[4][16] *= 0.012; qmatOffset[16][4] *= 0.012; qmatOffset[4][17] *= 0.062; qmatOffset[17][4] *= 0.062; qmatOffset[4][18] *= 0.053; qmatOffset[18][4] *= 0.053; 
	qmatOffset[4][19] *= 0.536; qmatOffset[19][4] *= 0.536; qmatOffset[5][6] *= 0.023; qmatOffset[6][5] *= 0.023; qmatOffset[5][7] *= 0.006; qmatOffset[7][5] *= 0.006; 
	qmatOffset[5][8] *= 0.027; qmatOffset[8][5] *= 0.027; qmatOffset[5][9] *= 0.006; qmatOffset[9][5] *= 0.006; qmatOffset[5][10] *= 0.014; qmatOffset[10][5] *= 0.014;
	qmatOffset[5][11] *= 0.081; qmatOffset[11][5] *= 0.081; qmatOffset[5][12] *= 0.024; qmatOffset[12][5] *= 0.024; qmatOffset[5][13] *= 0.026; qmatOffset[13][5] *= 0.026; 
	qmatOffset[5][14] *= 0.137; qmatOffset[14][5] *= 0.137; qmatOffset[5][15] *= 0.201; qmatOffset[15][5] *= 0.201; qmatOffset[5][16] *= 0.033; qmatOffset[16][5] *= 0.033; 
	qmatOffset[5][17] *= 0.047; qmatOffset[17][5] *= 0.047; qmatOffset[5][18] *= 0.055; qmatOffset[18][5] *= 0.055; qmatOffset[5][19] *= 0.008; qmatOffset[19][5] *= 0.008; 
	qmatOffset[6][7] *= 0.016; qmatOffset[7][6] *= 0.016; qmatOffset[6][8] *= 0.045; qmatOffset[8][6] *= 0.045; qmatOffset[6][9] *= 0.056; qmatOffset[9][6] *= 0.056; 
	qmatOffset[6][10] *= 0.033; qmatOffset[10][6] *= 0.033; qmatOffset[6][11] *= 0.391; qmatOffset[11][6] *= 0.391; qmatOffset[6][12] *= 0.115; qmatOffset[12][6] *= 0.115; 
	qmatOffset[6][13] *= 0.597; qmatOffset[13][6] *= 0.597; qmatOffset[6][14] *= 0.328; qmatOffset[14][6] *= 0.328; qmatOffset[6][15] *= 0.073; qmatOffset[15][6] *= 0.073; 
	qmatOffset[6][16] *= 0.046; qmatOffset[16][6] *= 0.046; qmatOffset[6][17] *= 0.011; qmatOffset[17][6] *= 0.011; qmatOffset[6][18] *= 0.008; qmatOffset[18][6] *= 0.008; 
	qmatOffset[6][19] *= 0.573; qmatOffset[19][6] *= 0.573; qmatOffset[7][8] *= 0.021; qmatOffset[8][7] *= 0.021; qmatOffset[7][9] *= 0.229; qmatOffset[9][7] *= 0.229; 
	qmatOffset[7][10] *= 0.479; qmatOffset[10][7] *= 0.479; qmatOffset[7][11] *= 0.047; qmatOffset[11][7] *= 0.047; qmatOffset[7][12] *= 0.01; qmatOffset[12][7] *= 0.01; 
	qmatOffset[7][13] *= 0.009; qmatOffset[13][7] *= 0.009; qmatOffset[7][14] *= 0.022; qmatOffset[14][7] *= 0.022; qmatOffset[7][15] *= 0.04; qmatOffset[15][7] *= 0.04; 
	qmatOffset[7][16] *= 0.245; qmatOffset[16][7] *= 0.245; qmatOffset[7][17] *= 0.961; qmatOffset[17][7] *= 0.961; qmatOffset[7][18] *= 0.009; qmatOffset[18][7] *= 0.009; 
	qmatOffset[7][19] *= 0.032; qmatOffset[19][7] *= 0.032; qmatOffset[8][9] *= 0.014; qmatOffset[9][8] *= 0.014; qmatOffset[8][10] *= 0.065; qmatOffset[10][8] *= 0.065; 
	qmatOffset[8][11] *= 0.263; qmatOffset[11][8] *= 0.263; qmatOffset[8][12] *= 0.021; qmatOffset[12][8] *= 0.021; qmatOffset[8][13] *= 0.292; qmatOffset[13][8] *= 0.292; 
	qmatOffset[8][14] *= 0.646; qmatOffset[14][8] *= 0.646; qmatOffset[8][15] *= 0.047; qmatOffset[15][8] *= 0.047; qmatOffset[8][16] *= 0.103; qmatOffset[16][8] *= 0.103; 
	qmatOffset[8][17] *= 0.014; qmatOffset[17][8] *= 0.014; qmatOffset[8][18] *= 0.01; qmatOffset[18][8] *= 0.01; qmatOffset[8][19] *= 0.008; qmatOffset[19][8] *= 0.008; 
	qmatOffset[9][10] *= 0.388; qmatOffset[10][9] *= 0.388; qmatOffset[9][11] *= 0.012; qmatOffset[11][9] *= 0.012; qmatOffset[9][12] *= 0.102; qmatOffset[12][9] *= 0.102; 
	qmatOffset[9][13] *= 0.072; qmatOffset[13][9] *= 0.072; qmatOffset[9][14] *= 0.038; qmatOffset[14][9] *= 0.038; qmatOffset[9][15] *= 0.059; qmatOffset[15][9] *= 0.059; 
	qmatOffset[9][16] *= 0.025; qmatOffset[16][9] *= 0.025; qmatOffset[9][17] *= 0.18; qmatOffset[17][9] *= 0.18; qmatOffset[9][18] *= 0.052; qmatOffset[18][9] *= 0.052; 
	qmatOffset[9][19] *= 0.024; qmatOffset[19][9] *= 0.024; qmatOffset[10][11] *= 0.03; qmatOffset[11][10] *= 0.03; qmatOffset[10][12] *= 0.016; qmatOffset[12][10] *= 0.016; 
	qmatOffset[10][13] *= 0.043; qmatOffset[13][10] *= 0.043; qmatOffset[10][14] *= 0.044; qmatOffset[14][10] *= 0.044; qmatOffset[10][15] *= 0.029; qmatOffset[15][10] *= 0.029; 
	qmatOffset[10][16] *= 0.226; qmatOffset[16][10] *= 0.226; qmatOffset[10][17] *= 0.323; qmatOffset[17][10] *= 0.323; qmatOffset[10][18] *= 0.024; qmatOffset[18][10] *= 0.024; 
	qmatOffset[10][19] *= 0.018; qmatOffset[19][10] *= 0.018; qmatOffset[11][12] *= 0.015; qmatOffset[12][11] *= 0.015; qmatOffset[11][13] *= 0.086; qmatOffset[13][11] *= 0.086; 
	qmatOffset[11][14] *= 0.045; qmatOffset[14][11] *= 0.045; qmatOffset[11][15] *= 0.503; qmatOffset[15][11] *= 0.503; qmatOffset[11][16] *= 0.232; qmatOffset[16][11] *= 0.232; 
	qmatOffset[11][17] *= 0.016; qmatOffset[17][11] *= 0.016; qmatOffset[11][18] *= 0.008; qmatOffset[18][11] *= 0.008; qmatOffset[11][19] *= 0.07; qmatOffset[19][11] *= 0.07; 
	qmatOffset[12][13] *= 0.164; qmatOffset[13][12] *= 0.164; qmatOffset[12][14] *= 0.074; qmatOffset[14][12] *= 0.074; qmatOffset[12][15] *= 0.285; qmatOffset[15][12] *= 0.285; 
	qmatOffset[12][16] *= 0.118; qmatOffset[16][12] *= 0.118; qmatOffset[12][17] *= 0.023; qmatOffset[17][12] *= 0.023; qmatOffset[12][18] *= 0.006; qmatOffset[18][12] *= 0.006; 
	qmatOffset[12][19] *= 0.01; qmatOffset[19][12] *= 0.01; qmatOffset[13][14] *= 0.31; qmatOffset[14][13] *= 0.31; qmatOffset[13][15] *= 0.053; qmatOffset[15][13] *= 0.053; 
	qmatOffset[13][16] *= 0.051; qmatOffset[16][13] *= 0.051; qmatOffset[13][17] *= 0.02; qmatOffset[17][13] *= 0.02; qmatOffset[13][18] *= 0.018; qmatOffset[18][13] *= 0.018; 
	qmatOffset[13][19] *= 0.024; qmatOffset[19][13] *= 0.024; qmatOffset[14][15] *= 0.101; qmatOffset[15][14] *= 0.101; qmatOffset[14][16] *= 0.064; qmatOffset[16][14] *= 0.064; 
	qmatOffset[14][17] *= 0.017; qmatOffset[17][14] *= 0.017; qmatOffset[14][18] *= 0.126; qmatOffset[18][14] *= 0.126; qmatOffset[14][19] *= 0.02; qmatOffset[19][14] *= 0.02; 
	qmatOffset[15][16] *= 0.477; qmatOffset[16][15] *= 0.477; qmatOffset[15][17] *= 0.038; qmatOffset[17][15] *= 0.038; qmatOffset[15][18] *= 0.035; qmatOffset[18][15] *= 0.035; 
	qmatOffset[15][19] *= 0.063; qmatOffset[19][15] *= 0.063; qmatOffset[16][17] *= 0.112; qmatOffset[17][16] *= 0.112; qmatOffset[16][18] *= 0.012; qmatOffset[18][16] *= 0.012; 
	qmatOffset[16][19] *= 0.021; qmatOffset[19][16] *= 0.021; qmatOffset[17][18] *= 0.025; qmatOffset[18][17] *= 0.025; qmatOffset[17][19] *= 0.016; qmatOffset[19][17] *= 0.016; 
	qmatOffset[18][19] *= 0.071; qmatOffset[19][18] *= 0.071;
	}

void Model::MultiplyByMtMamAAMatrix(){
	int modNum=0;
	MODEL_FLOAT **qmatOffset = qmat[modNum];

	qmatOffset [ 0 ][ 14 ] *= 0.0337 ; qmatOffset [ 14 ][ 0 ] *= 0.0337 ;
	qmatOffset [ 0 ][ 11 ] *= 0.0021 ; qmatOffset [ 11 ][ 0 ] *= 0.0021 ;
	qmatOffset [ 0 ][ 2 ] *= 0.0116 ; qmatOffset [ 2 ][ 0 ] *= 0.0116 ;
	qmatOffset [ 0 ][ 1 ] *= 0.0000 ; qmatOffset [ 1 ][ 0 ] *= 0.0000 ;
	qmatOffset [ 0 ][ 13 ] *= 0.0000 ; qmatOffset [ 13 ][ 0 ] *= 0.0000 ;
	qmatOffset [ 0 ][ 3 ] *= 0.0000 ; qmatOffset [ 3 ][ 0 ] *= 0.0000 ;
	qmatOffset [ 0 ][ 5 ] *= 0.0821 ; qmatOffset [ 5 ][ 0 ] *= 0.0821 ;
	qmatOffset [ 0 ][ 6 ] *= 0.0084 ; qmatOffset [ 6 ][ 0 ] *= 0.0084 ;
	qmatOffset [ 0 ][ 7 ] *= 0.0790 ; qmatOffset [ 7 ][ 0 ] *= 0.0790 ;
	qmatOffset [ 0 ][ 9 ] *= 0.0221 ; qmatOffset [ 9 ][ 0 ] *= 0.0221 ;
	qmatOffset [ 0 ][ 8 ] *= 0.0000 ; qmatOffset [ 8 ][ 0 ] *= 0.0000 ;
	qmatOffset [ 0 ][ 10 ] *= 0.0800 ; qmatOffset [ 10 ][ 0 ] *= 0.0800 ;
	qmatOffset [ 0 ][ 4 ] *= 0.0000 ; qmatOffset [ 4 ][ 0 ] *= 0.0000 ;
	qmatOffset [ 0 ][ 12 ] *= 0.0558 ; qmatOffset [ 12 ][ 0 ] *= 0.0558 ;
	qmatOffset [ 0 ][ 15 ] *= 0.3601 ; qmatOffset [ 15 ][ 0 ] *= 0.3601 ;
	qmatOffset [ 0 ][ 16 ] *= 0.7170 ; qmatOffset [ 16 ][ 0 ] *= 0.7170 ;
	qmatOffset [ 0 ][ 18 ] *= 0.0053 ; qmatOffset [ 18 ][ 0 ] *= 0.0053 ;
	qmatOffset [ 0 ][ 19 ] *= 0.0000 ; qmatOffset [ 19 ][ 0 ] *= 0.0000 ;
	qmatOffset [ 0 ][ 17 ] *= 0.4191 ; qmatOffset [ 17 ][ 0 ] *= 0.4191 ;
	qmatOffset [ 14 ][ 11 ] *= 0.0042 ; qmatOffset [ 11 ][ 14 ] *= 0.0042 ;
	qmatOffset [ 14 ][ 2 ] *= 0.0000 ; qmatOffset [ 2 ][ 14 ] *= 0.0000 ;
	qmatOffset [ 14 ][ 1 ] *= 0.1958 ; qmatOffset [ 1 ][ 14 ] *= 0.1958 ;
	qmatOffset [ 14 ][ 13 ] *= 0.2590 ; qmatOffset [ 13 ][ 14 ] *= 0.2590 ;
	qmatOffset [ 14 ][ 3 ] *= 0.0000 ; qmatOffset [ 3 ][ 14 ] *= 0.0000 ;
	qmatOffset [ 14 ][ 5 ] *= 0.0190 ; qmatOffset [ 5 ][ 14 ] *= 0.0190 ;
	qmatOffset [ 14 ][ 6 ] *= 0.2443 ; qmatOffset [ 6 ][ 14 ] *= 0.2443 ;
	qmatOffset [ 14 ][ 7 ] *= 0.0000 ; qmatOffset [ 7 ][ 14 ] *= 0.0000 ;
	qmatOffset [ 14 ][ 9 ] *= 0.0063 ; qmatOffset [ 9 ][ 14 ] *= 0.0063 ;
	qmatOffset [ 14 ][ 8 ] *= 0.0526 ; qmatOffset [ 8 ][ 14 ] *= 0.0526 ;
	qmatOffset [ 14 ][ 10 ] *= 0.0000 ; qmatOffset [ 10 ][ 14 ] *= 0.0000 ;
	qmatOffset [ 14 ][ 4 ] *= 0.0000 ; qmatOffset [ 4 ][ 14 ] *= 0.0000 ;
	qmatOffset [ 14 ][ 12 ] *= 0.0095 ; qmatOffset [ 12 ][ 14 ] *= 0.0095 ;
	qmatOffset [ 14 ][ 15 ] *= 0.0032 ; qmatOffset [ 15 ][ 14 ] *= 0.0032 ;
	qmatOffset [ 14 ][ 16 ] *= 0.0000 ; qmatOffset [ 16 ][ 14 ] *= 0.0000 ;
	qmatOffset [ 14 ][ 18 ] *= 0.0168 ; qmatOffset [ 18 ][ 14 ] *= 0.0168 ;
	qmatOffset [ 14 ][ 19 ] *= 0.0000 ; qmatOffset [ 19 ][ 14 ] *= 0.0000 ;
	qmatOffset [ 14 ][ 17 ] *= 0.0000 ; qmatOffset [ 17 ][ 14 ] *= 0.0000 ;
	qmatOffset [ 11 ][ 2 ] *= 0.9097 ; qmatOffset [ 2 ][ 11 ] *= 0.9097 ;
	qmatOffset [ 11 ][ 1 ] *= 0.0000 ; qmatOffset [ 1 ][ 11 ] *= 0.0000 ;
	qmatOffset [ 11 ][ 13 ] *= 0.0084 ; qmatOffset [ 13 ][ 11 ] *= 0.0084 ;
	qmatOffset [ 11 ][ 3 ] *= 0.0000 ; qmatOffset [ 3 ][ 11 ] *= 0.0000 ;
	qmatOffset [ 11 ][ 5 ] *= 0.0495 ; qmatOffset [ 5 ][ 11 ] *= 0.0495 ;
	qmatOffset [ 11 ][ 6 ] *= 0.4822 ; qmatOffset [ 6 ][ 11 ] *= 0.4822 ;
	qmatOffset [ 11 ][ 7 ] *= 0.0200 ; qmatOffset [ 7 ][ 11 ] *= 0.0200 ;
	qmatOffset [ 11 ][ 9 ] *= 0.0000 ; qmatOffset [ 9 ][ 11 ] *= 0.0000 ;
	qmatOffset [ 11 ][ 8 ] *= 0.4296 ; qmatOffset [ 8 ][ 11 ] *= 0.4296 ;
	qmatOffset [ 11 ][ 10 ] *= 0.0221 ; qmatOffset [ 10 ][ 11 ] *= 0.0221 ;
	qmatOffset [ 11 ][ 4 ] *= 0.0063 ; qmatOffset [ 4 ][ 11 ] *= 0.0063 ;
	qmatOffset [ 11 ][ 12 ] *= 0.0347 ; qmatOffset [ 12 ][ 11 ] *= 0.0347 ;
	qmatOffset [ 11 ][ 15 ] *= 0.4696 ; qmatOffset [ 15 ][ 11 ] *= 0.4696 ;
	qmatOffset [ 11 ][ 16 ] *= 0.1158 ; qmatOffset [ 16 ][ 11 ] *= 0.1158 ;
	qmatOffset [ 11 ][ 18 ] *= 0.0063 ; qmatOffset [ 18 ][ 11 ] *= 0.0063 ;
	qmatOffset [ 11 ][ 19 ] *= 0.1643 ; qmatOffset [ 19 ][ 11 ] *= 0.1643 ;
	qmatOffset [ 11 ][ 17 ] *= 0.0000 ; qmatOffset [ 17 ][ 11 ] *= 0.0000 ;
	qmatOffset [ 2 ][ 1 ] *= 0.0000 ; qmatOffset [ 1 ][ 2 ] *= 0.0000 ;
	qmatOffset [ 2 ][ 13 ] *= 0.0516 ; qmatOffset [ 13 ][ 2 ] *= 0.0516 ;
	qmatOffset [ 2 ][ 3 ] *= 0.5991 ; qmatOffset [ 3 ][ 2 ] *= 0.5991 ;
	qmatOffset [ 2 ][ 5 ] *= 0.0832 ; qmatOffset [ 5 ][ 2 ] *= 0.0832 ;
	qmatOffset [ 2 ][ 6 ] *= 0.0116 ; qmatOffset [ 6 ][ 2 ] *= 0.0116 ;
	qmatOffset [ 2 ][ 7 ] *= 0.0000 ; qmatOffset [ 7 ][ 2 ] *= 0.0000 ;
	qmatOffset [ 2 ][ 9 ] *= 0.0000 ; qmatOffset [ 9 ][ 2 ] *= 0.0000 ;
	qmatOffset [ 2 ][ 8 ] *= 0.0000 ; qmatOffset [ 8 ][ 2 ] *= 0.0000 ;
	qmatOffset [ 2 ][ 10 ] *= 0.0000 ; qmatOffset [ 10 ][ 2 ] *= 0.0000 ;
	qmatOffset [ 2 ][ 4 ] *= 0.0053 ; qmatOffset [ 4 ][ 2 ] *= 0.0053 ;
	qmatOffset [ 2 ][ 12 ] *= 0.0021 ; qmatOffset [ 12 ][ 2 ] *= 0.0021 ;
	qmatOffset [ 2 ][ 15 ] *= 0.0168 ; qmatOffset [ 15 ][ 2 ] *= 0.0168 ;
	qmatOffset [ 2 ][ 16 ] *= 0.0000 ; qmatOffset [ 16 ][ 2 ] *= 0.0000 ;
	qmatOffset [ 2 ][ 18 ] *= 0.0000 ; qmatOffset [ 18 ][ 2 ] *= 0.0000 ;
	qmatOffset [ 2 ][ 19 ] *= 0.0000 ; qmatOffset [ 19 ][ 2 ] *= 0.0000 ;
	qmatOffset [ 2 ][ 17 ] *= 0.0105 ; qmatOffset [ 17 ][ 2 ] *= 0.0105 ;
	qmatOffset [ 1 ][ 13 ] *= 0.0000 ; qmatOffset [ 13 ][ 1 ] *= 0.0000 ;
	qmatOffset [ 1 ][ 3 ] *= 0.0000 ; qmatOffset [ 3 ][ 1 ] *= 0.0000 ;
	qmatOffset [ 1 ][ 5 ] *= 0.0000 ; qmatOffset [ 5 ][ 1 ] *= 0.0000 ;
	qmatOffset [ 1 ][ 6 ] *= 0.3211 ; qmatOffset [ 6 ][ 1 ] *= 0.3211 ;
	qmatOffset [ 1 ][ 7 ] *= 0.0432 ; qmatOffset [ 7 ][ 1 ] *= 0.0432 ;
	qmatOffset [ 1 ][ 9 ] *= 0.0284 ; qmatOffset [ 9 ][ 1 ] *= 0.0284 ;
	qmatOffset [ 1 ][ 8 ] *= 0.0000 ; qmatOffset [ 8 ][ 1 ] *= 0.0000 ;
	qmatOffset [ 1 ][ 10 ] *= 0.0000 ; qmatOffset [ 10 ][ 1 ] *= 0.0000 ;
	qmatOffset [ 1 ][ 4 ] *= 0.0074 ; qmatOffset [ 4 ][ 1 ] *= 0.0074 ;
	qmatOffset [ 1 ][ 12 ] *= 0.0000 ; qmatOffset [ 12 ][ 1 ] *= 0.0000 ;
	qmatOffset [ 1 ][ 15 ] *= 0.3654 ; qmatOffset [ 15 ][ 1 ] *= 0.3654 ;
	qmatOffset [ 1 ][ 16 ] *= 0.1200 ; qmatOffset [ 16 ][ 1 ] *= 0.1200 ;
	qmatOffset [ 1 ][ 18 ] *= 0.0684 ; qmatOffset [ 18 ][ 1 ] *= 0.0684 ;
	qmatOffset [ 1 ][ 19 ] *= 0.5580 ; qmatOffset [ 19 ][ 1 ] *= 0.5580 ;
	qmatOffset [ 1 ][ 17 ] *= 0.0000 ; qmatOffset [ 17 ][ 1 ] *= 0.0000 ;
	qmatOffset [ 13 ][ 3 ] *= 0.2885 ; qmatOffset [ 3 ][ 13 ] *= 0.2885 ;
	qmatOffset [ 13 ][ 5 ] *= 0.0000 ; qmatOffset [ 5 ][ 13 ] *= 0.0000 ;
	qmatOffset [ 13 ][ 6 ] *= 0.5791 ; qmatOffset [ 6 ][ 13 ] *= 0.5791 ;
	qmatOffset [ 13 ][ 7 ] *= 0.0000 ; qmatOffset [ 7 ][ 13 ] *= 0.0000 ;
	qmatOffset [ 13 ][ 9 ] *= 0.0211 ; qmatOffset [ 9 ][ 13 ] *= 0.0211 ;
	qmatOffset [ 13 ][ 8 ] *= 0.2548 ; qmatOffset [ 8 ][ 13 ] *= 0.2548 ;
	qmatOffset [ 13 ][ 10 ] *= 0.0232 ; qmatOffset [ 10 ][ 13 ] *= 0.0232 ;
	qmatOffset [ 13 ][ 4 ] *= 0.0000 ; qmatOffset [ 4 ][ 13 ] *= 0.0000 ;
	qmatOffset [ 13 ][ 12 ] *= 0.0537 ; qmatOffset [ 12 ][ 13 ] *= 0.0537 ;
	qmatOffset [ 13 ][ 15 ] *= 0.0316 ; qmatOffset [ 15 ][ 13 ] *= 0.0316 ;
	qmatOffset [ 13 ][ 16 ] *= 0.0000 ; qmatOffset [ 16 ][ 13 ] *= 0.0000 ;
	qmatOffset [ 13 ][ 18 ] *= 0.0000 ; qmatOffset [ 18 ][ 13 ] *= 0.0000 ;
	qmatOffset [ 13 ][ 19 ] *= 0.0569 ; qmatOffset [ 19 ][ 13 ] *= 0.0569 ;
	qmatOffset [ 13 ][ 17 ] *= 0.0347 ; qmatOffset [ 17 ][ 13 ] *= 0.0347 ;
	qmatOffset [ 3 ][ 5 ] *= 0.0232 ; qmatOffset [ 5 ][ 3 ] *= 0.0232 ;
	qmatOffset [ 3 ][ 6 ] *= 0.0232 ; qmatOffset [ 6 ][ 3 ] *= 0.0232 ;
	qmatOffset [ 3 ][ 7 ] *= 0.0000 ; qmatOffset [ 7 ][ 3 ] *= 0.0000 ;
	qmatOffset [ 3 ][ 9 ] *= 0.0000 ; qmatOffset [ 9 ][ 3 ] *= 0.0000 ;
	qmatOffset [ 3 ][ 8 ] *= 0.2264 ; qmatOffset [ 8 ][ 3 ] *= 0.2264 ;
	qmatOffset [ 3 ][ 10 ] *= 0.0000 ; qmatOffset [ 10 ][ 3 ] *= 0.0000 ;
	qmatOffset [ 3 ][ 4 ] *= 0.0000 ; qmatOffset [ 4 ][ 3 ] *= 0.0000 ;
	qmatOffset [ 3 ][ 12 ] *= 0.0000 ; qmatOffset [ 12 ][ 3 ] *= 0.0000 ;
	qmatOffset [ 3 ][ 15 ] *= 0.0221 ; qmatOffset [ 15 ][ 3 ] *= 0.0221 ;
	qmatOffset [ 3 ][ 16 ] *= 0.0042 ; qmatOffset [ 16 ][ 3 ] *= 0.0042 ;
	qmatOffset [ 3 ][ 18 ] *= 0.0000 ; qmatOffset [ 18 ][ 3 ] *= 0.0000 ;
	qmatOffset [ 3 ][ 19 ] *= 0.0000 ; qmatOffset [ 19 ][ 3 ] *= 0.0000 ;
	qmatOffset [ 3 ][ 17 ] *= 0.0211 ; qmatOffset [ 17 ][ 3 ] *= 0.0211 ;
	qmatOffset [ 5 ][ 6 ] *= 0.0000 ; qmatOffset [ 6 ][ 5 ] *= 0.0000 ;
	qmatOffset [ 5 ][ 7 ] *= 0.0000 ; qmatOffset [ 7 ][ 5 ] *= 0.0000 ;
	qmatOffset [ 5 ][ 9 ] *= 0.0000 ; qmatOffset [ 9 ][ 5 ] *= 0.0000 ;
	qmatOffset [ 5 ][ 8 ] *= 0.0000 ; qmatOffset [ 8 ][ 5 ] *= 0.0000 ;
	qmatOffset [ 5 ][ 10 ] *= 0.0000 ; qmatOffset [ 10 ][ 5 ] *= 0.0000 ;
	qmatOffset [ 5 ][ 4 ] *= 0.0000 ; qmatOffset [ 4 ][ 5 ] *= 0.0000 ;
	qmatOffset [ 5 ][ 12 ] *= 0.0000 ; qmatOffset [ 12 ][ 5 ] *= 0.0000 ;
	qmatOffset [ 5 ][ 15 ] *= 0.1179 ; qmatOffset [ 15 ][ 5 ] *= 0.1179 ;
	qmatOffset [ 5 ][ 16 ] *= 0.0000 ; qmatOffset [ 16 ][ 5 ] *= 0.0000 ;
	qmatOffset [ 5 ][ 18 ] *= 0.0000 ; qmatOffset [ 18 ][ 5 ] *= 0.0000 ;
	qmatOffset [ 5 ][ 19 ] *= 0.0011 ; qmatOffset [ 19 ][ 5 ] *= 0.0011 ;
	qmatOffset [ 5 ][ 17 ] *= 0.0053 ; qmatOffset [ 17 ][ 5 ] *= 0.0053 ;
	qmatOffset [ 6 ][ 7 ] *= 0.0000 ; qmatOffset [ 7 ][ 6 ] *= 0.0000 ;
	qmatOffset [ 6 ][ 9 ] *= 0.0274 ; qmatOffset [ 9 ][ 6 ] *= 0.0274 ;
	qmatOffset [ 6 ][ 8 ] *= 0.0000 ; qmatOffset [ 8 ][ 6 ] *= 0.0000 ;
	qmatOffset [ 6 ][ 10 ] *= 0.0000 ; qmatOffset [ 10 ][ 6 ] *= 0.0000 ;
	qmatOffset [ 6 ][ 4 ] *= 0.0000 ; qmatOffset [ 4 ][ 6 ] *= 0.0000 ;
	qmatOffset [ 6 ][ 12 ] *= 0.0558 ; qmatOffset [ 12 ][ 6 ] *= 0.0558 ;
	qmatOffset [ 6 ][ 15 ] *= 0.0211 ; qmatOffset [ 15 ][ 6 ] *= 0.0211 ;
	qmatOffset [ 6 ][ 16 ] *= 0.0011 ; qmatOffset [ 16 ][ 6 ] *= 0.0011 ;
	qmatOffset [ 6 ][ 18 ] *= 0.0000 ; qmatOffset [ 18 ][ 6 ] *= 0.0000 ;
	qmatOffset [ 6 ][ 19 ] *= 1.6057 ; qmatOffset [ 19 ][ 6 ] *= 1.6057 ;
	qmatOffset [ 6 ][ 17 ] *= 0.0000 ; qmatOffset [ 17 ][ 6 ] *= 0.0000 ;
	qmatOffset [ 7 ][ 9 ] *= 0.2443 ; qmatOffset [ 9 ][ 7 ] *= 0.2443 ;
	qmatOffset [ 7 ][ 8 ] *= 0.0063 ; qmatOffset [ 8 ][ 7 ] *= 0.0063 ;
	qmatOffset [ 7 ][ 10 ] *= 0.3980 ; qmatOffset [ 10 ][ 7 ] *= 0.3980 ;
	qmatOffset [ 7 ][ 4 ] *= 0.0600 ; qmatOffset [ 4 ][ 7 ] *= 0.0600 ;
	qmatOffset [ 7 ][ 12 ] *= 0.0053 ; qmatOffset [ 12 ][ 7 ] *= 0.0053 ;
	qmatOffset [ 7 ][ 15 ] *= 0.0000 ; qmatOffset [ 15 ][ 7 ] *= 0.0000 ;
	qmatOffset [ 7 ][ 16 ] *= 0.3790 ; qmatOffset [ 16 ][ 7 ] *= 0.3790 ;
	qmatOffset [ 7 ][ 18 ] *= 0.0000 ; qmatOffset [ 18 ][ 7 ] *= 0.0000 ;
	qmatOffset [ 7 ][ 19 ] *= 0.0168 ; qmatOffset [ 19 ][ 7 ] *= 0.0168 ;
	qmatOffset [ 7 ][ 17 ] *= 2.3375 ; qmatOffset [ 17 ][ 7 ] *= 2.3375 ;
	qmatOffset [ 9 ][ 8 ] *= 0.0042 ; qmatOffset [ 8 ][ 9 ] *= 0.0042 ;
	qmatOffset [ 9 ][ 10 ] *= 0.6412 ; qmatOffset [ 10 ][ 9 ] *= 0.6412 ;
	qmatOffset [ 9 ][ 4 ] *= 0.2590 ; qmatOffset [ 4 ][ 9 ] *= 0.2590 ;
	qmatOffset [ 9 ][ 12 ] *= 0.0453 ; qmatOffset [ 12 ][ 9 ] *= 0.0453 ;
	qmatOffset [ 9 ][ 15 ] *= 0.0779 ; qmatOffset [ 15 ][ 9 ] *= 0.0779 ;
	qmatOffset [ 9 ][ 16 ] *= 0.0358 ; qmatOffset [ 16 ][ 9 ] *= 0.0358 ;
	qmatOffset [ 9 ][ 18 ] *= 0.0126 ; qmatOffset [ 18 ][ 9 ] *= 0.0126 ;
	qmatOffset [ 9 ][ 19 ] *= 0.0263 ; qmatOffset [ 19 ][ 9 ] *= 0.0263 ;
	qmatOffset [ 9 ][ 17 ] *= 0.1053 ; qmatOffset [ 17 ][ 9 ] *= 0.1053 ;
	qmatOffset [ 8 ][ 10 ] *= 0.0621 ; qmatOffset [ 10 ][ 8 ] *= 0.0621 ;
	qmatOffset [ 8 ][ 4 ] *= 0.0000 ; qmatOffset [ 4 ][ 8 ] *= 0.0000 ;
	qmatOffset [ 8 ][ 12 ] *= 0.0190 ; qmatOffset [ 12 ][ 8 ] *= 0.0190 ;
	qmatOffset [ 8 ][ 15 ] *= 0.0684 ; qmatOffset [ 15 ][ 8 ] *= 0.0684 ;
	qmatOffset [ 8 ][ 16 ] *= 0.0526 ; qmatOffset [ 16 ][ 8 ] *= 0.0526 ;
	qmatOffset [ 8 ][ 18 ] *= 0.0000 ; qmatOffset [ 18 ][ 8 ] *= 0.0000 ;
	qmatOffset [ 8 ][ 19 ] *= 0.0705 ; qmatOffset [ 19 ][ 8 ] *= 0.0705 ;
	qmatOffset [ 8 ][ 17 ] *= 0.0000 ; qmatOffset [ 17 ][ 8 ] *= 0.0000 ;
	qmatOffset [ 10 ][ 4 ] *= 0.0116 ; qmatOffset [ 4 ][ 10 ] *= 0.0116 ;
	qmatOffset [ 10 ][ 12 ] *= 0.0000 ; qmatOffset [ 12 ][ 10 ] *= 0.0000 ;
	qmatOffset [ 10 ][ 15 ] *= 0.0495 ; qmatOffset [ 15 ][ 10 ] *= 0.0495 ;
	qmatOffset [ 10 ][ 16 ] *= 0.7276 ; qmatOffset [ 16 ][ 10 ] *= 0.7276 ;
	qmatOffset [ 10 ][ 18 ] *= 0.0137 ; qmatOffset [ 18 ][ 10 ] *= 0.0137 ;
	qmatOffset [ 10 ][ 19 ] *= 0.0000 ; qmatOffset [ 19 ][ 10 ] *= 0.0000 ;
	qmatOffset [ 10 ][ 17 ] *= 0.8760 ; qmatOffset [ 17 ][ 10 ] *= 0.8760 ;
	qmatOffset [ 4 ][ 12 ] *= 0.0179 ; qmatOffset [ 12 ][ 4 ] *= 0.0179 ;
	qmatOffset [ 4 ][ 15 ] *= 0.0948 ; qmatOffset [ 15 ][ 4 ] *= 0.0948 ;
	qmatOffset [ 4 ][ 16 ] *= 0.0084 ; qmatOffset [ 16 ][ 4 ] *= 0.0084 ;
	qmatOffset [ 4 ][ 18 ] *= 0.0000 ; qmatOffset [ 18 ][ 4 ] *= 0.0000 ;
	qmatOffset [ 4 ][ 19 ] *= 0.7181 ; qmatOffset [ 19 ][ 4 ] *= 0.7181 ;
	qmatOffset [ 4 ][ 17 ] *= 0.0063 ; qmatOffset [ 17 ][ 4 ] *= 0.0063 ;
	qmatOffset [ 12 ][ 15 ] *= 0.2127 ; qmatOffset [ 15 ][ 12 ] *= 0.2127 ;
	qmatOffset [ 12 ][ 16 ] *= 0.0821 ; qmatOffset [ 16 ][ 12 ] *= 0.0821 ;
	qmatOffset [ 12 ][ 18 ] *= 0.0074 ; qmatOffset [ 18 ][ 12 ] *= 0.0074 ;
	qmatOffset [ 12 ][ 19 ] *= 0.0084 ; qmatOffset [ 19 ][ 12 ] *= 0.0084 ;
	qmatOffset [ 12 ][ 17 ] *= 0.0000 ; qmatOffset [ 17 ][ 12 ] *= 0.0000 ;
	qmatOffset [ 15 ][ 16 ] *= 0.6465 ; qmatOffset [ 16 ][ 15 ] *= 0.6465 ;
	qmatOffset [ 15 ][ 18 ] *= 0.0179 ; qmatOffset [ 18 ][ 15 ] *= 0.0179 ;
	qmatOffset [ 15 ][ 19 ] *= 0.1127 ; qmatOffset [ 19 ][ 15 ] *= 0.1127 ;
	qmatOffset [ 15 ][ 17 ] *= 0.0000 ; qmatOffset [ 17 ][ 15 ] *= 0.0000 ;
	qmatOffset [ 16 ][ 18 ] *= 0.0000 ; qmatOffset [ 18 ][ 16 ] *= 0.0000 ;
	qmatOffset [ 16 ][ 19 ] *= 0.0000 ; qmatOffset [ 19 ][ 16 ] *= 0.0000 ;
	qmatOffset [ 16 ][ 17 ] *= 0.2495 ; qmatOffset [ 17 ][ 16 ] *= 0.2495 ;
	qmatOffset [ 18 ][ 19 ] *= 0.0147 ; qmatOffset [ 19 ][ 18 ] *= 0.0147 ;
	qmatOffset [ 18 ][ 17 ] *= 0.0000 ; qmatOffset [ 17 ][ 18 ] *= 0.0000 ;
	qmatOffset [ 19 ][ 17 ] *= 0.0000 ; qmatOffset [ 17 ][ 19 ] *= 0.0000 ;
	}

void Model::MultiplyByMtRevAAMatrix(){
	int modNum=0;
	MODEL_FLOAT **qmatOffset = qmat[modNum];
	qmatOffset  [ 0 ][ 14 ] *= 23.18   ; qmatOffset [ 14 ][ 0 ] *= 23.18 ;
	qmatOffset  [ 0 ][ 11 ] *= 26.95   ; qmatOffset [ 11 ][ 0 ] *= 26.95 ;
	qmatOffset  [ 0 ][ 2 ] *= 17.67   ; qmatOffset [ 2 ][ 0 ] *= 17.67 ;
	qmatOffset  [ 0 ][ 1 ] *= 59.93   ; qmatOffset [ 1 ][ 0 ] *= 59.93 ;
	qmatOffset  [ 0 ][ 13 ] *= 1.9   ; qmatOffset [ 13 ][ 0 ] *= 1.9 ;
	qmatOffset  [ 0 ][ 3 ] *= 9.77   ; qmatOffset [ 3 ][ 0 ] *= 9.77 ;
	qmatOffset  [ 0 ][ 5 ] *= 120.71   ; qmatOffset [ 5 ][ 0 ] *= 120.71 ;
	qmatOffset  [ 0 ][ 6 ] *= 13.9   ; qmatOffset [ 6 ][ 0 ] *= 13.9 ;
	qmatOffset  [ 0 ][ 7 ] *= 96.49   ; qmatOffset [ 7 ][ 0 ] *= 96.49 ;
	qmatOffset  [ 0 ][ 9 ] *= 25.46   ; qmatOffset [ 9 ][ 0 ] *= 25.46 ;
	qmatOffset  [ 0 ][ 8 ] *= 8.36   ; qmatOffset [ 8 ][ 0 ] *= 8.36 ;
	qmatOffset  [ 0 ][ 10 ] *= 141.88   ; qmatOffset [ 10 ][ 0 ] *= 141.88 ;
	qmatOffset  [ 0 ][ 4 ] *= 6.37   ; qmatOffset [ 4 ][ 0 ] *= 6.37 ;
	qmatOffset  [ 0 ][ 12 ] *= 54.31   ; qmatOffset [ 12 ][ 0 ] *= 54.31 ;
	qmatOffset  [ 0 ][ 15 ] *= 387.86   ; qmatOffset [ 15 ][ 0 ] *= 387.86 ;
	qmatOffset  [ 0 ][ 16 ] *= 480.72   ; qmatOffset [ 16 ][ 0 ] *= 480.72 ;
	qmatOffset  [ 0 ][ 18 ] *= 1.9   ; qmatOffset [ 18 ][ 0 ] *= 1.9 ;
	qmatOffset  [ 0 ][ 19 ] *= 6.48   ; qmatOffset [ 19 ][ 0 ] *= 6.48 ;
	qmatOffset  [ 0 ][ 17 ] *= 195.06   ; qmatOffset [ 17 ][ 0 ] *= 195.06 ;
	qmatOffset  [ 14 ][ 11 ] *= 13.24   ; qmatOffset [ 11 ][ 14 ] *= 13.24 ;
	qmatOffset  [ 14 ][ 2 ] *= 1.9   ; qmatOffset [ 2 ][ 14 ] *= 1.9 ;
	qmatOffset  [ 14 ][ 1 ] *= 103.33   ; qmatOffset [ 1 ][ 14 ] *= 103.33 ;
	qmatOffset  [ 14 ][ 13 ] *= 220.99   ; qmatOffset [ 13 ][ 14 ] *= 220.99 ;
	qmatOffset  [ 14 ][ 3 ] *= 1.9   ; qmatOffset [ 3 ][ 14 ] *= 1.9 ;
	qmatOffset  [ 14 ][ 5 ] *= 23.03   ; qmatOffset [ 5 ][ 14 ] *= 23.03 ;
	qmatOffset  [ 14 ][ 6 ] *= 165.23   ; qmatOffset [ 6 ][ 14 ] *= 165.23 ;
	qmatOffset  [ 14 ][ 7 ] *= 1.9   ; qmatOffset [ 7 ][ 14 ] *= 1.9 ;
	qmatOffset  [ 14 ][ 9 ] *= 15.58   ; qmatOffset [ 9 ][ 14 ] *= 15.58 ;
	qmatOffset  [ 14 ][ 8 ] *= 141.4   ; qmatOffset [ 8 ][ 14 ] *= 141.4 ;
	qmatOffset  [ 14 ][ 10 ] *= 1.9   ; qmatOffset [ 10 ][ 14 ] *= 1.9 ;
	qmatOffset  [ 14 ][ 4 ] *= 4.69   ; qmatOffset [ 4 ][ 14 ] *= 4.69 ;
	qmatOffset  [ 14 ][ 12 ] *= 23.64   ; qmatOffset [ 12 ][ 14 ] *= 23.64 ;
	qmatOffset  [ 14 ][ 15 ] *= 6.04   ; qmatOffset [ 15 ][ 14 ] *= 6.04 ;
	qmatOffset  [ 14 ][ 16 ] *= 2.08   ; qmatOffset [ 16 ][ 14 ] *= 2.08 ;
	qmatOffset  [ 14 ][ 18 ] *= 21.95   ; qmatOffset [ 18 ][ 14 ] *= 21.95 ;
	qmatOffset  [ 14 ][ 19 ] *= 1.9   ; qmatOffset [ 19 ][ 14 ] *= 1.9 ;
	qmatOffset  [ 14 ][ 17 ] *= 7.64   ; qmatOffset [ 17 ][ 14 ] *= 7.64 ;
	qmatOffset  [ 11 ][ 2 ] *= 794.38   ; qmatOffset [ 2 ][ 11 ] *= 794.38 ;
	qmatOffset  [ 11 ][ 1 ] *= 58.94   ; qmatOffset [ 1 ][ 11 ] *= 58.94 ;
	qmatOffset  [ 11 ][ 13 ] *= 173.56   ; qmatOffset [ 13 ][ 11 ] *= 173.56 ;
	qmatOffset  [ 11 ][ 3 ] *= 63.05   ; qmatOffset [ 3 ][ 11 ] *= 63.05 ;
	qmatOffset  [ 11 ][ 5 ] *= 53.3   ; qmatOffset [ 5 ][ 11 ] *= 53.3 ;
	qmatOffset  [ 11 ][ 6 ] *= 496.13   ; qmatOffset [ 6 ][ 11 ] *= 496.13 ;
	qmatOffset  [ 11 ][ 7 ] *= 27.1   ; qmatOffset [ 7 ][ 11 ] *= 27.1 ;
	qmatOffset  [ 11 ][ 9 ] *= 15.16   ; qmatOffset [ 9 ][ 11 ] *= 15.16 ;
	qmatOffset  [ 11 ][ 8 ] *= 608.7   ; qmatOffset [ 8 ][ 11 ] *= 608.7 ;
	qmatOffset  [ 11 ][ 10 ] *= 65.41   ; qmatOffset [ 10 ][ 11 ] *= 65.41 ;
	qmatOffset  [ 11 ][ 4 ] *= 15.2   ; qmatOffset [ 4 ][ 11 ] *= 15.2 ;
	qmatOffset  [ 11 ][ 12 ] *= 73.31   ; qmatOffset [ 12 ][ 11 ] *= 73.31 ;
	qmatOffset  [ 11 ][ 15 ] *= 494.39   ; qmatOffset [ 15 ][ 11 ] *= 494.39 ;
	qmatOffset  [ 11 ][ 16 ] *= 238.46   ; qmatOffset [ 16 ][ 11 ] *= 238.46 ;
	qmatOffset  [ 11 ][ 18 ] *= 10.68   ; qmatOffset [ 18 ][ 11 ] *= 10.68 ;
	qmatOffset  [ 11 ][ 19 ] *= 191.36   ; qmatOffset [ 19 ][ 11 ] *= 191.36 ;
	qmatOffset  [ 11 ][ 17 ] *= 1.9   ; qmatOffset [ 17 ][ 11 ] *= 1.9 ;
	qmatOffset  [ 2 ][ 1 ] *= 1.9   ; qmatOffset [ 1 ][ 2 ] *= 1.9 ;
	qmatOffset  [ 2 ][ 13 ] *= 55.28   ; qmatOffset [ 13 ][ 2 ] *= 55.28 ;
	qmatOffset  [ 2 ][ 3 ] *= 583.55   ; qmatOffset [ 3 ][ 2 ] *= 583.55 ;
	qmatOffset  [ 2 ][ 5 ] *= 56.77   ; qmatOffset [ 5 ][ 2 ] *= 56.77 ;
	qmatOffset  [ 2 ][ 6 ] *= 113.99   ; qmatOffset [ 6 ][ 2 ] *= 113.99 ;
	qmatOffset  [ 2 ][ 7 ] *= 4.34   ; qmatOffset [ 7 ][ 2 ] *= 4.34 ;
	qmatOffset  [ 2 ][ 9 ] *= 1.9   ; qmatOffset [ 9 ][ 2 ] *= 1.9 ;
	qmatOffset  [ 2 ][ 8 ] *= 2.31   ; qmatOffset [ 8 ][ 2 ] *= 2.31 ;
	qmatOffset  [ 2 ][ 10 ] *= 1.9   ; qmatOffset [ 10 ][ 2 ] *= 1.9 ;
	qmatOffset  [ 2 ][ 4 ] *= 4.98   ; qmatOffset [ 4 ][ 2 ] *= 4.98 ;
	qmatOffset  [ 2 ][ 12 ] *= 13.43   ; qmatOffset [ 12 ][ 2 ] *= 13.43 ;
	qmatOffset  [ 2 ][ 15 ] *= 69.02   ; qmatOffset [ 15 ][ 2 ] *= 69.02 ;
	qmatOffset  [ 2 ][ 16 ] *= 28.01   ; qmatOffset [ 16 ][ 2 ] *= 28.01 ;
	qmatOffset  [ 2 ][ 18 ] *= 19.86   ; qmatOffset [ 18 ][ 2 ] *= 19.86 ;
	qmatOffset  [ 2 ][ 19 ] *= 21.21   ; qmatOffset [ 19 ][ 2 ] *= 21.21 ;
	qmatOffset  [ 2 ][ 17 ] *= 1.9   ; qmatOffset [ 17 ][ 2 ] *= 1.9 ;
	qmatOffset  [ 1 ][ 13 ] *= 75.24   ; qmatOffset [ 13 ][ 1 ] *= 75.24 ;
	qmatOffset  [ 1 ][ 3 ] *= 1.9   ; qmatOffset [ 3 ][ 1 ] *= 1.9 ;
	qmatOffset  [ 1 ][ 5 ] *= 30.71   ; qmatOffset [ 5 ][ 1 ] *= 30.71 ;
	qmatOffset  [ 1 ][ 6 ] *= 141.49   ; qmatOffset [ 6 ][ 1 ] *= 141.49 ;
	qmatOffset  [ 1 ][ 7 ] *= 62.73   ; qmatOffset [ 7 ][ 1 ] *= 62.73 ;
	qmatOffset  [ 1 ][ 9 ] *= 25.65   ; qmatOffset [ 9 ][ 1 ] *= 25.65 ;
	qmatOffset  [ 1 ][ 8 ] *= 1.9   ; qmatOffset [ 8 ][ 1 ] *= 1.9 ;
	qmatOffset  [ 1 ][ 10 ] *= 6.18   ; qmatOffset [ 10 ][ 1 ] *= 6.18 ;
	qmatOffset  [ 1 ][ 4 ] *= 70.8   ; qmatOffset [ 4 ][ 1 ] *= 70.8 ;
	qmatOffset  [ 1 ][ 12 ] *= 31.26   ; qmatOffset [ 12 ][ 1 ] *= 31.26 ;
	qmatOffset  [ 1 ][ 15 ] *= 277.05   ; qmatOffset [ 15 ][ 1 ] *= 277.05 ;
	qmatOffset  [ 1 ][ 16 ] *= 179.97   ; qmatOffset [ 16 ][ 1 ] *= 179.97 ;
	qmatOffset  [ 1 ][ 18 ] *= 33.6   ; qmatOffset [ 18 ][ 1 ] *= 33.6 ;
	qmatOffset  [ 1 ][ 19 ] *= 254.77   ; qmatOffset [ 19 ][ 1 ] *= 254.77 ;
	qmatOffset  [ 1 ][ 17 ] *= 1.9   ; qmatOffset [ 17 ][ 1 ] *= 1.9 ;
	qmatOffset  [ 13 ][ 3 ] *= 313.56   ; qmatOffset [ 3 ][ 13 ] *= 313.56 ;
	qmatOffset  [ 13 ][ 5 ] *= 6.75   ; qmatOffset [ 5 ][ 13 ] *= 6.75 ;
	qmatOffset  [ 13 ][ 6 ] *= 582.4   ; qmatOffset [ 6 ][ 13 ] *= 582.4 ;
	qmatOffset  [ 13 ][ 7 ] *= 8.34   ; qmatOffset [ 7 ][ 13 ] *= 8.34 ;
	qmatOffset  [ 13 ][ 9 ] *= 39.7   ; qmatOffset [ 9 ][ 13 ] *= 39.7 ;
	qmatOffset  [ 13 ][ 8 ] *= 465.58   ; qmatOffset [ 8 ][ 13 ] *= 465.58 ;
	qmatOffset  [ 13 ][ 10 ] *= 47.37   ; qmatOffset [ 10 ][ 13 ] *= 47.37 ;
	qmatOffset  [ 13 ][ 4 ] *= 19.11   ; qmatOffset [ 4 ][ 13 ] *= 19.11 ;
	qmatOffset  [ 13 ][ 12 ] *= 137.29   ; qmatOffset [ 12 ][ 13 ] *= 137.29 ;
	qmatOffset  [ 13 ][ 15 ] *= 54.11   ; qmatOffset [ 15 ][ 13 ] *= 54.11 ;
	qmatOffset  [ 13 ][ 16 ] *= 94.93   ; qmatOffset [ 16 ][ 13 ] *= 94.93 ;
	qmatOffset  [ 13 ][ 18 ] *= 1.9   ; qmatOffset [ 18 ][ 13 ] *= 1.9 ;
	qmatOffset  [ 13 ][ 19 ] *= 38.82   ; qmatOffset [ 19 ][ 13 ] *= 38.82 ;
	qmatOffset  [ 13 ][ 17 ] *= 19   ; qmatOffset [ 17 ][ 13 ] *= 19 ;
	qmatOffset  [ 3 ][ 5 ] *= 28.28   ; qmatOffset [ 5 ][ 3 ] *= 28.28 ;
	qmatOffset  [ 3 ][ 6 ] *= 49.12   ; qmatOffset [ 6 ][ 3 ] *= 49.12 ;
	qmatOffset  [ 3 ][ 7 ] *= 3.31   ; qmatOffset [ 7 ][ 3 ] *= 3.31 ;
	qmatOffset  [ 3 ][ 9 ] *= 1.9   ; qmatOffset [ 9 ][ 3 ] *= 1.9 ;
	qmatOffset  [ 3 ][ 8 ] *= 313.86   ; qmatOffset [ 8 ][ 3 ] *= 313.86 ;
	qmatOffset  [ 3 ][ 10 ] *= 1.9   ; qmatOffset [ 10 ][ 3 ] *= 1.9 ;
	qmatOffset  [ 3 ][ 4 ] *= 2.67   ; qmatOffset [ 4 ][ 3 ] *= 2.67 ;
	qmatOffset  [ 3 ][ 12 ] *= 12.83   ; qmatOffset [ 12 ][ 3 ] *= 12.83 ;
	qmatOffset  [ 3 ][ 15 ] *= 54.71   ; qmatOffset [ 15 ][ 3 ] *= 54.71 ;
	qmatOffset  [ 3 ][ 16 ] *= 14.82   ; qmatOffset [ 16 ][ 3 ] *= 14.82 ;
	qmatOffset  [ 3 ][ 18 ] *= 1.9   ; qmatOffset [ 18 ][ 3 ] *= 1.9 ;
	qmatOffset  [ 3 ][ 19 ] *= 13.12   ; qmatOffset [ 19 ][ 3 ] *= 13.12 ;
	qmatOffset  [ 3 ][ 17 ] *= 21.14   ; qmatOffset [ 17 ][ 3 ] *= 21.14 ;
	qmatOffset  [ 5 ][ 6 ] *= 1.9   ; qmatOffset [ 6 ][ 5 ] *= 1.9 ;
	qmatOffset  [ 5 ][ 7 ] *= 5.98   ; qmatOffset [ 7 ][ 5 ] *= 5.98 ;
	qmatOffset  [ 5 ][ 9 ] *= 2.41   ; qmatOffset [ 9 ][ 5 ] *= 2.41 ;
	qmatOffset  [ 5 ][ 8 ] *= 22.73   ; qmatOffset [ 8 ][ 5 ] *= 22.73 ;
	qmatOffset  [ 5 ][ 10 ] *= 1.9   ; qmatOffset [ 10 ][ 5 ] *= 1.9 ;
	qmatOffset  [ 5 ][ 4 ] *= 1.9   ; qmatOffset [ 4 ][ 5 ] *= 1.9 ;
	qmatOffset  [ 5 ][ 12 ] *= 1.9   ; qmatOffset [ 12 ][ 5 ] *= 1.9 ;
	qmatOffset  [ 5 ][ 15 ] *= 125.93   ; qmatOffset [ 15 ][ 5 ] *= 125.93 ;
	qmatOffset  [ 5 ][ 16 ] *= 11.17   ; qmatOffset [ 16 ][ 5 ] *= 11.17 ;
	qmatOffset  [ 5 ][ 18 ] *= 10.92   ; qmatOffset [ 18 ][ 5 ] *= 10.92 ;
	qmatOffset  [ 5 ][ 19 ] *= 3.21   ; qmatOffset [ 19 ][ 5 ] *= 3.21 ;
	qmatOffset  [ 5 ][ 17 ] *= 2.53   ; qmatOffset [ 17 ][ 5 ] *= 2.53 ;
	qmatOffset  [ 6 ][ 7 ] *= 12.26   ; qmatOffset [ 7 ][ 6 ] *= 12.26 ;
	qmatOffset  [ 6 ][ 9 ] *= 11.49   ; qmatOffset [ 9 ][ 6 ] *= 11.49 ;
	qmatOffset  [ 6 ][ 8 ] *= 127.67   ; qmatOffset [ 8 ][ 6 ] *= 127.67 ;
	qmatOffset  [ 6 ][ 10 ] *= 11.97   ; qmatOffset [ 10 ][ 6 ] *= 11.97 ;
	qmatOffset  [ 6 ][ 4 ] *= 48.16   ; qmatOffset [ 4 ][ 6 ] *= 48.16 ;
	qmatOffset  [ 6 ][ 12 ] *= 60.97   ; qmatOffset [ 12 ][ 6 ] *= 60.97 ;
	qmatOffset  [ 6 ][ 15 ] *= 77.46   ; qmatOffset [ 15 ][ 6 ] *= 77.46 ;
	qmatOffset  [ 6 ][ 16 ] *= 44.78   ; qmatOffset [ 16 ][ 6 ] *= 44.78 ;
	qmatOffset  [ 6 ][ 18 ] *= 7.08   ; qmatOffset [ 18 ][ 6 ] *= 7.08 ;
	qmatOffset  [ 6 ][ 19 ] *= 670.14   ; qmatOffset [ 19 ][ 6 ] *= 670.14 ;
	qmatOffset  [ 6 ][ 17 ] *= 1.9   ; qmatOffset [ 17 ][ 6 ] *= 1.9 ;
	qmatOffset  [ 7 ][ 9 ] *= 329.09   ; qmatOffset [ 9 ][ 7 ] *= 329.09 ;
	qmatOffset  [ 7 ][ 8 ] *= 19.57   ; qmatOffset [ 8 ][ 7 ] *= 19.57 ;
	qmatOffset  [ 7 ][ 10 ] *= 517.98   ; qmatOffset [ 10 ][ 7 ] *= 517.98 ;
	qmatOffset  [ 7 ][ 4 ] *= 84.67   ; qmatOffset [ 4 ][ 7 ] *= 84.67 ;
	qmatOffset  [ 7 ][ 12 ] *= 20.63   ; qmatOffset [ 12 ][ 7 ] *= 20.63 ;
	qmatOffset  [ 7 ][ 15 ] *= 47.7   ; qmatOffset [ 15 ][ 7 ] *= 47.7 ;
	qmatOffset  [ 7 ][ 16 ] *= 368.43   ; qmatOffset [ 16 ][ 7 ] *= 368.43 ;
	qmatOffset  [ 7 ][ 18 ] *= 1.9   ; qmatOffset [ 18 ][ 7 ] *= 1.9 ;
	qmatOffset  [ 7 ][ 19 ] *= 25.01   ; qmatOffset [ 19 ][ 7 ] *= 25.01 ;
	qmatOffset  [ 7 ][ 17 ] *= 1222.94   ; qmatOffset [ 17 ][ 7 ] *= 1222.94 ;
	qmatOffset  [ 9 ][ 8 ] *= 14.88   ; qmatOffset [ 8 ][ 9 ] *= 14.88 ;
	qmatOffset  [ 9 ][ 10 ] *= 537.53   ; qmatOffset [ 10 ][ 9 ] *= 537.53 ;
	qmatOffset  [ 9 ][ 4 ] *= 216.06   ; qmatOffset [ 4 ][ 9 ] *= 216.06 ;
	qmatOffset  [ 9 ][ 12 ] *= 40.1   ; qmatOffset [ 12 ][ 9 ] *= 40.1 ;
	qmatOffset  [ 9 ][ 15 ] *= 73.61   ; qmatOffset [ 15 ][ 9 ] *= 73.61 ;
	qmatOffset  [ 9 ][ 16 ] *= 126.4   ; qmatOffset [ 16 ][ 9 ] *= 126.4 ;
	qmatOffset  [ 9 ][ 18 ] *= 32.44   ; qmatOffset [ 18 ][ 9 ] *= 32.44 ;
	qmatOffset  [ 9 ][ 19 ] *= 44.15   ; qmatOffset [ 19 ][ 9 ] *= 44.15 ;
	qmatOffset  [ 9 ][ 17 ] *= 91.67   ; qmatOffset [ 17 ][ 9 ] *= 91.67 ;
	qmatOffset  [ 8 ][ 10 ] *= 91.37   ; qmatOffset [ 10 ][ 8 ] *= 91.37 ;
	qmatOffset  [ 8 ][ 4 ] *= 6.44   ; qmatOffset [ 4 ][ 8 ] *= 6.44 ;
	qmatOffset  [ 8 ][ 12 ] *= 50.1   ; qmatOffset [ 12 ][ 8 ] *= 50.1 ;
	qmatOffset  [ 8 ][ 15 ] *= 105.79   ; qmatOffset [ 15 ][ 8 ] *= 105.79 ;
	qmatOffset  [ 8 ][ 16 ] *= 136.33   ; qmatOffset [ 16 ][ 8 ] *= 136.33 ;
	qmatOffset  [ 8 ][ 18 ] *= 24   ; qmatOffset [ 18 ][ 8 ] *= 24 ;
	qmatOffset  [ 8 ][ 19 ] *= 51.17   ; qmatOffset [ 19 ][ 8 ] *= 51.17 ;
	qmatOffset  [ 8 ][ 17 ] *= 1.9   ; qmatOffset [ 17 ][ 8 ] *= 1.9 ;
	qmatOffset  [ 10 ][ 4 ] *= 90.82   ; qmatOffset [ 4 ][ 10 ] *= 90.82 ;
	qmatOffset  [ 10 ][ 12 ] *= 18.84   ; qmatOffset [ 12 ][ 10 ] *= 18.84 ;
	qmatOffset  [ 10 ][ 15 ] *= 111.16   ; qmatOffset [ 15 ][ 10 ] *= 111.16 ;
	qmatOffset  [ 10 ][ 16 ] *= 528.17   ; qmatOffset [ 16 ][ 10 ] *= 528.17 ;
	qmatOffset  [ 10 ][ 18 ] *= 21.71   ; qmatOffset [ 18 ][ 10 ] *= 21.71 ;
	qmatOffset  [ 10 ][ 19 ] *= 39.96   ; qmatOffset [ 19 ][ 10 ] *= 39.96 ;
	qmatOffset  [ 10 ][ 17 ] *= 387.54   ; qmatOffset [ 17 ][ 10 ] *= 387.54 ;
	qmatOffset  [ 4 ][ 12 ] *= 17.31   ; qmatOffset [ 12 ][ 4 ] *= 17.31 ;
	qmatOffset  [ 4 ][ 15 ] *= 64.29   ; qmatOffset [ 15 ][ 4 ] *= 64.29 ;
	qmatOffset  [ 4 ][ 16 ] *= 33.85   ; qmatOffset [ 16 ][ 4 ] *= 33.85 ;
	qmatOffset  [ 4 ][ 18 ] *= 7.84   ; qmatOffset [ 18 ][ 4 ] *= 7.84 ;
	qmatOffset  [ 4 ][ 19 ] *= 465.58   ; qmatOffset [ 19 ][ 4 ] *= 465.58 ;
	qmatOffset  [ 4 ][ 17 ] *= 6.35   ; qmatOffset [ 17 ][ 4 ] *= 6.35 ;
	qmatOffset  [ 12 ][ 15 ] *= 169.9   ; qmatOffset [ 15 ][ 12 ] *= 169.9 ;
	qmatOffset  [ 12 ][ 16 ] *= 128.22   ; qmatOffset [ 16 ][ 12 ] *= 128.22 ;
	qmatOffset  [ 12 ][ 18 ] *= 4.21   ; qmatOffset [ 18 ][ 12 ] *= 4.21 ;
	qmatOffset  [ 12 ][ 19 ] *= 16.21   ; qmatOffset [ 19 ][ 12 ] *= 16.21 ;
	qmatOffset  [ 12 ][ 17 ] *= 8.23   ; qmatOffset [ 17 ][ 12 ] *= 8.23 ;
	qmatOffset  [ 15 ][ 16 ] *= 597.21   ; qmatOffset [ 16 ][ 15 ] *= 597.21 ;
	qmatOffset  [ 15 ][ 18 ] *= 38.58   ; qmatOffset [ 18 ][ 15 ] *= 38.58 ;
	qmatOffset  [ 15 ][ 19 ] *= 64.92   ; qmatOffset [ 19 ][ 15 ] *= 64.92 ;
	qmatOffset  [ 15 ][ 17 ] *= 1.9   ; qmatOffset [ 17 ][ 15 ] *= 1.9 ;
	qmatOffset  [ 16 ][ 18 ] *= 9.99   ; qmatOffset [ 18 ][ 16 ] *= 9.99 ;
	qmatOffset  [ 16 ][ 19 ] *= 38.73   ; qmatOffset [ 19 ][ 16 ] *= 38.73 ;
	qmatOffset  [ 16 ][ 17 ] *= 204.54   ; qmatOffset [ 17 ][ 16 ] *= 204.54 ;
	qmatOffset  [ 18 ][ 19 ] *= 26.25   ; qmatOffset [ 19 ][ 18 ] *= 26.25 ;
	qmatOffset  [ 18 ][ 17 ] *= 5.37   ; qmatOffset [ 17 ][ 18 ] *= 5.37 ;
	qmatOffset  [ 19 ][ 17 ] *= 1.9   ; qmatOffset [ 17 ][ 19 ] *= 1.9 ;
	}

void Model::MultiplyByDayhoffAAMatrix(){
	int modNum=0;
	MODEL_FLOAT **qmatOffset = qmat[modNum];

	qmatOffset[0][1] *= 0.036; qmatOffset[0][2] *= 0.12; qmatOffset[0][3] *= 0.198; qmatOffset[0][4] *= 0.018; qmatOffset[0][5] *= 0.24; qmatOffset[0][6] *= 0.023;
	qmatOffset[0][7] *= 0.065; qmatOffset[0][8] *= 0.026; qmatOffset[0][9] *= 0.041; qmatOffset[0][10] *= 0.072; qmatOffset[0][11] *= 0.098; qmatOffset[0][12] *= 0.25;
	qmatOffset[0][13] *= 0.089; qmatOffset[0][14] *= 0.027; qmatOffset[0][15] *= 0.409; qmatOffset[0][16] *= 0.371; qmatOffset[0][17] *= 0.208; qmatOffset[0][19] *= 0.024;
	qmatOffset[1][0] *= 0.036; qmatOffset[1][5] *= 0.011; qmatOffset[1][6] *= 0.028; qmatOffset[1][7] *= 0.044; qmatOffset[1][12] *= 0.019; qmatOffset[1][14] *= 0.023;
	qmatOffset[1][15] *= 0.161; qmatOffset[1][16] *= 0.016; qmatOffset[1][17] *= 0.049; qmatOffset[1][19] *= 0.096; qmatOffset[2][0] *= 0.12; qmatOffset[2][3] *= 1.153; 
	qmatOffset[2][5] *= 0.125; qmatOffset[2][6] *= 0.086; qmatOffset[2][7] *= 0.024; qmatOffset[2][8] *= 0.071; qmatOffset[2][11] *= 0.905; qmatOffset[2][12] *= 0.013; 
	qmatOffset[2][13] *= 0.134; qmatOffset[2][15] *= 0.095; qmatOffset[2][16] *= 0.066; qmatOffset[2][17] *= 0.018; qmatOffset[3][0] *= 0.198; qmatOffset[3][2] *= 1.153; 
	qmatOffset[3][5] *= 0.081; qmatOffset[3][6] *= 0.043; qmatOffset[3][7] *= 0.061; qmatOffset[3][8] *= 0.083; qmatOffset[3][9] *= 0.011; qmatOffset[3][10] *= 0.03; 
	qmatOffset[3][11] *= 0.148; qmatOffset[3][12] *= 0.051; qmatOffset[3][13] *= 0.716; qmatOffset[3][14] *= 0.001; qmatOffset[3][15] *= 0.079; qmatOffset[3][16] *= 0.034; 
	qmatOffset[3][17] *= 0.037; qmatOffset[3][19] *= 0.022; qmatOffset[4][0] *= 0.018; qmatOffset[4][5] *= 0.015; qmatOffset[4][6] *= 0.048; qmatOffset[4][7] *= 0.196; 
	qmatOffset[4][9] *= 0.157; qmatOffset[4][10] *= 0.092; qmatOffset[4][11] *= 0.014; qmatOffset[4][12] *= 0.011; qmatOffset[4][14] *= 0.014; qmatOffset[4][15] *= 0.046; 
	qmatOffset[4][16] *= 0.013; qmatOffset[4][17] *= 0.012; qmatOffset[4][18] *= 0.076; qmatOffset[4][19] *= 0.698; qmatOffset[5][0] *= 0.24; qmatOffset[5][1] *= 0.011; 
	qmatOffset[5][2] *= 0.125; qmatOffset[5][3] *= 0.081; qmatOffset[5][4] *= 0.015; qmatOffset[5][6] *= 0.01; qmatOffset[5][8] *= 0.027; qmatOffset[5][9] *= 0.007; 
	qmatOffset[5][10] *= 0.017; qmatOffset[5][11] *= 0.139; qmatOffset[5][12] *= 0.034; qmatOffset[5][13] *= 0.028; qmatOffset[5][14] *= 0.009; qmatOffset[5][15] *= 0.234; 
	qmatOffset[5][16] *= 0.03; qmatOffset[5][17] *= 0.054; qmatOffset[6][0] *= 0.023; qmatOffset[6][1] *= 0.028; qmatOffset[6][2] *= 0.086; qmatOffset[6][3] *= 0.043; 
	qmatOffset[6][4] *= 0.048; qmatOffset[6][5] *= 0.01; qmatOffset[6][7] *= 0.007; qmatOffset[6][8] *= 0.026; qmatOffset[6][9] *= 0.044; qmatOffset[6][11] *= 0.535; 
	qmatOffset[6][12] *= 0.094; qmatOffset[6][13] *= 0.606; qmatOffset[6][14] *= 0.24; qmatOffset[6][15] *= 0.035; qmatOffset[6][16] *= 0.022; qmatOffset[6][17] *= 0.044; 
	qmatOffset[6][18] *= 0.027; qmatOffset[6][19] *= 0.127; qmatOffset[7][0] *= 0.065; qmatOffset[7][1] *= 0.044; qmatOffset[7][2] *= 0.024; qmatOffset[7][3] *= 0.061; 
	qmatOffset[7][4] *= 0.196; qmatOffset[7][6] *= 0.007; qmatOffset[7][8] *= 0.046; qmatOffset[7][9] *= 0.257; qmatOffset[7][10] *= 0.336; qmatOffset[7][11] *= 0.077; 
	qmatOffset[7][12] *= 0.012; qmatOffset[7][13] *= 0.018; qmatOffset[7][14] *= 0.064; qmatOffset[7][15] *= 0.024; qmatOffset[7][16] *= 0.192; qmatOffset[7][17] *= 0.889; 
	qmatOffset[7][19] *= 0.037; qmatOffset[8][0] *= 0.026; qmatOffset[8][2] *= 0.071; qmatOffset[8][3] *= 0.083; qmatOffset[8][5] *= 0.027; qmatOffset[8][6] *= 0.026; 
	qmatOffset[8][7] *= 0.046; qmatOffset[8][9] *= 0.018; qmatOffset[8][10] *= 0.243; qmatOffset[8][11] *= 0.318; qmatOffset[8][12] *= 0.033; qmatOffset[8][13] *= 0.153; 
	qmatOffset[8][14] *= 0.464; qmatOffset[8][15] *= 0.096; qmatOffset[8][16] *= 0.136; qmatOffset[8][17] *= 0.01; qmatOffset[8][19] *= 0.013; qmatOffset[9][0] *= 0.041; 
	qmatOffset[9][3] *= 0.011; qmatOffset[9][4] *= 0.157; qmatOffset[9][5] *= 0.007; qmatOffset[9][6] *= 0.044; qmatOffset[9][7] *= 0.257; qmatOffset[9][8] *= 0.018; 
	qmatOffset[9][10] *= 0.527; qmatOffset[9][11] *= 0.034; qmatOffset[9][12] *= 0.032; qmatOffset[9][13] *= 0.073; qmatOffset[9][14] *= 0.015; qmatOffset[9][15] *= 0.017; 
	qmatOffset[9][16] *= 0.033; qmatOffset[9][17] *= 0.175; qmatOffset[9][18] *= 0.046; qmatOffset[9][19] *= 0.028; qmatOffset[10][0] *= 0.072; qmatOffset[10][3] *= 0.03; 
	qmatOffset[10][4] *= 0.092; qmatOffset[10][5] *= 0.017; qmatOffset[10][7] *= 0.336; qmatOffset[10][8] *= 0.243; qmatOffset[10][9] *= 0.527; qmatOffset[10][11] *= 0.001; 
	qmatOffset[10][12] *= 0.017; qmatOffset[10][13] *= 0.114; qmatOffset[10][14] *= 0.09; qmatOffset[10][15] *= 0.062; qmatOffset[10][16] *= 0.104; qmatOffset[10][17] *= 0.258;
	qmatOffset[11][0] *= 0.098; qmatOffset[11][2] *= 0.905; qmatOffset[11][3] *= 0.148; qmatOffset[11][4] *= 0.014; qmatOffset[11][5] *= 0.139; qmatOffset[11][6] *= 0.535; 
	qmatOffset[11][7] *= 0.077; qmatOffset[11][8] *= 0.318; qmatOffset[11][9] *= 0.034; qmatOffset[11][10] *= 0.001; qmatOffset[11][12] *= 0.042; qmatOffset[11][13] *= 0.103;
	qmatOffset[11][14] *= 0.032; qmatOffset[11][15] *= 0.495; qmatOffset[11][16] *= 0.229; qmatOffset[11][17] *= 0.015; qmatOffset[11][18] *= 0.023; qmatOffset[11][19] *= 0.095;
	qmatOffset[12][0] *= 0.25; qmatOffset[12][1] *= 0.019; qmatOffset[12][2] *= 0.013; qmatOffset[12][3] *= 0.051; qmatOffset[12][4] *= 0.011; qmatOffset[12][5] *= 0.034;
	qmatOffset[12][6] *= 0.094; qmatOffset[12][7] *= 0.012; qmatOffset[12][8] *= 0.033; qmatOffset[12][9] *= 0.032; qmatOffset[12][10] *= 0.017; qmatOffset[12][11] *= 0.042;
	qmatOffset[12][13] *= 0.153; qmatOffset[12][14] *= 0.103; qmatOffset[12][15] *= 0.245; qmatOffset[12][16] *= 0.078; qmatOffset[12][17] *= 0.048; qmatOffset[13][0] *= 0.089;
	qmatOffset[13][2] *= 0.134; qmatOffset[13][3] *= 0.716; qmatOffset[13][5] *= 0.028; qmatOffset[13][6] *= 0.606; qmatOffset[13][7] *= 0.018; qmatOffset[13][8] *= 0.153;
	qmatOffset[13][9] *= 0.073; qmatOffset[13][10] *= 0.114; qmatOffset[13][11] *= 0.103; qmatOffset[13][12] *= 0.153; qmatOffset[13][14] *= 0.246; qmatOffset[13][15] *= 0.056;
	qmatOffset[13][16] *= 0.053; qmatOffset[13][17] *= 0.035; qmatOffset[14][0] *= 0.027; qmatOffset[14][1] *= 0.023; qmatOffset[14][3] *= 0.001; qmatOffset[14][4] *= 0.014;
	qmatOffset[14][5] *= 0.009; qmatOffset[14][6] *= 0.24; qmatOffset[14][7] *= 0.064; qmatOffset[14][8] *= 0.464; qmatOffset[14][9] *= 0.015; qmatOffset[14][10] *= 0.09;
	qmatOffset[14][11] *= 0.032; qmatOffset[14][12] *= 0.103; qmatOffset[14][13] *= 0.246; qmatOffset[14][15] *= 0.154; qmatOffset[14][16] *= 0.026; qmatOffset[14][17] *= 0.024;
	qmatOffset[14][18] *= 0.201; qmatOffset[14][19] *= 0.008; qmatOffset[15][0] *= 0.409; qmatOffset[15][1] *= 0.161; qmatOffset[15][2] *= 0.095; qmatOffset[15][3] *= 0.079;
	qmatOffset[15][4] *= 0.046; qmatOffset[15][5] *= 0.234; qmatOffset[15][6] *= 0.035; qmatOffset[15][7] *= 0.024; qmatOffset[15][8] *= 0.096; qmatOffset[15][9] *= 0.017;
	qmatOffset[15][10] *= 0.062; qmatOffset[15][11] *= 0.495; qmatOffset[15][12] *= 0.245; qmatOffset[15][13] *= 0.056; qmatOffset[15][14] *= 0.154; qmatOffset[15][16] *= 0.55;
	qmatOffset[15][17] *= 0.03; qmatOffset[15][18] *= 0.075; qmatOffset[15][19] *= 0.034; qmatOffset[16][0] *= 0.371; qmatOffset[16][1] *= 0.016; qmatOffset[16][2] *= 0.066;
	qmatOffset[16][3] *= 0.034; qmatOffset[16][4] *= 0.013; qmatOffset[16][5] *= 0.03; qmatOffset[16][6] *= 0.022; qmatOffset[16][7] *= 0.192; qmatOffset[16][8] *= 0.136;
	qmatOffset[16][9] *= 0.033; qmatOffset[16][10] *= 0.104; qmatOffset[16][11] *= 0.229; qmatOffset[16][12] *= 0.078; qmatOffset[16][13] *= 0.053; qmatOffset[16][14] *= 0.026;
	qmatOffset[16][15] *= 0.55; qmatOffset[16][17] *= 0.157; qmatOffset[16][19] *= 0.042; qmatOffset[17][0] *= 0.208; qmatOffset[17][1] *= 0.049; qmatOffset[17][2] *= 0.018;
	qmatOffset[17][3] *= 0.037; qmatOffset[17][4] *= 0.012; qmatOffset[17][5] *= 0.054; qmatOffset[17][6] *= 0.044; qmatOffset[17][7] *= 0.889; qmatOffset[17][8] *= 0.01;
	qmatOffset[17][9] *= 0.175; qmatOffset[17][10] *= 0.258; qmatOffset[17][11] *= 0.015; qmatOffset[17][12] *= 0.048; qmatOffset[17][13] *= 0.035; qmatOffset[17][14] *= 0.024;
	qmatOffset[17][15] *= 0.03; qmatOffset[17][16] *= 0.157; qmatOffset[17][19] *= 0.028; qmatOffset[18][4] *= 0.076; qmatOffset[18][6] *= 0.027; qmatOffset[18][9] *= 0.046;
	qmatOffset[18][11] *= 0.023; qmatOffset[18][14] *= 0.201; qmatOffset[18][15] *= 0.075; qmatOffset[18][19] *= 0.061; qmatOffset[19][0] *= 0.024; qmatOffset[19][1] *= 0.096;
	qmatOffset[19][3] *= 0.022; qmatOffset[19][4] *= 0.698; qmatOffset[19][6] *= 0.127; qmatOffset[19][7] *= 0.037; qmatOffset[19][8] *= 0.013; qmatOffset[19][9] *= 0.028;
	qmatOffset[19][11] *= 0.095; qmatOffset[19][14] *= 0.008; qmatOffset[19][15] *= 0.034; qmatOffset[19][16] *= 0.042; qmatOffset[19][17] *= 0.028; 
	qmatOffset[19][18] *= 0.061;
	//here are the zero entries
	qmatOffset[0][18]=qmatOffset[18][0]=0.0;
	qmatOffset[2][14]=qmatOffset[14][2]=0.0;
	qmatOffset[1][11]=qmatOffset[11][1]=0.0;
	qmatOffset[1][2]=qmatOffset[2][1]=0.0;
	qmatOffset[9][2]=qmatOffset[2][9]=0.0;
	qmatOffset[10][2]=qmatOffset[2][10]=0.0;
	qmatOffset[4][2]=qmatOffset[2][4]=0.0;
	qmatOffset[18][2]=qmatOffset[2][18]=0.0;
	qmatOffset[19][2]=qmatOffset[2][19]=0.0;
	qmatOffset[13][1]=qmatOffset[1][13]=0.0;
	qmatOffset[3][1]=qmatOffset[1][3]=0.0;
	qmatOffset[9][1]=qmatOffset[1][9]=0.0;
	qmatOffset[8][1]=qmatOffset[1][8]=0.0;
	qmatOffset[10][1]=qmatOffset[1][10]=0.0;
	qmatOffset[4][1]=qmatOffset[1][4]=0.0;
	qmatOffset[18][1]=qmatOffset[1][18]=0.0;
	qmatOffset[4][13]=qmatOffset[13][4]=0.0;
	qmatOffset[18][13]=qmatOffset[13][18]=0.0;
	qmatOffset[19][13]=qmatOffset[13][19]=0.0;
	qmatOffset[4][3]=qmatOffset[3][4]=0.0;
	qmatOffset[18][3]=qmatOffset[3][18]=0.0;
	qmatOffset[7][5]=qmatOffset[5][7]=0.0;
	qmatOffset[18][5]=qmatOffset[5][18]=0.0;
	qmatOffset[19][5]=qmatOffset[5][19]=0.0;
	qmatOffset[10][6]=qmatOffset[6][10]=0.0;
	qmatOffset[18][7]=qmatOffset[7][18]=0.0;
	qmatOffset[4][8]=qmatOffset[8][4]=0.0;
	qmatOffset[18][8]=qmatOffset[8][18]=0.0;
	qmatOffset[18][10]=qmatOffset[10][18]=0.0;
	qmatOffset[19][10]=qmatOffset[10][19]=0.0;
	qmatOffset[18][12]=qmatOffset[12][18]=0.0;
	qmatOffset[19][12]=qmatOffset[12][19]=0.0;
	qmatOffset[18][16]=qmatOffset[16][18]=0.0;
	qmatOffset[17][18]=qmatOffset[18][17]=0.0;
	}

void Model::MultiplyByWAGAAMatrix(){
	int modNum=0;
	MODEL_FLOAT **qmatOffset = qmat[modNum];

	qmatOffset[14][0] *= 1.75252;
	qmatOffset[0][14] *= 1.75252;
	qmatOffset[11][0] *= 1.61995;
	qmatOffset[0][11] *= 1.61995;
	qmatOffset[11][14] *= 2.0187;
	qmatOffset[14][11] *= 2.0187;
	qmatOffset[2][0] *= 2.34804;
	qmatOffset[0][2] *= 2.34804;
	qmatOffset[2][14] *= 0.468033;
	qmatOffset[14][2] *= 0.468033;
	qmatOffset[2][11] *= 17.251;
	qmatOffset[11][2] *= 17.251;
	qmatOffset[1][0] *= 3.26324;
	qmatOffset[0][1] *= 3.26324;
	qmatOffset[1][14] *= 1.67824;
	qmatOffset[14][1] *= 1.67824;
	qmatOffset[1][11] *= 0.842805;
	qmatOffset[11][1] *= 0.842805;
	qmatOffset[1][2] *= 0.0962568;
	qmatOffset[2][1] *= 0.0962568;
	qmatOffset[13][0] *= 2.88691;
	qmatOffset[0][13] *= 2.88691;
	qmatOffset[13][14] *= 9.64477;
	qmatOffset[14][13] *= 9.64477;
	qmatOffset[13][11] *= 4.90465;
	qmatOffset[11][13] *= 4.90465;
	qmatOffset[13][2] *= 1.95972;
	qmatOffset[2][13] *= 1.95972;
	qmatOffset[13][1] *= 0.313977;
	qmatOffset[1][13] *= 0.313977;
	qmatOffset[3][0] *= 5.02923;
	qmatOffset[0][3] *= 5.02923;
	qmatOffset[3][14] *= 1.39535;
	qmatOffset[14][3] *= 1.39535;
	qmatOffset[3][11] *= 3.00956;
	qmatOffset[11][3] *= 3.00956;
	qmatOffset[3][2] *= 19.6173;
	qmatOffset[2][3] *= 19.6173;
	qmatOffset[3][1] *= 0.0678423;
	qmatOffset[1][3] *= 0.0678423;
	qmatOffset[3][13] *= 17.3783;
	qmatOffset[13][3] *= 17.3783;
	qmatOffset[5][0] *= 4.50138;
	qmatOffset[0][5] *= 4.50138;
	qmatOffset[5][14] *= 1.85767;
	qmatOffset[14][5] *= 1.85767;
	qmatOffset[5][11] *= 3.57627;
	qmatOffset[11][5] *= 3.57627;
	qmatOffset[5][2] *= 2.75024;
	qmatOffset[2][5] *= 2.75024;
	qmatOffset[5][1] *= 0.974403;
	qmatOffset[1][5] *= 0.974403;
	qmatOffset[5][13] *= 1.04868;
	qmatOffset[13][5] *= 1.04868;
	qmatOffset[5][3] *= 1.80382;
	qmatOffset[3][5] *= 1.80382;
	qmatOffset[6][0] *= 1.00707;
	qmatOffset[0][6] *= 1.00707;
	qmatOffset[6][14] *= 6.79042;
	qmatOffset[14][6] *= 6.79042;
	qmatOffset[6][11] *= 12.5704;
	qmatOffset[11][6] *= 12.5704;
	qmatOffset[6][2] *= 2.95706;
	qmatOffset[2][6] *= 2.95706;
	qmatOffset[6][1] *= 0.791065;
	qmatOffset[1][6] *= 0.791065;
	qmatOffset[6][13] *= 13.6438;
	qmatOffset[13][6] *= 13.6438;
	qmatOffset[6][3] *= 1.81116;
	qmatOffset[3][6] *= 1.81116;
	qmatOffset[6][5] *= 0.792457;
	qmatOffset[5][6] *= 0.792457;
	qmatOffset[7][0] *= 0.614288;
	qmatOffset[0][7] *= 0.614288;
	qmatOffset[7][14] *= 0.594093;
	qmatOffset[14][7] *= 0.594093;
	qmatOffset[7][11] *= 1.76099;
	qmatOffset[11][7] *= 1.76099;
	qmatOffset[7][2] *= 0.125304;
	qmatOffset[2][7] *= 0.125304;
	qmatOffset[7][1] *= 0.540574;
	qmatOffset[1][7] *= 0.540574;
	qmatOffset[7][13] *= 0.361952;
	qmatOffset[13][7] *= 0.361952;
	qmatOffset[7][3] *= 0.404776;
	qmatOffset[3][7] *= 0.404776;
	qmatOffset[7][5] *= 0.0967499;
	qmatOffset[5][7] *= 0.0967499;
	qmatOffset[7][6] *= 0.439075;
	qmatOffset[6][7] *= 0.439075;
	qmatOffset[9][0] *= 1.26431;
	qmatOffset[0][9] *= 1.26431;
	qmatOffset[9][14] *= 1.58126;
	qmatOffset[14][9] *= 1.58126;
	qmatOffset[9][11] *= 0.417907;
	qmatOffset[11][9] *= 0.417907;
	qmatOffset[9][2] *= 0.269452;
	qmatOffset[2][9] *= 0.269452;
	qmatOffset[9][1] *= 1.22101;
	qmatOffset[1][9] *= 1.22101;
	qmatOffset[9][13] *= 2.76265;
	qmatOffset[13][9] *= 2.76265;
	qmatOffset[9][3] *= 0.490144;
	qmatOffset[3][9] *= 0.490144;
	qmatOffset[9][5] *= 0.194782;
	qmatOffset[5][9] *= 0.194782;
	qmatOffset[9][6] *= 1.58695;
	qmatOffset[6][9] *= 1.58695;
	qmatOffset[9][7] *= 10.0752;
	qmatOffset[7][9] *= 10.0752;
	qmatOffset[8][0] *= 2.8795;
	qmatOffset[0][8] *= 2.8795;
	qmatOffset[8][14] *= 17.0032;
	qmatOffset[14][8] *= 17.0032;
	qmatOffset[8][11] *= 9.57014;
	qmatOffset[11][8] *= 9.57014;
	qmatOffset[8][2] *= 1.52466;
	qmatOffset[2][8] *= 1.52466;
	qmatOffset[8][1] *= 0.23523;
	qmatOffset[1][8] *= 0.23523;
	qmatOffset[8][13] *= 12.3754;
	qmatOffset[13][8] *= 12.3754;
	qmatOffset[8][3] *= 8.21158;
	qmatOffset[3][8] *= 8.21158;
	qmatOffset[8][5] *= 1.18692;
	qmatOffset[5][8] *= 1.18692;
	qmatOffset[8][6] *= 2.82919;
	qmatOffset[6][8] *= 2.82919;
	qmatOffset[8][7] *= 1.02892;
	qmatOffset[7][8] *= 1.02892;
	qmatOffset[8][9] *= 0.818336;
	qmatOffset[9][8] *= 0.818336;
	qmatOffset[10][0] *= 2.83893;
	qmatOffset[0][10] *= 2.83893;
	qmatOffset[10][14] *= 2.17063;
	qmatOffset[14][10] *= 2.17063;
	qmatOffset[10][11] *= 0.629813;
	qmatOffset[11][10] *= 0.629813;
	qmatOffset[10][2] *= 0.32966;
	qmatOffset[2][10] *= 0.32966;
	qmatOffset[10][1] *= 1.24069;
	qmatOffset[1][10] *= 1.24069;
	qmatOffset[10][13] *= 4.9098;
	qmatOffset[13][10] *= 4.9098;
	qmatOffset[10][3] *= 1.00125;
	qmatOffset[3][10] *= 1.00125;
	qmatOffset[10][5] *= 0.553173;
	qmatOffset[5][10] *= 0.553173;
	qmatOffset[10][6] *= 1.28409;
	qmatOffset[6][10] *= 1.28409;
	qmatOffset[10][7] *= 13.5273;
	qmatOffset[7][10] *= 13.5273;
	qmatOffset[10][9] *= 15.4228;
	qmatOffset[9][10] *= 15.4228;
	qmatOffset[10][8] *= 2.9685;
	qmatOffset[8][10] *= 2.9685;
	qmatOffset[4][0] *= 0.668808;
	qmatOffset[0][4] *= 0.668808;
	qmatOffset[4][14] *= 0.326346;
	qmatOffset[14][4] *= 0.326346;
	qmatOffset[4][11] *= 0.305538;
	qmatOffset[11][4] *= 0.305538;
	qmatOffset[4][2] *= 0.148478;
	qmatOffset[2][4] *= 0.148478;
	qmatOffset[4][1] *= 1.26464;
	qmatOffset[1][4] *= 1.26464;
	qmatOffset[4][13] *= 0.317481;
	qmatOffset[13][4] *= 0.317481;
	qmatOffset[4][3] *= 0.257789;
	qmatOffset[3][4] *= 0.257789;
	qmatOffset[4][5] *= 0.158647;
	qmatOffset[5][4] *= 0.158647;
	qmatOffset[4][6] *= 2.15858;
	qmatOffset[6][4] *= 2.15858;
	qmatOffset[4][7] *= 3.36628;
	qmatOffset[7][4] *= 3.36628;
	qmatOffset[4][9] *= 6.72059;
	qmatOffset[9][4] *= 6.72059;
	qmatOffset[4][8] *= 0.282261;
	qmatOffset[8][4] *= 0.282261;
	qmatOffset[4][10] *= 3.78302;
	qmatOffset[10][4] *= 3.78302;
	qmatOffset[12][0] *= 4.57074;
	qmatOffset[0][12] *= 4.57074;
	qmatOffset[12][14] *= 2.15896;
	qmatOffset[14][12] *= 2.15896;
	qmatOffset[12][11] *= 0.619836;
	qmatOffset[11][12] *= 0.619836;
	qmatOffset[12][2] *= 1.34714;
	qmatOffset[2][12] *= 1.34714;
	qmatOffset[12][1] *= 0.347612;
	qmatOffset[1][12] *= 0.347612;
	qmatOffset[12][13] *= 2.96563;
	qmatOffset[13][12] *= 2.96563;
	qmatOffset[12][3] *= 2.16806;
	qmatOffset[3][12] *= 2.16806;
	qmatOffset[12][5] *= 0.773901;
	qmatOffset[5][12] *= 0.773901;
	qmatOffset[12][6] *= 2.21205;
	qmatOffset[6][12] *= 2.21205;
	qmatOffset[12][7] *= 0.317506;
	qmatOffset[7][12] *= 0.317506;
	qmatOffset[12][9] *= 1.32127;
	qmatOffset[9][12] *= 1.32127;
	qmatOffset[12][8] *= 1.76944;
	qmatOffset[8][12] *= 1.76944;
	qmatOffset[12][10] *= 0.544368;
	qmatOffset[10][12] *= 0.544368;
	qmatOffset[12][4] *= 0.51296;
	qmatOffset[4][12] *= 0.51296;
	qmatOffset[15][0] *= 10.7101;
	qmatOffset[0][15] *= 10.7101;
	qmatOffset[15][14] *= 3.88965;
	qmatOffset[14][15] *= 3.88965;
	qmatOffset[15][11] *= 12.6274;
	qmatOffset[11][15] *= 12.6274;
	qmatOffset[15][2] *= 3.40533;
	qmatOffset[2][15] *= 3.40533;
	qmatOffset[15][1] *= 4.4726;
	qmatOffset[1][15] *= 4.4726;
	qmatOffset[15][13] *= 3.26906;
	qmatOffset[13][15] *= 3.26906;
	qmatOffset[15][3] *= 2.23982;
	qmatOffset[3][15] *= 2.23982;
	qmatOffset[15][5] *= 4.2634;
	qmatOffset[5][15] *= 4.2634;
	qmatOffset[15][6] *= 2.35176;
	qmatOffset[6][15] *= 2.35176;
	qmatOffset[15][7] *= 1.01497;
	qmatOffset[7][15] *= 1.01497;
	qmatOffset[15][9] *= 1.09535;
	qmatOffset[9][15] *= 1.09535;
	qmatOffset[15][8] *= 3.07289;
	qmatOffset[8][15] *= 3.07289;
	qmatOffset[15][10] *= 1.5693;
	qmatOffset[10][15] *= 1.5693;
	qmatOffset[15][4] *= 1.7346;
	qmatOffset[4][15] *= 1.7346;
	qmatOffset[15][12] *= 5.12592;
	qmatOffset[12][15] *= 5.12592;
	qmatOffset[16][0] *= 6.73946;
	qmatOffset[0][16] *= 6.73946;
	qmatOffset[16][14] *= 1.76155;
	qmatOffset[14][16] *= 1.76155;
	qmatOffset[16][11] *= 6.45016;
	qmatOffset[11][16] *= 6.45016;
	qmatOffset[16][2] *= 1.19107;
	qmatOffset[2][16] *= 1.19107;
	qmatOffset[16][1] *= 1.62992;
	qmatOffset[1][16] *= 1.62992;
	qmatOffset[16][13] *= 2.72592;
	qmatOffset[13][16] *= 2.72592;
	qmatOffset[16][3] *= 2.61419;
	qmatOffset[3][16] *= 2.61419;
	qmatOffset[16][5] *= 0.717545;
	qmatOffset[5][16] *= 0.717545;
	qmatOffset[16][6] *= 1.50385;
	qmatOffset[6][16] *= 1.50385;
	qmatOffset[16][7] *= 4.63305;
	qmatOffset[7][16] *= 4.63305;
	qmatOffset[16][9] *= 1.03778;
	qmatOffset[9][16] *= 1.03778;
	qmatOffset[16][8] *= 4.40689;
	qmatOffset[8][16] *= 4.40689;
	qmatOffset[16][10] *= 4.81721;
	qmatOffset[10][16] *= 4.81721;
	qmatOffset[16][4] *= 0.546192;
	qmatOffset[4][16] *= 0.546192;
	qmatOffset[16][12] *= 2.52719;
	qmatOffset[12][16] *= 2.52719;
	qmatOffset[16][15] *= 13.9104;
	qmatOffset[15][16] *= 13.9104;
	qmatOffset[18][0] *= 0.35946;
	qmatOffset[0][18] *= 0.35946;
	qmatOffset[18][14] *= 3.69815;
	qmatOffset[14][18] *= 3.69815;
	qmatOffset[18][11] *= 0.228503;
	qmatOffset[11][18] *= 0.228503;
	qmatOffset[18][2] *= 0.412312;
	qmatOffset[2][18] *= 0.412312;
	qmatOffset[18][1] *= 2.27837;
	qmatOffset[1][18] *= 2.27837;
	qmatOffset[18][13] *= 0.685467;
	qmatOffset[13][18] *= 0.685467;
	qmatOffset[18][3] *= 0.497433;
	qmatOffset[3][18] *= 0.497433;
	qmatOffset[18][5] *= 1.07071;
	qmatOffset[5][18] *= 1.07071;
	qmatOffset[18][6] *= 0.834267;
	qmatOffset[6][18] *= 0.834267;
	qmatOffset[18][7] *= 0.675128;
	qmatOffset[7][18] *= 0.675128;
	qmatOffset[18][9] *= 2.1139;
	qmatOffset[9][18] *= 2.1139;
	qmatOffset[18][8] *= 0.436898;
	qmatOffset[8][18] *= 0.436898;
	qmatOffset[18][10] *= 1.63857;
	qmatOffset[10][18] *= 1.63857;
	qmatOffset[18][4] *= 4.86017;
	qmatOffset[4][18] *= 4.86017;
	qmatOffset[18][12] *= 0.442935;
	qmatOffset[12][18] *= 0.442935;
	qmatOffset[18][15] *= 1.6641;
	qmatOffset[15][18] *= 1.6641;
	qmatOffset[18][16] *= 0.352251;
	qmatOffset[16][18] *= 0.352251;
	qmatOffset[19][0] *= 0.764894;
	qmatOffset[0][19] *= 0.764894;
	qmatOffset[19][14] *= 1.21225;
	qmatOffset[14][19] *= 1.21225;
	qmatOffset[19][11] *= 3.45058;
	qmatOffset[11][19] *= 3.45058;
	qmatOffset[19][2] *= 1.03489;
	qmatOffset[2][19] *= 1.03489;
	qmatOffset[19][1] *= 1.72794;
	qmatOffset[1][19] *= 1.72794;
	qmatOffset[19][13] *= 0.723509;
	qmatOffset[13][19] *= 0.723509;
	qmatOffset[19][3] *= 0.623719;
	qmatOffset[3][19] *= 0.623719;
	qmatOffset[19][5] *= 0.329184;
	qmatOffset[5][19] *= 0.329184;
	qmatOffset[19][6] *= 12.3072;
	qmatOffset[6][19] *= 12.3072;
	qmatOffset[19][7] *= 1.33502;
	qmatOffset[7][19] *= 1.33502;
	qmatOffset[19][9] *= 1.26654;
	qmatOffset[9][19] *= 1.26654;
	qmatOffset[19][8] *= 0.423423;
	qmatOffset[8][19] *= 0.423423;
	qmatOffset[19][10] *= 1.36128;
	qmatOffset[10][19] *= 1.36128;
	qmatOffset[19][4] *= 20.5074;
	qmatOffset[4][19] *= 20.5074;
	qmatOffset[19][12] *= 0.686449;
	qmatOffset[12][19] *= 0.686449;
	qmatOffset[19][15] *= 2.50053;
	qmatOffset[15][19] *= 2.50053;
	qmatOffset[19][16] *= 0.925072;
	qmatOffset[16][19] *= 0.925072;
	qmatOffset[19][18] *= 7.8969;
	qmatOffset[18][19] *= 7.8969;
	qmatOffset[17][0] *= 6.37375;
	qmatOffset[0][17] *= 6.37375;
	qmatOffset[17][14] *= 0.800207;
	qmatOffset[14][17] *= 0.800207;
	qmatOffset[17][11] *= 0.623538;
	qmatOffset[11][17] *= 0.623538;
	qmatOffset[17][2] *= 0.484018;
	qmatOffset[2][17] *= 0.484018;
	qmatOffset[17][1] *= 3.18413;
	qmatOffset[1][17] *= 3.18413;
	qmatOffset[17][13] *= 0.957268;
	qmatOffset[13][17] *= 0.957268;
	qmatOffset[17][3] *= 1.87059;
	qmatOffset[3][17] *= 1.87059;
	qmatOffset[17][5] *= 0.594945;
	qmatOffset[5][17] *= 0.594945;
	qmatOffset[17][6] *= 0.376062;
	qmatOffset[6][17] *= 0.376062;
	qmatOffset[17][7] *= 24.8508;
	qmatOffset[7][17] *= 24.8508;
	qmatOffset[17][9] *= 5.72027;
	qmatOffset[9][17] *= 5.72027;
	qmatOffset[17][8] *= 0.970464;
	qmatOffset[8][17] *= 0.970464;
	qmatOffset[17][10] *= 6.54037;
	qmatOffset[10][17] *= 6.54037;
	qmatOffset[17][4] *= 2.06492;
	qmatOffset[4][17] *= 2.06492;
	qmatOffset[17][12] *= 1.0005;
	qmatOffset[12][17] *= 1.0005;
	qmatOffset[17][15] *= 0.739488;
	qmatOffset[15][17] *= 0.739488;
	qmatOffset[17][16] *= 4.41086;
	qmatOffset[16][17] *= 4.41086;
	qmatOffset[17][18] *= 1.1609;
	qmatOffset[18][17] *= 1.1609;
	qmatOffset[17][19] *= 1;
	qmatOffset[19][17] *= 1;
	}