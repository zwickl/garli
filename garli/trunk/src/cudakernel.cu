/*
 * cudakernel.cu
 *
 *  Created on: Nov 12, 2008
 *      Author: ayres
 */

#ifndef CUDAKERNELS
#define CUDAKERNELS

__global__ void GarliCLANucleotideNRate(FLOAT_TYPE* CLA, FLOAT_TYPE* Lpr,
		FLOAT_TYPE* Rpr, FLOAT_TYPE* LCL, FLOAT_TYPE* RCL, int nRateCats) {
	// has to be run on blocks of size RATES*4*16, one dimension

	// Block index pre-multiplied
	int bxpm = blockIdx.x * nRateCats * 64;

	// Thread index
	int tx = threadIdx.x;

	// CLAsubL and R are used to store the left and right elements of the block sub-matrix that are computed by the thread
	float CLAsubL = 0;
	float CLAsubR = 0;

	// Declaration of the shared memory array Lprs and Rprs used to store the sub-matrix of Lpr and Rpr
	__shared__ float Lprs[4*16];
	__shared__ float Rprs[4*16];

	// Declaration of the shared memory array LCLs and RCLs used to store the sub-matrix of LCL and RCL
	__shared__ float LCLs[4*4*16];
	__shared__ float RCLs[4*4*16];

	// Load the matrices from device memory to shared memory; each thread loads one element of each matrix
	if (tx < 64) {
		Lprs[tx] = Lpr[tx];
		Rprs[tx] = Rpr[tx];
	}

	// it's not clear to me why this needs to be here but the 1 rate case gives errors without it
	__syncthreads();

	LCLs[tx] = LCL[bxpm + tx];
	RCLs[tx] = RCL[bxpm + tx];

	// Synchronize to make sure the matrices are loaded
	__syncthreads();

	// Multiply the two matrices together, each thread computes one element of the block sub-matrix
	for (int k = 0; k < 4; ++k) {
		CLAsubL += Lprs[4 * (tx & (nRateCats * 4 - 1)) + k] * LCLs[tx - (tx
				& (3)) + k];
		CLAsubR += Rprs[4 * (tx & (nRateCats * 4 - 1)) + k] * RCLs[tx - (tx
				& (3)) + k];
	}

	// Write the block sub-matrix to device memory each thread writes one element
	CLA[bxpm + tx] = CLAsubL * CLAsubR;
}

__global__ void GarliCLAAminoAcidNRate(FLOAT_TYPE* CLA, FLOAT_TYPE* Lpr,
		FLOAT_TYPE* Rpr, FLOAT_TYPE* LCL, FLOAT_TYPE* RCL, int nRateCats) {
	// Block size has to be 20x20

	// Declaration of the shared memory array prs used to store the sub-matrix of Lpr and Rpr
	__shared__ float prs[20][20];

	// Declaration of the shared memory array CLs used to store the sub-matrix of LCL and RCL
	__shared__ float CLs[20][20];

	// Thread index
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	// Block index
	int CLBlockIdx = blockIdx.x * nRateCats * 400;

	for (int i = 0; i < nRateCats; ++i) {

		int prIdx = i * 400 + ty * 20 + tx;
		int CLIdx = i * 20 + nRateCats * 20 * ty + tx;

		float CLAsubL = 0;
		float CLAsubR = 0;

		__syncthreads();

		// Load the matrices from device memory to shared memory; each thread loads one element of each matrix
		prs[tx][ty] = Lpr[prIdx];
		CLs[ty][tx] = LCL[CLBlockIdx + CLIdx];

		// Synchronize to make sure the matrices are loaded
		__syncthreads();

		// Multiply the two matrices together, each thread computes one element of the block sub-matrix
		for (int k = 0; k < 20; ++k)
			CLAsubL += prs[k][tx] * CLs[ty][k];

		// Synchronize to make sure the L side is done
		__syncthreads();

		// Load the matrices from device memory to shared memory; each thread loads one element of each matrix
		prs[tx][ty] = Rpr[prIdx];
		CLs[ty][tx] = RCL[CLBlockIdx + CLIdx];

		// Synchronize to make sure the matrices are loaded
		__syncthreads();

		// Multiply the two matrices together, each thread computes one element of the block sub-matrix
		for (int k = 0; k < 20; ++k)
			CLAsubR += prs[k][tx] * CLs[ty][k];

		// Write the block sub-matrix to device memory each thread writes one element
		CLA[CLIdx + CLBlockIdx] = CLAsubL * CLAsubR;
	}
}

__global__ void GarliCLACodonNRate(FLOAT_TYPE* CLA, FLOAT_TYPE* Lpr,
		FLOAT_TYPE* Rpr, FLOAT_TYPE* LCL, FLOAT_TYPE* RCL, int nRateCats) {
	// Block size has to be 16x16

	// Declaration of the shared memory array prs used to store the sub-matrix of Lpr and Rpr
	__shared__ float prs[16][16];

	// Declaration of the shared memory array CLs used to store the sub-matrix of LCL and RCL
	__shared__ float CLs[16][16];

	// Thread index
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	// Block indexes
	int prBlockIdx = blockIdx.y * 61 * 16;
	int CLBlockIdx = blockIdx.x * nRateCats * 61 * 16;

	for (int i = 0; i < nRateCats; ++i) {

		int prIdx = i * 61 * 61 + ty * 61 + tx;
		int CLIdx = i * 61 + nRateCats * 61 * ty + tx;

		float CLAsubL = 0;
		float CLAsubR = 0;

		int j;
		for (j = 0; j < 48; j += 16) {

			// Synchronize before next round
			__syncthreads();

			// Load the matrices from device memory to shared memory; each thread loads one element of each matrix
			CLs[ty][tx] = LCL[CLBlockIdx + CLIdx + j];
			prs[tx][ty] = Lpr[prBlockIdx + prIdx + j];

			// Synchronize to make sure the matrices are loaded
			__syncthreads();

			// Multiply the two matrices together, each thread computes one element of the block sub-matrix
			for (int k = 0; k < 16; ++k) {
				CLAsubL += prs[k][tx] * CLs[ty][k];
			}

			// Synchronize to make sure the L side is done
			__syncthreads();

			// Load the matrices from device memory to shared memory; each thread loads one element of each matrix
			CLs[ty][tx] = RCL[CLBlockIdx + CLIdx + j];
			prs[tx][ty] = Rpr[prBlockIdx + prIdx + j];

			// Synchronize to make sure the matrices are loaded
			__syncthreads();

			// Multiply the two matrices together, each thread computes one element of the block sub-matrix
			for (int k = 0; k < 16; ++k) {
				CLAsubR += prs[k][tx] * CLs[ty][k];
			}
		}

		// For the last j step k only goes to 13 (61 states)
		__syncthreads();
		CLs[ty][tx] = LCL[CLBlockIdx + CLIdx + j];
		prs[tx][ty] = Lpr[prBlockIdx + prIdx + j];
		__syncthreads();
		for (int k = 0; k < 13; ++k) {
			CLAsubL += prs[k][tx] * CLs[ty][k];
		}
		__syncthreads();
		CLs[ty][tx] = RCL[CLBlockIdx + CLIdx + j];
		prs[tx][ty] = Rpr[prBlockIdx + prIdx + j];
		__syncthreads();
		for (int k = 0; k < 13; ++k) {
			CLAsubR += prs[k][tx] * CLs[ty][k];
		}

		// Write the block sub-matrix to device memory each thread writes one element
		if (blockIdx.y != 3 || tx < 13)
			CLA[CLIdx + (blockIdx.y * 16) + CLBlockIdx] = CLAsubL * CLAsubR;
	}
}

__global__ void GarliDerivNucleotideNRate(FLOAT_TYPE* d_partial,
		FLOAT_TYPE* d_CL1, int* d_partial_underflow_mult,
		int* d_CL1_underflow_mult, FLOAT_TYPE* d_prmat, FLOAT_TYPE* d_d1mat,
		FLOAT_TYPE* d_d2mat, FLOAT_TYPE* d_rateProb, FLOAT_TYPE* d_freqs,
		int* d_countit, int* d_conStates, FLOAT_TYPE* d_Tots_arr,
		int lastConst, bool NoPinvInModel, FLOAT_TYPE prI, int nRateCats) {
	__shared__ float prmat[4*4];
	__shared__ float d1mat[4*4];
	__shared__ float d2mat[4*4];
	__shared__ float CL1[4*128];
	__shared__ float freqs[4];
	__shared__ float rateProb[4];
	__shared__ float Tots[3*128];

	int tx = threadIdx.x;
	int bx = blockIdx.x;
	int dimx = blockDim.x;
	int nstates = 4;
	int i = bx * blockDim.x + tx;

	if (tx < nstates)
		freqs[tx] = d_freqs[tx];
	else if (tx < (nstates + nRateCats))
		rateProb[tx - nstates] = d_rateProb[tx - nstates];

	// calculate remaining chars

	float siteL = 0, siteD1 = 0, siteD2 = 0;
	float tempL, tempD1, tempD2;
	float rateL, rateD1, rateD2;

	d_partial += nstates * nRateCats * i;
	d_CL1 += nstates * nRateCats * i;

	for (int rate = 0; rate < nRateCats; rate++) {
		rateL = rateD1 = rateD2 = 0;
		int rateOffset = rate * nstates * nstates;

		__syncthreads();
		if (tx < 16)
			prmat[tx] = d_prmat[rateOffset + tx];
		else if (tx < 32)
			d1mat[tx - 16] = d_d1mat[rateOffset + tx - 16];
		else if (tx < 48)
			d2mat[tx - 32] = d_d2mat[rateOffset + tx - 32];

		int from = 0;
		tempL = tempD1 = tempD2 = 0;
		for (int to = 0; to < nstates; to++) {
			CL1[nstates * tx + to] = d_CL1[to];
			__syncthreads();
			tempL += prmat[to] * CL1[nstates * tx + to];
			tempD1 += d1mat[to] * CL1[nstates * tx + to];
			tempD2 += d2mat[to] * CL1[nstates * tx + to];
		}
		float partial = d_partial[from];
		rateL += tempL * partial * freqs[from];
		rateD1 += tempD1 * partial * freqs[from];
		rateD2 += tempD2 * partial * freqs[from];

		for (from = 1; from < nstates; from++) {
			tempL = tempD1 = tempD2 = 0;
			int offset = from * nstates;
			for (int to = 0; to < nstates; to++) {
				tempL += prmat[offset + to] * CL1[nstates * tx + to];
				tempD1 += d1mat[offset + to] * CL1[nstates * tx + to];
				tempD2 += d2mat[offset + to] * CL1[nstates * tx + to];
			}
			partial = d_partial[from];
			rateL += tempL * partial * freqs[from];
			rateD1 += tempD1 * partial * freqs[from];
			rateD2 += tempD2 * partial * freqs[from];
		}
		siteL += rateL * rateProb[rate];
		siteD1 += rateD1 * rateProb[rate];
		siteD2 += rateD2 * rateProb[rate];

		d_partial += nstates;
		d_CL1 += nstates;
	}

	float partial_underflow_mult = d_partial_underflow_mult[i];
	float CL1_underflow_mult = d_CL1_underflow_mult[i];

	if ((NoPinvInModel == false) && (i <= lastConst)) {
		float btot = 0.0f;
		int conStates = d_conStates[i];
		if (conStates & 1)
			btot += freqs[0];
		if (conStates & 2)
			btot += freqs[1];
		if (conStates & 4)
			btot += freqs[2];
		if (conStates & 8)
			btot += freqs[3];
		siteL += (prI * btot)
				* exp(partial_underflow_mult + CL1_underflow_mult);
	}

	int countit = d_countit[i];

	Tots[tx] = (log(siteL) - partial_underflow_mult - CL1_underflow_mult)
			* countit;
	siteD1 /= siteL;
	Tots[tx + dimx] = countit * siteD1;
	Tots[tx + dimx * 2] = countit * ((siteD2 / siteL) - siteD1 * siteD1);
	__syncthreads();
	for (unsigned int s = dimx / 2; s > 0; s >>= 1) {
		if (tx < s) {
			Tots[tx] += Tots[tx + s];
			Tots[tx + dimx] += Tots[tx + dimx + s];
			Tots[tx + dimx * 2] += Tots[tx + dimx * 2 + s];
		}
		__syncthreads();
	}
	if (tx == 0) {
		d_Tots_arr[bx] = Tots[0];
		d_Tots_arr[bx + gridDim.x] = Tots[0 + dimx];
		d_Tots_arr[bx + gridDim.x * 2] = Tots[0 + dimx * 2];
	}

}

__global__ void GarliDerivAminoAcidNRate(FLOAT_TYPE* d_partial,
		FLOAT_TYPE* d_CL1, int* d_partial_underflow_mult,
		int* d_CL1_underflow_mult, FLOAT_TYPE* d_prmat, FLOAT_TYPE* d_d1mat,
		FLOAT_TYPE* d_d2mat, FLOAT_TYPE* d_rateProb, FLOAT_TYPE* d_freqs,
		int* d_countit, int* d_conStates, FLOAT_TYPE* d_Tots_arr,
		int lastConst, bool NoPinvInModel, FLOAT_TYPE prI, int nRateCats) {
	__shared__ float prmat[20*20];
	__shared__ float d1mat[20*20];
	__shared__ float d2mat[20*20];
	__shared__ float CL1[20*128];
	__shared__ float freqs[20];
	__shared__ float rateProb[4];
	__shared__ float Tots[128];

	int tx = threadIdx.x;
	int bx = blockIdx.x;
	int dimx = blockDim.x;
	int nstates = 20;
	int i = bx * blockDim.x + tx;

	if (tx < nstates)
		freqs[tx] = d_freqs[tx];
	else if (tx < (nstates + nRateCats))
		rateProb[tx - nstates] = d_rateProb[tx - nstates];

	// calculate remaining chars

	float siteL = 0, siteD1 = 0, siteD2 = 0;
	float tempL, tempD1, tempD2;
	float rateL, rateD1, rateD2;

	d_partial += nstates * nRateCats * i;
	d_CL1 += nstates * nRateCats * i;

	for (int rate = 0; rate < nRateCats; rate++) {
		rateL = rateD1 = rateD2 = 0;
		int rateOffset = rate * nstates * nstates;

		__syncthreads();
		for (int j = 0; j < 3; j++) {
			prmat[tx + dimx * j] = d_prmat[rateOffset + tx + dimx * j];
			d1mat[tx + dimx * j] = d_d1mat[rateOffset + tx + dimx * j];
			d2mat[tx + dimx * j] = d_d2mat[rateOffset + tx + dimx * j];
		}
		if (tx < 16)
			prmat[tx + dimx * 3] = d_prmat[rateOffset + tx + dimx * 3];
		else if (tx < 32)
			d1mat[tx - 16 + dimx * 3]
					= d_d1mat[rateOffset + tx - 16 + dimx * 3];
		else if (tx < 48)
			d2mat[tx - 32 + dimx * 3]
					= d_d2mat[rateOffset + tx - 32 + dimx * 3];

		int from = 0;
		tempL = tempD1 = tempD2 = 0;
		for (int to = 0; to < nstates; to++) {
			CL1[nstates * tx + to] = d_CL1[to];
			__syncthreads();
			tempL += prmat[to] * CL1[nstates * tx + to];
			tempD1 += d1mat[to] * CL1[nstates * tx + to];
			tempD2 += d2mat[to] * CL1[nstates * tx + to];
		}
		float partial = d_partial[from];
		rateL += tempL * partial * freqs[from];
		rateD1 += tempD1 * partial * freqs[from];
		rateD2 += tempD2 * partial * freqs[from];

		for (from = 1; from < nstates; from++) {
			tempL = tempD1 = tempD2 = 0;
			int offset = from * nstates;
			for (int to = 0; to < nstates; to++) {
				tempL += prmat[offset + to] * CL1[nstates * tx + to];
				tempD1 += d1mat[offset + to] * CL1[nstates * tx + to];
				tempD2 += d2mat[offset + to] * CL1[nstates * tx + to];
			}
			partial = d_partial[from];
			rateL += tempL * partial * freqs[from];
			rateD1 += tempD1 * partial * freqs[from];
			rateD2 += tempD2 * partial * freqs[from];
		}
		siteL += rateL * rateProb[rate];
		siteD1 += rateD1 * rateProb[rate];
		siteD2 += rateD2 * rateProb[rate];

		d_partial += nstates;
		d_CL1 += nstates;
	}

	float partial_underflow_mult = d_partial_underflow_mult[i];
	float CL1_underflow_mult = d_CL1_underflow_mult[i];

	if ((NoPinvInModel == false) && (i <= lastConst))
		siteL += (prI * freqs[d_conStates[i]] * exp(partial_underflow_mult
				+ CL1_underflow_mult));

	int countit = d_countit[i];

	Tots[tx] = (log(siteL) - partial_underflow_mult - CL1_underflow_mult)
			* countit;
	__syncthreads();
	for (unsigned int s = dimx / 2; s > 0; s >>= 1) {
		if (tx < s)
			Tots[tx] += Tots[tx + s];
		__syncthreads();
	}
	if (tx == 0)
		d_Tots_arr[bx] = Tots[0];

	siteD1 /= siteL;
	Tots[tx] = countit * siteD1;
	__syncthreads();
	for (unsigned int s = dimx / 2; s > 0; s >>= 1) {
		if (tx < s)
			Tots[tx] += Tots[tx + s];
		__syncthreads();
	}
	if (tx == 0)
		d_Tots_arr[bx + gridDim.x] = Tots[0];

	Tots[tx] = countit * ((siteD2 / siteL) - siteD1 * siteD1);
	__syncthreads();
	for (unsigned int s = dimx / 2; s > 0; s >>= 1) {
		if (tx < s)
			Tots[tx] += Tots[tx + s];
		__syncthreads();
	}
	if (tx == 0)
		d_Tots_arr[bx + gridDim.x * 2] = Tots[0];

}

__global__ void GarliDerivCodonNRate(FLOAT_TYPE* d_partial, FLOAT_TYPE* d_CL1,
		int* d_partial_underflow_mult, int* d_CL1_underflow_mult,
		FLOAT_TYPE* d_prmat, FLOAT_TYPE* d_d1mat, FLOAT_TYPE* d_d2mat,
		FLOAT_TYPE* d_rateProb, FLOAT_TYPE* d_freqs, int* d_countit,
		int* d_conStates, FLOAT_TYPE* d_Tots_arr, int lastConst,
		bool NoPinvInModel, FLOAT_TYPE prI, int nRateCats) {
	__shared__ float freqs[61];
	__shared__ float rateProb[4];
	__shared__ float Tots[3*128];

	int tx = threadIdx.x;
	int bx = blockIdx.x;
	int dimx = blockDim.x;
	int nstates = 61;
	int i = bx * blockDim.x + tx;

	if (tx < nstates)
		freqs[tx] = d_freqs[tx];
	else if (tx < (nstates + nRateCats))
		rateProb[tx - nstates] = d_rateProb[tx - nstates];

	__syncthreads();

	// calculate remaining chars

	float siteL = 0, siteD1 = 0, siteD2 = 0;
	float tempL, tempD1, tempD2;
	float rateL, rateD1, rateD2;

	d_partial += nstates * nRateCats * i;
	d_CL1 += nstates * nRateCats * i;

	for (int rate = 0; rate < nRateCats; rate++) {
		rateL = rateD1 = rateD2 = 0;
		int rateOffset = rate * nstates * nstates;

		for (int from = 0; from < nstates; from++) {
			tempL = tempD1 = tempD2 = 0;
			int offset = from * nstates;
			for (int to = 0; to < nstates; to++) {
				float CL1 = d_CL1[to];
				tempL += d_prmat[rateOffset + offset + to] * CL1;
				tempD1 += d_d1mat[rateOffset + offset + to] * CL1;
				tempD2 += d_d2mat[rateOffset + offset + to] * CL1;
			}
			float partial = d_partial[from];
			rateL += tempL * partial * freqs[from];
			rateD1 += tempD1 * partial * freqs[from];
			rateD2 += tempD2 * partial * freqs[from];
		}
		siteL += rateL * rateProb[rate];
		siteD1 += rateD1 * rateProb[rate];
		siteD2 += rateD2 * rateProb[rate];

		d_partial += nstates;
		d_CL1 += nstates;
	}

	float partial_underflow_mult = d_partial_underflow_mult[i];
	float CL1_underflow_mult = d_CL1_underflow_mult[i];

	if ((NoPinvInModel == false) && (i <= lastConst))
		siteL += (prI * freqs[d_conStates[i]] * exp(partial_underflow_mult
				+ CL1_underflow_mult));

	int countit = d_countit[i];

	Tots[tx] = (log(siteL) - partial_underflow_mult - CL1_underflow_mult)
			* countit;
	siteD1 /= siteL;
	Tots[tx + dimx] = countit * siteD1;
	Tots[tx + dimx * 2] = countit * ((siteD2 / siteL) - siteD1 * siteD1);
	__syncthreads();
	for (unsigned int s = dimx / 2; s > 0; s >>= 1) {
		if (tx < s) {
			Tots[tx] += Tots[tx + s];
			Tots[tx + dimx] += Tots[tx + dimx + s];
			Tots[tx + dimx * 2] += Tots[tx + dimx * 2 + s];
		}
		__syncthreads();
	}
	if (tx == 0) {
		d_Tots_arr[bx] = Tots[0];
		d_Tots_arr[bx + gridDim.x] = Tots[0 + dimx];
		d_Tots_arr[bx + gridDim.x * 2] = Tots[0 + dimx * 2];
	}
}

#endif// #ifndef CUDAKERNELS
