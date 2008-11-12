/*
 * cudainterface.cu
 *
 *  Created on: Nov 12, 2008
 *      Author: ayres
 */

#include "cutil.h"
#include <iostream>

// includes, kernels
#include "cudakernel.cu"

extern "C"
void AllocatePinnedMemory(void** arr, unsigned int mem_size_bytes) {
	// allocate host pinned memory
	CUDA_SAFE_CALL(cudaMallocHost((void**)&(*arr), mem_size_bytes));

}

extern "C"
void AllocateGPU(void** arr, unsigned int mem_size_bytes) {
	// allocate device memory
	CUDA_SAFE_CALL(cudaMalloc((void**) &(*arr), mem_size_bytes));
}

extern "C"
void RunGPUCLA(FLOAT_TYPE* h_Lpr, FLOAT_TYPE* h_Rpr, FLOAT_TYPE* h_LCL, FLOAT_TYPE* h_RCL, FLOAT_TYPE* h_CLA,
		FLOAT_TYPE* d_Lpr, FLOAT_TYPE* d_Rpr, FLOAT_TYPE* d_LCL, FLOAT_TYPE* d_RCL, FLOAT_TYPE* d_CLA,
		unsigned int mem_size_pr, unsigned int mem_size_CL,
		int nstates, int nRateCats, int nchar, int ncharGPU, dim3 dimBlock, dim3 dimGrid) {
	// copy matrices to the device
	CUDA_SAFE_CALL(cudaMemcpy(d_Lpr, h_Lpr, mem_size_pr, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_Rpr, h_Rpr, mem_size_pr, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_LCL, h_LCL, mem_size_CL, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_RCL, h_RCL, mem_size_CL, cudaMemcpyHostToDevice));

	// run kernel
	if (nstates == 4)
		GarliCLANucleotideNRate<<< dimGrid, dimBlock >>>(d_CLA, d_Lpr, d_Rpr, d_LCL, d_RCL, nRateCats);
	else if (nstates == 20)
		GarliCLAAminoAcidNRate<<< dimGrid, dimBlock >>>(d_CLA, d_Lpr, d_Rpr, d_LCL, d_RCL, nRateCats);
	else
		GarliCLACodonNRate<<< dimGrid, dimBlock >>>(d_CLA, d_Lpr, d_Rpr, d_LCL, d_RCL, nRateCats);

	//check if kernel execution generated and error
	CUT_CHECK_ERROR("Kernel execution failed");

	// calculate remaining chars
	FLOAT_TYPE L1, R1;
	FLOAT_TYPE *dest = h_CLA + nstates * nRateCats * ncharGPU;
	h_LCL += nstates * nRateCats * ncharGPU;
	h_RCL += nstates * nRateCats * ncharGPU;
	for(unsigned int i=ncharGPU;i<nchar;i++) {
		for(unsigned int rate=0;rate<nRateCats;rate++) {
			for(unsigned int from=0;from<nstates;from++) {
				L1 = R1 = 0;
				for(unsigned int to=0;to<nstates;to++) {
					L1 += h_Lpr[rate*nstates*nstates + from*nstates + to] * h_LCL[to];
					R1 += h_Rpr[rate*nstates*nstates + from*nstates + to] * h_RCL[to];
				}
				dest[from] = L1 * R1;
			}
			h_LCL += nstates;
			h_RCL += nstates;
			dest += nstates;
		}
	}

	// copy result back to host
	CUDA_SAFE_CALL(cudaMemcpy(h_CLA, d_CLA, mem_size_CL, cudaMemcpyDeviceToHost));
}

extern "C"
void RunGPUDeriv(FLOAT_TYPE* h_partial, FLOAT_TYPE* h_CL1, int* h_partial_underflow_mult,
		int* h_CL1_underflow_mult, FLOAT_TYPE* h_prmat, FLOAT_TYPE* h_d1mat, FLOAT_TYPE* h_d2mat,
		FLOAT_TYPE* h_rateProb, FLOAT_TYPE* h_freqs, int* h_countit, int* h_conStates,
		FLOAT_TYPE* h_Tots, FLOAT_TYPE* h_Tots_arr,
		FLOAT_TYPE* d_partial, FLOAT_TYPE* d_CL1, int* d_partial_underflow_mult,
		int* d_CL1_underflow_mult, FLOAT_TYPE* d_prmat, FLOAT_TYPE* d_d1mat, FLOAT_TYPE* d_d2mat,
		FLOAT_TYPE* d_rateProb, FLOAT_TYPE* d_freqs, int* d_countit, int* d_conStates,
		FLOAT_TYPE* d_Tots, FLOAT_TYPE* d_Tots_arr,
		unsigned int mem_size_pr, unsigned int mem_size_CL, unsigned int mem_size_int_char,
		unsigned int mem_size_rates, unsigned int mem_size_states, unsigned int mem_size_Tots,
		unsigned int mem_size_Tots_arr, int lastConst, bool NoPinvInModel, FLOAT_TYPE prI,
		int nstates, int nRateCats, int nchar, int ncharGPU, dim3 dimBlock, dim3 dimGrid) {

	// copy matrices to the device
	CUDA_SAFE_CALL(cudaMemcpy(d_partial, h_partial, mem_size_CL, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_CL1, h_CL1, mem_size_CL, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_partial_underflow_mult, h_partial_underflow_mult, mem_size_int_char, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_CL1_underflow_mult, h_CL1_underflow_mult, mem_size_int_char, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_prmat, h_prmat, mem_size_pr, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_d1mat, h_d1mat, mem_size_pr, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_d2mat, h_d2mat, mem_size_pr, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_rateProb, h_rateProb, mem_size_rates, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_freqs, h_freqs, mem_size_states, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_countit, h_countit, mem_size_int_char, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_conStates, h_conStates, mem_size_int_char, cudaMemcpyHostToDevice));

	// run kernel
	if (nstates == 4)
		GarliDerivNucleotideNRate<<< dimGrid, dimBlock >>>(d_partial, d_CL1, d_partial_underflow_mult,
				d_CL1_underflow_mult, d_prmat, d_d1mat, d_d2mat,
				d_rateProb, d_freqs, d_countit, d_conStates, d_Tots_arr,
				lastConst, NoPinvInModel, prI, nRateCats);
	else if (nstates == 20)
		GarliDerivAminoAcidNRate<<< dimGrid, dimBlock >>>(d_partial, d_CL1, d_partial_underflow_mult,
				d_CL1_underflow_mult, d_prmat, d_d1mat, d_d2mat,
				d_rateProb, d_freqs, d_countit, d_conStates, d_Tots_arr,
				lastConst, NoPinvInModel, prI, nRateCats);
	else
		GarliDerivCodonNRate<<< dimGrid, dimBlock >>>(d_partial, d_CL1, d_partial_underflow_mult,
				d_CL1_underflow_mult, d_prmat, d_d1mat, d_d2mat,
				d_rateProb, d_freqs, d_countit, d_conStates, d_Tots_arr,
				lastConst, NoPinvInModel, prI, nRateCats);

	//check if kernel execution generated and error
	CUT_CHECK_ERROR("Kernel execution failed");

	// calculate remaining chars
	FLOAT_TYPE tot1=0, tot2=0, totL = 0;

	FLOAT_TYPE siteL, siteD1, siteD2;
	FLOAT_TYPE tempL, tempD1, tempD2;
	FLOAT_TYPE rateL, rateD1, rateD2;

	h_partial += nstates * nRateCats * ncharGPU;
	h_CL1 += nstates * nRateCats * ncharGPU;

	for(int i=ncharGPU;i<nchar;i++) {
		siteL = siteD1 = siteD2 = 0;
		for(int rate=0;rate<nRateCats;rate++) {
			rateL = rateD1 = rateD2 = 0;
			int rateOffset = rate*nstates*nstates;
			for(int from=0;from<nstates;from++) {
				tempL = tempD1 = tempD2 = 0;
				int offset = from * nstates;
				for(int to=0;to<nstates;to++) {
					tempL += h_prmat[rateOffset + offset + to]*h_CL1[to];
					tempD1 += h_d1mat[rateOffset + offset + to]*h_CL1[to];
					tempD2 += h_d2mat[rateOffset + offset + to]*h_CL1[to];
				}
				rateL += tempL * h_partial[from] * h_freqs[from];
				rateD1 += tempD1 * h_partial[from] * h_freqs[from];
				rateD2 += tempD2 * h_partial[from] * h_freqs[from];
			}
			siteL += rateL * h_rateProb[rate];
			siteD1 += rateD1 * h_rateProb[rate];
			siteD2 += rateD2 * h_rateProb[rate];
			h_partial += nstates;
			h_CL1 += nstates;
		}

		if((NoPinvInModel == false) && (i<=lastConst))
		siteL += (prI*h_freqs[h_conStates[i]] * exp((FLOAT_TYPE)h_partial_underflow_mult[i]) * exp((FLOAT_TYPE)h_CL1_underflow_mult[i]));

		totL += (log(siteL) - h_partial_underflow_mult[i] - h_CL1_underflow_mult[i]) * h_countit[i];
		siteD1 /= siteL;
		tot1 += h_countit[i] * siteD1;
		tot2 += h_countit[i] * ((siteD2 / siteL) - siteD1*siteD1);
	}

	// copy result back to host
	CUDA_SAFE_CALL(cudaMemcpy(h_Tots_arr, d_Tots_arr, mem_size_Tots_arr, cudaMemcpyDeviceToHost));

	// clear previous results
	h_Tots[0] = 0;
	h_Tots[1] = 0;
	h_Tots[2] = 0;

	for (int i=0;i<dimGrid.x;i++) {
		h_Tots[0] += h_Tots_arr[i];
		h_Tots[1] += h_Tots_arr[i+dimGrid.x];
		h_Tots[2] += h_Tots_arr[i+dimGrid.x*2];
	}

	// add up the results from the remaining chars
	h_Tots[0] += totL;
	h_Tots[1] += tot1;
	h_Tots[2] += tot2;

}

extern "C"
void FreeGPU(FLOAT_TYPE* arr) {
	// free device memory
	CUDA_SAFE_CALL(cudaFree(arr));
}

extern "C"
void FreePinnedMemory(void* arr) {
	// free host pinned memory
	CUDA_SAFE_CALL(cudaFreeHost(arr));

}
