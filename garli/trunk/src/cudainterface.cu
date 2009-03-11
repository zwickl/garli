/*
 * cudainterface.cu
 *
 *  Created on: Nov 12, 2008
 *      Author: ayres
 */


#include <cutil.h>

#include "defs.h"
#include "outputman.h"

// includes, kernels
#include "cudakernel.cu"

extern OutputManager outman;

extern "C"
bool CheckCuda() {
    int deviceCount;
    bool cuda_support = false;

    CUDA_SAFE_CALL(cudaGetDeviceCount(&deviceCount));
    if (deviceCount == 0)
    	outman.UserMessageNoCR("There is no device supporting CUDA\n");
    else
    	cuda_support = true;

    return cuda_support;
}

extern "C"
void DeviceQuery() {

	outman.UserMessageNoCR(
			"========================= GPU Device Query =========================\n");

    int deviceCount;
    CUDA_SAFE_CALL(cudaGetDeviceCount(&deviceCount));
    if (deviceCount == 0)
    	outman.UserMessageNoCR("There is no device supporting CUDA\n");
    int dev;
    for (dev = 0; dev < deviceCount; ++dev) {
        cudaDeviceProp deviceProp;
        CUDA_SAFE_CALL(cudaGetDeviceProperties(&deviceProp, dev));
        if (dev == 0) {
            if (deviceProp.major == 9999 && deviceProp.minor == 9999)
            	outman.UserMessageNoCR("There is no device supporting CUDA.\n");
            else if (deviceCount == 1)
            	outman.UserMessageNoCR("There is 1 device supporting CUDA\n");
            else
            	outman.UserMessageNoCR("There are %d devices supporting CUDA\n", deviceCount);
        }
        outman.UserMessageNoCR("\nDevice %d: \"%s\"\n", dev, deviceProp.name);
        outman.UserMessageNoCR("  Major revision number:                         %d\n",
               deviceProp.major);
        outman.UserMessageNoCR("  Minor revision number:                         %d\n",
               deviceProp.minor);
        outman.UserMessageNoCR("  Total amount of global memory:                 %u bytes\n",
               deviceProp.totalGlobalMem);
    #if CUDART_VERSION >= 2000
        outman.UserMessageNoCR("  Number of multiprocessors:                     %d\n",
               deviceProp.multiProcessorCount);
        outman.UserMessageNoCR("  Number of cores:                               %d\n",
               8 * deviceProp.multiProcessorCount);
    #endif
        outman.UserMessageNoCR("  Total amount of constant memory:               %u bytes\n",
               deviceProp.totalConstMem);
        outman.UserMessageNoCR("  Total amount of shared memory per block:       %u bytes\n",
               deviceProp.sharedMemPerBlock);
        outman.UserMessageNoCR("  Total number of registers available per block: %d\n",
               deviceProp.regsPerBlock);
        outman.UserMessageNoCR("  Warp size:                                     %d\n",
               deviceProp.warpSize);
        outman.UserMessageNoCR("  Maximum number of threads per block:           %d\n",
               deviceProp.maxThreadsPerBlock);
        outman.UserMessageNoCR("  Maximum sizes of each dimension of a block:    %d x %d x %d\n",
               deviceProp.maxThreadsDim[0],
               deviceProp.maxThreadsDim[1],
               deviceProp.maxThreadsDim[2]);
        outman.UserMessageNoCR("  Maximum sizes of each dimension of a grid:     %d x %d x %d\n",
               deviceProp.maxGridSize[0],
               deviceProp.maxGridSize[1],
               deviceProp.maxGridSize[2]);
        outman.UserMessageNoCR("  Maximum memory pitch:                          %u bytes\n",
               deviceProp.memPitch);
        outman.UserMessageNoCR("  Texture alignment:                             %u bytes\n",
               deviceProp.textureAlignment);
        outman.UserMessageNoCR("  Clock rate:                                    %.2f GHz\n",
               deviceProp.clockRate * 1e-6f);
    #if CUDART_VERSION >= 2000
        outman.UserMessageNoCR("  Concurrent copy and execution:                 %s\n",
               deviceProp.deviceOverlap ? "Yes" : "No");
    #endif
    }

	outman.UserMessageNoCR(
				"====================================================================\n\n");
}


extern "C"
void SetDevice(unsigned int device_number) {
	cudaSetDevice(device_number);
	    cudaDeviceProp deviceProp;
	    cudaGetDeviceProperties(&deviceProp, device_number);
	    outman.UserMessageNoCR ("Using GPU device %d: %s\n\n", device_number,deviceProp.name);
}

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
void CuComputeGPUCLA(FLOAT_TYPE* h_Lpr, FLOAT_TYPE* h_Rpr, FLOAT_TYPE* h_LCL, FLOAT_TYPE* h_RCL, FLOAT_TYPE* h_CLA,
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
void CuComputeGPUDeriv(const FLOAT_TYPE* h_partial, const FLOAT_TYPE* h_CL1, const int* h_partial_underflow_mult,
		const int* h_CL1_underflow_mult, const FLOAT_TYPE* h_prmat, const FLOAT_TYPE* h_d1mat, const FLOAT_TYPE* h_d2mat,
		const FLOAT_TYPE* h_rateProb, const FLOAT_TYPE* h_freqs, const int* h_countit, const int* h_conStates,
		FLOAT_TYPE* h_Tots, FLOAT_TYPE* h_Tots_arr, int* h_nchar_boot_index,
		FLOAT_TYPE* d_partial, FLOAT_TYPE* d_CL1, int* d_partial_underflow_mult,
		int* d_CL1_underflow_mult, FLOAT_TYPE* d_prmat, FLOAT_TYPE* d_d1mat, FLOAT_TYPE* d_d2mat,
		FLOAT_TYPE* d_rateProb, FLOAT_TYPE* d_freqs, int* d_countit, int* d_conStates,
		FLOAT_TYPE* d_Tots, FLOAT_TYPE* d_Tots_arr, int* d_nchar_boot_index,
		unsigned int mem_size_pr, unsigned int mem_size_CL, unsigned int mem_size_int_char,
		unsigned int mem_size_rates, unsigned int mem_size_states, unsigned int mem_size_Tots,
		unsigned int mem_size_Tots_arr, unsigned int mem_size_nchar_boot_index, int lastConst,
		bool NoPinvInModel, FLOAT_TYPE prI,
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
	CUDA_SAFE_CALL(cudaMemcpy(d_nchar_boot_index, h_nchar_boot_index, mem_size_nchar_boot_index, cudaMemcpyHostToDevice));

	// run kernel
	if (nstates == 4)
	GarliDerivNucleotideNRate<<< dimGrid, dimBlock >>>(d_partial, d_CL1, d_partial_underflow_mult,
			d_CL1_underflow_mult, d_prmat, d_d1mat, d_d2mat,
			d_rateProb, d_freqs, d_countit, d_conStates, d_Tots_arr, d_nchar_boot_index,
			lastConst, NoPinvInModel, prI, nRateCats);
	else if (nstates == 20)
	GarliDerivAminoAcidNRate<<< dimGrid, dimBlock >>>(d_partial, d_CL1, d_partial_underflow_mult,
			d_CL1_underflow_mult, d_prmat, d_d1mat, d_d2mat,
			d_rateProb, d_freqs, d_countit, d_conStates, d_Tots_arr, d_nchar_boot_index,
			lastConst, NoPinvInModel, prI, nRateCats);
	else
	GarliDerivCodonNRate<<< dimGrid, dimBlock >>>(d_partial, d_CL1, d_partial_underflow_mult,
			d_CL1_underflow_mult, d_prmat, d_d1mat, d_d2mat,
			d_rateProb, d_freqs, d_countit, d_conStates, d_Tots_arr, d_nchar_boot_index,
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

		if((NoPinvInModel == false) && (h_nchar_boot_index[i]<=lastConst)) {
//			if (nstates == 4) {
//				float btot = 0.0f;
//				if (h_conStates[h_nchar_boot_index[i]] & 1)
//				btot += h_freqs[0];
//				if (h_conStates[h_nchar_boot_index[i]] & 2)
//				btot += h_freqs[1];
//				if (h_conStates[h_nchar_boot_index[i]] & 4)
//				btot += h_freqs[2];
//				if (h_conStates[h_nchar_boot_index[i]] & 8)
//				btot += h_freqs[3];
//				siteL += (prI * btot) * exp(h_partial_underflow_mult[h_nchar_boot_index[i]]
//						+ h_CL1_underflow_mult[h_nchar_boot_index[i]]);
//			} else
			siteL += (prI*h_freqs[h_conStates[h_nchar_boot_index[i]]] * exp((FLOAT_TYPE)h_partial_underflow_mult[h_nchar_boot_index[i]]) * exp((FLOAT_TYPE)h_CL1_underflow_mult[h_nchar_boot_index[i]]));
		}

		totL += (log(siteL) - h_partial_underflow_mult[h_nchar_boot_index[i]] - h_CL1_underflow_mult[h_nchar_boot_index[i]]) * h_countit[h_nchar_boot_index[i]];
		siteD1 /= siteL;
		tot1 += h_countit[h_nchar_boot_index[i]] * siteD1;
		tot2 += h_countit[h_nchar_boot_index[i]] * ((siteD2 / siteL) - siteD1*siteD1);
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

