/*
 * cudaman.cpp
 *
 *  Created on: Nov 12, 2008
 *      Author: ayres
 */

using namespace std;

#include "cudaman.h"
#include "defs.h"
#include "outputman.h"
#include "cutil.h"
#include "vector_types.h"

extern OutputManager outman;

////////////////////////////////////////////////////////////////////////////////
// declaration, forward
extern "C" {
void AllocatePinnedMemory(void** arr, unsigned int mem_size_bytes);

void AllocateGPU(void** arr, unsigned int mem_size_bytes);

void FreeGPU(void* arr);

void FreePinnedMemory(void* arr);

void RunGPUCLA(FLOAT_TYPE* h_Lpr, FLOAT_TYPE* h_Rpr, FLOAT_TYPE* h_LCL,
		FLOAT_TYPE* h_RCL, FLOAT_TYPE* h_CLA, FLOAT_TYPE* d_Lpr,
		FLOAT_TYPE* d_Rpr, FLOAT_TYPE* d_LCL, FLOAT_TYPE* d_RCL,
		FLOAT_TYPE* d_CLA, unsigned int mem_size_pr, unsigned int mem_size_CL,
		int nstates, int nRateCats, int nchar, int ncharGPU, dim3 dimBlock,
		dim3 dimGrid);

void RunGPUDeriv(FLOAT_TYPE* h_partial, FLOAT_TYPE* h_CL1,
		int* h_partial_underflow_mult, int* h_CL1_underflow_mult,
		FLOAT_TYPE* h_prmat, FLOAT_TYPE* h_d1mat, FLOAT_TYPE* h_d2mat,
		FLOAT_TYPE* h_rateProb, FLOAT_TYPE* h_freqs, int* h_countit,
		int* h_conStates, FLOAT_TYPE* h_Tots, FLOAT_TYPE* h_Tots_arr,
		FLOAT_TYPE* d_partial, FLOAT_TYPE* d_CL1,
		int* d_partial_underflow_mult, int* d_CL1_underflow_mult,
		FLOAT_TYPE* d_prmat, FLOAT_TYPE* d_d1mat, FLOAT_TYPE* d_d2mat,
		FLOAT_TYPE* d_rateProb, FLOAT_TYPE* d_freqs, int* d_countit,
		int* d_conStates, FLOAT_TYPE* d_Tots, FLOAT_TYPE* d_Tots_arr,
		unsigned int mem_size_pr, unsigned int mem_size_CL,
		unsigned int mem_size_int_char, unsigned int mem_size_rates,
		unsigned int mem_size_states, unsigned int mem_size_Tots,
		unsigned int mem_size_Tots_arr, int lastConst, bool NoPinvInModel,
		FLOAT_TYPE prI, int nstates, int nRateCats, int nchar, int ncharGPU,
		dim3 dimBlock, dim3 dimGrid);
}
////////////////////////////////////////////////////////////////////////////////


CudaManager::CudaManager() {
}

CudaManager::CudaManager(int nstatesIn, int numRateCatsIn, int ncharIn) {
	nstates = nstatesIn;
	numRateCats = numRateCatsIn;
	nchar = ncharIn;
}

CudaManager::~CudaManager() {
}

void CudaManager::BenchmarkGPU() {
	BenchmarkCLA();
}

// Private methods
void CudaManager::RandomInit(FLOAT_TYPE* data, int size) {
	for (int i = 0; i < size; ++i)
		data[i] = rand() / (FLOAT_TYPE) RAND_MAX;
}

void CudaManager::ComputeGarliCLA(FLOAT_TYPE* dest, const FLOAT_TYPE* Lpr,
		const FLOAT_TYPE* Rpr, const FLOAT_TYPE* LCL, const FLOAT_TYPE* RCL) {
	FLOAT_TYPE L1, R1;
	for (unsigned int i = 0; i < nchar; i++) {
		for (unsigned int rate = 0; rate < numRateCats; rate++) {
			for (unsigned int from = 0; from < nstates; from++) {
				L1 = R1 = 0;
				for (unsigned int to = 0; to < nstates; to++) {
					L1 += Lpr[rate * nstates * nstates + from * nstates + to]
							* LCL[to];
					R1 += Rpr[rate * nstates * nstates + from * nstates + to]
							* RCL[to];
				}
				dest[from + (i * numRateCats * nstates + rate * nstates)] = L1
						* R1;
			}
			LCL += nstates;
			RCL += nstates;
		}
	}
}

void CudaManager::BenchmarkCLA() {
	printf("\ntest \t");
	printf("memory \t");
	printf("float \t");
	printf("states\t");
	printf("rates\t");
	printf("chars\t");
	printf("cpu   \t");
	printf("gpu   \t");
	printf("speedup\n");

	FLOAT_TYPE *h_Lpr, *h_Rpr, *h_LCL, *h_RCL, *h_CLA, *reference_CLA;
	FLOAT_TYPE *d_Lpr, *d_Rpr, *d_LCL, *d_RCL, *d_CLA;

	int ncharGPU;
	if (nstates == 20)
		ncharGPU = 20 * ((int) nchar / 20);
	else
		ncharGPU = 16 * ((int) nchar / 16);

	int size_pr = nstates * nstates * numRateCats;
	int size_CL = nstates * ncharGPU * numRateCats;
	int size_CL_host = nstates * nchar * numRateCats;
	unsigned int mem_size_pr = sizeof(FLOAT_TYPE) * size_pr;
	unsigned int mem_size_CL = sizeof(FLOAT_TYPE) * size_CL;
	unsigned int mem_size_CL_host = sizeof(FLOAT_TYPE) * size_CL_host;

	// setup block and grid size for gpu
	dim3 dimBlock;
	dim3 dimGrid;

	if (nstates == 4) {
		dimBlock = dim3(64 * numRateCats);
		dimGrid = dim3(ncharGPU / 16);
	} else if (nstates == 20) {
		dimBlock = dim3(20, 20);
		dimGrid = dim3(ncharGPU / 20);
	} else {
		dimBlock = dim3(16, 16);
		dimGrid = dim3(ncharGPU / 16, 4);
	}

	reference_CLA = (FLOAT_TYPE*) malloc(mem_size_CL_host);

	// allocate host memory
	AllocatePinnedMemory((void**) &h_Lpr, mem_size_pr);
	AllocatePinnedMemory((void**) &h_Rpr, mem_size_pr);
	AllocatePinnedMemory((void**) &h_LCL, mem_size_CL_host);
	AllocatePinnedMemory((void**) &h_RCL, mem_size_CL_host);
	AllocatePinnedMemory((void**) &h_CLA, mem_size_CL_host);

	// allocate memory on gpu
	if (nstates == 61) {
		unsigned int mem_size_pr_64 = sizeof(FLOAT_TYPE) * 64 * 64
				* numRateCats;
		unsigned int mem_size_CL_64 = sizeof(FLOAT_TYPE) * 64 * ncharGPU
				* numRateCats;

		AllocateGPU((void**) &d_Lpr, mem_size_pr_64);
		AllocateGPU((void**) &d_Rpr, mem_size_pr_64);
		AllocateGPU((void**) &d_LCL, mem_size_CL_64);
		AllocateGPU((void**) &d_RCL, mem_size_CL_64);
		AllocateGPU((void**) &d_CLA, mem_size_CL_64);
	} else {
		AllocateGPU((void**) &d_Lpr, mem_size_pr);
		AllocateGPU((void**) &d_Rpr, mem_size_pr);
		AllocateGPU((void**) &d_LCL, mem_size_CL);
		AllocateGPU((void**) &d_RCL, mem_size_CL);
		AllocateGPU((void**) &d_CLA, mem_size_CL);
	}

	// initialize host memory with random data
	srand(42);
	RandomInit(h_Lpr, size_pr);
	RandomInit(h_Rpr, size_pr);
	RandomInit(h_LCL, size_CL_host);
	RandomInit(h_RCL, size_CL_host);

	// create and start GPU timer
	unsigned int gpu_timer = 0;
	CUT_SAFE_CALL(cutCreateTimer(&gpu_timer));
	CUT_SAFE_CALL(cutStartTimer(gpu_timer));

	// run likelihood calculation on GPU (includes memory transfers)
	RunGPUCLA(h_Lpr, h_Rpr, h_LCL, h_RCL, h_CLA, d_Lpr, d_Rpr, d_LCL, d_RCL,
			d_CLA, mem_size_pr, mem_size_CL, nstates, numRateCats, nchar,
			ncharGPU, dimBlock, dimGrid);

	// stop GPU timer
	CUT_SAFE_CALL(cutStopTimer(gpu_timer));

	// create and start CPU timer
	unsigned int cpu_timer = 0;
	CUT_SAFE_CALL(cutCreateTimer(&cpu_timer));
	CUT_SAFE_CALL(cutStartTimer(cpu_timer));

	// run likelihood calculation on CPU as a reference solution
	ComputeGarliCLA(reference_CLA, h_Lpr, h_Rpr, h_LCL, h_RCL);

	// stop CPU timer
	CUT_SAFE_CALL(cutStopTimer(cpu_timer));

	// check result
	CUTBoolean res = cutCompareL2fe((const float*) reference_CLA,
			(const float*) h_CLA, size_CL_host, 1e-6f);
	printf("%s\t", (1 == res) ? "PASSED" : "FAILED");

	// free memory cuda allocated pinned memory on host
	FreePinnedMemory(&h_Lpr);
	FreePinnedMemory(&h_Rpr);
	FreePinnedMemory(&h_LCL);
	FreePinnedMemory(&h_RCL);
	FreePinnedMemory(&h_CLA);

	// free gpu memory
	FreeGPU(d_Lpr);
	FreeGPU(d_Rpr);
	FreeGPU(d_LCL);
	FreeGPU(d_RCL);
	FreeGPU(d_CLA);

	free(reference_CLA);

	printf("pinned\t");
	printf("%i bytes\t", sizeof(FLOAT_TYPE));
	printf("%i\t", nstates);
	printf("%i\t", numRateCats);
	printf("%i\t", nchar);
	printf("%4.2f\t", cutGetTimerValue(cpu_timer));
	printf("%4.4f\t", cutGetTimerValue(gpu_timer));
	printf("%4.4f\n", cutGetTimerValue(cpu_timer) / cutGetTimerValue(gpu_timer));

	CUT_SAFE_CALL(cutDeleteTimer(gpu_timer));
	CUT_SAFE_CALL(cutDeleteTimer(cpu_timer));

}
