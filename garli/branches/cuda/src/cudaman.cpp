/*
 * cudaman.cpp
 *
 *  Created on: Nov 12, 2008
 *      Author: ayres
 */

using namespace std;

#include <math.h>

#include "cudaman.h"
#include "defs.h"
#include "outputman.h"
#include "cutil.h"
#include "vector_types.h"

extern OutputManager outman;

CudaManager::CudaManager() {
	gpu_cla_enabled = gpu_deriv_enabled = false;
}

CudaManager::CudaManager(int nstates_in, int numRateCats_in, int nchar_in) {
	gpu_cla_enabled = gpu_deriv_enabled = true;
	nstates = nstates_in;
	numRateCats = numRateCats_in;
	nchar = nchar_in;
	test_iterations = 0;
	print_to_screen = false;

	SetGPUCLAParameters();
	SetGPUDerivParameters();
}

CudaManager::CudaManager(int nstates_in, int numRateCats_in, int nchar_in,
		int test_iterations_in, bool print_to_screen_in) {
	gpu_cla_enabled = gpu_deriv_enabled = true;
	nstates = nstates_in;
	numRateCats = numRateCats_in;
	nchar = nchar_in;
	test_iterations = test_iterations_in;
	print_to_screen = print_to_screen_in;

	if (test_iterations > 0)
		TestGPU();

	if (gpu_cla_enabled)
		SetGPUCLAParameters();

	if (gpu_deriv_enabled)
		SetGPUDerivParameters();
}

CudaManager::~CudaManager() {
	if (gpu_cla_enabled) {
		FreeGPU(cla_d_Lpr);
		FreeGPU(cla_d_Rpr);
		FreeGPU(cla_d_LCL);
		FreeGPU(cla_d_RCL);
		FreeGPU(cla_d_CLA);
	}

	if (gpu_deriv_enabled) {
		FreeGPU(deriv_d_countit);
		FreeGPU(deriv_d_partial);
		FreeGPU(deriv_d_CL1);
		FreeGPU(deriv_d_partial_underflow_mult);
		FreeGPU(deriv_d_CL1_underflow_mult);
		FreeGPU(deriv_d_prmat);
		FreeGPU(deriv_d_d1mat);
		FreeGPU(deriv_d_d2mat);
		FreeGPU(deriv_d_rateProb);
		FreeGPU(deriv_d_freqs);
		FreeGPU(deriv_d_Tots);
		FreeGPU(deriv_d_Tots_arr);
		FreeGPU(deriv_d_conStates);

#ifdef CUDA_MEMORY_PINNED
		FreePinnedMemory(&deriv_h_Tots);
		FreePinnedMemory(&deriv_h_Tots_arr);
#else
		free(deriv_h_Tots);
		free(deriv_h_Tots_arr);
#endif
	}
}

bool CudaManager::GetGPUCLAEnabled() {
	return gpu_cla_enabled;
}

bool CudaManager::GetGPUDerivEnabled() {
	return gpu_deriv_enabled;
}

void CudaManager::ComputeGPUCLA(const FLOAT_TYPE* h_Lpr,
		const FLOAT_TYPE* h_Rpr, const FLOAT_TYPE* h_LCL,
		const FLOAT_TYPE* h_RCL, const FLOAT_TYPE* h_CLA) {
	CuComputeGPUCLA(h_Lpr, h_Rpr, h_LCL, h_RCL, h_CLA, cla_d_Lpr, cla_d_Rpr,
			cla_d_LCL, cla_d_RCL, cla_d_CLA, cla_mem_size_pr, cla_mem_size_CL,
			nstates, numRateCats, nchar, cla_ncharGPU, cla_dimBlock,
			cla_dimGrid);
}

void CudaManager::ComputeGPUDeriv(FLOAT_TYPE* h_partial, FLOAT_TYPE* h_CL1,
		int* h_partial_underflow_mult, int* h_CL1_underflow_mult,
		const FLOAT_TYPE* h_prmat, const FLOAT_TYPE* h_d1mat,
		const FLOAT_TYPE* h_d2mat, const FLOAT_TYPE* h_rateProb,
		FLOAT_TYPE* h_freqs, const int* h_countit, const int* h_conStates,
		int lastConst, bool NoPinvInModel, FLOAT_TYPE prI) {
	CuComputeGPUDeriv(h_partial, h_CL1, h_partial_underflow_mult,
			h_CL1_underflow_mult, h_prmat, h_d1mat, h_d2mat, h_rateProb,
			h_freqs, h_countit, h_conStates, deriv_h_Tots, deriv_h_Tots_arr,
			deriv_d_partial, deriv_d_CL1, deriv_d_partial_underflow_mult,
			deriv_d_CL1_underflow_mult, deriv_d_prmat, deriv_d_d1mat,
			deriv_d_d2mat, deriv_d_rateProb, deriv_d_freqs, deriv_d_countit,
			deriv_d_conStates, deriv_d_Tots, deriv_d_Tots_arr,
			deriv_mem_size_pr, deriv_mem_size_CL, deriv_mem_size_int_char,
			deriv_mem_size_rates, deriv_mem_size_states, deriv_mem_size_Tots,
			deriv_mem_size_Tots_arr, lastConst, NoPinvInModel, prI, nstates,
			numRateCats, nchar, deriv_ncharGPU, deriv_dimBlock, deriv_dimGrid);
}

FLOAT_TYPE CudaManager::GetDerivTots(int index) {
	return deriv_h_Tots[index];
}

// Private methods
void CudaManager::TestGPU() {
	if (print_to_screen) {
		outman.UserMessageNoCR(
				"======================= GPU Tests =======================");

		outman.UserMessageNoCR("\nParameters:\n");
		outman.UserMessageNoCR("memory \t");
		outman.UserMessageNoCR("float \t");
		outman.UserMessageNoCR("states\t");
		outman.UserMessageNoCR("rates\t");
		outman.UserMessageNoCR("chars\n");
#ifdef CUDA_MEMORY_PINNED
		outman.UserMessageNoCR("pinned\t");
#else
		outman.UserMessageNoCR("paged\t");
#endif
		outman.UserMessageNoCR("%i bytes\t", sizeof(FLOAT_TYPE));
		outman.UserMessageNoCR("%i\t", nstates);
		outman.UserMessageNoCR("%i\t", numRateCats);
		outman.UserMessageNoCR("%i\n", nchar);

		outman.UserMessageNoCR(
				"\nBenchmarks: (best results out of %i iterations)\n",
				test_iterations);
		outman.UserMessageNoCR("function\t");
		outman.UserMessageNoCR("speedup\t");
		outman.UserMessageNoCR("cpu   \t");
		outman.UserMessageNoCR("gpu   \t");
		outman.UserMessageNoCR("test\t");
		outman.UserMessageNoCR("enabled\n");
	}

	TestGPUCLA();
	TestGPUDeriv();

	if (print_to_screen)
		outman.UserMessageNoCR(
				"=========================================================\n\n");
}

void CudaManager::SetGPUCLAParameters() {
	// setup block and grid size for gpu
	if (nstates == 4) {
		cla_dimBlock = dim3(64 * numRateCats);
		cla_dimGrid = dim3(nchar / 16);
		cla_ncharGPU = 16 * cla_dimGrid.x;
	} else if (nstates == 20) {
		cla_dimBlock = dim3(20, 20);
		cla_dimGrid = dim3(nchar / 20);
		cla_ncharGPU = 20 * cla_dimGrid.x;
	} else {
		cla_dimBlock = dim3(16, 16);
		cla_dimGrid = dim3(nchar / 16, 4);
		cla_ncharGPU = 16 * cla_dimGrid.x;
	}

	int cla_size_pr = nstates * nstates * numRateCats;
	int cla_size_CL = nstates * cla_ncharGPU * numRateCats;
	cla_mem_size_pr = sizeof(FLOAT_TYPE) * cla_size_pr;
	cla_mem_size_CL = sizeof(FLOAT_TYPE) * cla_size_CL;

	if (nstates == 61) {
		unsigned int cla_mem_size_pr_64 = sizeof(FLOAT_TYPE) * 64 * 64
				* numRateCats;
		unsigned int cla_mem_size_CL_64 = sizeof(FLOAT_TYPE) * 64
				* cla_ncharGPU * numRateCats;

		AllocateGPU((void**) &cla_d_Lpr, cla_mem_size_pr_64);
		AllocateGPU((void**) &cla_d_Rpr, cla_mem_size_pr_64);
		AllocateGPU((void**) &cla_d_LCL, cla_mem_size_CL_64);
		AllocateGPU((void**) &cla_d_RCL, cla_mem_size_CL_64);
		AllocateGPU((void**) &cla_d_CLA, cla_mem_size_CL_64);
	} else {
		AllocateGPU((void**) &cla_d_Lpr, cla_mem_size_pr);
		AllocateGPU((void**) &cla_d_Rpr, cla_mem_size_pr);
		AllocateGPU((void**) &cla_d_LCL, cla_mem_size_CL);
		AllocateGPU((void**) &cla_d_RCL, cla_mem_size_CL);
		AllocateGPU((void**) &cla_d_CLA, cla_mem_size_CL);
	}
}

void CudaManager::SetGPUDerivParameters() {
	int deriv_nBlock = 128;
	int deriv_aaBlock = 128;
	int deriv_cBlock = 128;

	if (nstates == 4) {
		deriv_dimBlock = dim3(deriv_nBlock);
		deriv_ncharGPU = deriv_nBlock * ((int) nchar / deriv_nBlock);
		deriv_dimGrid = dim3(deriv_ncharGPU / deriv_nBlock);
	} else if (nstates == 20) {
		deriv_dimBlock = dim3(deriv_aaBlock);
		deriv_ncharGPU = deriv_aaBlock * ((int) nchar / deriv_aaBlock);
		deriv_dimGrid = dim3(deriv_ncharGPU / deriv_aaBlock);
	} else {
		deriv_dimBlock = dim3(deriv_cBlock);
		deriv_ncharGPU = deriv_cBlock * ((int) nchar / deriv_cBlock);
		deriv_dimGrid = dim3(deriv_ncharGPU / deriv_cBlock);
	}

	int deriv_size_pr = nstates * nstates * numRateCats;
	int deriv_size_CL = nstates * deriv_ncharGPU * numRateCats;
	int deriv_size_CL_host = nstates * nchar * numRateCats;
	int deriv_size_Tots = 3;
	int deriv_size_Tots_arr = 3 * deriv_dimGrid.x;
	deriv_mem_size_pr = sizeof(FLOAT_TYPE) * deriv_size_pr;
	deriv_mem_size_CL = sizeof(FLOAT_TYPE) * deriv_size_CL;
	deriv_mem_size_int_char = sizeof(int) * deriv_ncharGPU;
	deriv_mem_size_rates = sizeof(int) * numRateCats;
	deriv_mem_size_states = sizeof(int) * nstates;
	deriv_mem_size_Tots = sizeof(FLOAT_TYPE) * deriv_size_Tots;
	deriv_mem_size_Tots_arr = sizeof(FLOAT_TYPE) * deriv_size_Tots_arr;

#ifdef CUDA_MEMORY_PINNED
	AllocatePinnedMemory((void**) &deriv_h_Tots, deriv_mem_size_Tots);
	AllocatePinnedMemory((void**) &deriv_h_Tots_arr, deriv_mem_size_Tots_arr);
#else
	deriv_h_Tots = (FLOAT_TYPE*) malloc(deriv_mem_size_Tots);
	deriv_h_Tots_arr = (FLOAT_TYPE*) malloc(deriv_mem_size_Tots_arr);
#endif

	// allocate memory on gpu
	AllocateGPU((void**) &deriv_d_countit, deriv_mem_size_int_char);
	AllocateGPU((void**) &deriv_d_partial, deriv_mem_size_CL);
	AllocateGPU((void**) &deriv_d_CL1, deriv_mem_size_CL);
	AllocateGPU((void**) &deriv_d_partial_underflow_mult,
			deriv_mem_size_int_char);
	AllocateGPU((void**) &deriv_d_CL1_underflow_mult, deriv_mem_size_int_char);
	AllocateGPU((void**) &deriv_d_prmat, deriv_mem_size_pr);
	AllocateGPU((void**) &deriv_d_d1mat, deriv_mem_size_pr);
	AllocateGPU((void**) &deriv_d_d2mat, deriv_mem_size_pr);
	AllocateGPU((void**) &deriv_d_rateProb, deriv_mem_size_rates);
	AllocateGPU((void**) &deriv_d_freqs, deriv_mem_size_states);
	AllocateGPU((void**) &deriv_d_Tots, deriv_mem_size_Tots);
	AllocateGPU((void**) &deriv_d_Tots_arr, deriv_mem_size_Tots_arr);
	AllocateGPU((void**) &deriv_d_conStates, deriv_mem_size_int_char);
}

void CudaManager::RandomInit(FLOAT_TYPE* data, int size) {
	for (int i = 0; i < size; ++i)
		data[i] = rand() / (FLOAT_TYPE) RAND_MAX;
}

void CudaManager::RandomIntInit(int* data, int size, int max) {
	for (int i = 0; i < size; ++i)
		data[i] = 1 + ((int) (((double) rand() / ((double) (RAND_MAX)
				+ (double) (1))) * max));
}

void CudaManager::ComputeRefCLA(FLOAT_TYPE* dest, const FLOAT_TYPE* Lpr,
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

void CudaManager::ComputeRefDeriv(const FLOAT_TYPE *partial,
		const FLOAT_TYPE *CL1, const int *partial_underflow_mult,
		const int *CL1_underflow_mult, const FLOAT_TYPE *prmat,
		const FLOAT_TYPE *d1mat, const FLOAT_TYPE *d2mat,
		const FLOAT_TYPE *rateProb, const FLOAT_TYPE *freqs,
		const int *countit, const int lastConst, const bool NoPinvInModel,
		const FLOAT_TYPE prI, const int *conStates, FLOAT_TYPE *reference_Tots) {
	FLOAT_TYPE tot1 = 0, tot2 = 0, totL = 0;

	FLOAT_TYPE siteL, siteD1, siteD2;
	FLOAT_TYPE tempL, tempD1, tempD2;
	FLOAT_TYPE rateL, rateD1, rateD2;

	for (int i = 0; i < nchar; i++) {
		siteL = siteD1 = siteD2 = 0;
		for (int rate = 0; rate < numRateCats; rate++) {
			rateL = rateD1 = rateD2 = 0;
			int rateOffset = rate * nstates * nstates;
			for (int from = 0; from < nstates; from++) {
				tempL = tempD1 = tempD2 = 0;
				int offset = from * nstates;
				for (int to = 0; to < nstates; to++) {
					tempL += prmat[rateOffset + offset + to] * CL1[to];
					tempD1 += d1mat[rateOffset + offset + to] * CL1[to];
					tempD2 += d2mat[rateOffset + offset + to] * CL1[to];
				}
				rateL += tempL * partial[from] * freqs[from];
				rateD1 += tempD1 * partial[from] * freqs[from];
				rateD2 += tempD2 * partial[from] * freqs[from];
			}
			siteL += rateL * rateProb[rate];
			siteD1 += rateD1 * rateProb[rate];
			siteD2 += rateD2 * rateProb[rate];
			partial += nstates;
			CL1 += nstates;
		}

		if ((NoPinvInModel == false) && (i <= lastConst)) {
			if (nstates == 4) {
				float btot = 0.0f;
				if (conStates[i] & 1)
					btot += freqs[0];
				if (conStates[i] & 2)
					btot += freqs[1];
				if (conStates[i] & 4)
					btot += freqs[2];
				if (conStates[i] & 8)
					btot += freqs[3];
				siteL += (prI * btot) * exp(partial_underflow_mult[i]
						+ CL1_underflow_mult[i]);
			} else
				siteL += (prI * freqs[conStates[i]] * exp(
						(FLOAT_TYPE) partial_underflow_mult[i]) * exp(
						(FLOAT_TYPE) CL1_underflow_mult[i]));
		}

		totL += countit[i] * (log(siteL) - partial_underflow_mult[i]
				- CL1_underflow_mult[i]);
		siteD1 /= siteL;
		tot1 += countit[i] * siteD1;
		tot2 += countit[i] * ((siteD2 / siteL) - siteD1 * siteD1);
	}

	reference_Tots[0] = totL;
	reference_Tots[1] = tot1;
	reference_Tots[2] = tot2;
}

void CudaManager::TestGPUCLA() {
	FLOAT_TYPE *h_Lpr, *h_Rpr, *h_LCL, *h_RCL, *h_CLA, *reference_CLA;
	FLOAT_TYPE * d_Lpr, *d_Rpr, *d_LCL, *d_RCL, *d_CLA;

	FLOAT_TYPE best_cpu = 0;
	FLOAT_TYPE best_gpu = 0;
	bool test_passed = true;

	for (int i = 0; i < test_iterations; i++) {
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
#ifdef CUDA_MEMORY_PINNED
		AllocatePinnedMemory((void**) &h_Lpr, mem_size_pr);
		AllocatePinnedMemory((void**) &h_Rpr, mem_size_pr);
		AllocatePinnedMemory((void**) &h_LCL, mem_size_CL_host);
		AllocatePinnedMemory((void**) &h_RCL, mem_size_CL_host);
		AllocatePinnedMemory((void**) &h_CLA, mem_size_CL_host);
#else
		h_Lpr = (FLOAT_TYPE*) malloc(mem_size_pr);
		h_Rpr = (FLOAT_TYPE*) malloc(mem_size_pr);
		h_LCL = (FLOAT_TYPE*) malloc(mem_size_CL_host);
		h_RCL = (FLOAT_TYPE*) malloc(mem_size_CL_host);
		h_CLA = (FLOAT_TYPE*) malloc(mem_size_CL_host);
#endif

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
		srand(i);
		RandomInit(h_Lpr, size_pr);
		RandomInit(h_Rpr, size_pr);
		RandomInit(h_LCL, size_CL_host);
		RandomInit(h_RCL, size_CL_host);

		// create and start GPU timer
		unsigned int gpu_timer = 0;
		CUT_SAFE_CALL(cutCreateTimer(&gpu_timer));
		CUT_SAFE_CALL(cutStartTimer(gpu_timer));

		// run likelihood calculation on GPU (includes memory transfers)
		CuComputeGPUCLA(h_Lpr, h_Rpr, h_LCL, h_RCL, h_CLA, d_Lpr, d_Rpr, d_LCL,
				d_RCL, d_CLA, mem_size_pr, mem_size_CL, nstates, numRateCats,
				nchar, ncharGPU, dimBlock, dimGrid);

		// stop GPU timer
		CUT_SAFE_CALL(cutStopTimer(gpu_timer));

		// create and start CPU timer
		unsigned int cpu_timer = 0;
		CUT_SAFE_CALL(cutCreateTimer(&cpu_timer));
		CUT_SAFE_CALL(cutStartTimer(cpu_timer));

		// run likelihood calculation on CPU as a reference solution
		ComputeRefCLA(reference_CLA, h_Lpr, h_Rpr, h_LCL, h_RCL);

		// stop CPU timer
		CUT_SAFE_CALL(cutStopTimer(cpu_timer));

		if (i == 0 | cutGetTimerValue(cpu_timer) < best_cpu)
			best_cpu = cutGetTimerValue(cpu_timer);

		if (i == 0 | cutGetTimerValue(gpu_timer) < best_gpu)
			best_gpu = cutGetTimerValue(gpu_timer);

		CUT_SAFE_CALL(cutDeleteTimer(gpu_timer));
		CUT_SAFE_CALL(cutDeleteTimer(cpu_timer));

		// check result
		CUTBoolean res = cutCompareL2fe((const float*) reference_CLA,
				(const float*) h_CLA, size_CL_host, 1e-6f);

		if (res == false)
			test_passed = false;

		// free memory cuda allocated pinned memory on host
#ifdef CUDA_MEMORY_PINNED
		FreePinnedMemory(&h_Lpr);
		FreePinnedMemory(&h_Rpr);
		FreePinnedMemory(&h_LCL);
		FreePinnedMemory(&h_RCL);
		FreePinnedMemory(&h_CLA);
#else
		free(h_Lpr);
		free(h_Rpr);
		free(h_LCL);
		free(h_RCL);
		free(h_CLA);
#endif

		// free gpu memory
		FreeGPU(d_Lpr);
		FreeGPU(d_Rpr);
		FreeGPU(d_LCL);
		FreeGPU(d_RCL);
		FreeGPU(d_CLA);

		free(reference_CLA);
	}

	if (!test_passed || best_cpu / best_gpu < 1)
		gpu_cla_enabled = false;

	if (print_to_screen) {
		outman.UserMessageNoCR("GPUCLA   \t");
		outman.UserMessageNoCR("%4.4f\t", best_cpu / best_gpu);
		outman.UserMessageNoCR("%4.2f\t", best_cpu);
		outman.UserMessageNoCR("%4.4f\t", best_gpu);
		outman.UserMessageNoCR("%s\t", (test_passed) ? "passed" : "failed");
		outman.UserMessageNoCR("%s\n", (gpu_cla_enabled) ? "yes" : "no");
	}
}

void CudaManager::TestGPUDeriv() {
	int *h_countit, *h_partial_underflow_mult, *h_CL1_underflow_mult,
			*h_conStates;

	FLOAT_TYPE best_cpu = 0;
	FLOAT_TYPE best_gpu = 0;
	bool test_passed = true;

	for (int i = 0; i < test_iterations; i++) {

		FLOAT_TYPE *h_partial, *h_CL1, *h_prmat, *h_d1mat, *h_d2mat,
				*h_rateProb, *h_freqs;
		FLOAT_TYPE *h_Tots_arr, *h_Tots;
		FLOAT_TYPE *reference_Tots;

		int lastConst;
		FLOAT_TYPE prI;

		bool NoPinvInModel = false;

		int *d_countit, *d_partial_underflow_mult, *d_CL1_underflow_mult,
				*d_conStates;
		FLOAT_TYPE *d_partial, *d_CL1, *d_prmat, *d_d1mat, *d_d2mat,
				*d_rateProb, *d_freqs;
		FLOAT_TYPE *d_Tots_arr, *d_Tots;

		int nBlock = 128;
		int aaBlock = 128;
		int cBlock = 128;

		int ncharGPU;
		if (nstates == 4)
			ncharGPU = nBlock * ((int) nchar / nBlock);
		else if (nstates == 20)
			ncharGPU = aaBlock * ((int) nchar / aaBlock);
		else
			ncharGPU = cBlock * ((int) nchar / cBlock);

		int size_pr = nstates * nstates * numRateCats;
		int size_CL = nstates * ncharGPU * numRateCats;
		int size_CL_host = nstates * nchar * numRateCats;
		int size_Tots = 3;
		unsigned int mem_size_pr = sizeof(FLOAT_TYPE) * size_pr;
		unsigned int mem_size_CL = sizeof(FLOAT_TYPE) * size_CL;
		unsigned int mem_size_CL_host = sizeof(FLOAT_TYPE) * size_CL_host;
		unsigned int mem_size_int_char = sizeof(int) * ncharGPU;
		unsigned int mem_size_int_char_host = sizeof(int) * nchar;
		unsigned int mem_size_rates = sizeof(int) * numRateCats;
		unsigned int mem_size_states = sizeof(int) * nstates;
		unsigned int mem_size_Tots = sizeof(FLOAT_TYPE) * size_Tots;

		// setup block and grid size for gpu
		dim3 dimBlock;
		dim3 dimGrid;

		if (nstates == 4) {
			dimBlock = dim3(nBlock);
			dimGrid = dim3(ncharGPU / nBlock);
		} else if (nstates == 20) {
			dimBlock = dim3(aaBlock);
			dimGrid = dim3(ncharGPU / aaBlock);
		} else {
			dimBlock = dim3(cBlock);
			dimGrid = dim3(ncharGPU / cBlock);
		}

		int size_Tots_arr = 3 * dimGrid.x;
		unsigned int mem_size_Tots_arr = sizeof(FLOAT_TYPE) * size_Tots_arr;

		// allocate host memory
#ifdef CUDA_MEMORY_PINNED
		AllocatePinnedMemory((void**) &h_countit, mem_size_int_char_host);
		AllocatePinnedMemory((void**) &h_partial, mem_size_CL_host);
		AllocatePinnedMemory((void**) &h_CL1, mem_size_CL_host);
		AllocatePinnedMemory((void**) &h_partial_underflow_mult,
				mem_size_int_char_host);
		AllocatePinnedMemory((void**) &h_CL1_underflow_mult,
				mem_size_int_char_host);
		AllocatePinnedMemory((void**) &h_prmat, mem_size_pr);
		AllocatePinnedMemory((void**) &h_d1mat, mem_size_pr);
		AllocatePinnedMemory((void**) &h_d2mat, mem_size_pr);
		AllocatePinnedMemory((void**) &h_rateProb, mem_size_rates);
		AllocatePinnedMemory((void**) &h_freqs, mem_size_states);
		AllocatePinnedMemory((void**) &h_Tots, mem_size_Tots);
		AllocatePinnedMemory((void**) &h_Tots_arr, mem_size_Tots_arr);
		AllocatePinnedMemory((void**) &h_conStates, mem_size_int_char_host);
#else
		h_countit = (int*) malloc(mem_size_int_char_host);
		h_partial = (FLOAT_TYPE*) malloc(mem_size_CL_host);
		h_CL1 = (FLOAT_TYPE*) malloc(mem_size_CL_host);
		h_partial_underflow_mult = (int*) malloc(mem_size_int_char_host);
		h_CL1_underflow_mult = (int*) malloc(mem_size_int_char_host);
		h_prmat = (FLOAT_TYPE*) malloc(mem_size_pr);
		h_d1mat = (FLOAT_TYPE*) malloc(mem_size_pr);
		h_d2mat = (FLOAT_TYPE*) malloc(mem_size_pr);
		h_rateProb = (FLOAT_TYPE*) malloc(mem_size_rates);
		h_freqs = (FLOAT_TYPE*) malloc(mem_size_states);
		h_Tots = (FLOAT_TYPE*) malloc(mem_size_Tots);
		h_Tots_arr = (FLOAT_TYPE*) malloc(mem_size_Tots_arr);
		h_conStates = (int*) malloc(mem_size_int_char_host);
#endif

		reference_Tots = (FLOAT_TYPE*) malloc(mem_size_Tots);

		// allocate memory on gpu
		AllocateGPU((void**) &d_countit, mem_size_int_char);
		AllocateGPU((void**) &d_partial, mem_size_CL);
		AllocateGPU((void**) &d_CL1, mem_size_CL);
		AllocateGPU((void**) &d_partial_underflow_mult, mem_size_int_char);
		AllocateGPU((void**) &d_CL1_underflow_mult, mem_size_int_char);
		AllocateGPU((void**) &d_prmat, mem_size_pr);
		AllocateGPU((void**) &d_d1mat, mem_size_pr);
		AllocateGPU((void**) &d_d2mat, mem_size_pr);
		AllocateGPU((void**) &d_rateProb, mem_size_rates);
		AllocateGPU((void**) &d_freqs, mem_size_states);
		AllocateGPU((void**) &d_Tots, mem_size_Tots);
		AllocateGPU((void**) &d_Tots_arr, mem_size_Tots_arr);
		AllocateGPU((void**) &d_conStates, mem_size_int_char);

		// initialize host memory with random data
		srand(i);
		RandomIntInit(&lastConst, 1, nchar);
		RandomInit(&prI, 1);
		RandomIntInit(h_countit, nchar, 20);
		RandomInit(h_partial, size_CL_host);
		RandomInit(h_CL1, size_CL_host);
		RandomIntInit(h_partial_underflow_mult, nchar, 20);
		RandomIntInit(h_CL1_underflow_mult, nchar, 20);
		RandomInit(h_prmat, size_pr);
		RandomInit(h_d1mat, size_pr);
		RandomInit(h_d2mat, size_pr);
		RandomInit(h_rateProb, numRateCats);
		RandomInit(h_freqs, nstates);
		RandomIntInit(h_conStates, nchar, nstates - 1);
		h_Tots[0] = h_Tots[1] = h_Tots[2] = 0;

		reference_Tots[0] = reference_Tots[1] = reference_Tots[2] = 0;

		// create and start GPU timer
		unsigned int gpu_timer = 0;
		CUT_SAFE_CALL(cutCreateTimer(&gpu_timer));
		CUT_SAFE_CALL(cutStartTimer(gpu_timer));

		// run likelihood calculation on GPU (includes memory transfers)
		CuComputeGPUDeriv(h_partial, h_CL1, h_partial_underflow_mult,
				h_CL1_underflow_mult, h_prmat, h_d1mat, h_d2mat, h_rateProb,
				h_freqs, h_countit, h_conStates, h_Tots, h_Tots_arr, d_partial,
				d_CL1, d_partial_underflow_mult, d_CL1_underflow_mult, d_prmat,
				d_d1mat, d_d2mat, d_rateProb, d_freqs, d_countit, d_conStates,
				d_Tots, d_Tots_arr, mem_size_pr, mem_size_CL,
				mem_size_int_char, mem_size_rates, mem_size_states,
				mem_size_Tots, mem_size_Tots_arr, lastConst, NoPinvInModel,
				prI, nstates, numRateCats, nchar, ncharGPU, dimBlock, dimGrid);

		// stop GPU timer
		CUT_SAFE_CALL(cutStopTimer(gpu_timer));

		// create and start CPU timer
		unsigned int cpu_timer = 0;
		CUT_SAFE_CALL(cutCreateTimer(&cpu_timer));
		CUT_SAFE_CALL(cutStartTimer(cpu_timer));

		// run likelihood calculation on CPU as a reference solution
		ComputeRefDeriv(h_partial, h_CL1, h_partial_underflow_mult,
				h_CL1_underflow_mult, h_prmat, h_d1mat, h_d2mat, h_rateProb,
				h_freqs, h_countit, lastConst, NoPinvInModel, prI, h_conStates,
				reference_Tots);

		// stop CPU timer
		CUT_SAFE_CALL(cutStopTimer(cpu_timer));

		if (i == 0 | cutGetTimerValue(cpu_timer) < best_cpu)
			best_cpu = cutGetTimerValue(cpu_timer);

		if (i == 0 | cutGetTimerValue(gpu_timer) < best_gpu)
			best_gpu = cutGetTimerValue(gpu_timer);

		CUT_SAFE_CALL(cutDeleteTimer(gpu_timer));
		CUT_SAFE_CALL(cutDeleteTimer(cpu_timer));

		// check result
		CUTBoolean res = cutCompareL2fe(reference_Tots, h_Tots, size_Tots,
				1e-5f);

		if (res == false)
			test_passed = false;

		// free memory cuda allocated memory on host
#ifdef CUDA_MEMORY_PINNED
		FreePinnedMemory(&h_countit);
		FreePinnedMemory(&h_partial);
		FreePinnedMemory(&h_CL1);
		FreePinnedMemory(&h_partial_underflow_mult);
		FreePinnedMemory(&h_CL1_underflow_mult);
		FreePinnedMemory(&h_prmat);
		FreePinnedMemory(&h_d1mat);
		FreePinnedMemory(&h_d2mat);
		FreePinnedMemory(&h_rateProb);
		FreePinnedMemory(&h_freqs);
		FreePinnedMemory(&h_Tots);
		FreePinnedMemory(&h_Tots_arr);
		FreePinnedMemory(&h_conStates);
#else
		free(h_countit);
		free(h_partial);
		free(h_CL1);
		free(h_partial_underflow_mult);
		free(h_CL1_underflow_mult);
		free(h_prmat);
		free(h_d1mat);
		free(h_d2mat);
		free(h_rateProb);
		free(h_freqs);
		free(h_Tots);
		free(h_Tots_arr);
		free(h_conStates);
#endif

		free(reference_Tots);

		// free gpu memory
		FreeGPU(d_countit);
		FreeGPU(d_partial);
		FreeGPU(d_CL1);
		FreeGPU(d_partial_underflow_mult);
		FreeGPU(d_CL1_underflow_mult);
		FreeGPU(d_prmat);
		FreeGPU(d_d1mat);
		FreeGPU(d_d2mat);
		FreeGPU(d_rateProb);
		FreeGPU(d_freqs);
		FreeGPU(d_Tots);
		FreeGPU(d_Tots_arr);
		FreeGPU(d_conStates);
	}

	if (!test_passed || best_cpu / best_gpu < 1)
		gpu_deriv_enabled = false;

	if (print_to_screen) {
		outman.UserMessageNoCR("GPUDeriv   \t");
		outman.UserMessageNoCR("%4.4f\t", best_cpu / best_gpu);
		outman.UserMessageNoCR("%4.2f\t", best_cpu);
		outman.UserMessageNoCR("%4.4f\t", best_gpu);
		outman.UserMessageNoCR("%s\t", (test_passed) ? "passed" : "failed");
		outman.UserMessageNoCR("%s\n", (gpu_deriv_enabled) ? "yes" : "no");
	}
}
