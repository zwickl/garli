/*
 * cudaman.h
 *
 *  Created on: Nov 11, 2008
 *      Author: ayres
 */

#ifndef CUDAMANAGER
#define CUDAMANAGER

using namespace std;

#include "defs.h"
#include "vector_types.h"

class CudaManager {

public:
	CudaManager();
	CudaManager(int nstatesIn, int numRateCatsIn, int ncharIn);
	~CudaManager();

	void BenchmarkGPU();

private:
	int nstates;
	int numRateCats;
	int nchar;

private:
	// Allocates a matrix with random FLOAT_TYPE entries.
	void RandomInit(FLOAT_TYPE* data, int size);

	// Allocates a matrix with random integer entries.
	void RandomIntInit(int* data, int size, int max);

	// Computes reference CLA result
	void ComputeRefCLA(FLOAT_TYPE* dest, const FLOAT_TYPE* Lpr,
					const FLOAT_TYPE* Rpr, const FLOAT_TYPE* LCL,
					const FLOAT_TYPE* RCL);

	// Computes reference deriv result
	void ComputeRefDeriv(const FLOAT_TYPE *partial,
			const FLOAT_TYPE *CL1, const int *partial_underflow_mult,
			const int *CL1_underflow_mult, const FLOAT_TYPE *prmat,
			const FLOAT_TYPE *d1mat, const FLOAT_TYPE *d2mat,
			const FLOAT_TYPE *rateProb, const FLOAT_TYPE *freqs,
			const int *countit, const int lastConst, const bool NoPinvInModel,
			const FLOAT_TYPE prI, const int *conStates, FLOAT_TYPE *reference_Tots);

	void BenchmarkGPUCLA();

	void BenchmarkGPUDeriv();
};

////////////////////////////////////////////////////////////////////////////////
// declaration, forward
extern "C" {
void AllocatePinnedMemory(void** arr, unsigned int mem_size_bytes);

void AllocateGPU(void** arr, unsigned int mem_size_bytes);

void FreeGPU(void* arr);

void FreePinnedMemory(void* arr);

void ComputeGPUCLA(FLOAT_TYPE* h_Lpr, FLOAT_TYPE* h_Rpr, FLOAT_TYPE* h_LCL,
		FLOAT_TYPE* h_RCL, FLOAT_TYPE* h_CLA, FLOAT_TYPE* d_Lpr,
		FLOAT_TYPE* d_Rpr, FLOAT_TYPE* d_LCL, FLOAT_TYPE* d_RCL,
		FLOAT_TYPE* d_CLA, unsigned int mem_size_pr, unsigned int mem_size_CL,
		int nstates, int nRateCats, int nchar, int ncharGPU, dim3 dimBlock,
		dim3 dimGrid);

void ComputeGPUDeriv(FLOAT_TYPE* h_partial, FLOAT_TYPE* h_CL1,
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


#endif

