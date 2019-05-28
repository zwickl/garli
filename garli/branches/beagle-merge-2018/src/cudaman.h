/*
 * cudaman.h
 *
 *  Created on: Nov 11, 2008
 *      Author: ayres
 */

#ifndef CUDAMANAGER
#define CUDAMANAGERe


#include <vector_types.h>

#include "defs.h"

class CudaManager {

public:
	CudaManager();
	CudaManager(int nstates_in, int numRateCats_in, int nchar_in);
	CudaManager(int nstates_in, int numRateCats_in, int nchar_in, int device_number,
			int test_iterations_in, bool print_to_screen_in, bool print_device_query);
	~CudaManager();

	void Initialization();

	void ChangeNChar(int nchar_boot_in, const int* counts_in);

	bool GetPinnedMemoryEnabled();

	bool GetGPUCLAEnabled();
	bool GetGPUDerivEnabled();

	void ComputeGPUCLA(const FLOAT_TYPE* h_Lpr, const FLOAT_TYPE* h_Rpr,
			const FLOAT_TYPE* h_LCL, const FLOAT_TYPE* h_RCL,
			const FLOAT_TYPE* h_CLA);

	void ComputeGPUDeriv(const FLOAT_TYPE* h_partial, const FLOAT_TYPE* h_CL1,
			const int* h_partial_underflow_mult, const int* h_CL1_underflow_mult,
			const FLOAT_TYPE* h_prmat, const FLOAT_TYPE* h_d1mat,
			const FLOAT_TYPE* h_d2mat, const FLOAT_TYPE* h_rateProb,
			const FLOAT_TYPE* h_freqs, const int* h_countit, const int* h_conStates,
			int lastConst, bool NoPinvInModel, FLOAT_TYPE prI);

	FLOAT_TYPE GetDerivTots(int index);

private:
	bool pinned_memory_enabled;
	bool gpu_cla_enabled, gpu_deriv_enabled;

	int test_iterations;
	bool print_gpu_tests;
	bool print_device_query;
	int device_number;

	int nstates, numRateCats, nchar;

	FLOAT_TYPE *cla_d_Lpr, *cla_d_Rpr, *cla_d_LCL, *cla_d_RCL, *cla_d_CLA;
	dim3 cla_dimBlock, cla_dimGrid;
	int cla_ncharGPU;
	unsigned int cla_mem_size_pr, cla_mem_size_CL;

	int *deriv_d_countit, *deriv_d_partial_underflow_mult,
			*deriv_d_CL1_underflow_mult, *deriv_d_conStates,
			*deriv_d_nchar_boot_index;
	FLOAT_TYPE *deriv_d_partial, *deriv_d_CL1, *deriv_d_prmat, *deriv_d_d1mat,
			*deriv_d_d2mat, *deriv_d_rateProb, *deriv_d_freqs,
			*deriv_d_Tots_arr, *deriv_d_Tots, *deriv_h_Tots_arr, *deriv_h_Tots;
	dim3 deriv_dimBlock, deriv_dimGrid;
	int deriv_ncharGPU;
	int deriv_nchar_full;
	int *deriv_h_nchar_boot_index;
	int *deriv_counts;
	unsigned int deriv_mem_size_pr, deriv_mem_size_CL, deriv_mem_size_int_char,
			deriv_mem_size_rates, deriv_mem_size_states, deriv_mem_size_Tots,
			deriv_mem_size_Tots_arr, deriv_mem_size_nchar_boot_index;

private:
	// Calculate all GPU parameters that don't change from call to call
	void SetGPUCLAParameters();
	void SetGPUDerivParameters();

	// Initialize nchar_full and nchar_boot_index
	void InitGPUDeriv(int nchar_full_in);

	// Free GPU allocated memory
	void FreeGPUCLA();
	void FreeGPUDeriv();

	// Allocates a matrix with random FLOAT_TYPE entries.
	void RandomInit(FLOAT_TYPE* data, int size);

	// Allocates a matrix with random integer entries.
	void RandomIntInit(int* data, int size, int max);

	// Computes reference CLA result
	void
			ComputeRefCLA(FLOAT_TYPE* dest, const FLOAT_TYPE* Lpr,
					const FLOAT_TYPE* Rpr, const FLOAT_TYPE* LCL,
					const FLOAT_TYPE* RCL);

	// Computes reference deriv result
	void ComputeRefDeriv(const FLOAT_TYPE *partial, const FLOAT_TYPE *CL1,
			const int *partial_underflow_mult, const int *CL1_underflow_mult,
			const FLOAT_TYPE *prmat, const FLOAT_TYPE *d1mat,
			const FLOAT_TYPE *d2mat, const FLOAT_TYPE *rateProb,
			const FLOAT_TYPE *freqs, const int *countit, const int lastConst,
			const bool NoPinvInModel, const FLOAT_TYPE prI,
			const int *conStates, FLOAT_TYPE *reference_Tots);

	// Test GPU performance and correctness and enable functions accordingly
	void TestGPU();
	void TestGPUCLA();
	void TestGPUDeriv();
};

////////////////////////////////////////////////////////////////////////////////
// declaration, forward
extern "C" {
bool CheckCuda(unsigned int device_number);

void DeviceQuery();

void SetDevice(unsigned int device_number);

void AllocatePinnedMemory(void** arr, unsigned int mem_size_bytes);

void AllocateGPU(void** arr, unsigned int mem_size_bytes);

void FreeGPU(void* arr);

void FreePinnedMemory(void* arr);

void CuComputeGPUCLA(const FLOAT_TYPE* h_Lpr, const FLOAT_TYPE* h_Rpr,
		const FLOAT_TYPE* h_LCL, const FLOAT_TYPE* h_RCL,
		const FLOAT_TYPE* h_CLA, const FLOAT_TYPE* d_Lpr, FLOAT_TYPE* d_Rpr,
		FLOAT_TYPE* d_LCL, FLOAT_TYPE* d_RCL, FLOAT_TYPE* d_CLA,
		unsigned int mem_size_pr, unsigned int mem_size_CL, int nstates,
		int nRateCats, int nchar, int ncharGPU, dim3 dimBlock, dim3 dimGrid);

void CuComputeGPUDeriv(const FLOAT_TYPE* h_partial, const FLOAT_TYPE* h_CL1,
		const int* h_partial_underflow_mult, const int* h_CL1_underflow_mult,
		const FLOAT_TYPE* h_prmat, const FLOAT_TYPE* h_d1mat,
		const FLOAT_TYPE* h_d2mat, const FLOAT_TYPE* h_rateProb,
		const FLOAT_TYPE* h_freqs, const int* h_countit, const int* h_conStates,
		FLOAT_TYPE* h_Tots, FLOAT_TYPE* h_Tots_arr, int* h_nchar_boot_index,
		FLOAT_TYPE* d_partial, FLOAT_TYPE* d_CL1,
		int* d_partial_underflow_mult, int* d_CL1_underflow_mult,
		FLOAT_TYPE* d_prmat, FLOAT_TYPE* d_d1mat, FLOAT_TYPE* d_d2mat,
		FLOAT_TYPE* d_rateProb, FLOAT_TYPE* d_freqs, int* d_countit,
		int* d_conStates, FLOAT_TYPE* d_Tots, FLOAT_TYPE* d_Tots_arr,
		int* d_nchar_boot_index, unsigned int mem_size_pr,
		unsigned int mem_size_CL, unsigned int mem_size_int_char,
		unsigned int mem_size_rates, unsigned int mem_size_states,
		unsigned int mem_size_Tots, unsigned int mem_size_Tots_arr,
		unsigned int mem_size_nchar_boot_index, int lastConst,
		bool NoPinvInModel, FLOAT_TYPE prI, int nstates, int nRateCats,
		int nchar, int ncharGPU, dim3 dimBlock, dim3 dimGrid);
}
////////////////////////////////////////////////////////////////////////////////


#endif // #ifndef CUDAMANAGER
