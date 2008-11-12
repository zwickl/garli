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

	// Computes reference result
	void ComputeGarliCLA(FLOAT_TYPE* dest, const FLOAT_TYPE* Lpr,
					const FLOAT_TYPE* Rpr, const FLOAT_TYPE* LCL,
					const FLOAT_TYPE* RCL);

	void BenchmarkCLA();
};

#endif

