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

#ifndef CUDAMANAGER
#define CUDAMANAGER

using namespace std;

class CudaManager {
	int nstates;
	int numRateCats;
	int nchar;

public:
	CudaManager() {
		nstates=0;
		numRateCats=0;
		nchar=0;
	}

	CudaManager(int nstatesIn, int numRateCatsIn, int ncharIn) {
		nstates = nstatesIn;
		numRateCats = numRateCatsIn;
		nchar = ncharIn;
	}

	~CudaManager() {
	}

	void BenchmarkGPU() {
		outman.UserMessage("nstates %i, numRateCats %i, nchar %i", nstates, numRateCats, nchar);
	}
};

#endif

