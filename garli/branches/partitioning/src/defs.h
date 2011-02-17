// GARLI version 2.0 source code
// Copyright 2005-2011 Derrick J. Zwickl
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
//

#ifndef DEFS
#define DEFS

#if defined(HAVE_CONFIG_H)
	#include "config.h"
#endif
//these will be defined by either the Microsoft compiler
//or the intel compiler when openmp support is turned on
//by compiling with /openmp (ms) or -openmp (icc)
//Nothing else should need to be defined anywhere to get 
//openMP working
#if defined (__OPENMP) || defined (_OPENMP)
	#include "omp.h"
	#define OPEN_MP
	#define OMP_INTINTCLA
	#define OMP_INTTERMCLA
	#define OMP_TERMDERIV
	#define OMP_INTDERIV
	
	#define OMP_INTINTCLA_NSTATE
	#define OMP_INTTERMCLA_NSTATE
	#define OMP_TERMDERIV_NSTATE
	#define OMP_INTDERIV_NSTATE
	#define OMP_INTSCORE_NSTATE
	#define OMP_TERMSCORE_NSTATE
#endif

/*
#ifndef NDEBUG
#undef NDEBUG
#endif
*/

#define USE_COUNTS_IN_BOOT

//#undef OPT_DEBUG

#define ONE_BRANCH_INS_DEL

//The ONLY thing that should need to be done to turn on memcheck leak detection
//should be defining MONITORING_ALLOCATION here
#undef MONITORING_ALLOCATION
#include "memchk.h"

#define ADAPTIVE_BOUNDED_OPT
#define ALT_NR_BAIL
#define PUSH_TO_MIN_BLEN
#define SUM_AA_REL_RATES
#define NEW_BUMPING
#define STOCHASTIC_STARTING_BLENS
#undef IGNORE_SMALL_TOPO_IMP
#undef INCLUDE_PERTURBATION
#undef SUBTREE_VERSION
//#undef ENABLE_CUSTOM_PROFILER
//#undef SINGLE_PRECISION_FLOATS
//#undef SWAP_BASED_TERMINATION

//#undef OUTPUT_UNIQUE_TREES
#undef VARIABLE_OPTIMIZATION

#undef INPUT_RECOMBINATION
#define NUM_INPUT 12

//#undef ALLOW_SINGLE_SITE

#undef EQUIV_CALCS

typedef double MODEL_FLOAT;

#ifdef SINGLE_PRECISION_FLOATS
	typedef float FLOAT_TYPE;
	#define ONE_POINT_ZERO 1.0f
	#define ZERO_POINT_FIVE 0.5f
	#define ZERO_POINT_ZERO 0.0f
	#define DEF_MIN_BRLEN 1e-8f
	#define DEF_MAX_BRLEN 100.0f
	#define DEF_STARTING_BRLEN 0.05f
	#define GARLI_FP_EPS FLT_EPSILON
	#define LUMP_LIKES
	#if !defined(LUMP_FREQ)
		#define LUMP_FREQ 400
	#endif
#else
	typedef double FLOAT_TYPE;
	#define ONE_POINT_ZERO 1.0
	#define ZERO_POINT_FIVE 0.5
	#define ZERO_POINT_ZERO 0.0
	#define DEF_MIN_BRLEN 1e-8
	#define DEF_MAX_BRLEN 100.0
	#define DEF_STARTING_BRLEN 0.05
	#define GARLI_FP_EPS DBL_EPSILON
	#if !defined(LUMP_FREQ)
		#define LUMP_FREQ 400
	#endif
#endif

#define MAXPATH   		256
#define DEF_PRECISION	8

#define MEM_DELETE_ARRAY(v)		{ delete [] v; v=NULL; }
#define MEM_NEW_ARRAY(a,t,n)	{ a = new t[n]; }

#ifdef BOINC
	#define WRITE_TO_FILE(ptr, size, count) write((void *) ptr, (size_t) size, (size_t) count)
	#define OUTPUT_CLASS MFILE

	#include "boinc_api.h"
	#include "filesys.h"
	#ifdef _WIN32
		#include "boinc_win.h"
	#else
		#include "config.h"
	#endif
#else
	#define WRITE_TO_FILE(ptr, size, count) write((const char *) ptr, (streamsize) size*count)
	#define OUTPUT_CLASS ofstream
#endif

//mpi message tags
#ifdef MPI_VERSION
#define TAG_PARAMS_SIZE			1
#define TAG_PARAMS				2
#define TAG_DATA_SIZE			3
#define TAG_DATA				4
#define TAG_TREE_STRINGS_COUNT	5
#define TAG_TREE_STRINGS_SIZE	6
#define TAG_TREE_STRINGS		7
#define TAG_CONFIG				8
#define TAG_QUIT				9
#define TAG_KAPPAS				10
#define TAG_NINDIVS				11
#define TAG_ACCEPT_COUNT		12
#define TAG_TREE_STRINGS_REQUEST 13
#define TAG_SCORE				14
#define TAG_PIS					15
#define TAG_MODEL				16
#define TAG_REMOTE_TYPE_SWITCH	17
#define TAG_SUBTREE_DEFINE		18
#define TAG_SUBTREE_ITERATION	19
#define TAG_PERTURB				20
#endif

#endif
