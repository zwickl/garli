// GARLI version 0.951 source code
// Copyright  2005-2006 by Derrick J. Zwickl
// All rights reserved.
//
// This code may be used and modified for non-commercial purposes
// but redistribution in any form requires written permission.
// Please contact:
//
//  Derrick Zwickl
//	National Evolutionary Synthesis Center
//	2024 W. Main Street, Suite A200
//	Durham, NC 27705
//  email: zwickl@nescent.org
//

#ifndef DEFS
#define DEFS

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

#undef OPT_DEBUG

//The ONLY thing that should need to be done to turn on memcheck leak detection
//should be defining MONITORING_ALLOCATION here
#undef MONITORING_ALLOCATION
#include "memchk.h"

#define CONSTRAINTS
#define STOCHASTIC_STARTING_BLENS
#undef IGNORE_SMALL_TOPO_IMP
#undef INCLUDE_PERTURBATION
#undef SUBTREE_VERSION
#undef ENABLE_CUSTOM_PROFILER
#undef SINGLE_PRECISION_FLOATS
#undef SWAP_BASED_TERMINATION

#undef OUTPUT_UNIQUE_TREES
#undef VARIABLE_OPTIMIZATION

#undef INPUT_RECOMBINATION
#define NUM_INPUT 12

#undef ALLOW_SINGLE_SITE

#undef EQUIV_CALCS

#ifdef SINGLE_PRECISION_FLOATS
	typedef float FLOAT_TYPE;
	#define ONE_POINT_ZERO 1.0f
	#define ZERO_POINT_FIVE 0.5f
	#define ZERO_POINT_ZERO 0.0f
	#define DEF_MIN_BRLEN 1e-8f
	#define DEF_MAX_BRLEN 100.0f
	#define DEF_STARTING_BRLEN 0.05f
#else
	typedef double FLOAT_TYPE;
	#define ONE_POINT_ZERO 1.0
	#define ZERO_POINT_FIVE 0.5
	#define ZERO_POINT_ZERO 0.0
	#define DEF_MIN_BRLEN 1e-8
	#define DEF_MAX_BRLEN 100.0
	#define DEF_STARTING_BRLEN 0.05
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
