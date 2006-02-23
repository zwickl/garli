#ifndef DEFS
#define DEFS

//unix is defined in the Makefile, so don't define
//another system type here
#ifndef UNIX
#undef WINDOWS
#undef MAC
#endif

#ifndef NDEBUG
#define NDEBUG
#endif

#undef INCLUDE_PERTURBATION
#undef SUBTREE_VERSION

#define DEF_MIN_BRLEN 1e-8
#define DEF_MAX_BRLEN 10.0
#define DEF_STARTING_BRLEN		0.05
#define MAXPATH   		256
#define DEF_PRECISION	10


#define MEM_DELETE_ARRAY(v)		{ delete [] v; v=NULL; }
#define MEM_NEW_ARRAY(a,t,n)	{ a = new t[n]; }

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
