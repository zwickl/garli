/* define this if pecr related code in files other than pecr.cpp has to be
 * compiled */
#undef GANESH

/* to debug the PECRMutate in pecr.cpp */
//#define DEBUG_PECR_MAIN

/* to debug the GenerateMatching routine in pecr.cpp */
//#define DEBUG_PECR_GEN_MATCHING

/* to debug MatchingToTree routine in pecr.cpp */
//#define DEBUG_PECR_MATCHING2TREE

/* to debug AlternateChoosePSubtreeGivenRoot routine in pecr.cpp */
//#define DEBUG_PECR_ALT_CHOOSE

/* to debug the SortMatching routine in pecr.cpp */
//#define DEBUG_PECR_SORT_MATCHING

/* to debug the allocation and initial values of the iedges array in the
 * IdentifyPSubtree routine */
//#define DEBUG_PECR_IEDGES

/* to debug the postorder traversal of the tree nodes in pecr.cpp */
//#define DEBUG_PECR_POSTORDER

/* This will give a sense of how often the fast p-subtree choosing routine
 * AlternateChoosePSubtreeGivenRoot succeeds */
//#define PECR_ALT_CHOOSE_SUCCESS_RATE

/* This will print out each mutation chosen, whether it is model or
 * topological */
//#define GANESH_CHECK_MOVE_PROPORTIONS

//#define DEBUG_PECR_TESTCHOOSE
//#define DEBUG_PECR_PARSIMONY_BRLEN
//efine PECR_SET_PARSIMONY_BRLEN
//#define PECR_SET_OLD_BRLEN


#define DEFAULT_P_VALUE     5 


/* TODO add PECR to all these error codes */
/* error and no-error codes */
#define ERROR_INSUFFICIENT_EDGES    -99
#define NOERROR_P_SUBTREE_LEAF_REACHED    1
#define NOERROR_BY_ALTCHOOSE 2
#define NOERROR_BY_CHOOSE 4
