//This source all written by Ganesh Ganapathy for his PECR operator
//largely deprecated at this point

#include "tree.h"
#include "funcs.h"

#ifdef GANESH

/* TODO this declaration must be places in funcs.h and funcs.h must be
 * included in Tree.h */
double CalculateHammingDistance(const char *str1, const char *str2, 
                                const int *col_count, int nchar);


/* what is the max index of allNodes?
   
   proceedure to reroot - is it correct?

   random number generator - how does it work?
*/

/* TODO since p_value is a static member, eliminate all passing of p */

void Tree::PECRMutate(rng& rnd, double optPrecision) {


    int p = p_value;

#ifdef DEBUG_PECR_MAIN
    cout << "entered PECR. p = " << p << "\n";
#endif

    

    /* the p_subtree is rooted, so that #tips = p+2, and #vertices = 2p+3
     * */
    n_p_subtree_tips  = p+2;
    int num_p_subtree_vertices = 2*n_p_subtree_tips - 1;

/* TODO send error  or noerror codes as return values from PECRMutate */
    if (numTipsTotal < 4) {
        /* There is only one unrooted tree on less than 4 leaves */
        cout << "Tree is trivial\n";
        return;
    }

    if (p > numTipsTotal-3) {
        cout << "Insuffiecient edges in input tree \n";
        return;
    }
    /* there is one extra space in pedges and p_subtree_leaves, which calls
     * for an explanation, the last spaces pedges[p] and
     * p_subtree_leaves[p+2] are used
     * thus:
       if (p_subtree_root == root) {
            pedges[p] = the edge (root->left->next, root)
            p_subtree_leaves[p+2] = root->left->next
       }
       else {
            pedges[p] = the edge (p_subtree_root, p_subtree_root->anc)
            p_subtree_leaves[p+2] = p_subtree_root->anc
       }
     * */
    pedges = new int[p+1];
    p_subtree_leaves = new TreeNode*[p+3];    

    /* in GAML the array allNodes contains pointers to each node in the
     * tree. allNodes[0] = root; allNodes[1]..allNodes[numTipsTotal]
     * contain pointers to all leaves (tips), and
     * allNodes[numTipsTotal+1]..allNodes[numNodesTotal-1] contain pointers
     * to all internal nodes except the root. numNodesTotal =
     * 2*numTipsTotal-2 */

    /* root tree at the highest-labeled leaf; 
     * This means that the
     * parent tree is rooted at the ancestor of leaf
     * leaf being its middle
     * descendant */
    RootWithMiddle(numTipsTotal);

    int ret_val = IdentifyPSubtree(p_subtree_leaves, pedges, p, rnd);
    assert (ret_val != ERROR_INSUFFICIENT_EDGES);
    if (ret_val == ERROR_INSUFFICIENT_EDGES) {
        cout << "insufficient edges. exitting from PECRMutate\n";
        PECRCleanUp(p, ret_val);
        return;
    }

    int z = 1; 

    /* z = (2*n_p_subtree_tips-3)!! */
    for (int i = 1; i <= 2*n_p_subtree_tips-3; i = i+2) {
        z = z * i;
    }

/* random_int(z) gives a random number in the range [0 z-1] */
    int matching_number = rnd.random_int(z);

    A = new int[2*n_p_subtree_tips-2];

    for (int i = 0; i <= 2*n_p_subtree_tips-3; i++) {
        A[i] = i;
    }

#ifdef DEBUG_PECR_GEN_MATCHING
    for (int i = 0; i <= 2*n_p_subtree_tips-3; i++) {
        cout << A[i] << " ";
    }
    cout << "\n";
#endif
    GenerateMatching(A, matching_number, z, num_p_subtree_vertices);
#ifdef DEBUG_PECR_GEN_MATCHING
    for (int i = 0; i <= 2*n_p_subtree_tips-3; i++) {
        cout << A[i] << " ";
    }
    cout << "\n";
#endif

	//DJZ- Think that these are the correct nodes to mark dirty
	for(int i=0;i<p+3;i++){
		SweepDirtynessOverTree(p_subtree_leaves[i]);
		}

    MatchingToTree(A, 2*n_p_subtree_tips-3, p_subtree_leaves, pedges);

#if 1
    //DJZ trying some optimization (although probably not optimally)
	for(int i=0;i<p+3;i++){
		if(p_subtree_leaves[i]->nodeNum > numTipsTotal) OptimizeBranchesWithinRadius(p_subtree_leaves[i], optPrecision, 0);
		}
#endif

#ifdef PECR_SET_PARSIMONY_BRLEN
#if 0
    int **W;
    W = new int*[4];
    for (int i = 0; i < 4; i++) {
        W[i] = new int[4];
    }
#else
    int W[4][4];
#endif

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (i==j) {
                W[i][j] = 0;
            }
            else {
                W[i][j] = 1;
            }
        }
    }

    ComputeParsimonyBrLen(pedges, p_subtree_leaves, W);
#if 0
    for (int i = 0; i < 4; i++) {
        free(W[i]);
    }
    free(W);
#endif
#endif
    PECRCleanUp(p, ret_val);
    return;
}

int Tree::IdentifyPSubtree(TreeNode **p_subtree_leaves, int *pedges, int p, rng& rnd) {

#if 0
    C = new int[p+1];
    ComputeRealCatalan(p);
#endif
    
/* iedges[i] is the number of internal edges under allNodes[i]; it is a
 * private variable */
    
    iedges = new int[numNodesTotal];
#ifdef DEBUG_PECR_IEDGES
    cout << "allocated iedges \n";  
#endif

    int middle_node;
    middle_node = root->left->next->nodeNum;
    iedges[middle_node] = -1;

    /* middle_node will not be (should not be) visited in the following dfs to compute
     * iedges */
    ComputeNumInternalEdges(root, iedges);
#ifdef DEBUG_PECR_IEDGES
    for (int i = 0; i < numNodesTotal; i++) {
        cout << "iedges[" << i << "] = " << iedges[i] << "\n";
    }
#endif

    /* if middle_node was visited then iedges[middle_node] >= 0 */
    assert(iedges[middle_node] == -1);
    iedges[middle_node] = 0;


    int ret_val = ERROR_INSUFFICIENT_EDGES;
    int j = 0;
    int n_trials = (numTipsTotal*p)/10;

#ifdef PECR_ALT_CHOOSE_SUCCESS_RATE
    cout << "Trying AltChoose. n p n_trials: " << numTipsTotal
         << " " << p << " " << n_trials << "\n";
#endif

    while ((ret_val == ERROR_INSUFFICIENT_EDGES) && (j < n_trials)) {
#ifdef DEBUG_PECR_ALT_CHOOSE
        cout << "insufficient edges. retrying... ret_val = " << ret_val << "\n";
#endif
        /* choose an internal node u.a.r. Note that tips occupy
        * positions 1 through numTipsTotal. So we need to choose one intereger
        * u.a.r from [numTipsTotal+1 numNodesTotal], if numNodesTotal is
        * chosen, we will choose the root, which is numbered 0. */
        int k = rnd.random_int(numNodesTotal - (numTipsTotal));
        k = k + numTipsTotal + 1;


        if (k == numNodesTotal) {
            k = 0;
        }
        p_subtree_root = allNodes[k];
    
        if (p_subtree_root == root) {
            p_subtree_leaves[p+2] = p_subtree_root->left->next;
            pedges[p] = p_subtree_root->left->next->nodeNum;
        }
        else {
            p_subtree_leaves[p+2] = p_subtree_root->anc;
#if 0
        pedges[p] = p_subtree_root->anc->nodeNum;
#else
        pedges[p] = p_subtree_root->nodeNum;
#endif
        }
    #ifdef DEBUG_PECR_ALT_CHOOSE 
        cout << "nodeNum of psubtree root = " << k << "\n";
        cout << "p_subtree_leaves[p+2] = " << p_subtree_leaves[p+2]->nodeNum <<
            "\n";
    #endif
        running_index = 0;
        pedges_index = 0;
        ret_val = AlternateChoosePSubtreeGivenRoot(p_subtree_leaves, pedges, p, 
                               p_subtree_root, rnd);
        j++;
    }
    if (ret_val == NOERROR_P_SUBTREE_LEAF_REACHED) {
#ifdef  DEBUG_PECR_ALT_CHOOSE
        for (int i = 0; i <= running_index; i++) {
            cout << p_subtree_leaves[i]->nodeNum << " ";
            cout << "\n";
        }
#endif
#ifdef PECR_ALT_CHOOSE_SUCCESS_RATE
        cout << "AltChoose Succeeded. Attempt " << j+1 << "\n"; 
#endif
        return NOERROR_BY_ALTCHOOSE;
    }

#ifdef PECR_ALT_CHOOSE_SUCCESS_RATE
    cout << "AltChoose Unsuccessful. Trying FullChoose... \n";
#endif

    /* treecatalan is a (p+1) x numNodesTotal array of integers */
    /* treecatalan[k, j] is the number of k-subtrees under the subtree rooted
       at node with nodeNum = j */
    treecatalan = new int*[p+1];
    for (int i=0; i < p+1; i++) {
        treecatalan[i] = new int[numNodesTotal];
        treecatalan[i][root->left->next->nodeNum] = -1;
    }

    ComputeTreeCatalan(p);

    /* initialize with #p_subtrees under the root */
    int total_psubtrees  = treecatalan[p][0];
    for (int i=numTipsTotal+1; i < numNodesTotal; i++) {
        total_psubtrees += treecatalan[p][i];
    }

    int n_internal_nodes = numNodesTotal - numTipsTotal;
    int *intervals = new int[n_internal_nodes];

    intervals[0] = treecatalan[p][0];
    for (int i = numTipsTotal+1; i < numNodesTotal; i++) {
       intervals[i-numTipsTotal] = intervals[i-numTipsTotal-1] + 
                                   treecatalan[p][i];  
    }

    double d = total_psubtrees * rnd.uniform();

    int k = -1;
    for (int j = 0; j < n_internal_nodes; j++) {
        if  (d < intervals[j]) {
            k = j;
            break;
        }
    }
    assert (k >= 0);

    int rand_internal_node;
    if (k == 0) {
        rand_internal_node = k;
    }
    else {
        rand_internal_node = k + numTipsTotal;
    }

    assert(rand_internal_node >= 0);
    assert(rand_internal_node < numNodesTotal);

    p_subtree_root = allNodes[rand_internal_node];

    if (p_subtree_root == root) {
        p_subtree_leaves[p+2] = p_subtree_root->left->next;
        pedges[p] = p_subtree_root->left->next->nodeNum;
    }
    else {
        p_subtree_leaves[p+2] = p_subtree_root->anc;
#if 0
        pedges[p] = p_subtree_root->anc->nodeNum;
#else
        pedges[p] = p_subtree_root->nodeNum;
#endif
    }
    running_index = 0;
    pedges_index = 0;
    ret_val = ChoosePSubtreeGivenRoot(p_subtree_leaves, pedges, p, p_subtree_root,
                                      rnd);
    assert(ret_val == NOERROR_P_SUBTREE_LEAF_REACHED);
    return NOERROR_BY_CHOOSE;

}

/* make node->anc the root, and node its middle child */
void Tree::RootWithMiddle(int node_number) {

    TreeNode *node = allNodes[node_number];

    /* This does not guarantee that node is the middle descendant */
    /* really, this check should not be necessary, but adding it for
     * debugging purposes TODO */
    if (node->anc != root) {
        /* These two lines have been suggested by Derrick */
        //assert(root->dlen == -1);
        root->dlen = 0.1;
        RerootHere(node->anc->nodeNum);
        
    }
    
    /* Now make node the middle descendant of the new root */
    
    /* check if the left child of root is the node */
    if (root->left == node) {
        root->left = node->next;
        root->left->prev = NULL;
        root->left->next = node;
        node->prev = root->left;
        node->next = root->right;
        root->right->prev = node;
        root->AdjustClasForReroot(UPLEFT);  //DJZ this will update the cla directions at the root
    }
    else if (root->right == node) {
        root->right = node->prev;
        root->right->next = NULL;
        root->right->prev = node;
        node->next  = root->right;
        node->prev  = root->left;
        root->left->next = node;
        root->AdjustClasForReroot(UPRIGHT); //DJZ this will update the cla directions at the root
    }
    else {
        /* in this case nothing needs to be done, since node is already the
           middle descendant of root */
//        MakeAllNodesDirty();
        return;
    }

/* conservatively make all nodes dirty so that all conditional likelihoods
 * are recalculated */
//    MakeAllNodesDirty();
    return;
}

void Tree::GenerateMatching(int *A, int matching_number, 
                            int num_matchings, int num_p_subtree_vertices) {
    
    int quotient[100000];
    
    int double_factorial_term = num_p_subtree_vertices - 2;

    for(int j = 0; j <= double_factorial_term; j++) {
        quotient[j] = 0;
    }

    int dividend = matching_number;
    int divisor  = num_matchings/double_factorial_term;
    int index = 0;
    int remainder = 0;

    while (dividend > 0) {
        quotient[index] = dividend/divisor;
        remainder = dividend - quotient[index]*divisor;
        dividend = remainder;
        double_factorial_term = double_factorial_term - 2;
        divisor = divisor/double_factorial_term;
        index++;
    }

    for (int i = 0; i < index; i++) {

        int temp0, temp1, fixed_swap_index;
        temp0 = quotient[i];
        fixed_swap_index = 2*i+1;
        temp1 = A[fixed_swap_index]; 
        A[fixed_swap_index] = A[fixed_swap_index + temp0];
        A[fixed_swap_index + temp0] = temp1;
    }
    return;
}

#if 0
void Tree::ComputeRealCatalan(int p) {
    
    C[0] = 1;
    C[1] = 2;

    for (int i = 2; i <= p; i++) {
        C[i] = C[i-1];
        for (int j = 1; j <= i-1; j++) {
            C[i] = C[i] + C[j-1]*C[i-1-j];
        }
        C[i] = C[i] + C[i-1];
    }
    return;
}
#endif

int  Tree::AlternateChoosePSubtreeGivenRoot(TreeNode **p_subtree_leaves,
                                            int *pedges,
                                            int p, TreeNode *subtree_root, 
                                            rng& rnd) {
    if (p > iedges[subtree_root->nodeNum]) {
        return ERROR_INSUFFICIENT_EDGES;
    }

    if (IsALeaf(subtree_root)) {
        /* subtree_root is a leaf */
#ifdef DEBUG_PECR_ALT_CHOOSE
        cout << "hit a leaf\n";
#endif
        return ERROR_INSUFFICIENT_EDGES;
    }

#ifdef DEBUG_PECR_ALT_CHOOSE
    cout << "p = " << p << "\n";
#endif

    if (p == 0) {
       p_subtree_leaves[running_index] = subtree_root->left; 
#ifdef DEBUG_PECR_ALT_CHOOSE
       cout << "running_index center = " << running_index << ", " <<
                                     p_subtree_leaves[running_index]->nodeNum << "\n";
#endif
       running_index++;
       p_subtree_leaves[running_index] = subtree_root->right; 
#ifdef DEBUG_PECR_ALT_CHOOSE
       cout << "running_index center = " << running_index << ", " <<
                                     p_subtree_leaves[running_index]->nodeNum << "\n";
#endif
       running_index++;
       return NOERROR_P_SUBTREE_LEAF_REACHED;
    }    

    double d = C[p] * rnd.uniform();
    int k = -1;

#if 0
    int *local_intervals;
    local_intervals = new int[p+1];

    local_intervals[0] = C[p-1];

    for (int j = 1; j <= p-1; j++) {
        local_intervals[j] = local_intervals[j-1] + (C[j-1] * C[p-j-1]);
    }

    local_intervals[p] = local_intervals[p-1] + C[p-1];


    for (int j = 0; j <= p; j++) {
        if  (d < local_intervals[j]) {
            k = j;
            break;
        }
    }
#else
    int floord = (int)d;
    assert ((floord >= 0) && (floord <= C[p]-1));
    k = inv_realcat_intervals[p][floord];
#endif
    assert ((k >= 0) && (k <= p));

#ifdef DEBUG_PECR_ALT_CHOOSE
    cout << "k = " << k << "\n";
#endif

    if (k == 0) {
        pedges[pedges_index] = subtree_root->right->nodeNum;
        assert (pedges[pedges_index] < numNodesTotal);
        pedges_index++;
        p_subtree_leaves[running_index] = subtree_root->left; 
#ifdef DEBUG_PECR_ALT_CHOOSE
        cout << "running_index left = " << running_index << ", " <<
                                      p_subtree_leaves[running_index]->nodeNum << "\n";
#endif
        running_index++;
        return AlternateChoosePSubtreeGivenRoot(p_subtree_leaves, pedges, p-1,
                                                subtree_root->right, rnd); 
    } 
    if (k == p) {
        pedges[pedges_index] = subtree_root->left->nodeNum;
        assert (pedges[pedges_index] < numNodesTotal);
        pedges_index++;
        p_subtree_leaves[running_index] = subtree_root->right; 
#ifdef DEBUG_PECR_ALT_CHOOSE
        cout << "running_index right = " << running_index << ", " <<
                                      p_subtree_leaves[running_index]->nodeNum << "\n";
#endif
        running_index++;
        return AlternateChoosePSubtreeGivenRoot(p_subtree_leaves, pedges, p-1,
                                                subtree_root->left, rnd); 
    }
    
    pedges[pedges_index] = subtree_root->right->nodeNum;
    pedges_index++;
    pedges[pedges_index] = subtree_root->left->nodeNum;
    pedges_index++;

    int retval = AlternateChoosePSubtreeGivenRoot(p_subtree_leaves, pedges, k-1,
                                                  subtree_root->left, rnd);
    if (retval == ERROR_INSUFFICIENT_EDGES) {
        return retval;
    }

    return AlternateChoosePSubtreeGivenRoot(p_subtree_leaves, pedges, p-k-1,
                                            subtree_root->right, rnd); 
}

void Tree::ComputeNumInternalEdges(TreeNode *node, int *iedges) {


    if (IsALeaf(node)) {
        /* if node is a leaf, the #internal edges under it is 0 */
        iedges[node->nodeNum] = 0;
        return; 
    }   

    ComputeNumInternalEdges(node->left, iedges); 
    ComputeNumInternalEdges(node->right, iedges); 

    if (IsALeaf(node->left) && IsALeaf(node->right)) {
        /* the left child and the right child are both leaves. so the
         * #internal edges under node is 0; */
        iedges[node->nodeNum] = 0;
        return;
    }

    if (IsALeaf(node->left)) {
        iedges[node->nodeNum] = 1 + iedges[node->right->nodeNum];
        return;
    }
    if (IsALeaf(node->right)) {
        iedges[node->nodeNum] = 1 + iedges[node->left->nodeNum];
        return;
    }

    iedges[node->nodeNum] = 2 + iedges[node->left->nodeNum] + 
                                iedges[node->right->nodeNum];
    return;
}    

bool Tree::IsALeaf(TreeNode *node) {
    return ((node->nodeNum > 0) && (node->nodeNum <= numTipsTotal)); 
}

bool Tree::IsALeaf(int nodeNum) {
    return ((nodeNum > 0) && (nodeNum <= numTipsTotal)); 
}

#ifdef PECR_SET_PARSIMONY_BRLEN
bool Tree::IsLeafMaskOn(TreeNode *node) {
    return (node->leaf_mask);
}

bool Tree::IsLeafMaskOn(int nodeNum) {
    return (allNodes[nodeNum]->leaf_mask);
}
#endif

void Tree::SortMatching(int *A, int max_index) {
    
    /* max_index = 2l-3, where l is #p_subtree tips and so should be odd */ 
    assert((max_index % 2) == 1);
    QSort(A, 0, (max_index-1)); 
    return;
}

void Tree::QSort(int *A, int left, int right) {

    /* left and right should be even. assert. TODO */
    assert ((left % 2) == 0);
    assert ((right % 2) == 0);

    if (right > left) {
        /* currently setting pivot_index to left; TODO this should be
         * randomized */
        int pivot_index = left;
        int pivot_value = Key(A, pivot_index);
        int pivot_new_index = Partition(A, left, right, pivot_index);
        assert ((pivot_new_index % 2) == 0);
        QSort(A, left, pivot_new_index-2);
        QSort(A, pivot_new_index+2, right);
    }
    return;
}

int Tree::Key(int *A, int i) {

    assert ((i % 2) == 0);
    /* i should be even  TODO */
    if (A[i] < A[i+1]) {
        return A[i];
    }
    return A[i+1];
}

int Tree::MaxLabelPair(int *A, int i) {

    assert ((i % 2) == 0);
    /* i should be even  TODO */
    if (A[i] > A[i+1]) {
        return A[i];
    }
    return A[i+1];
}

int Tree::Swap(int *A, int index1, int index2) {

    /* index1 and index2 should be even. assert. TODO  */
    
    assert((index1 % 2) == 0);
    assert((index2 % 2) == 0);

    int t1 = A[index1];
    int t2 = A[index1+1];

    A[index1] = A[index2];
    A[index1+1] = A[index2+1];

    A[index2] = t1;
    A[index2+1] = t2;
}

/* TODO: in comparing label pairs, initially I had used the *minimum* of
 * the two labels in the pair as the Key. The present Key() function
 * relects that. It turns out that what actually works as a valid key so
 * that the MatchingToTree works fine is the *maximum* of two the labels in
 * a pair. Either change the Key() function to relect this fact and
 * eliminate MaxLabelPair(), or retain MaxLabelPair() and eliminate Key().
 * The latter option will be more readable, I think. */ 
int Tree::CompareLabelPairs(int *A, int index1, int index2) {

    /* if (A[index1], A[index1+1]) > (A[index2], A[index2+1]) 
            return -1;
       if (A[index1], A[index1+1]) < (A[index2], A[index2+1]) 
            return 1;
       the two pairs compared cannot be equal, unless i = j */

    /* index1 and index2 should be even. assert. TODO  */
    
    assert((index1 % 2) == 0);
    assert((index2 % 2) == 0);
    
    if (index1 == index2) {
        return 1;
    }

#if 1
    int max1 = MaxLabelPair(A, index1);
    int max2 = MaxLabelPair(A, index2);
    if (max1 < max2) {
        return 1;
    }
    assert(max1 > max2);
    return -1;
#else
    if ((MaxLabelPair(A, index1) < n_p_subtree_tips) && 
        (MaxLabelPair(A, index2) < n_p_subtree_tips)) {
        /* both index1 and index2 denote leaf siblings */
        if (Key(A, index1) > Key(A, index2)) {
            return -1;
        }
#ifdef DEBUG_PECR_SORT_MATCHING
        cout << "Keys = " << Key(A, index1) << " " << Key(A, index2)
             <<  "\n";
        cout << "indices = " << index1 << " " << index2 << "\n";
#endif
        assert (Key(A, index1) < Key(A, index2));
        return 1;
    } 
    if ((MaxLabelPair(A, index1) < n_p_subtree_tips) && 
        (MaxLabelPair(A, index2) >= n_p_subtree_tips)) {
#ifdef DEBUG_PECR_SORT_MATCHING
        cout << "Keys = " << Key(A, index1) << " " << Key(A, index2)
             <<  "\n";
        cout << "indices = " << index1 << " " << index2 << "\n";
#endif
        return 1;
    } 
    if ((MaxLabelPair(A, index1) >= n_p_subtree_tips) && 
        (MaxLabelPair(A, index2) < n_p_subtree_tips)) {
#ifdef DEBUG_PECR_SORT_MATCHING
        cout << "Keys = " << Key(A, index1) << " " << Key(A, index2)
             <<  "\n";
        cout << "indices = " << index1 << " " << index2 << "\n";
#endif
        return -1;
    } 
    if ((MaxLabelPair(A, index1) >= n_p_subtree_tips) && 
        (MaxLabelPair(A, index2) >= n_p_subtree_tips)) {
        /* neither index1 nor index2 denote leaf siblings */
        if (Key(A, index1) > Key(A, index2)) {
            return -1;
        }
#ifdef DEBUG_PECR_SORT_MATCHING
        cout << "Keys = " << Key(A, index1) << " " << Key(A, index2)
             <<  "\n";
        cout << "indices = " << index1 << " " << index2 << "\n";
#endif
        assert (Key(A, index1) < Key(A, index2));
        return 1;
    } 
#endif
}

int Tree::Partition(int *A, int left, int right, int pivot_index) {
     
    int pivot_value = Key(A, pivot_index);

    Swap(A, pivot_index, right);

    int store_index = left;

    for (int i = left; i <= right-2; i = i+2) {
#if 0
        int k = Key(A, i);
        if (k <= pivot_value) {
            Swap(A, store_index, i); 
            store_index = store_index + 2;
        }
#else
        pivot_index = right;
        if (CompareLabelPairs(A, i, pivot_index) == 1) {
            Swap(A, store_index, i); 
            store_index = store_index + 2;
        }
#endif
    }

    Swap(A, right, store_index);
    assert ((store_index % 2) == 0);
    return store_index;
}


void Tree::MatchingToTree(int *A, int max_index, TreeNode **p_subtree_leaves, int *pedges)
{

#ifdef DEBUG_PECR_MATCHING2TREE
#if 0
    A[0] = 0;
    A[1] = 9;
    A[2] = 6;
    A[3] = 5;
    A[4] = 2;
    A[5] = 3;
    A[6] = 1;
    A[7] = 4;
    A[8] = 7;
    A[9] = 10;
    A[10] = 11;
    A[11] = 8;
    for (int i = 0; i <= max_index; i++) {
        cout << A[i] << " ";
    }
    cout << "\n";
#endif
#endif

    SortMatching(A, max_index);

/* variables l and n_p_subtree_tips denote the same quantity - the number
 * of tips in the identified p subtree. n_p_subtree_tips = l = p+2; I use l
 * in this routine for readability. */

    int l = (max_index + 3)/2;
    assert (l == n_p_subtree_tips);

    /* Translate numbers in A to nodeNum */

    /* translate[max_index+1] = p_subtree_root->nodeNum */
    /* remember: the root of the p_subtree does not change during a pECR
     * move, and the root itself is not mentioned in the
     * array A that holds the matching corresponding to the new p_subtree.
     * But when we build the new tree, in the last step as we set up the two
     * remaining nodes as the left and right children of the p_subtree
     * root, we will need p_subtree_root's nodeNum. */
    int *translate = new int[max_index+2];

#ifdef DEBUG_PECR_MATCHING2TREE
    for (int i = 0; i <= max_index; i++) {
        cout << A[i] << " ";
    }
    cout << "\n";
    for (int i = 0; i < l; i++) {
        cout << p_subtree_leaves[i]->nodeNum << " ";
        cout << "\n";
    }
    cout << "\n\n";
#endif

    for (int i = 0; i <= max_index; i++) {
        if (A[i] < l) {
            /* then A[i] denotes a leaf */
            translate[A[i]] = p_subtree_leaves[A[i]]->nodeNum;
        }
        else {
            translate[A[i]] = pedges[A[i]-l];
        }
#ifdef DEBUG_PECR_MATCHING2TREE
        cout << i << " " << A[i] << " " << translate[A[i]] << "\n";
#endif
    }
    translate[max_index+1] = p_subtree_root->nodeNum;
#ifdef DEBUG_PECR_MATCHING2TREE
        cout << max_index+1 << " " << "    " << translate[max_index+1] << "\n";
#endif

/* current_node_label is the label according to A; this has to be translated to a nodeNum
 * using translate */
    int current_node_label = l;
#ifdef DEBUG_PECR_MATCHING2TREE
    cout << "l = " << l << "\n";
    cout << "p_subtree_leaves[l] = " << p_subtree_leaves[l]->nodeNum <<
            "\n";
#endif

    for (int i = 0; i < max_index; i = i + 2) {
        TreeNode *parent = allNodes[translate[current_node_label]];        
#if 0
        parent->left =NULL;
        parent->right =NULL;
        parent->next =NULL;
        parent->prev =NULL;
        parent->anc = NULL;
#endif     
#ifdef DEBUG_PECR_MATCHING2TREE
        cout << "\n";
        cout << "parent = " << current_node_label << " " 
             << translate[current_node_label] << "\n";
        cout << "left = " << A[i] << " " 
             << translate[A[i]] << "\n";
        cout << "right = " << A[i+1] << " " 
             << translate[A[i+1]] << "\n";
        cout << "\n";
#endif
#if 1
        assert(A[i] < current_node_label);
        assert(A[i+1] < current_node_label);
#endif
        parent->left = allNodes[translate[A[i]]];
        parent->right = allNodes[translate[A[i+1]]];
#ifndef PECR_SET_OLD_BRLEN
        parent->left->dlen = parent->right->dlen = 0.05;
#endif
        if ((i == max_index - 1) && (parent == root)) {
            /* max_index is odd, and so this is the last iteration; plus the
             * root of our p_subtree (parent) is the root of the entire tree */ 
            parent->left->anc = parent;
            parent->right->anc = parent;
            assert(p_subtree_leaves[l]->anc == parent);
            parent->left->next = p_subtree_leaves[l];
            p_subtree_leaves[l]->next = parent->right;
            parent->right->next = NULL;
            parent->right->prev = p_subtree_leaves[l];
            p_subtree_leaves[l]->prev = parent->left;
            parent->left->prev = NULL;
        }
        else {        
            parent->left->anc = parent;
            parent->right->anc = parent;
            parent->left->next = parent->right;
            parent->left->prev = NULL;
            parent->right->prev = parent->left;
            parent->right->next = NULL;
        }
        current_node_label++;
    }

    free(translate);
    return;
}        

/* TODO follow naming conventions, if any, for identifiers */

void Tree::ComputeTreeCatalan(int p) {

    postorder = new int[numNodesTotal];
    postorder_reverse = new int[numNodesTotal]; 

    postorder_index = 0;
    PostOrderTraverse(root);

    /* the middle descendant of the root is arbitrarily assigned the last
     * position in the PO traversal */
    postorder[postorder_index] = root->left->next->nodeNum;

#ifdef DEBUG_PECR_POSTORDER
    cout << "PO index = " << postorder_index << "\n";
    assert(postorder_index == numNodesTotal-1);
    for (int i = 0; i < numNodesTotal; i++) {
        postorder_reverse[postorder[i]] = i;
    }
    for (int i = 0; i < numNodesTotal; i++) {
        int middle = root->left->next->nodeNum;
        if ((i != middle) && (IsALeaf(i) == false)) {
            assert(postorder_reverse[i] >
                   postorder_reverse[allNodes[i]->left->nodeNum]);
            assert(postorder_reverse[i] >
                   postorder_reverse[allNodes[i]->right->nodeNum]);
        }
    }
#endif


    for(int i = 0; i < numNodesTotal-1; i++) {
        if (IsALeaf(postorder[i])) {
           /* allNodes[postorder[i]] is a leaf */
            for (int j = 0; j <= p; j++) {
                treecatalan[j][postorder[i]] = 0;
            }    
            continue;
        }
        if (IsALeaf(allNodes[postorder[i]]->left) &&
            IsALeaf(allNodes[postorder[i]]->right)) {
            treecatalan[0][postorder[i]] = 1;
            for (int j = 1; j <= p; j++) {
                treecatalan[j][postorder[i]] = 0;
            }
            continue;
        }
        if (IsALeaf(allNodes[postorder[i]]->left)) {
            treecatalan[0][postorder[i]] = 1;
            treecatalan[1][postorder[i]] = 1;
            for (int j=2; j<=p; j++) {
                treecatalan[j][postorder[i]] =
                treecatalan[j-1][allNodes[postorder[i]]->right->nodeNum];
            }
            continue;
        }
        if (IsALeaf(allNodes[postorder[i]]->right)) {
            treecatalan[0][postorder[i]] = 1;
            treecatalan[1][postorder[i]] = 1;
            for (int j=2; j<=p; j++) {
                treecatalan[j][postorder[i]] =
                treecatalan[j-1][allNodes[postorder[i]]->left->nodeNum];
            }
            continue;
        }
        assert(IsALeaf(allNodes[postorder[i]]->left->nodeNum) == false);
        assert(IsALeaf(allNodes[postorder[i]]->right->nodeNum) == false);

        treecatalan[0][postorder[i]] = 1;
        treecatalan[1][postorder[i]] = 2;
        for(int j=2; j <= p; j++) {
            treecatalan[j][postorder[i]] =
                treecatalan[j-1][allNodes[postorder[i]]->left->nodeNum] +
                treecatalan[j-1][allNodes[postorder[i]]->right->nodeNum];
            int temp_sum = treecatalan[j][postorder[i]];
            for(int k=2; k <= j; k++) {
                temp_sum = temp_sum +
                           treecatalan[j-k][allNodes[postorder[i]]->left->nodeNum]* 
                           treecatalan[k-2][allNodes[postorder[i]]->right->nodeNum];
            }
            treecatalan[j][postorder[i]] = temp_sum;
        }
    } 
    return;
}

void Tree::PostOrderTraverse(TreeNode *node) {
    
    if (IsALeaf(node)) {
        postorder[postorder_index] = node->nodeNum;
        postorder_index++;
        return;
    } 
    assert(node->left != NULL);
    assert(node->right != NULL);
    PostOrderTraverse(node->left);
    PostOrderTraverse(node->right);
    postorder[postorder_index] = node->nodeNum;
    postorder_index++;
    return;
}


#ifdef PECR_SET_PARSIMONY_BRLEN
void Tree::MaskedPostOrderTraverse(TreeNode *node) {
    
#ifdef DEBUG_PECR_PARSIMONY_BRLEN
        cout << "node = " << node->nodeNum << "\n";
#endif
    /* note that leaves are already in the postorder array, even before
     * MaskedPostOrderTraverse is called for the first time. */
    if (IsLeafMaskOn(node)) {
        //node->leaf_mask = false;
        return;
    } 
    assert(node->left != NULL);
    assert(node->right != NULL);
    MaskedPostOrderTraverse(node->left);
    MaskedPostOrderTraverse(node->right);
    masked_postorder[masked_postorder_index] = node->nodeNum;
    masked_postorder_reverse[node->nodeNum] = masked_postorder_index;
    masked_postorder_index++;
    return;
}
#endif

int  Tree::ChoosePSubtreeGivenRoot(TreeNode **p_subtree_leaves,
                                   int *pedges,
                                   int p, TreeNode *subtree_root,  
                                   rng& rnd) {
    assert(treecatalan[p][subtree_root->nodeNum] > 0);

    if (p == 0) {
       p_subtree_leaves[running_index] = subtree_root->left; 
#ifdef DEBUG_PECR_CHOOSE
       cout << "running_index center = " << running_index << ", " <<
                                     p_subtree_leaves[running_index]->nodeNum << "\n";
#endif
       running_index++;
       p_subtree_leaves[running_index] = subtree_root->right; 
#ifdef DEBUG_PECR_CHOOSE
       cout << "running_index center = " << running_index << ", " <<
                                     p_subtree_leaves[running_index]->nodeNum << "\n";
#endif
       running_index++;
       return NOERROR_P_SUBTREE_LEAF_REACHED;
    }    

    assert(IsALeaf(subtree_root) == false);

    TreeNode *left_node = subtree_root->left;
    TreeNode *right_node = subtree_root->right;

    int left_num = left_node->nodeNum;
    int right_num = right_node->nodeNum;
    

    int *intervals;
    intervals = new int[p+1];

    intervals[0] = treecatalan[p-1][right_num];
    for (int j = 1; j <= p-1; j++) {
        intervals[j] = intervals[j-1] + (treecatalan[j-1][left_num] *
                                         treecatalan[p-j-1][right_num]);
    }
    intervals[p] = intervals[p-1] + treecatalan[p-1][left_num];

    double d = treecatalan[p][subtree_root->nodeNum] * rnd.uniform();

    int k = -1;
    for (int j = 0; j <= p; j++) {
        if  (d < intervals[j]) {
            k = j;
            break;
        }
    }
    assert (k >= 0);
#ifdef DEBUG_PECR_CHOOSE
    cout << "k = " << k << "\n";
#endif
    free(intervals);

    if (k == 0) {
        pedges[pedges_index] = subtree_root->right->nodeNum;
        assert (pedges[pedges_index] < numNodesTotal);
        pedges_index++;
        p_subtree_leaves[running_index] = subtree_root->left; 
#ifdef DEBUG_PECR_CHOOSE
        cout << "running_index left = " << running_index << ", " <<
                                      p_subtree_leaves[running_index]->nodeNum << "\n";
#endif
        running_index++;
        return ChoosePSubtreeGivenRoot(p_subtree_leaves, pedges, p-1,
                                                subtree_root->right, rnd); 
    } 
    if (k == p) {
        pedges[pedges_index] = subtree_root->left->nodeNum;
        assert (pedges[pedges_index] < numNodesTotal);
        pedges_index++;
        p_subtree_leaves[running_index] = subtree_root->right; 
#ifdef DEBUG_PECR_CHOOSE
        cout << "running_index right = " << running_index << ", " <<
                                      p_subtree_leaves[running_index]->nodeNum << "\n";
#endif
        running_index++;
        return ChoosePSubtreeGivenRoot(p_subtree_leaves, pedges, p-1,
                                                subtree_root->left, rnd); 
    }
    
    pedges[pedges_index] = subtree_root->right->nodeNum;
    pedges_index++;
    pedges[pedges_index] = subtree_root->left->nodeNum;
    pedges_index++;

    int retval = ChoosePSubtreeGivenRoot(p_subtree_leaves, pedges, k-1,
                                                  subtree_root->left, rnd);

    assert(retval == NOERROR_P_SUBTREE_LEAF_REACHED);

    return ChoosePSubtreeGivenRoot(p_subtree_leaves, pedges, p-k-1,
                                            subtree_root->right, rnd); 
}
void Tree::PECRCleanUp(int p, int code) {
    root->CheckTreeFormation();
    MakeAllNodesDirty();
    if (code == ERROR_INSUFFICIENT_EDGES) {
        free(pedges);
        free(p_subtree_leaves);
        free(iedges);
        if (random_p == true) {
            free(C);
            for(int a=0; a <= p; a++) {
                free(realcat_intervals[a]);
                free(inv_realcat_intervals[a]);
            }
        }
        return;
    }
    if (code == NOERROR_BY_ALTCHOOSE) {
        free(pedges);
        free(p_subtree_leaves);
        free(iedges);
        if (random_p == true) {
            free(C);
            for(int a=0; a <= p; a++) {
                free(realcat_intervals[a]);
                free(inv_realcat_intervals[a]);
            }
        }
        free(A);
        return;
    }
    if (code == NOERROR_BY_CHOOSE) {
        free(pedges);
        free(p_subtree_leaves);
        free(iedges);
        if (random_p == true) {
            free(C);
            for(int a=0; a <= p; a++) {
                free(realcat_intervals[a]);
                free(inv_realcat_intervals[a]);
            }
        }
        free(A);
        for (int i = 0; i < p+1; i++) {
            free(treecatalan[i]);
        }
        free(treecatalan);
        free(postorder);
        free(postorder_reverse);
        return;
    }
}

#ifdef PECR_SET_PARSIMONY_BRLEN
void Tree::ComputeParsimonyBrLen(int *pedges, TreeNode **p_subtree_leaves,
                                 int W[][4]) {
    
    int ***pscore;
    int p = p_value;
    int l = p+2;
    int n_psubtree_nodes = 2*p+3; 
    int n_psubtree_leaves = l;


    int nstates = mod->NStates();
    int nsites = data->NChar();

/* the  number of times a column appears in the data */
    const int *col_count = data->GetCounts();
#if 1
    for (int c = 0; c < nsites; c++) {
        cout << col_count[c] << " ";
    }
    cout << "\n";
#endif

#ifdef DEBUG_PECR_PARSIMONY_BRLEN
    char *string = new char[nsites];

    for (int i = 0; i < n_psubtree_leaves; i++) {
        if ((p_subtree_leaves[i]->nodeNum > 0) &&
           (p_subtree_leaves[i]->nodeNum <= numTipsTotal)) {
               GetInternalStateString(string, p_subtree_leaves[i]->nodeNum);    
        }
    }
    exit(0);

    for (int i = 0; i < p; i++) {
        cout << "node:anc:left:right " << pedges[i] << ":"
                                       << allNodes[pedges[i]]->anc->nodeNum << ":"
                                       << allNodes[pedges[i]]->left->nodeNum << ":"
                                       << allNodes[pedges[i]]->right->nodeNum << "\n";
    }
    cout << "node " << pedges[p] << "\n";
#endif

    


    pscore = new int**[n_psubtree_nodes];
    for (int i = 0; i < n_psubtree_nodes; i++) {
        pscore[i] = new int*[nsites];
        for (int j = 0; j < nsites; j++) {
            pscore[i][j] = new int[nstates];
        }
    }


    masked_postorder_index = 0;
    masked_postorder = new int[n_psubtree_nodes];
    masked_postorder_reverse = new int[numNodesTotal];

    for (int i = 0; i < numNodesTotal; i++) {
        masked_postorder_reverse[i] = -1;
        allNodes[i]->leaf_mask = false;
    }

    /* we want a post-order the p-subtree nodes where all the leaves occur
     * before any internal node. The following code reflects that. */

    for (int i = 0; i < n_psubtree_leaves; i++) {
        masked_postorder[i] = p_subtree_leaves[i]->nodeNum;
        masked_postorder_reverse[p_subtree_leaves[i]->nodeNum] = i;
        masked_postorder_index++;
    }

    /* remember that p_subtree_leaves[p+2] holds the middle descendant of
     * the p_subtree_root */
    TreeNode *p_subtree_root;
    if (pedges[p] == p_subtree_leaves[p+2]->nodeNum) {
        p_subtree_root = p_subtree_leaves[p+2]->anc;
#ifdef DEBUG_PECR_PARSIMONY_BRLEN
        cout << "subtree root is main root\n";
#endif
    }
    else {
        assert(p_subtree_leaves[p+2] == allNodes[pedges[p]]->anc);
        p_subtree_root = allNodes[pedges[p]];
#ifdef DEBUG_PECR_PARSIMONY_BRLEN
        cout << "subtree root is not main root\n";
#endif
    }

#ifdef DEBUG_PECR_PARSIMONY_BRLEN
        cout << "root = " << p_subtree_root->nodeNum << "\n";
#endif

    /* the mask is not set for p_subtree_leaves[p+2], the middle descendant
     * of the p_subtree_root */
    for (int i = 0; i < n_psubtree_leaves; i++) {
        p_subtree_leaves[i]->leaf_mask = true;
    }

#ifdef DEBUG_PECR_PARSIMONY_BRLEN
    for (int i = 0; i < n_psubtree_leaves; i++) {
        cout << p_subtree_leaves[i]->nodeNum << ":" <<
                p_subtree_leaves[i]->anc->nodeNum << ":" << 
                p_subtree_leaves[i]->leaf_mask << "\n";
    }
    cout << "last leaf " << p_subtree_leaves[n_psubtree_leaves]->nodeNum <<"\n";
#endif
    MaskedPostOrderTraverse(p_subtree_root);

    /* remove the masks */
    for (int i = 0; i < n_psubtree_leaves; i++) {
        p_subtree_leaves[i]->leaf_mask = false;
    }

    /* inf should be greater than the upperbound for pscore[i][j][k]. Now
     * pscore[i][j][k] can't be greater than the #edges in the p-subtree as
     * the parsimony score can increase by at most 1 in each edge. So we
     * set inf = 10*p*max(W), where max(W) denotes the maximum entry in W.
     * */
    /* TODO max_weight = 1 only for now. Generalize this. */
    int max_weight = 1;
    int inf = 10*p*max_weight;


    char **internal_string;
    internal_string = new char*[n_psubtree_nodes];

    for (int i = 0; i < n_psubtree_nodes; i++) {
       internal_string[i] = new char[nsites]; 
    }

    /* for each leaf of the p-subtree*/
    for (int i = 0; i < n_psubtree_leaves; i++) {
        if ((p_subtree_leaves[i]->nodeNum > 0) &&
           (p_subtree_leaves[i]->nodeNum <= numTipsTotal)) {
            TranslateTipString(p_subtree_leaves[i]->tipData,
                               internal_string[i], nsites);
        }
        else {
            GetInternalStateString(internal_string[i], p_subtree_leaves[i]->nodeNum);
        }
        /* for each site */
        for (int j = 0; j < nsites; j++) {
            /* for each base */
            for (int k = 0; k < 4; k++) {
                if (internal_string[i][j] == k) {
                    pscore[i][j][k] = 0;
                }
                else {
                    pscore[i][j][k] = inf;
                }
            }
        }
    }

#ifdef DEBUG_PECR_PARSIMONY_BRLEN
    cout << "masked_postorder_index = " << masked_postorder_index << "\n";
#endif
    /* for each internal node in the subtree, in post-order */
    for (int i = l; i < n_psubtree_nodes; i++) {
        int left  = masked_postorder_reverse[allNodes[masked_postorder[i]]->left->nodeNum];
        int right = masked_postorder_reverse[allNodes[masked_postorder[i]]->right->nodeNum];
        /* for each site */
        for (int j = 0; j < nsites; j++) {
            /* for each state */
            for (int k = 0; k < nstates; k++) {
                int minleft = inf;
                int minright = inf;
                for (int r = 0; r < nstates; r++) {
                    if (W[k][r] + pscore[left][j][r] < minleft) {
                        minleft = W[k][r] + pscore[left][j][r];
                    }
                    if (W[k][r] + pscore[right][j][r] < minright) {
                        minright = W[k][r] + pscore[right][j][r];
                    }
                }
                pscore[i][j][k] = minleft + minright;
            }
        }
    }

    /* Now assign the states in pre-order. Root is handled separately */    
    
    /* for each site at the root */
    for (int j = 0; j < nsites; j++) {
        /* compare costs of setting root's j'th site to r */
        int minstate = -1;
        int mincost = inf;
        for (int r = 0; r < nstates; r++) {
            /* TODO note that the state assignment changes depending on
             * whether we use < or <= here, without affecting the total cost */
            if (pscore[n_psubtree_nodes-1][j][r] < mincost) {
                mincost = pscore[n_psubtree_nodes-1][j][r];
                minstate = r;
            } 
        }
        internal_string[n_psubtree_nodes-1][j] = minstate;
    }
        
    /* for the rest of the internal nodes in pre order */
    for (int i = n_psubtree_nodes-2; i >= l; i--) {
        int anc = masked_postorder_reverse[allNodes[masked_postorder[i]]->anc->nodeNum]; 
        /* for each state */
        for (int j = 0; j < nsites; j++) {
            int minstate = -1;
            int mincost = inf;
            for (int r = 0; r < nstates; r++) {
                /* TODO note that the state assignment changes depending on
                * whether we use < or <= here, without affecting the total cost */
                if (W[internal_string[anc][j]][r]+pscore[i][j][r] < mincost) {
                    mincost = W[internal_string[anc][j]][r]+pscore[i][j][r];
                    minstate = r;
                }
            }
            internal_string[i][j] = minstate;
        }
    }

    /* now assign branch lengths to the edges of the p-subtree.
     * we will have to reassign brlens for only the internal branches
     * (those not attached to a leaf). Here internal/external should be
     * understood with respect to the p-subtree. There are p edges,
     * therefor for which we need to reassign branch lengths. */

    for (int i = 0; i < p; i++) {
        int desc = masked_postorder_reverse[pedges[i]]; 
        int anc  = masked_postorder_reverse[allNodes[pedges[i]]->anc->nodeNum];
        double dlen = CalculateHammingDistance(internal_string[anc],
                                               internal_string[desc],
                                               col_count, nsites);
        allNodes[pedges[i]]->dlen = (dlen > min_brlen ? (dlen < max_brlen ? dlen : max_brlen) : min_brlen);
//        SweepDirtynessOverTree(allNodes[pedges[i]]);
  //      allNodes[pedges[i]]->dlen = dlen;
#if 1
        cout << pedges[i] << " " << dlen << "\n";
#endif
    }

    for (int i = 0; i < n_psubtree_nodes; i++) {
        for (int j = 0; j < nsites; j++) {
            free(pscore[i][j]); 
        }
        free(pscore[i]); 
    }
    free(pscore);

    free(masked_postorder);
    free(masked_postorder_reverse);

    for (int i = 0; i < n_psubtree_nodes; i++) {
        free(internal_string[i]);
    }
    free(internal_string);
    return;
}

void  Tree::TranslateTipString(const char *tip, char *string, int nsites) {
    int t = 0;
    for (int i = 0; i < nsites; i++) {
#if 0
        if(tip[t] == '-') {
            cout << tip[t];
            t++;
            int n_ambiguity = tip[t];
            int r = RandomInt(1, n_ambiguity);
            string[i] = tip[t+r];
            t = t+r+1;
        }
#else
        if(tip[t] < 0) {
            if (tip[t] == -4) {
                string[i] = tip[t];
                t++;
            }
            else {
                int n_ambiguity = -tip[t];
                int r = RandomInt(1, n_ambiguity);
                string[i] = tip[t+r];
                t = t+r+1;
            }
        }
#endif
        else {
            string[i] = tip[t];
            t++;
        }
    }
    return;
}

#endif
#endif
