
// GARLI version 0.94 source code
// Copyright  2005 by Derrick J. Zwickl
// All rights reserved.
//
// This code may be used and modified for non-commercial purposes
// but redistribution in any form requires written permission.
// Please contact:
//
//  Derrick Zwickl
//	Integrative Biology, UT
//	1 University Station, C0930
//	Austin, TX  78712
//  email: garli.support@gmail.com
//
//	Note: In 2006  moving to NESCENT (The National
//	Evolutionary Synthesis Center) for a postdoc

// all of the mpi related code appears here or in mpifuncs.cpp

#ifdef MPI_VERSION

#include "defs.h"
#include "threaddcls.h"
#include "mpifuncs.h"
#include "individual.h"

// local vars
transferred_data_t *node_results;
pthread_mutex_t lock_pm;
pthread_mutex_t lock_pop;
pthread_cond_t cond_pm;
bool g_quit_time;
bool g_processing_message;
int remote_types[32];

#define AMR  1
#define SM   2
#define SW	 3

extern int calcCount;

void *thread_func2(void *varg)	{
	int who, size, tag, quits = 0;
	thread_arg_t *targs = (thread_arg_t*)varg;
	
	MasterGamlConfig *conf = targs->conf;
	char *buf;

	bool poo=true;
//	while(poo);
	
	//initialize the array that shows what type of remote each node is
	int method = 0, nprocs = targs->nprocs;
	
	for(who=0;who<nprocs;who++){
		//make all remotes SubtreeWorkers
		remote_types[who]=SW;
		
		/*if (conf->gc.method == "sm" || (conf->gc.method == "hybrid" && who <= (int) (conf->gc.hybridpercent*(nprocs -1))))
			remote_types[who]=SM;
		else if (conf->gc.method == "amr" || (conf->gc.method == "hybrid" && who > conf->gc.hybridpercent*(nprocs -1)))
			remote_types[who]=AMR;
		else
			debug_mpi("ERROR: can't determine method proper type of remote, remote #%d", who);
		*/
		}
	timespec sleepTime;
	sleepTime.tv_nsec=50000000;
	sleepTime.tv_sec=0;
	int nextStart=1, nextRemote;
	int maxNumReceives=2, numReceives;
	while (quits < targs->nprocs-1)	{
		nextRemote=nextStart;
		bool checkedAll=false, received=false;
		do{
			who=nextRemote;
			received=RecvMPIMessage(&buf, &size, who, &tag, false);
			if(received==true) assert(tag==TAG_TREE_STRINGS || tag==TAG_QUIT);
			nextRemote=(nextRemote<(targs->nprocs-1) ? nextRemote+1 : 1);
			if(nextRemote==nextStart) checkedAll=true;
			}while(received==false && checkedAll==false);
		nextStart=(nextStart<(targs->nprocs-1) ? nextStart+1 : 1);
		if(checkedAll==true) nanosleep(&sleepTime, NULL);
		if(received==true){
			if (tag == TAG_QUIT){
				++quits;
				debug_mpi("received quit message from %d.  %d quits received.", who, quits);
				}
			else	{
				pthread_mutex_lock(&lock_pm);
					g_processing_message = true;
				pthread_mutex_unlock(&lock_pm);
			
				pthread_mutex_lock(&lock_pop);
					quits += process_message(buf, size, who, tag, targs);
				pthread_mutex_unlock(&lock_pop);
			
//				if(numReceives == maxNumReceives){
					pthread_mutex_lock(&lock_pm);
						g_processing_message = false;
						pthread_cond_signal(&cond_pm);
					pthread_mutex_unlock(&lock_pm);
//					numReceives=0;
//					}
//				else numReceives++;
				}
			delete []buf;
			numReceives++;
			}
		
		if(g_quit_time == true){
			send_quit_messages(nprocs);
			break;
			}

/*		if((checkedAll==true && g_processing_message==true) || (received==true && numReceives >=  maxNumReceives)){
			pthread_mutex_lock(&lock_pm);
				g_processing_message = false;
				pthread_cond_signal(&cond_pm);
			pthread_mutex_unlock(&lock_pm);
			numReceives=0;
			}*/
//		else numReceives++;
		}
	debug_mpi("thread terminating");
	g_quit_time = true;
}

void send_quit_messages(int np){
	for(int i=1;i<np;i++){
		SendMPIMessage(NULL, 0, i, TAG_QUIT);
		}
	}

int process_message(char *buf, int size, int who, int tag, thread_arg_t *targ)	{
	MasterGamlConfig *conf = targ->conf;
	Population *pop = targ->pop;
	int method = 0, nprocs = targ->nprocs;
	int foundQuit=0;

	if(remote_types[who]==SM) DoMasterSM(buf, size, who, tag, targ);
	else if(remote_types[who]==AMR) DoMasterAMR(buf, size, who, tag, targ);

	else if(remote_types[who]==SW) foundQuit=DoMasterSW(buf, size, who, tag, targ);

	else debug_mpi("ERROR: can't determine method to use in process_message(), remote #%d", who);

	return foundQuit;
}

void DoMasterSM(char *buf, int size, int who, int tag, thread_arg_t *targ)	{
	MasterGamlConfig *conf = targ->conf;
	Population *pop = targ->pop;
	int count, start_shield, *which = new int;//[conf->gc.numshields];
	char *tree_strings;
	double *models;

	assert(tag == TAG_TREE_STRINGS);
	tree_strings = buf;
	count = CountTreeStrings(tree_strings);
	
	RecvMPIMessage(&buf, &size, who, &tag, true);
	assert(tag == TAG_MODEL);
	models=(double*)buf;

	// print some messages
	debug_mpi("SYNCHRONOUS COMMUNICATION (%d, SM)", who);
	debug_mpi("\trecv: %d tree strings", count);
	debug_mpi("\trecv: %d models", count);
	
	bool poo=true;
	//while(poo);
	
	*which = start_shield = pop->params->nindivs + (who-1);
//	for (int i = 0; i < conf->gc.numshields; ++i)
//		which[i] = start_shield++;
	pop->ReplaceSpecifiedIndividuals(count, which, tree_strings, models);
	pop->CalcAverageFitness();
	
	delete [] which;
	//this will be deleted back where the initial call to RecvMPIMessage was made
	//delete [] tree_strings;
	delete [] (char*)models;
	
	//under certain conditions, tell the remote to become an AMR node
	if(((double)rand()/RAND_MAX)<.002){
		debug_mpi("sent message to change to AMR to remote %d", who);
		SendMPIMessage(NULL, 0, who, TAG_REMOTE_TYPE_SWITCH);
		//get rid of any SM messages that might be waiting
		while(RecvMPIMessage(&buf, &size, who, &tag, false)==true);
		//update the remote_types array
		remote_types[who]=AMR;
	}
}

void DoMasterAMR(char *buf, int size, int who, int tag, thread_arg_t *targ)	{
	MasterGamlConfig *conf = targ->conf;
	Population *pop = targ->pop;
	int count, start_shield, which;
	char *tree_strings, *model_buf;
	double *models, score;

	assert(tag == TAG_SCORE);
	memcpy(&score, buf, sizeof(double));
	//delete [] buf;

	// print some messages
	debug_mpi("SYNCHRONOUS COMM (%d, AMR, %f)", who, pop->bestFitness);
	debug_mpi("\trecv: score %f", score);

	if (score > pop->bestFitness)	{
		SendMPIMessage(NULL, 0, who, TAG_TREE_STRINGS_REQUEST);
		debug_mpi("\tsent: TAG_TREE_STRINGS_REQUEST");
		RecvMPIMessage(&buf, &size, who, &tag, true);
		assert(tag == TAG_TREE_STRINGS);
		tree_strings = buf;
		debug_mpi("\trecv: %d tree strings", CountTreeStrings(tree_strings));

		RecvMPIMessage(&buf, &size, who, &tag, true);
		assert(tag == TAG_MODEL);
		models=(double*)buf;
		
		which = pop->total_size-1;
		pop->ReplaceSpecifiedIndividuals(1, &which, tree_strings, models);
		pop->CalcAverageFitness();
		debug_mpi("score sent: %f, score calced: %f", score, pop->IndivFitness(which));
		//assert(abs(score -  pop->IndivFitness(which))<.00001);
		}
	else	{
		which = (int)pop->cumfit[pop->total_size-1][0];
		pop->GetSpecifiedTreeStrings(&tree_strings, 1, &which);
		int model_size=pop->GetSpecifiedModels(&models, 1, &which);

		SendMPIMessage(tree_strings, strlen2(tree_strings)+2, who, TAG_TREE_STRINGS);
		debug_mpi("\tsent: %d tree strings", CountTreeStrings(tree_strings));

		model_buf = (char*)models;
		SendMPIMessage(model_buf, sizeof(double)*model_size, who, TAG_MODEL);
		debug_mpi("\tsent: %d models", 1);
		}
	delete [] tree_strings;
	delete [] models;
	}

int DoMasterSW(char *buf, int size, int who, int tag, thread_arg_t *targ)	{
	MasterGamlConfig *conf = targ->conf;
	Population *pop = targ->pop;
	int count, start_shield, *which = new int;//[conf->gc.numshields];
	char *tree_strings, *model_buf;
	char *out_tree_strings;
	char *defbuf, *subbuf;
	double *out_models;
	double *models;
	int remoteSubtreeDef, remoteSubtreeNode;
	ParallelManager *paraMan=(pop->paraMan);

//first get the tree and model from the remote and include it in
//the master population.  If there are multiple messages, chuck 
//earlier ones and just get the most recent
	bool firstmessage=true;
	do{
		debug_mpi("Remote %d", who);
		assert(tag == TAG_TREE_STRINGS || tag == TAG_QUIT);
		if(tag==TAG_QUIT) return 1;
		tree_strings = buf;
		count = CountTreeStrings(tree_strings);
		
//		debug_mpi("about to get model strings...");

		RecvMPIMessage(&buf, &size, who, &tag, true);
		assert(tag == TAG_MODEL);
		models=(double*)buf;
		
//		debug_mpi("about to get subdef strings...");

		//determine what the remote was doing when it sent this tree
		RecvMPIMessage(&defbuf, &size, who, &tag, true);
		assert(tag==TAG_SUBTREE_ITERATION);
		remoteSubtreeDef=atoi(defbuf);
		
		if(remoteSubtreeDef>0){
			RecvMPIMessage(&subbuf, &size, who, &tag, true);
			assert(tag==TAG_SUBTREE_DEFINE);
			remoteSubtreeNode=atoi(subbuf);
			if(remoteSubtreeDef==paraMan->subtreeDefNumber)
				paraMan->localSubtreeAssign[who]=remoteSubtreeNode;
			else paraMan->localSubtreeAssign[who]=0;
			delete []subbuf;
			}
		//DJZ 5-18-05
		else {
			paraMan->localSubtreeAssign[who]=0;
			remoteSubtreeNode=0;
			}
		
		double score;
		char *scoreBuf;
		RecvMPIMessage(&scoreBuf, &size, who, &tag, true);
		assert(tag==TAG_SCORE);
		memcpy(&score, scoreBuf, sizeof(double));
//		debug_mpi("recieved score of %f", score);
		delete []scoreBuf;
		
		if(firstmessage==false) debug_mpi("\tfound another tree from remote %d", who);

		*which = start_shield = pop->params->nindivs + (who-1);
		pop->ReplaceSpecifiedIndividuals(count, which, tree_strings, models);
	
		pop->indiv[*which].SetFitness(score);
		if(firstmessage==false) delete []tree_strings;
		delete [](char*)models;
		delete []defbuf;

		firstmessage=false;
		}while(RecvMPIMessage(&buf, &size, who, &tag, false)==true);
	
	bool subtreesCurrent = ((remoteSubtreeDef == paraMan->subtreeDefNumber) && remoteSubtreeDef > 0);
	
	if(paraMan->subtreeModeActive==false || subtreesCurrent==false){
		pop->indiv[*which].accurateSubtrees=false;
		pop->newindiv[*which].accurateSubtrees=false;
		}
	else {
		pop->indiv[*which].accurateSubtrees=true;
		pop->newindiv[*which].accurateSubtrees=true;
		}

//	debug_mpi("about to CalcFitness...");
	double prevBestScore=pop->BestFitness();
//	pop->indiv[*which].CalcFitness(0);
	pop->indiv[*which].treeStruct->calcs=calcCount;
	pop->CalcAverageFitness();

	//reclaim clas if the new tree has essentially no chance of reproducing
	if(((pop->indiv[*which].Fitness() - pop->indiv[pop->bestIndiv].Fitness()) < (-11.5/pop->params->selectionIntensity))){
//		debug_mpi("about to reclaim...");
		pop->indiv[*which].treeStruct->ReclaimUniqueClas();
		}
	

	//Now, take a look at what we got from the remote and decide what to do
	double inscore=pop->indiv[*which].Fitness();
	double scorediff=prevBestScore - inscore;
	debug_mpi("\tnew ind - def %d - node %d - lnL: %f", remoteSubtreeDef, remoteSubtreeNode, inscore);
	if(scorediff < 0) debug_mpi("\tPrev Best=%f, diff=%f (new best)", prevBestScore, scorediff);
	else debug_mpi("\tPrev Best=%f, diff=%f", prevBestScore, scorediff);
//	debug_mpi("\tbest=%d, bestAc=%d, bestlnL=%f, bestAcclnL=%f", pop->bestIndiv, pop->bestAccurateIndiv, pop->BestFitness(), pop->indiv[pop->bestAccurateIndiv].Fitness());

	
	bool recalcSubtrees=false;
	if(scorediff < -0.01){
		pop->LogNewBestFromRemote(-scorediff, *which);
		}

	int subtreeNum;
	bool send=false;
	
	//there are really 8 possible cases here
	//1. Subtree mode active, 	got accurate tree,	score good -> do nothing
	//2. 											score bad  -> send best accurate tree
	//3. 						inaccurate tree,	score good -> recalc subtrees, send?
	//4.											score bad  -> send best accurate tree
	//5. Subtree mode inactive, got accurate tree, 	score good -> send best tree
	//6.											score bad  -> send best tree
	//7.						inaccurate tree, 	score good -> do nothing
	//8.											score bad  -> send best tree
	//so, 2 "do nothings" 3 "send best", 2 "send best accurate" and 1 "subtree recalc"
	
//if subtree mode isn't active, send the remote our best tree if the
//tree we got from it is worse by some amount, or if it is still working
//on a subtree
	
	double updateThresh=paraMan->updateThresh;
	if(paraMan->subtreeModeActive==false){
		if((paraMan->perturbModeActive==false && (scorediff > updateThresh || remoteSubtreeDef>0))/* || (paraMan->needToSend[who]==true)*/){
			debug_mpi("\tupdate thresh = %f, send indiv", updateThresh);
			//cases 5, 6 and 8
			*which = (int)pop->cumfit[pop->total_size-1][0];
			subtreeNum=0;
			send=true;
			}
		else debug_mpi("\tupdate thresh = %f", updateThresh);
		}
		
	else if(paraMan->subtreeModeActive==true){
		//cases 1-4
		if((scorediff > updateThresh) || (subtreesCurrent==false)/* || paraMan->perturb==true*/){
			//cases 2 and 4.  send the best accurate tree
			*which=pop->bestAccurateIndiv;	
			if(paraMan->remoteSubtreeAssign[who] != 0) subtreeNum=paraMan->remoteSubtreeAssign[who];
			else subtreeNum=paraMan->ChooseSubtree();
			debug_mpi("\tsend best accurate ind, %f (best=%f)", pop->indiv[*which].Fitness(), pop->bestFitness);
//			debug_mpi("\tperturb=%d, bestFit=%f, indFit=%f", paraMan->perturb, pop->bestFitness, pop->indiv[*which].Fitness());
			send=true;
			}
		else if(recalcSubtrees==true && subtreesCurrent==false){
			//case 3
			//if the new inaccurate tree that came in is better than what we have,
			//recalcuate the subtrees, and send the same tree back, but with a 
			//subtree asignment
			pop->StartSubtreeMode();
			debug_mpi("Recalculating subtrees");
			subtreeNum=paraMan->ChooseSubtree();
			send=true;
			}
		}

	if(paraMan->needToSend[who]){
		char pertbuf[5];
		int perttype = (pop->pertMan->pertType > 0 ? pop->pertMan->pertType : (int)(rnd.uniform() * 2 + 1));
		sprintf(pertbuf, "%d", perttype);
		SendMPIMessage(pertbuf, strlen(pertbuf)+2, who, TAG_PERTURB);
		debug_mpi("sending pertub message to %d, type %d", who, perttype);
		paraMan->needToSend[who]=false;
		}


	if(send==true){
		pop->GetSpecifiedTreeStrings(&out_tree_strings, 1, which);
		
		assert(*out_tree_strings == '(');
		int model_size=pop->GetSpecifiedModels(&out_models, 1, which);
		
		SendMPIMessage(out_tree_strings, strlen2(out_tree_strings)+2, who, TAG_TREE_STRINGS);
		SendMPIMessage((char*)out_models, sizeof(double)*model_size, who, TAG_MODEL);

/*		if(paraMan->needToSend[who]){
			char pertbuf[5];
			int perttype = (pop->pertMan->pertType > 0 ? pop->pertMan->pertType : (rnd.uniform() * 2 + 1));
			sprintf(pertbuf, "%d", subtreeNum);
			SendMPIMessage(NULL, 0, who, TAG_PERTURB);
			debug_mpi("sending pertub message to %d, type %d", who, perttype);
			paraMan->needToSend[who]=false;
			}
*/			
//		else{
			char stn[5];
			sprintf(stn, "%d", subtreeNum);
			SendMPIMessage(stn, strlen(stn)+2, who, TAG_SUBTREE_DEFINE);
			debug_mpi("\tsent ind %d, lnL %f", *which, pop->indiv[*which].Fitness());

			if(subtreeNum > 0){
				//if this node was already assigned a subtree, be sure to subtract the old one from the assigned array
				sprintf(stn, "%d", paraMan->subtreeDefNumber);
				debug_mpi("\tsubdef %d, node %d", paraMan->subtreeDefNumber, subtreeNum);
				SendMPIMessage(stn, strlen(stn)+2, who, TAG_SUBTREE_ITERATION);
				}
			
//			}
		paraMan->remoteSubtreeAssign[who]=subtreeNum;
		
		delete []out_models;
		delete []out_tree_strings;		
		}
#ifndef NDEBUG
	if(paraMan->subtreeModeActive && paraMan->subtreeDefNumber==remoteSubtreeDef){
        //if we think that this remote gave us a tree with accurate subtrees, check
        paraMan->CheckSubtreeAccuracy(pop->indiv[which[0]].treeStruct);
		}
#endif
	
	//the tree_strings that were passed in will be deleted back
	//where the initial call to RecvMPIMessage was made
	delete [] which;

	pop->CalcAverageFitness();
	return 0;
}

void purge_results(transferred_data_t *r)	{
	if (r->tree_strings)
		delete [] r->tree_strings;
	if (r->kappas)
		delete [] r->kappas;
	memset(r, 0, sizeof(transferred_data_t));
}

void copy_results(transferred_data_t *lhs, transferred_data_t rhs)	{
	memcpy(lhs, &rhs, sizeof(transferred_data_t));
	if (rhs.tree_strings)	{
		lhs->tree_strings = new char[rhs.ts_size];
		memcpy(lhs->tree_strings, rhs.tree_strings, rhs.ts_size);
	}
	if (rhs.kappas)	{
		lhs->kappas = new double[rhs.k_size/sizeof(double)];
		memcpy(lhs->kappas, rhs.kappas, rhs.k_size);
	}
}

bool valid_results(transferred_data_t r)	{
	if (r.tree_strings && r.kappas)
		return true;
	if (r.tag == TAG_SCORE)
		return true;
	return false;
}

#endif