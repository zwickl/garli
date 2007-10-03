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
//  email: zwickl@mail.utexas.edu
//
//	Note: In 2006  moving to NESCENT (The National
//	Evolutionary Synthesis Center) for a postdoc

//	NOTE: Portions of this source adapted from GAML source, written by Paul O. Lewis

// all of the mpi related code appears here or in threadfuncs.cpp


#ifdef MPI_VERSION

#include <stdarg.h>
#include <mpi.h>
#include <string.h>

#include "mpifuncs.h"
#include "defs.h"
#include "population.h"
#include "individual.h"
#include "configoptions.h"
#include "configreader.h"
#include "funcs.h"
#include "stopwatch.h"
#include "threaddcls.h"
#include "adaptation.h"
#include "errorexception.h"

// globals
Stopwatch *g_sw=NULL;
long int* g_gen = NULL;
extern rng rnd;

FILE *fhandle;

int MPIMain(int argc, char** argv)	{

	MPI_Init(&argc, &argv);
	
	int rank, nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	bool poo=true;
	//if(rank==0) while (poo) ;

	try{
	if (rank == 0)	{
        	MasterGamlConfig conf;
	        int err=conf.Read("garli.conf", true);
		if(err != 0){
			send_quit_messages(nprocs);
			throw ErrorException("Error in config file (Master)...aborting.");
			}
		
		LogConfig(conf);

	        // Create the data object
        	HKYData data;
        	ReadData(conf.datafname.c_str(), &data);

	        // start the remote nodes going...
        	StartProcs(conf, data);

	        // start yourself!
        	MasterMaster(conf, data);

		}
	else	{ // rank != 0
	
		int from, tag, size;
		char* buf;
		HKYData data;
		GeneralGamlConfig conf;
		if (conf.Read("garli.conf") < 0)	{
			throw ErrorException("Error in config file (Remote)...aborting.");
			}

		//no longer sending the conf, just letting the remote read it from file
//		for (int i = 0; i < 2; ++i)	{
			RecvMPIMessage(&buf, &size, &from, &tag, true);
			assert(from == 0); // sanity check
			if (tag == TAG_DATA)	{
				data.Deserialize(buf, size);
				debug_mpi("receieved data from %d", from);
				}
/*			else if (tag == TAG_CONFIG)	{
				conf.Deserialize(buf, size);
				debug_mpi("received conf from %d", from);
				}
*/			else	{
				debug_mpi("ERROR: received unexpected message from %d with tag %d", from, tag);
				debug_mpi("aborting from MPIMain()");
				}
			delete [] buf;
//		}
		
//		LogConfig(conf);
		RemoteMaster(conf, data);
		}
	}catch(ErrorException err){
		err.Print(cout);
		}
	// time to kill some global vars
	delete [] node_results;
	
	MPI_Finalize();
	return 0;
}

int StartProcs(const GeneralGamlConfig& conf, HKYData& data)	{

//	debug_mpi("entering StartProcs()");

	char* conf_buf, *data_buf;
	int conf_size, data_size;
	GeneralGamlConfig ctest;
//	HKYData dtest;

	data.Serialize(&data_buf, &data_size);
//	conf.Serialize(&conf_buf, &conf_size);
	
	// sanity check: make sure the serialization code works
/*	dtest.Deserialize(data_buf, data_size);
	assert(data == dtest);
	ctest.Deserialize(conf_buf, conf_size);
	assert(conf == ctest);
*/	
	int nprocs;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	for (int i = 1; i < nprocs; ++i)	{
		SendMPIMessage(data_buf, data_size, i, TAG_DATA);
		debug_mpi("sent data to node %d", i);
//		SendMPIMessage(conf_buf, conf_size, i, TAG_CONFIG);
//		debug_mpi("sent conf to node %d", i);
	}
	
	delete [] data_buf;
//	delete [] conf_buf;
	
	debug_mpi("leaving StartProcs()");
	
	return 0;
}

/* threaded MasterMaster */
/* i would prefer that the thread initialization code happen in MPIMain(), but
 * it needs Population pop, which is declared here */
int MasterMaster(MasterGamlConfig& conf, HKYData& data)	{
	Parameters params;
	params.SetParams(conf, data);
	LogParams(params);

	int nprocs;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	bool poo=true;
//	while(poo);
		
	Tree::alpha = params.gammaShapeBrlen;
	Tree::meanBrlenMuts = params.meanBrlenMuts;

	Population pop;
//	debug_mpi("about to setup");
	pop.Setup(params, &conf, nprocs, 0);
	g_sw=&pop.stopwatch;
//	debug_mpi("setup");
	g_gen = &pop.gen;

	pop.CalcAverageFitness();

	// start the thread
	pthread_t thread;
	thread_arg_t targ;
	pthread_mutex_init(&lock_pm, NULL);
	pthread_mutex_init(&lock_pop, NULL);
	pthread_cond_init(&cond_pm, NULL);
	g_quit_time = false;
	g_processing_message = false;
	targ.conf = &const_cast<MasterGamlConfig&>(conf);
	targ.pop = &pop;
	targ.nprocs = nprocs;
	
	pthread_create(&thread, NULL, thread_func2, (void*)&targ);
	
	cout << "Master running..." << endl;

	pop.gen=0;
	while (!g_quit_time){
		pthread_mutex_lock(&lock_pop);
			pop.keepTrack();
			pop.OutputFate();
			if (pop.gen % conf.logevery == 0) pop.OutputLog();
			++pop.gen;
			pop.NextGeneration();
			if(pop.gen % pop.params->saveEvery == 0) pop.CreateTreeFile( pop.params->treefname );
			if(pop.gen % pop.adap->intervalLength == 0){
	            	bool reduced=false;
	    	        if(pop.gen-pop.lastTopoImprove >= pop.adap->intervalsToStore*pop.adap->intervalLength){
	         	           reduced=pop.adap->ReducePrecision();
	            	        }
	    	        if(reduced){
	            	        pop.lastTopoImprove=pop.gen;
	                	   pop.indiv[pop.bestIndiv].treeStruct->OptimizeAllBranches(pop.adap->branchOptPrecision);
	                   	   pop.indiv[pop.bestIndiv].SetDirty();
	                   	   pop.CalcAverageFitness();
							//DJZ 2/20/06
							//reducing parallel remote update thresh based on same criteria as opt precision
							double prev=pop.paraMan->updateThresh;
							pop.paraMan->ReduceUpdateThresh();
							debug_mpi("Remote update threshold reduced from %f to %f", prev, pop.paraMan->updateThresh);	                    	
							}
/*	           		 else if(!(pop.gen%(pop.adap->intervalLength*pop.adap->intervalsToStore))){
	                            pop.indiv[pop.bestIndiv].treeStruct->OptimizeAllBranches(pop.adap->branchOptPrecision);
	                            pop.indiv[pop.bestIndiv].SetDirty();
	                            pop.CalcAverageFitness();
	                            }
*/
				if(pop.enforceTermConditions == true
					&& pop.gen-pop.lastTopoImprove > pop.lastTopoImproveThresh 
					&& pop.adap->improveOverStoredIntervals < pop.improveOverStoredIntervalsThresh
					&& pop.adap->branchOptPrecision == pop.adap->minOptPrecision){
			//		&& pop.paraMan->updateThresh == pop.paraMan->minUpdateThresh){
					cout << "Reached termination condition!\nlast topological improvement at gen " << pop.lastTopoImprove << endl;
					cout << "Improvement over last " << pop.adap->intervalsToStore*pop.adap->intervalLength << " gen = " << pop.adap->improveOverStoredIntervals << endl;
					g_quit_time=true;
					break;
					}
				pop.CheckSubtrees();

#ifdef INCLUDE_PERTURBATION
				pop.CheckPerturbParallel();
#endif
				}
		pthread_mutex_unlock(&lock_pop);
		pthread_mutex_lock(&lock_pm);
			while (g_processing_message)
				pthread_cond_wait(&cond_pm, &lock_pm);
		pthread_mutex_unlock(&lock_pm);
		}
	
	//DJZ 3-1-06 Need to give control back to the thread one more time so that it can deal with any final messages
	pthread_mutex_unlock(&lock_pop);
	pthread_mutex_lock(&lock_pm);
		while (g_processing_message)
			pthread_cond_wait(&cond_pm, &lock_pm);
	pthread_mutex_unlock(&lock_pm);
		
	pop.FinalOptimization();
	pop.FinalizeOutputStreams();
	pthread_join(thread, NULL);
	return 0;
}

/* old MasterMaster 
int MasterMaster(const GamlConfig& conf, HKYData& data)	{
	Parameters params;
	params.SetParams(conf, data, true);
	LogParams(params);

	Tree::brlen_mu       = params.brlenMutProb;
	Tree::mu             = params.topoMutProb;
	Tree::lambda         = params.crossoverProb;
	Tree::alpha          = params.gammaShape;
	
	int nprocs;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	
	Population pop;
	pop.Setup(params, conf, nprocs, 0);
	g_gen = &pop.gen;
	
	for (int i = 1; i < nprocs; ++i)	{
		if (i < conf.hybridp.nt*nprocs)
			debug_mpi("node %d is shielded migrants", i);
		else
			debug_mpi("node %d is alpha male replication", i);
	}
	
	if (conf.method == "fde")
		MasterFullDuplexExchange(pop, conf);
	else if (conf.method == "sm")
		MasterShieldedMigrants(pop, conf);
	else if (conf.method == "amr")
		MasterAlphaMaleReplication(pop, conf);
	else if (conf.method == "hybrid")
		MasterHybrid(pop, conf);
	else	{
		debug_mpi("ERROR: unknown method (GamlConfig::General::method): %s", conf.method.c_str());
		MPI_Abort(MPI_COMM_WORLD, -1);
	}

	return 0;
}
*/


int RemoteMaster(GeneralGamlConfig& conf, HKYData& data)	{

	debug_mpi("starting RemoteMaster()...");
	
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	Parameters params;
	params.SetParams(conf, data);
	//the seed has already been set in SetParams above, but we need to
	//modify it so that all remotes are not the same
	rnd.set_seed((rank+1) * rnd.init_seed());
	LogParams(params);


	Tree::alpha = params.gammaShapeBrlen;
	Tree::meanBrlenMuts = params.meanBrlenMuts;
	
	int nprocs;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	
	bool poo=true;
//	while(poo);
	
	Population pop;
	pop.Setup(params, &conf, nprocs, rank);
	g_sw=&pop.stopwatch;
	g_gen = &pop.gen;

	
	//for now all nodes will be SW
	debug_mpi("doing remote subtree worker");
	RemoteSubtreeWorker(pop, conf);
	/*
	if (conf.method == "fde")
		assert(0);
//		RemoteFullDuplexExchange(pop, conf);
	//DJZ changed this 9/10/03 to avoid occasional disparity in the acutal # of sm nodes and the number of shields allocated
	else if (conf.method == "sm" || (conf.method == "hybrid" && rank <= (int) (conf.hybridp.nt*(nprocs-1)))){
//	else if (conf.method == "sm" || (conf.method == "hybrid" && rank < conf.hybridp.nt*nprocs)){
		debug_mpi("doing remote shielded migrants");
		RemoteShieldedMigrants(pop, conf);
	}
	else if (conf.method == "amr" || conf.method == "hybrid")	{
		debug_mpi("doing remote alpha male replication");
		RemoteAlphaMaleReplication(pop, conf);
	}
	else	{
		debug_mpi("ERROR: unknown method (GamlConfig::General::method): %s", conf.method.c_str());
		MPI_Abort(MPI_COMM_WORLD, -1);
	}
	*/
	return 0;
}

/*
int RemoteMaster(const GamlConfig& conf, HKYData& data)	{

	debug_mpi("starting RemoteMaster()...");
	
	Parameters params;
	params.SetParams(conf, data, false);
	LogParams(params);
	
	Tree::brlen_mu       = params.brlenMutProb;
	Tree::mu             = params.topoMutProb;
	Tree::lambda         = params.crossoverProb;
	Tree::alpha          = params.gammaShape;

	Population pop;
	pop.Setup(params, conf.max_nindivs);
	
	char* tree_strings_in;
	char* tree_strings_out;
	int trans_count = 0;

	// tree's must be scored before calling next generation so call it here so it doesn't crash on gen == 1
	pop.CalcAverageFitness();
	
	pop.gen = 1;
	while (true)	{
		if (pop.gen % conf.interval == 0)	{
		
			debug_mpi("send interval reached.");

        	// sanity check: make sure param's nindivs == pop's current/original size
      		assert( (params.nindivs == pop.current_size) && (pop.current_size == pop.original_size) );

        	// see if there is a quit message
        	int flag;
			MPI_Status status;
			MPI_Iprobe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
			if (flag)	{
				if (status.MPI_TAG == TAG_QUIT)	{
					debug_mpi("quit message received.");
					break;
				}
				else // sanity check: should not get here
					assert(status.MPI_TAG == TAG_QUIT);
			}
			
			// save old scores for TransLog()
			pop.CalcAverageFitness();
			double* old_scores = new double[params.nindivs];
			for (int j = 0; j < params.nindivs; ++j)
				old_scores[j] = pop.indiv[j].Fitness();

			// send tree strings to node
			pop.ShrinkPopulation(params.nindivs, &tree_strings_out); // allocates buf on heap
			SendResultsToNode(0, pop.original_size, tree_strings_out);
			debug_mpi("sent %d tree strings to master...", pop.original_size);
			char* p = tree_strings_out;
			for (int i = 0; i < pop.original_size; ++i)	{
				debug_mpi("%s", p);
				p += strlen(p) +1;
			}

			// wait for tree's to be sent back.
			// must probe this incase a quit message is sent
			int nindivs;
			MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			if (status.MPI_SOURCE != 0)
				assert(status.MPI_SOURCE == 0); // should not receive messages from any node except the master
			if (status.MPI_TAG == TAG_QUIT)	{
				debug_mpi("quit message received while waiting for tree strings.");
				break;
			}
			else if (status.MPI_TAG == TAG_TREE_STRINGS_COUNT)	{
			
				GetResultsFromNode(0, &nindivs, &tree_strings_in);
				debug_mpi("received %d tree strings from master...", nindivs);
				char* p = tree_strings_in;
				for (int i = 0; i < pop.original_size; ++i)	{
					debug_mpi("%s", p);
					p += strlen(p) +1;
				}
	
				// sanity check: make sure we're getting back as many as we sent out
				assert(nindivs == pop.original_size);

				pop.ExtendPopulation(nindivs, tree_strings_in);
				pop.CalcAverageFitness();
				
				// save new scores for TransLog()
				double* new_scores = new double[params.nindivs];
				for (int j = 0; j < params.nindivs; ++j)
					new_scores[j] = pop.indiv[j].Fitness();
				
				// remote nodes send all indivs and replace all indivs
				int* temp = new int[params.nindivs+1];
				for (int j = 0; j < params.nindivs; ++j)
					temp[j] = j;
					
				TransLog(trans_count++, params.nindivs, nindivs, tree_strings_out,
						 nindivs, tree_strings_in,
						 temp, temp,
						 old_scores, new_scores);
					
				delete [] tree_strings_in;
				delete [] tree_strings_out;
				delete [] new_scores;
				delete [] old_scores;
				delete [] temp;
				
			}
			else // sanity check: should not get here!
				assert(status.MPI_TAG != TAG_QUIT || status.MPI_TAG != TAG_TREE_STRINGS_COUNT);
		}
		else	{
			pop.NextGeneration();
			++(pop.gen);
		}
	} // end while (true)

	return 0;
}
*/
/*	returns:	how many nodes have results.
	pre condition:	nodes[] must be of size nprocs-1
	post condition:	nodes[] will hold >which< nodes have results.
	example:	if remote nodes 1, 2 and 5 have results, then
				return == 3
				nodes[] == {1, 2, 5, (garbage)...}
*/
int PollForResults(int nodes[])	{
	int nprocs, flag, count = 0;
	MPI_Status status;

	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	for (int i = 2 ; i < nprocs; ++i)	{
		MPI_Iprobe(i, TAG_TREE_STRINGS_COUNT, MPI_COMM_WORLD, &flag, &status);
		if (flag)
			nodes[count++] = status.MPI_SOURCE;
	}

	return count;
}

bool PollForResults(int n)	{
	int flag;
	MPI_Status status;
	MPI_Iprobe(n, TAG_TREE_STRINGS_COUNT, MPI_COMM_WORLD, &flag, &status);
	if (flag)
		return 1;
	else
		return 0;
}
	

/*	results format:	4 bytes (int) - nindivs
					4 bytes (int) - size of tree strings including null terminators
					n bytes (char) - tree strings seperated by NULLs, terminated by a double NULL
	returns:		which node it received from
	side affects:	if returns 0, then tree_strings is allocated on the heap
*/
int GetResultsFromNode(int node, int* n_, char** tree_strings_)	{
	int& n = *n_;
	char*& tree_strings = *tree_strings_;
	int buf_size = 0;
	MPI_Status status;

	// sanity check: make sure node_num is a valid node
	int nprocs;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	assert( ((node >= 0) && (node < nprocs)) || (node == MPI_ANY_SOURCE) );
	
	MPI_Recv(&n, 1, MPI_INT, node, TAG_TREE_STRINGS_COUNT, MPI_COMM_WORLD, &status);
	assert(n > 0);
	if (node == MPI_ANY_SOURCE)
		node = status.MPI_SOURCE;
	MPI_Recv(&buf_size, 1, MPI_INT, node, TAG_TREE_STRINGS_SIZE, MPI_COMM_WORLD, &status);
	assert(buf_size > 0);
	tree_strings = new char[buf_size];
	MPI_Recv(tree_strings, buf_size, MPI_CHAR, node, TAG_TREE_STRINGS, MPI_COMM_WORLD, &status);

	// sanity check: the number of tree strings in tree_strings should == nindivs
	int count = 0;
	char* p = tree_strings;
	while (*p)	{
		p += strlen(p) + 1;
		++count;
	}
	assert(count == n);
	
	// sanity check: make sure p ended up at the end of the string
	++p;
	assert(p-tree_strings == buf_size);
	
	return status.MPI_SOURCE;
}

/*	results format:	4 bytes (int) - nindivs
					4 bytes (int) - size of tree strings including null terminators
					n bytes (char) - tree strings seperated by NULLs, terminated by a double NULL
	pre conditions:	0 < node < nprocs
					n > 0
					tree_strings must have n tree_strings seperated by NULLs
*/
int SendResultsToNode(int node, int n, char* tree_strings)	{

	// sanity check: make sure node is a valid number
	int nprocs;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	assert( (node >= 0) && (node < nprocs) );

	// sanity check: make sure the actual number of tree strings in tree_strings is equal to n
	int count = 0;
	char* p = tree_strings;
	while (*p)	{
		p += strlen(p) + 1;
		++count;
	}
	assert(count == n);

	int buf_size = (p+1)-tree_strings;

	// sanity check: make sure i calculated buf_size correctly
	p = tree_strings;
	count = 0;
	while (*p || *(p+1))	{	// this is basically strlen() for string terminated by two NULLs instead of one
		++count;
		++p;
	}
	assert(count+2 == buf_size);

	MPI_Status status;
	MPI_Send(&n, 1, MPI_INT, node, TAG_TREE_STRINGS_COUNT, MPI_COMM_WORLD);
	MPI_Send(&buf_size, 1, MPI_INT, node, TAG_TREE_STRINGS_SIZE, MPI_COMM_WORLD);
	MPI_Send(tree_strings, buf_size, MPI_CHAR, node, TAG_TREE_STRINGS, MPI_COMM_WORLD);

	return 0;
}

int ReceiveParams(Parameters* params_, int node)	{
	Parameters& params = *params_;

    int size;
    MPI_Status status;
	
	MPI_Recv(&size, 1, MPI_INT, node, TAG_PARAMS_SIZE, MPI_COMM_WORLD, &status);
	char* buf = new char[size];
	MPI_Recv(buf, size, MPI_CHAR, node, TAG_PARAMS, MPI_COMM_WORLD, &status);

	params.Deserialize(buf);

	delete [] buf;

	return 0;
}

int ReceiveData(HKYData* data_, int node)	{
	HKYData& data = *data_;

    int size;
    MPI_Status status;

	MPI_Recv(&size, 1, MPI_INT, node, TAG_DATA_SIZE, MPI_COMM_WORLD, &status);
	char* buf = new char[size];
	MPI_Recv(buf, size, MPI_CHAR, node, TAG_DATA, MPI_COMM_WORLD, &status);

	data.Deserialize(buf, size);

	delete [] buf;

	return 0;
}

int debug_mpi(const char* fmt, ...)	{
	static bool first_call = true;
	static int rank;
	static char fname[13];

	if (first_call)	{
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if (rank < 10)
			sprintf(fname, "node0%d.log", rank);
		else
			sprintf(fname, "node%d.log", rank);
		fhandle = fopen(fname, "w");
		first_call = false;
	}
	if(g_sw != NULL)
		fprintf(fhandle, "g %d, t %d: ", (g_gen && *g_gen != 0 ? *g_gen : -1), g_sw->SplitTime());
	else 
		fprintf(fhandle, "g %d, t %d: ", (g_gen && *g_gen != 0 ? *g_gen : -1), 0);
	va_list vl;
	va_start(vl, fmt);
	vfprintf(fhandle, fmt, vl);
	va_end(vl);
	fprintf(fhandle, "\n");

	fflush(fhandle);

	return 0;
}

int TransLog(int to_who_or_count, int nindivs, int n, char* str_out, int m, char* str_in, int* to_send, int* to_replace, double* old_scores, double* new_scores)	{

	char fname[64];
	char temp_buf[64];
	int rank;
	static int* counts;
	static bool first_call = true;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	if (first_call && rank == 0)	{
		int nprocs;
		MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
		counts = new int[nprocs];  // yeah yeah, so i'm never deallocating this.  it'll be free'd when the process ends....=)
		memset(counts, 0, sizeof(int)*nprocs);
	}
	
	if (rank < 10)
		sprintf(fname, "trans0%d.log", rank);
	else
		sprintf(fname, "trans%d.log", rank);
	
	FILE* translog;
	if (first_call)	{
		translog = fopen(fname, "w");
		first_call = false;
	}
	else
		translog = fopen(fname, "a");
	assert(translog);
	
	if (rank == 0)
		sprintf(temp_buf, "%d-%d", to_who_or_count, counts[to_who_or_count]++);
	else
		sprintf(temp_buf, "%d-%d", rank, to_who_or_count);
		
	fprintf(translog, "transmission %s\n", temp_buf);
	
	fprintf(translog, "send count = %d\n", n);
	for (int j = 0; j < n; ++j)	{
		fprintf(translog, "%s\n", str_out);
		str_out += strlen(str_out) + 1;
	}
	
	fprintf(translog, "recv count = %d\n", m);
	for (int j = 0; j < m; ++j)	{
		fprintf(translog, "%s\n", str_in);
		str_in += strlen(str_in) + 1;
	}
	
	fprintf(translog, "indivs sent =");
	for (int j = 0; j < n; ++j)	{
		fprintf(translog, " %d", to_send[j]);
	}
	fprintf(translog, "\n");
	
	fprintf(translog, "indivs replaced =");
	for (int j = 0; j < m; ++j)	{
		fprintf(translog, " %d", to_replace[j]);
	}
	fprintf(translog, "\n");
	
	fprintf(translog, "old scores\n");
	for (int j = 0; j < nindivs; ++j)
		fprintf(translog, "%f\n", old_scores[j]);
		
	fprintf(translog, "new scores\n");
	for (int j = 0; j < nindivs; ++j)
		fprintf(translog, "%f\n", new_scores[j]);
	
	fprintf(translog, "\n");
	fflush(translog);
	fclose(translog);
	return 0;
}

int RecvMPIMessage(char** buf_, int* size_, int* who_, int* tag_, bool block)	{
	int& who = *who_;
	int& tag = *tag_;
	int& size = *size_;
	char*& buf = *buf_;
	
	int flag;
	MPI_Status status;
	
	MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
	
	if (flag == 0)	{
		if (block == false)
			return 0;
		else
			MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	}

	MPI_Get_count(&status, MPI_CHAR, &size);
	buf = new char[size];
	MPI_Recv(buf, size, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
	who = status.MPI_SOURCE;
	tag = status.MPI_TAG;
	
	return 1;
}

int RecvMPIMessage(char** buf_, int* size_, int who, int* tag_, bool block)	{
	int& tag = *tag_;
	int& size = *size_;
	char*& buf = *buf_;
	
	int flag;
	MPI_Status status;
	
	MPI_Iprobe(who, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
	
	if (flag == 0)	{
		if (block == false)
			return 0;
		else
			MPI_Probe(who, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	}

	MPI_Get_count(&status, MPI_CHAR, &size);
	buf = new char[size];
	MPI_Recv(buf, size, MPI_CHAR, who, status.MPI_TAG, MPI_COMM_WORLD, &status);
	tag = status.MPI_TAG;
	
	return 1;
}

int RecvMPIMessage(char** buf_, int* size_, int* who_, int tag, bool block)	{
	int& who = *who_;
	int& size = *size_;
	char*& buf = *buf_;
	
	int flag;
	MPI_Status status;
	
	MPI_Iprobe(MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &flag, &status);
	
	if (flag == 0)	{
		if (block == false)
			return 0;
		else
			MPI_Probe(MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
	}

	MPI_Get_count(&status, MPI_CHAR, &size);
	buf = new char[size];
	MPI_Recv(buf, size, MPI_CHAR, status.MPI_SOURCE, tag, MPI_COMM_WORLD, &status);
	who = status.MPI_SOURCE;
	
	return 1;
}

int RecvMPIMessage(char** buf_, int* size_, int who, int tag, bool block)	{
	int& size = *size_;
	char*& buf = *buf_;
	
	int flag;
	MPI_Status status;
	
	MPI_Iprobe(who, tag, MPI_COMM_WORLD, &flag, &status);
	
	if (flag == 0)	{
		if (block == false)
			return 0;
		else
			MPI_Probe(who, tag, MPI_COMM_WORLD, &status);
	}

	MPI_Get_count(&status, MPI_CHAR, &size);
	buf = new char[size];
	MPI_Recv(buf, size, MPI_CHAR, who, tag, MPI_COMM_WORLD, &status);
	
	return 1;
}

int SendMPIMessage(char* buf, int size, int who, int tag)	{
	return MPI_Send(buf, size, MPI_CHAR, who, tag, MPI_COMM_WORLD);
}

int LogConfig(const GeneralGamlConfig& conf)	{
	debug_mpi("logging GamlConfig structure...");
	
	debug_mpi("[general]");
	debug_mpi("logevery = %d", conf.logevery);
	debug_mpi("saveevery = %d", conf.saveevery);
	debug_mpi("datafname = %s", conf.datafname.c_str());
	debug_mpi("streefname = %s", conf.streefname.c_str());
	debug_mpi("ofprefix = %s", conf.ofprefix.c_str());
	
	debug_mpi("[master]");
	debug_mpi("holdover = %d", conf.holdover);
	debug_mpi("nindivs = %d %d", conf.min_nindivs, conf.max_nindivs);
	debug_mpi("stopgen = %d", conf.stopgen);
	
	debug_mpi("[remote]");
	debug_mpi("holdover = %d", conf.holdover);
	debug_mpi("nindivs = %d %d", conf.min_nindivs, conf.max_nindivs);
	debug_mpi("stopgen = %d", conf.stopgen);
	
	debug_mpi("done logging GamlConfig structure");
	return 0;
}

int LogParams(const Parameters& params)	{
	debug_mpi("logging Parameters (partial) structure...");
	debug_mpi("holdover = %d", params.holdover);
	debug_mpi("nindivs = %d", params.nindivs);
	debug_mpi("randomSeed = %d", params.randomSeed);
	return 0;
}

int LogTreeStrings(const char* tree_strings)	{
	int count = 0;
	const char* p = tree_strings;
	while (*p)	{
		p += strlen(p) + 1;
		++count;
	}
	debug_mpi("%d tree strings:", count);
	p = tree_strings;
	while (*p)	{
		debug_mpi("%s", p);
		p += strlen(p) + 1;
	}
	return 0;
}

int LogKappas(const double* kappa_probs, const int count)	{
	debug_mpi("%d kappa probs:", count);
	for (int i = 0; i < count; ++i)
		debug_mpi("%f", kappa_probs[i]);
	return 0;
}

int LogPis(const double* pis, const int count)	{
	debug_mpi("%d pis:", count);
	for (int i = 0; i < count; ++i)
		for(int b=0;b<4;b++)
			debug_mpi("%f, ", pis[i*4+b]);
	return 0;
}

int CountTreeStrings(char* p)	{
	int count = 0;
	while (*p)	{
		p += strlen(p) + 1;
		++count;
	}
	return count;
}

// string length for a string that is terminated by a double null
int strlen2(char* p)	{
	int count = 0;
	while (*p || *(p+1))	{
		++count;
		++p;
	}
	return count;
}

int CalcMaxIndivs(const HKYData& data, int mem)	{
	const int KB = 1024;
	const int MB = KB*KB;
	int sizeof_treenode = 4*data.NChar()*sizeof(double);
	debug_mpi("sizeof_treenode = %d", sizeof_treenode);
	int num_internal_treenode_per_indiv = data.NTax()-2;
	int size_of_terminal_indivs=sizeof(int) * data.NChar();
	debug_mpi("num_internal_treenode_per_indiv = %d", num_internal_treenode_per_indiv);
	int sizeof_indiv = sizeof_treenode*num_internal_treenode_per_indiv*2;
	sizeof_indiv+=size_of_terminal_indivs*data.NTax();
	debug_mpi("sizeof_indiv = %d", sizeof_indiv);
	debug_mpi("retval = %d", mem*MB / sizeof_indiv);
	return mem*MB / sizeof_indiv;
}

int RemoteShieldedMigrants(Population& pop, const GeneralGamlConfig& conf)	{
/*	int size, rank, nprocs, count = conf.numshields, restart_count = 0;
	int *which = new int[count];
	char *tree_strings, fname[64];
	double *models, old_score;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	if (rank < 10)	sprintf(fname, "log0%d.log", rank);
	else	sprintf(fname, "log%d.log", rank);
	
	pop.CalcAverageFitness();
	old_score = pop.bestFitness;
	
	//start the tree log file
	pop.CreateTreeLog(rank);
	
	for (pop.gen = 1; pop.gen < conf.stopgen; ++pop.gen)	{
		if (pop.gen % conf.repeatthresh == 0)	{
			if (pop.bestFitness - old_score <= conf.scorethresh)	{
				debug_mpi("repeat thresh=%d", conf.repeatthresh); 
				debug_mpi("score threshold exceeded, restarting population (%f, %f)", old_score, pop.bestFitness);
				pop.Restart(1, rank, nprocs, restart_count++);	
				pop.CalcAverageFitness();
			}
			old_score = pop.bestFitness;
		}
		if (pop.gen % conf.logevery == 0)
			pop.Log(fname, 0.0);
		if (pop.gen % conf.interval == 0)	{

			//Check for a message to change remote type
			int tag;
			RecvMPIMessage(&tree_strings, &size, 0, &tag, false);
			if(tag == TAG_REMOTE_TYPE_SWITCH){
				//call RemoteAlphaMaleReplication
				debug_mpi("\tchanging remote type to AMR...\n");
				RemoteAlphaMaleReplication(pop, conf);
				return 0;//if we return from RemoteAMR, everything must be done, so return
				}
		
			debug_mpi("STARTING SYNCHRONOUS COMMUNICATION (node %d)", 0);
			pop.GetNBestIndivIndices(&which, count);
			//pop.GetNRandomIndivIndices(&which, conf.numshieldedpernode);  // alternatively
			pop.GetSpecifiedTreeStrings(&tree_strings, count, which);
			SendMPIMessage(tree_strings, strlen2(tree_strings)+2, 0, TAG_TREE_STRINGS);
			
			int model_size=pop.GetSpecifiedModels(&models, count, which);
			SendMPIMessage((char*) models, count*sizeof(double)*model_size, 0, TAG_MODEL);
						
			debug_mpi("\tsent: %d tree strings", count);
			debug_mpi("\tsent: %d models, size: %d", count, model_size);
						
			delete [] tree_strings;
			delete [] models;
		
		}
		//adding this to make log files of trees for each population
//		if(pop.gen % pop.params->saveEvery == 0) pop.AppendTreeToTreeLog( rank );
		pop.NextGeneration();
		pop.OutputFate();
	}
	SendMPIMessage(NULL, 0, 0, TAG_QUIT);
	debug_mpi("STARTING SYNCHRONOUS COMMUNICATION (node %d)", 0);
	debug_mpi("\tsent: quit message");
	delete [] which;
	// TODO what to do on quit?
	debug_mpi("quitting");
*/	return 0;
}
/*
int RemoteAlphaMaleReplication(Population& pop, const GeneralGamlConfig& conf)	{
	int which, *all, size, rank, tag;
	char *tree_strings, *buf;
	double score, *models;
	char fname[32];
	
	all = new int[pop.params->nindivs];
	for (int i = 0; i < pop.params->nindivs; ++i)
		all[i] = i;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank < 10) sprintf(fname, "log0%d.log", rank);
	else sprintf(fname, "log%d.log", rank);
		
	pop.CalcAverageFitness();

	for (pop.gen = 1; pop.gen < conf.stopgen; ++pop.gen)	{
		if (pop.gen % conf.logevery == 0)
			pop.OutputLog();
		if (pop.gen % conf.interval == 0)	{
			debug_mpi("STARTING SYNCHRONOUS COMMUNICATION (node 0)");
			score = pop.bestFitness;
			SendMPIMessage((char*)&score, sizeof(double), 0, TAG_SCORE);
			debug_mpi("\tsend: score = %f", score);
			RecvMPIMessage(&tree_strings, &size, 0, &tag);
			if (tag == TAG_TREE_STRINGS_REQUEST)	{
				debug_mpi("\trecv: tree strings request");
				which = (int)pop.cumfit[pop.total_size-1][0];
				pop.GetSpecifiedTreeStrings(&tree_strings, 1, &which);
				SendMPIMessage(tree_strings, strlen2(tree_strings)+2, 0, TAG_TREE_STRINGS);
				
				int model_size=pop.GetSpecifiedModels(&models, 1, &which);
				SendMPIMessage((char*)models, sizeof(double)*model_size, 0, TAG_MODEL);
				
				debug_mpi("\tsend: %d tree string", 1);
				debug_mpi("\tsend: %d models, size:%d", 1, model_size*sizeof(double));	
	
				delete [] tree_strings;
				delete [] models;
			}
			else if (tag == TAG_TREE_STRINGS)	{
				debug_mpi("\trecv: %d tree strings", CountTreeStrings(tree_strings));
				RecvMPIMessage(&buf, &size, 0, TAG_MODEL);
				models=(double*) buf;

				pop.ReplicateSpecifiedIndividuals(pop.total_size, all, tree_strings, models);
				pop.CalcAverageFitness();
				delete [] tree_strings;
				delete [] models;
			}
			else	{
				debug_mpi("alpha male replication recved bad message tag");
				assert(false);
			}
		}
		//adding this to make log files of trees for each population
//		if(pop.gen % pop.params->saveEvery == 0) pop.AppendTreeToTreeLog( rank );
	//	if(pop.gen % 100) pop.NNIoptimization();
		pop.NextGeneration();
		pop.OutputFate();
	}
	debug_mpi("sending quit message");
	SendMPIMessage(NULL, 0, 0, TAG_QUIT);
	delete [] all;
	return 0;
}
*/

void RemoteSendBestTree(Population& pop){
	int *which=new int;
	char *tree_strings;
	double *models;
	
	pop.GetNBestIndivIndices(&which, 1);
	pop.GetSpecifiedTreeStrings(&tree_strings, 1, which);
	int size=strlen2(tree_strings)+2;
//	debug_mpi("about to send treestrings...");
	SendMPIMessage(tree_strings, size, 0, TAG_TREE_STRINGS);
	debug_mpi("\tsent ind %d, lnL %f", *which, pop.indiv[*which].Fitness());
//	debug_mpi("about to send modelstrings...");
	int model_size=pop.GetSpecifiedModels(&models, 1, which);
	SendMPIMessage((char*) models, sizeof(double)*model_size, 0, TAG_MODEL);
//	debug_mpi("about to send subdef...");
	char std[5];
	sprintf(std, "%d", pop.subtreeDefNumber);
	SendMPIMessage(std, strlen(std)+2, 0, TAG_SUBTREE_ITERATION);
	
	if(pop.subtreeDefNumber!=0){
		char stn[10];
		sprintf(stn, "%d", pop.subtreeNode);
		SendMPIMessage(stn, strlen(stn)+2, 0, TAG_SUBTREE_DEFINE);
		debug_mpi("\tvalid for subtree def %d, node %d", pop.subtreeDefNumber, pop.subtreeNode);
		}
	else
		debug_mpi("\tno defined subtree");
		
	//finally, send the score
	double score = pop.indiv[*which].Fitness();
	SendMPIMessage((char*)&score, sizeof(double), 0, TAG_SCORE);	
	
	delete which;
	delete [] tree_strings;
	delete [] models;
	}
	


int RemoteSubtreeWorker(Population& pop, const GeneralGamlConfig& conf){
	int *which, size, rank, tag;
	char *tree_strings, *buf;
	double score, *models;
	bool perturb;
	
	which=new int[5];
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		
	pop.CalcAverageFitness();

	int lastSend=g_sw->SplitTime();	

	cout << "Remote number " << rank << " running..." << endl;

	for (pop.gen = 1; pop.gen < conf.stopgen;){

		pop.keepTrack();
		pop.OutputFate();
		if (pop.gen % conf.logevery == 0)
			pop.OutputLog();
		++pop.gen;
		pop.NextGeneration();

	           if(pop.gen % pop.adap->intervalLength == 0){
	                bool reduced=false;
	                if(pop.gen-pop.lastTopoImprove >= pop.adap->intervalsToStore*pop.adap->intervalLength){
	                   reduced=pop.adap->ReducePrecision();
	                   }
	                if(reduced){
	                   pop.lastTopoImprove=pop.gen;
	                   pop.indiv[pop.bestIndiv].treeStruct->OptimizeAllBranches(pop.adap->branchOptPrecision);
	                   pop.indiv[pop.bestIndiv].SetDirty();
	                   pop.CalcAverageFitness();
	                        }
/*	                 else if(!(pop.gen%(pop.adap->intervalLength*pop.adap->intervalsToStore))){
	                        pop.indiv[pop.bestIndiv].treeStruct->OptimizeAllBranches(pop.adap->branchOptPrecision);
	                        pop.indiv[pop.bestIndiv].SetDirty();
	                        pop.CalcAverageFitness();
	                        }
*/					}
			


		if(g_sw->SplitTime() - lastSend > conf.sendInterval){
			debug_mpi("SYNCH COMM (node 0)");
			//send our best individual to the master
			RemoteSendBestTree(pop);		
			lastSend=g_sw->SplitTime();
			if(pop.params->stoptime - g_sw->SplitTime() < 0){
				debug_mpi("time limit of %d seconds reached...", pop.params->stoptime);
				break;
				}
			}
		//Check for a new tree from the master
		bool firstmessage=true;
		bool gotmessage=false;
		int subtreeNode;
		while(RecvMPIMessage(&tree_strings, &size, 0, &tag, false)==true){
			//check for a quit message
			if(tag == TAG_QUIT) {
				debug_mpi("\trecv: quit message");
				delete [] which;
				debug_mpi("quitting");
				return 0;	
				}
			//
			bool gotNewIndiv=false;
			int recievedDefNumber;
			debug_mpi("SYNCH COMM (node 0)");
			gotmessage=true;
			assert(tag == TAG_TREE_STRINGS || tag==TAG_PERTURB);
			if(firstmessage==false) debug_mpi("\tfound a newer message...");
			if(tag != TAG_PERTURB){
				gotNewIndiv=true;
				RecvMPIMessage(&buf, &size, 0, &tag, true);
				assert(tag == TAG_MODEL);
				models=(double*) buf;
			
				debug_mpi("\tgot new ind" );
				RecvMPIMessage(&buf, &size, 0, &tag, true);
		
	//		if(tag != TAG_PERTURB){
				perturb=false;
				assert(tag == TAG_SUBTREE_DEFINE);
				subtreeNode=atoi(buf);
				if(subtreeNode!=0){
					delete []buf;
					RecvMPIMessage(&buf, &size, 0, &tag, true);
					assert(tag == TAG_SUBTREE_ITERATION);
					recievedDefNumber=atoi(buf);
					debug_mpi("\tworking on subtree def %d, node %d", recievedDefNumber, subtreeNode);
					}
				else recievedDefNumber=0;
				}
			else{
				pop.pertMan->pertType=atoi(tree_strings);
				perturb=true;
				}	
			
			//if the current best and the new tree are either both accurate for the same subtree def or both
			//inaccurate for subtrees, just replace the worst individual, rather than the
			// whole pop, that way if the tree is old and worse that what the remote
			// already has it won't matter
			if(gotNewIndiv){
				*which=(int)pop.cumfit[0][0];
				debug_mpi("\treplacing indiv %d", *which);
				pop.ReplaceSpecifiedIndividuals(1, which, tree_strings, models);
				if(recievedDefNumber!=pop.subtreeDefNumber || (pop.subtreeNode!=0 && subtreeNode!=0)){
					pop.AssignSubtree(subtreeNode, *which);
					pop.CalcAverageFitness();
					debug_mpi("\tfilling pop with clones of %d", *which);
					pop.SetNewBestIndiv(*which);
					pop.FillPopWithClonesOfBest();
					pop.subtreeDefNumber=recievedDefNumber;
					}

				delete [] models;
				delete [] buf;
				}
#ifdef INCLUDE_PERTURBATION
			if(perturb==true){
				pop.CalcAverageFitness();
				if(pop.pertMan->pertType==1){
					debug_mpi("peforming NNI perturbation...");
					int toReplace=(pop.bestIndiv == 0 ? 1 : 0);
					pop.AppendTreeToTreeLog(-1, pop.bestIndiv);
					pop.NNIPerturbation(pop.bestIndiv, toReplace);
					pop.SetNewBestIndiv(toReplace);
					pop.FillPopWithClonesOfBest();
					pop.AppendTreeToTreeLog(-1, pop.bestIndiv);
					}
				else if(pop.pertMan->pertType==2){
					debug_mpi("peforming SPR perturbation...");
					int toReplace=(pop.bestIndiv == 0 ? 1 : 0);
					pop.AppendTreeToTreeLog(-1, pop.bestIndiv);
					pop.SPRPerturbation(pop.bestIndiv, toReplace);
					pop.SetNewBestIndiv(toReplace);
					pop.FillPopWithClonesOfBest();
					pop.AppendTreeToTreeLog(-1, pop.bestIndiv);
					}
				else assert(0);
				}				
#endif

			delete [] tree_strings;
			tag=0;
			firstmessage=false;
			}
		if(gotmessage==true){
//			if(pop.subtreeNode != subtreeNode) pop.AssignSubtree(subtreeNode);
			pop.CalcAverageFitness();
			debug_mpi("\tbest score= %f", pop.indiv[*which].Fitness());
			pop.AppendTreeToTreeLog(-1, *which);
			}
		}

	SendMPIMessage(NULL, 0, 0, TAG_QUIT);
	debug_mpi("\tsent: quit message");
	delete [] which;
	pop.FinalizeOutputStreams();
	debug_mpi("quitting");
	return 0;
	}



#endif // #ifdef MPI_VERSION

