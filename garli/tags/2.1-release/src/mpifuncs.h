// GARLI version 0.93 source code
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

#ifdef MPI_VERSION

#ifndef MPIFUNCS_H
#define MPIFUNCS_H

#include "configoptions.h"
#include "sequencedata.h"
#include "parameters.h"
#include "population.h"
#include "threaddcls.h"

int MPIMain(int arc, char** argv);

int StartProcs(const GeneralGamlConfig&, NucleotideData&);

int MasterMaster(MasterGamlConfig&, NucleotideData&);
int RemoteMaster(GeneralGamlConfig&, NucleotideData&);
int MasterFullDuplexExchange(Population& pop, const MasterGamlConfig& conf);
int RemoteFullDuplexExchange(Population& pop, const GeneralGamlConfig& conf);
int MasterShieldedMigrants(Population& pop, const MasterGamlConfig& conf);
int RemoteShieldedMigrants(Population& pop, const GeneralGamlConfig& conf);
int RemoteAlphaMaleReplication(Population& pop, const GeneralGamlConfig& conf);


int RemoteSubtreeWorker(Population& pop, const GeneralGamlConfig& conf);
int MasterAlphaMaleReplication(Population& pop, const MasterGamlConfig& conf);
int MasterHybrid(Population& pop, const MasterGamlConfig& conf);	
int MasterLastCall(Population& pop, int master_mem);

//int DoMasterSM(Population& pop, const GamlConfig& conf, int who, transferred_data_t results);
//int DoMasterAMR(Population& pop, const GamlConfig& conf, int who, transferred_data_t results);

int CalcMaxIndivs(const NucleotideData&, int);

// buf, size in bytes, from, tag, blocking
int RecvMPIMessage(char**, int*, int*, int*, bool block = true);
int RecvMPIMessage(char**, int*, int, int*, bool block = true);
int RecvMPIMessage(char**, int*, int*, int, bool block = true);
int RecvMPIMessage(char**, int*, int, int, bool block = true);
// buf, size, to, tag
int SendMPIMessage(char*, int, int, int);

int PollForResults(int nodes[]);
bool PollForResults(int n);

int GetResultsFromNode(int node_num, int* nindivs_, char** tree_strings_);
int SendResultsToNode(int node, int n, char* tree_strings);

int ReceiveParams(Parameters* params_, int node);
int ReceiveData(NucleotideData* data_, int node);

int debug_mpi(const char* fmt, ...);
int LogConfig(const GeneralGamlConfig&);
int LogParams(const Parameters& params);
int LogTreeStrings(const char* tree_strings);
int LogKappas(const double* kappa_probs, const int count);
int LogPis(const double* pis, const int count);
int TransLog(int count, int nindivs, int n, char* str_out, int m, char* str_in, int* to_send, int* to_replace, double* old_scores, double* new_scores);

int strlen2(char* p);
int CountTreeStrings(char* p);

#endif

#endif
