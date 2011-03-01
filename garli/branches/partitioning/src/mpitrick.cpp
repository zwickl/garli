// GARLI version 1.00 source code
// Copyright 2005-2010 Derrick J. Zwickl
// email zwickl@nescent.org
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

#if defined(SUBROUTINE_GARLI) || defined(OLD_SUBROUTINE_GARLI)

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include "defs.h"
#include "string.h"
#include "mpi.h"
#include <time.h>
#include "funcs.h"
#include "outputman.h"

using namespace std;

int SubGarliMain(int);

void UsageMessage(char *execName);

extern OutputManager outman;

//old (parallel batch) and new (parallel replicates) mpi behavior now rolled into a single function

#if(1)
int jobloop(int, int, MPI_Comm, int numJobs);

string MyFormattedTime(){
	time_t rawtime;
	struct tm * timeinfo;
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	string s = asctime(timeinfo);
	s.erase(s.end()-1, s.end());
	return s;
	}

int main(int argc,char **argv){

  if(argc == 2){
	if(!strcmp(argv[1], "--help") || !strcmp(argv[1], "--h") || !strcmp(argv[1], "-help") || !strcmp(argv[1], "-h")){
		UsageMessage(argv[0]);
		return 0;
		}
	}

  int rc = MPI_Init(&argc,&argv); 
  if(rc != MPI_SUCCESS){
	outman.SetLogFile("mpi_messages.log");
    	outman.UserMessage("Error starting MPI.  Terminating.");
	MPI_Abort(MPI_COMM_WORLD, rc);
	}

  MPI_Comm comm,mycomm; 
  int nproc, rank;
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm,&nproc);
  MPI_Comm_rank(comm,&rank);

  timespec wait;
  int numJobsTotal = 0;
  
  if(rank == 0){
    	outman.SetLogFile("mpi_messages.log");
	outman.UserMessageNoCR("MPI Garli started with command line: ");
	for(int i=0;i<argc;i++) outman.UserMessageNoCR("%s ", argv[i]);
        outman.UserMessage("\n");
#ifdef OLD_SUBROUTINE_GARLI
        outman.UserMessage("This is the original batch MPI GARLI version.  It expects a series of configuration");
        outman.UserMessage("files named \"run0.conf\", \"run1.conf\", etc.  If no number is passed on the command");
        outman.UserMessage("line after the executable name, then it assumes that the # configs = # processors.");
        outman.UserMessage("Otherwise it looks for the specified number of configs.\n");

	if(argc > 1){
		if(! isdigit(argv[1][0])){
			outman.UserMessage("***ERROR***:GARLI is expecting <exe> <total # configs>\n\tor\n\t<exe> <nothing>\n\tGot <exe> %s", argv[1]);
			UsageMessage(argv[0]);
			MPI_Finalize();
			return 1;
			}
		else numJobsTotal = atoi(argv[1]);
		}
	else numJobsTotal = nproc;
#else
	if(argc == 1 || (argv[1][0] != '-' && !isdigit(argv[1][0]))){
		outman.UserMessage("***ERROR***:Garli is expecting the number of jobs to be run to follow\n\tthe executable name on the command line\n");
		UsageMessage(argv[0]);
		MPI_Finalize();
		return 1;
		}
	else{
		if(argv[1][0] == '-') numJobsTotal = atoi(&argv[1][1]);
		else numJobsTotal = atoi(&argv[1][0]);
		}
#endif
  	outman.UserMessage("#####%d total executions of the config file were requested######", numJobsTotal);
	}
  else{//wait a moment for proc 0 to output to the messages file, then attach to the stream
	wait.tv_sec = 1;
	wait.tv_nsec=0;
	nanosleep(&wait, NULL);
	outman.SetLogFileForAppend("mpi_messages.log");
	}
 
 
 //These barriers really shouldn't be necessary, but adding them seemed to resolve a weird issue that Jelesko
 //was having where processes with rank >= 8 were hanging in the Bcast until some of the first 8 were completely 
 //finished and returned from the job loop
  //outman.UserMessage("#####Process %d approaching barrier 1 at %s######", rank, MyFormattedTime().c_str());
  MPI_Barrier(comm);
  //outman.UserMessage("#####Process %d passed barrier 1 at %s######", rank, MyFormattedTime().c_str());
 
//send all of the processors the number of jobs total
  MPI_Bcast(&numJobsTotal, 1, MPI_INT, 0, comm);

//  outman.UserMessage("#####Process %d passed broadcast, approaching barrier 2 at %s######", rank, MyFormattedTime().c_str());
  MPI_Barrier(comm);
//  outman.UserMessage("#####Process %d passed barrier 2 at %s######", rank, MyFormattedTime().c_str());

//DEBUG
  //if startjob lockfiles exist at this point that must mean that a previous run bailed.  Remove them.
  //Then processes will start those same runs and possibly restart from checkpoint (if any were being 
  //written and restart=1 was specified).  Otherwise they will just start them over.
  if(rank == 0){
	for(int j = 0;j < numJobsTotal;j++){
		char startfile[100], donefile[100];
		sprintf(startfile, ".s-lock%d", j);
		sprintf(donefile, ".d-lock%d", j);
		if(FileExists(donefile)){
			outman.UserMessage("It appears that run %d was completed in a previous MPI invocation.\n\tRun %d will not be re-run unless the hidden file \"%s\" and any checkpoint files for this run (if present) are deleted from this directory.", j, j, donefile);
			}
		else if(FileExists(startfile)){
			outman.UserMessage("It appears that run %d was started but not completed in a previous MPI invocation.\n\tRun %d will either be re-run or restarted from a checkpoint (if restart = 1 was specified in the GARLI config file).", j, j);
			remove(startfile);
			}
  		}
	}

  int jobsCompleted = jobloop(rank,nproc,mycomm,numJobsTotal);
  outman.SetLogFileForAppend("mpi_messages.log");
  if(jobsCompleted > -1){
#ifdef OLD_SUBROUTINE_GARLI
 	outman.UserMessage("process %d finished, did %d run(s), no more configs to execute at %s. Waiting for other procs...", rank, jobsCompleted,  MyFormattedTime().c_str());
#else
 	outman.UserMessage("process %d finished, did %d run(s), no further runs to do at %s. Waiting for other procs...", rank, jobsCompleted,  MyFormattedTime().c_str());
#endif
	}

  MPI_Barrier(comm);
  if(rank == 0)  outman.UserMessage("all processes completed at %s", MyFormattedTime().c_str());
  else nanosleep(&wait, NULL);//this is just to keep proper ordering in the output file
 
  outman.UserMessage("process %d terminating", rank);
  
//Not sure if deleting lock files should or should not be done.
/*  if(rank == 0){
	char temp[100];
	for(int i=0;i<numJobsTotal;i++){
		sprintf(temp, ".lock%d", i);
		remove(temp);
		}
	}
  */
  MPI_Finalize();
  return 0;
}

int jobloop(int mytid,int ntids,MPI_Comm comm, int numJobs){

	//to start off with each process takes the run = to its tid
	int jobNum=mytid;	
	char temp[100];
	int jobsCompleted = 0;

#ifdef SUBROUTINE_GARLI
	//this wait ensures that two runs don't start with the same seed
	timespec wait;
	wait.tv_sec = mytid * 2;
	wait.tv_nsec=0;
	nanosleep(&wait, NULL);
#endif
	while(jobNum < numJobs){
		//DEBUG
		//sprintf(temp, ".lock%d", jobNum);
		ofstream lock;
		char startfile[100], donefile[100];
		sprintf(startfile, ".s-lock%d", jobNum);
		sprintf(donefile, ".d-lock%d", jobNum);
		//DEBUG
		
		if(FileExists(donefile) || FileExists(startfile)) jobNum++;
		else{
			lock.open(startfile);
			lock.close();
			outman.SetLogFileForAppend("mpi_messages.log");
			outman.UserMessage("process %d starting run %d at %s", mytid, jobNum, MyFormattedTime().c_str());
			int err = SubGarliMain(jobNum);
			if(err){
				outman.SetLogFileForAppend("mpi_messages.log");
				outman.UserMessage("***process %d aborted run %d at %s", mytid, jobNum, MyFormattedTime().c_str());
				outman.UserMessage("\tsee the <filename>.screen.log files for details on what went wrong");
				return -1;
				}
			jobsCompleted++;
			lock.open(donefile);
			lock.close();
			remove(startfile);
			jobNum++;
			}
		}
	return jobsCompleted;
}

#elif defined(__OLD_SUBROUTINE_GARLI)
void jobloop(int, int, MPI_Comm, int numJobs=-1);

int main(int argc,char **argv) {

	MPI_Comm comm,mycomm; int ntids,mytid;
	MPI_Init(&argc,&argv); comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm,&ntids); MPI_Comm_rank(comm,&mytid);
	MPI_Comm_split(comm,mytid,mytid,&mycomm);
	
	if(mytid == 0){
		outman.SetLogFile("mpi_messages.log");
		outman.UserMessageNoCR("MPI Garli started with command line: ");
		for(int i=0;i<argc;i++) outman.UserMessageNoCR("%s ", argv[i]);
		outman.UserMessage("\n");
		outman.UserMessage("This is the original batch MPI GARLI version.  It expects a series of configuration");
		outman.UserMessage("files named \"run0.conf\", \"run1.conf\", etc.  If no number is passed on the command");
		outman.UserMessage("line after the executable name, then it assumes that the # configs = # processors.");
		outman.UserMessage("Otherwise it looks for the specified number of configs.");
		}
	else{
	 	timespec wait;
		wait.tv_sec = 2;
		wait.tv_nsec=0;
		nanosleep(&wait, NULL);
		outman.SetLogFileForAppend("mpi_messages.log");
		}

	if(argc > 1){
		if(! isdigit(argv[1][0])){
			cout << "I'm confused.\nExpecting <exe> <total # jobs>\n  or\n<exe> <nothing>\nGot <exe> " << argv[1] << endl;
			MPI_Finalize();
			return 1;
			}
		jobloop(mytid,ntids,mycomm, atoi(argv[1]));
		}
	else jobloop(mytid, ntids, mycomm);
	
	MPI_Finalize();
	return 0;
	}

void jobloop(int mytid,int ntids,MPI_Comm comm, int numJobs /*=-1*/){

	if(numJobs < 0) numJobs=ntids;
	int jobNum=mytid;	
	cout << "total proc: " << ntids << " total jobs: " << numJobs << endl;
	while(jobNum < numJobs){
		timespec wait;
		wait.tv_sec = mytid * 2;
		wait.tv_nsec=0;
		nanosleep(&wait, NULL);
		SubGarliMain(jobNum);
		jobNum += ntids;
		}
	return;
	}
#endif
#endif
