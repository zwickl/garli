// GARLI version 0.96b8 source code
// Copyright 2005-2008 Derrick J. Zwickl
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

#ifdef SUBROUTINE_GARLI

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include "defs.h"
#include "string.h"
#include "mpi.h"
#include <time.h>
#include "funcs.h"

using namespace std;

void SubGarliMain(int);
void jobloop(int, int, MPI_Comm, int numJobs);
void UsageMessage(char *execName);

int main(int argc,char **argv){

  int rc = MPI_Init(&argc,&argv); 
  if(rc != MPI_SUCCESS){
	cout << "Error starting MPI.  Terminating." << endl;
	MPI_Abort(MPI_COMM_WORLD, rc);
	}

  MPI_Comm comm,mycomm; 
  int nproc, rank;
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm,&nproc); 
  MPI_Comm_rank(comm,&rank);

//  MPI_Comm_split(comm,rank,rank,&mycomm);
  int numJobsTotal = 0;
  
  if(rank == 0){
  	cout << "Garli started with command line:" << endl;
	for(int i=0;i<argc;i++) cout << argv[i] << " ";
	cout << endl;  

	if(argc == 1 || (argv[1][0] != '-' && !isdigit(argv[1][0]))){
		cout << "ERROR:\nGarli is expecting the number of jobs to be run to follow\nthe executable name on the command line" << endl;
		UsageMessage(argv[0]);
		exit(1);
		}
	else{
		if(argv[1][0] == '-') numJobsTotal = atoi(&argv[1][1]);
		else numJobsTotal = atoi(&argv[1][0]);
/*		for(int i=1;i<nproc;i++){
			rc = MPI_SSend(&numJobsTotal, 1, MPI_INT, i, 0, comm);
			if(rc != MPI_SUCCESS){
				cout << "Error communicating with process " << i << ".  Terminating." << endl;
				MPI_Abort(MPI_COMM_WORLD, rc);
				}
			}
*/		}

	}
 /*  else{
	MPI_Status stat;
	rc = MPI_Recv(&numJobsTotal, 1, MPI_INT, 0, MPI_ANY_TAG, comm, &stat);
	if(rc != MPI_SUCCESS){
		cout << "Error communicating with process " << 0 << ".  Terminating." << endl;
		MPI_Abort(MPI_COMM_WORLD, rc);
		}
	}
*/
  cout << "process " << rank << " about to broadcast" << endl;
  MPI_Bcast(&numJobsTotal, 1, MPI_INT, 0, comm);

  jobloop(rank,nproc,mycomm,numJobsTotal);

  cout << "process " << rank << " retunred from jobloop" << endl;

  MPI_Barrier(comm);
  cout << "process " << rank << " passed barrier" << endl;
  if(rank == 0){
	char temp[100];
	for(int i=0;i<numJobsTotal;i++){
		sprintf(temp, ".lock%d", i);
		remove(temp);
		}
	}
  
  MPI_Finalize();
  return 0;
}

void jobloop(int mytid,int ntids,MPI_Comm comm, int numJobs)
{
	//to start off with each process takes the run = to its tid
	int jobNum=mytid;	
	char temp[100];

	//this wait ensures that two runs don't start with the same seed
	timespec wait;
	wait.tv_sec = mytid * 2;
	wait.tv_nsec=0;
	nanosleep(&wait, NULL);

	sprintf(temp, ".lock%d", jobNum);
	ofstream lock(temp);
	lock.close();

	while(jobNum < numJobs){
		cout << "process " << mytid << " starting run #" << jobNum << endl;
		SubGarliMain(jobNum);
		jobNum++;
		while(jobNum < numJobs){
			sprintf(temp, ".lock%d", jobNum);
			if(FileExists(temp)) jobNum++;
			else{
				ofstream lock(temp);
				lock.close();
				break;
				}
			}
		}
	return;
}
#endif
