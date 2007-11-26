#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include "string.h"
#include "mpi.h"
#include <time.h>
#include "funcs.h"

using namespace std;

void SubGarliMain(int);
void jobloop(int, int, MPI_Comm, int numJobs);

int main(int argc,char **argv){

  cout << "Garli started with command line:" << endl;
  for(int i=0;i<argc;i++) cout << argv[i] << " ";
  cout << endl;  

  if(argc == 1 || !isdigit(argv[1][0])){
	cout << "ERROR:\nGarli is expecting the number of jobs to be run to follow\nthe executable name on the command line" << endl;
	return 1;
	}

  MPI_Comm comm,mycomm; int ntids,mytid;
  MPI_Init(&argc,&argv); comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm,&ntids); MPI_Comm_rank(comm,&mytid);
  MPI_Comm_split(comm,mytid,mytid,&mycomm);

  int numJobs = atoi(argv[1]);

  jobloop(mytid,ntids,mycomm,numJobs);
	
  MPI_Finalize();

  char temp[100];
  for(int i=0;i<numJobs;i++){
	sprintf(temp, ".lock%d", i);
	remove(temp);
	}

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
		do{
			jobNum++;
			sprintf(temp, ".lock%d", jobNum);
			}while(FileExists(temp) && (jobNum < numJobs));
		ofstream lock(temp);
		lock.close();
		}
	return;
}
