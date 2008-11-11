#!/bin/bash

#S -S /bin/bash -cwd
#$ -o out.internalOMP.$JOB_ID -j y
#$ -M zwickl@duke.edu
#$ -pe threaded 2
#$ -v OMP_NUM_THREADS=2
cd ../internalOMP
./doRuns ~/bin/working.OMP -b


