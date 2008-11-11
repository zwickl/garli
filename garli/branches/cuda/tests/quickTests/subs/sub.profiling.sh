#!/bin/bash

#S -S /bin/bash -cwd
#$ -o out.profiling.$JOB_ID -j y
#$ -M zwickl@duke.edu
cd ../profiling
./doRuns ~/bin/working.iccProfile -b


