#!/bin/bash

#S -S /bin/bash -cwd
#$ -o out.memcheck.$JOB_ID -j y
#$ -M zwickl@duke.edu
cd ../memcheck
./doRuns ~/bin/working.iccMemcheck -b


