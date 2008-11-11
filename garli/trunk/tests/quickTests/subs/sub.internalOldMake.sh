#!/bin/bash

#S -S /bin/bash -cwd
#$ -o out.internal.$JOB_ID -j y
#$ -M zwickl@duke.edu
cd ../internalOldMake
./doRuns ~/bin/working.oldMake -b


