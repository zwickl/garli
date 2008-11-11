#!/bin/bash

#S -S /bin/bash -cwd
#$ -o out.internalConst.$JOB_ID -j y
#$ -M zwickl@duke.edu
cd ../internalConstraints
./doRuns ~/bin/working.norm -b


