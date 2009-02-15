#!/bin/bash

#S -S /bin/bash -cwd
#$ -o out.internal.$JOB_ID -j y
#$ -M zwickl@duke.edu
cd ../internal
./doRuns ~/bin/working.norm -b


