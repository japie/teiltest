#!/bin/bash
# Usage:
#  This script takes no argument and runs only on variable defined here:
# The total number of runs
#runs=200
# The number of runs to be executed in a job, it shall divide $runs
#group_size=2
#prog=../ecga
datadir=./
#  It goes through these steps:
#    1. It executes the program "../ecga" and takes arguments from
#       "inputfile" in background, multi-process simultaneously.
#    2. The output files is set to "raw/<num>.out", where <num> is the
#       number of that run. The standard output of ecga is piped to
#       /dev/null.
#    3. Compress and backup raw outputs to raw.tar.bz2
#    4. Then all outputs in "raw" are clarified and averaged to
#       standard out.

###### step 1. and 2. ######

cd $datadir

message="echo Working directory is \$PBS_O_WORKDIR;\
cd \$PBS_O_WORKDIR;\
echo Running on host \`hostname\`;\
echo Time is \`date\`;\
echo Directory is \`pwd\`;"
#jobname="Online_Modified_DSMGA-$$"

# generate input files
#for((i=0; i<$runs; i++)); do
#    rnd=$(($RANDOM%10000))
#    rnd_frac=$((rnd/1000))$((rnd%1000/100))$((rnd%100/10))$((rnd%10))
#    sed -e 's/\(seed *\)[0-9\.]*/\10.'$rnd_frac'/' ../inputfile > $i.in
#done
# submit jobs

idlower=1000
idupper=2000

for((i = ${idlower}; i <= ${idupper}; i += 1)); do
  qdel $i
done

# wait for all task to stop
#while [ "$(qstat | grep $jobname)" ]; do sleep 3; done
# remove input files and outputs from qsub
#sleep 5
#rm *.in $jobname.*

#cd ..

###### step 3. ######
#tar cjf $datadir.tar.bz2 $datadir

###### step 4. ######
#./comp_bit.sh $runs $datadir
