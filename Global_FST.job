#!/bin/bash
#PBS -l select=1:ncpus=1
#PBS -l walltime=5000:00:00
#PBS -N Global_FST
#PBS -j oe
#PBS -m ae
#PBS -M lmoreira@amnh.org
#PBS -k oe

# change to the working directory
cd $PBS_O_WORKDIR
echo "pbsworkdir"
echo $PBS_O_WORKDIR
EXECDIR=`pwd`
export PATH=./:$PATH
echo $PATH

date
time

for i in *.fst.idx; do echo "Get global estimate for" $i; realSFS fst stats $i; done > Global_pairwise.fst
