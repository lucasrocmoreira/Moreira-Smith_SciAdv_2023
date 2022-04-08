#!/bin/bash
#PBS -l select=1:ncpus=16
#PBS -l walltime=5000:00:00
#PBS -N ANGSD_SFS
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

module load parallel-20171122

ref=/home/lmoreira/reference/Picoides_pubescens_ref-pseudochromosome.v2/pseudochromosomes-subset_sorted.fasta

# for FST, saf should not be folded. the folding is done when producing the sfs 

date
time

ls pop-*.bamlist | parallel angsd -GL 2 \
-dosaf 1 \
-anc $ref \
-C 50 \
-ref $ref \
-baq 1 \
-minInd 6 \
-minMapQ 20 \
-minQ 20 \
-nThreads 16 \
-bam {} \
-out $PBS_O_WORKDIR/SAF/{}

