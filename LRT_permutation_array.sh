#!/bin/bash
#PBS -l select=1:ncpus=64:mem=32Gb
#PBS -N LRT_randomization
#PBS -J 1-200:1

# change to the working directory
cd $PBS_O_WORKDIR
echo "pbsworkdir"
echo $PBS_O_WORKDIR
EXECDIR=`pwd`
export PATH=./:$PATH
echo $PATH

date
time

# Run separately:
#for i in {1..200}; do mkdir -p $i; cat PC.env | shuf > $i/PC_shuf.env; done

angsd -beagle PP.beagleOut.complete_filtered.gprobs.gz \
-fai pseudochromosomes-subset_sorted.fasta.fai \
-P 5 -doMaf 4 -doAsso 4 -yQuant ${PBS_ARRAY_INDEX}/PC_shuf.env -cov PCA.cov -out ${PBS_ARRAY_INDEX}/PC.res

exit

# After running:
#for i in {1..200}; do echo $i; cd $i; zcat PC.res.lrt0.gz | awk '{print $7}' > LRT.lrt0; cd ..; done
#paste {1..200}/LRT.lrt0 > permutation_LRT_lrt0.tsv

# Then in R:
#library('data.table')
#permutations <- fread("permutation_LRT_lrt0.tsv")
#threshold_0.00001 <- apply(permutations, 1, quantile, probs=c(0.99999))
#write(threshold_0.00001,file = "LRT_threshold_0.00001.tsv",sep="\n")
