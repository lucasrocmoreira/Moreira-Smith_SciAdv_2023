#!/bin/bash
#PBS -l select=1:ncpus=16
#PBS -N ANGSD_asso

# change to the working directory
cd $PBS_O_WORKDIR
echo "pbsworkdir"
echo $PBS_O_WORKDIR
EXECDIR=`pwd`
export PATH=./:$PATH
echo $PATH

date
time

# Calculate genotype likelihoods and subselect SNPs
angsd/angsd -GL 2 -bam Picoides_pubescens.bamlist -doGeno 3 -doGlf 2 -doCounts 1 -dumpCounts 2 \
-doMaf 1 -doMajorMinor 1 -doPost 1 -minInd 50 -minMaf 0.05 -minMapQ 30 -minQ 20 -rf regions.txt \
-doPlink 2 -nThreads 16 -SNP_pval 0.000001 -out Picoides-pubescens

# Split gl file
zcat Picoides-pubescens.beagle.gz | split -l 1000000 - Picoides-pubescens.beagle.part

gzip -f -c Picoides-pubescens.beagle.partaa > Picoides-pubescens.beagle.partaa.ready.gz

for i in Picoides-pubescens.beagle.parta[a-z]; do cat beagle.header $i > $i.ready; done
gzip -f *.ready

# Impute genotypes
for i in Picoides-pubescens.*.ready.gz;
do 
    java -Xmx16g -jar ./beagle.jar like=$i out=PP.beagleOut
done

# Merge results for all files
first=1
for f in PP.beagleOut.Picoides-pubescens*.gprobs.gz
do
    if [ "$first" ]
    then
        zcat "$f"
        first=
    else
        zcat "$f"| tail -n +2
    fi
done > PP.beagleOut.complete.gprobs
gzip -f PP.beagleOut.complete.gprobs

# Run association
# Make sure the .fai file has no chromosomes with '_'
angsd -beagle PP.beagleOut.complete.gprobs.gz \
-fai pseudochromosomes-subset_sorted.fasta.fai \
-P 5 -doMaf 4 -doAsso 4 -yQuant PC.env -cov PCA.cov -out PC.res -nThreads 8

# Run association with no structure covariants
angsd -beagle PP.beagleOut.complete.gprobs.gz \
-fai pseudochromosomes-subset_sorted.fasta.fai \
-P 5 -doMaf 4 -doAsso 4 -yQuant PC.env -out PC.res.nostruct -nThreads 8

