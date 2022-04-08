#!/bin/bash
#PBS -l select=1:ncpus=1
#PBS -l walltime=5000:00:00
#PBS -N ANGSD_FST
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

for i in PP-*.sfs;
do

	name=`echo $i | cut -d '.' -f1`
	pop1=`echo $name | cut -d '-' -f2`
	pop2=`echo $name | cut -d '-' -f3`
	echo $name
	
	#echo "Sum across each 10Mb bin"
	
	#awk '{for(i=1;i<=NF;i++)$i=(a[i]+=$i)}END{print}' $i > $name.total.sfs
	
	echo "Prep for window analysis"

	realSFS fst index \
	SAF/pop-Picoides_pubescens-$pop1.bamlist.saf.idx \
	SAF/pop-Picoides_pubescens-$pop2.bamlist.saf.idx \
	-sfs $name.total.sfs \
	-fold 1 \
	-fstout $name

	#echo "Get global estimate"

	#realSFS fst stats \
	#$name.fst.idx 

	echo "Sliding Windows"

	realSFS fst stats2 \
	$name.fst.idx \
	-win 50000 -step 10000 > $name.fst_slidingwindow

done
