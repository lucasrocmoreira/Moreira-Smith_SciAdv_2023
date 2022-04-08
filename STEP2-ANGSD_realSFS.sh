#!/bin/bash
#PBS -l select=1:ncpus=1
#PBS -l walltime=5000:00:00
#PBS -N ANGSD_realSFS
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

#echo PP-AK-SE.sfs

realSFS SAF/pop-Picoides_pubescens-AK.bamlist.saf.idx SAF/pop-Picoides_pubescens-SE.bamlist.saf.idx -fold 1 -nSites 10000000 > PP-AK-SE.sfs2

#echo PP-AK-NE.sfs

#realSFS SAF/pop-Picoides_pubescens-AK.bamlist.saf.idx SAF/pop-Picoides_pubescens-NE.bamlist.saf.idx -fold 1 -nSites 10000000 > PP-AK-NE.sfs

#echo PP-AK-MW.sfs

#realSFS SAF/pop-Picoides_pubescens-AK.bamlist.saf.idx SAF/pop-Picoides_pubescens-MW.bamlist.saf.idx -fold 1 -nSites 10000000 > PP-AK-MW.sfs

#echo PP-AK-SR.sfs

#realSFS SAF/pop-Picoides_pubescens-AK.bamlist.saf.idx SAF/pop-Picoides_pubescens-SR.bamlist.saf.idx -fold 1 -nSites 10000000 > PP-AK-SR.sfs

#echo PP-AK-NR.sfs

#realSFS SAF/pop-Picoides_pubescens-AK.bamlist.saf.idx SAF/pop-Picoides_pubescens-NR.bamlist.saf.idx -fold 1 -nSites 10000000 > PP-AK-NR.sfs

#echo PP-AK-NW.sfs

#realSFS SAF/pop-Picoides_pubescens-AK.bamlist.saf.idx SAF/pop-Picoides_pubescens-NW.bamlist.saf.idx -fold 1 -nSites 10000000 > PP-AK-NW.sfs

#echo PP-NE-SE.sfs

#realSFS SAF/pop-Picoides_pubescens-NE.bamlist.saf.idx SAF/pop-Picoides_pubescens-SE.bamlist.saf.idx -fold 1 -nSites 10000000 > PP-NE-SE.sfs

#echo PP-NE-MW.sfs

#realSFS SAF/pop-Picoides_pubescens-NE.bamlist.saf.idx SAF/pop-Picoides_pubescens-MW.bamlist.saf.idx -fold 1 -nSites 10000000 > PP-NE-MW.sfs

#echo PP-NE-SR.sfs

#realSFS SAF/pop-Picoides_pubescens-NE.bamlist.saf.idx SAF/pop-Picoides_pubescens-SR.bamlist.saf.idx -fold 1 -nSites 10000000 > PP-NE-SR.sfs

#echo PP-NE-NR.sfs

#realSFS SAF/pop-Picoides_pubescens-NE.bamlist.saf.idx SAF/pop-Picoides_pubescens-NR.bamlist.saf.idx -fold 1 -nSites 10000000 > PP-NE-NR.sfs

#echo PP-NE-NW.sfs

#realSFS SAF/pop-Picoides_pubescens-NE.bamlist.saf.idx SAF/pop-Picoides_pubescens-NW.bamlist.saf.idx -fold 1 -nSites 10000000 > PP-NE-NW.sfs

#echo PP-MW-SE.sfs

#realSFS SAF/pop-Picoides_pubescens-MW.bamlist.saf.idx SAF/pop-Picoides_pubescens-SE.bamlist.saf.idx -fold 1 -nSites 10000000 > PP-MW-SE.sfs

#echo PP-MW-SR.sfs

#realSFS SAF/pop-Picoides_pubescens-MW.bamlist.saf.idx SAF/pop-Picoides_pubescens-SR.bamlist.saf.idx -fold 1 -nSites 10000000 > PP-MW-SR.sfs

#echo PP-MW-NR.sfs

#realSFS SAF/pop-Picoides_pubescens-MW.bamlist.saf.idx SAF/pop-Picoides_pubescens-NR.bamlist.saf.idx -fold 1 -nSites 10000000 > PP-MW-NR.sfs

#echo PP-MW-NW.sfs

#realSFS SAF/pop-Picoides_pubescens-MW.bamlist.saf.idx SAF/pop-Picoides_pubescens-NW.bamlist.saf.idx -fold 1 -nSites 10000000 > PP-MW-NW.sfs

#echo PP-SE-SR.sfs

#realSFS SAF/pop-Picoides_pubescens-SE.bamlist.saf.idx SAF/pop-Picoides_pubescens-SR.bamlist.saf.idx -fold 1 -nSites 10000000 > PP-SE-SR.sfs

#echo PP-SE-NR.sfs

#realSFS SAF/pop-Picoides_pubescens-SE.bamlist.saf.idx SAF/pop-Picoides_pubescens-NR.bamlist.saf.idx -fold 1 -nSites 10000000 > PP-SE-NR.sfs

#echo PP-SE-NW.sfs

#realSFS SAF/pop-Picoides_pubescens-SE.bamlist.saf.idx SAF/pop-Picoides_pubescens-NW.bamlist.saf.idx -fold 1 -nSites 10000000 > PP-SE-NW.sfs

#echo PP-SR-NR.sfs

#realSFS SAF/pop-Picoides_pubescens-SR.bamlist.saf.idx SAF/pop-Picoides_pubescens-NR.bamlist.saf.idx -fold 1 -nSites 10000000 > PP-SR-NR.sfs

#echo PP-SR-NW.sfs

#realSFS SAF/pop-Picoides_pubescens-SR.bamlist.saf.idx SAF/pop-Picoides_pubescens-NW.bamlist.saf.idx -fold 1 -nSites 10000000 > PP-SR-NW.sfs

#echo PP-NR-NW.sfs

#realSFS SAF/pop-Picoides_pubescens-NR.bamlist.saf.idx SAF/pop-Picoides_pubescens-NW.bamlist.saf.idx -fold 1 -nSites 10000000 > PP-NR-NW.sfs

