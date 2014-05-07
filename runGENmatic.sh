#!/bin/bash
RSCRIPT_FILE=$SCRATCH/$1/$2.sh
echo "#! /bin/bash" > $RSCRIPT_FILE
echo "module load intel R/3.0.1" >> $RSCRIPT_FILE
echo "Rscript $SCRATCH/$1/$2.R $SCRATCH/ $1/ $3" >> $RSCRIPT_FILE

qsub -v foldername=$1 -q debug -l nodes=1:ppn=1,walltime=2:00:00 $RSCRIPT_FILE

#chmod 711 ./Script.sh
#./Script.sh GENmatictest GENmatic.R
