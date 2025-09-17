#!/bin/bash
#$ -M pguaita@nd.edu
#$ -m abe
#$ -pe smp 6
#$ -q long@@crippa_d32cepyc
#$ -N matlabjob

export MATLABPATH=/data/pguaita/downscaling/matlab_code/visualization/

cd /data/pguaita/downscaling/matlab_code/visualization/

module load matlab

matlab -nodisplay -nosplash -r "plot_split; exit;"
