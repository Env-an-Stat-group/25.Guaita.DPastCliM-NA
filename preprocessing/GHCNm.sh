#!/bin/bash
#$ -M pguaita@nd.edu
#$ -m abe
#$ -pe smp 24
#$ -q long@@crippa_d32cepyc
#$ -N matlabjob

export MATLABPATH=/data/pguaita/downscaling/matlab_code/preprocessing/

cd /data/pguaita/downscaling/matlab_code/preprocessing/

module load matlab

matlab -nodisplay -nosplash -r "GHCNm; exit;"
