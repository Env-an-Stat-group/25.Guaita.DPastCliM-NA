#!/bin/bash
#$ -M pguaita@nd.edu
#$ -m abe
#$ -pe smp 8
#$ -q long@@crippa_d32cepyc
#$ -N matlabjob

export MATLABPATH=/data/pguaita/downscaling/matlab_code/downscaling/

cd /data/pguaita/downscaling/matlab_code/downscaling/

module load matlab

matlab -nodisplay -nosplash -r "PCR_downscaling_gridding_v1_LR; exit;"
