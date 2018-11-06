#!/bin/bash
#SBATCH --account=quinlan-kp
#SBATCH --partition=quinlan-kp
#SBATCH -o ./%j-%N.out
#SBATCH -e ./%j-%N.err
#SBATCH --time=24:00:00

# CCR for gnomAD:
bash regions.sh -c -s -w -v gnomAD10x.5syn -d 10 -d 0.5 -x data/segmental.bed.gz -x data/self-chains.id90.bed.gz -q X -q Y -g -u

# CCR for ExAC v1:
bash regions.sh -c -s -w -v ExACv1syn -d 10 -d 0.5 -x data/segmental.bed.gz -x data/self-chains.id90.bed.gz -q X -q Y -e -u

# X-CCR for gnomAD:
bash regions.sh -c -s -w -v Xchromonly -d 10 -d 0.5 -x data/segmental.bed.gz -x data/self-chains.id90.bed.gz -r X -g -u
