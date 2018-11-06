#!/bin/bash
#SBATCH --account=quinlan-kp
#SBATCH --partition=quinlan-kp
#SBATCH -o ./%j-%N.out
#SBATCH -e ./%j-%N.err
#SBATCH --time=24:00:00

# CCR for gnomAD:
bash newccrs.sh -c -w -s -v newgnomAD10x.5 -d 10 -d 0.5 -x data/segmental.bed.gz -x data/self-chains.id90.bed.gz -q X -q Y -u

# X-CCR for gnomAD:
bash newccrs.sh -c -w -s -v newXchrom -d 10 -d 0.5 -x data/segmental.bed.gz -x data/self-chains.id90.bed.gz -r X -u
