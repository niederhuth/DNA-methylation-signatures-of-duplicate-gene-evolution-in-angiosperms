#!/bin/bash --login
#SBATCH --time=3:49:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#SBATCH --job-name GC_content
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR
export PATH="$HOME/miniconda3/envs/gene-duplication/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/miniconda3/envs/gene-duplication/lib:$LD_LIBRARY_PATH"

#variables
sample=$(pwd | sed s/.*data\\/// | sed s/\\/.*//)

#get total weighted mC
echo "Get GC content of $sample"
python ../../scripts/analyses_py/GC_content.py $sample

