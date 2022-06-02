#!/bin/bash --login
#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10GB
#SBATCH --job-name classify_genes
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR

export PATH="$HOME/miniconda3/envs/gene-duplication/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/miniconda3/envs/gene-duplication/lib:$LD_LIBRARY_PATH"

#Set variables
sample=$(pwd | sed s/^.*\\///)

#Classify genes
cd methylpy
python ../../../scripts/methylC/py/classify_genes.py "$sample"

cut -f30 results/"$sample"_classified_genes.tsv | sort | uniq -c 
