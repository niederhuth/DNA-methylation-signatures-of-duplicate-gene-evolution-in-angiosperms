#!/bin/bash --login
#SBATCH --time=3:55:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=60GB
#SBATCH --job-name lineage_specific_diamond
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR
export PATH="$HOME/miniconda3/envs/gene-duplication/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/miniconda3/envs/gene-duplication/lib:$LD_LIBRARY_PATH"

#Set Variables
threads=10
evalue=0.00001
evalue2=0.0000000001

#Run diamond
species=$(pwd | sed s/^.*\\///)
cd diamond2/lineage_specific

echo "BLASTP against 1kp proteins"
diamond blastp \
	--threads $threads \
	--db ../../../../misc/1kp.dmnd \
	--query "$species"-specific-protein.fa \
	--out 1kp.m8 \
	--un 1kp-un.fa \
	--more-sensitive \
	--masking 0 \
	--evalue $evalue \
	--unal 0

echo "BLASTP against trep proteins"
diamond blastp \
        --threads $threads \
        --db ../../../../misc/trep_proteins.dmnd \
        --query "$species"-specific-protein.fa \
        --out trep.m8 \
        --un trep-un.fa \
        --more-sensitive \
	--masking 0 \
        --evalue $evalue \
        --unal 0

echo "BLASTN against mipsREdat nucleotides"
diamond blastp \
        --threads $threads \
        --db ../../../../misc/mipsREdat.dmnd \
        --query "$species"-specific-cds.fa \
        --out mipsREdat.m8 \
        --un mipsREdat-un.fa \
        --more-sensitive \
	--masking 0 \
        --evalue $evalue2 \
        --unal 0


