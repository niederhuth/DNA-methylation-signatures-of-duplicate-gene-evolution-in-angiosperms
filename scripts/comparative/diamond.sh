#!/bin/bash --login
#SBATCH --time=3:55:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=100GB
#SBATCH --job-name diamond
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR
export PATH="$HOME/miniconda3/envs/gene-duplication/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/miniconda3/envs/gene-duplication/lib:$LD_LIBRARY_PATH"

#Set Variables
threads=50
evalue=0.00001
max_target_seqs=0

#Run diamond
x=$(pwd | sed s/^.*\\///)
species=$(cut -d ',' -f1 ../../misc/genomes.csv | sed '1d' | tr '\n' ' ')
mkdir diamond
cd diamond

for i in ${species}
do

	echo "${x} ${i} blastp"
	diamond blastp \
		--threads ${threads} \
		--db ../../${i}/ref/mcscanx/${i}.dmnd \
		--query ../ref/mcscanx/${x}-protein.fa \
		--out ${x}-${i}.m8 \
		--un ${x}-${i}-un.fa \
 		--more-sensitive \
		--evalue ${evalue} \
		--max-target-seqs ${max_target_seqs} \
		--unal 0

done

