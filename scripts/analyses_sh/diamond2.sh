#!/bin/bash --login
#SBATCH --time=3:55:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=50GB
#SBATCH --job-name diamond2
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR
export PATH="$HOME/miniconda3/envs/gene-duplication/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/miniconda3/envs/gene-duplication/lib:$LD_LIBRARY_PATH"

#Set Variables
threads=10
evalue=0.00001

#Run diamond
x=$(pwd | sed s/^.*\\///)
species=$(cut -d ',' -f1 ../../misc/genomes.csv | sed '1d' | tr '\n' ' ')
mkdir diamond2
cd diamond2

for i in $species
do

	echo "$x $i blastp"
	diamond blastp \
		--threads $threads \
		--db ../../$i/ref/mcscanx/"$i".dmnd \
		--query ../ref/mcscanx/"$x"-protein.fa \
		--out "$x"-"$i".m8 \
		--un "$x"-"$i"-un.fa \
 		--more-sensitive \
		--evalue $evalue \
		--unal 0

done

