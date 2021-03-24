#!/bin/bash --login
#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#SBATCH --job-name job_reports/build_table
#SBATCH --output=%x-%j.SLURMout

cd $PBS_O_WORKDIR

export PATH="$HOME/miniconda3/envs/gene-duplication/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/miniconda3/envs/gene-duplication/lib:$LD_LIBRARY_PATH"

#Build Table
cd variation

cut -f1,30 ../methylpy/results/Athaliana_classified_genes.tsv | sed s/Classification/Reference/ > At_variation.tsv

for i in GSM*
do
	cut -f30 "$i"/"$i"_classified_genes.tsv | sed s/Classification/"$i"/ > tmp
	paste -d '\t' At_variation.tsv tmp > tmp2
	mv tmp2 At_variation.tsv
	rm tmp
done 

mv At_variation.tsv ../methylpy/results

echo "Done"

