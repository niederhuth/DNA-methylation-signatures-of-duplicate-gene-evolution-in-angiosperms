#!/bin/bash --login
#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10GB
#SBATCH --job-name get_baseline
#SBATCH --output=%x-%j.SLURMout

cd $PBS_O_WORKDIR

export PATH="$HOME/miniconda3/envs/gene-duplication/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/miniconda3/envs/gene-duplication/lib:$LD_LIBRARY_PATH"

#Get baseline values for DNA methylation contexts
samples=$(awk -v FS="," '$5=="yes"' ../misc/genomes.csv | cut -d ',' -f 1 | tr '\n' ' ')
header=$(echo $samples | cut -d ' ' -f 1)
head -1 $header/methylpy/results/CDS_methylation.tsv > tmp
for i in $samples
do
	sed '1d' $i/methylpy/results/CDS_methylation.tsv >> tmp
done

python ../scripts/analyses_py/get_baseline.py

rm tmp tmp2

