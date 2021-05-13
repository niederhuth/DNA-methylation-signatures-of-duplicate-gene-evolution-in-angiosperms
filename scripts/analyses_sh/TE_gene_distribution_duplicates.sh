#!/bin/bash --login
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB
#SBATCH --job-name TE_gene_distribution_duplicates
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR
export PATH="$HOME/miniconda3/envs/gene-duplication/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/miniconda3/envs/gene-duplication/lib:$LD_LIBRARY_PATH"

#Set tmp directories
export TMPDIR=$PBS_O_WORKDIR
export TMP=$PBS_O_WORKDIR
export TEMP=$PBS_O_WORKDIR

#variables
sample=$(pwd | sed s/.*data\\/// | sed s/\\/.*//)

#Filter methylation data to keep only duplicate genes
awk '$2!="unclassified"' dupgen/results-unique/classified_genes.tsv | cut -f1 > tmp
fgrep -w -f tmp methylpy/results/${sample}_classified_genes.tsv > methylpy/results/${sample}_classified_genes_duplicates.tsv
wc -l methylpy/results/${sample}_classified_genes.tsv
wc -l methylpy/results/${sample}_classified_genes_duplicates.tsv
wc -l tmp
rm tmp

#Get TE & gene distributions
echo "Get TE & gene distributions for $sample"
python ../../scripts/analyses_py/TE_gene_distribution_duplicates.py $sample

echo "Done"
