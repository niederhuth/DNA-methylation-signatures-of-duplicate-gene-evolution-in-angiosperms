#!/bin/bash --login
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB
#SBATCH --job-name At_CDS_methylation
#SBATCH --output=%x-%j.SLURMout

cd ${PBS_O_WORKDIR}
export PATH="${HOME}/miniconda3/envs/gene-duplication/bin:${PATH}"
export LD_LIBRARY_PATH="${HOME}/miniconda3/envs/gene-duplication/lib:${LD_LIBRARY_PATH}"

#Set tmp directories
export TMPDIR=${PBS_O_WORKDIR}
export TMP=${PBS_O_WORKDIR}
export TEMP=${PBS_O_WORKDIR}

#Set variables
sample=$(pwd | sed s/.*GSM/GSM/) 

#Fix header & chrs
header=$(zcat ${sample}*tsv.gz | head -1 | cut -f1)
if [ ${header} == "chrom" ]
then
	echo "Header line present, removing header line"
	zcat ${sample}*tsv.gz | sed '1d' | sed s/^/Chr/ | gzip -c > tmp.gz
else
	echo "No header line"
	zcat ${sample}*tsv.gz | sed s/^/Chr/ | gzip -c > tmp.gz
fi
	
#get total weighted mC
echo "Get gene CDS methylation data"
python ../../../../scripts/analyses_py/At_CDS_methylation_variation.py

rm tmp.gz allc.bed

echo "done" 

