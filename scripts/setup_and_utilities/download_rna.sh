#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100GB
#SBATCH --job-name download_rna
#SBATCH --output=job_reports/%x-%j.SLURMout

cd ${PBS_O_WORKDIR}
export PATH="${HOME}/miniconda3/envs/gene-duplication/bin:${PATH}"
export LD_LIBRARY_PATH="${HOME}/miniconda3/envs/gene-duplication/lib:${LD_LIBRARY_PATH}"

#Set tmp directories
export TMPDIR=${PBS_O_WORKDIR}
export TMP=${PBS_O_WORKDIR}
export TEMP=${PBS_O_WORKDIR}

#List Variables
species=$(pwd | sed s/^.*\\///)
path1=$(pwd | sed s/data.*// | sed s/$/misc/)
samples=$(cut -f 5 ${path1}/${species}_rna.tsv)

#Create working directory
if [ -d rna ]
then
	echo "Working directory exists"
	cd rna
else
	mkdir rna
	cd rna
fi

#Download rna data
for i in ${samples}
do
	if ls ${i}/*_1.fastq.gz >/dev/null 2>&1
	then
		echo "${i} data found."
		echo "If this is a mistake, please delete ${i} and resubmit"
	else
		mkdir ${i}
		cd ${i}
		echo "Downloading SRA files for ${i}"
		sra_list=$(awk -v a=${i} '$5==a' | cut -f4)
		for sra in ${sra_list}
		do
			echo "Downloading ${sra}"
			prefetch --max-size 100000000 ${sra}
			mv ${sra}/*sra ./
			rmdir ${sra}
			echo "Converting ${sra}"
			fastq-dump --gzip --split-3 ${sra}.sra
			rm ${sra}.sra
		done
		cd ../
	fi
done

echo "Done"