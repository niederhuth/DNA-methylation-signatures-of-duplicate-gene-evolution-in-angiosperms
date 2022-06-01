#!/bin/bash --login
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#SBATCH --job-name pain_in_my_ass
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR
cd dupgen/results/kaks_results

cut -f1 ../../results-unique/classified_genes.tsv > tmp


for i in Fxananassa.*.kaks.KKC.format
do
	a=$(echo $i | sed s/\.KKC.format//)
	head -1 $a > tmp2
	while read line
	do
		grep ^$line\\- $i | sed "s/$line\-/$line,/" | cut -f1,3,4,5,6 | tr ',' '\t' >> tmp2
	done < tmp
	mv tmp2 $a
done

