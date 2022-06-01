#!/bin/bash --login
#SBATCH --time=3:55:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#SBATCH --job-name MCScanX
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR
SCRIPTS=$(pwd | sed s/methylation.*/methylation\\/scripts/)
export PATH="$HOME/miniconda3/envs/gene-duplication/bin:$PATH"
export PATH="$SCRIPTS/comparative/DupGen_finder:$PATH"
export LD_LIBRARY_PATH="$HOME/miniconda3/envs/gene-duplication/lib:$LD_LIBRARY_PATH"

#Set Variables
sample=$(pwd | sed s/^.*\\///)
outgroup=$(awk -v FS="," -v a="$sample" '$1 == a' ../../misc/genomes.csv | cut -d ',' -f 6)

#Run MCScanX
cd dupgen
mkdir collinearity

for i in $outgroup
do
	echo "Outgroup: $i"
	cat ../ref/mcscanx/"$sample".gff ../../"$i"/ref/mcscanx/"$i".gff > data/"$sample"_"$i".gff
	cp ../diamond/"$sample"-"$i".m8 data/"$sample"_"$i".blast

	MCScanX data/"$sample"_"$i" collinearity/"$sample"_"$i"
done

cd collinearity
echo "$sample $outgroup" | tr ' ' ',' > collinearity.csv
cat ../data/"$sample".gff | cut -f2 | while read line
do
	a="$line"
	for i in $outgroup
	do
		if grep -q $line "$sample"_"$i".collinearity
		then
			a=$(echo "$a,$i")
		else
			a=$(echo "$a,")
		fi
	done
	echo $a >> collinearity.csv
done

echo "Done"
