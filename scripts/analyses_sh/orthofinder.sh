#!/bin/bash --login
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=128GB
#SBATCH --job-name orthofinder
#SBATCH --output=%x-%j.SLURMout

cd $PBS_O_WORKDIR
export PATH="$HOME/miniconda3/envs/orthofinder/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/miniconda3/envs/orthofinder/lib:$LD_LIBRARY_PATH"

#Set Variables
species=$(cut -d ',' -f1 ../../misc/genomes.csv | sed '1d' | tr '\n' ' ')
threads=128
threads2=4

#Run OrthoFinder
echo "Copying sequence files"
mkdir seqs
for i in $species
do
	cp ../"$i"/ref/mcscanx/"$i"-protein.fa  seqs/"$i".fa
done

echo "Running OrthoFinder"
orthofinder \
	-t $threads \
	-a $threads2 \
	-M dendroblast \
	-S diamond \
	-o orthofinder \
	-f seqs/

echo "Create orthogroup list"
while read line
do
	og=$(echo $line | cut -d ' ' -f1 | sed s/\\:$//)
	echo $line | tr ' ' '\n' | sed '1d' | awk -v OFS='\t' -v a=$og '{print $0,a}' >> orthogroup_list.tsv
done < orthofinder/*/Orthogroups/Orthogroups.txt

echo "Get orthogroups for each species & gene"
for i in $species
do
	echo $i
	cut -f2 ../"$i"/ref/mcscanx/"$i".gff > tmp
	fgrep -f tmp orthogroup_list.tsv > ../"$i"/ref/mcscanx/"$i"_orthogroups.tsv
done
rm tmp

echo "Done"

