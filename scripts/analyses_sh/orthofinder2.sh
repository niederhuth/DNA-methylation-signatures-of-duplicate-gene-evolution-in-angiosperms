#!/bin/bash --login
#SBATCH --time=168:00:00
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
threads=128
threads2=3
species=$(cut -d ',' -f1 ../../misc/genomes.csv | sed '1d' | tr '\n' ' ')
search="diamond"
msa="mafft"
tree="fasttree"
inflation=1.5
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
	-M msa \
	-S $search \
	-A $msa \
	-T $tree \
	-I $inflation \
	-o orthofinder \
	-f seqs/

echo "Create orthogroup list"
while read line
do
	og=$(echo $line | cut -d ' ' -f1 | sed s/\\:$//)
	echo $line | tr ' ' '\n' | sed '1d' | awk -v OFS='\t' -v a=$og '{print $0,a}' >> orthogroup_list.tsv
done < */Orthogroups/Orthogroups.txt

echo "Done"

