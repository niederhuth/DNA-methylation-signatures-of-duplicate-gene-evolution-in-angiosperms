#!/bin/bash --login
#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#SBATCH --job-name lineage_specific
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR
export PATH="$HOME/miniconda3/envs/gene-duplication/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/miniconda3/envs/gene-duplication/lib:$LD_LIBRARY_PATH"

species=$(pwd | sed s/.*data\\/// | sed s/\\/.*//)
genomes=$(sed '1d' ../../misc/genomes.csv | cut -d ',' -f1 | tr '\n' ' ')
cd diamond2
mkdir lineage_specific
for i in *fa
do
	name=$(echo $i | sed s/"$species"-// | sed s/\-un.fa//)
	grep \> $i | sed s/^\>// | sed s/\ .*$// | awk -v a=$name -v OFS="\t" '{print $0,a}' > lineage_specific/"$name"-unmapped.txt
done
cd lineage_specific
number=$(expr $(ls *txt | wc -l) - 1)

echo "Feature,Missing,$(echo $genomes | tr ' ' ',')" > unmapped_genes.csv

for i in $(cat *txt | cut -f1 | sort | uniq)
do
	a=$(grep $i *txt | wc -l)
	b=""
	for x in $genomes
	do
		if grep -q $i $x-unmapped.txt
		then
			b="$b,"
		else
			b="$b,$x"
		fi
	done
	echo $i $a $b | tr ' ' ',' >> unmapped_genes.csv
done

#Report on results
echo "Total genes: " $(sed '1d' unmapped_genes.csv | wc -l)
echo "Species-specific genes: " $(awk -v a=$number -v FS="," '$2>=a' unmapped_genes.csv | grep "$species" | wc -l)

#Make list of species specific genes
awk -v a=$number -v FS="," '$2>=a' unmapped_genes.csv | cut -d ',' -f1 > tmp

#Get protein sequences
fasta_formatter -w 0 -i ../../ref/mcscanx/"$species"-protein.fa -o tmp2
fgrep -A 1 -f tmp tmp2 | grep -v \\-\\- > "$species"-specific-protein.fa

#Get CDS sequences
fasta_formatter -w 0 -i ../../ref/mcscanx/"$species"-cds.fa -o tmp2
fgrep -A 1 -f tmp tmp2 | grep -v \\-\\- > "$species"-specific-cds.fa
rm tmp tmp2

echo "Done"






