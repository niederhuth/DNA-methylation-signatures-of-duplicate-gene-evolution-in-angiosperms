#!/bin/bash --login
#SBATCH --time=3:55:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50GB
#SBATCH --job-name MCScanX-transposed
#SBATCH --output=job_reports/%x-%j.SLURMout

cd ${PBS_O_WORKDIR}
export PATH="$HOME/miniconda3/envs/gene-duplication/bin:${PATH}"
export LD_LIBRARY_PATH="$HOME/miniconda3/envs/gene-duplication/lib:${LD_LIBRARY_PATH}"

#Set Variables
sample=$(pwd | sed s/^.*\\///)
species=$(awk -v FS="," -v a=${sample} '$1 == a' ../../misc/genomes.csv | cut -d ',' -f 6)
outgroup_list=$(echo ${species} | tr ' ' ',')
epochs=$(awk -v FS="," -v a=${sample} '$1 == a' ../../misc/genomes.csv | cut -d ',' -f 7)

#Say a bunch of stuff
echo "Sample: ${sample}"
echo "Outgroups: ${species}"
echo "Epochs: ${epochs}"
echo "Copying ${sample} data"

#Create output dirs
mkdir mcscanx
mkdir mcscanx/data

#Copy over data
cp ref/mcscanx/${sample}.gff mcscanx/data/${sample}.gff
cp diamond/${sample}-${sample}.m8 mcscanx/data/${sample}_unfiltered.blast
cp ref/mcscanx/${sample}_orthogroups.tsv dupgen/data/
for outgroup in $species
do
	echo "Copying ${outgroup} data"
	cat ref/mcscanx/${sample}.gff ../${outgroup}/ref/mcscanx/${outgroup}.gff > mcscanx/data/${sample}_${outgroup}.gff
	cat diamond/${sample}-${outgroup}.m8 ../${outgroup}/diamond/${outgroup}-${sample}.m8 > mcscanx/data/${sample}_${outgroup}_unfiltered.blast
	cp ../${outgroup}/ref/mcscanx/${outgroup}_orthogroups.tsv dupgen/data/mcscanx/data
done

#Filter blast hits based on orthogroup
cd mcscanx/data
join -1 2 -2 1 <(sort -k2,2 ${sample}_unfiltered.blast) <(sort -k1,1 ${sample}_orthogroups.tsv) | \
tr ' ' '\t' > tmp
join -1 2 -2 1 <(sort -k2,2 tmp) <(sort -k1,1 ${sample}_orthogroups.tsv) | tr ' ' '\t' | \
awk -v OFS="\t" '{if ($13==$14) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' | sort | uniq > ${sample}.blast
rm tmp
for outgroup in ${species}
do
	cat ${sample}_orthogroups.tsv ${outgroup}_orthogroups.tsv > tmp
	join -1 2 -2 1 <(sort -k2,2 ${sample}_${outgroup}_unfiltered.blast) <(sort -k1,1 tmp) | \
        tr ' ' '\t' > tmp2
        join -1 2 -2 1 <(sort -k2,2 tmp2) <(sort -k1,1 tmp) | tr ' ' '\t' | \
        awk -v OFS="\t" '{if ($13==$14) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' | sort | uniq > ${sample}_${outgroup}.blast
        rm tmp tmp2
done
cd ..

#Run MCScanX-transposed
echo "Running MCScanX-transposed"
cp ../../../scripts/MCScanX-transposed.tar.gz ./
tar -xjvf MCScanX-transposed.tar.gz
mv MCScanX-transposed/* ./
rmdir MCScanX-transposed
perl MCScanX-transposed.pl \
	-i data \
	-t "$sample" \
	-c "$outgroup_list" \
	-o results \
	-x "$epochs"

rm *

cd results
head -1 "$sample".transposed.pairs | \
awk -v OFS="\t" '{print $0,"epoch","epoch_species"}' >"$sample".transposed_epoch.pairs
count=1
for i in $species
do
	if [ $count -lt 2 ]
	then
		sed '1d' "$sample".transposed_after_"$i".pairs | \
		awk -v b=$count -v c=$i -v OFS="\t" '{print $0,b,c}' >> "$sample".transposed_epoch.pairs
		count=`expr $count + 1`
		a=$i
	else
				
		sed '1d' "$sample".transposed_between_"$a"_"$i".pairs | \
		awk -v b=$count -v c=$i -v OFS="\t" '{print $0,b,c}' >> "$sample".transposed_epoch.pairs
		count=`expr $count + 1`
		a=$i
	fi
done

echo "Done"
