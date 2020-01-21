#!/bin/bash --login
#SBATCH --time=3:55:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50GB
#SBATCH --job-name MCScanX-transposed
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR
export PATH="$HOME/miniconda3/envs/gene-duplication/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/miniconda3/envs/gene-duplication/lib:$LD_LIBRARY_PATH"

#Set Variables
sample=$(pwd | sed s/^.*\\///)
species=$(awk -v FS="," -v a="$sample" '$1 == a' ../../misc/genomes.csv | cut -d ',' -f 6)
outgroup_list=$(echo $species | tr ' ' ',')
epochs=$(awk -v FS="," -v a="$sample" '$1 == a' ../../misc/genomes.csv | cut -d ',' -f 7)

#Run MCScanX-transposed
mkdir mcscanx
mkdir mcscanx/data
echo "Sample: $sample"
echo "Outgroups: $species"
echo "Epochs: $epochs"
echo "Copying $sample data"
cp ref/mcscanx/"$sample".gff mcscanx/data/"$sample".gff
cp diamond/"$sample"-"$sample".m8 mcscanx/data/"$sample".blast
for i in $species
do
	echo "Copying $i data"
	cat ref/mcscanx/"$sample".gff ../"$i"/ref/mcscanx/"$i".gff > mcscanx/data/"$sample"_"$i".gff
	cat diamond/"$sample"-"$i".m8 ../"$i"/diamond/"$i"-"$sample".m8 > mcscanx/data/"$sample"_"$i".blast
done
cd mcscanx
echo "Running MCScanX"
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
