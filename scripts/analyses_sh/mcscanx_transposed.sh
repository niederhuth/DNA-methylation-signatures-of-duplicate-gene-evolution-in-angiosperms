#!/bin/bash --login
#SBATCH --time=3:55:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50GB
#SBATCH --job-name MCScanX-transposed
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/atac/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/atac/lib:${LD_LIBRARY_PATH}"

#Other variables, these should not have to be changed, but should set automatically
path1=$(pwd | sed s/data.*// | sed s/$/scripts/)
sample=$(pwd | sed s/^.*\\///)
species=$(awk -v FS="," -v a=${sample} '$1 == a' ../../misc/genomes.csv | cut -d ',' -f 6)
outgroup_list=$(echo ${species} | tr ' ' ',')
epochs=$(awk -v FS="," -v a=${sample} '$1 == a' ../../misc/genomes.csv | cut -d ',' -f 7)

#Say a bunch of stuff
echo "Sample: ${sample}"
echo "Outgroups: ${species}"
echo "Epochs: ${epochs}"
echo "Copying ${sample} data"

#Create directories for MCScanX-transposed
mkdir mcscanx mcscanx/data

#Copy over BLAST, gff, and orthogroup files for the sample 
cp ref/mcscanx/${sample}.gff mcscanx/data/${sample}.gff
cp diamond/${sample}-${sample}.m8 mcscanx/data/${sample}_unfiltered.blast
cp ref/mcscanx/${sample}_orthogroups.tsv mcscanx/data/
#Iterate over the list of species and copy over BLAST, gff, and orthogroup files
for outgroup in $species
do
	echo "Copying ${outgroup} data"
	#Concatanate the gff files
	cat ref/mcscanx/${sample}.gff ../${outgroup}/ref/mcscanx/${outgroup}.gff \
	> mcscanx/data/${sample}_${outgroup}.gff
	#Concatanate the BLAST files
	cat diamond/${sample}-${outgroup}.m8 ../${outgroup}/diamond/${outgroup}-${sample}.m8 \
	> mcscanx/data/${sample}_${outgroup}_unfiltered.blast
	cp ../${outgroup}/ref/mcscanx/${outgroup}_orthogroups.tsv mcscanx/data/mcscanx/data
done

#Filter BLAST hits based on orthogroup for the sample data
cd mcscanx/data
join -1 2 -2 1 <(sort -k2,2 ${sample}_unfiltered.blast) <(sort -k1,1 ${sample}_orthogroups.tsv) | \
tr ' ' '\t' > tmp
join -1 2 -2 1 <(sort -k2,2 tmp) <(sort -k1,1 ${sample}_orthogroups.tsv) | tr ' ' '\t' | \
awk -v OFS="\t" '{if ($13==$14) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' | sort | uniq > ${sample}.blast
rm tmp
#For each outgroup species, filter BLAST hits based on orthogroup
for outgroup in ${species}
do
	#Cat the sample and outgroup orthogroup lists together. 
	#This is necessary since MCScanX-transposed requires reciprocal BLAST hits
	cat ${sample}_orthogroups.tsv ${outgroup}_orthogroups.tsv > tmp
	join -1 2 -2 1 <(sort -k2,2 ${sample}_${outgroup}_unfiltered.blast) <(sort -k1,1 tmp) | \
        tr ' ' '\t' > tmp2
        join -1 2 -2 1 <(sort -k2,2 tmp2) <(sort -k1,1 tmp) | tr ' ' '\t' | \
        awk -v OFS="\t" '{if ($13==$14) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' | sort | \
        uniq > ${sample}_${outgroup}.blast
        rm tmp tmp2
done
#Make sure to move up a directory
cd ..

#Copy the MCScanX-transposed scripts and untar/gz them
#This is necessary since I have been too lazy to modify the scripts myself so that this isn't necessary
cp ${path1}/MCScanX-transposed.tar.gz ./
tar -xjvf MCScanX-transposed.tar.gz
mv MCScanX-transposed/* ./
rmdir MCScanX-transposed

#Run MCScanX-transposed
echo "Running MCScanX-transposed"
perl MCScanX-transposed.pl \
	-i data \
	-t "${sample}" \
	-c "${outgroup_list}" \
	-o results \
	-x "${epochs}"
#Remove all the junk...yeah, I know this step is super dangerous, sue me
#rm *

#Now we are going to create a table of transposed pairs and which "epoch" they transposed in
cd results
#Create the header
head -1 ${sample}.transposed.pairs | \
awk -v OFS="\t" '{print $0,"epoch","epoch_species"}' > ${sample}.transposed_epoch.pairs
#The count is used to define the "epoch", starting with the most recent split
count=1
#Iterate over each of the outgroup species
for i in ${species}
do
	#Start with the "transposed_after" pairs, the ones not present in any, but the target species
	#There is probably a more elegant way to do this, but it works, since the files are not named 
	#The same. So this first "if" ensures that the first file we look for is the "transposed_after"
	if [ ${count} -lt 2 ]
	then
		#Add those to our list, with "1" as the epoch
		sed '1d' ${sample}.transposed_after_${i}.pairs | \
		awk -v b=${count} -v c=${i} -v OFS="\t" '{print $0,b,c}' >> ${sample}.transposed_epoch.pairs
		#Once that is done, add one to our count, so that we can move on
		count=`expr ${count} + 1`
		#Remember which species we just looked at, so we can use it in next iteration
		a=${i}
	else
		#Now we look at those events that occured between the last species in our list and our current species
		sed '1d' ${sample}.transposed_between_${a}_${i}.pairs | \
		awk -v b=${count} -v c=${i} -v OFS="\t" '{print $0,b,c}' >> ${sample}.transposed_epoch.pairs
		count=`expr ${count} + 1`
		a=${i}
	fi
done

echo "Done"
