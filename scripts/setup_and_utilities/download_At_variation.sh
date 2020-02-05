#!/bin/bash --login
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=3GB
#SBATCH --job-name download_At_variation
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR
mkdir variation
cd variation

#Download datasets
echo "Downloading data"
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE43nnn/GSE43857/suppl/GSE43857_RAW.tar

#Untar data
echo "Unpack tar files"
tar -xvf GSE43857_RAW.tar

#make sample directories
echo "Creating sample directories"
for i in GSM*tsv.gz
do
	name=$(echo $i | sed s/_.*//)
	mkdir $name
	mv $i $name
done

#Remove tar files
#rm GSE43857_RAW.tar

echo "Done"


