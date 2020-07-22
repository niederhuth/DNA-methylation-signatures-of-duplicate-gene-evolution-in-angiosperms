#!/bin/bash --login
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10GB
#SBATCH --job-name download_1kp
#SBATCH --output=%x-%j.SLURMout

cd $PBS_O_WORKDIR
export PATH="$HOME/miniconda3/envs/gene-duplication/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/miniconda3/envs/gene-duplication/lib:$LD_LIBRARY_PATH"

#Download 1kp Data
mkdir 1kp_tmp
cd 1kp_tmp
echo "Downloading 1kp data"
cat ../1kp_list.txt | while read line
do
	a=$(echo $line | sed s/\-.*//)
	wget ftp://parrot.genomics.cn/gigadb/pub/10.5524/100001_101000/100627/assemblies/"$line"/"$a"-translated-protein.fa.gz
done
zcat *-translated-protein.fa.gz > 1kp.fa

#Make diamond DB
echo "Making diamond database"
diamond makedb --in 1kp.fa --db ../1kp.dmnd --masking 0

#Remove tmp data
cd ../
rm -R 1kp_tmp

echo "Done"

