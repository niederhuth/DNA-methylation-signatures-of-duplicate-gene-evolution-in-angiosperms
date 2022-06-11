#!/bin/bash --login
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100GB
#SBATCH --job-name=../job_reports/filter_genes
#SBATCH --output=%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=20

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/gene-duplication/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/gene-duplication/lib:$LD_LIBRARY_PATH"

#Set temporary directories for large memory operations
export TMPDIR=$(pwd)
export TMP=$(pwd)
export TEMP=$(pwd)

#The following shouldn't need to be changed, but should set automatically
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
path1=$(pwd | sed s/data.*/misc/)
path2=$(pwd | sed s/data.*/scripts/)
path3="gene_filtering"

#Make & cd to directory
if [ -d ${path3} ]
then
	cd ${path3}
else
	mkdir ${path3}
	cd ${path3}
fi

#Set temporary directories for large memory operations
export TMPDIR=$(pwd)
export TMP=$(pwd)
export TEMP=$(pwd)

#Generate a gene list from gff file
echo "Generating Gene List"
perl ${path2}/transposons/pl/create_gene_list.pl \
	--input_gff ../annotations/${species}.gff \
	--output_file input_gene_list.txt

#Check for Pfam A
echo "Downloading Pfam-A"
wget -q http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz

#Prepare Pfam A hmm-database
echo "Preparing Pfam A hmm database"
hmmpress Pfam-A.hmm

#Search Pfam A hmm domains
echo "Searching against Pfam A hmm profiles"
hmmscan \
	--domE 1e-5 \
	-E 1e-5 \
	--cpu ${threads} \
	-o pfam_alignments.out \
	--tblout prot_domains.out \
	Pfam-A.hmm \
	../mcscanx/${species}-protein.fa

#Download and make Transposase blast DB
echo "Downloading Tpases020812"
wget -q http://www.hrt.msu.edu/uploads/535/78637/Tpases020812.gz
gunzip Tpases020812.gz
#remove trailing white space in the Tpases020812 file
sed -i 's/[ \t]*$//' Tpases020812

echo "Making Transposase diamond blast DB"
diamond makedb \
	--threads ${threads} \
	--in Tpases020812 \
	--db Tpases020812.dmnd

#Run diamond blastp against Transposases 
echo "Running diamond blastp on Transposases"
diamond blastp \
	--threads ${threads} \
	--db Tpases020812.dmnd \
	--query ../mcscanx/${species}-protein.fa \
	--out TE_blast.out\
	--evalue 1e-10 \
	--outfmt 6

#Still running this, even though we will not use for the final analysis
#Download Gypsy DB hmm files and format the hmm database
echo "Downloading GyDB_collection"
wget -q https://gydb.org/extensions/Collection/collection/db/GyDB_collection.zip
unzip GyDB_collection.zip
echo "Combining GyDB hmm profiles"
cat GyDB_collection/profiles/*hmm > all_gypsy.hmm
echo "Formatting all_gypsy.hmm database"
hmmpress all_gypsy.hmm

#Search gypsy hmm profiles
echo "Searching against gypsy hmm profiles"
hmmscan \
	--domE 1e-5 \
	-E 1e-5 \
	--cpu ${threads} \
	-o gypsy_alignments.out \
	--tblout gypsyHMM_analysis.out \
	all_gypsy.hmm \
	../mcscanx/${species}-protein.fa

#Still running this, even though we will not use for the final analysis
#Create a genelist with no TEs
echo "Creating gene list with TEs removed"
python ${path2}/transposons/py/create_no_TE_genelist.py \
	--input_geneList input_gene_list.txt \
	--pfamhmm prot_domains.out \
	--TEpfam_list ${path1}/TE_Pfam_domains.txt \
	--TEblast TE_blast.out \
	--output_file noTE_gene_list.txt
	#Removed TEhmm against gypsy genes as this resulted in a high number of false positives
	#--TEhmm gypsyHMM_analysis.out

echo "Done"

