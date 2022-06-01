#! /usr/bin/python

#Original by Tiffany Liu 11/17/2016
#Modified by Chad Niederhuth 6/1/2022

import argparse
import os.path
import sys

parser = argparse.ArgumentParser(description="This program creates a gene list that does not contain TEs.")
parser.add_argument("--input_file_TEpfam", help="Path to input file: list of TE-related Pfam IDs.", required = True)
parser.add_argument("--input_file_maxPfam", help="Path to input file: hmmscan output (e.g. maxPfam.txt).", required = True)
parser.add_argument("--input_file_geneList_toKeep", help="Path to input file: genes to keep (e.g. maker standard gene list).", required = True)
parser.add_argument("--input_file_TEhmm", help="Path to input file: hmmscan output of TE-related genes (e.g. gypsy hmm)")
parser.add_argument("--input_file_TEblast", help="Path to input file: TE-related blast results (e.g. blastp with transposases).")
parser.add_argument("--output_file", help="Path to the output file.", required=True)
args = parser.parse_args()

input_file_TEpfam = args.input_file_TEpfam
if os.path.isfile(input_file_TEpfam) == False:
    print('\n\nThe file ' + input_file_TEpfam + ' does not exist.\n')
    sys.exit()
input_file_TEhmm = args.input_file_TEhmm
if os.path.isfile(input_file_TEhmm) == False:
    print('\n\nThe file ' + input_file_TEhmm + ' does not exist.\n')
    sys.exit()
input_file_maxPfam = args.input_file_maxPfam
if os.path.isfile(input_file_maxPfam) == False:
    print('\n\nThe file ' + input_file_maxPfam + ' does not exist.\n')
    sys.exit()
input_file_geneList_toKeep = args.input_file_geneList_toKeep
if os.path.isfile(input_file_geneList_toKeep) == False:
    print('\n\nThe file ' + input_file_geneList_toKeep + ' does not exist.\n')
    sys.exit()
input_file_TEblast = args.input_file_TEblast
if os.path.isfile(input_file_TEblast) == False:
    print('\n\nThe file ' + input_file_TEblast + ' does not exist.\n')
    sys.exit()

output_file = args.output_file
if os.path.isfile(output_file):
    print('\n\nThe file ' + output_file + ' already exists.\n')
    sys.exit()
output_fh = open(output_file, 'w')

#This loop adds all TE-related PfamIDs into a list
TEpfamList = []
TEpfamSet = set(TEpfamList)
with open(input_file_TEpfam) as input_fh_TEpfam:
    for each_line in input_fh_TEpfam:
        if 'Pfam' == each_line[0:4]:
            pass #Ignoring header
        else:
            line_string = each_line.strip()
            (Pfam,Domain,Description) = line_string.split('\t')
            TEpfamSet.add(Pfam)
print("Total number of TE-related pfam domains: ", len(TEpfamSet))

#This loop compares PfamIDs between the TE pfam list and the maxPfam file,
#resulting in a list of TE-related PfamIDs that were present in the maxPfam file
pfamHits = []
pfamSet = set(pfamHits)
TEgeneList = []
TEgeneSet = set(TEgeneList)
with open(input_file_maxPfam) as input_fh_maxPfam:
    for each_line in input_fh_maxPfam:
        if "#" not in each_line[0]: #ignoring header
            line_string = each_line.strip()
            line_tuple = line_string.split()
            pfam = line_tuple[1]
            pfamSplit = pfam.split(".")
            pfamID = pfamSplit[0]
            geneID = line_tuple[2]
            if pfamID in TEpfamSet:
                TEgeneSet.add(geneID)
                pfamSet.add(geneID)
output_pfam = open("pfam_filtered_genes.txt", 'w') 
for each_element in pfamSet:
    output_pfam.write("%s\n"%(each_element))
output_pfam.close()    
print("Number of TE-related genes from pfam hmmscan: ", len(pfamSet))

#Add the genes identified from the gypsy hmm
gypsyHits = []
gypsySet = set(gypsyHits)
with open(input_file_TEhmm) as input_fh_TEhmm:
    for each_line in input_fh_TEhmm:
        if "#" not in each_line[0]:
            line_string = each_line.strip()
            line_tuple = line_string.split()
            gene = line_tuple[2]
            TEgeneSet.add(gene)
            gypsySet.add(gene)
output_gypsy = open("gypsy_filtered_genes.txt", 'w') 
for each_element in gypsySet:
    output_gypsy.write("%s\n"%(each_element))
output_gypsy.close()    
print("Number of TE-related genes from gypsy hmmscan: ", len(gyspySet))            

#Add blast results (outfmt 6 (table))
blastHits = []
blastSet = set(gypsyHits)
with open(input_file_TEblast) as input_fh_TEblast:
    for each_line in input_fh_TEblast:
        if "#" not in each_line[0]:
            line_string = each_line.strip()
            line_tuple = line_string.split()
            gene = line_tuple[0]
            TEgeneSet.add(gene)
            blastSet.add(gene)
output_blast = open("blast_filtered_genes.txt", 'w') 
for each_element in blastSet:
    blast.write("%s\n"%(each_element))
output_blast.close()    
print("Number of TE-related genes from blast: ", len(gyspySet)) 
print("Number of TE-related genes filtered: ", len(TEgeneSet))

#This loop adds all genes from the gene list file to a list
oldGeneList = []
oldGeneSet = set(oldGeneList)
with open(input_file_geneList_toKeep) as input_fh_geneList_toKeep:
    for each_line in input_fh_geneList_toKeep:
        line_string = each_line.strip()
        oldGeneSet.add(line_string)

#Below: what is present in oldGeneListSet, absent in TEgeneSet
newGeneSet = oldGeneSet.difference(TEgeneSet)
print("Length of old gene set: ", len(oldGeneSet))
print("Length of new gene set: ", len(newGeneSet))
for each_element in newGeneSet:
    output_fh.write("%s\n"%(each_element))

output_fh.close()








