import sys
import pandas as pd
import pybedtools as pbt
from os import remove
from scipy import stats
from numpy import float64
from numpy import minimum
from subprocess import call
from io import TextIOWrapper
from itertools import product
from gzip import open as gzopen
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.SeqUtils import GC123
from Bio import motifs
from Bio.Seq import Seq

#split large files into temporary files of smaller size
def split_file(input,line_number=10000000):
	#set starting count values
	a = 1
	b = 0
	c = line_number
	#open first output file and add to list
	files = ['split' + str(a) + '.tmp']
	out=open('split' + str(a) + '.tmp','w')
	#open input file
	with TextIOWrapper(gzopen(input,'rb')) as e:
		#iterate over each line, adding an index
		for index, line in enumerate(e):
			#test if index fits between upper and lower limits and write to file
			if index <= c:
				if index > b:
					out.write(str(line))
			else:
				#close last ouput
				out.close()
				#reset count values
				a += 1
				b = c
				c += line_number
				#open new output and add to list
				filesappend('split' + str(a) + '.tmp')
				out=open('split' + str(a) + '.tmp','w')
				#output line
				out.write(str(line))
		#close last ouptut
		out.close()
	#return number of temporary files for use in other functions
	return(files)

#function for filtering annotation files based on feature (gene, exon, mRNA, etc)
def feature_filter(x,feature):
	if feature:
		return x[2] == feature
	else:
		return x

#function for filtering annotation files based on strand
def strand_filter(x,strand):
	return x.strand == strand

#function for filtering annotation files based on chromosome
def chr_filter(x,chr):
	return x.chrom in chr

#interpret sequence context, taken from methylpy utilities
def expand_nucleotide_code(mc_type=['C']):
	iub_dict = {'N':['A','C','G','T','N'],
				'H':['A','C','T','H'],
				'D':['A','G','T','D'],
				'B':['C','G','T','B'],
				'A':['A','C','G','A'],
				'R':['A','G','R'],
				'Y':['C','T','Y'],
				'K':['G','T','K'],
				'M':['A','C','M'],
				'S':['G','C','S'],
				'W':['A','T','W'],
				'C':['C'],'G':['G'],'T':['T'],'A':['A']}
	mc_class = list(mc_type) # copy
	if 'C' in mc_type:
		mc_class.extend(['CGN', 'CHG', 'CHH','CNN'])
	elif 'CG' in mc_type:
		mc_class.extend(['CGN'])
	mc_class_final = []
	for motif in mc_class:
		mc_class_final.extend([''.join(i) for i in product(*[iub_dict[nuc] for nuc in motif])])
	return(set(mc_class_final))

#Read allc file and convert to bedfile
def allc2bed(allc,return_bed=True):
	#get first line
	if allc.endswith('gz'):
		header = gzopen(allc).readline().rstrip()
	else:
		header = open(allc).readline().rstrip()
	#check if first line contains header
	if header == 'chr\tpos\tstrand\tmc_class\tmc_count\ttotal\tmethylated':
		#read in allc file to pandas dataframe
		a = pd.read_csv(allc,
			dtype={'chr':str,'pos':int,'strand':str,'mc_class':str,'mc_count':int,'total':int,'methylated':int},
			sep="\t")
	else:
		#read in allc file to pandas dataframe, add header if does not have one
		a = pd.read_csv(allc,names=['chr','pos','strand','mc_class','mc_count','total','methylated'],
			dtype={'chr':str,'pos':int,'strand':str,'mc_class':str,'mc_count':int,'total':int,'methylated':int},
			sep="\t")
	#if return_bed = True, convert to bedfile
	if return_bed is True:
		#add new columns
		a['pos2'] = a['pos']
		a['name'] = a.index
		a['score'] = '.'
		#reorder columns
		a = a[['chr','pos','pos2','name','score','strand','mc_class','mc_count','total','methylated']]
		#create pybedtools object
		a = pbt.BedTool.from_dataframe(a)
	#return a
	return a

#Collect mC data for a context
#If site_cutoff_only is set to "True", then the cutoff will only apply to tallying number of 
#sites & methylated reads called by methylpy
def get_mC_data(a,mc_type='C',cutoff=0,site_cutoff_only=False):
	#expand nucleotide list for a given context
	b = expand_nucleotide_code(mc_type=[mc_type])
	d1 = d2 = d3 = d4 = 0
	#iterate over each line
	for c in a.itertuples():
		#check if line is correct context
		if c[4] in b:
			#check if meets minimum cutoff for read depth
			if int(c[6]) >= int(cutoff):
				#add up number of sites
				d1 = d1 + 1
				#add up number of sites called methylated by methylpy
				d2 = d2 + int(c[7])
			#if site_cutoff_only is false, then apply cutoff to methylated reads
			if site_cutoff_only is False:
				if int(c[6]) >= int(cutoff):
					#add up total reads covering a site
					d3 = d3 + int(c[6])
					#add up total methylated reads covering a site
					d4 = d4 + int(c[5])
			#if site_cutoff_only is true, then do not cutoff to methylated reads
			else:
				#add up total reads covering a site
				d3 = d3 + int(c[6])
				#add up total methylated reads covering a site
				d4 = d4 + int(c[5])
	#if no sites for that context, set to 'NA'
	if d1 == 0:
		d1 = d2 = d3 = d4 = 'NA'
	#create list
	e = [mc_type,d1,d2,d3,d4]
	#return that list
	return e

#Collect total methylation data for genome or sets of chromosomes
def total_weighted_mC(allc,output=(),mc_type=['CG','CHG','CHH'],cutoff=0,chrs=[],site_cutoff_only=False):
	#read allc file
	a = allc2bed(allc,return_bed=False)
	#filter chromosome sequences
	if chrs:
		a = a[a.chr.isin(chrs)]
	#create data frame
	columns=['Context','Total_sites','Methylated_sites','Total_reads','Methylated_reads','Weighted_mC']
	b = pd.DataFrame(columns=columns)
	#iterate over each mC type and run get_mC_data
	for c in mc_type:
		d = get_mC_data(a,mc_type=c,cutoff=cutoff,site_cutoff_only=site_cutoff_only)
		#calculate weighted methylation
		d += [(float64(d[4])/float64(d[3]))]
		#add to data frame
		b = b.append(pd.DataFrame([d],columns=columns), ignore_index=True)
	#output results
	if output:
		b.to_csv(output, sep='\t', index=False)
	else:
		return b

#for calculating methylation levels in windows across the genome
def genome_window_methylation(allc,genome_file,output=(),mc_type=['CG','CHG','CHH'],
		window_size=100000,stepsize=50000,cutoff=0,chrs=[],site_cutoff_only=False):
	#read in allc file
	print("Reading allc file")
	a = allc2bed(allc)
	#create output data frame
	print("Creating output dataframe")
	c = ['Chr','Window']
	columns=['Total_sites','Methylated_sites','Total_reads','Methylated_reads','Weighted_mC']
	for d in mc_type:
		for e in columns:
			c = c + [d + '_' + e]
	b = pd.DataFrame(columns=c)
	#make windows
	print("Making genome windows")
	w_bed = pbt.bedtool.BedTool.window_maker(pbt.BedTool(genome_file),g=genome_file,
		w=window_size,s=stepsize,i='srcwinnum').filter(chr_filter,chrs)
	#intersect bedfiles with pybedtools
	print("Mapping DNA methylation data to windows")
	mapping = pbt.bedtool.BedTool.intersect(a,w_bed,wa=True,wb=True)
	del(w_bed,a)
	#convert to pandas dataframe
	print("Converting to pandas dataframe")
	m = pd.read_csv(mapping.fn,header=None,usecols=[13,6,7,8,9],sep="\t")
	del(mapping)
	#split srcwinnum
	print("Formatting names")
	f = m[13].str.rsplit('_', n = 1, expand = True)
	#make new columns from srcwinnum
	m['Chr'] = f[0]
	m['Window'] = f[1]
	del(f)
	#reorder data frame
	m = m[['Chr','Window',13,6,7,8,9]]
	#iterate over each chromosome
	print("Calculating window methylation data")
	for g in chrs:
		#get windows for that specific chromosome
		windows = list(m[m['Chr'].isin([str(g)])]['Window'].drop_duplicates())
		#iterate over each window in each chromosome
		for h in windows:
			#filter for rows matching specific chr & window number
			i = m[m['Chr'].isin([str(g)]) & m['Window'].isin([str(h)])]
			#make list for methylation data
			j = [g,h]
			#iterate over each mC type and run get_mC_data
			for k in mc_type:
				l = get_mC_data(i,mc_type=k,cutoff=cutoff,site_cutoff_only=site_cutoff_only)
				#delete first line of list
				del(l[0])
				#Calculate weighted methylation and add this to list of data for other mc_types
				j += l + [(float64(l[3])/float64(l[2]))]
			#append the results for that window to the dataframe
			b = b.append(pd.DataFrame([j],columns=c),ignore_index=True)
	#output results
	print("Outputting results")
	if output:
		b.to_csv(output, sep='\t', index=False)
	else:
		return b

#create filtered allc of sites mapping to annotations
def allc_annotation_filter(allc,annotations,genome_file,output=(),updown_stream=2000,
		first_feature=(),second_feature=(),chrs=[]):
	#read in annotations and filter by first feature, typically a something like 'gene'
	#this is used solely to accurately create flanking regions
	print("Reading annotations")
	bed = pbt.BedTool(annotations).filter(feature_filter,first_feature).filter(chr_filter,chrs)
	#create bedfile of flanking regions (if specified)
	print("Getting flanking regions")
	flank_bed = pbt.bedtool.BedTool.flank(bed,g=genome_file,l=updown_stream,
		r=updown_stream,s=True).saveas('f_bed.tmp')
	#correct sites where the start is greater than the end
	command="awk -v FS='\t' -v OFS='\t' '{$5=($4>$5&&$4==1?$4:$5)}; \
			{$4=($4>$5&&$4!=1?$5:$4)} 1' f_bed.tmp > tmp; mv tmp f_bed.tmp"
	call(command, shell=True)
	#read in annotations and filter by second annotation, typically a something like coding sequences 'CDS'
	#this is the annotation used to first filter the data
	print("Filtering second feature annotations")
	cds_bed = pbt.BedTool(annotations).filter(feature_filter,second_feature).filter(chr_filter,
		chrs).saveas('c_bed.tmp')
	#combine flanking regions with second feature
	print("Combining second feature and flanking region bed files")
	combined_bed = cds_bed.cat(flank_bed, postmerge=False)
	#read in allc file and map to annotations
	print("Reading allc file")
	a = allc2bed(allc)
	print("Mapping sites to annotations")
	mapping = pbt.bedtool.BedTool.intersect(a,combined_bed,wa=True)
	print("Converting mapped sites to table")
	m = pd.read_csv(mapping.fn, header=None, usecols = [0,1,5,6,7,8,9],sep="\t")
	#create new filtered allc file of sites mapping to regions
	print("Reformat data and drop duplicate sites")
	m.columns = ['chr','pos','strand','mc_class','mc_count','total','methylated']
	m = m.drop_duplicates().sort_values(['chr','pos'],ascending=[True,True])
	#output results	#output results
	if output:
		print("Outputing results")
		m.to_csv(output, sep='\t', index=False)
	else:
		return m
	#remove temporary files created
	print("Removing temporary files")
	tmp=['f_bed.tmp','c_bed.tmp']
	for b in tmp:
		remove(b)

#output methylation data for making metaplots of features (e.g. genes, CDS, transposons), 
#will not filter out data from introns, etc...use gene_metaplot for that
def metaplot(allc,annotations,genome_file,output=(),mc_type=['CG','CHG','CHH'],
		window_number=60,updown_stream=2000,feature=(),cutoff=0,chrs=[],site_cutoff_only=False):
	#read in allc file
	print("Reading allc file")
	a = allc2bed(allc)
	#create output data frame
	print("Create output dataframe")
	c = ['Window']
	columns=['Weighted_mC']
	for d in mc_type:
		for e in columns:
			c = c + [d + '_' + e]
	b = pd.DataFrame(columns=c)
	#read annotation file and filter on feature
	print("Reading annotations")
	f_bed = pbt.BedTool(annotations).filter(feature_filter,feature).filter(chr_filter,chrs).saveas('f_bed.tmp')
	#if updown_stream set to 0, only region to look at is the feature itself, e.g. f_bed
	if updown_stream == 0:
		print("No flanking regions specified. Analyze annotations only")
		regions=[f_bed]
	#if updown_stream specified, create bed files for upstream regions (u_bed) and down stream regions (d_bed)
	else:
		print("Get flanking regions")
		#upstream regions
		u_bed = pbt.bedtool.BedTool.flank(f_bed,g=genome_file,l=updown_stream,r=0,s=True).saveas('u_bed.tmp')
		#correct sites where the start is greater than the end
		command="awk -v FS='\t' -v OFS='\t' '{$5=($4>$5&&$4==1?$4:$5)}; \
				{$4=($4>$5&&$4!=1?$5:$4)} 1' u_bed.tmp > tmp; mv tmp u_bed.tmp"
		call(command, shell=True)
		#downstream regions
		d_bed = pbt.bedtool.BedTool.flank(f_bed,g=genome_file,l=0,r=updown_stream,s=True).saveas('d_bed.tmp')
		#correct sites where the start is greater than the end
		command="awk -v FS='\t' -v OFS='\t' '{$5=($4>$5&&$4==1?$4:$5)}; {$4=($4>$5&&$4!=1?$5:$4)} 1' \
				d_bed.tmp > tmp; mv tmp d_bed.tmp"
		call(command, shell=True)
		regions=[u_bed,f_bed,d_bed]
	#set window number to 1
	window = 1
	#iterate over each region and collect methylation data
	print("Mapping sites to annotations and collecting methylation data")
	for f in regions:
		#filter bed files based on strand
		p_bed = f.filter(strand_filter,strand='+').saveas('p_bed.tmp')
		n_bed = f.filter(strand_filter,strand='-').saveas('n_bed.tmp')
		#make windows
		pw_bed = pbt.bedtool.BedTool.window_maker(p_bed,b=p_bed,n=int(window_number/3),i='srcwinnum')
		nw_bed = pbt.bedtool.BedTool.window_maker(n_bed,b=n_bed,n=int(window_number/3),i='srcwinnum',
			reverse=True)
		#combine windows
		w_bed = pw_bed.cat(nw_bed, postmerge=False)
		#intersect bedfiles with pybedtools
		mapping = pbt.bedtool.BedTool.intersect(a,w_bed,wa=True,wb=True)
		del(w_bed,pw_bed,nw_bed)
		#convert to pandas dataframe
		m = pd.read_csv(mapping.fn,header=None,usecols=[13,6,7,8,9],sep="\t")
		del(mapping)
		#split srcwinnum
		g = m[13].str.rsplit('_', n = 1, expand = True)
		#make new columns from srcwinnum
		m['Name'] = g[0]
		m['Window'] = g[1]
		del(g)
		#reorder data frame
		m = m[['Name','Window',13,6,7,8,9]]
		#iterate list of window numbers
		for h in list(range(1,int(window_number/3)+1)):
			#filter for rows matching specific window number
			i = m[m['Window'].isin([str(h)])]
			#make list for methylation data
			j = [window]
			#iterate over each mC type and run get_mC_data
			for k in mc_type:
				#check if i is empty
				if i.empty:
					#if empty, add 0 column
					j += ['NA','NA','NA','NA','NA']
				#else if not empty
				else:
					l = get_mC_data(i,mc_type=k,cutoff=cutoff,site_cutoff_only=site_cutoff_only)
					#Check if there are any sites
					if l[1] == 'NA':
						#If no sites, output 'NA'
						j += ['NA','NA','NA','NA','NA']
					else:
						#Calculate weighted methylation and add this to list of data for other mc_types
						j += [l[1],l[2],l[3],l[4],(float64(l[4])/float64(l[3]))]
			#append the results for that window to the dataframe
			b = b.append(pd.DataFrame([j],columns=c),ignore_index=True)
			#count windows
			window += 1
	#output results
	if output:
		print("Outputing results")
		b.to_csv(output, sep='\t', index=False)
	else:
		return b
	#remove temporary files created
	print("Removing temporary files")
	tmp=['f_bed.tmp','u_bed.tmp','d_bed.tmp','p_bed.tmp','n_bed.tmp']
	for n in tmp:
		remove(n)

#output methylation data for making metaplots of features (e.g. genes, CDS, transposons), 
#for use when need to first filter data from another feature, such as intron sequences
def gene_metaplot(allc,annotations,genome_file,output=(),mc_type=['CG','CHG','CHH'],
		window_number=60,updown_stream=2000,cutoff=0,first_feature=(),second_feature=(),
		chrs=[],filtered_data_output='annotation_filtered_allc.tmp',remove_tmp=True):
	#prefilter allc file based on annotations
	print("Filtering data for sites in annotations")
	allc_annotation_filter(allc,annotations,genome_file,output=filtered_data_output,
		updown_stream=updown_stream,first_feature=first_feature,second_feature=second_feature,chrs=chrs)
	#collect methylation data
	print("Collecting metaplot data")
	metaplot(filtered_data_output,annotations,genome_file,output=output,mc_type=mc_type,
		window_number=window_number,updown_stream=updown_stream,feature=first_feature,cutoff=0,chrs=chrs)
	#remove annotation filtered allc file, if set to false, this will be kept
	if remove_tmp:
		remove(filtered_data_output)

#Calculate methylation levels for features
def feature_methylation(allc,annotations,genome_file,output=(),mc_type=['CG','CHG','CHH'],
		updown_stream=0,feature=(),cutoff=0,chrs=[],site_cutoff_only=False):
	#read in allc file
	a = allc2bed(allc)
	#create output data frame
	c = ['Feature']
	columns=['Total_C','Methylated_C','Total_Reads','Methylated_Reads','Weighted_mC']
	for d in mc_type:
		for e in columns:
			c = c + [d + '_' + e]
	b = pd.DataFrame(columns=c)
	#read annotation file and filter on feature
	f_bed = pbt.BedTool(annotations).filter(feature_filter,feature).filter(chr_filter,chrs).saveas('f_bed.tmp')
	#if updown_stream set to 0, only region to look at is the feature itself, e.g. f_bed
	if updown_stream == 0:
		regions=['f_bed.tmp']
	#if updown_stream specified, create bed files for upstream regions (u_bed) and down stream regions (d_bed)
	else:
		#upstream regions
		u_bed = pbt.bedtool.BedTool.flank(f_bed,g=genome_file,l=updown_stream,r=0,s=True).saveas('u_bed.tmp')
		#correct sites where the start is greater than the end
		command="awk -v FS='\t' -v OFS='\t' '{$5=($4>$5&&$4==1?$4:$5)}; \
				{$4=($4>$5&&$4!=1?$5:$4)} 1' u_bed.tmp > tmp; mv tmp u_bed.tmp"
		call(command, shell=True)
		#downstream regions
		d_bed = pbt.bedtool.BedTool.flank(f_bed,g=genome_file,l=0,r=updown_stream,s=True).saveas('d_bed.tmp')
		#correct sites where the start is greater than the end
		command="awk -v FS='\t' -v OFS='\t' '{$5=($4>$5&&$4==1?$4:$5)}; \
				{$4=($4>$5&&$4!=1?$5:$4)} 1' d_bed.tmp > tmp; mv tmp d_bed.tmp"
		call(command, shell=True)
		regions=['u_bed.tmp','d_bed.tmp']
	#iterate over each region and collect methylation data
	for f in regions:
		#intersect bedfiles with pybedtools
		mapping = pbt.bedtool.BedTool.intersect(a,f,wa=True,wb=True)
		del(a)
		#convert to pandas dataframe
		m = pd.read_csv(mapping.fn,header=None,usecols=[18,6,7,8,9],sep="\t")
		del(mapping)
		#split column 18
		#g = m[18].str.split(';', n = 1, expand = True)
		#g = g[0].str.split('=', n = 1, expand = True)
		#make new columns from column 18
		m['Name'] = m[18].str.split(';', n = 1, expand = True)[0].str.split('=', n = 1, expand = True)[1]
		m['Window'] = '.'
		#reorder m
		m = m[['Name','Window',18,6,7,8,9]]
		#get gene names
		g = pd.read_csv(f_bed.fn,header=None,usecols=[2,8],sep="\t")
		g = g[g[2]=="gene"][8].str.split(';', n = 1, expand = True)[0].str.split('=', n = 1, expand = True)[1]
		#iterate list of gene names
		for h in list(g):
			#filter for rows matching specific gene
			i = m[m['Name'].isin([str(h)])]
			#make list for methylation data
			j = [h]
			#iterate over each mC type and run get_mC_data
			for k in mc_type:
				#check if i is empty
				if i.empty:
					#if empty, add 0 column
					j += ['NA','NA','NA','NA','NA']
				#else if not empty
				else:
					l = get_mC_data(i,mc_type=k,cutoff=cutoff,site_cutoff_only=site_cutoff_only)
					#Check if there are any sites
					if l[1] == 'NA':
						#If no sites, output 'NA'
						j += ['NA','NA','NA','NA','NA']
					else:
						#Calculate weighted methylation and add this to list of data for other mc_types
						j += [(float64(l[4])/float64(l[3]))]
			#append the results for that window to the dataframe
			b = b.append(pd.DataFrame([j],columns=c),ignore_index=True)
	#output results
	if output:
		b.to_csv(output, sep='\t', index=False)
	else:
		return b
	#remove temporary files created
	for n in regions:
		remove(n)

#filter DMR output from methylpy based on minimum cutoffs
def filter_dmr(dmr_file,output=(),min_dms=5,min_mC_diff=0.1):
	a = pd.read_csv(dmr_file,sep='\t')
	list = []
	for b in a.itertuples():
		if b[4] >= min_dms:
			if max(b[7:]) - min(b[7:]) >= min_mC_diff:
				list.append(b[0])
	a = a.ix[list]
	if output:
		a.to_csv(output, sep='\t', index=False)
	else:
		return a

#Perform FDR correction on a column in a pandas data frame.
#Returns the data frame with the added column with the FDR correction
def FDR(df,column,new_col):
	#read the column from the data frame, exclude rows with NA and sort in descending order
	a = df[df[column] != "NaN"][column].sort_values(ascending=False)
	#Get the number of rows in the data frame
	b = len(a)
	#create a list in reverse order for the range of 1 to the length of the data frame
	c = list(reversed(range(1,b+1)))
	#Perform FDR, I can't remember why this works, but it does
	d = minimum.accumulate([b/x*y for x,y in zip(c,a)])
	#If a value is great than 1, change it to 1, else keep it the same
	e = [x if x < 1.0 else 1.0 for x in d]
	#create a data frame out of the new FDR column
	f = pd.DataFrame(e,columns=[new_col])
	#match the index for the new column to the original data frame
	f.index=a.index
	#Ad the new column to the data frame
	df = pd.concat([df,f],axis=1)
	#Return the new data frame
	return df

#Perform binomial test on number of methylated sites in a gene
def gene_binom_test(df,output=(),mc_type=['CG','CHG','CHH'],baseline={},min_sites=20,calc_baseline=True):
	#read the table
	a = pd.read_csv(df,sep='\t')
	#iterate through each listed methylation context
	for b in mc_type:
		#Check to calculate the baseline methylation
		if calc_baseline:
			#If true, then calculate the baseline for that methylation context
			BL = pd.DataFrame.sum(a[b+'_Methylated_C'])/pd.DataFrame.sum(a[b+'_Total_C'])
			#print the baseline to standard output
			print(b+' baseline value is: '+str(BL))
		#If not calculating baseline, check to see if one is provided
		elif not baseline[b]:
			#give warning and quit
			print("Baseline not given")
			return
		#If a baseline is supplied
		else:
			#Get the supplied baseline 
			BL = baseline[b]
		#Perform the binomial test for each row and and return the p-value
		a[b+'_pvalue'] = stats.binom.sf(a[b+'_Methylated_C']-1,a[b+'_Total_C'],BL)
		#If not sufficient data to perform the test, return NaN
		a[b+'_pvalue'] = [x if y >= min_sites else "NaN" for x,y in zip(a[b+'_pvalue'],a[b+'_Total_C'])]
		#Perform FDR correction
		a = FDR(a,column=b+'_pvalue',new_col=b+'_qvalue')
	#If output file specified, write to output
	if output:
		a.to_csv(output, sep='\t', index=False)
	#else return the data frame
	else:
		return a

#Subscript for classifying genes based on binomial test results
#min_sites is minimum number of sites with read coverage that are needed to consider a gene for 
#classification. qvalue is the adjusted p-value used as a cutoff for gene classificaiton
#For unmethylated genes, you can simply set a cutoff "uM_cutoff" on number of methylated reads permitted to 
#call a gene as unmethylated. The default value is 0, but this may be too harsh. However, if this value is
#set too high, then you run the risk of forcing many methylated genes as unmethylated. I recommend setting 
#this in the range of 0-2. As an alternative to classify additional genes as unmethylated, while retaining
#some rigor, you can also use the weighted_mC level for the gene as additional evidence in addition to the
#hard read cutoff. This is activated by setting a max value for the uM_weighted_mC_cutoff. By default this 
#is set to False and turned off. If you use it, it is used in addition to the uM_cutoff. In this case I would 
#recommend setting the uM_cutoff more conservatively to 0-1 and then putting a threshhold of ~0-0.02. This 
#setting works by calculating the total weighted methylation for all contexts (CG + CHG + CHH) and then
#classifying any gene with a total weighted methylation below this cutoff as unmethylated. So if the cutoff
#is set at 0.01, that means any gene with a weighted methylation above 0.01 (1%) will not be considered 
#unmethylated. The 'use_CH' argument determine whether or not to use CH in classifying genes, while, 
#'use_subCH' determines whether to use CHG & CHH separately for classification. Both options can be set 
#to True at the same time.
def classify(row,min_sites=20,qvalue=0.05,uM_cutoff=0,uM_weighted_mC_cutoff=False,use_subCH=True,use_CH=False):
	mC_sites=(row['CG_Methylated_C']+row['CHG_Methylated_C']+row['CHH_Methylated_C'])
	total_sites=(row['CG_Total_C']+row['CHG_Total_C']+row['CHH_Total_C'])
	mC_reads=(row['CG_Methylated_Reads']+row['CHG_Methylated_Reads']+row['CHH_Methylated_Reads'])
	total_reads=(row['CG_Total_Reads']+row['CHG_Total_Reads']+row['CHH_Total_Reads'])
	#Classify genes based on qvalue
	if use_CH and row['CH_qvalue'] <= qvalue and row['CH_Total_C'] >= min_sites:
		return 'TE-like'
	elif use_subCH and row['CHH_qvalue'] <= qvalue and row['CHH_Total_C'] >= min_sites:
		return 'TE-like'
	elif use_subCH and row['CHG_qvalue'] <= qvalue and row['CHG_Total_C'] >= min_sites:
		return 'TE-like'
	elif row['CG_qvalue'] <= qvalue and row['CG_Total_C'] >= min_sites:
		return 'gbM'
	#Check if number of methylated sites is at or below the cutoff to call unmethylated genes
	elif mC_sites <= uM_cutoff and total_sites >= min_sites:
		return 'Unmethylated'
	#Setting a hard cutoff on unmethylated sites may be too harsh
	elif (mC_reads/total_reads) < uM_weighted_mC_cutoff and total_sites >= min_sites:
			return 'Unmethylated'
	elif total_sites >= min_sites:
		return 'Unclassified'
	else:
		return 'NA'

#Script for calling "classify" function and applying to table
#See "classify" function for details of arguments
def classify_genes(df,output=(),min_sites=20,qvalue=0.05,uM_cutoff=0,uM_weighted_mC_cutoff=False,
		use_subCH=True,use_CH=False):
	a = pd.read_table(df,sep="\t")
	a['Classification'] = a.apply(classify,axis=1,min_sites=min_sites,qvalue=qvalue,uM_cutoff=uM_cutoff,
					uM_weighted_mC_cutoff=uM_weighted_mC_cutoff,use_subCH=use_subCH,use_CH=use_CH)
	if output:
		a.to_csv(output, sep='\t', index=False)
	else:
		return a

#Prepares allc files as input for HOME DMR calling
def allc2HOME(allc,output=()):
	#get first line
	if allc.endswith('gz'):
		#assign first line as "header"
		header = gzopen(allc).readline().rstrip()
		#open allc in text mode
		a = gzopen(allc,'rt')
	else:
		#assign first line as "header"
		header = open(allc).readline().rstrip()
		#open allc
		a = open(allc)
	#open output file
	b = open(output, 'w')
	#iterate over allc file
	with a:
		#check if there is a header line, if so, skip that line
		if header == 'chr\tpos\tstrand\tmc_class\tmc_count\ttotal\tmethylated':
			next(a)
		#iterate over each line
		for c in a:
			#split the line by tabs
			d = c.split('\t')
			#check if column 4 is CG, CHG, CHH, or CN, set variable e to that
			if d[3].startswith("CG"):
				e = 'CG'
			elif d[3].endswith("G"):
				e = 'CHG'
			elif d[3].startswith("CN") or c[3].endswith("N"):
				e = 'CNN'
			else:
				e = 'CHH'
			#output needed lines
			b.write(d[0]+'\t'+d[1]+'\t'+d[2]+'\t'+e+'\t'+d[4]+'\t'+d[5]+'\n')
	#cloese the files
	b.close()
	a.close()

#Returns GC content and GC for each codon for a fasta file
def gc123(fasta,output=()):
	#create header
	a=[['Transcript','GC','GC1','GC2','GC3']]
	#parse and read over fasta file
	for b in SeqIO.parse(open(fasta, "r"), "fasta"):
		#use Biopython GC123 to calculate GC content and codon GC content
		c = GC123(b.seq)
		#for each sequence, add these values to the array
		a = a + [[b.id,str(c[0]),str(c[1]),str(c[2]),str(c[3])]]
	#output results
	if output:
		with open(output, 'w') as d:
			d.writelines('\t'.join(e) + '\n' for e in a)
	else:
		return(a)

#Function for creating a bed file of motifs found in fasta
#Example use: motif2bed("ATCG","test.fa",reverse_strand=True,output="test.bed")
def motif2bed(motif,fasta,reverse_strand=True,output=()):
	motif_bed = []
	m = motifs.create([Seq(motif)])
	#iterate over fasta file and search for motifs
	for record in SeqIO.parse(open(fasta, "r"), "fasta"):
		chrom = record.id
		#Check forward strand
		for pos, seq in m.instances.search(str(record.seq)):
			start_pos = pos
			stop_pos = pos + len(seq)
			motif_bed = motif_bed + [[str(chrom),str(start_pos),str(stop_pos),str(motif),".","+"]]
		#Check reverse strand
		if reverse_strand:
			for pos, seq in m.reverse_complement().instances.search(str(record.seq)):
				start_pos = pos
				stop_pos = pos + len(seq)
				motif_bed = motif_bed + [[str(chrom),str(start_pos),str(stop_pos),str(motif),".","-"]]
	if output:
		with open(output, 'w') as file:
			file.writelines('\t'.join(i) + '\n' for i in motif_bed)
	else:
		return(motif_bed)

#Function for adding columns with combined CHG & CHH methylation data
def add_CH(df,output):
	a = pd.read_csv(df,sep='\t')
	a['CH_Total_C'] = a['CHG_Total_C'] + a['CHH_Total_C']
	a['CH_Methylated_C'] = a['CHG_Methylated_C'] + a['CHH_Methylated_C']
	a['CH_Total_Reads'] = a['CHG_Total_Reads'] + a['CHH_Total_Reads']
	a['CH_Methylated_Reads'] = a['CHG_Methylated_Reads'] + a['CHH_Methylated_Reads']
	a['CH_Weighted_mC'] = a['CHG_Methylated_Reads']/a['CHH_Total_Reads']
	if output:
		a.to_csv(output, sep='\t', index=False)
	else:
		return(a)



