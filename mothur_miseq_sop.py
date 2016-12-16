import os
import subprocess
import argparse
from argparse import RawTextHelpFormatter
__author__ = 'Andrew Thomas'
 
parser = argparse.ArgumentParser(description='This is a script to analyze V4 16S amplicon data using the mothur MiSeq SOP by Andrew Thomas.\nExample usage:\npython mothur_miseq_sop.py -fq /work/data/ -wd /work/ -s samples.txt -oligos Oligos.txt\nThe sample-file mapping should be as follows:\nOct_12\tZLAM0113IQMSIPE01R01R1.fastq\tZLAM0113IQMSIPE01R01R2.fastq\n', formatter_class=RawTextHelpFormatter)
parser.add_argument('-fq','--fastq_dir', help="Absolute fastq input file directory",required=True)
parser.add_argument('-wd','--working_dir',help='Absolute working directory', required=True)
parser.add_argument('-s','--sample_mapping',help='Mapping of sample names to file names', required=True)
parser.add_argument('-oligos','--oligos_file',help='File with the 16s primer sequences', required=True)
args = parser.parse_args()

sample_files = {}
sample_mapping = args.sample_mapping
sample_mapping1 = open(args.sample_mapping, "r")

for line in sample_mapping1.readlines():
                line = line.rstrip()
                line = line.split("\t")
                sample_files[line[0]] = [line[1], line[2]]
                sample_mapping1.close()

# Edit directories
fastq_dir = args.fastq_dir
working_directory = args.working_dir
oligos = args.oligos_file

print "Making Contigs..."

# Make contigs
subprocess.call("mkdir Contigs", shell = True)
os.chdir("Contigs/")
for key in sample_files:
		subprocess.call("mothur '# make.contigs(ffastq=" + fastq_dir + sample_files[key][0] + ", rfastq=" + fastq_dir + sample_files[key][1] + ", oligos=" + working_directory + oligos + ", pdiffs=2, processors=8)'", shell = True)

subprocess.call("mv " + fastq_dir + "*.trim.contigs.fasta " + working_directory + "Contigs/", shell = True)

## Remove raw data
subprocess.call("rm " + fastq_dir + "*scrap.contigs*", shell = True)
subprocess.call("rm " + fastq_dir + "*contigs.report", shell = True)
subprocess.call("rm " + fastq_dir + "*trim.contigs.qual", shell = True)

## Create run maps
os.chdir(working_directory)

run_dic = {}
for fasta in os.listdir(str(working_directory + "Contigs/")):
        if fasta.endswith("trim.contigs.fasta"):
                contig = fasta.split(".trim")[0]
                sample = subprocess.check_output("grep " + contig + " " + str(sample_mapping) + " | cut -f1", shell = True).rstrip()
                run_dic[sample] = fasta

## Make groups file
fastas = ""
for fasta in run_dic.values():
        fastas += "Contigs/" + str(fasta) + "-"

fastas = fastas[:-1]

samples = ""
for sample in run_dic.keys():
        samples += str(sample) + "-"

samples = samples[:-1]

print fastas
print samples

## MiSEQ SOP 
print "MiSeq SOP started"
subprocess.call("mothur '# make.group(fasta=" + fastas + ", groups=" + samples + ")'", shell = True)

subprocess.call("mothur '# merge.files(input=" + fastas + ", output=merged.fasta)'", shell = True)

subprocess.call("mothur '# summary.seqs(fasta=merged.fasta)'", shell = True)

subprocess.call("mv Contigs/mergegroups " + working_directory, shell = True)

subprocess.call("mothur '# screen.seqs(fasta=merged.fasta, group=mergegroups, maxambig=0, maxlength=275)'", shell = True)

subprocess.call("mothur '# unique.seqs(fasta=merged.good.fasta)'", shell = True)

subprocess.call("mothur '# count.seqs(name=merged.good.names, group=mergegroupsgood, processors=8)'", shell = True)

subprocess.call("mothur '# summary.seqs(count=merged.good.count_table, fasta=merged.good.unique.fasta)'", shell = True)

subprocess.call("mothur '# pcr.seqs(fasta=silva.bacteria.fasta, start=11894, end=25319, keepdots=F, processors=8)'", shell = True)

subprocess.call("mothur '# system(mv silva.bacteria.pcr.fasta silva.v4.fasta)'", shell = True)

subprocess.call("mothur '# align.seqs(fasta=merged.good.unique.fasta, reference=silva.v4.fasta, processors=8)'", shell = True)

subprocess.call("mothur '# summary.seqs(fasta=merged.good.unique.align, count=merged.good.count_table)'", shell = True)

subprocess.call("mothur '# screen.seqs(fasta=merged.good.unique.align, count=merged.good.count_table, start=1968, end=11550, maxhomop=8, processors=8)'", shell = True)

subprocess.call("mothur '# filter.seqs(fasta=merged.good.unique.good.align, vertical=T, processors=8, trump=.)'", shell = True)

subprocess.call("mothur '# unique.seqs(fasta=merged.good.unique.good.filter.fasta, count=merged.good.good.count_table)'", shell = True)

print "Clustering sequences..."
subprocess.call("mothur '# pre.cluster(fasta=merged.good.unique.good.filter.unique.fasta, count=merged.good.unique.good.filter.count_table, diffs=2, processors=8)'", shell = True)

subprocess.call("mothur '# chimera.uchime(fasta=merged.good.unique.good.filter.unique.precluster.fasta, count=merged.good.unique.good.filter.unique.precluster.count_table, dereplicate=t, processors=8)'", shell = True)

subprocess.call("mothur '# remove.seqs(fasta=merged.good.unique.good.filter.unique.precluster.fasta, accnos=merged.good.unique.good.filter.unique.precluster.denovo.uchime.accnos)'", shell = True)

print "Classifying sequences..."
subprocess.call("mothur '# classify.seqs(fasta=merged.good.unique.good.filter.unique.precluster.pick.fasta, count=merged.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table, reference=trainset9_032012.pds.fasta, taxonomy=trainset9_032012.pds.tax, cutoff=80, processors=8)'", shell = True)

subprocess.call("mothur '# remove.lineage(fasta=merged.good.unique.good.filter.unique.precluster.pick.fasta, count=merged.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table, taxonomy=merged.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)'", shell = True)

subprocess.call("mothur '# cluster.split(fasta=merged.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=merged.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table, taxonomy=merged.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, splitmethod=classify, taxlevel=4, cutoff=0.15, processors=8)'", shell = True)

subprocess.call("mothur '# make.shared(list=merged.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.list, count=merged.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table, label=0.03)'", shell = True)

subprocess.call("mothur '# classify.otu(list=merged.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.list, count=merged.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table, taxonomy=merged.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, label=0.03)'", shell = True)

subprocess.call("mothur '# phylotype(taxonomy=merged.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy)'", shell = True)

subprocess.call("mothur '# make.shared(list=merged.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.list, count=merged.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table, label=1)'", shell = True)

subprocess.call("mothur '# classify.otu(list=merged.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.list, count=merged.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table, taxonomy=merged.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, label=1)'", shell = True)

subprocess.call("mothur '# make.biom(shared=merged.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared, constaxonomy=merged.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.taxonomy)'", shell = True)

subprocess.call("mv merged.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.biom mothur.OTUbased.table.biom", shell = True)

print "Making OTU table..."
subprocess.call("mothur '# make.biom(shared=merged.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.shared, constaxonomy=merged.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.1.cons.taxonomy)'", shell = True)

subprocess.call("mv merged.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.1.biom mothur.phylotypebased.table.biom", shell = True)

subprocess.call("biom summarize-table -i mothur.OTUbased.table.biom -o mothur.OTUbased.table.summary", shell = True)

print "Making phylogenetic tree.."
subprocess.call("mothur '# dist.seqs(fasta=merged.good.unique.good.filter.unique.precluster.pick.pick.fasta, output=lt, processors=8)'", shell = True)

subprocess.call("mothur '# clearcut(phylip=merged.good.unique.good.filter.unique.precluster.pick.pick.phylip.dist)'", shell = True)

subprocess.call("mothur '# mothur '# unifrac.weighted(tree=merged.good.unique.good.filter.unique.precluster.pick.pick.phylip.tre, count=merged.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table, subsample=T, distance=square)'", shell = True)

subprocess.call("mothur '# mothur '# unifrac.unweighted(tree=merged.good.unique.good.filter.unique.precluster.pick.pick.phylip.tre, count=merged.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table, subsample=T, distance=square)'", shell = True)


