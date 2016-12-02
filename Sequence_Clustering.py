import os
import subprocess
import argparse
from argparse import RawTextHelpFormatter
__author__ = 'Andrew Thomas'
 
parser = argparse.ArgumentParser(description='This is a sequence clustering script using UPARSE for filtered reads by Andrew Thomas.\nExample usage:\npython sequence_clustering.py -fasta_dir /work/Pre-Processing/Sequence-Filtering/ -wd /work/Pre-Processing/ -sc /work/16S_Analysis/ -p HMP', formatter_class=RawTextHelpFormatter)
parser.add_argument('-fasta_dir','--fasta_dir', help="Absolute path of filtered reads",required=True)
parser.add_argument('-wd','--working_dir',help='Absolute working directory', required=True)
parser.add_argument('-sc','--script_dir', help="Absolute file directory for this repository",required=True)
parser.add_argument('-p','--project_name',help='Project name', required=True)
args = parser.parse_args()


# Edit directories
filtered_seqs = args.fasta_dir
working_directory = args.working_dir
scripts_directory = args.script_dir

project = args.project_name

os.chdir(working_directory)
subprocess.call("mkdir Sequence-Clustering", shell = True)
subprocess.call("cat " + filtered_seqs + "*.fna > " + working_directory + "Sequence-Clustering/All.filtered.fna", shell = True)
os.chdir(working_directory + "Sequence-Clustering/")
subprocess.call("perl " + scripts_directory + "Qiime2Uparse.pl -i All.filtered.fna -o All.filtered.uparse.fna", shell = True)
p1 = subprocess.Popen(["usearch", "-derep_fulllength", "All.filtered.uparse.fna", "-fastaout", "All.filtered.uparse.derep.fna", "-sizeout"])
p1.wait()
p2 = subprocess.Popen(["usearch", "-sortbysize", "All.filtered.uparse.derep.fna", "-fastaout", "All.filtered.uparse.derep.sorted.fna", "-minsize", "2"])
p2.wait()
p3 = subprocess.Popen(["usearch", "-cluster_otus", "All.filtered.uparse.derep.sorted.fna", "-otus", "All.filtered.uparse.derep.sorted.otus1.fna"])
p3.wait()
p4 = subprocess.Popen(["usearch", "-uchime_ref", "All.filtered.uparse.derep.sorted.otus1.fna", "-db", scripts_directory + "gold.fa", "-strand", "plus", "-nonchimeras", "All.filtered.uparse.derep.sorted.otus2.fna"])
p4.wait()
subprocess.call("fasta_formatter -i All.filtered.uparse.derep.sorted.otus2.fna -o All.filtered.uparse.derep.sorted.otus2.renamed.fna", shell = True)
subprocess.call("perl " + scripts_directory + "otuName.pl -i All.filtered.uparse.derep.sorted.otus2.renamed.fna -o All.filtered.rep.set.otus.fna", shell = True)
p5 = subprocess.Popen(["usearch", "-usearch_global", "All.filtered.uparse.fna", "-db", "All.filtered.rep.set.otus.fna", "-strand", "plus", "-id", "0.97", "-uc", "All.filtered.uc"])
p5.wait()
p6 = subprocess.Popen(["/usr/bin/qiime", "assign_taxonomy.py", "-i", "All.filtered.rep.set.otus.fna", "-o", "All.filtered.rep.set.otus.Taxonomy", "-m", "rdp" ])
p6.wait()
p7 = subprocess.Popen(["/usr/bin/qiime", "align_seqs.py", "-i", "All.filtered.rep.set.otus.fna", "-o", "All.filtered.rep.set.otus.Aligned"])
p7.wait()
p8 = subprocess.Popen(["/usr/bin/qiime", "filter_alignment.py", "-i", "All.filtered.rep.set.otus.Aligned/All.filtered.rep.set.otus_aligned.fasta", "-o", "All.filtered.rep.set.otus.Filtered.Aligned", "-m", scripts_directory + "lanemask_in_1s_and_0s"])
p8.wait()
p9 = subprocess.Popen(["/usr/bin/qiime", "make_phylogeny.py", "-i", "All.filtered.rep.set.otus.Filtered.Aligned/All.filtered.rep.set.otus_aligned_pfiltered.fasta", "-o", "All.filtered.rep.set.otus.tree"])
p9.wait()
subprocess.call("python " + scripts_directory + "map2qiime.py All.filtered.uc > All.filtered.txt", shell = True)
p10 = subprocess.Popen(["/usr/bin/qiime", "make_otu_table.py", "-i", "All.filtered.txt", "-t", "All.filtered.rep.set.otus.Taxonomy/All.filtered.rep.set.otus.taxonomy.txt", "-o", project + ".otu.table.biom"])
p10.wait()
p11 = subprocess.Popen(["/usr/bin/qiime", "filter_otus_from_otu_table.py", "-i", project + ".otu.table.biom", "-s", "3", "-o", project + ".otu.table.filtered.biom"])
p11.wait()
p12 = subprocess.Popen(["/usr/bin/qiime", "sort_otu_table.py", "-i", project + ".otu.table.filtered.biom", "-o", project + ".otu.table.filtered.sorted.biom"])
p12.wait()
subprocess.call("biom summarize-table -i " + project + ".otu.table.filtered.sorted.biom -o " + project + ".otu.table.filtered.sorted.summary", shell = True)

