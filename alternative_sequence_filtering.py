import os
import subprocess
import argparse
from argparse import RawTextHelpFormatter
__author__ = 'Andrew Thomas'
 
parser = argparse.ArgumentParser(description='This is a sequence filtering script for Illumina V4-V5 paired-end reads by Andrew Thomas.\nExample usage:\npython sequence_filtering.py -fq /work/data/ -wd /work/ -sc /home/16S_Analysis/ -s samples.txt\nThe sample-file mapping should be as follows:\nOct_12\tZLAM0113IQMSIPE01R01R1.fastq\tZLAM0113IQMSIPE01R01R2.fastq\n', formatter_class=RawTextHelpFormatter)
parser.add_argument('-fq','--fastq_dir', help="Absolute fastq input file directory",required=True)
parser.add_argument('-wd','--working_dir',help='Absolute working directory', required=True)
parser.add_argument('-sc','--script_dir', help="Absolute file directory for this repository",required=True)
parser.add_argument('-s','--sample_mapping',help='Mapping of sample names to file names', required=True)
args = parser.parse_args()

sample_files = {}
sample_mapping = open(args.sample_mapping, "r")

for line in sample_mapping.readlines():
    line = line.rstrip()
    line = line.split("\t")
    sample_files[line[0]] = [line[1], line[2]]
sample_mapping.close()

# Edit directories
fastq_dir = args.fastq_dir
working_directory = args.working_dir
script_dir = args.script_dir

# Make directories
subprocess.call("mkdir Pre-Processing/", shell = True)
os.chdir(working_directory + "Pre-Processing")
subprocess.call("mkdir Sequence-Filtering/", shell = True)
os.chdir(working_directory + "Pre-Processing/Sequence-Filtering")


# Filter sequences + Join Paired-End reads
subprocess.call("echo 'sample\tTotal_Reads\tQuality_Filtered_Reads\tMerged_Reads\tPrimer_Trimmed_Reads\tPercentage_of_Total' > Sequence.Filtering.Report", shell = True)
for key in sample_files:
	subprocess.call("mkdir " + key, shell = True)
	os.chdir(key)
	total_reads = int(subprocess.check_output("wc -l " + fastq_dir + sample_files[key][0], shell = True).rstrip().split()[0]) / 4
	
	print "Filtering reads for " + key
	subprocess.call("mkdir Filtered_Reads/", shell = True)
	p1 = subprocess.Popen("prinseq-lite -fastq " + fastq_dir + sample_files[key][0] + " -fastq2 " +  fastq_dir + sample_files[key][1] + " -out_good Filtered_Reads/" + key + ".quality.20.Good -out_bad null -min_qual_mean 20 -ns_max_n 1", shell = True)
	p1.wait()
	subprocess.call("rm Filtered_Reads/*_singletons.fastq", shell = True)
	p2 = subprocess.Popen("python " + script_dir + "Utilities/combine_pairs.py Filtered_Reads/" + key + ".quality.20.Good_1.fastq Filtered_Reads/" + key + ".quality.20.Good_2.fastq", shell = True)
	p2.wait()
	print
	quality_filtered_reads = int(subprocess.check_output("wc -l Filtered_Reads/" + key + ".quality.20.Good_1.fastq", shell = True).rstrip().split()[0]) / 4
	subprocess.call("rm Filtered_Reads/" + key + ".quality.20.Good_1.fastq Filtered_Reads/" + key + ".quality.20.Good_2.fastq", shell = True)
	subprocess.call("rm Filtered_Reads/*_singles.fastq", shell = True)
	subprocess.call("mv Filtered_Reads/" + key + ".quality.20.Good_1.fastq_pairs_R1.fastq Filtered_Reads/" + key + ".quality.20.filtered.R1.fastq", shell = True)
	subprocess.call("mv Filtered_Reads/" + key + ".quality.20.Good_2.fastq_pairs_R2.fastq Filtered_Reads/" + key + ".quality.20.filtered.R2.fastq", shell = True)
	
	print "Joining paired-end reads for " + key
	subprocess.call("mkdir PEAR/", shell = True)
	p4 = subprocess.Popen("pear-merger -f Filtered_Reads/" + key + ".quality.20.filtered.R1.fastq -r Filtered_Reads/" + key + ".quality.20.filtered.R2.fastq -o PEAR/" + key, shell = True)
	p4.wait()
	merged_reads = int(subprocess.check_output("wc -l PEAR/" + key + ".assembled.fastq", shell = True).rstrip().split()[0]) / 4

	print "Trimming primers for " + key
	p3 = subprocess.Popen("cutadapt -g ^CCTACGGGRSGCAGCAG --discard-untrimmed -e 0.12 PEAR/" + key + ".assembled.fastq > " + key + ".Ftrim.fastq 2> " + key + ".Ftrim.report.txt", shell = True)
	p3.wait()
	p3 = subprocess.Popen("cutadapt -a ATTAGAWACCCDBGTAGTCC$ --minimum-length 300 --discard-untrimmed -e 0.1 " + key + ".Ftrim.fastq > " + key + ".Primer.Trim.fastq 2> " + key + ".Rtrim.report.txt", shell = True)
	p3.wait()
	subprocess.call("rm *.Ftrim.fastq", shell = True)
	primer_trimmed_reads = int(subprocess.check_output("wc -l " + key + ".Primer.Trim.fastq", shell = True).rstrip().split()[0]) / 4
	remaining_reads = (primer_trimmed_reads / float(total_reads)) * 100.0
	subprocess.call("fastq_to_fasta -i " +  key + ".Primer.Trim.fastq -o " + key + ".Primer.Trim.fna", shell = True) 
	subprocess.call("python " + script_dir + "Utilities/fasta_number.py " + key + ".Primer.Trim.fna " + key + "_ > " + key + "-Filtered.renamed.fna", shell = True)
	subprocess.call("mv " + key + "-Filtered.renamed.fna " + working_directory + "Pre-Processing/Sequence-Filtering/", shell = True)
	os.chdir(str(working_directory + "Pre-Processing/Sequence-Filtering"))
	subprocess.call("echo '" + key + "\t" + str(total_reads) + "\t" + str(quality_filtered_reads) + "\t" + str(merged_reads) + "\t" + str(primer_trimmed_reads) + "\t" + str(remaining_reads) + "' >> Sequence.Filtering.Report", shell = True)

