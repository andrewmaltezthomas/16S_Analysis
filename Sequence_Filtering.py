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
scripts_directory = args.script_dir

def write_report(sample):
        total_seqs = int(subprocess.check_output("wc -l " + fastq_dir + sample_files[key][0], shell = True).rstrip().split()[0]) / 4
        remaining_seqs = int(subprocess.check_output("wc -l " + key + "-Filtered.renamed.fna", shell = True).rstrip().split()[0]) / 2
        joined_seqs = int(subprocess.check_output("wc -l " + key + "/" + key + "-join.fastq | cut -f1", shell = True).rstrip().split()[0]) / 4
        unjoined_seqs = int(subprocess.check_output("wc -l " + key + "/" + key + "-un1 | cut -f1", shell = True).rstrip().split()[0]) / 4
        fprimer_filtered = int(subprocess.check_output("cat " + key + "/" + key + "-Split-Libraries/split_library_log.txt | grep 'Num mismatches in primer exceeds' | cut -f2 -d':'", shell = True).lstrip())
        size_filtered = int(subprocess.check_output("cat " + key + "/" + key + "-Split-Libraries/split_library_log.txt | grep 'Length outside bounds' | cut -f2", shell = True).lstrip())
        ambiguous_bases = int(subprocess.check_output("cat " + key + "/" + key + "-Split-Libraries/split_library_log.txt | grep 'Num ambiguous bases exceeds' | cut -f2", shell = True).lstrip())
        quality_score = int(subprocess.check_output("cat " + key + "/" + key + "-Split-Libraries/split_library_log.txt | grep 'Number of sequences where a low quality score' | cut -f2 -d':'", shell = True).lstrip())
        homopolymer = int(subprocess.check_output("cat " + key + "/" + key + "-Split-Libraries/split_library_log.txt | grep 'Max homopolymer' | cut -f2", shell = True).lstrip())
        rprimer_filtered = int(subprocess.check_output("cat " + key + "/" + key + "-Reverse-Primer-Split/rev_primer_truncation.log | grep 'Reverse primers not found' | cut -f2 -d':'", shell = True).lstrip())
        percent_remain = float(remaining_seqs) / total_seqs * 100.0
        subprocess.call("echo -e '" + key + "\t" + str(total_seqs) + "\t" + str(joined_seqs) + "\t" + str(unjoined_seqs) + "\t" + str(fprimer_filtered) + "\t" + str(rprimer_filtered) + "\t" + str(size_filtered) + "\t" + str(ambiguous_bases) + "\t" + str(quality_score) + "\t" + str(homopolymer) + "\t" + str(remaining_seqs) + "\t" + str(percent_remain) + "' >> Sequence.Filtering.Report", shell = True)

# Make directories
subprocess.call("mkdir Pre-Processing/", shell = True)
os.chdir(working_directory + "Pre-Processing")
subprocess.call("mkdir Sequence-Filtering/", shell = True)
os.chdir(working_directory + "Pre-Processing/Sequence-Filtering")
 
# Join paired reads + filter sequences
subprocess.call("echo -e 'Parameters for filtering sequences\nForward Primer mismatches: 2\nReverse Primer mismatches: 2\nAmbiguous Bases: 2\nMinimum length: 100\nMaximum length: 1000\nMax homopolymer: 10\nMin qual score: 20\nWindow size: 50\n' > Sequence.Filtering.Report", shell = True)
subprocess.call("echo -e 'Sample\tTotal_Reads\tJoined_Reads\tUnjoined_Reads\tForward_Mismatches\tReverse_Mismatches\tSize_Filtered\tAmbiguous_Bases\tQuality_Filtered\tHomopolymer_Filtered\tRemaining_Reads\tPercent_of_total' >> Sequence.Filtering.Report", shell = True)
for key in sample_files:
	subprocess.call("mkdir " + key, shell = True)
	os.chdir(key)
	print "Joining paired-end reads for " + key
	subprocess.call("fastq-join " + fastq_dir + sample_files[key][0] + " " + fastq_dir + sample_files[key][1] + " -o " + key + "-", shell = True)
	subprocess.call("mv " + key + "-join " + key + "-join.fastq", shell = True)
	p1 = subprocess.Popen(["/usr/bin/qiime", "convert_fastaqual_fastq.py", "-f", key + "-join.fastq", "-c", "fastq_to_fastaqual", "-o", "Fasta-Qual/"])
	p1.wait()
	print "Filtering sequences for " + key
	p2 = subprocess.Popen(["/usr/bin/qiime", "split_libraries.py", "-f", "Fasta-Qual/*.fna", "-q", "Fasta-Qual/*.qual", "-l", "100", "-s", "20", "-a", "2", "-H", "10", "-w", "50", "-g", "-M", "1", "-b", "8", "-e", "1", "-o", key + "-Split-Libraries/", "-m", script_dir + "Utilities/Amplicon.map"])
	p2.wait()
	p3 = subprocess.Popen(["/usr/bin/qiime", "truncate_reverse_primer.py", "-f", key + "-Split-Libraries/seqs.fna", "-m", script_dir + "Utilities/Amplicon.map", "-z", "truncate_remove", "-M", "2", "-o", key + "-Reverse-Primer-Split/"])
	p3.wait()
	subprocess.call("mv " + key + "-Reverse-Primer-Split/seqs_rev_primer_truncated.fna " + key + "-Filtered.fna", shell = True)
	subprocess.call("python " + script_dir + "Utilities/fasta_number.py " + key + "-Filtered.fna " + key + "_ > " + key + "-Filtered.renamed.fna", shell = True)
	subprocess.call("rm " + key + "-Filtered.fna", shell = True)
	subprocess.call("mv " + key + "-Filtered.renamed.fna " + working_directory + "Pre-Processing/Sequence-Filtering/" + key + "-Filtered.renamed.fna", shell = True)
	os.chdir(str(working_directory + "Pre-Processing/Sequence-Filtering"))
	write_report(key)

subprocess.call("sed -i 's/-e //g' Sequence.Filtering.Report", shell = True)

