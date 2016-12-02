# 16S Analysis using V4-V5 MiSeq and UPARSE

#### Dependencies:
* Qiime (>= version 1.7.0)
* USEARCH (>= version 9.0)
* Fastq-join (>= version )

#### Notes:
* Assumes that Qiime executable is in "/usr/bin/qiime". If not, change this in the script. 
* Assumes that the primer pair is Forward: CCTACGGGNGGCWGCAG and Reverse: GACTACHVGGGTATCTAATCC. If not, change this in the file Utilities/Amplicon.map. 

Sequence_Filtering.py - Will join V4-V5 16S paired-end MiSeq reads using fastq-join, filter by quality, primer mismatches, homopolymers and size. A report is generated at the end showing numbers for each filter. 

Sequence_Clustering.py - Will concatenate all quality filtered reads generated by the Sequece_Filtering.py script and use them as input for *de novo* OTU clustering using UPARSE at 97% identity. Will also classify representative seed sequences using RDP classifier and a 80% confidence threshold. 

