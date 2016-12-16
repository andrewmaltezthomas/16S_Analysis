# 16S rRNA pre-processing of Illumina MiSeq reads

## Sequence filtering and clustering using *Qiime* and *UPARSE* 

#### Dependencies:
* Qiime (version >= 1.7.0)
* USEARCH (version >= 9.0)
* Fastq-join (version >= 1.1.2)
* Python (version 2.7)
* Perl (version >= 5.18.2)
* Fastx Toolkit (version >= 0.0.14-1)

#### Notes:
* Assumes that Qiime executable is in "/usr/bin/qiime". **If not, change this in the script.** 
* Assumes that the primer pair is Forward: CCTACGGGNGGCWGCAG and Reverse: GACTACHVGGGTATCTAATCC. 
**If not, change this in the file Utilities/Amplicon.map.**

**Sequence_Filtering.py** - Will join V4-V5 16S paired-end MiSeq reads using fastq-join, filter by quality, primer mismatches, homopolymers and size. A report is generated at the end showing numbers for each filter. 

**Sequence_Clustering.py** - Will concatenate all quality filtered reads generated by the **Sequence_Filtering.py** script and use them as input for *de novo* OTU clustering using UPARSE at 97% identity. Will also classify representative seed sequences using RDP classifier and a 80% confidence threshold. 

For more information on the methods, see [Antunes et al. 2016](http://www.nature.com/articles/srep38915)

## Alternative sequence filtering pipeline

#### Dependencies:
* Prinseq-lite (version )
* PEAR (Paired-end reAd mergeR - version )
* cutadapt (version )

### Notes:
* Assumes that pear is called as "pear-merger". **If not, change this in the script.**
* By default, reads are quality filtered using an average phred-like quality of 20. **If not, change this in the script.**
* Assumes that the primer pair is Forward: CCTACGGGNGGCWGCAG and Reverse: GACTACHVGGGTATCTAATCC. 
**If not, change this in the script.**

**alternative_sequence_filtering.py** - Will filter reads by average quality (default 20) using prinseq, merge paired-end reads using PEAR, trim the 16S primers using cutadapt (allowing a maximum of 2 mismatches) and filter by size (minimum 300nt). A report is generated at the end showing numbers for each filter.

## Mothur V4 MiSeq SOP script

#### Dependencies:
* mothur (version >= )

### Notes:
* Download the [Silva database](https://www.mothur.org/w/images/9/98/Silva.bacteria.zip) and the [RDP training set](https://www.mothur.org/w/images/5/59/Trainset9_032012.pds.zip). 
* Assumes that 16S primers haven't been trimmed. 

For more information on the methods and specific command calls see the [SOP](https://www.mothur.org/wiki/MiSeq_SOP). 


