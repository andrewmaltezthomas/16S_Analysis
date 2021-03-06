# 16S rRNA pre-processing of Illumina MiSeq reads

## Alternative 1 - Sequence filtering and clustering using *Qiime* and *UPARSE* 

#### Dependencies (tested on):
* Qiime (version 1.7.0)
* USEARCH (version v8.1.1812)
* Fastq-join (version 1.1.2)
* Python (2.7)
* Perl (version 5.18.2)
* Fastx Toolkit (version 0.0.14)
* [Biom](http://biom-format.org/index.html)

#### Notes:
* Assumes that Qiime executable is in "/usr/bin/qiime". **If not, change this in the script.** 
* Assumes that the primer pair is Forward: CCTACGGGNGGCWGCAG and Reverse: GACTACHVGGGTATCTAATCC. 
**If not, change this in the file Utilities/Amplicon.map.**

**Sequence_Filtering.py** - Will join V3-V4 16S paired-end MiSeq reads using fastq-join, filter by quality, primer mismatches, homopolymers and size. A report is generated at the end showing numbers for each filter. 

### Alternative sequence filtering pipeline

#### Dependencies (tested on):
* Prinseq-lite (version 0.20.4)
* PEAR (Paired-end reAd mergeR - version 0.9.10)
* cutadapt (version 1.12)

### Notes:
* Assumes that pear is called as "pear-merger". **If not, change this in the script.**
* By default, reads are quality filtered using an average *phred-like* quality of 20. **If not, change this in the script.**
* Assumes that the primer pair is Forward: CCTACGGGRSGCAGCAG and Reverse: ATTAGAWACCCVHGTAGTCC (reverse complement). 
**If not, change this in the script.**

**alternative_sequence_filtering.py** - Will filter reads by average quality (default 20) using prinseq, merge paired-end reads using PEAR, trim the 16S primers using cutadapt (allowing a maximum of 2 mismatches) and filter by size (minimum 300nt). A report is generated at the end showing numbers for each filter.

**Sequence_Clustering.py** - Will concatenate all quality filtered reads generated by the **Sequence_Filtering.py** script and use them as input for *de novo* OTU clustering using UPARSE at 97% identity. Will also classify representative seed sequences using RDP classifier and a 80% confidence threshold. 

For more information on the methods, see [Antunes et al. 2016](http://www.nature.com/articles/srep38915)


## Alternative 2 - Mothur V4 MiSeq SOP script (accessed on Jan-2017)

#### Dependencies (tested on):
* mothur (version 1.39.4)

### Notes:
* Download the [Silva database](https://www.mothur.org/w/images/9/98/Silva.bacteria.zip) and the [RDP training set](https://www.mothur.org/w/images/5/59/Trainset9_032012.pds.zip). 
* Assumes that 16S primers haven't been trimmed. 
* This script runs the cluster.split heuristic as opposed to the more traditional dist.seqs. 
"In this approach, we use the taxonomic information to split the sequences into bins and then cluster within each bin. In our testing, the MCC values when splitting the datasets at the class and genus levels were within 98.0 and 93.0%, respectively, of the MCC values obtained from the entire test dataset. These decreases in MCC value resulted in the formation of as many as 4.7 and 22.5% more OTUs, respectively, than were observed from the entire dataset. The use of the cluster splitting heuristic was probably not worth the loss in clustering quality. However, as datasets become larger, it may be necessary to use the heuristic to clustering the data into OTUs. The advantage of the cluster.split approach is that it should be faster, use less memory, and can be run on multiple processors. In an ideal world we would prefer the traditional route because "Trad is rad", but we also think that kind of humor is funny.... In this command we use taxlevel=4, which corresponds to the level of Order."

For more information on the methods and specific command calls see the [SOP](https://www.mothur.org/wiki/MiSeq_SOP). 

## Alternative 3 - dada2 pipeline in snakemake

COMING SOON

#### Dependencies (tested on):
* 
