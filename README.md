# RNA Seq Preprocessing Using Snakemake
 

# Table of contents
1. [Introduction](#introduction)
2. [Data](#data)
3. [Pipeline](#pipeline)
4. [Usage](#usage)
    1. [Input](#input)
    2. [Output](#output)
5. [Results](#results)
6. [References](#references)


## RNA Sequencing Preprocessing <a name="introduction"></a>

RNA Sequencing quantifies the amount of RNA in a biological sample at a given moment. Thus, between 2 conditions (eg. before/after treatment or healthy/disease case) RNA sequencing gives information about the transcriptome of the samples quantitatively. Comparing the RNA levels in these 2 conditions, it is possible to make interpretation about the differentially expressed genes. However, before diving into DGE analysis, it is important to do the preprocessing. Here, I have a simple RNA Seq preprocessing workflow to do quality check and perform alignment. 

The output of the RNA Sequencing is the fastq files, containing the read sequences and some information about the read, such as the sequencing run, cluster and the quality score of the bases. This information can be used to clean the low quality reads from the fastq files. For this first I will assess the read quality using [Fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). During the alignment, if a read has high number of mismatches it will not be aligned. So in order to reduce the number of mismatches between the query sequence and reference data, [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) will be used to trim the adaptor sequencing, remove the low quality bases and get rid of any contamination. However, this step is optinal as many modern alignmnet tools can account for these. For the alignment step [STAR Aigner]( is used to map the reads to the reference genome. The output alignment 


## Data <a name="data"></a>

- The *E. Coli* data [this paper](https://www.nature.com/articles/s41598-019-51046-7#Sec10)
- The reference [*E. Coli* genome](https://www.ncbi.nlm.nih.gov/assembly/GCF_002848225.1/)

The study aims to find the effects of Wi-Fi radiofrequency radiation of 2.4 GHz on global gene expression. For this the data is an RNA seq from control and radio frequency exposed bacteria. I just perform the preprocessing on the raw sequence data.

## Pipeline <a name="pipeline"></a>

![alt text](https://github.com/iremycl/AlignmentWorkflow_Snakemake/blob/main/dag.svg)

1. QC Check : `FastQC`
2. Adaptor Trimming: `Trimmomatic`
3. Alignment to reference genome:  `STAR` 
4. Alignment Quality Check: `samtools`


## Usage <a name="usage"></a>
Usage

### Input <a name="input"></a>

The input data consists of:

|File|Description|
|---|---|
|GCF_002848225.1_ASM284822v1_genomic.fna.gz| *E. Coli* reference genome fasta file|
|GCF_002848225.1_ASM284822v1_genomic.gtf.gz| *E. Coli*  genome GTF file|
|SRR8580014_1.fastq.gz|Control *E. Coli* fastq file Read1|
|SRR8580014_2.fastq.gz|Control *E. Coli* fastq file Read2|
| TruSeq3-PE.fa| Adaptor sequences for Trimming |


### Output <a name="output"></a>



## References <a name="references"></a>
1. Said-Salman, I.H., Jebaii, F.A., Yusef, H.H. et al. Global gene expression analysis of Escherichia coli K-12 DH5α after exposure to 2.4 GHz wireless fidelity radiation. Sci Rep 9, 14425 (2019). https://doi.org/10.1038/s41598-019-51046-7
2. Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170.
3. FastQC: A Quality Control Tool for High Throughput Sequence Data [Online]. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ (2015), "FastQC," https://qubeshub.org/resources/fastqc.
4. Alexander Dobin, Carrie A. Davis, Felix Schlesinger, Jorg Drenkow, Chris Zaleski, Sonali Jha, Philippe Batut, Mark Chaisson, Thomas R. Gingeras, STAR: ultrafast universal RNA-seq aligner, Bioinformatics, Volume 29, Issue 1, January 2013, Pages 15–21, https://doi.org/10.1093/bioinformatics/bts635
5. Heng Li, Bob Handsaker, Alec Wysoker, Tim Fennell, Jue Ruan, Nils Homer, Gabor Marth, Goncalo Abecasis, Richard Durbin, 1000 Genome Project Data Processing Subgroup, The Sequence Alignment/Map format and SAMtools, Bioinformatics, Volume 25, Issue 16, 15 August 2009, Pages 2078–2079, https://doi.org/10.1093/bioinformatics/btp352

