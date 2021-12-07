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

The output of the RNA Sequencing is the fastq files, containing the read sequences and some information about the read, such as the sequencing run, cluster and the quality score of the bases. This information can be used to clean the low quality reads from the fastq files. For this first I will assess the read quality using [Fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). During the alignment, if a read has high number of mismatches it will not be aligned. So in order to reduce the number of mismatches between the query sequence and reference data, [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) will be used to trim the adaptor sequencing, remove the low quality bases and get rid of any contamination. However, this step is optinal as many modern alignmnet tools can account for these. For the alignment step [STAR Aigner](https://github.com/alexdobin/STAR) is used to map the reads to the reference genome. Then [samtools](https://github.com/samtools/samtools) is used to quality check the output alignment file. 

This is a very general workflow and there are many more tools that can be used for each step. In this pipeline I used the most commonly used tools for each step, so that this can be used by many user as the first step of their analysis. After this step, many further analysis can be done, such as variant calling, DGEA, or GSEA. 

## Data <a name="data"></a>

The study aims to find the effects of Wi-Fi radiofrequency radiation of 2.4 GHz on global gene expression. For this the data is an RNA seq from control and radio frequency exposed bacteria. I just perform the preprocessing on the raw sequence data. There are RNA Seq data from 3 controls and 3 exposed *E. Coli*, but I will only use 1 of the control samples for my workflow. 

- The *E. Coli* data from [this paper](https://www.nature.com/articles/s41598-019-51046-7#Sec10)
- The reference [*E. Coli* genome](https://www.ncbi.nlm.nih.gov/assembly/GCF_002848225.1/)


## Pipeline <a name="pipeline"></a>

![alt text](https://github.com/iremycl/AlignmentWorkflow_Snakemake/blob/main/dag.svg)

The main steps of this workflow consists of:

1. QC Checkusing `FastQC`
2. Adaptor Trimming with `Trimmomatic`
3. Alignment to reference genome using  `STAR` 
4. Alignment Quality Check with `samtools`

In the workflow the versions are:

  - python =3.8
  - FastQC=v0.11.9
  - trimmomatic=0.39
  - STAR=2.7.9a
  - samtools=1.14
  - 
## Usage <a name="usage"></a>

To use in OSx:

To use the same workflow, just clone the repository:
```
git clone https://github.com/iremycl/AlignmentWorkflow_Snakemake.git
cd AlignmentWorkflow_Snakemake
```
and create the conda environment using the environment.yaml :

```
conda env create --file environment.yml

```
Then the pipeline can be run in the same directory with the Snakemake file:

```
snakemake --cores 2

```


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

The output data of the first run is the html files showing the base qualities for the reads in the input Fastq files. Then the trimming part removes the reads with high mismatches, such as adaptors and contaminants Then the trimmed reads are used along with the indexed reference genome for the alignment. The final output is from samtools, reporting the alignment statistics.

## References <a name="references"></a>
1. Said-Salman, I.H., Jebaii, F.A., Yusef, H.H. et al. Global gene expression analysis of Escherichia coli K-12 DH5α after exposure to 2.4 GHz wireless fidelity radiation. Sci Rep 9, 14425 (2019). https://doi.org/10.1038/s41598-019-51046-7
2. Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170.
3. FastQC: A Quality Control Tool for High Throughput Sequence Data [Online]. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ (2015), "FastQC," https://qubeshub.org/resources/fastqc.
4. Alexander Dobin, Carrie A. Davis, Felix Schlesinger, Jorg Drenkow, Chris Zaleski, Sonali Jha, Philippe Batut, Mark Chaisson, Thomas R. Gingeras, STAR: ultrafast universal RNA-seq aligner, Bioinformatics, Volume 29, Issue 1, January 2013, Pages 15–21, https://doi.org/10.1093/bioinformatics/bts635
5. Heng Li, Bob Handsaker, Alec Wysoker, Tim Fennell, Jue Ruan, Nils Homer, Gabor Marth, Goncalo Abecasis, Richard Durbin, 1000 Genome Project Data Processing Subgroup, The Sequence Alignment/Map format and SAMtools, Bioinformatics, Volume 25, Issue 16, 15 August 2009, Pages 2078–2079, https://doi.org/10.1093/bioinformatics/btp352

