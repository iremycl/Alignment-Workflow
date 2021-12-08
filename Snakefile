# configfile: "environment.yaml"

SAMPLE = ["SRR8580014"]

rule all:
    input:
        "data/samples/raw/SRR8580014_fastqc.html",
        "data/samples/raw/SRR8580014_fastqc.zip",
        # "data/samples/trimmed/SRR8580014_1P.fastq.gz",
        # "data/samples/trimmed/SRR8580014_2P.fastq.gz",
        # "data/samples/trimmed/SRR8580014_1U.fastq.gz",
        # "data/samples/trimmed/SRR8580014_2U.fastq.gz",
        "data/reference/index",
        "data/samples/mapped/SRR8580014_Aligned.sortedByCoord.out.bam",
        "results/SRR8580014_Aligned.sortedByCoord.out.bam.flagstat"


rule fastqc:
    input:
        "data/samples/raw/SRR8580014_1.fastq.gz",
        "data/samples/raw/SRR8580014_2.fastq.gz"
    output:
        "data/samples/raw/SRR8580014_fastqc.html",
        "data/samples/raw/SRR8580014_fastqc.zip"
    shell:
        "fastqc {input} -o {output}"

# rule trimming:
#     input:
#         s1r1="data/samples/raw/SRR8580014_1.fastq.gz",
#         s1r2="data/samples/raw/SRR8580014_2.fastq.gz"
#     output:
#         s1r1_paired="data/samples/trimmed/SRR8580014_1P.fastq.gz",
#         s1r2_paired="data/samples/trimmed/SRR8580014_2P.fastq.gz",
#         s1r1_unpaired="data/samples/trimmed/SRR8580014_1U.fastq.gz",
#         s1r2_unpaired="data/samples/trimmed/SRR8580014_2U.fastq.gz"
#     shell:
#         "trimmomatic PE -threads 4 {input.s1r1} {input.s1r2} {output.s1r1_paired} {output.s1r1_unpaired} {output.s1r2_paired} {output.s1r2_unpaired} ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:15"



rule generate_genome_index:
    input:
        fasta="data/reference/GCF_002848225.1_ASM284822v1_genomic.fna",
        gtf="data/reference/GCF_002848225.1_ASM284822v1_genomic.gtf"
    output:
        "data/reference/index/chrLength.txt",
        "data/reference/index/chrNameLength.txt",
        "data/reference/index/chrName.txt",
        "data/reference/index/chrStart.txt",
        "data/reference/index/exonGeTrInfo.tab",
        "data/reference/index/exonInfo.tab",
        "data/reference/index/geneInfo.tab",
        "data/reference/index/Genome",
        "data/reference/index/genomeParameters.txt",
        "data/reference/index/Log.out",
        "data/reference/index/SA",
        "data/reference/index/SAindex",
        "data/reference/index/sjdbInfo.txt",
        "data/reference/index/sjdbList.fromGTF.out.tab",
        "data/reference/index/sjdbList.out.tab",
        "data/reference/index/transcriptInfo.tab"
    shell:
        "STAR --runThreadN {threads} "
        "--runMode genomeGenerate "
        "--genomeDir data/reference/index "
        "--genomeFastaFiles {input.fasta} "
        "--sjdbGTFfile {input.gtf} "
        "--sjdbOverhang 100 "
        "--genomeSAindexNbases 10 "

rule unzip_1:
    input:
        s1r1="data/samples/raw/SRR8580014_1.fastq.gz"
    output:
        s1r1="data/samples/raw/SRR8580014_1.fastq"
    shell:
        "gzip -dc {input} > {output}"


rule unzip_2:
    input:
        s1r2="data/samples/raw/SRR8580014_2.fastq.gz"
    output:
        s1r2="data/samples/raw/SRR8580014_2.fastq"
    shell:
        "gzip -dc {input} > {output}"

rule star_map:
    input:
        "data/reference/index/chrLength.txt",
        "data/reference/index/chrNameLength.txt",
        "data/reference/index/chrName.txt",
        "data/reference/index/chrStart.txt",
        "data/reference/index/exonGeTrInfo.tab",
        "data/reference/index/exonInfo.tab",
        "data/reference/index/geneInfo.tab",
        "data/reference/index/Genome",
        "data/reference/index/genomeParameters.txt",
        "data/reference/index/Log.out",
        "data/reference/index/SA",
        "data/reference/index/SAindex",
        "data/reference/index/sjdbInfo.txt",
        "data/reference/index/sjdbList.fromGTF.out.tab",
        "data/reference/index/sjdbList.out.tab",
        "data/reference/index/transcriptInfo.tab",
        s1r1="data/samples/raw/SRR8580014_1.fastq",
        s1r2="data/samples/raw/SRR8580014_2.fastq"
    output:
        "data/samples/mapped/SRR8580014_Aligned.sortedByCoord.out.bam"
    shell:
        "STAR --runThreadN {threads} "
        " --readFilesIn {input.s1r1} {input.s1r2} "
        " --genomeDir data/reference/index "
        " --outFileNamePrefix data/samples/mapped/{SAMPLE}_ "
        " --outSAMtype BAM SortedByCoordinate "
        " --outSAMunmapped Within "
        " --outSAMattributes Standard "
        " --alignIntronMax 20000 "

rule samtools_alignment_qc:
    input:
        "data/samples/mapped/SRR8580014_Aligned.sortedByCoord.out.bam"
    output:
        "results/SRR8580014_Aligned.sortedByCoord.out.bam.flagstat"
    shell:
        "samtools flagstat {input} > {output}"
