# Optimized denatured RiboMeth-seq (dRMS) for small RNAs
# 1. Overview
Despite recent advancements, Nm site quantification remains challenging. The commonly used RiboMeth-seq (RMS) hydrolyze RNA at high pH, where 2’-O-methylation protects the phosphodiester bond, which is detected by sequencing. However, stable RNA structures and dense modifications, especially for tRNAs, cause strong biases in the fragmentation and reverse transcription. Here, we developed a denatured RMS (dRMS) method by combining alkaline pH, high temperature and DMSO, which increased the efficiency and uniformity of RNA fragmentation.

To avoid the false positive results caused by other modifications on small RNAs, the ribosome methylation score (RiboMethScore) was computed with four neighboring nucleotides (±2 nucleotide window) from 5’-end and 3’-end reads at nucleotide resolution. The empirical rule with 1.5 standard deviation (µ ± 1.5*σ, the probability of 1.5 standard deviations is about 86.6%) was also used to filter the discrete values. If there are observations that fall outside the 1.5 standard deviation, only two neighboring nucleotides (±1 nucleotide window) from 5’-end and 3’-end reads were used to calculate RiboMethScore.

The RiboMethScore calling script is written in Python and available as source code that could be downloaded and run directly on local machine (no compiling needed). The dRMS analysis pipeline consists of five steps:

* Step1: STAR mapping
* Step2: Calculate coverage of 5' and 3' positions of mapped reads
* Step3: RiboMethScore calling

# 2. System and environment requirements
Download the scripts and save it to a known path/location. No special installation is needed, but the python package dependencies need to be properly resolved before use. Additional tools used for mapping and general processing of high throughput sequencing data, including STAR, samtools and bedtools. For visualization of the results, we recommend IGV.

• STAR 

• samtools

• bedtools

• IGV v2.4+


# 3. dRMS analysis
## 3.1 STAR mapping

The input fastq files need to be mapped to related genomes. Here, STAR was shown as an example for the reads alignment.

    STAR --runMode alignReads --runThreadN 8 --genomeDir staridx_2.7.9 --genomeLoad NoSharedMemory --readFilesIn Sample.fastq --outFileNamePrefix Outprefix --outStd Log --outReadsUnmapped Fastx --outSAMtype BAM Unsorted SortedByCoordinate --outSAMmode Full --outSAMattributes All --outFilterType BySJout  --outFilterMultimapNmax 80 --outFilterMultimapScoreRange 1 --outFilterScoreMin 1 --outFilterScoreMinOverLread 0.1 --outFilterMatchNmin 15 --outFilterMatchNminOverLread 0.1 --outFilterMismatchNmax 5 --outFilterMismatchNoverLmax 0.1 --outFilterMismatchNoverReadLmax 0.1 --alignIntronMin 20 --alignIntronMax 1000000

## 3.2 Calculate coverage of 5' and 3' positions of mapped reads

Bedtools genomecov was used to summarize the feature coverage of 5' and 3' positions of mapped reads.

    bedtools genomecov -ibam Sample.bam -5 -bg > Sample_5end.bedgraph

    bedtools genomecov -ibam Sample.bam -3 -bg > Sample_3end.bedgraph

## 3.3 RiboMethScore calling

The script RiboMethScore_calling.py could be used to analyze the ribosome methylation score (RiboMethScore) with four neighboring nucleotides (±2 nucleotide window) from 5’-end and 3’-end reads at nucleotide resolution. Only the RiboMethScore of target sites with more than 20 mapped reads will be outputted. 

•**Gene bed** file should contain six required BED fields: Chrome, ChromStart, ChromEnd, Name, Score, Strand.
  
    $ python RiboMethScore_calling.py  Sample_5end.bedgraph  Sample_3end.bedgraph  gene.bed  5end  Outputprefix
  
• The output file contains four columns:

• Column1: Genomic position of target sites.

• Column2: Gene name.

• Column3: Relative genetic position of target sites.

• Column4: The coverage of 5' positions of mapped reads. 

• Column5: The RiboMethScore calculated only by 5'end coverage.

    $ python RiboMethScore_calling.py  Sample_5end.bedgraph  Sample_3end.bedgraph  gene.bed  3end  Outputprefix
  
• The output file contains four columns:

• Column1: Genomic position of target sites.

• Column2: Gene name.

• Column3: Relative genetic position of target sites.

• Column4: The coverage of 3' positions of mapped reads. 

• Column5: The RiboMethScore calculated only by 3'end coverage.


    $ python RiboMethScore_calling.py  Sample_5end.bedgraph  Sample_3end.bedgraph  gene.bed  both  Outputprefix
  
• The output file contains four columns:

• Column1: Genomic position of target sites.

• Column2: Gene name.

• Column3: Relative genetic position of target sites.

• Column4: The coverage of 5' positions of mapped reads. 

• Column5: The coverage of 3' positions of mapped reads. 

• Column6: The sum coverage of 5' and 3' positions of mapped reads.

• Column7: The RiboMethScore calculated by 5' + 3' end coverage.


