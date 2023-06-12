# Optimized denatured RiboMeth-seq (dRMS) for small RNAs
# 1.	Overview
Despite recent advancements, Nm site quantification remains challenging. The commonly used RiboMeth-seq (RMS) hydrolyze RNA at high pH, where 2’-O-methylation protects the phosphodiester bond, which is detected by sequencing. However, stable RNA structures and dense modifications, especially for tRNAs, cause strong biases in the fragmentation and reverse transcription. Here, we developed a denatured RMS (dRMS) method by combining alkaline pH, high temperature and DMSO, which increased the efficiency and uniformity of RNA fragmentation.

To avoid the false positive results caused by other modifications on small RNAs, the ribosome methylation score (RiboMethScore) was computed with four neighboring nucleotides (±2 nucleotide window) from 5’-end and 3’-end reads at nucleotide resolution. The empirical rule with 1.5 standard deviation (µ ± 1.5*σ, the probability of 1.5 standard deviations is about 86.6%) was also used to filter the discrete values. If there are observations that fall outside the 1.5 standard deviation, only two neighboring nucleotides (±1 nucleotide window) from 5’-end and 3’-end reads were used to calculate RiboMethScore.

The RiboMethScore calling script is written in Python and available as source code that could be downloaded and run directly on local machine (no compiling needed). The dRMS analysis pipeline consists of five steps:

* STAR mapping
* Calculate coverage of 5' and 3' positions of mapped reads
* RiboMethScore calling

# 2.	Overview
