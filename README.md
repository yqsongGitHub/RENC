
<!-- README.md is generated from README.Rmd. Please edit that file -->

</div>

<h2 align="center"> RENC: Ranking of target genes based on ENhancer-promoter Contacts </h2>

<p align="center">
<img align="center" src="https://github.com/yqsongGitHub/RENC/blob/main/Images/Schematic.png">
</p>   

------------------------------------------------------------------------

<!-- badges: start -->
<!-- badges: end -->
## Introduction 
We developed a methodology called Ranking of target genes based on ENhancer-promoter Contacts (RENC), which utilizes HiChIP data to prioritize target genes associated with putative enhancers located within duplication hotspots.

For enhancers in each duplication hotspot, we use HiChIP paired-end tags (PETs) to measure the amount of enhancer activity delivered to each gene’s promoter via significant HiChIP loops (hotspot-delivered enhancer activity). Then, for each of the candidate genes, we measure the relative contribution of the hotspot-delivered enhancer activity as percentage of activity delivered from all enhancers (inside and outside of the hotspot) interacting with the gene’s promoter. We multiply these two factors to prioritize target genes for the enhancers in the duplication hotspots.

The formula is as follows：
<p align="center">
<img align="center" src="https://github.com/yqsongGitHub/RENC/blob/main/Images/Formula.png" alt="Formula" width="50%" height="50%">
</p>   

For a given duplication hotspot, $`E_{in}`$ is the sum of the PETs from all significant HiChIP loops (PETs >= 3 and FDR < 0.05) connecting the enhancers within the hotspot to a candidate gene’s promoter ($`{EP}_{in}`$): $`E_{in}=\ \sum_{i\in{EP}_{in}}{PET}_i`$. $`E_{out}`$ is the sum of the PETs  from all significant HiChIP loops connecting the enhancers outside the hotspot to the candidate gene’s promoter ($`{EP}_{out}`$):  $`E_{out}=\ \sum_{i\in{EP}_{out}}{PET}_i`$. Therefore, $E_{in}$ represents the absolute enhancer activity delivered from the hotspot to each gene promoter. $`\frac{E_{in}}{E_{in}+E_{out}}`$ represents the relative contribution of hotspot-delivered enhancer activity as percentage of activity delivered from all its enhancers (including the enhancers inside or outside the duplication hotspots). The multiplication of the two factors: the absolute enhancer activity delivered (PETs) and their relative contribution to each gene promoter (%), is used for ranking the target genes for each hotspot. 

For genes that have multiple promoters, each of the promoters will be treated separately and the promoter with the highest RENC score will be used to represent the gene. Genes not linked to the duplication hotspot through a significant enhancer-promoter HiChIP loop are given a RENC score of zero. For cancer types that have HiChIP data available from more than one cell line, we selected the prioritized gene promoters that are shared. As the RENC score is calculated per duplication hotspot, genes linked to different hotspots are not compared to each other.

The RENC methodology provides a framework to prioritize target genes for enhancers in duplication hotspots. RENC is designed with reference to [BEDTools](http://quinlanlab.org/tutorials/bedtools/bedtools.html) for command-line style programming. We provide a shell script and user guide in this github page.

------
## Requirements
For each cell-type, the inputs to the RENC methodology are:

* Required Inputs
 	* 1. BED files for the genomic coordinates (e.g. duplication hotspot regions) 
 	* 2. BEDPE files derived from HiChIP results
 	* 3. BED files for transcription start site（TSS）of protein coding genes
 	* 4. BED files for gene body of protein coding genes
* Optional Inputs
 	* 5. specific genes or the genomic coordinates of duplication hotspots to search; e.g.,"C4BPB","chr1:201970000-202085000"or "C4BPB;chr1:201970000-202085000"  

------

## Basic Usage and Quick Guide
### Example and reference data background introduction

Example data for testing is available at [RENC/example](https://github.com/yqsongGitHub/RENC/tree/main/example/input). 
 * 1. BED file for the genomic coordinates of duplication hotspots of stomach adenocarcinoma (STAD), derived from Pan-Cancer Atlas of Whole Genomes (PCAWG) data.
 * 2. BEDPE file of significant loops derived from HiChIP results for chromosome 1 in the STAD cell line AGS.

Reference data is available at [RENC/reference](https://github.com/yqsongGitHub/RENC/tree/main/reference/hg19).
 * 3. BED file for transcription start site（TSS）of protein coding genes downloaded from  [UCSC](https://genome.ucsc.edu/)
 * 4. BED file for gene body of protein coding genes downloaded from  [UCSC](https://genome.ucsc.edu/)

### Installation
#### 1. Create Conda environment
By creating a Conda environment, you can have separate environments for different projects or purposes. Ensure that the packages and dependencies within each environment do not conflict with each other.
```
conda create -n RENC
```
#### 2. Activate your Conda environment
Once the environment is created, you can activate it and start installing packages or running programs within that environment. Ensure that the installed packages and dependencies are isolated from your system's global environment.
```
conda activate RENC
```
#### 3. Install bedtools using the Conda package manager
Install bedtools in your Conda environment.
```
conda install -c bioconda bedtools
# If the installation is successful, the script below should run without any issues
bedtools -h
```
#### 4. Download RENC script via GitHub
It is recommended to download using ***git clone***. However, in case of any network issues, you can download directly by clicking the "Download ZIP" button on the GitHub website. By downloading the script, you acquire a local copy of the file that you can use and modify as needed on your own machine.
```
git clone https://github.com/yqsongGitHub/RENC.git
cd RENC
```
Then, you can use the following example command lines to test the functionality of RENC after installation.

### Routine analysis 1: Get all the enhancer-gene regulatory relationships from the input BEDPE file and BED file.
```
bash ./src/RENC.sh -r ./src/RENC.R  \
-d ./example/input/Duplication_STAD.bed \
-c ./example/input/AGS_STAD_chr1.bedpe  \
-t ./reference/hg19/RefSeq_proteinCoding.tss.bed \
-g ./reference/hg19/RefSeq_proteinCoding.body.bed 
```
Please note, in RENC, all files using the reference genome must be consistent, with the example using hg19. If there are files using the hg38 reference genome, the [LiftOver](https://liftover.broadinstitute.org/) tool can be used to convert the hg38 genome to hg19.

The output is a .RENC.txt file with annotation of information as follows.  

column | name | explanation
------ | ---- | ------------
1th | chr1 | chromosomal for the loop first anchor
2th | start1 | genomic coordinate of the start site for the first anchor
3th | end1 | genomic coordinate of the end site for the first anchor
4th | chr2 | chromosomal for the loop second anchor
5th | start2 | genomic coordinate of the start site for the second anchor
6th | end2 | genomic coordinate of the end site for the second anchor
7th | pets | observed paired-end tags number linking the two anchors
8th | chr | chromosomal for the duplication
9th | start | genomic coordinate of the start site for the duplication
10th | end | genomic coordinate of the end site for the duplication
11th | dup_sum_pets | the amount of duplicated enhancer activity delivered to the candidate gene
12th | gene_ID_sum_pets | the amount of all enhancers activity delivered to the candidate gene
13th | percentage | the relative contribution of the duplicated enhancer activity as percentage of activity delivered from all enhancers interacting with the candidate gene
14th | score | RENC score for the enhancer interacting with the gene is calculated by multiplying the amount of duplicated enhancer activity (dup_sum_pets) by the relative contribution of the duplicated enhancer activity (percentage)
15th | rank_gene | the rank of the target genes for the enhancers in the duplication hotspot
16th | gene | the specific symbol or identifier associated with the candidate gene
17th | ID | the transcript ID or identifier associated with the candidate gene
18th | rank_enhancer | the rank of the enhancers in the duplication hotspot delivered to the candidate gene
19th | enhancer | genomic coordinates of the linked enhancers


### Routine analysis 2: Get the target genes from input genomic coordinates of specific duplication hotspots.

If using multiple regions as input, please separate them with a semicolon (;). If no results are found, return a placeholder (.).

```
bash ./src/RENC.sh -r ./src/RENC.R  \
-d ./example/input/Duplication_STAD.bed \
-c ./example/input/AGS_STAD_chr1.bedpe  \
-t ./reference/hg19/RefSeq_proteinCoding.tss.bed \
-g ./reference/hg19/RefSeq_proteinCoding.body.bed \
-s  "chr1:201970000-202085000;chr1:33135000-33370000"
```

The output is a .RENC.search.txt file with annotation of information as follows.  
column | name | explanation
------ | ---- | ------------
1th | enhancer | the input genomic coordinates of specific duplication hotspots
2th | gene | the target gene. If no results are found, return a placeholder (.)

### Routine analysis 3: Get genomic coordinates of the duplicated enhancers associated with input target genes.

If using multiple regions as input, please separate them with a semicolon (;). If no results are found, return a placeholder (.).

```
bash ./src/RENC.sh -r ./src/RENC.R  \
-d ./example/input/Duplication_STAD.bed \
-c ./example/input/AGS_STAD_chr1.bedpe  \
-t ./reference/hg19/RefSeq_proteinCoding.tss.bed \
-g ./reference/hg19/RefSeq_proteinCoding.body.bed \
-s  "C4BPB;ELF3"
```
The output is a .RENC.search.txt file with annotation of information as follows.  
column | name | explanation
------ | ---- | ------------
1th | enhancer | genomic coordinate of duplicated enhancers associated with input target genes. If no results are found, return a placeholder (.)
2th | gene | the input target gene

------
## RENC Main Functions

Run ***bash RENC.sh*** or ***bash RENC.sh -h*** can show the main functions of RENC.sh with short descriptions and examples.     

```
Usage:   
  bash RENC.sh [-r <RSCRIPT>][-d <DUPLICATION_FILE>] [-c <HICHIP_FILE>][-t <TSS_FILE>][-g <GENEBODY_FILE>]
  -r Rscript          Rscript file [required].
  -d DUPLICATION_FILE Duplication or region file [required].
  -c HICHIP_FILE      HiChIP file from Hichipper pipeline or other bedpe file [required].
  -t TSS_FILE         TSS file for protein coding genes or other TSS file [required].
  -g GENEBODY_FILE    Gene body file for protein coding genes [required].
  -s GENE_ENH         Gene or region in Duplication to search; ALL gene and region by default. e.g.,"C4BPB","chr1:201970000-202085000" or "C4BPB:chr1:201970000-202085000"

Example:     
  bash ./src/RENC.sh -r ./src/RENC.R  
  -d ./example/input/Duplication_STAD.bed  
  -c ./example/input/AGS_STAD_chr1.bedpe  
  -t ./reference/hg19/RefSeq_proteinCoding.tss.bed  
  -g ./reference/hg19/RefSeq_proteinCoding.body.bed 

```

------
## Example: KLF5 is the top candidate target in squamous cell carcinoma
### Step1:  Get all the enhancer-gene regulatory relationships from the input BEDPE file and BED file. 
The ***Duplication_Squamous.bed*** file contains the genomic coordinates of duplication hotspots of squamous cell carcinoma. The ***SqCC_BICR31.bedpe*** file, derived from HiChIP experiments, is mapped to hg19 in the BICR31 (squamous cell carcinoma) cell line. The ***output*** provides all the enhancer-gene regulatory relationships and the RENC score for each candidate gene, including KLF5.
```
bash ./src/RENC.sh -r ./src/RENC.R  \
-d ./example/input/Duplication_Squamous.bed \
-c ./example/input/SqCC_BICR31.bedpe  \
-t ./reference/hg19/RefSeq_proteinCoding.tss.bed \
-g ./reference/hg19/RefSeq_proteinCoding.body.bed
```

### Step2: Visualize the region near KLF5.
#### Step2.1：Set global variables and load R packages.
We used R version 4.2.2 and the [gTrack](https://github.com/mskilab-org/gTrack) package for visualization.
We also need to load some functions.
```
library(gTrack)
library(gUtils)
library(rtracklayer)
source("./src/gTrack.R")
source("./src/function.R")
options(scipen = 100)
options(warn = -1)
```

#### Step2.2：Load the required files.
Here we load some necessary files, including the ***SqCC_BICR31_RENC.txt*** and ***SqCC_BICR31_tss.bed*** files generated in step 1.
```
gt.ge <- readRDS("./data/gt.ge.rds")
genes <- read.table("./reference/hg19/RefSeq_proteinCoding.body.bed")
tss <- read.table("./reference/hg19/RefSeq_proteinCoding.tss.bed")
chip <- './data/GSM2356644_H3K27ac-BICR31.bw'
ins <- read.table("./data/BICR31_KLF5.bedpe")
dup <- read.table("./data/dup.Squamous.merged.1mb.sort.genomecov")
hs <- read.table("./example/input/Duplication_Squamous.bed")
pps <- read.table("./example/output/SqCC_BICR31_RENC.txt",header = T)
bp <- read.table("./example/output/SqCC_BICR31_tss.bed")
```

#### Step2.3：Visualize the region near KLF5
We selected the region from 73,205,000 to 74,345,000 on chromosome 13 for presentation, which includes H3K27ac HiChIP signal, positions of duplication hotspots, duplication events observed in squamous cancer, H3K27ac ChIP-seq signal, RENC scores prioritizing the KLF5 promoter that is more likely to be activated by enhancers within the left duplication hotspot presented in the window, and the number of PETs contributed by each interacting enhancer to KLF5. 

Here is the plot of KLF5 that we created. Further editing and enhancement using [Adobe Illustrator](https://www.adobe.com/products/illustrator/free-trial-download.html) is required.
<p align="center">
<img align="center" src="https://github.com/yqsongGitHub/RENC/blob/main/Images/BICR31_KLF5.png" width="80%" height="80%">
</p>  

```
Gene <- "KLF5"
Chr <- "chr13"
Chr.start <-  73205000
Chr.end <- 74345000
Height <- 3
Ygap <- 0.8

hs <- hs[hs$V1=="chr13" & hs$V2>73210000 & hs$V3 < 74340000,]
hs1 <- hs[1,] 
pps1 <- pps[pps$chr == hs1$V1 & pps$start == hs1$V2 & pps$end == hs1$V3,]
bp <- find_enh_anchor(bp)
hara <- draw_gtrack(ins,chr=Chr,chr.start = Chr.start,chr.end=Chr.end)
dat.gt <- hara[[1]]
wins <- hara[[2]]
ind <- which(names(gt.ge@data[[1]]) %in% genes$V4)
gt.ge@data[[1]] <- gt.ge@data[[1]][ind] 
gt.ge@formatting$height <- Height+2
gt.ge@formatting$ygap <- Ygap

# H3K27ac ChIP-seq signal 
chip.gr.score <- gTrack(chip,col = '#858480', name = 'chip',bars = TRUE,height = Height+0.5, ygap = Ygap)

# Duplication hotspot
hs.gr <- GRanges(seqnames = Rle(hs$V1), 
                 ranges = IRanges(start = hs$V2,
                                  end = hs$V3))
hs.gr <- gTrack(hs.gr,col = '#198330', name = '',bars = TRUE,height = Height/2, ygap =Ygap/2)

# Duplication events 
dup.gr <- GRanges(seqnames = Rle(dup$V1), 
                  ranges = IRanges(start = dup$V2,end = dup$V3),
                  score =dup$V4)
dup.gr <- gTrack(dup.gr, y.field = 'score',
                 col = '#858480', name = 'dup',bars = TRUE,height = Height+0.5, ygap =Ygap, y0=1)

# RENC scores 
pps1 <- pps1[,c("ID","score")]
pps1 <- merge(tss,pps1,by.x="V5",by.y="ID",all.x=T)
pps1$score[is.na(pps1$score)] <- 0
pps1 <- pps1[!(duplicated(pps1$V4) & pps1$score==0),]
pps.gr1 <- GRanges(seqnames = Rle(pps1$V1), 
                   ranges = IRanges(start = pps1$V2,end = pps1$V3),
                   score =pps1$score)
pps.gr1 <- gTrack(pps.gr1, y.field = 'score',
                  col = '#858480', name = 'RENC',circles = TRUE,height = Height+0.5, ygap = Ygap)

# Enhancer contributions
bp <- bp[bp$gene==Gene,]
enr.gr <- GRanges(seqnames = Rle(bp$chr), 
                  ranges = IRanges(start = bp$start,end = bp$end),
                  pets =bp$pets)
enr.gr <- gTrack(enr.gr, y.field = 'pets',
                 col = '#858480', name = 'contrib.',circles = TRUE, height = Height+1,ygap = Ygap)

pdf("./Images/BICR31_KLF5.pdf",height = 12,width = 8)
plot(c(enr.gr,pps.gr1,chip.gr.score,dup.gr,hs.gr,gt.ge,dat.gt),wins)
dev.off()
```

### Step3 : Enhance the image using Adobe Illustrator.
In step 2, we obtained the plot of KLF5. Then, we used [Adobe Illustrator](https://www.adobe.com/products/illustrator/free-trial-download.html) to enhance the representation of protein-coding genes and annotated the previously published functional enhancers for KLF5.

<p align="center">
<img align="center" src="https://github.com/yqsongGitHub/RENC/blob/main/Images/BICR31_KLF5_ai.png" width="80%" height="80%">
</p>  

--------
## Questions & Answers
Please submit a github issue with any questions or if you experience any issues/bugs.
- Create a new [issue](https://github.com/yqsongGitHub/RENC/issues).

-------
## Warning   
No tested bugs needed warning. 

--------
## Citation
Coming soon