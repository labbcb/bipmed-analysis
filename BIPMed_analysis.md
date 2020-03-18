# BIPMed data analysis

Author: Rodrigo Secolin  
Date: September 25th, 2019

## Introduction

Here we generate a descriptive statistics from two BIPMed VCF files: one including WES and one including SNP array genotype data, and compare allele frequencies with [1000 Genome dataset](https://doi.org/10.1038/nature15393). BIPMed data were evaluated based on hg19 genome built. We changed directory and file names to maintain privacy. Commands are performed by `R software v3.6.3` and `PLINK v1.9`.

## Analysis

### Load require package for plot

```r
library(ggplot2)
```

### Merging SNP array and WES datasets

Remove samples from WES data according heterogeneity, IBS, and relatedness parameters from [Secolin et al., 2019](https://doi.org/10.1038/s41598-019-50362-2), and convert in PLINK format.
```r
system('plink --vcf wes.vcf.gz --make-bed --out wes_release --double-id --biallelic-only strict --remove remove_sample_ids.txt')
# change population definition
wes.fam<-read.table('wes_release.fam')
wes.fam$V1<-"BRS"
write.table(wes.fam,file='wes_release.fam',quote = F,row.names = F,col.names = F,sep="\t")
```

Change code "." to mutation nomenclature in WES *.bim* file. 
```r
# Change code "." to mutation nomenclature in WES data from bim file.
wes.bim<-read.table('wes_release.bim',stringsAsFactors = F)
for(i in 1:nrow(wes.bim)){
  if(wes.bim$V2[i]=="."){
    wes.bim$V2[i]<-paste0("chr",wes.bim$V1[i],"_",wes.bim$V4[i],wes.bim$V5[i],"_",wes.bim$V6[i])
  }
}
write.table(wes.bim,file='wes_release.bim',quote = F,row.names = F,col.names = F,sep="\t")
```
Remove samples from SNP array data according parameters from Secolin et al., 2019 (doi:10.1038/s41598-019-50362-2), and removing SNPs with MAF < 0.01 due to error rate of SNP array genotyping.
```r
system('plink --vcf array_raw.vcf.gz --remove remove_sample_ids.txt --maf 0.01 --make-bed --out array_id_filtered_maf')
# change population definition
array.fam<-read.table('array_id_filtered_maf.fam')
array.fam$V1<-"BRS"
write.table(array.fam,file='array_id_filtered_maf.fam',quote = F,row.names = F,col.names = F,sep="\t")

```
Merge SNP array and WES datasets
```r
system('plink --bfile array_id_filtered_maf --bmerge wes_release --make-bed --out array_wes')

## flip WES
system('plink --bfile wes_release --flip array_wes-merge.missnp --make-bed --out wes_release_flip')

## merging again
system('plink --bfile array_id_filtered_maf --bmerge wes_release_flip --make-bed --out array_wes')
```

### Merge BIPMed dataset and 1000 Genome dataset

Remove variants with more than 5% of missing data, ambiguous variants, variants with Hardy-Weinberg disequilibrium and non-autossomal variants.
```r
#Remove ambiguous variants in R
bsp.aut<-read.table('array_wes.bim')
bsp.aut<-cbind(bsp.aut, gen=paste0(bsp.aut$V5,bsp.aut$V6))
transitions<-c("GC","CG","AT","TA")
ambiguous<-subset(bsp.aut,bsp.aut$gen %in% transitions)
ambiguous<-droplevels(ambiguous)
write.table(as.matrix(ambiguous$V2),file='ambiguous.snplist',quote=F,row.names=F,col.names=F)
#filter variants
system('plink --bfile array_wes_aut --mind 0.05 --geno 0.05 --extract ambiguous.snplist --hwe 0.01 --make-bed --out array_wes_filtered')
```

Overlap BIPMed and 1000 Genome variant list.
```r
# create BIPMed SNP list and extract them from 1000 genome
system('plink --bfile array_wes_filtered --write-snplist --out array_wes_filtered')
# extract variants from 1000 genome dataset
for (i in 1:22){
system(paste0("plink --vcf ALL.chr",i,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --make-bed --out thousand_genome_chr",i,"_extract --biallelic-only strict --allow-extra-chr"))
}
system("ls thousand_genome_chr*log | cut -d'.' -f1 > mergelist.txt")
system("plink --merge-list mergelist.txt --out thousand_genome_raw")
# change population information in 1000 Genome .fam file
thousandgenome.panel<-read.table('integrated_call_samples_v3.20130502.ALL.panel',header=T)
thousandgenome.fam<-read.table('thousand_genome_raw.fam')
thousandgenome.fam$pop<-thousand.genome.ids$pop
write.table(thousandgenome.fam[,c(7,2:6)],file='thousand_genome_raw.fam',quote = F,row.names = F,col.names = F,sep="\t")

system('plink --bfile thousand_genome_raw --extract array_wes_filtered.snplist --make-bed --out thousand_genome_extract')
# extract the 1000 genome SNP list from BIPMed, creating only overlap SNPs to further comparisons
system('plink --bfile thousand_genome_extract --write-snplist --out thousand_genome_extract')
system('plink --bfile array_wes_filtered --extract thousand_genome_extract.snplist --make-bed --out array_wes_extract')
```

Merge datasets
```r
system('plink --bfile array_wes_extract --bmerge thousand_genome_extract --make-bed --out array_wes_thousandgenome')
## flip SNPs
system('plink --bfile array_wes_extract --flip array_wes_thousandgenome-merge.missnp --make-bed --out array_wes_flip')

## merge again
system('plink --bfile array_wes_flip --bmerge thousand_genome_extract --make-bed --out array_wes_thousandgenome')

## exclude SNPs that did not match even after flipping
system('plink --bfile array_wes_flip --exclude array_wes_thousandgenome-merge.missnp --make-bed --out array_wes_exclude')
system('plink --bfile  thousand_genome_extract --exclude array_wes_thousandgenome-merge.missnp --make-bed --out thousand_genome_exclude')

## merge again
system('plink --bfile array_wes_exclude --bmerge thousand_genome_exclude --make-bed --out array_wes_thousandgenome')
```


## Compare allele frequencies between BIPMed and 1000 Genome datasets

Calculate MAF from the data separating by continental population
```r
# create a fam file combining EUR, AFR and AMR populations 
array.wes.thousandgenome<-read.table('array_wes_thousandgenome.fam')
thousandgenome.ids<-read.table('integrated_call_samples_v3.20130502.ALL.panel',header=T)
array.wes.thousandgenome$V1<-recode(array.wes.thousandgenome$V1,"c('ESN','GWD','LWK','MSL','YRI')='AFR';c('CLM','MXL','PEL','PUR')='AMR';c('CDX','CHB','CHS','JPT','KHV')='EAS';c('CEU','GBR','IBS','TSI')='EUR'")
write.table(array.wes.thousandgenome,file='array_wes_thousandgenome_continental.fam',quote = F,row.names = F,col.names = F,sep="\t")

# calculate MAF separating by continental populations
system('plink --bed array_wes_thousandgenome.bed --bim array_wes_thousandgenome.bim --fam array_wes_thousandgenome_continental.fam --freq --family --out array_wes_thousandgenome')

# load frequency estimates file
all.freq<-read.table('/home/secolin/populationstructure/data/array_wes_thousandgenome.frq.strat',header = T)

# separate by population
all.freq.split<-split(all.freq,all.freq$CLST)

# include MAF columns from EUR, AFR and AMR to BRS data.frame
all.freq.split$BRS$MAF.EUR<-all.freq.split$EUR$MAF
all.freq.split$BRS$MAF.AFR<-all.freq.split$AFR$MAF
all.freq.split$BRS$MAF.AMR<-all.freq.split$AMR$MAF
```

### MAF comparison between Europeans and BIPMed
```r
## estimates count
addmargins(table(all.freq.split$BRS$MAF.EUR<0.01,all.freq.split$BRS$MAF<0.01))
## estimates percentage
addmargins(round(prop.table(table(all.freq.split$BRS$MAF.EUR<0.01,all.freq.split$BRS$MAF<0.01)),digits = 3))
```

### MAF comparison between Africans and BIPMed
```r
addmargins(table(all.freq.split$BRS$MAF.AFR<0.01,all.freq.split$BRS$MAF<0.01))
addmargins(round(prop.table(table(all.freq.split$BRS$MAF.AFR<0.01,all.freq.split$BRS$MAF<0.01)),digits = 3))
```

### MAF comparison between admixed Americans and BIPMed
```r
addmargins(table(all.freq.split$BRS$MAF.AMR<0.01,all.freq.split$BRS$MAF<0.01))
addmargins(round(prop.table(table(all.freq.split$BRS$MAF.AMR<0.01,all.freq.split$BRS$MAF<0.01)),digits = 3))
```

# R Package version information
```
R version 3.6.3 (2020-02-29)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 16.04.6 LTS

Matrix products: default
BLAS:   /usr/lib/libblas/libblas.so.3.6.0
LAPACK: /usr/lib/lapack/liblapack.so.3.6.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] ggplot2_3.2.1

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.2       knitr_1.25       magrittr_1.5     tidyselect_0.2.5
 [5] munsell_0.5.0    colorspace_1.4-1 R6_2.4.0         rlang_0.4.0     
 [9] stringr_1.4.0    dplyr_0.8.3      tools_3.6.3      grid_3.6.3      
[13] gtable_0.3.0     xfun_0.10        withr_2.1.2      htmltools_0.4.0 
[17] yaml_2.2.0       lazyeval_0.2.2   digest_0.6.21    assertthat_0.2.1
[21] tibble_2.1.3     crayon_1.3.4     purrr_0.3.3      glue_1.3.1      
[25] evaluate_0.14    rmarkdown_1.16   stringi_1.4.3    compiler_3.6.3  
[29] pillar_1.4.2     scales_1.0.0     pkgconfig_2.0.3
```
