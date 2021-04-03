# GWA tutorial script 1

## figure out what this means:
# Before main variant filters, 112 founders and 53 nonfounders present.


# first get the repository for the tutorial:
git clone https://github.com/MareesAT/GWA_tutorial.git
cd /scratch/jmoggrid/project4_gwa/GWA_tutorial/


# now, setup an interactive session
# srun --pty --account {my-account-name} –mem=16G -N 1 -n 6 -t 0-01:30 /bin/bash
srun --pty --account def-nricker --mem=16G -N 1 -n 6 -t 0-01:30 /bin/bash

# load modules plink and r (with dependencies on graham)
module load nixpkgs/16.09 gcc/7.3.0 r/4.0.2 plink/1.9b_4.1-x86_64 

# unzip first tutorial
unzip 1_QC_GWAS.zip 
cd 1_QC_GWAS

## these are the data files that we start with; 
## they are linked, ie. they have data for the same individuals

# HapMap_3_r3_1.bed # binary something file
# HapMap_3_r3_1.bim # this has our snpnames and genotypes
# HapMap_3_r3_1.fam # has relationship info


# Investigate missingness per individual and per SNP and make histograms.
plink --bfile HapMap_3_r3_1 --missing
# output: plink.imiss and plink.lmiss, these files show respectively 
# the proportion of missing SNPs per individual and the proportion of 
# missing individuals per SNP.

# useful lines in the output:
# > 165 people (80 males, 85 females) loaded from .fam.
# > 112 phenotype values loaded from .fam.
# > Total genotyping rate is 0.997378.


# Generate plots to visualize the missingness results.
Rscript --no-save hist_miss.R

# ``` R code for: histimiss.pdf  histlmiss.pdf
# indmiss <- read.table(file = "plink.imiss", header = TRUE)
# snpmiss<-read.table(file = "plink.lmiss", header = TRUE)
# # read data into R 
# 

## histogram of missingness for individuals
# pdf("histimiss.pdf") #indicates pdf format and gives title to file
# hist(indmiss[,6], main="Histogram individual missingness") 

## histogram of missingness for loci 
# pdf("histlmiss.pdf") 
# hist(snpmiss[,5],main="Histogram SNP missingness")  
# dev.off() 
# ```

### Step 1: Missing data? ###

# Delete SNPs and individuals with high levels of
# missingness, explanation of this and all following steps 
# can be found in box 1 and table 1 of the article mentioned 
# in the comments of this script.
# The following two QC commands will not remove any SNPs or 
# individuals. However, it is good practice to start the QC with
# these non-stringent thresholds.  

# Create new bed file after each step (HapMap_3_r3_2,...,_5)
# Delete SNPs with missingness >0.2. 
plink --bfile HapMap_3_r3_1 --geno 0.2 --make-bed --out HapMap_3_r3_2
# note --geno 0.2 doesn't actually filter anything

# Delete individuals with missingness >0.2.
plink --bfile HapMap_3_r3_2 --mind 0.2 --make-bed --out HapMap_3_r3_3
# you'd have to drop --mind 0.01 to lose more than 1 individual

# Delete SNPs with missingness >0.02.
plink --bfile HapMap_3_r3_3 --geno 0.02 --make-bed --out HapMap_3_r3_4
# > 21579 variants removed due to missing genotype data (--geno).
# note that we dropped some of the data we were getting warnings about
# > Warning: 179 het. haploid genotypes present (see HapMap_3_r3_5.hh ); many
# > commands treat these as missing.

# Delete individuals with missingness >0.02.
plink --bfile HapMap_3_r3_4 --mind 0.02 --make-bed --out HapMap_3_r3_5

# > Total genotyping rate is 0.997899.
# > 1430443 variants and 165 people pass filters and QC.
# > Among remaining phenotypes, 56 are cases and 56 are controls.  
# > (53 phenotypes are missing.)


### Step2: sex discrepency ####

# Check for sex discrepancy.
# Subjects who were a priori determined as females must have a 
# F value of <0.2, and subjects who were a priori determined as
# males must have a F value >0.8. This F value is based on the X
# chromosome inbreeding (homozygosity) estimate.
# Subjects who do not fulfil these requirements are flagged "PROBLEM"
# by PLINK.

plink --bfile HapMap_3_r3_5 --check-sex 

grep "PROBLEM" plink.sexcheck
# we see that NA10854 (female) is flagged 'Problem' because has F=0.99
# pedsex is 2 but snpsex is 1
# >1349   NA10854     2     1     PROBLEM   0.99

# Generate plots to visualize the sex-check results.
Rscript gender_check.R
# > Gender_check.pdf Men_check.pdf Women_check.pdf
# ```
# gender <- read.table("plink.sexcheck", header=T,as.is=T)
# 
# pdf("Gender_check.pdf")
# hist(gender[,6], main = "Gender", xlab = "F")
# dev.off()
# 
# pdf("Men_check.pdf")
# male = subset(gender, gender$PEDSEX==1)
# hist(male[, 6], main = "Men", xlab = "F")
# dev.off()
# 
# pdf("Women_check.pdf")
# female = subset(gender, gender$PEDSEX == 2)
# hist(female[, 6] ,main = "Women", xlab = "F")
# dev.off()
# ```


# These checks indicate that there is one woman with a 
# sex discrepancy, F value of 0.99. (When using other datasets 
# often a few discrepancies will be found). 

# The following two scripts can be used to deal with 
# individuals with a sex discrepancy.
# Note, please use one of the two options below to 
# generate the bfile hapmap_r23a_6, this file we will use 
# in the next step of this tutorial.

# 1) Delete individuals with sex discrepancy.
grep "PROBLEM" plink.sexcheck| awk '{print$1,$2}'> sex_discrepancy.txt

# This command generates a list of individuals with the status ?PROBLEM?.
plink --bfile HapMap_3_r3_5 --remove sex_discrepancy.txt --make-bed --out HapMap_3_r3_6 
# > 165 people (80 males, 85 females) loaded from .fam.
# > 112 phenotype values loaded from .fam.
# > --remove: 164 people remaining.


# This command removes the list of individuals with the status ?PROBLEM?.

# 2) impute-sex.
plink --bfile HapMap_3_r3_5 --impute-sex --make-bed --out HapMap_3_r3_6
# This imputes the sex based on the genotype information into 
# your data set.
# > --impute-sex: 23424 Xchr and 0 Ychr variant(s) scanned, all sexes imputed.




### Step 3 ### 


# Generate a bfile with autosomal SNPs only and delete SNPs 
# with a low minor allele frequency (MAF).

# extract all snps that are on chromosomes 1-22, not 25 (mitochondrial)
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' HapMap_3_r3_6.bim > snp_1_22.txt

# Select autosomal SNPs only (i.e., from chromosomes 1 to 22).
plink --bfile HapMap_3_r3_6 --extract snp_1_22.txt --make-bed --out HapMap_3_r3_7
# > 1430443 variants loaded from .bim file.
# > --extract: 1398544 variants remaining.

# Generate a plot of the MAF distribution. (MAF_distribution.pdf)
plink --bfile HapMap_3_r3_7 --freq --out MAF_check
# > --freq: Allele frequencies (founders only) written to MAF_check.frq .

Rscript --no-save MAF_check.R

# Remove SNPs with a low MAF frequency. (5%)
plink --bfile HapMap_3_r3_7 --maf 0.05 --make-bed --out HapMap_3_r3_8
# > 325318 variants removed due to minor allele threshold(s)

# More that a million (1,073,226) SNPs are left
# A conventional MAF threshold for a regular GWAS is between 
# 0.01 or 0.05, depending on sample size.



### Step 4 ###

# Delete SNPs which are not in Hardy-Weinberg equilibrium (HWE).
# Check the distribution of HWE p-values of all SNPs. (plink.hwe output)
plink --bfile HapMap_3_r3_8 --hardy
# less plink.hwe

# Selecting SNPs with HWE p-value below 0.00001, required for one
# of the two plot generated by the next Rscript, allows to zoom in
# on strongly deviating SNPs.
awk '{ if ($9 <0.00001) print $0 }' plink.hwe > plinkzoomhwe.hwe
# > CHR  SNP     TEST A1  A2 GENO      O(HET) E(HET)    P 
# > 3 rs7623291   ALL  T  C  22/28/62  0.25   0.4362    8.938e-06

Rscript --no-save hwe.R
# ```
# hwe<-read.table (file="plink.hwe", header=TRUE)
# pdf("histhwe.pdf")
# hist(hwe[,9],main="Histogram HWE")
# dev.off()
# 
# hwe_zoom<-read.table (file="plinkzoomhwe.hwe", header=TRUE)
# pdf("histhwe_below_theshold.pdf")
# hist(hwe_zoom[,9],main="Histogram HWE: strongly deviating SNPs only")
# dev.off()
# ```

# By default the --hwe option in plink only filters for controls.
# Therefore, we use two steps, first we use a stringent HWE threshold for controls,
# followed by a less stringent threshold for the case data.
plink --bfile HapMap_3_r3_8 --hwe 1e-8 --make-bed --out HapMap_hwe_filter_step1
# > --hwe: 0 variants removed due to Hardy-Weinberg exact test.

# -------- UP TO HERE ----------- ~~~~~~

# The HWE threshold for the cases filters out only SNPs which 
# deviate extremely from HWE. 
# This second HWE step only focusses on cases because in the
# controls all SNPs with a HWE p-value < hwe 1e-6 were already
# removed
plink --bfile HapMap_hwe_filter_step1 --hwe 1e-10 --hwe-all --make-bed --out HapMap_3_r3_9

# Theoretical background for this step is given in our accompanying article:
# https://www.ncbi.nlm.nih.gov/pubmed/29484742 .

############################################################
### step 5 ###

# Generate a plot of the distribution of the heterozygosity rate of your subjects.
# And remove individuals with a heterozygosity rate
# deviating more than 3 sd from the mean.

# Checks for heterozygosity are performed on a set of SNPs
# which are not highly correlated.
# Therefore, to generate a list of non-(highly)correlated SNPs,
# we exclude high inversion regions (inversion.txt [High LD regions])
# and prune the SNPs using the command --indep-pairwise?.
# The parameters `50 5 0.2` stand respectively for: 
#  1 the window size,
#  2 the number of SNPs to shift the window at each step, 
#  3 and the multiple # correlation coefficient for a SNP being
#    regressed on all other SNPs simultaneously.

plink --bfile HapMap_3_r3_9 --exclude inversion.txt --range --indep-pairwise 50 5 0.2 --out indepSNP
# Note, don't delete the file indepSNP.prune.in, we will use this file in later steps of the tutorial.

plink --bfile HapMap_3_r3_9 --extract indepSNP.prune.in --het --out R_check
# This file contains your pruned data set.

# Plot of the heterozygosity rate distribution
Rscript --no-save check_heterozygosity_rate.R

# The following code generates a list of individuals who deviate more than 3 standard deviations from the heterozygosity rate mean.
# For data manipulation we recommend using UNIX. However, when performing statistical calculations R might be more convenient, hence the use of the Rscript for this step:
Rscript --no-save heterozygosity_outliers_list.R

# Output of the command above: fail-het-qc.txt .
# When using our example data/the HapMap data this list contains 2 individuals (i.e., two individuals have a heterozygosity rate deviating more than 3 SD's from the mean).
# Adapt this file to make it compatible for PLINK, by removing all quotation marks from the file and selecting only the first two columns.
sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt

# Remove heterozygosity rate outliers.
plink --bfile HapMap_3_r3_9 --remove het_fail_ind.txt --make-bed --out HapMap_3_r3_10


############################################################
### step 6 ###

# It is essential to check datasets you analyse for cryptic relatedness.
# Assuming a random population sample we are going to exclude all individuals
# above the pihat threshold of 0.2 in this tutorial.

# Check for relationships between individuals with a pihat > 0.2.
plink --bfile HapMap_3_r3_10 --extract indepSNP.prune.in --genome \
  --min 0.2 --out pihat_min0.2

# The HapMap dataset is known to contain parent-offspring relations. 
# The following commands will visualize specifically these parent-offspring relations, using the z values. 
awk '{ if ($8 >0.9) print $0 }' pihat_min0.2.genome>zoom_pihat.genome

# Generate a plot to assess the type of relationship.
Rscript --no-save Relatedness.R

# The generated plots show a considerable amount of related individuals (explentation plot; PO = parent-offspring, UN = unrelated individuals) in the Hapmap data, this is expected since the dataset was constructed as such.
# Normally, family based data should be analyzed using specific family based methods. In this tutorial, for demonstrative purposes, we treat the relatedness as cryptic relatedness in a random population sample.
# In this tutorial, we aim to remove all 'relatedness' from our dataset.
# To demonstrate that the majority of the relatedness was due to parent-offspring we only include founders (individuals without parents in the dataset).

plink --bfile HapMap_3_r3_10 --filter-founders --make-bed --out HapMap_3_r3_11

# Now we will look again for individuals with a pihat >0.2.
plink --bfile HapMap_3_r3_11 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2_in_founders
# The file 'pihat_min0.2_in_founders.genome' shows that, after exclusion of all non-founders, only 1 individual pair with a pihat greater than 0.2 remains in the HapMap data.
# This is likely to be a full sib or DZ twin pair based on the Z values. Noteworthy, they were not given the same family identity (FID) in the HapMap data.

# For each pair of 'related' individuals with a pihat > 0.2, we recommend to remove the individual with the lowest call rate. 
plink --bfile HapMap_3_r3_11 --missing
# Use an UNIX text editor (e.g., vi(m) ) to check which individual has the highest call rate in the 'related pair'. 

# Generate a list of FID and IID of the individual(s) with a Pihat above 0.2, to check who had the lower call rate of the pair.
# In our dataset the individual 13291  NA07045 had the lower call rate.
vi 0.2_low_call_rate_pihat.txt
i 
13291  NA07045
# Press esc on keyboard!
:x
# Press enter on keyboard
# In case of multiple 'related' pairs, the list generated above can be extended using the same method as for our lone 'related' pair.

# Delete the individuals with the lowest call rate in 'related' pairs with a pihat > 0.2 
plink --bfile HapMap_3_r3_11 --remove 0.2_low_call_rate_pihat.txt --make-bed --out HapMap_3_r3_12

################################################################################################################################

# CONGRATULATIONS!! You've just succesfully completed the first tutorial! You are now able to conduct a proper genetic QC. 

# For the next tutorial, using the script: 2_Main_script_MDS.txt, you need the following files:
# - The bfile HapMap_3_r3_12 (i.e., HapMap_3_r3_12.fam,HapMap_3_r3_12.bed, and HapMap_3_r3_12.bim
# - indepSNP.prune.in







# ------
# 1000 genomes data
# ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp

# Hapmap data
# https://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2010-05_phaseIII/plink_format/

# family relationship metadata
# https://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2010-05_phaseIII/relationships_w_pops_041510.txt


