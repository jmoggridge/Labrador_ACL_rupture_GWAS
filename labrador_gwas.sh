# Labrador retriever GWAS for ACL rupture (air-bud syndrome)

# Article Source: Genome-wide association analysis in dogs implicates 99 
# loci as risk variants for anterior cruciate ligament rupture
# Baker LA, Kirkpatrick B, Rosa GJM, Gianola D, Valente B, et al. (2017) 
# PLOS ONE 12(4): e0173810. https://doi.org/10.1371/journal.pone.0173810

# The data is available on data dryad (curl download below)
# There is better-explained info on Plink here:

# use scp to grab pdf files in you terminal and check them as you go!

### Do this on a login node    #########################################

curl -L http://datadryad.org/api/v2/datasets/doi%253A10.5061%252Fdryad.8kk06/download --output labrador_download
unzip labrador_download
rm labrador_download

# note that we get the cr237_dryad.bed, .bim., and .fam files 
# .bed is a binary file with genotypes
# .bim is a binary file with info about variant loci
# .fam has the phenotype info and genders and is text

# The .fam file has no header, so it's useful to know the columns:
# 1 Family ID ('FID')
# 2 Within-family ID ('IID'; cannot be '0')
# 3 Within-family ID of father ('0' if father isn't in dataset)
# 4 Within-family ID of mother ('0' if mother isn't in dataset)
# 5 Sex code ('1' = male, '2' = female, '0' = unknown)
# 6 Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)

head cr237_dryad.fam

# We don't have any relatives so the FIDs are all the same as the IIDs
# and columns 3 and 4 are all zeros. Our outcome of interest is column 6.

# see here for more info about file formats (other input formats are possible):
# https://www.cog-genomics.org/plink/1.9/formats

### Getting plinky with it: Basic Plink usage primer

# plink --(b)file {my_data} --dog --options --commands --out {new_name} 

# we use the --bfile flag because we have bed/bim/fam files
# we use the --make-bed flag to get the same type of files back
# we always have to use the --dog flag because 38 chromosomes


#### Start your interactive session now ################################

srun --pty --account def-nricker --mem=4G -N 1 -n 4 -t 0-01:30 /bin/bash

# modules to load on graham
module load nixpkgs/16.09 gcc/7.3.0 r/4.0.2 plink/1.9b_4.1-x86_64 

### Filtering done in Baker et al.

# "Genome-wide SNP genotyping was performed in 98 cases and 139 controls
# using the Illumina CanineHD BeadChip, which genotypes 173,662 SNPs evenly
# spaced across the genome. Data underwent quality control filtering using
# PLINK [32]. All samples had a genotyping call rate of ≥95%. 49,859 SNPs
# were excluded because minor allele frequency (MAF) was ≤0.05 and 7,468 
# SNPs were excluded because of a low genotyping rate (≤95%). 153 SNPs were
# excluded because of deviation from Hardy-Weinberg equilibrium at P<1E-07.
# 118,992 SNPs were used for further analysis."

#################################################################
## PART 1 - QC & FILTERING SNPS AND SAMPLES
#################################################################


## 1: Data Missingness ############################

# check missingness with --missing
# because we have non-human data, we always specify the organism --dog
plink --bfile cr237_dryad --missing --dog 
# other you'll get an error like this:
# "Error: Invalid chromosome code '27' on line 128594 of .bim file."

# files .imiss .lmiss show the proportion of missing data per loci and indiv.
# plot the missingness results.
Rscript --no-save ./R/hist_miss.R

# Delete SNPs and individuals with missingness >0.05.
plink --bfile cr237_dryad --geno 0.05 --dog --make-bed --out cr237_dryad_2
plink --bfile cr237_dryad_2 --mind 0.05 --dog --make-bed --out cr237_dryad_3

# you can see in the plink output how many variants were dropped:
# 7468 variants removed due to missing genotype data (--geno).
# 0 dogs removed due to missing genotype data (--mind).

## 2: Sex Discrepency  ##########################################

# ignore snps in the canine pseudoautosomal region of X-chromosome
plink --bfile cr237_dryad_3 --dog --split-x 6630000 126883977 \
  --make-bed --out cr237_dryad_3b
# 568 chromosome codes changed; this is consistent with the 565 in the PAR here:
# https://www.genetics.org/content/184/2/595#T3
# 25567 het. haploid genotypes decreases to 5284 het. haploid genotypes

# Checking the F values, which are related to X chromosome homozygosity
# / inbreeding,  and informs us if any individuals are sexed incorrectly
# or ambiguously. Default to pass is Male F > 0.8, Female F < 0.2. 
# Males should clump near 1, females more widely around 0.

plink --bfile cr237_dryad_3b --check-sex 0.5 0.7 --dog
grep "PROBLEM" plink.sexcheck| awk '{print$1,$2}'> sex_discrepancy.txt
wc -l sex_discrepancy.txt 
cat sex_discrepancy.txt

# Generate plots to visualize the sex-check results.
Rscript ./R/gender_check.R 
# Spits out Gender_check.pdf

# We see a number of discrepencies, but the authors didn't filter these out
# according to their methods section.
# We could impute them from genotypes, and see the changes
plink --bfile cr237_dryad_3b --impute-sex 0.5 0.7 --dog --make-bed \
  --out cr237_dryad_3b_imputed
plink --bfile cr237_dryad_3b_imputed --check-sex 0.5 0.7 --dog
grep "PROBLEM" plink.sexcheck

## This command removes individuals with PROBLEM status
# plink --bfile cr237_dryad_3 --remove sex_discrepancy.txt --dog --make-bed --out cr237_dryad_4 

# To be honest, I think it's a bit beyond this tutorial to worry
# about this ambiguities. I hope we can just trust the original data
# but if we encounter something weird, it might be worth looking into deeper.
rm plink.sex_check cr237_dryad_3b_imputed cr237_dryad_3b

#### 3: Autosomal SNPs only   ##################################

# Filter .bed file to autosomal SNPs only and delete SNPs with low MAF.
# For this GWAS we aren't interested in sex-linked traits.

# extract only snps that are on chromosomes 1-38 (autosomal only)
awk '{ if ($1 >= 1 && $1 <= 38) print $2 }' cr237_dryad_3.bim > snp_1_38.txt
head snp_1_38.txt 
plink --bfile cr237_dryad_3 --extract snp_1_38.txt --dog --make-bed \
  --out cr237_dryad_4

# 160929 variants and 237 dogs pass filters and QC.


#### 4: MAF   ###################################################

# Generally, we want a threshold of MAF > 0.05, 
# or possibly lower if our sample is very large.
# Generate a plot of the MAF distribution. (MAF_distribution.pdf)
plink --bfile cr237_dryad_4 --dog --freq --out MAF_check
Rscript --no-save ./R/MAF_check.R

# Remove SNPs with MAF frequency <5%
plink --bfile cr237_dryad_4 --maf 0.05 --dog --make-bed \
  --out cr237_dryad_5
# 45090 variants removed due to minor allele threshold(s)
# 115839 variants and 237 dogs pass filters and QC. (similar to Baker et al.)


#### 5: Hardy-Weinberg Equilibrium ############################################

# Far more heterozygous calls than expected? Probably systemic calling error.
# Fewer heterozygotes than expected? Could be population stratification.
# (Migration, natural selection, and/or non-random mating violations).
# GWAS methods often assume HW-equilibrium.

# We want to delete SNPs which are far from Hardy-Weinberg equilibrium (HWE)
# based on some threshold for p-value (eg. 10e-7 in Baker et al.).
# The --hwe default in plink only filters the controls, so we use two steps, 
# first on controls, then cases. The thresholds can be looser on cases 
# than controls (more significant).

# compute HWE
plink --bfile cr237_dryad_5 --hardy --dog
head plink.hwe

# filter to significant HWE violations only for plot
awk '{ if ($9 <0.00001) print $0 }' plink.hwe > plinkzoomhwe.hwe
head plinkzoomhwe.hwe

# Plot the distribution of HWE p-values across SNPs from plink.hwe output
Rscript --no-save ./R/hwe.R

# --hwe filter for controls with HWE p-value > 1e-7
plink --bfile cr237_dryad_5 --hwe 1e-7 --dog  \
  --make-bed --out cr237_hwe_filter_step1
# --hwe: 43 variants removed due to Hardy-Weinberg exact test.

# Do --hwe 1e-7 and --hwe-all to filter the case SNPs
plink --bfile cr237_hwe_filter_step1 --hwe 1e-7 --dog \
  --hwe-all --make-bed --out cr237_dryad_6
# --hwe: 93 variants removed due to Hardy-Weinberg exact test.
# 115703 variants and 237 dogs pass filters and QC.


### 5: Linkage Disequilibrium & Pruning ###########################

# 'We will 'prune' the set of SNPs with `--indep-pairwise`
# to get a subset consisting of loci that are
# in approximate linkage equilibrium. This removes a bunch of SNPs
# that are highly correlated.
# The --indep-pairwise parameters `50 5 0.2` stand respectively for:  
#  1 the window size 50 kbp, 
#  2 the number of SNPs to shift the window at each step,  
#  3 and the multiple correlation coefficient for a SNP being 
#    regressed on all other SNPs simultaneously. 
# Any pairs of SNPs that are too correlated will have one removed
# Then we filter any individuals with heterozygosity >3sd from mean 
# across the pruned set of SNPs.

# Get a list of SNPs with high correlation to 'Prune' 
plink --bfile cr237_dryad_6 --dog \
  --indep-pairwise 50kb 5 0.2 --out indepSNP

# Subset of pruned SNPs with --extract
# Check the heterozygosity of these snps with --het
plink --bfile cr237_dryad_6 --dog \
  --extract indepSNP.prune.in --het --out R_check

# Plot of the heterozygosity rate distribution
Rscript --no-save ./R/check_heterozygosity_rate.R
head R_check.het

# Create a list of the outliers (>3sd from mean); remove quot marks
Rscript --no-save ./R/heterozygosity_outliers_list.R
sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt
cat het_fail_ind.txt

# --remove heterozygosity rate outliers using that list
plink --bfile cr237_dryad_6 --dog \
  --remove het_fail_ind.txt \
  --make-bed --out cr237_dryad_7


### 6: Cryptic Relatedness Test ##########################################

# Want to exclude relatives from analysis (independant sample assumed).
# Examine data for individuals with high relatedness.
# Normally remove any individuals with 'cryptic relatedness',
# where `pi_hat` is greater than some threshold

# Here we use the --genome flag to investigate relation between
# pairs of individuals (identity by descent).
# The --min flag sets minimum PI_HAT for inclusion, where PI_HAT
# 1 is identical, 0.5 is parent-child or siblings, 0.25 is aunt-neice.
# Our choice 0.2 is somewhere between double first cousins and first 
# cousins.

# Check for relationships between individuals with a pihat > 0.2.
plink --bfile cr237_dryad_7 --dog \
  --extract indepSNP.prune.in --genome --min 0.2 \
  --out pihat_min0.2
head pihat_min0.2.genome

# visualize parent-offspring relations using the z values. 
awk '{ if ($8 >0.9) print $0 }' pihat_min0.2.genome > zoom_pihat.genome
cat zoom_pihat.genome
# these could be dog siblings based on the pi_hat value and Z1 value

# Generate a plot to assess the type of relationship.
Rscript --no-save ./R/relatedness_tidy.R
# (explentation plot; PO = parent-offspring, UN = unrelated individuals) 

# The file 'pihat_min0.2_in_founders.genome' shows that,
# 82 pairs of dogs have pihat >0.2

# For each pair of 'related' individuals with a pihat > 0.2,
# we recommend to remove the individual with the lowest call rate. 
plink --bfile cr237_dryad_8 --dog --missing
head plink.imiss 

## ok so we could venture to identify which to remove based on 
## missingness but it doesn't seem that the authors did this,
## but rather made up for it in their modelling 

# Generally, drop the individual with the lower call rate in the 'related pair'. 

# Generate a list of FID and IID of the individual(s) with
# a Pihat above 0.2, to check who had the lower call rate of the pair.
# In our dataset the individual 13291  NA07045 had the lower call rate.
# nano 0.2_low_call_rate_pihat.txt 
# eg. put a line for each individual to remove with the id columsn
# : 13291  NA07045


## In case of multiple 'related' pairs, the list generated 
## above can be extended using the same method as for our lone
## 'related' pair.
## Delete the individuals with the lowest call rate in 'related' 
## pairs with a pihat > 0.2 

# plink --bfile cr237_dryad_8 --remove 0.2_low_call_rate_pihat.txt --make-bed --out HapMap_3_r3_12


# The end of QC filtering steps


# download all the plots and data of interest 
# tidy up workspace
for i in 2 3 4 5 6 7; do rm -i *dryad_$i.*; done
rm *_dryad_3b*
rm *_hwe_filter*




#################################################################
# PART 2: DEALING WITH POP STRATIFICATION
#################################################################

## If we have resources to place our study subjects within their 
## 'ethnic' groups (eg. 1000 genomes SNPs data for human ethnicity),
## then we can 'anchor' our samples with this data to find out where
## the samples lie in relation to these groups. This uses MDS to
## show distance between individuals, such that we can identify any
## individuals that are outliers from their groups. Because we want 
## to do our study on a homogenous group, we remove any samples 
## that deviate from the group based on a selected threshold.
## Since we are dealing with dogs, (and there is no 1000 dog genomes)
## we'll just perform MDS to check for any outliers in the labrador
## retrievers. We also take the MDS components as covariates to control
## for any stratification in our sample

## For more details, see
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/, which
## provides a worked example with human SNP data anchored to the 
## 1000 genomes data.

## For more info on plink clustering options: 
# https://www.cog-genomics.org/plink/1.9/strat

# Get the pairwise IBS metrics
plink --bfile cr237_dryad_8 --dog --extract indepSNP.prune.in \
  --genome --out cr237_dryad_8
head cr237_dryad_8.genome

# Do MDS clustering of individuals based on autosomal SNP genotypes 
# We pass the IBS data along with --read-genome.
# --cluster does complete linkage hclust by default,
# so we pass the --mds-plot flag to get MDS and tell it to give
# the first 10 dimensions

plink --bfile cr237_dryad_8 --dog --read-genome cr237_dryad_8.genome \
  --cluster --mds-plot 10 --out cr237_dryad_8

# download these two files into your local project directory:
# cr237_dryad_8.mds
# cr237_dryad_8.fam
# plot the MDS in Rstudio by running: source("./R/MDS_presentation.R")

# In the first two MDS components we can vaguely see two groups,
# but no individuals are extremely far from the cloud.

# We dont need to, but we could exclude the couple controls that
# are around -0.1 in coefficient 2 as an example

# Using awk to exlcude them from the list
awk '{ if ($5 >-0.07) print $1,$2 }' cr237_dryad_8.mds \
  > cr237_dryad_8_strat

# To only keep the individuals we want
plink --bfile cr237_dryad_8 --dog --keep cr237_dryad_8_strat \
  --make-bed --out cr237_dryad_8_strat

# You can see that this removed those three distant dogs.
# However we don't need to do this because we don't really have 
# any other data to reference against here.
# Removing individuals decreases our power (and otherwise
# we'd have to redo the MDS).
rm cr237_dryad_8_stra*

# Speaking of MDS, we want to reformat the .mds file to a
# plink covariate file for the association analysis.
# This will help to control for any population stratification.

awk '{print$1, $2, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' \
  cr237_dryad_8.mds > covar_mds.txt
head covar_mds.txt 


#################################################################
## PART 3 - Associations
#################################################################

## Association testing for binary traits
# (--assoc) without correction for MDS covariates (bad)
plink --bfile cr237_dryad_8 --dog \
  --extract indepSNP.prune.in \
  --assoc --out result1
head result1.assoc 


# Logistic regression (--logistic) for binary trait with 
# 10 principal components. We include MDS components by providing
# our list in covar_mds.txt to the `--covar` parameter
# from. Option --hide-covar shows only additive results (not the
# tests against all the covariates)
# --extract is providing the pruned SNP subset indices
plink --bfile cr237_dryad_8 --dog \
  --extract indepSNP.prune.in \
  --covar covar_mds.txt \
  --logistic \
  --hide-covar \
  --out result2
head result2.assoc.logistic 
wc -l result2.assoc.logistic 

# Check if there are any NA values in the results
awk '/'NA'/' result2.assoc.logistic | wc -l 
# Remove any NA values from results
awk '!/'NA'/' result2.assoc.logistic  > result2.assoc.logistic




