# Labrador retriever GWAS for ACL rupture (air-bud syndrome)

# Article Source: Genome-wide association analysis in dogs implicates 99 
# loci as risk variants for anterior cruciate ligament rupture
# Baker LA, Kirkpatrick B, Rosa GJM, Gianola D, Valente B, et al. (2017) 
# PLOS ONE 12(4): e0173810. https://doi.org/10.1371/journal.pone.0173810

# The data is available on data dryad (curl download below)
# There is better-explained info on Plink here:


### Do this on a login node    #########################################

curl -L http://datadryad.org/api/v2/datasets/doi%253A10.5061%252Fdryad.8kk06/download --output labrador_download
unzip labrador_download
rm labrador_download

# note that we get the .bed, .bim., and .fam files 
# (genotypes, locus info, phenotypes)

#### Start your interactive session now ################################

srun --pty --account def-nricker --mem=4G -N 1 -n 4 -t 0-01:30 /bin/bash

# modules to load on graham
module load nixpkgs/16.09 gcc/7.3.0 r/4.0.2 plink/1.9b_4.1-x86_64 

## 1: Filter by genotype missingness ####################################

# check missingness with --missing
# because we have non-human data, we always specify the organism --dog
plink --bfile cr237_dryad --missing --dog 

# files .imiss .lmiss show the proportion of missing data per loci and indiv.
# plot the missingness results.
Rscript --no-save ./R/hist_miss.R

# Delete SNPs with missingness >0.02.
plink --bfile cr237_dryad --geno 0.02 --make-bed --out cr237_dryad_2

# Delete individuals with missingness >0.02.
plink --bfile cr237_dryad_2 --mind 0.02 --make-bed --out cr237_dryad_3

 # Total genotyping rate is 0.997899.
 # 1430443 variants and 165 people pass filters and QC.
 # Among remaining phenotypes, 56 are cases and 56 are controls.  
 # (53 phenotypes are missing.) --> # (53 phenotypes are missing.)


We dropped 21,579 of 1,457,897 variants due to missing genotype data (--geno) and none from the --mind filter. Note that we dropped some of the data we were getting warnings about:

 > Warning: 179 het. haploid genotypes present (see HapMap_3_r3_5.hh ); many
 > commands treat these as missing.

### Step2: Sex discrepency 

Filtering based on F value which is related to X chromosome homozygosity / inbreeding, and informs us if any individuals are sexed incorrectly or ambiguously. Male > 0.8, F < 0.2 to pass. Uses plink `--check-sex` flag

<!-- # chromosome inbreeding (homozygosity) estimate
<!-- # Check for sex discrepancy. -->
<!-- # Subjects who were a priori determined as females must have a  -->
<!-- # F value of <0.2, and subjects who were a priori determined as -->
<!-- # males must have a F value >0.8. This F value is based on the X -->
<!-- # chromosome inbreeding (homozygosity) estimate. -->
<!-- # Subjects who do not fulfil these requirements are flagged "PROBLEM" -->
<!-- # by PLINK. -->

```{r, engine='bash'}
plink --bfile HapMap_3_r3_3 --check-sex 
grep "PROBLEM" plink.sexcheck
```

 we see that NA10854 (female) is flagged 'Problem' because has F=0.99
 pedsex is 2 but snpsex is 1
 >1349   NA10854     2     1     PROBLEM   0.99

This is flagged for removal later.

```{r, engine='bash'}
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
```


<!-- # These checks indicate that there is one woman with a  -->
<!-- # sex discrepancy, F value of 0.99. (When using other datasets  -->
<!-- # often a few discrepancies will be found).  -->

<!-- # The following two scripts can be used to deal with  -->
<!-- # individuals with a sex discrepancy. -->
<!-- # Note, please use one of the two options below to  -->
<!-- # generate the bfile hapmap_r23a_6, this file we will use  -->
<!-- # in the next step of this tutorial. -->

```{r, engine='bash'}
# 1) Delete individuals with sex discrepancy.
# This command generates a list of individuals with the status ?PROBLEM?.
grep "PROBLEM" plink.sexcheck| awk '{print$1,$2}'> sex_discrepancy.txt
# This command removes the list of individuals with the status ?PROBLEM?.
plink --bfile HapMap_3_r3_3 --remove sex_discrepancy.txt --make-bed --out HapMap_3_r3_4 

# 2) OR  - impute-sex.
# This imputes the sex based on the genotype information into 
# your data set.
plink --bfile HapMap_3_r3_3 --impute-sex --make-bed --out HapMap_3_r3_4
```

<!-- # > 165 people (80 males, 85 females) loaded from .fam. -->
<!-- # > 112 phenotype values loaded from .fam. -->
<!-- # > --remove: 164 people remaining. -->
<!--  or -->

<!-- # > --impute-sex: 23424 Xchr and 0 Ychr variant(s) scanned, all sexes imputed. -->

### Step 3 

Filter .bed file to autosomal SNPs only and delete SNPs with low MAF.
For this GWAS we aren't interested in sex-linked traits.
Generally, we want a threshold of MAF > 0.05, or possibly lower if our sample is large.

```{r, engine='bash'}
# extract all snps that are on chromosomes 1-22, not 25 (mitochondrial)
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' HapMap_3_r3_6.bim > snp_1_22.txt

# Select autosomal SNPs only (i.e., from chromosomes 1 to 22).
plink --bfile HapMap_3_r3_4 --extract snp_1_22.txt --make-bed --out HapMap_3_r3_5

# Generate a plot of the MAF distribution. (MAF_distribution.pdf)
plink --bfile HapMap_3_r3_5 --freq --out MAF_check
Rscript --no-save MAF_check.R

# Remove SNPs with a low MAF frequency. (5%)
plink --bfile HapMap_3_r3_5 --maf 0.05 --make-bed --out HapMap_3_r3_6
```
<!-- # More that a million (1,073,226) SNPs are left -->
<!-- # A conventional MAF threshold for a regular GWAS is between  -->
<!-- # 0.01 or 0.05, depending on sample size. -->

<!-- # > 325318 variants removed due to minor allele threshold(s) -->


### Step 4 ###

Delete SNPs which are not in Hardy-Weinberg equilibrium (HWE) based on some threshold for p-value. Threshold is higher for controls than cases, so we do two steps, where 2nd has an added flag.

Plot the distribution of HWE p-values for SNPs. (plink.hwe output)

<!-- # By default the --hwe option in plink only filters for controls. Therefore, we use two steps, first we use a stringent HWE threshold for controls, followed by a less stringent threshold for the case data. -->
<!--  threshold of p<1-e6 for controls -->
<!--  threshold of p<1-e10 for cases -->

<!-- # Selecting SNPs with HWE p-value below 0.00001, required for one -->
<!-- # of the two plot generated by the next Rscript, allows to zoom in -->
<!-- # on strongly deviating SNPs. -->

<!-- # The HWE threshold for the cases filters out only SNPs which  -->
<!-- # deviate extremely from HWE.  -->
<!-- # This second HWE step only focusses on cases because in the -->
<!-- # controls all SNPs with a HWE p-value < hwe 1e-6 were already -->
<!-- # removed -->
<!-- # Theoretical background for this step is given in our accompanying article: -->
<!-- # https://www.ncbi.nlm.nih.gov/pubmed/29484742 . -->

```{r, engine='bash'}
# compute HWE
plink --bfile HapMap_3_r3_6 --hardy
# output: plink.hwe

# filter to significant HWE violations only for plot
awk '{ if ($9 <0.00001) print $0 }' plink.hwe > plinkzoomhwe.hwe
Rscript --no-save hwe.R

# -hwe filters only Controls (default) with HWE p-value > 1e-6
plink --bfile HapMap_3_r3_6 --hwe 1e-6 --make-bed --out HapMap_hwe_filter_step1
# to filter Case SNPs, we add hwe-all flag
plink --bfile HapMap_hwe_filter_step1 --hwe 1e-10 --hwe-all --make-bed --out HapMap_3_r3_7
```
<!-- --hwe: 0 variants removed due to Hardy-Weinberg exact test. -->
<!-- --hwe: 10 variants removed due to Hardy-Weinberg exact test. -->


<!-- ```{r, engine='bash'} -->
<!-- # ``` -->
<!-- # hwe<-read.table (file="plink.hwe", header=TRUE) -->
<!-- # pdf("histhwe.pdf") -->
<!-- # hist(hwe[,9],main="Histogram HWE") -->
<!-- # dev.off() -->
<!-- #  -->
<!-- # hwe_zoom<-read.table (file="plinkzoomhwe.hwe", header=TRUE) -->
<!-- # pdf("histhwe_below_theshold.pdf") -->
<!-- # hist(hwe_zoom[,9],main="Histogram HWE: strongly deviating SNPs only") -->
<!-- # dev.off() -->
<!-- # ``` -->
<!-- ``` -->


### Step 5: Heterozygosity filter.

Prune SNPs to only consider uncorrelated loci (requires an input file of LD regions to exclude).

<!-- # Generate a plot of the distribution of the heterozygosity rate of your subjects. -->
<!-- # And remove individuals with a heterozygosity rate -->
<!-- # deviating more than 3 sd from the mean. -->
<!-- # Checks for heterozygosity are performed on a set of SNPs -->
<!-- # which are not highly correlated. -->
<!-- # Therefore, to generate a list of non-(highly)correlated SNPs, -->
<!-- # we exclude high inversion regions (inversion.txt [High LD regions]) -->
<!-- # and prune the SNPs using the command --indep-pairwise?. -->
<!-- # The parameters `50 5 0.2` stand respectively for:  -->
<!-- #  1 the window size, -->
<!-- #  2 the number of SNPs to shift the window at each step,  -->
<!-- #  3 and the multiple # correlation coefficient for a SNP being -->
<!-- #    regressed on all other SNPs simultaneously. -->

'Pruning' SNPs ...

Filtering individuals with heterozygosity that deviates from the mean.


```{r, engine='bash', eval = FALSE}
cat inversion.txt
# Note, don't delete the file indepSNP.prune.in, we will use this file in later steps of the tutorial.
plink --bfile HapMap_3_r3_7 --exclude inversion.txt --range --indep-pairwise 50 5 0.2 --out indepSNP
# gives .in files to use for --extract. (and .out to check removed loci)
plink --bfile HapMap_3_r3_7 --extract indepSNP.prune.in --het --out R_check
# This file contains your pruned data set.

# Plot of the heterozygosity rate distribution
Rscript --no-save check_heterozygosity_rate.R

# The following code generates a list of individuals who deviate more than 3 standard deviations from the heterozygosity rate mean.
# For data manipulation we recommend using UNIX. However, when performing statistical calculations R might be more convenient, hence the use of the Rscript for this step:
Rscript --no-save heterozygosity_outliers_list.R
# output: fail-het-qc.txt .

# When using our example data/the HapMap data this list contains 2 individuals (i.e., two individuals have a heterozygosity rate deviating more than 3 SD's from the mean).
# Adapt this file to make it compatible for PLINK, by removing all quotation marks from the file and selecting only the first two columns.
sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt
cat het_fail_ind.txt

# --remove: heterozygosity rate outliers from index file.
plink --bfile HapMap_3_r3_9 --remove het_fail_ind.txt --make-bed --out HapMap_3_r3_10

```


### step 6

<!-- It is essential to check datasets you analyse for cryptic relatedness. -->
<!-- Assuming a random population sample we are going to exclude all individuals above the pihat threshold of 0.2 in this tutorial. -->

- Examine individuals with high relatedness.
- Want to exclude relatives from analysis (indep. sample).
- Keep 'founders'; remove any individuals whose parents are in the study.
- Normally remove any individuals with 'cryptic relatedness'.



```{r, engine='bash', eval = FALSE}
# Check for relationships between individuals with a pihat > 0.2.
plink --bfile HapMap_3_r3_10 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2

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
# < !-- vi 0.2_low_call_rate_pihat.txt --> vi blows
nano 0.2_low_call_rate_pihat.txt 
# put in this line: 13291  NA07045
```

In case of multiple 'related' pairs, the list generated above can be extended using the same method as for our lone 'related' pair.

```{r, engine='bash', eval = FALSE}
# Delete the individuals with the lowest call rate in 'related' pairs with a pihat > 0.2 
plink --bfile HapMap_3_r3_11 --remove 0.2_low_call_rate_pihat.txt --make-bed --out HapMap_3_r3_12
```

The end of filtering steps


# Part 2:



















# My dataset: Poplar GWAS for biofuel characteristics

Full dataset contains 882 Poplar trees, 

### Using Globus to transfer poplar data to computecanada  

- 1 Go to Globus help page on compute canada https://docs.computecanada.ca/wiki/Globus.
- 2 Follow the link to log on to globus. (You use can use a gmail or orcid to log on to globus; or computecanada credentials - but that didn't work for me.

- 3 Once into the system, set up the two panel view (top left buttons)

- 4 In the left panel, put:
  - Collection: OLCF DOI-DOWNLOADS
  - Path: /OLCF/201712/10.13139_OLCF_1411410/

- 5 In the right panel, put (eg. for graham):
  - Collection: computecanada#graham-dtn 
  - Path: < ~/scratch/where_you_want

- Click the transfer arrow on the left to transfer the poplar data to your chosen directory. (Not the one on the right - that does the opposite transfer, which I tried a few times *lol*).

This is all a bit of a pain, but the download is very fast.

### 

- Unzip the files: `tar -vxzf 55.tar.gz `
 
- What is there?

```

