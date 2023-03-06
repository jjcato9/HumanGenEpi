---
editor_options: 
  markdown: 
    wrap: 72
---

# Polygenic risk score manual

In this manual you will learn how to compute polygenic risck scores
using **PRSice** (<https://www.prsice.info/>)

**You will need to download the following files for this practical**

``` bash
# In your home directory create a folder called PRS and download the files into it

mkdir PRS
cd PRS

wget https://github.com/WCSCourses/HumanGenEpi/raw/main/manuals/Polygenic_risk_scores/ASN.bed.xz

wget https://github.com/WCSCourses/HumanGenEpi/raw/main/manuals/Polygenic_risk_scores/ASN.pheno

wget https://github.com/WCSCourses/HumanGenEpi/raw/main/manuals/Polygenic_risk_scores/ASN.bim

wget https://github.com/WCSCourses/HumanGenEpi/raw/main/manuals/Polygenic_risk_scores/test

wget https://github.com/WCSCourses/HumanGenEpi/raw/main/manuals/Polygenic_risk_scores/validate

wget https://github.com/WCSCourses/HumanGenEpi/raw/main/manuals/Polygenic_risk_scores/ASN.gwas.txt.xz

wget https://github.com/WCSCourses/HumanGenEpi/raw/main/manuals/Polygenic_risk_scores/ASN.cov

wget https://github.com/WCSCourses/HumanGenEpi/raw/main/manuals/Polygenic_risk_scores/ASN.fam

wget https://github.com/WCSCourses/HumanGenEpi/raw/main/manuals/Polygenic_risk_scores/PRSice_linux

wget https://github.com/WCSCourses/HumanGenEpi/raw/main/manuals/Polygenic_risk_scores/PRSice.R

#Unzip the following files as shown below

xz -dv ASN.gwas.txt.xz
xz -dv ASN.bed.xz
```

## 1. Training the PRS to find the best predictive one using the validation target data set

***Its key to do some QC on your discovery data set***

``` bash
Rscript PRSice.R --prsice PRSice_mac --base ASN.gwas.txt --target ASN --binary F --keep validate --pheno ASN.pheno --cov ASN.cov --out trial1
```

Check out the trial1.log file which contains the errors shown above. So,
we see that initially the all the SNPs in the base file are read.
Ambiguous variants are removed to avoid strand errors, however, other
strand flips are automatically detected and accounted for. You can
filter the variants to include based on the quality of imputation and if
duplicate variants are available the latest version of PRSice will stop
as shown above. You can use the --extract trial1.valid however this
option will not show you the detailed breakdown of SNPs included so we
will remove these duplicates using the code below.

``` bash
uniq ASN.gwas.txt > ASN.gwasqc.txt
```

[**This has removed the duplicate files, so we only have unique SNP
files and no repeat files which would inappropriately bias the
data.**]{.underline}

***Preliminary analysis using default settings***

``` bash
Rscript PRSice.R --prsice PRSice_mac --base ASN.gwasqc.txt --target ASN --keep validate --pheno ASN.pheno --binary F --cov ASN.cov --out Prelim
```

Let's look at the .summary file and the plots and ensure you understand
them. What is the PRS R2 and how many SNPs are in the best preforming
PRS ?

[**Look at the .prsice file --\> Look at the highest threshold and you
can see only around 1400 SNPs are performing here using
GWAS.**]{.underline}

[**Looking at summary file, threshold is 0.456???**]{.underline}

[**Look at .log file --\> there are initially 529502 variants for 503
people (240 males, 263 females).**]{.underline}

[**After clumping of the vaitents, 131186 remain. With this, there are
only 240 samples included in the final linear regression (i.e. these are
the only ones with complete datasets, which are needed for linear
regression).**]{.underline}

***Optimise computation to get the most predictive PRS***

-   clump-kb 500 clump-r2 0.1

``` bash
Rscript PRSice.R --prsice PRSice_mac --base ASN.gwasqc.txt --target ASN --binary F --keep validate --pheno ASN.pheno --cov ASN.cov --clump-kb 500 --clump-r2 0.1 --base-info INFO:0.4 --out Opt500_0.1
```

-   clump-kb 250 clump-r2 0.3

``` bash
Rscript PRSice.R --prsice PRSice_mac --base ASN.gwasqc.txt --binary F --target ASN --keep validate --pheno ASN.pheno --cov ASN.cov --clump-kb 250 --clump-r2 0.3 --base-info INFO:0.4 --out Opt250_0.5
```

***Which parametes give the best predictive PRS ? Lets print out the
SNPs of the best predictive PRS***

We need to compare the variance...

``` bash
Rscript PRSice.R --prsice PRSice_mac --base ASN.gwasqc.txt --target ASN --keep validate --pheno ASN.pheno --binary F --cov ASN.cov --clump-kb 500 --clump-r2 0.1 --base-info INFO:0.4 --print-snp --out validation
```

```{bash}
(base) Jessicas-MacBook-Pro-7:Polygenic_risk_scores JJCaterson$ head Opt500_0.1.summary
Phenotype	Set	Threshold	PRS.R2	Full.R2	Null.R2	Prevalence	Coefficient	Standard.ErrorNum_SNP
-	Base	0.46875	0.275608	0.422932	0.203377	-	67987.8	7159.75	2.31876e-18	133903
(base) Jessicas-MacBook-Pro-7:Polygenic_risk_scores JJCaterson$ head Opt250_0.5.summary
Phenotype	Set	Threshold	PRS.R2	Full.R2	Null.R2	Prevalence	Coefficient	Standard.ErrorNum_SNP
-	Base	1	0.221685	0.379976	0.203377	-	76962.1	9367.25	1.31142e-14	306903
```

Looking at PRS.R2 - this is R-squared - the higher the R-squared, the
better the fit of the data to the model.

## 2. Apply the best PRS in the testing data set

Look for the validation.snp file and filter SNPs that are working best
at the optimal p value threshold of **0.46875** indicated in the
**validation.summary** file

This is to get the SNPs that overcome the threshold p-value defined

``` bash
awk '$4 <= 0.46875' validation.snp | awk '{print $2}' > PRS_snps
```

**Lets make a discovery data set of the SNPs for the best predictive
PRS**\*

``` bash
instead:
  
awk 'NR==FNR ? a[$1] : $3 in a' PRS_snps ASN.gwasqc.txt > temp2
  
( echo -e "CHR\tBP\tSNP\tA1\tA2\tN\tSE\tP\tOR\tINFO\tMAF"; cat temp2) > Bestprs_disc
```

***Finally lets run this best PRS in the test dataset and create a
decile plot***

``` bash
Rscript PRSice.R --prsice PRSice_mac --base Bestprs_disc --target ASN --keep test --pheno ASN.pheno --binary F --cov ASN.cov --no-clump --keep-ambig --fastscore --bar-levels 1 --base-info INFO:0.4 --quantile 10 --quant-break 1,2,3,4,5,6,7,8,9,10 --quant-ref 1 --out test
```

Outputs from (0,1] (1st decile), the difference then in each of the
deciles of genetic variation from the PRS calvulation. Lines around dots
are the SD around the mean.

The PRS stratifies into groups (in deciles).

The bar chart is showing that we are using all the SNPs because in the
data here, we have already selected for the SNPs which have exceeded are
P-value threshold, we aren't doing any clumping, and we are not removing
an ambiguous SNPs any more. (Basically, ignore this plot lol).
