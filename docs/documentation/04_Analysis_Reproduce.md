---
layout: page
title: Analysis Reproduce
---

We provided the source codes to reproduce all analyses both in simulations and real data applications. Note that, the LD matrices were computed from individual-level genotypes from UK Biobank in our study, which requires a separate data access application to UK Biobank directly. While uploading and sharing LD matrices is impractical, we hope the code and the workflow will help users reproduce our results or apply the software to their own data. Additionally, users may compute their own LD matrices using the 1,000 Genomes Project or other appropriate reference datasets.

Simulations
-------------------
We performed comprehensive simulations across 40 different settings to highlight the advantage of sdSuSiE. Compared methods include univariate sex-dimorphic analysis, stepwise regression, SuSiE-modify, MESuSiE-modify, and sdSuSiE. Specifically, we randomly selected [100 genomic regions](https://github.com/yuanzhongshang/sdSuSiE/blob/main/reproduce/simulation/simulation_region.txt) from the 376 regions analyzed in the real data applications. The length of each region ranged from 1Mb to 1.78 Mb (mean = 1.12 Mb; median = 1.03 Mb) and the number of SNPs per region ranged from 611 to 6,345 (mean = 2,772 SNPs; median = 2,758 SNPs). You can find the source codes for each setting:
- [different causal SNP numbers](https://github.com/yuanzhongshang/sdSuSiE/blob/main/reproduce/simulation/simulation_ncausal.sh)
- [different heritability sizes of common and sex-dimorphic effect](https://github.com/yuanzhongshang/sdSuSiE/blob/main/reproduce/simulation/simulation_homoh2.sh)
- [different heritability ratios of common effect to sex-dimorphic effect](https://github.com/yuanzhongshang/sdSuSiE/blob/main/reproduce/simulation/simulation_heteh2.sh)
- [different correlations between common and sex-dimorphic effect](https://github.com/yuanzhongshang/sdSuSiE/blob/main/reproduce/simulation/simulation_rhocd.sh)
- [different male-to-female sample size ratios](https://github.com/yuanzhongshang/sdSuSiE/blob/main/reproduce/simulation/simulation_samplesize.sh)
- [different MAFs of causal SNPs](https://github.com/yuanzhongshang/sdSuSiE/blob/main/reproduce/simulation/simulation_maf.sh)
- [extreme settings](https://github.com/yuanzhongshang/sdSuSiE/blob/main/reproduce/simulation/simulation_extreme.sh)
- different reference panels: [in-sample LD](https://github.com/yuanzhongshang/sdSuSiE/blob/main/reproduce/simulation/simulation_LD.sh)
[out-of-sample LD](https://github.com/yuanzhongshang/sdSuSiE/blob/main/reproduce/simulation/simulation_LD1000G.sh)

Real data applications
-------------------
We performed sex-dimorphic fine-mapping analysis on 626 traits (312 quantitative traits and 314 binary traits) across 50 phenotypic categories from three different data sources in the UK Biobank. 
## Step 1:  Sex-stratified GWAS summary statistics download, quality control and univariate sex-dimorphic analysis
We used Slurm to process 626 traits in parallel. This step includes three sub-steps: data download, quality control, and univariate sex-dimorphic analysis. The analyses were submitted to the Slurm scheduler using the following command:
```bash
sbatch 1_sde.sh
```
Here, sex-stratified GWAS summary statistics for males and females were obtained from the [Neale Lab’s round 2 (imputed-v3) GWAS analysis](https://www.nealelab.is/uk-biobank). 
For male- or female-only GWAS summary statistics, we filtered out SNPs with a Hardy-Weinberg equilibrium p-value < 10<sup>−6</sup>, MAF < 0.001, genotype call rate < 95%, strand ambiguous or multi-allelic, and SNPs in the human leukocyte antigen (HLA) regions (chr6: 25Mb - 36Mb). 
Next, we conducted the univariate sex-dimorphic analysis on the remaining SNPs.

## Step 2: Candidate region determination
We used Slurm to process 626 traits in parallel. For each trait in turn, we examined one SNP at a time and identified SNPs that display sex-dimorphic effects at a genome-wide significance threshold (p-value < 5×10<sup>−8</sup>). We then created a 1-Mb window centered on each significant sdSNP (500 Kb upstream and 500 Kb downstream) and merged the overlapped genomic regions. 
```bash
sbatch 2_region.sh
```

## Step 3:  Sex-dimorphic analysis
We uploaded analysis codes for an exemplary trait, diastolic blood pressure (Data-Field 4079_irnt), , allowing users to follow the same procedure to reproduce results for all traits. For this trait, we prepared the inputs and applied sdSuSiE, stepwise regression, SuSiE-modify, and MESuSiE-modify for sex-dimorphic analysis.
```bash
sbatch 3_analysis.R
```
Since we had already constructed different folders for each trait, we used the following bash script to process the remaining traits in parallel.
```bash
sbatch 3_analysis_run.sh
```
All files can be find at (here)[https://github.com/yuanzhongshang/sdSuSiE/blob/main/reproduce/real_data].