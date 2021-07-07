# Analysis for PLN-(arg14del) in four groups of inidviduals

Created by: Esteban Lopera.\
Created on: 07-07-2021. \
Contact information: ealopera@gmail.com


### 1. Independent analysis for multiple phenotypes
The script Phenotype_analysis_lmekin will calculate the mixed model \
 Phenotype~ age +gender + age^2 + outcome  \
where outcome is defined as 0="healthy non-carriers" 1="healthy carriers"

Requirements: \
-Phenotype(s) file: File containing the columns "id","fam","mother","father","group", "age","gender". \
-list of selected SNPs for GRM  (GDS or plink format). \
-list of covariates (or model). 

### INSTRUCTIONS
If you are going for a quick analysis we recomend to use the GDS version of SAIGE. it uses the genotype (GT) field of the files and each chromosome might take up to 20 min for 40k samples. You will need all your files to be in GDS format. Otherwise, if you need a more precise, but time consuming (~7h per chormosome in 40k samples) anaylsis use normal SAIGE and prepare vcf of bgen files with their respective index files and make sure they have the dosage (DS field). \
Detailed instructions are found in the respective folders.
