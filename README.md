# Analysis for PLN-(arg14del) in four groups of inidviduals

Created by: Esteban Lopera.\
Created on: 07-07-2021. \
Contact information: ealopera@gmail.com


### 1. Independent analysis for quantitative phenotypes
The script Phenotype_analysis_lmekin will calculate the mixed model \
 Phenotype~ age +gender + age^2 + outcome + family  \
where:
  "outcome" is defined as 0="healthy non-carriers" 1="healthy carriers" \
  and "family" is a random intercept based on the familial information

Requirements: \
-Phenotype(s) file: File containing the columns "id", "fam", "mother", "father", "group", "age", "gender". \
  where the colums "id","fam","mother","father" consist of unique identifiers for each individual, listing its family ID, its mother and father, respectively. Family ID connects individuals with blood relationships, as well spouses if there are offsprings \
  The column "group" classifies the individuals according to the cardiomiopathy symtoms and their state of PLN mutation in four groups, namely: "asympt_ncarr" "sympt_ncarr", "asympt_carr", "sympt_carr". \
-phenotypes list (optional, defaults to all columns of the phenotypes file except the ones mentioned above) to be evaluated in a file with a phenotype per line without header. \
-output directory.

### 2. Permutation analysis for quantitative phenotypes
the script Phenotype_permutations will add the permutation_FDR column to the independent models above. this is made in a separate script because it can take a lot of memory and time otherwise, but if you feel like your system can take it you can merge both scripts.

Requirements:
-Phenotype(s) file: the same phenotype described above. \
-independent models result: the output from step one. \
-number of permutations (optional, defaults to 100). \
-output directory. \

