################### calculate permutation for phenotypes ##################
### author:EALM
### version:4.0
### date: 07-07-2021
##################
### new

####load libraries ###
library(tidyverse)
library(data.table)
library(kinship2)
library(coxme)
###########functions ########

##extract coxme table as a result table
extract_coxme_table <- function (mod){
  beta <- mod$coefficients$fixed
  nvar <- length(beta)
  nfrail <- nrow(mod$var) - nvar
  se <- sqrt(diag(mod$var)[nfrail + 1:nvar])
  z<- round(beta/se, 2)
  p<- signif(1 - pchisq((beta/se)^2, 1), 2)
  table=data.frame(cbind(beta,se,z,p))
  return(table)
}

##function "not in"
`%!in%`<-Negate(`%in%`)

## function to calculate q-value for each p-value (x)
findq<-function(x,data){
  n_less<-length(data[ data$pvalue<=x,"pvalue" ])
  data$preq<-n_less/nrow(data)
}

### function to re-organize q-values from lesser to higher  
FDR_calc<-function(phenodb, permdb,cases){
  pheno1<-phenodb[phenodb$Cases==cases,]
  data1<-permdb[permdb$Cases==cases,]
  #data1<-permdb[permdb$Cases=="asymp_carrier",]
  q_pre<-unlist(lapply(pheno1$pvalue,findq,data1))
  pheno1[,10]<-q_pre
  pheno1<-pheno1[order(pheno1$pvalue),]
  head(pheno1)
  pheno1[,11]<-(c(pheno1[1,10],rep(NA,nrow(pheno1)-1)))
  for (n in 2:nrow(pheno1)){
    if (pheno1[n,10]< pheno1[n-1,10]){
      pheno1[n,11]<-pheno1[n-1,10]} 
    else
    {pheno1[n,11]<-pheno1[n,10]}
    #pheno1[,10]!=pheno1[,11]
  }
  return(pheno1)
}
###############################load data########################################
option_list = list(
  
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="full path to the inputfile containing the columns 'id','fam','mother',father','group','age','gender', and all phenptypes to be analysed", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="output directory", metavar="character"),
  make_option(c("-i", "--independent"), type="character", default=NULL,
              help="complete path to the independent models results", metavar="character"),
  make_option(c("-n", "--nperm"), type="character", default=NULL,
              help="number of permutations", metavar="character")
  
)

# Parse arguments 
parser  <- OptionParser(option_list=option_list)
opt <- parse_args(parser)

########################### setup variables #####################################
input<-opt$input
outdir<-opt$output
nperm<-opt$nperm
indep<-opt$independent

### check variables
if(!is.null(input)){
  ##load data
  pheno_complete2<-fread(input,data.table=F)
  if("id" %!in% colnames(pheno_complete2))
  {print("Error: No id column in input!");} 
  if("fam" %!in% colnames(pheno_complete2))
  {print("Error: No fam column in input!");} 
  if("mother" %!in% colnames(pheno_complete2))
  {print("Error: No mother column in input!");} 
  if("father" %!in% colnames(pheno_complete2))
  {print("Error: No father column in input!");} 
  if("group" %!in% colnames(pheno_complete2))
  {print("Error: No group column in input!");} 
  #### list of phenotypes to permute
  to_use<-which( names(pheno_complete2) %!in% c("id","fam","mother","father","group", "age","gender") )
  phenolist<-names(pheno_complete2)[to_use]
} else{
  print("Error: No input files!")
}

if(!is.null(outdir)){
  ## create output dir
  dir.create(outdir,recursive = T)
}else{
  print("Error: No output directory!")
}

if(!is.null(nperm)){
  #### assign number of permutations
  nperm<-as.numeric(nperm)
}else{
  print("Warning: No number of permutations specified, evaluating 100 by default")
  nperm<-100
}

if(!is.null(indep)){
  #### make list of phenotypes to analize, remove the unwanted
  phenoresults<-fread(indep,data.table=F)
}else{
  print("Error: No independent file results !!")
  
}

################################# main #########################################


################################# main #########################################

### make familial matrix
Kmatrix<-makekinship (famid=pheno_complete2$fam, id=pheno_complete2$id, 
                      father.id=pheno_complete2$father,mother.id=pheno_complete2$mother)
###factorize group variable
pheno_complete2$group<-factor(pheno_complete2$group,levels=c("asympt_ncarr",
                                                             "sympt_ncarr","asympt_carr","sympt_carr" ))
###create outcome variable
pheno_complete2$outcome<-ifelse(pheno_complete2$group=="asympt_carr",1,
                                ifelse(pheno_complete2$group=="asympt_ncarr",0,NA))
pheno_complete2$outcome<-factor(pheno_complete2$outcome,levels=c(0,1))
### models to create permutation results

phenoresults_permuted<-c()
for (cname in phenolist){
  v<-which(names(pheno_complete2)==cname)
  ncase=length(which(pheno_complete2$outcome==1 & !is.na(pheno_complete2[,v])))
  ncontrol=length(which(pheno_complete2$outcome==0 & !is.na(pheno_complete2[,v])))
  if (ncase>25){
  for (n in 1:nperm) {
  #### permutations
  ## I choose id as the genotype id, because Kmatrix is using this one
  pheno_complete2$id2<- pheno_complete2$id
  pheno_only_cols<-which(names(pheno_complete2) %in% c("id","outcome","group") ) 
  pheno_only<-pheno_complete2[,-pheno_only_cols] ##make phenotype only dataset
  geno_only<-pheno_complete2[,pheno_only_cols] ###make genotype only set
  pheno_id_permuted<-data.frame(cbind(sample(pheno_complete2$id2),pheno_complete2$id)) ### make shuffling
  ##format shuffled data
  names(pheno_id_permuted)<-c("id2","id")
  pheno_id_permuted$id<-as.factor(pheno_id_permuted$id)
  geno_only$id<-as.factor(geno_only$id)
  pheno_perm<-left_join(pheno_id_permuted,pheno_only,by="id2") ##shuffle phenotypes together
  pheno_perm<-left_join(pheno_perm,geno_only,by="id") ### join original genotypes with shuffled phenotypes
  v<-which(names(pheno_perm)==cname)
  modelperm1<-lmekin( pheno_perm[,v]  ~ age+gender+age^2 +outcome +(1|id),
                      data=pheno_perm, varlist=list(Kmatrix),na.action=na.omit )
  #store results
  mperm1<-extract_coxme_table(modelperm1)
  phenoresults_permuted<-rbind(phenoresults_permuted,c(names(pheno_complete2)[v],"asymp_carrier",ncase,ncontrol,mperm1[4,]))
  }
  }
}

# transform into a proper matrix
phenoresults_permuted<-as.data.frame(phenoresults_permuted)
# transform Ncase,beta, se and p in numeric vectors
phenoresults_permuted[,1]=as.character(phenoresults_permuted[,1])
phenoresults_permuted[,2]=as.character(phenoresults_permuted[,2])
phenoresults_permuted[,3]=as.numeric(as.character(phenoresults_permuted[,3]))
phenoresults_permuted[,4]=as.numeric(as.character(phenoresults_permuted[,4]))
phenoresults_permuted[,5]=as.numeric(as.character(phenoresults_permuted[,5]))
phenoresults_permuted[,6]=as.numeric(as.character(phenoresults_permuted[,6]))
phenoresults_permuted[,7]=as.numeric(as.character(phenoresults_permuted[,7]))
phenoresults_permuted[,8]=as.numeric(as.character(phenoresults_permuted[,8]))

names(phenoresults_permuted)<-c("Trait","Cases","Ncase","Ncontrols","Effect","StdErr","zvalue","pvalue")

################################################################################
##################calculate q-value for permutations
################################################################################

C1<-FDR_calc(phenodb=phenoresults, permdb=phenoresults_permuted,cases="asymp_carrier")
new_pheno_results<-C1
names(new_pheno_results)[11]<-"permutation_FDR"
new_pheno_results<-new_pheno_results[,-10]
write.table(new_pheno_results,paste0(outdir,"/phenotypes_age_square_with_permFDR.txt"),
            quote=F,sep="\t",row.names = F)

#### done ######