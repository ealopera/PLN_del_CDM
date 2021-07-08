################### calculate phenotype models  ##################
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
library(optparse)

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

###############################load data########################################
option_list = list(
  
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="full path to the inputfile containing the columns 'id','fam','mother',father','group','age','gender', and all phenptypes to be analysed", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="output directory", metavar="character"),
  make_option(c("-p", "--phenofile"), type="character", default=NULL,
              help="list of phenotypes to be analised (one line for each phenotype with no header), if left empty, defaults to all other columns", metavar="character")
  
)

# Parse arguments 
parser  <- OptionParser(option_list=option_list)
opt <- parse_args(parser)

########################### setup variables #####################################
input<-opt$input
outdir<-opt$output
phenofile<-opt$phenofile

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
} else{
  print("Error: No input files!")
}

if(!is.null(outdir)){
  ## create output dir
  dir.create(outdir,recursive = T)
}else{
  print("Error: No output directory!")
}

if(!is.null(phenofile)){
  #### make list of phenotypes to analize, remove the unwanted
  phenolist<-fread(phenofile,data.table=F,header=F)
  phenolist<-phenolist$V1
}else{
  print("Warning: No list of phenotypes, evaluating all non-identifier columns")
  to_use<-which( names(pheno_complete2) %!in% c("id","fam","mother","father","group", "age","gender") )
  phenolist<-names(pheno_complete2)[to_use]
}

################################# main #########################################

### make familial matrix
Kmatrix<-makekinship (famid=pheno_complete2$fam, id=pheno_complete2$id, 
                      father.id=pheno_complete2$father,mother.id=pheno_complete2$mother)
###factorize group variable
pheno_complete2$group<-factor(pheno_complete2$group,levels=c("asympt_carr","sympt_carr" ))
###create outcome variable
pheno_complete2$outcome<-ifelse(pheno_complete2$group=="asympt_carr",1,
                                 ifelse(pheno_complete2$group=="sympt_carr",0,NA))
pheno_complete2$outcome<-factor(pheno_complete2$outcome,levels=c(0,1))

### models to create independent results
phenoresults<-c()
  for (cname in phenolist){
    v<-which(names(pheno_complete2)==cname)
    ncase=length(which(pheno_complete2$outcome==1 & !is.na(pheno_complete2[,v])))
    ncontrol=length(which(pheno_complete2$outcome==0 & !is.na(pheno_complete2[,v])))
    if (ncase>25){
      ## lmekin model
      model1<-lmekin( pheno_complete2[,v]  ~ age+gender+age^2+outcome +(1|id),
                      data=pheno_complete2, varlist=list(Kmatrix),na.action=na.omit )
      #store results
      m1<-extract_coxme_table(model1)
      phenoresults<-rbind(phenoresults,c(names(pheno_complete2)[v],"asympt_carrier",ncase,ncontrol,m1[4,]))
    }
  }
# transform into a proper matrix
phenoresults<-as.data.frame(phenoresults)
# transform Ncase,beta, se and p in numeric vectors
phenoresults[,1]=as.character(phenoresults[,1])
phenoresults[,2]=as.character(phenoresults[,2])
phenoresults[,3]=as.numeric(as.character(phenoresults[,3]))
phenoresults[,4]=as.numeric(as.character(phenoresults[,4]))
phenoresults[,5]=as.numeric(as.character(phenoresults[,5]))
phenoresults[,6]=as.numeric(as.character(phenoresults[,6]))
phenoresults[,7]=as.numeric(as.character(phenoresults[,7]))
phenoresults[,8]=as.numeric(as.character(phenoresults[,8]))

names(phenoresults)<-c("Trait","Cases","Ncase","Ncontrols","Effect","StdErr","zvalue","pvalue")
phenoresults<-phenoresults%>% group_by(Cases)%>% mutate(FDR=p.adjust(pvalue, method = "fdr"))

write.table(phenoresults, 
            paste0(outdir,"/phenotypes_age_square.txt"),
            quote=F,row.names=F,sep="\t",col.names=T)

####done ####
