################### calculate and plot age between groups #######################
### author:EALM
### version:3.0
### date: 08-07-2021
##################
### new

library(ggpubr)
library(tidyverse)
library(data.table)
library(gridExtra)

###############################load data########################################
option_list = list(
  
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="full path to the inputfile containing the columns 'id','fam','mother',father','group','age','gender', and all phenptypes to be analysed", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="output directory", metavar="character")
)

# Parse arguments 
parser  <- OptionParser(option_list=option_list)
opt <- parse_args(parser)

########################### setup variables #####################################
input<-opt$input
outdir<-opt$output

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
} else{
  print("Error: No input files!")
}

if(!is.null(outdir)){
  ## create output dir
  dir.create(outdir,recursive = T)
}else{
  print("Error: No output directory!")
}

################################## main ########################################

trai="age"
sumstats<-pheno_complete2%>%group_by(group)%>% summarise(
  size =sum(!is.na(age)),
  mean = mean(age, na.rm = TRUE),
  sd = sd(age, na.rm = TRUE),
  median = median(age, na.rm = TRUE),
  IQR = IQR(age, na.rm = TRUE),
  min = min(age, na.rm = TRUE),
  max=  max(age, na.rm = TRUE)
  
)
xlabels<-paste0(sumstats$group,"\n","(n=",sumstats$size,")")

equat<-paste("age","~ group")
stat.test<-compare_means(as.formula(equat),data=pheno_complete2,method = "t.test",ref.group = "asympt_carr")
maxy<-max(sumstats$max)
miny<-min(sumstats$min)
frac<-(abs(maxy)-miny)/32
stat.test <- stat.test %>%
  mutate(y.position = maxy-2*frac)

ptfactor=0.5/2.141959
p<-ggplot(data=pheno_complete2,aes_string(x="group",y= "age"  ))+
  geom_boxplot(outlier.shape=21,width=0.6,outlier.size = 1.6,lwd=ptfactor,fatten=3.5/2.141959,
               outlier.stroke =ptfactor)+
  geom_violin(trim=TRUE, alpha=0.4,fill="lightgray",width=0.6,size=ptfactor)+
  ylab("age")+
  stat_pvalue_manual(stat.test,tip.length = 0.005,size=2.8, label ="p.format")+
  scale_x_discrete(labels=xlabels)+
  theme_bw()

ggsave(paste0(outdir,"age.box.pdf"),p,dpi=300)
