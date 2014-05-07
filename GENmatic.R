#Do not change these lines unless you know what you are doing.
args=commandArgs(trailingOnly=TRUE)
require(survival)
require(gap)
scratch=args[1]
folder=args[2]
datafile=args[3]
wd=paste0(scratch,folder)
setwd(wd)
#This Assumes your Fitting.R funciton is in $SCRATCH/GENmatic/Fitting.R
source(paste0(scratch,"GENmatic/Fitting.R"))
data=read.csv(paste0(wd,datafile))

######################################
##### START TO EDIT THE FILE HERE#####
######################################

#Name of directory in scinet with plink (use trailing /)
pd="/home/w/wxu/oespinga/software/plink/plink-1.07-x86_64/"

#Enter your calls to GENfit here. You can use some sort of apply if you want.
#Make sure you set pd=pd and wd=wd
GENfit(data[,c(1,1)],data[,c("SvRfs","Rfs")],data$PC1,data$SEX,"coxph","additive",
       "thinned","gwastest",qq=T,manhattan=T,pd=pd,wd=wd,
       topn=10,topprop=0.1,topcut=0.05)
GENfit(data[,c(1,1)],data$SEX,data$PC1,NULL,"logistic","additive",
       "thinned","logistictest",qq=T,manhattan=T,pd=pd,wd=wd,
       topn=10,topprop=0.1,topcut=0.05)
GENfit(data[,c(1,1)],data$SvRfs,data$PC1,NULL,"linear","additive",
       "thinned","lineartest",qq=T,manhattan=T,pd=pd,wd=wd,
       topn=10,topprop=0.1,topcut=0.05)
