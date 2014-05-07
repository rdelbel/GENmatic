ossystem=function(x){
 if(.Platform$OS.type=="windows") shell(x)
 else system(x)
}

px<-function(x)
  print.xtable(xtable(x),include.rownames=F,table.placement="H",sanitize.text.function=identity)

#'Make Rplink.txt file
#'
#'Makes Rplink.txt file for fash coxph gwas. 
#'Any element of the outcome,strata, or null vectors will correspond to either 
#'the column index (not including the 2 id columns) or column name of the 
#'corresponding covaraite in the covar.txt file. Any element of these vectors will 
#'be treated as an index if it consists of the digits 0-9 and will be considered
#'a column name otherwise. Even if you input an element as "101" (with quotes) it will
#'be treated as an index to the 101th column not the column name 101. Only one stratification
#'variable can be used. If you want to use more than one, create a variable that incorporates all of
#'the strata in one column. 
#'
#'@param outcome A Vector or Matrix/Dataframe corresponding to the outcome. If survival outcome then
#'the first column is the time, and the second column is the status.
#'@param covariates vector/Dataframe corresponding to confounding factors. Default
#'(no covariates) is NULL
#'@param strata vector corresponding to stratification variable. This is only for coxph models. Default
#'(no strata) is NULL. 
#'@param kind additive/dominant/recessive
#'@param wd string corresponding to the directory you want to output the 
#'Rplink.R file to. Default is to the current working directory
#'@param cname the name of the plink covariate file (no extension) default is covar
#'
#'@export
#'@import xtable stringr
makeRplink<-function(outcome,covariates=NULL,strata=NULL,kind="additive",wd="",cname="covar"){  
  makeRplinkinput<-function(x){
    if(is.null(x)) return(NULL)
    x<-as.character(x)
    sapply(x,function(y){
      z<-gsub("[0-9]","",y)
      if(nchar(z)==0)
        return(y)
      return(which(covarnames==y))           
    })}
  if(kind=="additive"){
    snptxt="snp<-as.double(snp[sorted])"
  }else if(kind=="dominant"){
    snptxt="snp<-ifelse(as.double(snp[sorted])==0,0,1)"
  }else if(kind=="recessive"){
    snptxt="snp<-ifelse(as.double(snp[sorted])==2,1,0)"    
  }
  covarnames<-names(read.table(paste0(wd,cname,".txt"),header=T))[-c(1,2)]
  outcome<-makeRplinkinput(outcome)
  strata<-makeRplinkinput(strata)
  covariates<-makeRplinkinput(covariates)
  snpstext=ifelse(!is.null(covariates),'cbind(snp,X)','snp')  
  subsnpstext=ifelse(!is.null(covariates),'cbind(snp[nm],X[nm,])','snp[nm]')    
  subnewstratatext=ifelse(!is.null(strata),'as.integer(c(1 * (diff(as.numeric(strata[nm])) !=0), 1))','newstrata[nm]')
  outapply= "apply(GENO,2,f1)"
  a4<-''
  a1<-paste0('require(survival)
             Rplink<-function(PHENO,GENO,CLUSTER,COVAR){
             time<-COVAR[,',outcome[1],']
             nn=length(time)')
  if(!is.null(strata)){
    a2<-paste0('strata<-COVAR[,',strata,']
               sorted <- order(strata, time)
               strata <- strata[sorted]
               newstrata<-as.integer(c(1 * (diff(as.numeric(strata)) !=0), 1))')
  }else{
    a2<-'sorted<-order(time)
    newstrata<-as.integer(rep(0,nn))'
  }
  a3<-paste0('time<-time[sorted]
             status<-as.integer(COVAR[sorted,',outcome[2],'])')
  if(!is.null(covariates)){
    a4<-paste0('X<-COVAR[sorted,',paste0("c(",paste(covariates,collapse=","),")"),',drop=F]
               init=rep(0,ncol(X)+1)')
  }
  a5<-paste0('alpha=qnorm(.975)
             nn=length(time)
             offset<-rep(0,nn)
             weights<-rep(1,nn)')  
  
  a6<-paste0('f1<-function(snp)
{
             ',snptxt,'
             nm<-which(!is.na(snp))    
             if(length(nm)==nn){
             model<-.Call("coxfit6", 20L, time, status,
             ',snpstext,',
             offset,weights,
             newstrata,
             1L,1e-09, 1.818989e-12,
             init,
             1L)
             }else{   
             model<-.Call("coxfit6", 20L, time[nm], status[nm],
             ',subsnpstext,',
             offset[nm],weights[nm],
             ',subnewstratatext,',
             1L,1e-09, 1.818989e-12,
             init,
             1L)
             }
             coef=model$coef[1]
             var=model$imat[1]
             if(any(model$flag==1000,abs(coef)>10,coef==0)){
             r<-rep(NA,4)
             }else{      
             r<-c(exp(c(coef,coef-alpha*sqrt(var),coef+alpha*sqrt(var))),
             1-pchisq(coef^2/var,1))
             }
             c(4,r)
}
             ',outapply,'
  }')
a<-paste(a1,a2,a3,a4,a5,a6,sep="\n")


  sink(paste0(wd,"Rplink.R"))
  cat(a)
  sink()

}

#'Extract data from plink bfile into a Data Frame
#'
#'Extracts data from a plink bfile into a Data Frame for use in R.
#'Can extract all data, or a subset of specified snps
#'
#'@param ifile String corresponding to the name of the plink bfile (no extension)
#'@param ofile String corresponding to the name of the resulting file which can be read into a dataframe. 
#'It will be called ofile.raw. By default ofile=ifile
#'@param snps A vector of snps that will be in the dataframe. By default all snps will be selected
#'@param pd directory of plink. Use only if plink is not in your path
#'@param wd directory that all files are located in. Is by default set to the current R working directory
#'
#'@export

plinktoR<-function(ifile,ofile=NULL,snps=NULL,pd="",wd=""){
  if(is.null(ofile)) ofile=ifile
  if(is.null(snps)){
    ossystem(paste0(pd,"plink --noweb  --recodeA --bfile ",wd,ifile," --out ",wd,ofile))
  }else{
    write.table(snps,"plinktoRsnps.txt",quote=F,row.names=F,col.names=F)
    ossystem(paste0(pd,"plink --noweb --extract plinktoRsnps.txt --recodeA --bfile ",wd,ifile," --out ",wd,ofile))
  }
  table=read.table(paste0(wd,ofile,".raw"),head=T)[,-c(3,4,5,6)] 
  colnames(table)[-c(1,2)]=sapply(colnames(table)[-c(1,2)],function(name)substr(name,1,nchar(name)-2))
  return(table)
}

#' Fit coxph genetic analysis on windows. 
#' 
#' This function uses R so probably not a good idea to use it on a GWAS. It calls coxph directly from C
#' so it will be an order of magnitute faster than calling coxph from R. 
#' 
#'@param outcome A Dataframe/matrix. the first column is the time, and the second column is the status.
#'@param covariates vector/Dataframe corresponding to confounding factors. Default
#'(no covariates) is NULL
#'@param strata vector corresponding to stratification variable.Default is no strata.
#'@param ifile String corresponding to the name of the plink bfile (no extension)
#'@param ofile String corresponding to the name of the resulting files
#'@param pd directory of plink. Use only if plink is not in your path
#'@param wd directory that all files are located in. Is by default set to the current R working directory
#'@export

GENwincox<-function(outcome,covariates,strata,kind,ifile,ofile,pd="",wd=""){
  geno=plinktoR(ifile,ofile,NULL,pd,wd)[,-c(1,2)] 
  
  #colnames(geno)=sapply(colnames(geno),function(x) substr(x,1,nchar(x)-2))
  if(kind=="dominant"){
    geno=sapply(geno,function(column) ifelse(column!=0,1,0))
  }else if(kind=="recessive"){
    geno=sapply(geno,function(column) ifelse(column==2,1,0))        
  }
  
  if(is.null(strata)){
    sorted=order(outcome[,1])
  }else{
    sorted=order(strata,outcome[,1])
    strata <- strata[sorted]
    newstrata<-as.integer(c(1 * (diff(as.numeric(strata)) !=0), 1))        
  }      
  
  time<-outcome[sorted,1]
  status<-as.integer(outcome[sorted,2])
  if(!is.null(covariates)){
    covariates<-as.matrix(covariates)[sorted,,drop=F]
  }       
  init=rep(0,ncol(covariates)+1)
  alpha=qnorm(.975)
  nn=length(time)
  offset<-rep(0,nn)
  weights<-rep(1,nn)
  fitted=data.frame(t(sapply(geno,function(snp)
  {
    snp<-as.double(snp[sorted])
    nm<-which(!is.na(snp))    
    if(length(nm)==nn){
      model<-.Call("coxfit6", 20L, time, status,
                   cbind(snp,covariates),
                   offset,weights,
                   newstrata,
                   1L,1e-09, 1.818989e-12,
                   init,
                   1L)
    }else{   
      model<-.Call("coxfit6", 20L, time[nm], status[nm],
                   cbind(snp[nm],covariates[nm,]),
                   offset[nm],weights[nm],
                   as.integer(c(1 * (diff(as.numeric(strata[nm])) !=0), 1)),
                   1L,1e-09, 1.818989e-12,
                   init,
                   1L)
    }
    coef=model$coef[1]
    var=model$imat[1]
    if(any(model$flag==1000,abs(coef)>10,coef==0)){
      r<-rep(NA,4)
    }else{      
      r<-c(exp(c(coef,coef-alpha*sqrt(var),coef+alpha*sqrt(var))),
           1-pchisq(coef^2/var,1))
    }
    r
  })))
  fitted=data.frame(rownames(fitted),fitted)
  colnames(fitted)=c("SNP","HR","L95","U95","P")
  bpchr=read.table(paste0(wd,ifile,".bim"),stringsAsFactors=F)[,c(1,2,4)]
  colnames(bpchr)=c("CHR","SNP","BP")
  return(merge(bpchr,fitted,by="SNP"))
}

#' Use plink to fit models
#' 
#' Use plink to fit models. If covariates are to be used assumes that pheno and covar
#'  files exist and are in the form ofile_pheno.txt and ofile_covar.txt. If coxph is
#'  to be used assumes ofile_Rplink.R exists and the system is unix
#'  
#'@param ifile String corresponding to the name of bed/bim/fam plink files. 
#'@param ofile String corresponding to the name of the resulting files. Default is ifile
#'@param type type of analysis (coxph, linear, logistic)
#'@param pd Directory that plink is located in. Only specify if it is not in your path
#'@param wd Directory that all input/output files live in. By default the current R working directory is used.
#'@export

plinkfit<-function(ifile,ofile=NULL,type,covar,pd,wd){
  if(is.null(ofile)) ofile=ifile
  if(type=="coxph"){
    if(.Platform$OS.type=="windows") stop("can not use plink to fix coxph on windows")
    require(Rserve)
    Rserve(args="--no-save")    
    ossystem(paste0(pd,"plink --noweb --bfile ",wd,ifile,"  --covar ",wd,ofile,"_covar.txt  --pheno ",wd,ofile,"_pheno.txt --R ",wd,"Rplink.R --out ",wd,ofile))
    ossystem("kill -9 `ps -eo pid,args | grep Rserve | grep -v grep | awk '{print $1}'`")
    filename=paste0(wd,ofile,".auto.R")    
    fitted=read.table(filename,stringsAsFactors=F)
    colnames(fitted)=c("CHR","SNP","BP","A1","HR","L95","U95","P")
    write.table(fitted,filename,quote=F,row.names=F)
  }
  else if(type=="logistic"&!covar){
    ossystem(paste0(pd,"plink --noweb --bfile ",wd,ifile,"  --pheno ",wd,ofile,"_pheno.txt --assoc --ci 0.95  --out ",wd,ofile))
  }
  else if(type=="logistic"&covar){
    ossystem(paste0(pd,"plink --noweb --bfile ",wd,ifile,"  --covar ",wd,ofile,"_covar.txt  --pheno ",wd,ofile,"_pheno.txt --logistic --ci 0.95 --hide-covar --out ",wd,ofile))
  }else if(type=="linear"&!covar){
    ossystem(paste0(pd,"plink --noweb --bfile ",wd,ifile,"  --pheno ",wd,ofile,"_pheno.txt --assoc --ci 0.95  --out ",wd,ofile))
    
  }else if(type=="linear"&covar){
    ossystem(paste0(pd,"plink --noweb --bfile ",wd,ifile,"  --covar ",wd,ofile,"_covar.txt   --pheno ",wd,ofile,"_pheno.txt --linear --ci 0.95 --hide-covar --out ",wd,ofile))
    
  }
}

#'Fit genetic models
#'
#'Fit genetic models. Currently can only fit coxph models
#'
#'@param id dataframe with 2 columns corresponding to the id in the .fam plink file
#'@param outcome dataframe/vector of outcomes. If coxph model first column is time, second column is status
#'@param covariates dataframe of covaraites. Default is NULL for no covariates
#'@param strata vector of statification membership status. Works only with coxph.
#'@param type Type of analysis to do. Must be one of "linear", "logistic", "coxph".
#'@param kind string corresponding to additive/dominant/recessive model. Currently only
#'coxph can do a non-additive analysis
#'@param ifile String corresponding to the name of bed/bim/fam plink files. Default is GENmatic
#'@param ofile String corresponding to the name of the resulting files. Default is GENmatic
#'@param qq Boolean indicating if a QQ plot of the p-values will be created. Default is F
#'@param manhattan Boolean indicating if a manhattan plot of the p-values will be created. Default is F
#'@param topn If not equal to 0, will create a file with the top n number of snps. Default is 10.
#'@param topprop If not equal to 0, will create a file with the top n proportion of snps (0<=n<=1)
#'@param topcut If not equal to 0, will create a file with all snps with a p-value <= n
#'@param pd Directory that plink is located in. Only specify if it is not in your path
#'@param wd Directory that all input/output files live in. By default the current R working directory is used.
#'@export
#'@import survival
GENfit<-function(id,outcome,covariates=NULL,strata=NULL,type,kind="additive",
                 ifile="GENmatic",ofile="GENmatic",qq=F,manhattan=F,topn=10,topprop=0,topcut=0,pd="",wd=""){
  if(!type%in%c("coxph","logistic","linear"))
    stop("type must be coxph, logistic, or linear")
  if(type!="coxph"&!is.null(strata))
    stop("strata is only defined when using coxph")
  if(type!="coxph"&kind!="additive")
    stop("kind can only not be additive when using coxph")
  if(!kind%in%c("additive","dominant","recessive"))
    stop("kind must be one of additive, dominant, recessive")
  cname=paste0(ofile,"_covar")
  pname=paste0(ofile,"_pheno")
  if(type=="coxph"){
    if(.Platform$OS.type=="unix"){      
    numbers=make_covar_pheno(type,id,outcome,covariates,strata,cname,pname,wd,call=T)   
    makeRplink(1:2,numbers[[1]],numbers[[2]],kind,wd,cname)    
    plinkfit(ifile,ofile,type,covar=T,pd,wd)    
    filename=paste0(wd,ofile,".auto.R")
    }else{        
      fitted=GENwincox(outcome,covariates,strata,kind,ifile,ofile,pd,wd)
      filename=paste0(wd,ofile,".windows.R")
      write.table(fitted,filename,row.names=F,quote=F)
    }
  }else if(type=="logistic"){
    make_covar_pheno(type,id,outcome,covariates,strata,cname,pname,wd)
    if(is.null(covariates)){
      plinkfit(ifile,ofile,type,covar=F,pd,wd)
      filename=paste0(wd,ofile,".assoc")
    }    
    else{
      plinkfit(ifile,ofile,type,covar=T,pd,wd)
      filename=paste0(wd,ofile,".assoc.logistic")
    }
  }else if(type=="linear"){
    make_covar_pheno(type,id,outcome,covariates,strata,cname,pname,wd) 
    if(is.null(covariates)){
      plinkfit(ifile,ofile,type,covar=F,pd,wd)
      
      filename=paste0(wd,ofile,".qassoc")
    }
    else{
      plinkfit(ifile,ofile,type,covar=T,pd,wd) 
      filename=paste0(wd,ofile,".assoc.linear")
    }
  }   
  fitted=read.table(filename,stringsAsFactors=F,head=T)  
  si=snpinfo(ifile,pd,wd)
  a=merge(si,fitted[,!colnames(fitted)%in%c("A1","A2","TEST")],by="SNP",check.names = FALSE)
  t1=which(colnames(a)=="CHR")
  t2=which(colnames(a)=="BP")
  a=a[,c(t1,1,t2,setdiff(1:ncol(a),c(1,t1,t2)))]  
  finalfile=paste0(wd,ofile,"_results.txt")
  write.table(a,finalfile,row.names=F,quote=F)   
    if(manhattan) manhattan(finalfile)
    if(qq) qq(finalfile)    
    topsnps(finalfile,ofile,topn=topn,topprop=topprop,topcut=topcut,wd=wd)  
 unlink(c(paste0(wd,ifile,".frq"),paste0(wd,ifile,".nosex"),paste0(wd,ifile,".hwe"),paste0(wd,ifile,".log")))
}


#'Make covariate file and phenotype file
#'
#'Make covariate file and phenotype file. These files can be made seperate if necessary, but will always be made via GENfit call
#'
#'@param type What kind of model you want to make the files for. Currently can use coxph, linear, logistic
#'@param id datframe corresponding to the 2 column id used in the .fam file
#'@param outcomedataframe of covariates wanted in the model. Do not include intercept term. Make sure all factors are actually factors
#'@param strata vector of stratification status. Only applicable for coxph
#'@param cname name of outputed covar file. Default is covar
#'@param pname name of outputed phenotype file. Default is pheno
#'@param wd name of directory (including trailing slash) where the outputed files will go. By default files go to the current working directory
#'@export
make_covar_pheno<-function(type,id,outcome,covar=NULL,strata=NULL,cname="covar",pname="pheno",wd="",call=F){
    if(!type%in%c("coxph","logistic","linear"))
    stop("Type must be coxph, logistic, linear")
  if(type!="coxph"&&!is.null(strata))
    stop("strata can only be defined for coxph")  
  C=is.null(covar)
  S=is.null(strata)
  nstrat=NULL
  numcov=NULL
  if(length(which(outcome==-9)!=0))
    stop("Data can not have any -9 values when using plink")
  if(!C){
    if(length(which(covar==-9)!=0))
      stop("Data can not have any -9 values when using plink")
    covar=data.frame(covar)
  }
    if(type=="coxph"){
    pheno=cbind(id,1)
    if(!C&&!S){
      if(length(which(strata==-9)!=0))
        stop("Data can not have any -9 values when using plink")
      cov=model.matrix(~.,model.frame(~.,covar,na.action=NULL))[,-1,drop=F]
      numcov=3:(ncol(cov)+2)
      nstrat=ncol(cov)+3
      data=cbind(id,outcome,cov,strata)
    }else if(C){      
      nstrat=3
      data=cbind(id,outcome,strata)
    }else if(S){
      cov=model.matrix(~.,model.frame(~.,covar,na.action=NULL))[,-1,drop=F]
      numcov=3:(ncol(cov)+2)
      data=cbind(id,outcome,cov)
    }else{
      data=cbind(id,outcome)
    }
    data[is.na(data)]=-9
    write.table(data,paste0(wd,cname,".txt"),quote=F,row.names=F)
    write.table(pheno,paste0(wd,pname,".txt"),quote=F,row.names=F,col.names=F)
    if(call){
      return(list(numcov,nstrat))
    }
    }
  else if(type=="logistic"){
    if(length(unique(outcome))!=2) stop("logistic must have exactly 2 outcomes")
    ma=which(outcome==max(outcome))
    mi=setdiff(setdiff(1:length(outcome),ma),which(outcome==-9))    
    outcome[ma]=2
    outcome[mi]=1
    pheno=cbind(id,outcome)
    if(!C){
      cov=model.matrix(~.,model.frame(~.,covar,na.action=NULL))[,-1,drop=F]
      write.table(data.frame(id,cov),paste0(wd,cname,".txt"),quote=F,row.names=F)    
    }
    write.table(pheno,paste0(wd,pname,".txt"),quote=F,row.names=F,col.names=F)
  }else{
    if(length(which(outcome==-9))!=0)
      stop("Can not use outcome with -9 in plink")
    if(!C){
      cov=model.matrix(~.,model.frame(~.,covar,na.action=NULL))[,-1,drop=F]
      write.table(data.frame(id,cov),paste0(wd,cname,".txt"),quote=F,row.names=F)    
    }
    pheno=cbind(id,outcome)
    write.table(pheno,paste0(wd,pname,".txt"),quote=F,row.names=F,col.names=F)
  }
} 

#' Create a qq plot of genetic analysis p-values
#' 
#' Create a qq plot of genetic analysis p-values
#' @param ifile the name of the output file of the genetic analysis
#' @param title title of plot
#' @param ... arguments to go into plot
#' @export
qq<-function(ifile,title="qq plot of p-values",...){
  results=read.table(ifile,head=T)
  obs=sort(results[,"P"],decreasing=F)
ept=c(1:length(obs))/(length(obs))
plot(-log10(ept),-log10(obs),col=4,xlab="Expected -log10(pvalue)",ylab="Observed -log10(pvalue)",main=title,...)
abline(0,1,col="red")
}
#' Create a qq plot of genetic analysis p-values
#' 
#' Create a qq plot of genetic analysis p-values
#' @param ifile the name of the output file of the genetic analysis
#' @param title title of plot
#' @param ... arguments to go into gap::mhtplot
#' @export
#' @import gap
manhattan<-function(ifile,title="Manhattan Plot",...){
gwas<-read.table(ifile,header=T)
d<-gwas[complete.cases(gwas),c("CHR", "BP","P")]
ord <- order(d$CHR,d$BP)
d <- d[ord,]
top=ceiling(-log10(min(d[,"P"])))
colors <- c(rep(c("blue","red"),15),"red")
oldpar<-par()
par(cex=0.6)
mhtplot(d,control=mht.control(colors=colors,gap=1000),pch=19,srt=0,ylim=c(0,top),...)
axis(2,cex.axis=2)
title(title)
#sline<--log10(3.60036E-05)
#gline<--log10(1.8E-06)
#abline(h=sline,col="blue")
#abline(h=gline,col="green")
abline(h=0)
}

#' Returns a dataframe with descriptive information about each snp
#'
#' Returns a dataframe with descriptive information about each snp
#'@param ifile String corresponding to the name of the plink bfile you want to summarize (no extension)
#'@param pd directory of plink. Use only if plink is not in your path
#'@param wd directory that all files are located in. Is by default set to the current R working directory
#'@export
snpinfo<-function(ifile,pd,wd){
  ossystem(paste0(pd,"plink --noweb --bfile ",wd,ifile," --freq --out ",wd,ifile))
  ossystem(paste0(pd,"plink --noweb --bfile ",wd,ifile," --hardy --out ",wd,ifile))
  f=read.table(paste0(wd,ifile,".frq"),stringsAsFactors=F,head=T)
  h=read.table(paste0(wd,ifile,".hwe"),stringsAsFactors=F,head=T)
  h=h[h$TEST=="ALL",!colnames(h)%in%c("A1","A2","CHR")]
  fh=merge(f,h,"SNP")
  colnames(fh)[which(colnames(fh)=="P")]="HWE.P"
  fh=fh[,!colnames(fh)%in%c("TEST","CHR")]
}


pts<-function(data){
  xt=xtable(data)
  digits(xt)=c(rep(2,ncol(data)),8)
  print.xtable(xt,include.rownames=F,table.placement="H")
}
#' Summarize top snp information
#' 
#' Summarize top snp information. Can get top n snps, a proportion of top snps, or all snps with p-values
#' less than a specified cutoff. Each option will write to a unique file
#'@param ifile String corresponding to the name of the results file you want to get the top snps from
#'@param ofile String corresponding to the name of the resulting files (by default same as ifile)
#'@param topn Specify how many snps you want to retreive
#'@param topprop specifiy what proportion of snps you want to retreive
#'@param topcut specify the p-value cutoff you want to retreive
#'@param wd directory that all files are located in. Is by default set to the current R working directory

#'@export
topsnps<-function(ifile,ofile=NULL,topn=0,topprop=0,topcut=0,wd=""){
  if(topprop>1 | topprop<0) stop("topprop must be between 0 and 1")
  if(is.null(ofile)) ofile=ifile
  data=read.table(ifile,header=T,stringsAsFactors=F)
  ord=order(data[,"P"],decreasing=F)
  temp=which(colnames(data)=="P")
  data=data[ord,c(setdiff(1:ncol(data),temp),temp)]
  if(topn!=0){
    write.csv(data[1:topn,],paste0(wd,ofile,"_ntopsnps.csv"),row.names=F,quote=F)    
  }
  if(topprop!=0){
    write.csv(data[1:ceiling(nrow(data)*topprop),],paste0(wd,ofile,"_proptopsnps.csv"),row.names=F,quote=F)
  }
  if(topcut!=0){
    write.csv(data[which(data[,"P"]<=topcut),],paste0(wd,ofile,"_cutofftopsnps.csv"),row.names=F,quote=F)
  }
}
