#'Make Rplink.txt file
#'
#'Makes Rplink.txt file for fash coxph gwas. 
#'#'Any element of the outcome,strata, or null vectors will correspond to either 
#'the column index (not including the 2 id columns) or column name of the 
#'corresponding covaraite in the covar.txt file. Any element of these vectors will 
#'be treated as an index if it consists of the digits 0-9 and will be considered
#'a column name otherwise. Even if you input an element as "101" (with quotes) it will
#'be treated as an index to the 101th column not the column name 101. Only one stratification
#'variable can be used. If you want to use more than one, create a variable that incorporates all of
#'the strata in one column. 
#'
#'@param outcome length 2 vector corresponding to time and status
#'@param covariates vector corresponding to confounding factors. Default
#'(no covariates) is NULL
#'@param strata vector corresponding to stratification variable. Default
#'(no strata) is NULL. 
#'@param mcore Dont use this
#'@param fout string corresponding to the directory you want to output the 
#'Rplink.R file. Default (write to working directory) is NULL.
#'@export
#'@import xtable stringr
makeRplink<-function(outcome,covariates=NULL,strata=NULL,mcore=F,fout=NULL){  
  makeRplinkinput<-function(x){
    if(is.null(x)) return(NULL)
    x<-as.character(x)
    sapply(x,function(y){
      z<-gsub("[0-9]","",y)
      if(nchar(z)==0)
        return(y)
      return(which(covarnames==y))           
  })}
  covarnames<-names(read.table("covar.txt",header=T))[-c(1,2)]
  outcome<-makeRplinkinput(outcome)
  strata<-makeRplinkinput(strata)
  covariates<-makeRplinkinput(covariates)
  snpstext=ifelse(!is.null(covariates),'cbind(snp,X)','snp')  
  subsnpstext=ifelse(!is.null(covariates),'cbind(snp[nm],X[nm,])','snp[nm]')    
  subnewstratatext=ifelse(!is.null(strata),'as.integer(c(1 * (diff(as.numeric(strata[nm])) !=0), 1))','newstrata[nm]')
  outapply=ifelse(mcore,"t(do.call(rbind,mclapply(as.list(data.frame(GENO)),f1)))",
                  "apply(GENO,2,f1)")
  a4<-''
  a1<-paste0('require(survival)
             require(multicore)
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
    newstrata<-rep(0,nn)'
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
             snp<-as.double(snp[sorted])
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
  
  if(!is.null(fout)){
    oldwd<-getwd()
    setwd(fout)
    sink("Rplink.R")
    cat(a)
    sink()
    setwd(oldwd)
  }else{
    sink("Rplink.R")
    cat(a)
    sink()
  }  
}


#'Fit genetic models
#'
#'Fit genetic models. Currently can only fit coxph models
#'
#'@param outcome length 2 vector corresponding to time and status
#'@param covariates vector corresponding to confounding factors. Default
#'(no covariates) is NULL
#'@param strata vector corresponding to stratification variable. Default
#'(no strata) is NULL.
#'@param type Type of analysis to do. Currently can only do coxph.
#'@param mcore Dont use this
#'@param fout string corresponding to the directory you want to output the 
#'Rplink.R file. Default (write to working directory) is NULL.
GENfit<-function(outcome,covariates=NULL,strata=NULL,type="coxph",mcore=F,fout=NULL){
  if(type=="coxph"){
    makeRplink(outcome,covariates,strata,mcore,fout)
    if(mcore)require(multicore)
    require(Rserve)
    Rserve(args="--no-save")
    system("plink --noweb --bfile GENmatic --covar covar.txt --pheno pheno.txt --R Rplink.R --out gwas")
  }}