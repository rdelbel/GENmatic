parse_GENmatic.log<-function(geno=F,mind=F,maf=F,hh=F){
  con  <- file("GENmatic.log", open = "r")
  result<-c(0,0,0,0)
  while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
    if(length(out<-str_split(oneLine," markers to be included from")[[1]])>1)
      result[4]<-as.numeric(out[1])
    if(geno){
      if(length(out<-str_split(oneLine," SNPs failed missingness test")[[1]])>1)
        result[3]<-as.numeric(out[1])
    }else if(maf){
      if(length(out<-str_split(oneLine," SNPs failed frequency test")[[1]])>1)
        result[3]<-as.numeric(out[1])
    }else if (hh){
      if(length(out<-str_split(oneLine," heterozygous haploid genotypes; set to missing")[[1]])>1)
        result[3]<-as.numeric(out[1])
    }else{      
      if(length(out<-str_split(oneLine," snps removed with --")[[1]])>1)
        result[3]<-as.numeric(out[1])
    }
    if(length(out<-str_split(oneLine," individuals read from")[[1]])>1)
      result[2]<-as.numeric(out[1])
    if(!mind){
      if(length(out<-str_split(oneLine," individuals removed with --")[[1]])>1)
        result[1]<-as.numeric(out[1])
    }else{
      if(length(out<-str_split(oneLine,"individuals removed for low genotyping")[[1]])>1)
        result[1]<-as.numeric(str_split(out[1],fixed(" "))[[1]][1])
    }
  }  
  result[2]<-result[2]-result[1]
  result[4]<-result[4]-result[3]
  close(con)
  return(result)
}

#'Removes all snps with MAF below cutoff 
#'
#'Removes all snps with MAF below cutoff
#'
#'@param cutoff numeric indicating the MAF cutoff to use
#' @param pd directory where plink is. Do not need to specify if plink is in your path
#' or working directory
#'@export
remove_maf<-function(cutoff=0.05,pd=""){
  GENmaticGWASQCcount<<-GENmaticGWASQCcount+1
  cat(paste0("\\section*{",GENmaticGWASQCcount,
             ": Identification of SNPs with low Minor Allele Frequency}"),
      "All SNPs with MAF$\\leq$",cutoff," were removed. ")
  ossystem(paste0(pd,"plink --noweb --bfile GENmatic --freq --out GENmatic"))
  maf <- read.table("GENmatic.frq",h=T,comment.char="`",stringsAsFactors=F);
  maf_remove = maf[maf$MAF<0.05,2]
  if(length(maf_remove)>0){
    cat(length(maf_remove),"SNPs were removed this way.")
    
    write.table(maf_remove,"maf_remove.txt",quote=FALSE,col.names=F,row.names=F)  
    ossystem(paste0(pd,"plink --noweb --bfile GENmatic --maf ",cutoff," --out GENmatic --make-bed"))
  }else{
    cat("There were no such SNPs.")
  }
  QCsummary<<-rbind(QCsummary,c("Low MAF",parse_GENmatic.log(maf=T)))
  
  hist(maf$MAF,main="MAF distribution",xlim=c(0,0.5),ylim=c(0,250000),axes=F,xlab="Minor Allele Frequency",ylab="Frequency",col="blue")
  axis(1,at=c(0,0.1,0.2,0.3,0.4,0.5),labels=c(0,0.1,0.2,0.3,0.4,0.5));
  axis(2,at=c(0,50000,100000,150000,200000,250000),labels=c("0","50K","100K","150K","200K","250K"),las=2)
  abline(v=cutoff, col="black", lty=2);
  legend(x=cutoff,y=200000,legend=paste0(cutoff*100,"% MAF threshold"))                    
}

#' Initialize GWAS Quality Control
#' 
#' Looks in current directory for ifile.bed, ifile.bim, ifile.fam.
#' Reads these files to get initial data. Saves new files as
#' GENmatic.bed/fam/bin so that we do not modify origional data
#' Sets up 2 secret global variables for the rest of the QC process.
#'
#' @param pd directory where plink is. Do not need to specify if plink is in your path
#' or working directory
#' @keywords GWAS
#' @export
initiate_QC<-function(ifile,pd=""){
  #cat("\\section*{1: Identification of Heterozygous Haploid Genotypes} All Heterozygous Haploid Genotypes were removed. ")      
  ossystem(paste0(pd,"plink --noweb --bfile ",ifile," --set-hh-missing --out GENmatic  --make-bed"))
  result <- t(data.frame(parse_GENmatic.log(), stringsAsFactors = F))
  result <- data.frame("Start", result, stringsAsFactors = F)
  colnames(result) <- c("Step", "Individual Removed", "Individual Remaining", 
                        "SNP Removed", "SNP Remaining")
  rownames(result) <- NULL
  result=result[-1,]
  #result<-rbind(result,c("Heterozygous Haploid",parse_GENmatic.log(hh=T)))
  #QCsummary <<- result
  GENmaticGWASQCcount <<- 0
  
  #cat(ifelse(QCsummary[2,3]=="0","No",QCsummary[2,4]), "snps were removed this way")
}


# remove_missing_pheno<-function(pd=""){
#   
#   filenames<-system("ls",intern=T)
#   filenames<-filenames[grepl(".",filenames,fixed=T)]
#   filenames<-split_filenames(filenames)
#   pheno<-read.csv(paste0(get_filename(filenames,"csv"),".csv"))
#   fam=get_filename(filenames,"fam",T)
#   findID<-matrix(unlist(str_split(system(paste("head -2", fam),intern=T),fixed(" "))),nrow=2,byrow=T)[,c(1,2)]
#   ID<-which(apply(findID,2,function(x)length(unique(x)))>1)
#   notID<-setdiff(c(1,2),ID)
#   notIDvalue<-findID[1,notID]
#   if(length(ID)==2){
#     print(ID)
#     if(length(unique(ID))==1){
#       write.table(cbind(as.character(pheno[,1]),notIDvalue),"keep.txt",quote=F,col.names=F,row.names=F)
#     }else
#       stop("fam file ID makes no sense")
#   }
#   else if(ID==1){
#     write.table(cbind(as.character(pheno[,1]),notIDvalue),"keep.txt",quote=F,col.names=F,row.names=F)
#   }else
#     write.table(cbind(notIDvalue,as.character(pheno[,1])),"keep.txt",quote=F,col.names=F,row.names=F)
#   plinkdataname<-get_filename(filenames,"bed")
#   ossystem(paste0(pd,"plink --noweb --bfile GENmatic --keep keep.txt --make-bed --out GENmatic"))
#   QCsummary<<-rbind(QCsummary,c("Missing Phenotype",parse_GENmatic.log()))
# }

#' Remove people with sex problems
#' 
#' Remove people with sex problems.
#' @param pd directory where plink is. Do not need to specify if plink is in your path
#' or working directory
#' @keywords GWAS
#' @export
remove_sex_problems<-function(pd=""){
  GENmaticGWASQCcount<<-GENmaticGWASQCcount+1
  cat(paste0("\\section*{",GENmaticGWASQCcount,": Identification of individuals with discordant sex information}"))  
  
  ossystem(paste0(pd,"plink --noweb --bfile GENmatic --check-sex --out GENmatic"))
  sexcheck<-read.table("GENmatic.sexcheck",header=T)
  SexProblems<-which(sexcheck$STATUS=="PROBLEM")
  if(length(SexProblems)==0){
    log<-parse_GENmatic.log()
    QCsummary<<-rbind(QCsummary,c("Sex Problems",c(0,log[2],0,log[4])))
  }else if(length(SexProblems)==nrow(sexcheck)){
    warning("No sex information in genotype")
    log<-parse_GENmatic.log()
    QCsummary<<-rbind(QCsummary,c("Sex Problems",c(log[2],0,0,log[4])))
  }else{
    write.table(sexcheck[SexProblems,c(1,2)],"SexProblems.txt",quote=F,row.names=F,col.names=F)
    ossystem(paste0(pd,"plink --noweb --bfile GENmatic --remove SexProblems.txt --make-bed --out GENmatic"))
    QCsummary<<-rbind(QCsummary,c("Sex Problems",parse_GENmatic.log()))
  }
  cat("This option uses X chromosome data to determine sex (i.e. based on heterozygosity rates) and flags individuals for whom the reported sex in the PED file does not match the estimated sex (given genomic data).\\\\")
  if(length(SexProblems)!=0){
    print.xtable(xtable(sexcheck[sexcheck$STATUS=="PROBLEM",]),table.placement="H",include.rownames=F)
  }else{
    cat("No individuals had discordant sex information")
  }
}
#' Remove snps and people with too many missing
#' 
#' First removed 100% missing snps, then missing snps, then
#' issing idividuals.
#' @param snp.cutoff numeric corresponding to cutoff for missing snps.
#'  Default is 0.05
#'  @param individual.cutoff numeric corresponding to cutoff for missing
#'  individuals. Default is 0.05
#' @param pd directory where plink is. Do not need to specify if plink is in your path
#' or working directory
#'  @keywords GWAS
#'  @export 
remove_missing<-function(snp.cutoff=0.05,individual.cutoff=0.05,pd=""){
  GENmaticGWASQCcount<<-GENmaticGWASQCcount+1
  cat(paste0("\\section*{",GENmaticGWASQCcount,": Identification of individuals/markers with a high missing rate}"))  
  
  ossystem(paste0(pd,"plink --noweb --bfile GENmatic --missing --out GENmatic"))
  lmiss<-read.table("GENmatic.lmiss",stringsAsFactors=F,colClasses=c("character","character","numeric","numeric","numeric"),header=T,comment.char = "`")
  fully_missing_snps<-length(lmiss[,5][lmiss[,5]==1])
  missing_snps<-length(lmiss[,5][lmiss[,5]>=snp.cutoff])-fully_missing_snps
  log<-parse_GENmatic.log()
  QCsummary<<-rbind(QCsummary,c("Fully Missing Snps",c(0,log[2],fully_missing_snps,log[4]-fully_missing_snps)))
  QCsummary<<-rbind(QCsummary,c(paste0(snp.cutoff*100,"% Missing Snps"),c(0,log[2],missing_snps,log[4]-fully_missing_snps-missing_snps)))
  ossystem(paste0(pd,"plink --noweb --bfile GENmatic --geno ",snp.cutoff," --out GENmatic --make-bed"))  
  ossystem(paste0(pd,"plink --noweb --bfile GENmatic --mind ",individual.cutoff," --out GENmatic --make-bed"))
  log<-parse_GENmatic.log(mind=T)
  QCsummary<<-rbind(QCsummary,c(paste0(individual.cutoff*100,"% missing individuals"),log))
  cat(paste0(ifelse(fully_missing_snps==0,"No",fully_missing_snps),
             " SNPs in the dataset are 100\\% missing",ifelse(fully_missing_snps==0,".",", which were removed from the data first."), 
             " Next SNPs with more than ",
             snp.cutoff*100,"\\% of individuals missing them were removed. ",ifelse(missing_snps==0,"No",missing_snps)," SNPs' call rate were below this cutoff.",
             " Next Individuals with more than ",individual.cutoff*100," \\% SNPs missing were removed. ",ifelse(log[1]==0,"No",log[1]),
             " individuals were removed under this criterion.\\\\"))
  
  ossystem(paste0(pd,"plink --noweb --bfile GENmatic --missing --out GENmatic"))
  imiss<-read.table("GENmatic.imiss",stringsAsFactors=F,header=T,comment.char = "`")
  
  lmiss$logF_MISS = log10(lmiss$F_MISS);
  ## frequency plot
  
  hist(lmiss$logF_MISS,main="SNP missing rate",xlim=c(-3,0),ylim=c(0,100000),axes=F,xlab="Missing rate per SNP",ylab="Frequency",col="blue")
  axis(1,at=c(-3,-2,-1,0),labels=c(0.001,0.01,0.1,1));
  axis(2,at=c(0,20000,40000,60000,80000,100000),labels=c("0","20K","40K","60K","80K","100K"),las=2)
  abline(v=log10(snp.cutoff), col="black", lty=2);
  legend(x=log10(snp.cutoff),y=80000,legend=paste0(snp.cutoff*100,"% missing rate"))
  
  imiss$logF_MISS = log10(imiss$F_MISS);
  head(imiss)
  ## frequency plot
  
  hist(imiss$logF_MISS,main="Missing rate distribution for individuals",xlim=c(-4,0),ylim=c(0,100),axes=F,xlab="Missing rate per individual",ylab="Frequency",col="blue")
  axis(1,at=c(-4,-3,-2,-1,0),labels=c(0.0001,0.001,0.01,0.1,1));
  axis(2,at=c(0,25,50,75,100))
  abline(v=log10(individual.cutoff), col="black", lty=2);
  legend(x=log10(individual.cutoff),y=75,legend=paste0(individual.cutoff*100,"% completion rate"))
}
#' Remove relatives
#' 
#' Removes the larger ID number for all pairs of individual over 
#' specified IBS
#' @param cutoff numeric corresponding to cutoff ibs.
#' Default is 0.25.
#' @param pd directory where plink is. Do not need to specify if plink is in your path
#' or working directory
#' @param skip boolean indicating if you want to skip getting the GENmatic.genome file. Use TRUE
#' only if you want to redo a previous QC in the same sequence and already have the file.
#' @keywords GWAS
#' @export 
remove_relatives<-function(cutoff=0.25,pd="",skip=F){
  GENmaticGWASQCcount<<-GENmaticGWASQCcount+1
  cat(paste0("\\section*{",GENmaticGWASQCcount,": Identification of duplicated or related individuals}"),
      "To detect duplicated or related individuals, identity by state is calculated for each pair of individuals. If two individuals are identified as relatives, as defined by having IBS$\\geq$",cutoff," the one with the higher sample number was removed.") 
  if(!skip){
  ossystem(paste0(pd,"plink --noweb --bfile GENmatic --indep 50 5 2 --out GENmatic"))
  ossystem(paste0(pd,"plink --noweb --bfile GENmatic --extract GENmatic.prune.in --genome --out GENmatic"))
  }
  ibs = read.table("GENmatic.genome",header=T,stringsAsFactors=F)
  ibs<-ibs[ibs$PI_HAT >= cutoff,c(1,2,3,4,10),drop=F]
  if(nrow(ibs)!=0){
    first_degree_relatives<-ibs[,c(3,4)]
    write.table(first_degree_relatives,"first_degree_relatives.txt",quote=FALSE,col.names=F,row.names=F)
    ossystem(paste0(pd,"plink --noweb --bfile GENmatic  --remove first_degree_relatives.txt --out GENmatic --make-bed"))
    cat("People Removed in this step are in bold and listed in the following table.")
    colnames(ibs)
    colnames(ibs)<-sanitizestr(colnames(ibs))
    ibs<-sapply(ibs,function(cols)sanitizestr(cols))
    ibs[,c(3,4)]<-sapply(ibs[,c(3,4)],lbld)
    print.xtable(xtable(ibs),include.rownames=F,table.placement="H",sanitize.text.function=identity)
    #QCsummary[7,-1]<<-parse_GENmatic.log()
    QCsummary<<-rbind(QCsummary,c("Related Individuals",parse_GENmatic.log()))
    
  }else{
    log<-parse_GENmatic.log()
    #QCsummary[7,-1]<<-c(0,log[2],0,log[4])
    QCsummary<<-rbind(QCsummary,c("Related Individuals",c(0,log[2],0,log[4])))
    cat("No indivudals were identified as relatives") 
  }
  
}
#' Remove indviduals with large heterozygocity
#' 
#' Removes all individuals with heterozogocity more than
#' cutoff SD away form the mean.
#' @param SD numeric corresponding to the numer of SD away
#' from the mean we will choose as cutoff. Default is 6.
#' @keywords GWAS
#' @param pd directory where plink is. Do not need to specify if plink is in your path
#' or working directory
#' @export 
remove_heterozygocity<-function(SD=6,pd=""){
  GENmaticGWASQCcount<<-GENmaticGWASQCcount+1
  cat(paste0("\\section*{",GENmaticGWASQCcount,": Identification of individuals with outlying heterozygosity rate}"))  
  
  cat("The outliers in the heterozygosity rate might be indicative of DNA sample contamination or inbreeding.
      The following plot shows the heterozygosity rate of the individuals in our CRC data, and the horizontal dashed lines represent",SD,"standard deviations from the mean.\\\\")
  ossystem(paste0(pd,"plink --noweb --bfile GENmatic --het --out GENmatic"))
  ossystem(paste0(pd,"plink --noweb --bfile GENmatic --missing --out GENmatic"))
  
  
  imiss <- read.table("GENmatic.imiss",h=T);
  imiss$logF_MISS = log10(imiss$F_MISS);
  het=read.table("GENmatic.het",h=T);
  
  het$meanHet = (het$N.NM. - het$O.HOM.)/het$N.NM.;
  
  colors  <- densCols(imiss$logF_MISS,het$meanHet);
  ## plot for heterozygocity rate
  plot(imiss$logF_MISS,het$meanHet, col=colors, xlim=c(-4,0),ylim=c(0,0.5),pch=20, xlab="Proportion of missing genotypes", 
       ylab="Heterozygosity rate",axes=F,main="Genotype failure rate vs. heterozygosity of all individuals");
  axis(2,at=c(0,0.05,0.10,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5),tick=T);
  axis(1,at=c(-4,-3,-2,-1,0),labels=c(0.0001,0.001,0.01,0.1,1));
  abline(h=mean(het$meanHet)-(SD*sd(het$meanHet)),col="RED",lty=2);
  abline(h=mean(het$meanHet)+(SD*sd(het$meanHet)),col="RED",lty=2);
  abline(v=log10(0.05), col="BLUE", lty=2);
  legend(x=log10(0.05),y=0.15,legend="5% completion rate");
  
  low=mean(het$meanHet)-(SD*sd(het$meanHet));
  high=mean(het$meanHet)+(SD*sd(het$meanHet));
  hetero_remove<-het[het$meanHet<low | het$meanHet>high, c("FID","IID")]
  if(nrow(hetero_remove)>0){
    write.table(hetero_remove,"hetero_remove.txt",quote=FALSE,col.names=F,row.names=F);
    ossystem(paste0(pd,"plink --noweb --bfile GENmatic --remove hetero_remove.txt --out GENmatic --make-bed"))
    #QCsummary[7,-1]<<-parse_GENmatic.log()  
    QCsummary<<-rbind(QCsummary,c(paste0("Heterozygocity Rate (-/+ ",SD,"SD)"),parse_GENmatic.log()))
    cat("\n","The following people were removed for having heterozygosity rate",SD,"standard deviations form the mean:")
    print.xtable(xtable(hetero_remove),table.placement="H",include.rownames=F)
  }else{
    log<-parse_GENmatic.log()
    QCsummary<<-rbind(QCsummary,c(paste0("Heterozygocity Rate (-/+ ",SD,"SD)"),0,log[2],0,log[4]))
    cat("No patients had heterozygosity rate",SD,"standard deviations form the mean")
  }
  
}

#' End GWAS QC
#' 
#' Outputs summary of QC process. There are no paramaters.
#' @export
end_QC<-function(){
  GENmaticGWASQCcount<<-GENmaticGWASQCcount+1
  cat(paste0("\\section*{",GENmaticGWASQCcount,": Summary of QC}"))  
  cat("After all QC we have",QCsummary[nrow(QCsummary),3], "individuals, and",QCsummary[nrow(QCsummary),5], "SNPs remaining.\\\\")
  #OQC<<-QCsummary
  temp3<-rep(F,nrow(QCsummary))
  temp5<-rep(F,nrow(QCsummary))  
  for(i in 2:nrow(QCsummary)){
    if(QCsummary[i,3]==QCsummary[i-1,3])
      temp3[i]=T
    if(QCsummary[i,5]==QCsummary[i-1,5])
      temp5[i]=T          
    QCsummary[temp3,3]<<-""
    QCsummary[temp5,5]<<-""  } 
  QCsummary[,c(2,4)][QCsummary[,c(2,4)]=="0"]<-""
  print.xtable(xtable(QCsummary),table.placement="H",include.rownames=F)
}