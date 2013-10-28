parse_reportRx.log<-function(geno=F,mind=F,maf=F,hh=F){
  con  <- file("reportRx.log", open = "r")
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
#'@export
remove_maf<-function(cutoff=0.05){
  reportRxGWASQCcount<<-reportRxGWASQCcount+1
  cat(paste0("\\section*{",reportRxGWASQCcount,
             ": Identification of SNPs with low Minor Allele Frequency}"),
      "All SNPs with MAF$\\leq$",cutoff," were removed. ")
  system("plink --noweb --bfile reportRx --freq --out reportRx")
  maf <- read.table("reportRx.frq",h=T,comment.char="`",stringsAsFactors=F);
  maf_remove = maf[maf$MAF<0.05,2]
  if(length(maf_remove)>0){
    cat(length(maf_remove),"SNPs were removed this way.")
    
    write.table(maf_remove,"maf_remove.txt",quote=FALSE,col.names=F,row.names=F)  
    system(paste("plink --noweb --bfile reportRx --maf",cutoff,"--out reportRx --make-bed"))
  }else{
    cat("There were no such SNPs.")
  }
  QCsummary<<-rbind(QCsummary,c("Low MAF",parse_reportRx.log(maf=T)))
  
  hist(maf$MAF,main="MAF distribution",xlim=c(0,0.5),ylim=c(0,250000),axes=F,xlab="Minor Allele Frequency",ylab="Frequency",col="blue")
  axis(1,at=c(0,0.1,0.2,0.3,0.4,0.5),labels=c(0,0.1,0.2,0.3,0.4,0.5));
  axis(2,at=c(0,50000,100000,150000,200000,250000),labels=c("0","50K","100K","150K","200K","250K"),las=2)
  abline(v=cutoff, col="black", lty=2);
  legend(x=cutoff,y=200000,legend=paste0(cutoff*100,"% MAF threshold"))                    
}

#' Initialize GWAS Quality Control
#' 
#' Looks in current directory for dirty.bed dirty.fam dirty.bin.
#' Reads these files to get initial data. Saves new files as
#' reportRx.bed/fam/bin so that we do not modify origional data
#' Sets up 2 secret global variables for the rest of the QC process.
#' This function takes no paramaters
#' @keywords GWAS
#' @export
initiate_QC<-function(){
  cat("\\section*{1: Identification of Heterozygous Haploid Genotypes} All Heterozygous Haploid Genotypes were removed.")      
  system("plink --noweb --bfile dirty --set-hh-missing --out reportRx  --make-bed")
  result <- t(data.frame(parse_reportRx.log(), stringsAsFactors = F))
  result <- data.frame("Start", result, stringsAsFactors = F)
  colnames(result) <- c("Step", "Individual Removed", "Individual Remaining", 
                        "SNP Removed", "SNP Remaining")
  rownames(result) <- NULL
  result<-rbind(result,c("Heterozygous Haploid",parse_reportRx.log(hh=T)))
  QCsummary <<- result
  reportRxGWASQCcount <<- 1
  
  cat(ifelse(QCsummary[2,3]=="0","No",QCsummary[2,3]), "snps were removed this way")
}

remove_missing_pheno<-function(){
  
  filenames<-system("ls",intern=T)
  filenames<-filenames[grepl(".",filenames,fixed=T)]
  filenames<-split_filenames(filenames)
  pheno<-read.csv(paste0(get_filename(filenames,"csv"),".csv"))
  fam=get_filename(filenames,"fam",T)
  findID<-matrix(unlist(str_split(system(paste("head -2", fam),intern=T),fixed(" "))),nrow=2,byrow=T)[,c(1,2)]
  ID<-which(apply(findID,2,function(x)length(unique(x)))>1)
  notID<-setdiff(c(1,2),ID)
  notIDvalue<-findID[1,notID]
  if(length(ID)==2){
    print(ID)
    if(length(unique(ID))==1){
      write.table(cbind(as.character(pheno[,1]),notIDvalue),"keep.txt",quote=F,col.names=F,row.names=F)
    }else
      stop("fam file ID makes no sense")
  }
  else if(ID==1){
    write.table(cbind(as.character(pheno[,1]),notIDvalue),"keep.txt",quote=F,col.names=F,row.names=F)
  }else
    write.table(cbind(notIDvalue,as.character(pheno[,1])),"keep.txt",quote=F,col.names=F,row.names=F)
  plinkdataname<-get_filename(filenames,"bed")
  system("plink --noweb --bfile reportRx --keep keep.txt --make-bed --out reportRx")
  QCsummary<<-rbind(QCsummary,c("Missing Phenotype",parse_reportRx.log()))
}
#' Remove people with sex problems
#' 
#' Remove people with sex problems. This function takes no paramaters.
#' @keywords GWAS
#' @export
remove_sex_problems<-function(){
  reportRxGWASQCcount<<-reportRxGWASQCcount+1
  cat(paste0("\\section*{",reportRxGWASQCcount,": Identification of individuals with discordant sex information}"))  
  
  system("plink --noweb --bfile reportRx --check-sex --out reportRx")  
  sexcheck<-read.table("reportRx.sexcheck",header=T)
  SexProblems<-which(sexcheck$STATUS=="PROBLEM")
  if(length(SexProblems)==0){
    log<-parse_reportRx.log()
    QCsummary<<-rbind(QCsummary,c("Sex Problems",c(0,log[2],0,log[4])))
  }else if(length(SexProblems)==nrow(sexcheck)){
    warning("No sex information in genotype")
    log<-parse_reportRx.log()
    QCsummary<<-rbind(QCsummary,c("Sex Problems",c(log[2],0,0,log[4])))
  }else{
    write.table(sexcheck[SexProblems,c(1,2)],"SexProblems.txt",quote=F,row.names=F,col.names=F)
    system("plink --noweb --bfile reportRx --remove SexProblems.txt --make-bed --out reportRx")
    QCsummary<<-rbind(QCsummary,c("Sex Problems",parse_reportRx.log()))
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
#'  @keywords GWAS
#'  @export 
remove_missing<-function(snp.cutoff=0.05,individual.cutoff=0.05){
  reportRxGWASQCcount<<-reportRxGWASQCcount+1
  cat(paste0("\\section*{",reportRxGWASQCcount,": Identification of individuals/markers with a high missing rate}"))  
  
  system("plink --noweb --bfile reportRx --missing --out reportRx")
  lmiss<-read.table("reportRx.lmiss",stringsAsFactors=F,colClasses=c("character","character","numeric","numeric","numeric"),header=T,comment.char = "`")
  fully_missing_snps<-length(lmiss[,5][lmiss[,5]==1])
  missing_snps<-length(lmiss[,5][lmiss[,5]>=0.05])-fully_missing_snps
  log<-parse_reportRx.log()
  QCsummary<<-rbind(QCsummary,c("Fully Missing Snps",c(0,log[2],fully_missing_snps,log[4]-fully_missing_snps)))
  QCsummary<<-rbind(QCsummary,c(paste0(snp.cutoff*100,"% Missing Snps"),c(0,log[2],missing_snps,log[4]-fully_missing_snps-missing_snps)))
  system(paste("plink --noweb --bfile reportRx --geno",snp.cutoff,"--out reportRx --make-bed"))  
  #   system("plink --noweb --bfile reportRx --geno 0.05 --out reportRx")
  #   system("sed '303342q;d' reportRx.lmiss")
  #   system("sed '303343q;d' reportRx.lmiss")
  #   system("sed '303344q;d' reportRx.lmiss")
  #   system("sed '303345q;d' reportRx.lmiss")
  #   system("sed '303346q;d' reportRx.lmiss")
  #   system("awk 'NR==303342q{print;exit}' reportRx.lmiss")
  system(paste("plink --noweb --bfile reportRx --mind",individual.cutoff,"--out reportRx --make-bed"))
  log<-parse_reportRx.log(mind=T)
  QCsummary<<-rbind(QCsummary,c(paste0(individual.cutoff*100,"% missing individuals"),log))
  cat(paste0(ifelse(fully_missing_snps==0,"No",fully_missing_snps),
             " SNPs in the dataset are 100\\% missing",ifelse(fully_missing_snps==0,".",", which were removed from the data first."), 
             " Next SNPs with more than ",
             snp.cutoff*100,"\\% of individuals missing them were removed. ",ifelse(missing_snps==0,"No",missing_snps)," SNPs' call rate were below this cutoff.",
             " Next Individuals with more than ",individual.cutoff*100," \\% SNPs missing were removed. ",ifelse(log[1]==0,"No",log[1]),
             " individuals were removed under this criterion.\\\\"))
  
  system("plink --noweb --bfile reportRx --missing --out reportRx")
  imiss<-read.table("reportRx.imiss",stringsAsFactors=F,header=T,comment.char = "`")
  
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
#' @keywords GWAS
#' @export 
remove_relatives<-function(cutoff=0.25){
  reportRxGWASQCcount<<-reportRxGWASQCcount+1
  cat(paste0("\\section*{",reportRxGWASQCcount,": Identification of duplicated or related individuals}"),
      "To detect duplicated or related individuals, identity by state is calculated for each pair of individuals. If two individuals are identified as relatives, as defined by having IBS$\\geq$",cutoff," the one with the higher sample number was removed.") 
  
  #system("plink --noweb --bfile reportRx --genome --out reportRx")
  ibs = read.table("reportRx.genome",header=T,stringsAsFactors=F)
  ibs<-ibs[ibs$PI_HAT >= cutoff,c(1,2,3,4,10),drop=F]
  if(nrow(ibs)!=0){
    first_degree_relatives<-ibs[,c(3,4)]
    write.table(first_degree_relatives,"first_degree_relatives.txt",quote=FALSE,col.names=F,row.names=F)
    system("plink --noweb --bfile reportRx  --remove first_degree_relatives.txt --out reportRx --make-bed")
    cat("People Removed in this step are in bold and listed in the following table.")
    colnames(ibs)
    colnames(ibs)<-sanitizestr(colnames(ibs))
    ibs<-sapply(ibs,function(cols)sanitizestr(cols))
    ibs[,c(3,4)]<-sapply(ibs[,c(3,4)],lbld)
    print.xtable(xtable(ibs),include.rownames=F,table.placement="H",sanitize.text.function=identity)
    #QCsummary[7,-1]<<-parse_reportRx.log()
    QCsummary<<-rbind(QCsummary,c("Related Individuals",parse_reportRx.log()))
    
  }else{
    log<-parse_reportRx.log()
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
#' @export 
remove_heterozygocity<-function(SD=6){
  reportRxGWASQCcount<<-reportRxGWASQCcount+1
  cat(paste0("\\section*{",reportRxGWASQCcount,": Identification of individuals with outlying heterozygosity rate}"))  
  
  cat("The outliers in the heterozygosity rate might be indicative of DNA sample contamination or inbreeding.
      The following plot shows the heterozygosity rate of the individuals in our CRC data, and the horizontal dashed lines represent",SD,"standard deviations from the mean.\\\\")
  system("plink --noweb --bfile reportRx --het --out reportRx")
  system("plink --noweb --bfile reportRx --missing --out reportRx")
  
  
  imiss <- read.table("reportRx.imiss",h=T);
  imiss$logF_MISS = log10(imiss$F_MISS);
  het=read.table("reportRx.het",h=T);
  
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
    system("plink --noweb --bfile reportRx --remove hetero_remove.txt --out reportRx --make-bed")
    #QCsummary[7,-1]<<-parse_reportRx.log()  
    QCsummary<<-rbind(QCsummary,c(paste0("Heterozygocity Rate (-/+ ",SD,"SD)"),parse_reportRx.log()))
    cat("\n","The following people were removed for having heterozygosity rate",SD,"standard deviations form the mean:")
    print.xtable(xtable(hetero_remove),table.placement="H",include.rownames=F)
  }else{
    log<-parse_reportRx.log()
    QCsummary<<-rbind(QCsummary,c(paste0("Heterozygocity Rate (-/+ ",SD,"SD)"),0,log[2],0,log[4]))
    cat("No patients had heterozygosity rate",SD,"standard deviations form the mean")
  }
  
}
prepare_hapmap<-function(){
  #   system("plink --noweb --bfile reportRx --indep-pairwise 50 5 0.2 --out reportRxHM")
  # system("plink --noweb --bfile reportRxHM --extract reportRxHM.prune.in --recode --out reportRxHM") 
  #   system("plink --noweb --file hapmap3 --filter-founders --recode --make-bed --out hapmap")
  #   system("plink --noweb --bfile hapmap --extract reportRxHM.prune.in --make-bed --out hapmap")
  #   system("awk '{print $2}' hapmap.bim > commonsnp.txt")
  #   system("plink --noweb --file reportRxHM --extract commonsnp.txt --make-bed --out reportRxHM")
  #   system("plink --noweb --bfile reportRxHM --bmerge hapmap.bed hapmap.bim hapmap.fam --make-bed --out hapmap")
  #   system("plink --noweb --bfile hapmap1 --flip gwashapmap.missnp --make-bed --out hapmap2")
  # /plink --noweb --bfile hapmap2 --freq --out hapmapfreq
  # /plink --noweb --bfile gwas10 --freq --out gwasfreq
  
  system("plink --noweb --bfile reportRx --indep-pairwise 50 5 0.2 --out gwas8")
  system("plink --noweb --bfile reportRx --extract gwas8.prune.in --recode --out gwas9")
  system("plink --noweb --file hapmaporigional --filter-founders --recode --make-bed --out hapmap")
  #1457897 markers, 1198 founders 
  #extract gwas independent SNPs: 160500 SNPS totally
  system("plink --noweb --bfile ./hapmap --extract gwas8.prune.in --make-bed --out hapmap1")
  # there are 105865 SNPs
  system("awk '{print $2}' ./hapmap1.bim > ./commonsnp.txt")
  system("plink --noweb --file gwas9 --extract commonsnp.txt --make-bed --out gwas10")
  system("plink --noweb --bfile gwas10 --bmerge hapmap1.bed hapmap1.bim hapmap1.fam --make-bed --out gwashapmap")
  system("plink --noweb --bfile hapmap1 --flip gwashapmap.missnp --make-bed --out hapmap2")
  system("plink --noweb --bfile hapmap2 --freq --out hapmapfreq")
  system("plink --noweb --bfile gwas10 --freq --out gwasfreq")
  system("join -1 2 -2 2 hapmapfreq.frq gwasfreq.frq > mergefrq.txt")
  system("awk '{print $1,$5,$10}' ./mergefrq.txt > ./mafreq.txt")
  system("cat ./mafreq.txt | awk 'NR>1 {print $0 \" \" ($2-$3)}' ./mafreq.txt > mafcomp.txt")
  system("cat ./mafcomp.txt | awk '($4 > 0.2) || ($4 < -0.2) {print $0}' ./mafcomp.txt > ./snp_bigdif.txt")
  system("wc -l ./snp_bigdif.txt")
  # 8475 SNPs
  system("awk '{print $1}' ./snp_bigdif.txt > ./snpid_bigdif.txt")
  system("cat ./mafcomp.txt | awk '($2 > 0.5) {print $0}' ./mafcomp.txt > ./snp_nomaf1.txt")
  system("cat ./mafcomp.txt | awk '($3 > 0.5) {print $0}' ./mafcomp.txt > ./snp_nomaf2.txt")
  system("plink --noweb --bfile gwas10 --bmerge hapmap2.bed hapmap2.bim hapmap2.fam --exclude snpid_bigdif.txt --make-bed --out gwashapmap1")
  # There's still one SNP: rs2862907 not able to be merged because of mis-matching strand, 
  system("plink --noweb --bfile gwas10 --exclude gwashapmap1.missnp --make-bed --out gwas11")
  system("plink --noweb --bfile gwas11 --exclude snpid_bigdif.txt --make-bed --out gwas12")
  system("plink --noweb --bfile hapmap2 --exclude gwashapmap1.missnp --make-bed --out hapmap3")
  system("plink --noweb --bfile hapmap3 --exclude snpid_bigdif.txt --make-bed --out hapmap4")
  system("plink --noweb --bfile gwas12 --bmerge hapmap4.bed hapmap4.bim hapmap4.fam --make-bed --out gwashapmap2")
  # 97389SNPs
  system("plink --noweb --file gwas9 --extract commonsnp.txt --make-bed --out ")
}
#' End GWAS QC
#' 
#' Outputs summary of QC process. There are no paramaters.
#' @export
end_QC<-function(){
  reportRxGWASQCcount<<-reportRxGWASQCcount+1
  cat(paste0("\\section*{",reportRxGWASQCcount,": Summary of QC}"))  
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