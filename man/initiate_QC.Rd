\name{initiate_QC}
\alias{initiate_QC}
\title{Initialize GWAS Quality Control}
\usage{
  initiate_QC()
}
\description{
  Looks in current directory for dirty.bed dirty.fam
  dirty.bin. Reads these files to get initial data. Saves
  new files as reportRx.bed/fam/bin so that we do not
  modify origional data Sets up 2 secret global variables
  for the rest of the QC process. This function takes no
  paramaters
}
\keyword{GWAS}

