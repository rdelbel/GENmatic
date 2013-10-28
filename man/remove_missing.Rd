\name{remove_missing}
\alias{remove_missing}
\title{Remove snps and people with too many missing}
\usage{
  remove_missing(snp.cutoff = 0.05,
    individual.cutoff = 0.05)
}
\arguments{
  \item{snp.cutoff}{numeric corresponding to cutoff for
  missing snps.  Default is 0.05}

  \item{individual.cutoff}{numeric corresponding to cutoff
  for missing individuals. Default is 0.05}
}
\description{
  First removed 100% missing snps, then missing snps, then
  issing idividuals.
}
\keyword{GWAS}

