\name{remove_heterozygocity}
\alias{remove_heterozygocity}
\title{Remove indviduals with large heterozygocity}
\usage{
  remove_heterozygocity(SD = 6)
}
\arguments{
  \item{SD}{numeric corresponding to the numer of SD away
  from the mean we will choose as cutoff. Default is 6.}
}
\description{
  Removes all individuals with heterozogocity more than
  cutoff SD away form the mean.
}
\keyword{GWAS}

