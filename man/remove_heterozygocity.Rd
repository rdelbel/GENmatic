% Generated by roxygen2 (4.0.0): do not edit by hand
\name{remove_heterozygocity}
\alias{remove_heterozygocity}
\title{Remove indviduals with large heterozygocity}
\usage{
remove_heterozygocity(SD = 6, pd = "")
}
\arguments{
\item{SD}{numeric corresponding to the numer of SD away
from the mean we will choose as cutoff. Default is 6.}

\item{pd}{directory where plink is. Do not need to specify if plink is in your path
or working directory}
}
\description{
Removes all individuals with heterozogocity more than
cutoff SD away form the mean.
}
\keyword{GWAS}

