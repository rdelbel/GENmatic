\name{makeRplink}
\alias{makeRplink}
\title{Make Rplink.txt file}
\usage{
  makeRplink(outcome, covariates = NULL, strata = NULL,
    mcore = F, fout = NULL)
}
\arguments{
  \item{outcome}{length 2 vector corresponding to time and
  status}

  \item{covariates}{vector corresponding to confounding
  factors. Default (no covariates) is NULL}

  \item{strata}{vector corresponding to stratification
  variable. Default (no strata) is NULL.}

  \item{mcore}{Dont use this}

  \item{fout}{string corresponding to the directory you
  want to output the Rplink.R file. Default (write to
  working directory) is NULL.}
}
\description{
  Makes Rplink.txt file for fash coxph gwas. #'Any element
  of the outcome,strata, or null vectors will correspond to
  either the column index (not including the 2 id columns)
  or column name of the corresponding covaraite in the
  covar.txt file. Any element of these vectors will be
  treated as an index if it consists of the digits 0-9 and
  will be considered a column name otherwise. Even if you
  input an element as "101" (with quotes) it will be
  treated as an index to the 101th column not the column
  name 101. Only one stratification variable can be used.
  If you want to use more than one, create a variable that
  incorporates all of the strata in one column.
}

