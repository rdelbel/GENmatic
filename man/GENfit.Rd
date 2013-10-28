\name{GENfit}
\alias{GENfit}
\title{Fit genetic models}
\usage{
  GENfit(outcome, covariates = NULL, strata = NULL,
    type = "coxph", mcore = F, fout = NULL,
    bfilename = "GENmatic")
}
\arguments{
  \item{outcome}{length 2 vector corresponding to time and
  status}

  \item{covariates}{vector corresponding to confounding
  factors. Default (no covariates) is NULL}

  \item{strata}{vector corresponding to stratification
  variable. Default (no strata) is NULL.}

  \item{type}{Type of analysis to do. Currently can only do
  coxph.}

  \item{mcore}{Dont use this}

  \item{fout}{string corresponding to the directory you
  want to output the Rplink.R file. Default (write to
  working directory) is NULL.}
}
\description{
  Fit genetic models. Currently can only fit coxph models
}

