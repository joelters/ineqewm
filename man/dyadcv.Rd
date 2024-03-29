% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dyadcv.R
\name{dyadcv}
\alias{dyadcv}
\title{dyadic ML cross-validation}
\usage{
dyadcv(
  X,
  Y,
  f,
  ML = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB", "Logit_lasso"),
  Kcv = 5
)
}
\arguments{
\item{X}{is a dataframe containing all the features on which the model was estimated}

\item{Y}{is a vector containing the labels for which the model was estimated}

\item{f}{function of Yi and Yj defining the dependent variable}

\item{ML}{string vector specifying which machine learners to use}

\item{Kcv}{number of folds}
}
\value{
list containing ML attaining minimum RMSE and RMSE
}
\description{
dyadic ML cross-validation
}
