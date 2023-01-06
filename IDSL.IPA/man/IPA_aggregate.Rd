\name{IPA_aggregate}
\alias{IPA_aggregate}
\title{
aggregation method for the IDSL.IPA modules
}
\description{
This module is to optimize the `indexVec` variable by removing elements that have redundant `idVec` numbers.
}
\usage{
IPA_aggregate(idVec, variableVec, indexVec, targetVar)
}
\arguments{
  \item{idVec}{
a vector of id numbers. Repeated id numbers are allowed.
}
  \item{variableVec}{
a vector of variable of the interest such as RT, m/z, etc.
}
  \item{indexVec}{
a vector of indices
}
  \item{targetVar}{
the targeted value in `variableVec`
}
}
\value{
a clean indexVec after removing repeated `idVec`.
}
