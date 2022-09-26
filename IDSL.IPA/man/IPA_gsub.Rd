\name{IPA_gsub}
\alias{IPA_gsub}
\title{
IPA gsub
}
\description{
This function is mostly `gsub` command of R, but can replace a vector of strings.
}
\usage{
IPA_gsub(pattern, replacement, x, ignore.case = FALSE, perl = FALSE,
fixed = FALSE, useBytes = FALSE)
}
\arguments{
  \item{pattern}{
pattern which can be a vector of strings. `gsub` command of R can only take one pattern string
}
  \item{replacement}{
replacement
}
  \item{x}{
a string
}
  \item{ignore.case}{
ignore.case
}
  \item{perl}{
perl
}
  \item{fixed}{
fixed
}
  \item{useBytes}{
useBytes
}
}
\value{
A substituted string
}
