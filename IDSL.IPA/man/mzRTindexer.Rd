\name{mzRTindexer}
\alias{mzRTindexer}
\title{
m/z - RT Indexer
}
\description{
This function locate the closeset pair of a reference (m/z - RT) pair in a 2-D grid of `m/z` and `RT` vectors.
}
\usage{
mzRTindexer(MZvec, RTvec, MZref, RTref, MZtolerance, RTtolerance)
}
\arguments{
  \item{MZvec}{
m/z vector
}
  \item{RTvec}{
RT vector
}
  \item{MZref}{
a reference m/z
}
  \item{RTref}{
a reference RT
}
  \item{MZtolerance}{
m/z tolerance
}
  \item{RTtolerance}{
RT tolerance
}
}
\note{
This function return NULL in case no match is detected.
}
\value{
index of closest pair to the reference (m/z - RT) pair
}
