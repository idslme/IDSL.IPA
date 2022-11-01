\name{mzRTindexer}
\alias{mzRTindexer}
\title{
m/z - RT Indexer
}
\description{
This function locate the closest pair of a reference (m/z - RT) pair in a 2-D grid of `m/z` and `RT` vectors.
}
\usage{
mzRTindexer(MZvec, RTvec, MZref, RTref, massAccuracy, RTtolerance)
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
  \item{massAccuracy}{
m/z tolerance
}
  \item{RTtolerance}{
RT tolerance
}
}
\note{
This function returns NULL in case no match is detected.
}
\value{
index of closest pair to the reference (m/z - RT) pair
}
