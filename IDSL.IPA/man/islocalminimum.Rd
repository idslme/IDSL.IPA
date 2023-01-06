\name{islocalminimum}
\alias{islocalminimum}
\title{
islocalminimum
}
\description{
This function returns indices of local minimum points on a curve.
}
\usage{
islocalminimum(y)
}
\arguments{
  \item{y}{
is a vector of y values.
}
}
\value{
A vector in the same size of the vector 'y'. Local minimum arrays represented by -1.
}
\examples{
data(chromatogramMatrix)
int <- chromatogramMatrix$smoothChromatogram
islocalminimum(int)
}
