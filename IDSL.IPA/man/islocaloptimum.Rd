\name{islocaloptimum}
\alias{islocaloptimum}
\title{
islocaloptimum
}
\description{
This function returns indices of local minimum and maximum points on a curve.
}
\usage{
islocaloptimum(y)
}
\arguments{
  \item{y}{
is a vector of y values.
}
}
\value{
A vector in the same size of the vector 'y'. Local minimum and maximum arrays represented by -1 and +1, respectively.
}
\examples{
data(chromatogram_builder)
int <- chromatogram_builder$SmoothedChromatogram
islocaloptimum(int)
}