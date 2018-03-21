\name{gradient}
\alias{gradient}
\title{
Function to do Gradient Descent Optimization
}
\description{
A bracketing function to peform gradient descent algorithm.
}
\usage{
gradient(f, start, gamma = 0.01, tol = 0.0001, maxiter = 1000,
                     showiter = FALSE, makinfo = TRUE)
}
\arguments{
  \item{f}{The function to be optimized}
  \item{start}{numeric. Starting reference point.}
  \item{gamma}{numeric. Step size multiplier.}
  \item{tol}{Desired accuracy level.}
  \item{maxiter}{Maximum iteration taken to run the bracketing function. Should
              be defined to prevent the case of infinite iteration.}
  \item{showiter}{logical. If TRUE give description or information for each iteration.}
  \item{makinfo}{logical. If TRUE give general information on function's results.}
}
\details{
}
\value{
A list of minimum (or maximum), total number of iteration, and objective
which give the location of the minimum (or maximum). The description of
result is also given.
}
\references{

}
\author{
Insan R. Adiwibowo
}
\note{

}

\seealso{

}
\examples{
func <- function(x) x^4 - 8*x^2 + 2*x
# See the result of the function, for example, with starting point -5
# to find the nearest minima from the steep. We also want to see the result for
# each iteration.
gradient(f = func, start = -5, showiter = TRUE)
}

