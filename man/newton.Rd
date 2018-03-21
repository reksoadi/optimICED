\name{newton}
\alias{newton}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to do Newton optimization
}
\description{
An iterative function to peform Newton iteration algorithm
}
\usage{
newton(f, start, tol = 0.0001, maxiter = 1000, showiter = FALSE,
                   makinfo = TRUE)
}
\arguments{
  \item{f}{The function to be optimized}
  \item{start}{numeric. Starting reference point.}
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
# See the result of the function, for example, with starting point -5. We also
# want to see the result for each iteration.
newton(f = func, start = -5, showiter = TRUE)
}

