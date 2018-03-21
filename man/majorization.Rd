\name{majorization}
\alias{majorization}
\title{
Function to do Iterative Majorization
}
\description{
An iterative function to peform iterative majorization algorithm.
}
\usage{
majorization(f, start, tol = 0.0001, maxiter = 1000,
                         showiter = FALSE, makinfo = TRUE)
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
# Before doing majorization technique, define the majorization function g(x, y)
# and its minimum point.
func <- function(x, option = 0) {
  # First, define our basic function h(x).
  hx <- function(x) x^4 - 8*x^2 + 2*x

  # Second function constructed by using the fact:
  # (x - y)^2 =< 0 <=> -x^2 = y^2 - 2*x*y.
  # So the g(x, y) function is resulted as follow:
  gxy <- function(x, y) x^4 + 8*y^2 - 16*x*y + 2*x

  # Third function finds the minimum of g(x, y) by doing partial derivation on x
  # in which g'(x, y) = 0, resulting x = (4*x - 0.5)^1/3.
  min.gxy <- function(x) sign(4*x - 0.5)*abs(4*x - 0.5)^(1/3)

  # Define the choice
  if (option == 0) return(hx(x))
  if (option == 1) return(gxy(x, y))
  if (option == 2) return(min.gxy(x))
}

# See the result of the function, for example, with starting point -5. We also
# want to see the result for each iteration.
majorization(f = func, start = -5, showiter = TRUE)
}

