\name{gss}
\alias{gss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to do Golden Section Search
}
\description{
A bracketing function to peform golden-section search.
}
\usage{
gss(f, lower, upper, tol = 0.0001, maxiter = 1000,
                showiter = FALSE, makinfo = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{the function to be optimized}
  \item{lower}{the lower end point of the interval to be searched.}
  \item{upper}{the upper end point of the interval to be searched.}
  \item{tol}{the desired accuracy.}
  \item{maxiter}{maximum iteration taken to run the bracketing function. Should
              be defined to prevent the case of infinite iteration.}
  \item{showiter}{give description or information for each iteration.}
  \item{makinfo}{give general information on function's results.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 A list of minimum (or maximum), total number of iteration, and objective
 which give the location of the minimum (or maximum). The description of
 result is also given.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Insan R. Adiwibowo
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
func <- function(x) x^4 - 8*x^2 + 2*x
res <- gss(f = func, lower = -5, upper = 5, showiter = TRUE)
res
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
