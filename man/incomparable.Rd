\name{incomparable}
\Rdversion{1.1}
\alias{incomparable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Checks if a set of sets is incomparable
}
\description{
Given a set of states of a knowledge structure the function checks if
these states are incomparable, i.e. no state is a subset of another.
}
\usage{
incomparable(KS)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{KS}{
Subset of a knowledge structure in '\link[sets]{set} of sets' or 'matrix' representation (cf. \link[skills]{skills-package}).
}
}
\details{
Given a set of states \code{incomparable} checks whether or not any state is contained in another. The empty state is not considered.
}
\value{
logical
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Angela Haidinger \email{angela.ulrike.haidinger@student.uni-augsburg.de},\cr
Ali Uenlue \email{uenlue@statistik.tu-dortmund.de}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
set1 = set(set(), set(1), set(3), set(2,5))
incomparable(set1)

set2 = set(set(), set(1), set(3), set(2,3,5))
incomparable(set2)

mat = matrix(c(0,0,0,0,1,0,1,0,0,1,1,0,1,1,0,1), nrow = 4)
rownames(mat) = 1:4
mat
incomparable(mat)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{incomparable}
\keyword{subset}
