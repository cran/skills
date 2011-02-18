\name{permute}
\Rdversion{1.1}
\alias{permute}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computation of all permutations of a vector
}
\description{
\code{permute} computes all permutations of a given part of a vector.
}
\usage{
permute(v, beg = 1, last = length(v))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v}{
vector
}
  \item{beg}{
First digit of that part of \code{v} which should be permuted.
}
  \item{last}{
Last digit of that part of \code{v} which should be permuted.
}
}
\details{
All vectors which contain a permutation of that part of \code{v} between \code{beg} and
\code{last} are computed.
}
\value{
A matrix whose columns contain the computed permutations.
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
v = 1:4
permute(v)

w = 1:5
permute(w, beg = 2, last = 4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{permutation}
