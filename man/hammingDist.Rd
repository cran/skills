\name{hammingDist}
\Rdversion{1.1}
\alias{hammingDist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hamming distance between two states
}
\description{
\code{hammingDist} computes the hamming distance between two states
of a knowledge structure. This distance counts the itmes which are
contained in exactly one of those states.
}
\usage{
hammingDist(s1, s2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s1}{
State of a knowledge structure either as 0-1-vector (with the same length as s2)
or as \link[sets]{set}
}
  \item{s2}{
State of a knowledge structure either as 0-1-vector (with the same length as s1)
or as \link[sets]{set}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
integer, counting the different items
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
s1 = c(1,1,1,0)
s2 = c(0,1,1,1)
hammingDist(s1, s2)

s1 = set(1,2,4,6)
s2 = set(2,3,4,7)
hammingDist(s1, s2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hamming}
\keyword{set}
\keyword{distance}
