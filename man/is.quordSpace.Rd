\name{is.quordSpace}
\Rdversion{1.1}
\alias{is.quordSpace}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Checks if knowledge structure is a quasi-ordinal space
}
\description{
The function \code{is.quordSpace} checks if a knowledge structure is a quasi-ordinal space,
that means if it is closed under union as well as under intersection.
}
\usage{
is.quordSpace(KS, domain = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{KS}{
A knowledge structure in 'set', 'matrix' or 'gset' representation (cf. \link[skills]{skills-package})
}
  \item{domain}{
Domain of the knowledge structure. If \code{NULL} it is the union of all items contained in \code{KS}.
}
}
\details{
\code{is.quordSpace} checks if a knowledge structure is a quasi-ordinal space,
that means if it is closed under union and under intersection. Therefore both closures
are derived and have to be identical to the knowledge structure itself.
}
\value{
logical;\cr
attribute "missing": states which are contained in one of the closures but not in KS
}
\references{
Doignon, J.-P., Falmagne, J.-C. (1999), \emph{Knowledge Spaces}. Berlin, Heidelberg and New York: Springer-Verlag.
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
Both \code{\link[skills]{is.KSpace}} and \code{\link[skills]{is.CSpace}} have to hold for this check.
}
\examples{
KS = set(set(), set("a"), set("b"), set("a","b"), set("a","c"),
    set("a","b","c"), set("a","b","c","d"), set("a","b","c","e"),
    set("a","b","c","d","e"), set("a","b","c","d","e","f"))
is.quordSpace(KS)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{quasi ordinal space}
\keyword{knowledge structure}
