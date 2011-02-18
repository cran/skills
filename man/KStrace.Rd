\name{KStrace}
\alias{KStrace}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Trace of a knowledge structure
}
\description{
\code{KStrace} computes the trace of a knowledge structure \code{KS}.
}
\usage{
KStrace(KS, subs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{KS}{
Knowledge structure in 'set' or 'matrix' representation (cf. \link[skills]{skills-package}); 'gset' for empirical knowledge structures
}
  \item{subs}{
Subset of the domain for which the trace should be computed
}
}
\details{
The trace of a knowledge structure \code{KS} is gained by restricting \code{KS} to a subset of items (\code{subs}). If the
passed parameter \code{subs} contains the domain of \code{KS} and thus consists of all items contained in the domain or even
additional items the returned knowledge structure is equal to the passed \code{KS}. If \code{subs} does not contain any item
contained in the domain the empty set is returned.\cr
If an empirical knowledge structure is passed the according frequencies of the states are summed up as well.
}
\value{
Knowledge structure (cf. \link[skills]{skills-package})
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
s1 = set(set(), set(1,3), set(2,3), set(1,2,3), set(1,2,4), set(1,2,3,4))
s1
KStrace(s1, subs = set(2,3))
KStrace(s1, subs = set(2))
KStrace(s1, subs = set(1,2,3,4,5))

eKS = gset(s1, memberships = 1:6)
eKS
KStrace(eKS, subs = set(2,3))
KStrace(eKS, subs = set(2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{trace}
\keyword{knowledge structure}
