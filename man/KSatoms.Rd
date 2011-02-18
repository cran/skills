\name{KSatoms}
\alias{KSatoms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Atoms of a knowledge structure
}
\description{
\code{KSatoms} computes the atoms of a knowledge structure.
}
\usage{
KSatoms(KS)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{KS}{
Knowledge structure in 'set' or 'matrix' representation (cf. \link[skills]{skills-package})
}
}
\details{
For each item the atoms are the (with respect to subset) minimal elements of the knowledge structure which contain the item itself.
}
\value{
list: set of atoms for every item in the domain (cf. \link[skills]{KSdomain})
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
KSatoms(as.KS(s1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{atom}
\keyword{knowledge structure}
