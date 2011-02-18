\name{is.KS}
\Rdversion{1.1}
\alias{is.KS}
\alias{is.KS.set}
\alias{is.KS.matrix}
\alias{is.KS.gset}
\title{
Check if KS is a knowledge structure
}
\description{
Checks if KS contains the empty set and the whole domain and thus satisfies the definition of a knowledge structure.
}
\usage{
is.KS(x, domain = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Knowledge structure in set, gset or matrix representation (cf. \link[skills]{skills-package})
}
  \item{domain}{
Domain of interested knowledge structure (as \link[sets]{set} of all items)
}
  \item{\dots}{
further arguments
}
}
\details{
If a \code{domain} is passed which is smaller than the union of all occuring items the returned value declares if the object
 defined by the trace on this domain is a knowledge structure (contains the empty set and the whole domain as set).
}
\value{
logical
}
\references{
Doignon, J.-P. and Falmagne, J.-C. (1999), \emph{Knowledge Spaces}. Berlin, Heidelberg and New York: Springer-Verlag.
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
KS1 = set(set(), set(1,3), set(2,3), set(1,2,3), set(1,2,4), set(1,2,3,4))
is.KS(KS1)
is.KS(KS1, domain = set(1,2,3))

KS2 = KS1 - set(set())
is.KS(KS2)

KS3 = matrix(c(0,0,0,0,1,0,1,0,0,1,1,0,1,1,1,0,1,1,0,1,1,1,1,1), ncol = 6) 
# equals KS1
is.KS(KS3)
KS3
is.KS(KS3, domain = set("a", "b", "c", "d"))
# TRUE because the convert function is applied which names the items
# according to the first letters if not specified different

KS4 = gset(KS1, memberships = 1:6)
is.KS(KS4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{knowledge structure}
