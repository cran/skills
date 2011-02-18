\name{is.KSpace}
\Rdversion{1.1}
\alias{is.KSpace}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Checks if knowledge structure is a knowledge space
}
\description{
The function \code{is.KSpace} checks whether a knowledge structure is a knowledge space or not,
that means if it is closed under union.
}
\usage{
is.KSpace(KS, domain = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{KS}{
A knowledge structure in 'set', 'matrix' or 'gset' representation (cf. \link[skills]{skills-package})
}
  \item{domain}{
Domain of the knowledge structure. If \code{NULL} it is set to the union of all items contained in \code{KS}.
}
}
\details{
\code{is.KSpace} checks if a knowledge structure is a knowledge space,
that means if it is closed under union. Therefore the closure under union
is derived and has to be identical to the knowledge structure itself.
}
\value{
logical;\cr
attribute "missing": states which are contained in the closure but not in KS
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
\code{\link[sets]{closure}} with parameter \code{operation = "union"} computes the relevant closure
}
\examples{
KS1 = set(set(), set(1,3), set(2,3), set(1,2,3), set(1,2,4), set(1,2,3,4))
is.KSpace(KS1)

KS2 = dualKS(KS1)
is.KSpace(KS2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{union}
\keyword{knowledge space}
\keyword{knowledge structure}
