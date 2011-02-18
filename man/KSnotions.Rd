\name{KSnotions}
\alias{KSnotions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Notions of a knowledge structure
}
\description{
\code{KSnotions} computes the notions of a knowledge structure.
}
\usage{
KSnotions(KS)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{KS}{
Knowledge structure in 'set' or 'matrix' representation (cf. \link[skills]{skills-package})
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
list of different notions
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
The disciminative reduction which is a knowledge structure which corresponds to the knowledge structure \code{KS}
but is reduced by its equally informative items, the notions, is computed by \code{\link[skills]{discred}}.
}
\examples{
s = set(set(), set(1,2), set(1,2,3), set(1,2,4), set(1,2,3,4))
KSnotions(as.KS(s))

re = endorelation(graph = set(tuple(1,1), tuple(1,2), tuple(2,1),
    tuple(2,2), tuple(3,3)))
relation_incidence(re)
KSnotions(as.KS(re))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{notion}
\keyword{knowledge structure}
