\name{KSbase}
\alias{KSbase}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Base of a knowledge structure
}
\description{
\code{KSbase} computes the base of a knowledge structure \code{KS} if \code{KS} is a knowledge space (\link[skills]{is.KSpace}).
}
\usage{
KSbase(KS)
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
set of states in \code{KS} which form the base of \code{KS}
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
\code{\link[skills]{is.KSpace}} checks if \code{KS} is a knowledge space (closed under union) and hence has a base.
If so, the base is formed by the collection of the atoms (\link[skills]{KSatoms}) of the knowledge space.
}
\examples{
s1 = set(set(), set(1,3), set(2,3), set(1,2,3), set(1,2,4), set(1,2,3,4))
KSbase(as.KS(s1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{base}
\keyword{knowledge structure}
\keyword{knowledge space}
