\name{KSdomain}
\alias{KSdomain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Domain of a knowledge structure
}
\description{
\code{KSdomain} computes the domain of a knowledge structure.
}
\usage{
KSdomain(KS)
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
sorted set of all items
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
KSdomain(as.KS(s1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{domain}
\keyword{knowledge structure}
