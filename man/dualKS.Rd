\name{dualKS}
\Rdversion{1.1}
\alias{dualKS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dual knowledge structure
}
\description{
\code{dualKS} computes the dual knowledge structure of the given \code{KS}.
}
\usage{
dualKS(KS)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{KS}{
Knowledge structure (cf. \link[skills]{skills-package})
}
}
\details{
To each knowledge structure there exists the dual knowledge structure which consists of the following states:
To each state S of the knowledge structure the associated dual knowledge state contains exactly those items of Q (the set of all items)
which are not in S.
}
\value{
\code{dualKS} returns the dual knowledge structure in the same representation as the given \code{KS}.
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
KS1 = set(set(), set(3), set(5), set(2,5), set(1,3,4,5), set(1,2,3,4,5))
KS1
dualKS(KS1)

KS2=matrix(c(0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,1,0,0,1,1,0,1,1,1,1,1,1,1,1),
    nrow = 5)
rownames(KS2) = 1:5
KS2
dualKS(KS2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{knowledge structure}
