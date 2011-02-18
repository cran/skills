\name{discred}
\alias{discred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Discriminative reduction of a knowledge structure
}
\description{
Computes the discriminative reduction of a knowledge structure.
}
\usage{
discred(KS, concat = ":")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{KS}{
Knowledge structure (cf. \link[skills]{skills-package})
}
  \item{concat}{
symbol for the concatenation of notions
}
}
\details{
The discriminative reduction of a knowledge structure consists of the same information as the knowledge structure itself
only the notions (containing the equal informative items) are combined to one item.
}
\value{
Knowledge structure (cf. \link[skills]{skills-package}) of the same type as the given \code{KS}
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
\code{\link[skills]{KSnotions}} computes the notions which are joined together such that the knowledge structure is discriminative,\cr
\code{\link[skills]{plothasse}} calls this function since the displayed quasi order has to be a partial order.
}
\examples{
KSs = set(set(), set(3), set(5), set(2,5), set(1,3,4,5), set(1,2,3,4,5))
discred(KSs)

KSm = convert(KSs)
discred(KSm)

KSgs = gset(KSs, memberships = 1:6)
discred(KSgs)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{notion}
\keyword{knowledge structure}
