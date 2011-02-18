\name{as.KS}
\Rdversion{1.1}
\alias{as.KS}
\alias{as.KS.set}
\alias{as.KS.matrix}
\alias{as.KS.relation}
\alias{as.KS.gset}
\title{
Creating knowledge structures
}
\description{
Creates an object of class KS
}
\usage{
as.KS(x, type = "set", space = TRUE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Knowledge structure in one of the described forms (\link[skills]{skills-package})
}
  \item{type}{
Type of knowledge structure which should be returned ("set" or "matrix")
}
  \item{space}{
logical; only has an effect if if \code{x} is a relation; then if \code{TRUE} the quasi ordinal space defined by this relation
is returned, otherwise the minimal knowledge structure.
}
  \item{\dots}{
further arguments
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Knowledge structure of class KS
}
\references{
%% ~put references to the literature/web site here ~
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

}
\examples{
KS1 = as.KS(set(set(), set(1,3), set(2,3), set(1,2,3), set(1,2,4),
                              set(1,2,3,4)))
class(KS1)

KS2=matrix(c(0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,1,0,0,1,1,0,1,1,1,1,1,1,1,1),
    nrow = 5)
as.KS(KS2, type = "matrix")
as.KS(KS2, type = "set")

KS3 = endorelation(graph = set(tuple(1,1), tuple(2,2), tuple(3,3),
    tuple(4,4), tuple(1,2), tuple(2,3), tuple(1,3), tuple(1,4)))
as.KS(KS3)
as.KS(KS3, space = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{knowledge structure}
