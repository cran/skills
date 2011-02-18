\name{permtest}
\Rdversion{1.1}
\alias{permtest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Permutation test of random assignment of items to skill sets
}
\description{
The function \code{permtest} computes the p-value of a permutation test.
}
\usage{
permtest(eKS, tKS, model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eKS}{
Empirical knowledge structure (cf. \link[skills]{skills-package})
in "\link[sets]{gset} of sets" - representation with
memberships equal to observed frequencies.
}
  \item{tKS}{
Theoretical knowledge structure (cf. \link[skills]{skills-package})
in "\link[sets]{set} of sets" - represtentation
}
  \item{model}{
Model for Skill Assignment (either "disjunctive" or "conjunctive")
}
}
\details{
The function \code{permtest} computes the p-value of a permutation test with the associated null hypothesis:
\code{w_cind(eKS,tKS)} is the result of a random assignment of problems to skill sets.
}
\value{
p-value
}
\references{
Duentsch, I., Gediga, G. (2002), \emph{Skill Set Analysis in Knowledge Structures}. British Journal of Mathematical and
Statistical Psychology, 55(2), 361 - 384.
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
The single weighted consistency indeces are computed by \code{\link[skills]{w_cind}}.
}
\examples{
tKS1 = set(set(), set(2), set(1,3,4), set(1,2,3,4), set(1,2,4,5),
    set(1,2,3,4,5))
eKS1 = gset(set(set(), set(2,3), set(1,3,4), set(1,2,3,4), set(1,2,4,5),
    set(1,2,3,4,5)), memberships = c(1,2,3,4,5,6))
permtest(eKS1, tKS1, "disjunctive")

tKS2 = set(set(), set(3), set(5), set(2,5), set(1,3,4,5), set(1,2,3,4,5))
eKS2 = gset(set(set(), set(3), set(3,5), set(2,4,5), set(1,3,5),
    set(1,3,4,5), set(1,2,3,4,5)),
    memberships = c(80,96,4,7,14,74,55))
permtest(eKS2, tKS2, "conjunctive")

eKS3 = as.KS(as.gset(pisa))
tKS3 = as.KS(as.relation(eKS3, empirical = TRUE, v = 1))
permtest(eKS3, tKS3, "disjunctive")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{permutation}
\keyword{empirical}
\keyword{knowledge structure}
