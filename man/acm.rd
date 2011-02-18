\name{acm}
\Rdversion{1.1}
\alias{acm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adjusted consistency measure
}
\description{
\code{acm} adjusts \link[skills]{w_cind} by considering the expected value of a random assignment of items to skill sets.
}
\usage{
acm(eKS, tKS, model)
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
in "\link[sets]{set} of sets" - representation
}
  \item{model}{
Model, either "disjunctive" or "conjunctive",
for the computation of the skill assignment (cf. \link[skills]{SkillAss})
}
}
\details{
\code{acm} is an adjusted measurement of consistency. The function computes the expected value E of the weighted consistency indices
for the permutations of skill sets to items (this corresponds to a permutation of the item labels, the structure of the knowledge
structure remains the same) and then adjusts the weighted consistency index \link[skills]{w_cind} as follows: acm = (w-E)/(1-E).\cr
Returned values near or even under 0 mean that the fit of \code{tKS} is only as well as or even worse than the mean of
all permutations.
}
\value{
numeric value in [-Inf,1]
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
\code{\link[skills]{SkillAss}} for computing the skill assignment underlying the theoretical knowledge structure,\cr
\code{\link[skills]{w_cind}} for the unadjusted weighted consistency index.
}
\examples{
tKS = set(set(), set(3), set(5), set(2,5), set(1,3,4,5), set(1,2,3,4,5))
eKS = gset(set(set(), set(3), set(3,5), set(2,4,5), set(1,3,5),
    set(1,3,4,5), set(1,2,3,4,5)),
    memberships = c(80,96,4,7,14,74,55))
acm(eKS, tKS, "conjunctive")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{consistency}
\keyword{empirical}
\keyword{knowledge structure}
