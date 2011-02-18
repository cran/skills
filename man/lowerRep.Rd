\name{lowerRep}
\Rdversion{1.1}
\alias{lowerRep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Coefficient of reproducibility for Guttman-scales (lower bounds)
}
\description{
\code{lowerRep} computes a coefficient of reproducibility for
Guttman-scales. This coefficient is a measure for the overall divergence
of empirical states to their lower bound (\link[skills]{lb_eKS})
in the theoretical knowledge structure.
}
\usage{
lowerRep(eKS, tKS, model)
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
Model for Skill Assignment (either "disjunctive" or "conjunctive")
}
}
\details{
\code{lowerRep} counts the number of element exchages which are necessary to gain from each empirical state
its lower bound (cf. \link[skills]{lb_eKS}) and then returns the proportion of this number in the number of possible element
exchanges at all.
}
\value{
proportion (between 0 and 1)
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
\code{\link[skills]{lowerDist}} computes the necessary single distances.
}
\examples{
tKS_d = set(set(),set(2), set(1,3,4), set(1,2,3,4), set(1,2,4,5),
    set(1,2,3,4,5))
eKS_d = gset(set(set(), set(2,4), set(1,3), set(1,2,3,5), set(1,2,3,4,5)),
    memberships = c(1,2,3,4,5))
lowerRep(eKS_d, tKS_d, "disjunctive")

tKS_c = set(set(), set(3), set(5), set(2,5), set(1,3,4,5), set(1,2,3,4,5))
eKS_c = gset(set(set(), set(3), set(3,5), set(2,4,5), set(1,3,5),
    set(1,2,3,4,5)), memberships = c(1,2,3,4,5,6))
lowerRep(eKS_c, tKS_c, "conjunctive")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{empirical}
\keyword{knowledge structure}
