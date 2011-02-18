\name{item_cind}
\Rdversion{1.1}
\alias{item_cind}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Item consistency index
}
\description{
\code{item_cind} computes the consistency indices for omitting
single items.
}
\usage{
item_cind(eKS, tKS, model)
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
Given a skill assignment, specified through \code{tKS} and \code{model}, \code{item_cind} omits one item at a time in the skill
assignment, defining the theoretical knowledge structure, as well as in the empirical knowledge structure and computes the weighted
consistency index (cf. \link[skills]{w_cind}) therefore.\cr
All values will be higher than the original one because both knowledge structures, \code{eKS} and \code{tKS}, are restricted to
traces and thus the consistency is higher.
If the resulting value differs much from the original one this is an indicator for a bad model fitting with respect to this item.
}
\value{
vector of weighted consistency indices, names according to omitted items
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
\code{\link[skills]{w_cind}} computes the single weighted consistency indeces.
}
\examples{
tKS_d = set(set(), set(2), set(1,3,4), set(1,2,3,4), set(1,2,4,5),
    set(1,2,3,4,5))
eKS_d = gset(set(set(), set(2,4), set(1,3), set(1,2,3,5), set(1,2,3,4,5)),
    memberships = c(1,2,3,4,5))
item_cind(eKS_d, tKS_d, "disjunctive")

eKS = as.KS(as.gset(pisa))
tKS = as.KS(as.relation(eKS, empirical = TRUE, v = 1))
item_cind(eKS, tKS, "disjunctive")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{consistency}
