\name{skill_cind}
\Rdversion{1.1}
\alias{skill_cind}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Consistency index for skill reduct
}
\description{
\code{skill_cind} computes the consistency indices for omitting single or sets of skills.
}
\usage{
skill_cind(eKS, tKS, model, allsubs = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eKS}{
Empirical knowledge structure (cf. \link[skills]{skills-package})
in "\link[sets]{gset} of sets" - represtentation with
memberships equal to observed frequencies.
}
  \item{tKS}{
Theoretical knowledge structure (cf. \link[skills]{skills-package})
in "\link[sets]{set} of sets" - represtentation
}
    \item{model}{
Model for Skill Assignment (either "disjunctive" or "conjunctive")   
}
    \item{allsubs}{
if TRUE the weighted consistency index is computed for all knowledge structures
which can be obtained by leaving any subset of skills out of the skill assignment;
if FALSE only subsets containing one skill are used.
}
}
\details{
Given a skill assignment, specified through \code{tKS} and \code{model}, \code{skill_cind} omits one skill at a time in the skill
assignment defining the theoretical knowledge structure, and computes the weighted
consistency index (cf. \link[skills]{w_cind}) for the modified theoretical structure and the given \code{eKS}.\cr
All values will be lower than without leaving skills because the modified structure is a substructure of the original one.
If the resulting value barely differs from the original one this is an indicator for a low influence of this skill on the 
theoretical knowledge structure.
}
\value{
vector of consistency indices, names according to left skills
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
skill_cind(eKS_d, tKS_d, "disjunctive", allsubs = FALSE)

eKS = as.KS(as.gset(pisa))
tKS = as.KS(as.relation(eKS, empirical = TRUE, v = 1))
skill_cind(eKS, tKS, "disjunctive", allsubs = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{consistency}
\keyword{empirical}
\keyword{knowledge structure}
