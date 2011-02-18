\name{lb_eKS}
\Rdversion{1.1}
\alias{lb_eKS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Lower bound of an empirical knowledge state in a theoretical structure
}
\description{
\code{lb_eKS} assignes a theoretical state as lower bound to each empirical state.
}
\usage{
lb_eKS(eKS, tKS, model)
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
In almost every study there will appear knowledge states which are not
contained in the knowledge structure modelled by theory due to lucky guesses
and careless errors. In case \code{tKS} is a closure space (closed under intersection,
formed by a conjunctive model) or a knowledge space (closed under union, formed
by a disjunctive model) Duentsch and Gediga suggest a way to find out a theoretical state
which is assumed to be real state of knowledge to each observed state.
}
\value{
A list, the same length as \code{eKS}, with states in \link[sets]{set}-representation. The j-th entry of the list is the upper bound
of the j-th state in \code{eKS}.
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
\code{\link[skills]{SkillAss}} computes the necessary skill assignment.
}
\examples{
# theoretical knowledge structure closed under union
tKS_d = set(set(), set(2), set(1,3,4), set(1,2,3,4), set(1,2,4,5),
    set(1,2,3,4,5))
eKS_d = gset(set(set(), set(2,4), set(1,3), set(1,2,3,5), set(1,2,3,4,5)),
    memberships = c(1,2,3,4,5))
as.list(eKS_d)
lb_eKS(eKS_d, tKS_d, "disjunctive")

# theoretical knowledge structure closed under intersection
tKS_c = set(set(), set(3), set(5), set(2,5), set(1,3,4,5), set(1,2,3,4,5))
eKS_c = gset(set(set(), set(3), set(3,5), set(2,4,5), set(1,3,5),
    set(1,2,3,4,5)), memberships = c(1,2,3,4,5,6))
as.list(eKS_c)
lb_eKS(eKS_c, tKS_c, "conjunctive")

# real dataset; the quasi order (closed under union and intersection)
# computed by iita (cf.as.relation.KS)
eKS = as.KS(as.gset(pisa))
tKS = as.KS(as.relation(eKS, empirical = TRUE, v = 1))
tKS
as.list(eKS)
lb_eKS(eKS, tKS, "disjunctive")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{bound}
\keyword{empirical}
\keyword{knowledge structure}
