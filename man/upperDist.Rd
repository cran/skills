\name{upperDist}
\Rdversion{1.1}
\alias{upperDist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hamming distances between upper bounds and empirical states
}
\description{
\code{upperDist} computes the hamming distances (\link[skills]{hammingDist})
between the upper bound (\link[skills]{lb_eKS}) and the empirical state itself for each
knowlegde state in \code{eKS}.
}
\usage{
upperDist(eKS, tKS, model)
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
%%  ~~ If necessary, more details than the description above ~~
}
\value{
list of distances each of which between an empirical state and its upper bound (cf. \link[skills]{ub_eKS}). The position corresponds
to that of the empirical state in the knowledge structure \code{eKS}.
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
\code{\link[skills]{SkillAss}} computes the necessary skill assignment,\cr
\code{\link[skills]{lb_eKS}} coontains the correspondend theoretical states,\cr
and \code{\link[skills]{hammingDist}} computes the distances.
}
\examples{
tKS_d = set(set(), set(2), set(1,3,4), set(1,2,3,4), set(1,2,4,5),
    set(1,2,3,4,5))
eKS_d = gset(set(set(), set(2,4), set(1,3), set(1,2,3,5), set(1,2,3,4,5)),
    memberships = c(1,2,3,4,5))
upperDist(eKS_d, tKS_d, "disjunctive")

tKS_c = set(set(), set(3), set(5), set(2,5), set(1,3,4,5), set(1,2,3,4,5))
eKS_c = gset(set(set(), set(3), set(3,5), set(2,4,5), set(1,3,5),
    set(1,2,3,4,5)), memberships = c(1,2,3,4,5,6))
upperDist(eKS_c, tKS_c, "conjunctive")

eKS = as.KS(as.gset(pisa))
tKS = as.KS(as.relation(eKS, empirical = TRUE, v = 1))
upperDist(eKS, tKS, "disjunctive")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{empirical}
\keyword{knowledge structure}
