\name{probfun}
\Rdversion{1.1}
\alias{probfun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computation of the problem function
}
\description{
Given an item-skill-assignment and a model (disjunctive, conjunctive, skill multimap)
the problem function computes the consequential knowledge structure.
}
\usage{
probfun(SA, model = "disjunctive", alldata = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SA}{
item-skill-matrix which row names name the items
}
  \item{model}{
either disjunctive, conjunctive or skill multimap
}
  \item{alldata}{
logical (see details)
}
}
\details{
An item-skill-assignment is a 0-1-matrix which ones indicate which skills are assigned to which items.
The row names of \code{mat} will be used as names of the items.
If method equals disjunctive (default) or conjunctive the rownames have to be unique otherwise \code{probfun}
changes to skill multimap automatically. If the parameter \code{alldata} is set to TRUE additionally two new attributes are returned:
a knowledge structure with nonunique states and the powerset of the skill set, comparing the columns
of the latter two specifying which subset of skills led to which state.
}
\value{
A knowledge structure in matrix representation (cf. \link[skills]{skills-package}) or a list as described in the details.
}
\references{
Duentsch, I., Gediga, G. (2002), \emph{Skill Set Analysis in Knowledge Structures}. British Journal of Mathematical and
Statistical Psychology, 55(2), 361 - 384.,\cr
Doignon, J.-P. and Falmagne, J.-C. (1999) \emph{Knowledge Spaces}. Berlin, Heidelberg and New York: Springer-Verlag,\cr
Doignon, J.-P. and Falmagne, J.-C. (1994) \emph{Knowledge spaces and Skill assignments.} In G.H.Fischer, D.Laming, editors,
\emph{Recent research in psychology: Contributions to mathematical psychology, psychometrics, and methodology},
chapter 8, pages 111--121, Berlin, Heidelberg and New York: Springer-Verlag, 1994.
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
\code{\link[skills]{SkillAss}} for the contrary direction.
}
\examples{
SA = matrix(c(1,1,0,1,1,1,0,1,1,0,0,1,0,0,0), nrow = 5)
rownames(SA) = 1:5
KS1 = probfun(SA, "conjunctive")
KS2 = probfun(SA, "disjunctive")
KS3 = probfun(SA, "disjunctive", alldata = TRUE)

SA = matrix(c(1,1,0,1,1,1,0,1,1,0,0,1,0,0,0), nrow = 5)
rownames(SA) = c(1,1,2,2,3)
KS4 = probfun(SA, "skill multimap")

\dontrun{
KS5 = probfun(SA, "disjunctive")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{skill assignment}
\keyword{knowledge structure}
