\name{minSA}
\Rdversion{1.1}
\alias{minSA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Minimize skill assignment
}
\description{
Minimization of the skill set and restriction of the skill-assignment-matrix
such that the knowledge structure resulting from the problem function (cf. \link[skills]{probfun}) remains the same.
}
\usage{
minSA(SA, model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SA}{
Skill assignment matrix which should be minimized
}
  \item{model}{
Either disjunctive, conjunctive or skill multimap (cf. \link[skills]{probfun})
}
}
\details{
In order to find the most economically solution of a skill assignment,
that means a skill assignment with as few as possible skills, \code{minSA}
checks if leaving skills out results in the same knowledge structure
using the given model. If no skill is redundant \code{minSA} returns
the initial skill assignment.
}
\value{
A matrix, similar to \code{SA}, with just as many rows (items retain unchanged)
but possibly fewer columns (skill set gets minimized).
}
\references{
Doignon, J.-P., Falmagne, J.-C. (1999), \emph{Knowledge Spaces}. Berlin, Heidelberg and New York: Springer-Verlag,\cr
Doignon, J.-P., Falmagne, J.-C. (1994), \emph{Knowledge spaces and Skill assignments.} In G.H.Fischer, D.Laming, editors,
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
SA1 = matrix(c(0,1,0,0,0,1,0,1,1,0,1,1,1,1,0,1,1,0,1,1,1,1,1,1,1),
    nrow = 5)
rownames(SA1) = 1:5
colnames(SA1) = letters[1:5]
SA1
minSA1 = minSA(SA1, "disjunctive")
minSA1

minSA(minSA1, "disjunctive")    
# remains the same because minSA1 is already minimal

# s1 is a knowledge space
s1 = set(set(), set(1,3), set(2,3), set(1,2,3), set(1,2,4), set(1,2,3,4)) 

SAdg = SkillAss(s1, method="DG", model="disjunctive", minimize=FALSE)
SAdg
SAdgmin = minSA(SAdg, "disjunctive") 
# remains the same -> SAdg is already minimal

SAdf = SkillAss(s1, method="DF", model="disjunctive", minimize=FALSE)
SAdf
SAdfmin = minSA(SAdf, "disjunctive") 
SAdfmin
# minimized and now equal to SAdg resp. SAdgmin
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{knowledge structure}
\keyword{skill assignment}
