\name{SkillAss}
\Rdversion{1.1}
\alias{SkillAss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Skill Assignment to items
}
\description{
The function \code{SkillAss} computes a skill assignment matrix to a given
knowledge structure. 
}
\usage{
SkillAss(KS, method = "DG", minimize = F, model = "disjunctive")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{KS}{
Theoretical knowledge structure in 'set' or 'matrix' representation (cf. \link[skills]{skills-package})
}
  \item{method}{
Two algorithms possible, either according to Doignon and Falmagne ("DF") or to Duentsch and Gediga ("DG")
}
  \item{minimize}{
logical: Should the resulting skill assignment be minimized?
}
  \item{model}{
either disjunctive, conjunctive or skill multimap
}
}
\details{
The function \code{SkillAss} derives a skill assignment matrix to a given
knowledge structure. Generally two different methods are possible: According
to Doignon and Falmagne ("DF") or to Duentsch and Gediga ("DG"). Both offer different
algorithms for the three models "disjunctive", "conjunctive" and "skill multimap".
If requested \code{minimize} attempts to minimize the skill set and thus the skill
assignment. By default \code{SkillAss} attempts to compute a disjunctive skill assignment
(\code{model = "disjunctive"}) according to Duentsch and Gediga ("DG") which is already minimal itself.
}
\value{
A skill assignment matrix with two attributes.\cr
The matrix has row names according to the items (and column names according to the skills; therefore the first letters of the ABC
are used). The attributes "method" and "model" specify which method and which model
the function have been used for the computed skill assignment.
}
\references{
Duentsch, I., Gediga, G. (1995), \emph{Skills and knowledge structures}. British Journal of Mathematical and Statistical
Psychology 48 (1995), 9 - 27,\cr
Duentsch, I., Gediga, G. (2002), \emph{Skill Set Analysis in Knowledge Structures}. British Journal of Mathematical and
Statistical Psychology, 55(2), 361 - 384.\cr
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
\code{\link[skills]{minSA}} for minimizing the returned skill assignment.
}
\examples{
# s1 is KnowledgeSpace
s1 = set(set(), set(1,3), set(2,3), set(1,2,3), set(1,2,4), set(1,2,3,4))
SkillAss(s1, method = "DG", model = "disjunctive")
SkillAss(s1, method = "DG", model = "skill multimap")
SkillAss(s1, method = "DF", model = "disjunctive")
SkillAss(s1, method = "DF", minimize = TRUE, model = "disjunctive")

\dontrun{
SkillAss(s1, method = "DG", model = "conjunctive")   
# error: conjuctive model not suitable for knowledge spaces}

# s2 is ClosureSpace
s2 = set(set(), set(3), set(5), set(2,5), set(1,3,4,5), set(1,2,3,4,5))
SkillAss(s2, model = "conjunctive")
SkillAss(s2, method = "DG", model = "skill multimap")

# s3 is neither nor
s3 = set(set(), set(2), set(3), set(3,4), set(1,2,3), set(1,3,4),
    set(2,3,4), set(1,2,3,4))
SkillAss(s3, model = "skill multimap")

# s4 is both
s4 = set(set(), set(1), set(2), set(3), set(1,2), set(1,3), set(2,3),
    set(1,2,3))
SkillAss(s4, model = "disjunctive")
SkillAss(s4, model = "conjunctive")


#s5 equals s2 only with letters
s5 = set(set(), set("c"), set("f"), set("b","f"), set("a","c","d","f"),
    set("a","b","c","d","f"))
SkillAss(s5, model = "conjunctive")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{skill assignment}
\keyword{knowledge space theory}
\keyword{skills}
\keyword{knowledge}
\keyword{kst}
