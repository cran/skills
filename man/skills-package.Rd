\name{skills-package}
\Rdversion{1.1}
\alias{skills-package}
\docType{package}
\title{
Skills in knowledge space theory
}
\description{
Knowledge space theory analyses the knowledge of persons in so called knowledge states meaning the
information which questions of a test (called items) a person was able to master.\cr
Whilst methods for the analysis of knowledge structures themselves are already available in different packages such as
\pkg{DAKS} by A.Sargin and A.Uenlue (2009) the interest of research often lies in the skills behind these items.
For this purpose the package \pkg{skills} provides functions to detect skills behind a knowledge structure
(cf. \link[skills]{SkillAss}) or to compute a knowledge structure from a skill assignment
(cf. \link[skills]{probfun}). In addition the package contains functions for general knowledge space theory as well as
functions which allow a comparison between a theoretical (constructed) knowledge structure and an empirical (observed) one.\cr
The functions in this package are in accordance with papers from Duentsch and Gediga and Doignon and Falmagne.
}
\details{
\tabular{ll}{
Package: \tab skills\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-02-18\cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab yes\cr
}
A theoretical knowledge structure can either be represented as a \link[sets]{set} (of sets), as an item-state-matrix
or as an \link[relations]{endorelation}. The second one is a 0-1-matrix whose column vectors represent the states where 1
indicates that the correponding row item is contained. The item labels are available as the rownames-attribute of the matrix.
The endorelation refers to the surmise relations of knowledge structures.\cr
An empirical knowledge structure is represented as a \link[sets]{gset} of sets with memberships equal to
the observed frequencies.
}
\author{
Angela Pilhoefer, University of Augsburg, Germany

Maintainer: \cr
Angela Pilhoefer \email{angela.pilhoefer@gmx.de}

}
\references{
Duentsch, I., Gediga, G. (2002), \emph{Skill Set Analysis in Knowledge Structures}. British Journal of Mathematical and
Statistical Psychology, 55(2), 361 - 384.,\cr
Duentsch, I., Gediga, G. (1995), \emph{Skills and knowledge structures}. British Journal of Mathematical and Statistical
Psychology 48 (1995), 9 - 27,\cr
Doignon, J.-P., Falmagne, J.-C. (1999), \emph{Knowledge Spaces}. Berlin, Heidelberg and New York: Springer-Verlag,\cr
Doignon, J.-P., Falmagne, J.-C. (1994), \emph{Knowledge spaces and Skill assignments.} In G.H.Fischer, D.Laming, editors,
\emph{Recent research in psychology: Contributions to mathematical psychology, psychometrics, and methodology},
chapter 8, pages 111--121, Berlin, Heidelberg and New York: Springer-Verlag, 1994,\cr
Sargin, A., Uenlue, A. (in press), \emph{The R package DAKS: Basic functions and complex algorithms in knowledge space theory.}
In H. Locarek-Junge, C. Weihs, editors, \emph{Studies in Classification, Data Analysis, and Knowledge Organization},
Berlin-Heidelberg-New York: Springer,\cr
Haidinger, A., Uenlue, A. (in preparation), \emph{Skills and Competencies Models in Knowledge Space Theory: An Implementation in R.}.
}

\keyword{ package }
\seealso{

}
\examples{
# convert: Conversion between set and matrix representation of knowledge
# structures
KS1 = set(set(), set("c"), set("f"), set("b","f"), set("a","c","d","f"),
    set("a","b","c","d","f"))
KS1
convert(KS1)
KS2 = gset(set(set(), set(3), set(3,5), set(2,4,5), set(1,3,5),
    set(1,2,3,4,5)), memberships = c(1,2,3,4,5,6))
KS2
convert(KS2, return.dataframe = TRUE)

KS3=matrix(c(0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,1,0,0,1,1,0,1,1,1,1,1,1,1,1),
    nrow = 5)
KS3
convert(KS3)

# probfun: Computation of the problem function
SA = matrix(c(1,1,0,1,1,1,0,1,1,0,0,1,0,0,0), nrow = 5)
rownames(SA) = 1:5
SA
probfun(SA, "conjunctive")
convert(probfun(SA, "conjunctive"))

# SkillAss: Skill Assignment to items
# s1 is KnowledgeSpace
s1 = set(set(), set(1,3), set(2,3), set(1,2,3), set(1,2,4), set(1,2,3,4))
s1
SkillAss(s1, method = "DG", model = "disjunctive")
}
