\name{as.relation.KS}
\Rdversion{1.1}
\alias{as.relation.KS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes the surmise relation of knowledge structures
}
\description{
Computation of the surmise relation for any knowledge structure.
}
\usage{
as.relation.KS(x, empirical = FALSE, v = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Knowledge structure as set, matrix or gset (for empirical structures) (cf. \link[skills]{skills-package}) of class \code{KS}
}
  \item{empirical}{
logical
}
  \item{v}{
parameter (numeric: 1,2,3) specifying the \link[DAKS]{iita} algorithm for computing the quasi order (only necessary for empirical knowledge strcutres)
}
  \item{\dots}{
further arguments
}
}
\details{
If \code{empirical = FALSE} the surmise relation of the given knowledge structure \code{x} is computed directly: The prerequisites for item q
are given as the intersection of those states which contain this item.\cr
If the parameter \code{empirical} is set to \code{TRUE} and \code{x} is given as \link[sets]{gset} the function \code{iita}
of package \pkg{DAKS} is applied which computes the best fitting quasi order behind empirical knowledge structures.
This relation is returned as \link[relations]{endorelation} and can be displayed by \link[relations]{relation_incidence}.
}
\value{
\link[relations]{endorelation} (in the displaying incidence matrix: row items as prerequisites of column items)
}
\references{
Sargin, A. and Uenlue, A. (in press) \emph{The R package DAKS: Basic functions and complex algorithms in knowledge space theory.}
In H. Locarek-Junge, C. Weihs, editors, \emph{Studies in Classification, Data Analysis, and Knowledge Organization},
Berlin-Heidelberg-New York: Springer.
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
\code{\link[skills]{as.KS}} for the converse direction, i.e. converting from relation to 'set'-/'matrix'-representation,\cr
\code{\link[relations]{relation}} for general information to relations,\cr
\code{\link[skills]{plothasse}} for plotting a Hasse diagram of a given partial order,\cr
\code{\link[DAKS]{iita}} is called to compute the relation for empirical knowledge structures.
}
\examples{
KS=set(set(),set(1),set(2),set(3),set(1,2),set(1,3),set(2,3),set(1,2,3))
KSr=as.relation(as.KS(KS))
relation_incidence(KSr)

eKS_c = gset(set(set(), set(3), set(3,5), set(2,4,5), set(1,3,5),
    set(1,2,3,4,5)), memberships = c(1,2,3,4,5,6))
eKS_r = as.relation(as.KS(eKS_c), empirical = TRUE, v = 1)
relation_incidence(eKS_r)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{surmise}
\keyword{knowledge structure}
\keyword{empirical}
