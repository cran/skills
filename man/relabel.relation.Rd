\name{relabel.relation}
\Rdversion{1.1}
\alias{relabel.relation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Relabelling the items in a relation
}
\description{
This function relabels the items in a relation
}
\usage{
relabel.relation(er, domain = NULL, data = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{er}{
a relation
}
  \item{domain}{
vector with new names of the items
}
  \item{data}{
dataframe whose names should be the new names of the items
}
}
\details{
Many functions using the package \pkg{relations} return a relation with item labelled by
"1L", "2L", ... . The function \code{relabel.relation} allows to relabel the items by given names.
}
\value{
an endorelation with changed names 
}
\references{
%% ~put references to the literature/web site here ~
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
gr = iita(pisa, v = 1) $ implications
re = endorelation(graph = gr)
relation_incidence(re)
re.new = relabel.relation(re, data = pisa)
relation_incidence(re.new)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{relabel}
\keyword{relation}
