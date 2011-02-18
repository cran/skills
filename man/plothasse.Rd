\name{plothasse}
\alias{plothasse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hasse diagram for knowledge structures
}
\description{
Plots a hasse diagram for a knowledge structure based on the surmise relation.
}
\usage{
plothasse(KS, empirical = FALSE, v = 1, concat = ":", lab.cex = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{KS}{
Theoretical ('set' or 'matrix' representation) or empirical ('gset') knowledge structure (\link[skills]{skills-package})
}
  \item{empirical}{
logical, if TRUE first a quasi order is computed by \link[DAKS]{iita} (package DAKS) which takes careless errors and lucky
guesses into account, if FALSE the quasi order is computed by \emph{r<q iff r in Intersection(K_q)}
}
  \item{v}{
parameter (numeric: 1,2,3) specifying the \link[DAKS]{iita} algorithm for computing the quasi order
}
  \item{concat}{
symbol for concatenation of equally informative items (cf. \link[skills]{discred})
}
  \item{lab.cex}{
mulipier for font size
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
\code{\link[skills]{as.relation.KS}} computes the relation underlying the given knowledge structure,\cr
\code{\link[skills]{discred}} computes the discriminative reduction such that the quasi order is even a partial
order which is necessary for plotting a Hasse diagram.
}
\examples{
KS = set(set(), set(1), set(2), set(3), set(1,2), set(1,3), set(2,3),
    set(1,2,3), set(1,2,4), set(1,2,3,4))
plothasse(KS)

z = as.gset(pisa)
plothasse(as.KS(z), empirical = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hasse}
\keyword{surmise}
\keyword{knowledge structure}
