\name{summary.KS}
\alias{summary.KS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary Method for Objects of Class KS
}
\description{
S3 method to summarize objects of class KS.
}
\usage{
summary.KS(object, npattern = 5, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Knowledge structure (cf. \link[skills]{skills-package}) of class KS
}
  \item{npattern}{
for empirical knowledge structures (represented as \link[sets]{gset}): how many of the most frequent states should be reported
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
list
\item{domain}{cf.\link[skills]{KSdomain}}
\item{Freq}{For a theoretical knowledge structure: table of frequency for every item; for a empirical knowledge structure: table of observed frequencies}
\item{N}{only available for empirical knowledge structures: overall number of observations}
\item{pattern}{only availabe for empirical knowledge structures: most frequent response patterns as data frame}
\item{notions}{cf. \link[skills]{KSnotions}}
\item{Nstates}{Number of states}
\item{atoms}{cf. \link[skills]{KSatoms}}
\item{KSpace}{cf. \link[skills]{is.KSpace}}
\item{base}{only if KSpace == TRUE: \link[skills]{KSbase}; else not available}
}
\references{
Doignon, J.-P., Falmagne, J.-C. (1999), \emph{Knowledge Spaces}. Berlin, Heidelberg and New York: Springer-Verlag.
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
s1 = set(set(), set(1,3), set(2,3), set(1,2,3), set(1,2,4), set(1,2,3,4))
summary(as.KS(s1))

s2 = as.KS(as.gset(pisa))
summary(s2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{knowledge structure}
