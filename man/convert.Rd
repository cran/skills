\name{convert}
\Rdversion{1.1}
\alias{convert}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conversion between set and matrix representation of knowledge structures
}
\description{
A knowledge structure can either be represented as a \link[sets]{set} (of sets) or as an item-state-matrix.
The \code{convert} function automatically recognizes the type of the given knowledge structure and changes it to the other.
}
\usage{
convert(KS, domain = NULL, return.dataframe = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{KS}{
Theoretical or empirical knowledge structure (cf. \link[skills]{skills-package})}
  \item{domain}{
if declared: a \link[sets]{set} containing at least all items of \code{KS}}
  \item{return.dataframe}{
logical. If TRUE an empirical \code{KS} is returned as a dataframe of observed states containing a frequency variable \code{Freq}. 
}
}
\details{
The \code{convert} function automatically recognizes the type of the given knowledge structure and changes it to the other.
The item-state-matrix representation is a 0-1-matrix which rows stand for the items, the columns for the states and a indicates
which items are contained in which states. The row names of the matrix name the items.
If \code{KS} is given as such a matrix which does not have row names then the first letters of the alphabet are chosen to be the
names of the items.\cr
Please note: After conversion the type of the rownames as well as of the items of the states is character.\cr
If \code{KS} is an empirical knowledge structure as passed as \link[sets]{gset} the convert function either returns a 
matrix with columns equal to the states (if parameter \code{return.dataframe = FALSE}) and omits the frequencies of the states or
a data frame with a frequency vector according to the memberships of the gset (if parameter \code{return.dataframe = TRUE}).
}
\value{
Set of sets, a matrix or data frame (see details).
If \code{KS} is neither a matrix nor a set of sets nor a gset no conversion will be applied and \code{KS} will be returned unchanged.
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
\code{\link[skills]{as.gset.data.frame}} for converting a data frame into a \link[sets]{gset}
}
\examples{
KS1 = set(set(), set(1,3), set(2,3), set(1,2,3), set(1,2,4), set(1,2,3,4))
convert(KS1)

KS2=set(set(), set("c"), set("f"), set("b","f"), set("a","c","d","f"),
    set("a","b","c","d","f"))
convert(KS2)

KS3=matrix(c(0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,1,0,0,1,1,0,1,1,1,1,1,1,1,1),
    nrow = 5)
KS3
convert(KS3)

KS4=KS3
rownames(KS4) = c("as", "df", "gh", "hj", "jk")
convert(KS4)

eKS = gset(set(set(), set(3), set(3,5), set(2,4,5), set(1,3,5),
    set(1,2,3,4,5)), memberships = c(1,2,3,4,5,6))
convert(eKS, return.dataframe = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{convert}
%\keyword{knowledge structure}
