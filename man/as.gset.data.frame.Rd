\name{as.gset.data.frame}
\alias{as.gset.data.frame}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Converts a dataframe into a gset
}
\description{
The function \code{as.gset.data.frame} is based on the generic function \link[sets]{as.gset}.\cr
It takes a data frame of response patterns (=observed states) and converts the data into a \link[sets]{gset} with support equal to the states
and memberships equal to the frequency.

}
\usage{
as.gset.data.frame(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A raw dataframe or one with a frequency vector (named "Freq")
}
}
\details{
Similar to the matrix representation of a knowledge structure a data frame can be seen as a representation of
an empirical knowledge structure if its column names specify the items and the single rows represent the different response patterns where a 1 in the j-th
column indicates that the j-th item was answered correctly and a 0 that if it was not solved. The data frame can either be a raw data frame
(without frequency vector) and thus contains the response pattern for each person in a single row or a frequency table with unique response pattern
but summed up frequency in an additional column called "Freq".
}
\value{
\link[sets]{gset} of sets containing the different observed states with memberships equal to their frequency.
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
\code{\link[skills]{convert}} with \code{return.dataframe = TRUE} for the conversion in the contrary direction,\cr
\code{\link[sets]{gset}} for general information to the class 'gset'.
}
\examples{
z = as.gset(pisa)
z

mat=matrix(c(0,0,0,1,0,1,0,0,0,0,1,1,0,1,1,1,0,1,0,0,0,0,1,1,0,0,1,1,1,1),
    ncol = 5)
mat = cbind(mat, 1:6)
colnames(mat) = c(1:5, "Freq")
dat = as.data.frame(mat)
g = as.gset(dat)
g
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{gset}
\keyword{dataframe}
