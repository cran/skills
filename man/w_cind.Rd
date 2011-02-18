\name{w_cind}
\Rdversion{1.1}
\alias{w_cind}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Weigthed index for consistency between empirical and theoretical knowledge structures
}
\description{
The weigthed consistency index \code{w_cind} indicates the proportion of people whose
empirical knowledge state is contained in the theoretical knowledge structure.
}
\usage{
w_cind(eKS, tKS)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eKS}{
Empirical knowledge structure (cf. \link[skills]{skills-package})
in "\link[sets]{gset} of sets" - representation with
memberships equal to observed frequencies.
}
  \item{tKS}{
Theoretical knowledge structure (cf. \link[skills]{skills-package})
in "\link[sets]{set} of sets" - representation
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
proportion
}
\references{
Duentsch, I., Gediga, G. (2002), \emph{Skill Set Analysis in Knowledge Structures}. British Journal of Mathematical and
Statistical Psychology, 55(2), 361 - 384.
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
\code{\link[skills]{cind}} only counts the proportion of empirical states which are contained in the theoretical structure, whereas\cr
\code{\link[skills]{acm}} is an improvement which also takes the expected value of arbitrary assignments of items to skills sets
into account.
}
\examples{
tKS = set(set(), set(2), set(1,3,4), set(1,2,3,4), set(1,2,4,5),
    set(1,2,3,4,5))
eKS = gset(set(set(), set(2), set(1,3), set(1,2,4,5), set(1,2,3,4,5)),
    memberships = c(1,2,3,4,5))
w_cind(eKS, tKS)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{consistency}
