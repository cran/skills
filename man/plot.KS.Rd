\name{plot.KS}
\alias{plot.KS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plotKS
}
\description{
Plots a hasse diagram for the knowledge states.
}
\usage{
plot.KS(x, y = NULL, ord = NULL, max.scale = "rel", lab.cex=1.5,
    border = "black", bg = c("grey70", "skyblue", "tomato"),
    fill = c("tomato", "dodgerblue"), shape = "rect", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Theoretical or empirical knowledge structure (cf. \link[skills]{skills-package})
}
  \item{y}{
Theoretical or empirical knowledge structure (cf. \link[skills]{skills-package})
}
  \item{ord}{
order for the items
}
  \item{max.scale}{
either "abs", "rel" or a numeric value in (0,1]; "abs" plots the absolute frequencies whereas both others plot
relative frequencies of the states, for "rel" the maximum of the scale is the maximal relative frequency otherwise it is equal to the
specified numeric
}
  \item{lab.cex}{
mulipier for font size
}
  \item{border}{
colour of the border
}
  \item{bg}{
colour of the background
}
  \item{fill}{
colour specifying the relative frequencies of the states
}
  \item{shape}{
shape of the states
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
For a single theoretical knowledge structure \code{plot} graphs a hasse-diagram of all states (case I). Given two theoretical
knowledge structures \code{x} and \code{y} \code{plot(x,y)} returns a hasse-diagram containing all states of both structures.
They are distinguishable in terms of colors and the states of each structure are connected by lines only to those of the same
structure (case II).\cr
A single empirical structure is graphed similar to case I with only one difference that the boxes symbolizing the states are
colored according to the relative frequency of observations of this state measured by the number of the most frequently state.
In case of one empirical and one theoretical knowledge structure the graph is similar to case II again with the difference that
the boxes of the empirical states are only colored according to the relative frequency as above. Note that here the first
parameter \code{x} has to be the empirical structure and the second \code{y} the theoretical one. If both knowledge structures
are empirical the underlying graph again corresponds to case II but the coloring of the states depends on the difference of the
observed frequencies, i.e. if one state is contained in only one of both structures it is colored like before, but elsewise that
part of the box is colored which lies between the different numbers of observations for this state and namely in the color of the
predominant structure.\cr
If at least one empirical knowledge structure is displayed the parameter \code{max.scale} can be specified either by \code{abs},
\code{rel} or a numeric value in (0,1]. \code{abs} plots the absolute frequencies whereas both others plot relative frequencies
of the states, for "rel" the maximum of the scale is the maximal relative frequency otherwise it is equal to the specified numeric.
This is in particular a good option if two empirical structures are displayed. If the sample size, i.e. the number of testees who
took place in the test, differ from each other, the distribution on the knowledge states can be compared thereby.
}
\value{

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

KS = set(set(), set("a"), set("b"), set("a","b"), set("a","c"),
     set("a","b","c"), set("a","b","c","d"), set("a","b","c","e"),
     set("a","b","c","d","e"), set("a","b","c","d","e","f"))
plot(as.KS(KS))

gpisa = as.KS(as.gset(pisa))
plot(gpisa)

pisaTh1 = as.KS(as.relation(gpisa, empirical = TRUE, v = 1))
plot(pisaTh1)
pisaTh3 = as.KS(as.relation(gpisa, empirical = TRUE, v = 3))
plot(pisaTh3)
plot(pisaTh1, pisaTh3)
# gpisaTh results from pisa by omitting lucky guesses and thus assigning
# each empirical state that one which characterizes
# its upper bound (cf. ub_eKS).
# The according frequencies are summed up.
gpisaTh = as.KS(gset(set(set(), set("a"), set("a","b"), set("a","b","c"),
    set("a","b","c","d"), set("a","b","c","e"), set("a","b","c","d","e")),
    memberships = c(20,41,72,94,66,34,13)))
plot(gpisa, gpisaTh)

# the max.scale = "abs" option
# it is also possibe to use absolute frequencies for the scales
# gpisaThx is 3 times gpisaTh
gpisaThx = as.KS(gset(set(set(), set("a"), set("a","b"), set("a","b","c"),
    set("a","b","c","d"), set("a","b","c","e"), set("a","b","c","d","e")),
    memberships = c(20,41,72,94,66,34,13) * 3))
plot(gpisa, gpisaThx)
plot(gpisa, gpisaThx, max.scale = "rel")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{knowledge structure}
\keyword{hasse}
