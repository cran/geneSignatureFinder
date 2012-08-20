\name{goodAndPoorClassification}
\alias{goodAndPoorClassification}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to classify the unsupervised clusters found through a
classification procedure in good and poor groups of samples.
}
\description{This function label the groups found by the unsupervised classification in good and poor.
}
\usage{
goodAndPoorClassification(clustering)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clustering}{factor of so many elements as the number of samples containing the classification in two groups. Missing values are allowed.
%%     ~~Describe \code{clustering} here~~
}
  
}
\details{The labelling of the two groups in good and poor is obtained through the comparion of the expected and the observed failures computed in the test statistics of the log-rank test. The group for which the expected are more than the observed are classified are good prognosis. The variable stData has to be defined in the environment.
}
\value{A factor with levels "good" and "poor".
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{Stefano M. Pagnotta, Michele Ceccarelli and Pietro Zoppoli}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~
%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
data(geNSCLC)
geData <- geNSCLC

data(stNSCLC)
stData <- stNSCLC
######################
randomClustering  <- runif(nrow(stData)) > .5
survdiff(stData ~ randomClustering)
###############
clusteringByPrognosis <- goodAndPoorClassification(randomClustering)
survdiff(stData ~ clusteringByPrognosis)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
