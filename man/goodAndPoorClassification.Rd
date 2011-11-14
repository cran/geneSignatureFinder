\name{goodAndPoorClassification}
\alias{goodAndPoorClassification}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to classify the unsupervised clusters found through a
classification procedure in good and poor groups of samples.
}
\description{This function label the groups found by the unsupervised classification in good and poor.
}
\usage{
goodAndPoorClassification(clustering, percentile = 0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clustering}{list of so many elements as the number of samples containing the classification in two groups. Missing values are allowed.
%%     ~~Describe \code{clustering} here~~
}
  \item{percentile}{ real value between 0.0+ and 0.50- (see details)
}
}
\details{The labelling of the two groups in good and poor is obtained through a
robust interpolation of the two survival curves. 
For each of the two groups the percentiles at level percentile and 1-percentile of the survival times are computed. The straight line through these two percentile is then evaluated.
The group for which the line has the lowest slope is labelled as poor.  
%%  ~~ If necessary, more details than the description above ~~
}
\value{A factor with levels "good" and "poor".
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{Stefano M. Pagnotta and Michele Ceccarelli
%%  ~~who you are~~
}
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
# set the dimension of the cpu's cluster (use a value different from 3 depending on the number of cpu available)
aNCPUS <-NCPUS(nchips = 2)
# set the starting gene to STX1A
geneSeed <- which(colnames(geData) == "STX1A")
# run ...
aSignatureFinder <- signatureFinder(geneSeed, logFilePrefix = "test", cpuCluster = aNCPUS)
cclassification <- aSignatureFinder$classification

par(mfrow = c(1,2))
sf <- survfit(stData ~ cclassification)
plot(sf, main = paste("Survivals for signature starting from:", aSignatureFinder$startingSignature),
  xlab = paste("tValue(Log-Rank test) = ", round(aSignatureFinder$tValue, 3)),
  col = c("green", "red"))
  
cclassification <- goodAndPoorClassification(cclassification, percentile = .1)
sf <- survfit(stData ~ cclassification)
plot(sf, main = paste("Survivals for signature starting from:", aSignatureFinder$startingSignature),
  xlab = paste("tValue(Log-Rank test) = ", round(aSignatureFinder$tValue, 3)),
  col = c("green", "red"))
  
# if necessary to check all the signatures ...
#K <- length(searchResults)
#for(k in 1:K) {
#  aSignatureFinder <- searchResults[[k]]
#  cclassification <- aSignatureFinder$classification
#  cclassification <- goodAndPoorClassification(cclassification, pct = .1)
#  searchResults[[k]]$classification <- cclassification
#}    
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
