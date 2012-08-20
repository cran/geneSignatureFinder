\name{seedsFinder}
\alias{seedsFinder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Evaluate some statistics on all genes in order to select those that can
be used as seeds for searching the signatures.
%%  ~~function to do ... ~~
}
\description{This function works on each column (gene expression level) of the geData and returns the test-value and p-value of the log-rank test, the bayesian information criterion value under the hypothesis that tha data are drawn from a single gaussian (bic1) and a mixture of two gaussians (bic2); at the end the clustering of the samples is added. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
seedsFinder(cutoff = 1.95, evaluateBICs = TRUE, cpuCluster = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cutoff}{argument passed to the BICs() function.}
  \item{evaluateBICs}{flag to force the computation of the  bayesian information criteria.}
  \item{cpuCluster}{If a parallel search is necessary, this variable has to be set to the output of NCPUS() function.}
}
\details{For each gene expression levels data an unbiased classification is performed resulting into two clusters coded by the values 0 and 1. The samples classified by 0 are those for which the mean is lower than that of the samples classified with 1. The classification method is the partitioning around medoids algorithm linked to the a leave-one-out re-classification strategy (see the pamUmbiased() function for further details). From the clusters two survival curves are estimated with the stData data and then tested for the null hypothesis of no difference among them (see the survdiff() function for further details) providing the tValue. The correponding pValue is given by 1-pchisq(tValue, df = 1). Two more indexes are computed, the bayesian information criteria under the hypotheses 1) the gene levels are from a univarite gaussians (bic1) and 2) the gene levels are from a mixture of two gaussians (bic2) (see the BICs() function for further details). The mixing coefficient is estimated from the classification as the fraction of samples classified as 1. The parameters of the gassians are robustly estimated.
%%  ~~ If necessary, more details than the description above ~~
}
\value{The result of the function is a matrix having so many rows as ncol(geData) and 4+nrow(geData) rows.
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{column no.1: tValue}{test-value of the log-rank test statistic under the null hypothesis that the two survival curves are equal (see details)}
\item{column no.2: pValue}{p-value corresponding to the test-value in column no.1; actuallly is 1-pchisq(tValue, df = 1)}
\item{column no.3: bic1}{value of the bayesian information criterion computed under the hypothesis that the data are drawn from a single gaussian}
\item{column no.4: bic1}{value of the bayesian information criterion computed under the hypothesis that the data are drawn from mixture of two gaussians}
\item{columns from no.5 to no.4+nrow(geData)}{result of the unbiased classification (see details)}
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

\seealso{\code{\link{pam}}, \code{\link{survdiff}}, \code{\link{BICs}}.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(geNSCLC)
geData <- geNSCLC

data(stNSCLC)
stData <- stNSCLC

# here few genes and samples are considered to speed up the timing of the example.
# please, try 
# genesToUse <- which(apply(!is.na(geData), 2, sum)/nrow(geData) > 0.75)
# geData <- geData[, genesToUse]
# and comment stData <- stData[1:50, ]
genesToUse <- which(apply(!is.na(geData), 2, sum) == nrow(geData))
geData <- geData[, genesToUse]
geData <- geData[1:50, ]
stData <- stData[1:50, ]
dim(geData)

aMakeCluster <- makeCluster(2)
aSeedsFinder <- seedsFinder(cpuCluster = aMakeCluster)
head(aSeedsFinder)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
