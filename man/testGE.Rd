\name{testGE}
\alias{testGE}
\title{Test the differential expression of the the genes in a signature with respect to the good and poor prognosis groups.}
\description{Given the clustering of the samples in good and poor prognosis associated to the signature, for each gene in the signature the test for the null hypothesis of equality of the expression levels is performed. Additional statistics are provided.
}
\usage{
testGE(aSignatureFinder, permutationReplications = 1000, cpuCluster = NULL, stopCpuCluster = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aSignatureFinder}{(structure) as results from the function signatureFinder().}
  \item{permutationReplications}{(integer) number of replications of the permutation test (default: 1000).}
  \item{cpuCluster}{structure as result from the NCPUS() function}
  \item{stopCpuCluster}{flag to control if the channel to the cpu-cluster has to be closed.}
}

\details{The t-test for testing the differential expression of the genes in the
signature is performed according to the procedure of Dudoit et al. (2002). The
test statistics is the Welch's one and the null distribution is obtained through
a permutation scheme. }


\value{
The function returns the same variable in the input aSignatureFinder structure and 
	\item{groupMedian}{real matrix with as many rows as length(aSignatureFinder$signature) and two columns containing the medians of each gene with respect to the good and poor prognosis group}
	\item{medianAbsDifference}{a list of real with as many elements as length(aSignatureFinder$signature) where each entry is the absolute difference of the medians computed in each group}
	\item{groupMean}{real matrix with as many rows as length(signature$signature) and two columns containing the means of each gene with respect to the good and poor prognosis group}
	\item{meanAbsDifference}{a list of real with as many elements as length(aSignatureFinder$signature) where each entry is the absolute difference of the means computed in each group}
	\item{meanDifferenceTValue}{a list of real with as many elements as length(aSignatureFinder$signature) where each entry is the value of the test statistic}
	\item{meanDifferencePValue}{a list of real with as many elements as length(aSignatureFinder$signature) where each entry is the p-value of the test statistic}
}
\references{Dudoit S.  {   et   al.}:   Statistical methods for identifying
differentially expressed genes in replicated cDNA microarray experiments
{Statistica Sinica}, 12, pp. 111--139, 2002.}
\author{Stefano M. Pagnotta and Michele Ceccarelli}
%%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{signatureFinder}},
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# find the signature starting from the gene SELP for  the Non Small Cell Lung Cancer 
# set the working data 
data(geNSCLC)
geData <- geNSCLC
data(stNSCLC)
stData <- stNSCLC
# set the dimension of the cpu's cluster 
# (use a value different from 3 depending on the number of cpu avalilable)
aNCPUS <-NCPUS(nchips = 2)
# set the starting gene to SELP
geneSeed <- which(colnames(geData) == "SELP")
# run ...
ans <- signatureFinder(geneSeed, logFilePrefix = "test", cpuCluster = aNCPUS)
print(ans)
ans <- testGE(ans)
print(ans)
barplot2(t(ans$groupMean), beside = TRUE,
                   main = paste("Signature starting from:", ans$startingSignature),
                   legend = paste(colnames(ans$groupMedian), "prognosis group"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
