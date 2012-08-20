\name{importance}
\alias{importance}

\title{Compute an index of importance of the genes in the signature.
}

\description{An index measuring the importance of each gene in the signature is computed. The more the index is close to 1, the most the corresponding gene is important to separate the survival curves. When the value is about 0, the gene has a low contribution in separating the survival curves. If it happens that the index is negative, the deletion of the gene from the signature improves the distance of the survival curves when all but it are in the signature.
}

\usage{
importance(aSignatureFinder, deep = FALSE, cpuCluster = NULL, stopCpuCluster = TRUE)
}


\arguments{
  \item{aSignatureFinder}{Structure as results from the function signatureFinder()}
  \item{deep}{flag for future use}
  \item{cpuCluster}{structure as result from the NCPUS() function}
  \item{stopCpuCluster}{flag to control if the channel to the cpu-cluster has to be closed.}
}

\details{The importance index associated to each gene of a signature is computed according to the leave one out strategy. Given a gene, it is removed from the signature and with the other genes in the  sequence first a classification is performed and then thet test-value of the log-rank test is evaluated. These test-values are collected in the list L1GeneOutTV. All these values are compared to test-value of the signature: importance = 1 - L1GeneOutTV/aSignatureFinder$tValue.
}
\value{The function returns the same variable in the input signature structure and 
  \item{L1GeneOutTV}{a list of real values with the test-values of the log-rank
test applied to all but one gene in the signature (leave-one-out strategy)}
  \item{importance}{A list of real values measuring the importance of each gene in the signature}
}

%%\references{
%%}

\author{Stefano M. Pagnotta and Michele Ceccarelli}

%%\note{
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{signatureFinder}}.
}

\examples{
data(geNSCLC)
geData <- geNSCLC
data(stNSCLC)
stData <- stNSCLC
aMakeCluster <- makeCluster(2)
geneSeed <- which(colnames(geData) == "SELP")
signature <- signatureFinder(geneSeed, 
  cpuCluster = aMakeCluster, stopCpuCluster = FALSE)
signature <- importance(signature, cpuCluster = aMakeCluster)
signature
signature$importance
#################################
barplot(signature$importance, main = "Importance based on L1GeneOut",
        sub = paste("Signature starting from:", signature$startingSignature))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
