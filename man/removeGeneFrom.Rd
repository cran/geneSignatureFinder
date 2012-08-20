\name{removeGeneFrom}
\alias{removeGeneFrom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Remove the gene with the lowest importance from a signature}
\description{This function implements a pruning algorithm devoted to remove the
gene with the lowest importance from the signature. The importances have
to be computed before calling this function. The gene with the
lowest importance is removed and the set of importances are computed again. The
function is designed to be used iteratively so that all the genes with
importance below a cutoff value are removed. If the signature has the results
of the testGEL() function just computed, they are removed.
}
\usage{removeGeneFrom(aSignatureFinder, cutoff = 0)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aSignatureFinder}{Structure as results from the function
signatureFinder()}
  \item{cutoff}{is a real value used when the function works iteratively}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{The function return NULL when no gene is removed, otherwise return a
structure as results from the importance() function with the additional slots 
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{removedGene}{list of genes (as string) removed from the signature}
  \item{originalSignature}{same as aSignatureFinder$signature}
  \item{originalClassification}{same as aSignatureFinder$classification}
  \item{originalTValue}{same as aSignatureFinder$tValue}
  \item{originalPValue}{same as aSignatureFinder$pValue}
}

%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{Stefano M. Pagnotta and Michele Ceccarelli}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{signatureFinder}}, \code{\link{importance}}.
}
\examples{
data(geNSCLC)
geData <- geNSCLC
data(stNSCLC)
stData <- stNSCLC
# develop the signature
aMakeCluster <- makeCluster(2)
geneSeed <- which(colnames(geData) == "RNF5")
signature <- signatureFinder(geneSeed, cpuCluster = aMakeCluster,
                             stopCpuCluster = FALSE)
signature <- importance(signature, cpuCluster = aMakeCluster)
print(signature$importance)

################################
# 1 step pruning
removeGeneFrom(signature, cutoff = 0.3)$signature

##################################
# iterative pruning of the signature
signature
repeat {
  tmp <- removeGeneFrom(signature, cutoff = 0.3)
  if(is.null(tmp)) break
  signature <- tmp
}  
signature
signature$importance
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
