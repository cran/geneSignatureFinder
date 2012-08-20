\name{signatureSummaryTable}
\alias{signatureSummaryTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Build a data frame populated with the basic information of a signature.
%%  ~~function to do ... ~~
}
\description{A table with as many rows as the number of genes belonging to the
signature. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
signatureSummaryTable(aSignatureFinder)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aSignatureFinder}{Structure as results from the 
signatureFinder() function}
}

%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{A data frame with the following columns
  \item{column no.1}{genes in the signature}
  \item{column no.2}{percentage of missing data per gene}
  \item{column no.3}{importance (if computed)}
  \item{column no.4}{median level in good prognosis group (if computed)}
  \item{column no.5}{median level in poor prognosis group (if computed)}
  \item{column no.6}{mean level in good prognosis group (if computed)}
  \item{column no.7}{mean level in poor prognosis group (if computed)}
  \item{column no.8}{absolute difference of the mean levels (if computed)}
  \item{column no.9}{Welch's test-Value (if computed)}
  \item{column no.10}{p-value associated to t-value (if computed; see
testGE() function)}
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{Stefano M. Pagnotta and Michele Ceccarelli}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{signatureFinder}},
\code{\link{importance}}, \code{\link{testGE}}
}
\examples{
data(geNSCLC)
geData <- geNSCLC
data(stNSCLC)
stData <- stNSCLC
# set the dimension of the cpu's cluster (use a value different from 2 
# depending on the number of cpu avalilable)
aMakeCluster <- makeCluster(2)
signature <- signatureFinder(which(colnames(geData) == "HIF1a"),
  cpuCluster = aMakeCluster)
signature <- importance(signature)
signature <- testGE(signature)
signatureTable <- signatureSummaryTable(signature)
print(signatureTable)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
