\name{ensembleTable}
\alias{ensembleTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Build a data frame populated with statistical indexes for each gene of each signature.
}
\description{Its a summary table from which the ensemble signature is obtained.
}
\usage{
ensembleTable(aSearchResults)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aSearchResults}{ a vector of list
%%     ~~Describe \code{aSearchResults} here~~
}
}
\details{aSearchResults is a vector of list having as many items as the
number of signature to be developed. Each entry of the vector is a structure as
results from the signatureFinder() function. An example of such a vector can be
found in the source of the puttingAllTogether() function.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 \item{column no.1}{probe/gene found at least one time in some signature.}
  \item{column no.2}{counts of signature in which the probe/gene is present.}
  \item{column no.3}{mean importance of the probe/gene}
  \item{column no.4}{weighted mean importance}
  \item{column no.5 and more}{importance of the probe/gene with respect to the
signature where it is found}
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{Stefano M. Pagnotta and Michele Ceccarelli}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{signatureFinder}}, \code{\link{puttingAllTogether}},
}
\examples{
## see the source code of the puttingAllTogether() function.
get("puttingAllTogether")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
