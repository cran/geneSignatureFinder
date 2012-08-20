\name{BICs}
\alias{BICs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute the bayesian information criteria
%%  ~~function to do ... ~~
}
\description{This function computes the bayesian information criteria (bic) under two hypotheses: 1) the data are drawn from a single gaussian; 2) the data are drawn from a mixture of two gaussians.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BICs(ddata, clusters, cutoff = 2.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ddata}{dataset on which the bayesian information criteria have to be computed
%%     ~~Describe \code{data} here~~
}
  \item{clusters}{a list 0-1 of the same length of dataset indicating to which cluster a value belongs
%%     ~~Describe \code{cclusters} here~~
}
  \item{cutoff}{real value that controls the RLS steps for the robust estimates of the parameters of the gaussians
%%     ~~Describe \code{cutoff} here~~
}
}
\details{The estimates of location and scale of the gaussians are obtained by carrying out reweighted least squares (RLS) steps on starting robust estimates. If center0 and deviation0 are the initial robust estimates of the gaussian, to each value in the dataset is assigned weigth = 1 if ((value - center0)/deviation0) < cutoff^2, weigth = 0 otherwise. The center and deviation estimates are updated with those value with weigth = 1. 

The paremeter that controls the mixing of the two gaussians under the second hypothesis is given by mean(clusters).
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{bics}{list of two values: the first element is the bic under the first hypothesis, the second element is the bic under the secod hypothesis}
  \item{bic1Parameters}{list of two values that are the estimates of center and deviation under the first hypothesis}
  \item{bic2Parameters}{list of five values that are the estimates of center and deviation for the first cluster (clusters == 0), center and deviation for the second cluster (clusters == 1), and the estimate of the mixing parameter of the two gaussians.}
%% ...
}
\author{Stefano M. Pagnotta and Michele Ceccarelli
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(geNSCLC)
ans <- classify(geNSCLC[, "STX1A"])
BICs(geNSCLC[, "STX1A"], ans$clusters)

# an example with missing values
data(geNSCLC)
ans <- classify(geNSCLC[, "CALCA"])
BICs(geNSCLC[!ans$missing, "CALCA"], ans$clusters[!ans$missing])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
