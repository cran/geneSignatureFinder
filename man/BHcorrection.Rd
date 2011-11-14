\name{BHcorrection}
\alias{BHcorrection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Benjamini & Hochberg (1995) method for p-values correction
%%  ~~function to do ... ~~
}
\description{This function evaluate the p-values correction according to the Benjamini & Hochberg (1995) method. 
}
\usage{
BHcorrection(pvs, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pvs}{list of p-values
%%     ~~Describe \code{pvs} here~~
}
  \item{alpha}{level of significance of the test
%%     ~~Describe \code{alpha} here~~
}
}
\value{list of corrected p-values
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Benjamini, Y. and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. J. Roy. Statist. Soc. Ser. B 57 289-300. 
%% ~put references to the literature/web site here ~
}
\author{Stefano M. Pagnotta and Michele Ceccarelli
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
###################
data <- matrix(rnorm(50000), 10000, 5)
pvs <- apply(data, 1, function(xx) t.test(xx, alternative = "less")$p.value)
qvs <- BHcorrection(pvs)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
