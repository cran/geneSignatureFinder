\name{NCPUS}
\alias{NCPUS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Set up a cluster of cpu's on the local machine.
%%  ~~function to do ... ~~
}
\description{The function set up a cluster of 'nchips' cpu's on the local host.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
NCPUS(nchips = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nchips}{is the number of cpu's on the local pc
%%     ~~Describe \code{nchips} here~~
}
}
\value{The function returns a structure useful to establish a connection to the
cluster of cpu's.
}
\references{http://sekhon.berkeley.edu/rgenoud/test/AutoCluster3.R
%% ~put references to the literature/web site here ~
}
\author{Jasjeet S. Sekhon - UC Berkeley
%%  ~~who you are~~
}

\note{
The definition and connection to the cpu cluster can be established also by
using the native function of the snow package. 
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as 
\code{%\link{snow}, 
\link{makeCluster}, \link{detectCores}}
}
\examples{
# to count how many cpu's are on the localhost, run
detectCores()

# we suggest to set a cluster of maximal dimension -1 to leave a cpu to control the system
aNCPUS <- NCPUS(nchips = 2)
stopCluster(aNCPUS)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
