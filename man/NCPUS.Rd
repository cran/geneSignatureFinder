\name{NCPUS}
\alias{NCPUS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Wrapper function to both makeCluster() and detecCore() functions of the 'parallel' library.
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
\details{The use of this function is depreceted. Use makeCluster() instead.
The function is here for compatibility with the version 2013.08.20 of the library.
}
\author{Stefano M. Pagnotta and Michele Ceccarelli
%%  ~~who you are~~
}

\note{
The definition and connection to the cpu cluster can be established also by
using the makeCluster() native function of the parallel package. 
%%  ~~further notes~~
}

\seealso{
%% ~~objects to See Also as 
\code{%\link{snow}, 
\link{makeCluster}, \link{detectCores}}
}
\examples{
# to count how many cpu's are on the localhost, run
detectCores()

# we suggest to set a cluster of maximal dimension -1 to leave a cpu to control the system
aMakeCluster <- makeCluster(2)
stopCluster(aMakeCluster)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
