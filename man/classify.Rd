\name{classify}
\alias{classify}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Wrapper function to the classification method
%%  ~~function to do ... ~~
}
\description{This function links the classification method to the procedures to
find the seed genes and the signatures. 
If necessary this function can be rewritten in order to use a different classification method. 
Actually the classification method linked is the partitioning around medoids (see pam() function for details)
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
classify(ddata)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ddata}{can be either a list of real values (gene expression levels) or a data-matrix where the rows are the samples and the columns are the genes.
%%     ~~Describe \code{data} here~~
}
}
\details{The function provides two clusters.
%%  ~~ If necessary, more details than the description above ~~
}
\value{If ddata is a list of real values (1 gene expression levels) the function calls an unbiased version of pam (see pamUnbiased() function for details). In this case two elements are returned
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{clusters}{list of 1-2 indicators of the two clusters.}
  \item{missing}{list of T-F logical values where T labels the values set to missing by pamUnbiased().}
%% ...

If ddata is a real matrix (more than 1 expression levels) the function calls pam() and returns an object of pam class with an additional slot
  \item{clusters}{list of 1-2 indicators of the two clusters.}
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

\seealso{ \code{\link{pam}},  \code{\link{pam.object}}.
%% ~~objects to See Also as \code{\link{pam}}, ~~~
}
\examples{
# univariate classification
data(geNSCLC)
sum(is.na(geNSCLC[, "SELP"]))
ans <- classify(geNSCLC[, "SELP"]) 
table(ans$clusters)
sum(ans$missing)
mean(geNSCLC[which(ans$clusters == 1), "SELP"], na.rm = TRUE)
mean(geNSCLC[which(ans$clusters == 2), "SELP"], na.rm = TRUE)

# multivariate classification
data(geNSCLC)
ddata <- geNSCLC[, c("STX1A", "FADD", "STC1", "RNF5")]
ans <- classify(ddata)$clusters
table(ans)
rbind(apply(ddata[ans == 1, ], 2, mean, na.rm = TRUE),
      apply(ddata[ans == 2, ], 2, mean, na.rm = TRUE))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
