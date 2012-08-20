\name{puttingAllTogether}
\alias{puttingAllTogether}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Example of strategy to find signature from different seeds at the same
time and to generate plots and summary table.
%%  ~~function to do ... ~~
}
\description{This function has been designed to show how to put the functions of
the package in sequence in order to analyse microarray data for signatures.
Given this function is an example, and can be run for small microarray data, the
user is invited to a reverse engineering task (start by typing
"puttingAllTogether") in order to define a procedure for more complex analysis
sessions. 

}
\usage{
puttingAllTogether(workingFile = "", nchips = 2, alpha = 0.05,
                 saveSurvivalCurvesPlot = FALSE,
                 saveIndividualSignature = FALSE, 
                 saveImportancePlot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{workingFile}{prefix string for the log file and other output}
  \item{nchips}{see the NCPUS() function for details; it has to be set to the number of cpu's to use if useCpuCluster = TRUE}
  \item{alpha}{see the BHcorrection() function}
  \item{saveSurvivalCurvesPlot}{if TRUE a plot of the two survival curves of the signature are stored in a .pdf file in the working directory}
  \item{saveIndividualSignature}{if TRUE a .RData file is stored in the working directory with the result of the signatureFinder() function}
  \item{saveImportancePlot}{if TRUE a plot of the importances of the genes in the signature are stored in a .pdf file in the working directory}
}
\details{The plotting of the figures may not work if a graphical device is not
instantiated. }
\value{Different files can be found in the getwd() directory.}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{Stefano M. Pagnotta and Michele Ceccarelli}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
get("puttingAllTogether")

# uncomment the following lines to performe an analysis
# of the non small cells lung cancer (NSCLC) data 

##################
#data(geNSCLC)
#geData <- geNSCLC
#data(stNSCLC)
#stData <- stNSCLC

#puttingAllTogether(
#	workingFile = "NSCLC",
######### modify the number of cpu according to your machine
#	nchips = 2,
#	saveIndividualSignature = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
