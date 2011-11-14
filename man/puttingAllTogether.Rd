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
puttingAllTogether(workingFile = "", nchips = NULL, useCpuCluster = FALSE, alpha = 0.05, 
saveIndividualSignature = FALSE, saveSurvivalCurvesPlot = FALSE, saveImportancePlot = FALSE, saveRegulationPlot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{workingFile}{prefix string for the log file and other output}
  \item{nchips}{see the NCPUS() function for details; it has to be set to the number of cpu's to use if useCpuCluster = TRUE}
  \item{useCpuCluster}{flag to force the parallel computation}
  \item{alpha}{see the BHcorrection() function}
  \item{saveIndividualSignature}{if TRUE a .RData file is stored in the working directory with the result of the signatureFinder() function}
  \item{saveSurvivalCurvesPlot}{if TRUE a plot of the two survival curves of the signature are stored in a .pdf file in the working directory}
  \item{saveImportancePlot}{if TRUE a plot of the importances of the genes in the signature are stored in a .pdf file in the working directory}
  \item{saveRegulationPlot}{if TRUE a plot of the regulation (up/down) of the genes in the signature are stored in a .pdf file in the working directory}
}
\details{The plotting of the figures may not work if a graphical device is not
instantiated. If the package "WriteXLS" is installed and loaded before running
the function, a number of .xls table are stored in the working directory.}
\value{NULL}
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

# if the library "WriteXLS" is installed,
# some XLS tables will be store in the working directory;
# otherwise they are print to screen.

#puttingAllTogether(
#	workingFile = "NSCLC",
######### modify the number of cpu according to your machine
#	nchips = 2,
#	useCpuCluster = TRUE,
#	saveIndividualSignature = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
