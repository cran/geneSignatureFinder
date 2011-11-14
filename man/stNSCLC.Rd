\name{stNSCLC}
\alias{stNSCLC}
\docType{data}
\title{Survival times of 147 samples affected by Non Small Cell Lung Cancer
disease.
}
\description{
This data are of the class "Surv" (defined in the library "survival") and are
built by considering the columns "OS" and "STATUS" of the "Supplementary
Table S1" at the reference below.
The data are defined as stNSCLC <- Surv(OS, STATUS == 1). This dataset is
connected to the geNSCLC data.
}
\usage{data(stNSCLC)}
%%\details{
%%}
%%\source{
%%}
\references{Lau, S.K., P. C. Boutros, M. Pintilie, F. H. Blackhall, C.-Q. Zhu, D. Strumpf, M. R. Johnston, G. Darling, S. Keshavjee, T. K. Waddell, N. Liu, D. Lau, L. Z. Penn, F. A. Shepherd, I. Jurisica, S. D. Der, M.-S. Tsao. A three-gene prognostic classifier for early-stage non-small cell lung cancer. J Clinical Oncology, 25(35): 5562-5569, 2007. 

%%  ~~ possibly secondary sources and usages ~~
}
\seealso{\code{\link{Surv}},\code{\link{geNSCLC}}}

\examples{
data(stNSCLC)
library(survival)
plot(survfit(stNSCLC ~ 1))
}

\keyword{datasets}
