\name{geNSCLC}
\alias{geNSCLC}
\docType{data}
\title{Gene expression levels on 147 samples affected by Non Small Cell Lung Cancer disease.
%%   ~~ data name/kind ... ~~
}
\description{This is a matrix is built from the table "Normalized RT qPCR data" available at the source below. Each row is a sample and each column is a gene expression levels.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(geNSCLC)}
\format{The columns are the 158 genes and the 147 rows are the samples. This dataset is connected to the stNSCLC data.}
%%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%%}
\source{http://www.cs.toronto.edu/~juris/data/JCO07/
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{Lau, S.K., P. C. Boutros, M. Pintilie, F. H. Blackhall, C.-Q. Zhu, D. Strumpf, M. R. Johnston, G. Darling, S. Keshavjee, T. K. Waddell, N. Liu, D. Lau, L. Z. Penn, F. A. Shepherd, I. Jurisica, S. D. Der, M.-S. Tsao. A three-gene prognostic classifier for early-stage non-small cell lung cancer. J Clinical Oncology, 25(35): 5562-5569, 2007. 

%%  ~~ possibly secondary sources and usages ~~
}
\seealso{\code{\link{stNSCLC}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(geNSCLC)
dim(geNSCLC)
head(colnames(geNSCLC))
}
\keyword{datasets}
