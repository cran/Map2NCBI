\name{Map2NCBI-package}
\alias{Map2NCBI-package}
\alias{Map2NCBI}
\docType{package}
\title{
Mapping Markers to the Nearest Genomic Feature
}
\description{
Map2NCBI is an R package that allows the user to generate a list of features (gene, pseudo, RNA, CDS, and/or UTR) directly from NCBI database for any species with a current build available. Option to save downloaded and formatted files is available, and the user can prioritize the feature list based on type and assembly builds present in the current build used. The user can then use the list of features generated or provide a list to map a set of markers (designed for SNP markers with a single base pair position available) to the closest feature based on the map build. This function does require map positions of the markers to be provided and the positions should be based on the build being queried through NCBI.
}
\details{
\tabular{ll}{
Package: \tab Map2NCBI\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-09-24\cr
License: \tab GPL (>= 2)\cr
}
This package can be used as a two part process with the "GetGeneList" function followed by the "MapMarkers" function. See individual function documentation for more information.
}
\author{
Lauren L. Hulsman Hanna and David G. Riley

Maintainer: Lauren Hanna <Lauren.Hanna@ndsu.edu>

}
\references{
Publication of journal article on the Map2NCBI package is pending.
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
\seealso{
  Functions: GetGeneList & MapMarkers

}
\examples{
#See individual function documentation for applied examples.
}
