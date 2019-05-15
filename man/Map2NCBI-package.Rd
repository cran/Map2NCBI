\name{Map2NCBI-package}
\alias{Map2NCBI-package}
\alias{Map2NCBI}
\docType{package}
\title{
Mapping Markers to the Nearest Genomic Feature
}
\description{
Allows the user to generate a list of features (gene, pseudo, RNA, CDS, and/or UTR) directly from NCBI database for any species with a current build available. Option to save downloaded and formatted files is available, and the user can prioritize the feature list based on type and assembly builds present in the current build used. Due to changes in NCBI ftp site organization, two functions to generate the list of features are available. The original GetGeneList function (now called \code{\link{GetGeneList_v11}}) will query the NCBI for genome builds released prior to 2018. The new \code{\link{GetGeneList}} function will query the NCBI for genome builds released in 2018 or later (e.g., UCD 1.2 for /emph{Bos taurus}). The user can then use the list of features generated or provide a list to map a set of markers (designed for SNP markers with a single base pair position available) to the closest feature based on the map build. This function does require map positions of the markers to be provided and the positions should be based on the build being queried through NCBI.
}
\details{
\tabular{ll}{
Package: \tab Map2NCBI\cr
Type: \tab Package\cr
Version: \tab 1.3\cr
Date: \tab 2019-05-15\cr
License: \tab GPL (>= 2)\cr
}
This package can be used as a two part process with the \code{\link{GetGeneList}} function followed by the \code{\link{MapMarkers}} function. See individual function documentation for more information.
}
\author{
Lauren L. Hulsman Hanna and David G. Riley

Maintainer: Lauren Hanna <\email{Lauren.Hanna@ndsu.edu}>

}
\references{
Hulsman Hanna, L. L., and D. G. Riley. 2014. Mapping genomic markers to closest feature using the R package Map2NCBI. Livest. Sci. 162:59-65. doi:\href{https://doi.org/10.1016/j.livsci.2014.01.019}{10.1016/j.livsci.2014.01.019}

National Center for Biotechnology Information. 2018. Latest assembly version 'README' file, last updated 26 February 2018. Available at: \url{https://ftp.ncbi.nlm.nih.gov/genomes/refseq/README.txt} (Accessed 14 April 2019).
}
\keyword{ package }
\seealso{
  Functions: \code{\link{GetGeneList}} & \code{\link{MapMarkers}}

}
\examples{
#See individual function documentation for applied examples.
}
