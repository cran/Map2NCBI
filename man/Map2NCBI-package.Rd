\name{Map2NCBI-package}
\alias{Map2NCBI-package}
\alias{Map2NCBI}
\docType{package}
\title{
Mapping Markers to the Nearest Genomic Feature
}
\description{
Allows the user to generate a list of features (gene, pseudo, RNA, CDS, and/or UTR) directly from NCBI database for any species with a build available. Option to save downloaded and formatted files is available, and the user can prioritize the feature list based on type and assembly builds present in the build used. The GetGeneList function can now handle query of the NCBI for genome builds released prior to 2018 as well as the latest build for that species. The user can then use the list of features generated or provide a list to map a set of markers (designed for SNP markers with a single base pair position available) to the closest feature based on the map build. This function does require map positions of the markers to be provided and the positions should be based on the build being queried through NCBI.
}
\details{
\tabular{ll}{
Package: \tab Map2NCBI\cr
Type: \tab Package\cr
Version: \tab 1.4\cr
Date: \tab 2020-01-23\cr
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

National Center for Biotechnology Information. 2018. Latest assembly version 'README' file, last updated 26 February 2018. Available at: \url{https://ftp.ncbi.nlm.nih.gov/genomes/refseq/README.txt} (Accessed 23 Jan 2020).
}

\keyword{ package }
\seealso{
  Functions: \code{\link{GetGeneList}} & \code{\link{MapMarkers}}

}
\examples{
#See individual function documentation for applied examples.
}
