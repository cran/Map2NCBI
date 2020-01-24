\name{MapMarkers}
\alias{MapMarkers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mapping SNP Markers to Closest Genomic Feature
}
\description{
\code{MapMarkers} allows the user to map the supplied DNA markers (primarily designed for SNP markers) to the genomic feature in closest proximity based on the feature list generated using the \code{\link{GetGeneList}} function or a properly formated feature list (see \code{Values} section).
}
\usage{
MapMarkers(features, markers, nAut, other = c("X"), savefiles = TRUE, destfile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{features}{
This is the table or matrix in the current R session that will be used to map the marker list to. If using the \code{GetGeneList} function, the name given to the output should be supplied here (e.g., "GeneList" from the example provided in the \code{\link{GetGeneList}} documentation file). Note: If using a feature list generated using version 1.1 of this package, please review column names below for changes in format to ensure the function runs properly.
}
  \item{markers}{
This is the table or matrix in the current R session that will provide marker map information to use for the function. See \code{Values} for format of the marker file needed.
}
  \item{nAut}{
The number of autosomes in the species. This should reflect the total number of autosomes in the species, not the number of autosomes in the marker file.
}
  \item{other}{
The sex chromosomes or other genomic information available (e.g., for eukaryotes this could include mitochondrial DNA). These must be specified inside quotation marks. If sex chromosomes or other genomic information is not provided in the marker file, set other=FALSE.
}
  \item{savefiles}{
Default is TRUE. This term allows you to save the final marker file with genomic feature information in the destfile location as "MappedMarkers.txt" format. Any markers that cannot be mapped due to lack of feature information are saved as "NotMapped.txt". Options: Must be either TRUE or FALSE.
}
  \item{destfile}{
This is the pathway to the folder in which files will be saved and must be specified using quotation marks (e.g., \code{"C:/Temp/"} or \code{paste(getwd(),"/",sep="")}).
}
}
\details{
The \code{MapMarkers} function processes each chromosome individually to search for features that fall closest to the markers provided based on the map information included. Map positions of the markers must match the assembly being used in the feature list. Once the closest feature has been found, the marker and feature information are saved together and take the format of binding the marker map file (which include at a minimum 3 columns) with the feature list columns provided (20 columns if using the \code{GetGeneList} function or a minimum of 4 columns if formatting yourself). The function also adds 2 additional columns described in \code{Value} section to identify the distance the marker is from the feature and a category to group the marker's proximity to the feature by.

}
\value{
1) Format for feature list if not generated using the \code{GetGeneList} function:
    \item{FeatureName}{The name of the feature provided. Column heading name can be changed, but should be included to identify the feature once the \code{MapMarkers} function is completed.}
    \item{chromosome}{The chromosome in which the genomic feature is located on. The column heading name must be given this name. If including sex chromosomes or other genomic information, label based on letters or abbreviation (e.g., "X").}
    \item{start}{The start position of the genomic feature based on the build used. The column heading name must be given this name. This used to be called "chr_start" in version 1.1 of this package.}
    \item{end}{The end or stop position of the genomic feature based on the build used. The column heading name must be given this name. This used to be called "chr_stop" in version 1.1 of this package.}


2) Format for the marker map file:
    \item{Marker}{Name of the marker. Be aware of R language and its restrictions. The name of this column heading can be changed to something else.}
    \item{chromosome}{The chromosome in which the marker is mapped to. The name of this column is required and must be exact. This must be numeric. If including sex chromosomes or other genomic information, assign numbers to each. Number the sex chromosomes or other genomic information in the order that matches the order listed in the other=c() statement (e.g., X and Y chromosomes are labeled 30 and 31, respectively, so other=c("X","Y") to follow that order). The function will automatically align the letter with the correct number as long as they are included in the order specified.}
    \item{position}{The base pair position of the marker based on the map build used. This build must also match the build in which you generated genomic feature from using the \code{GetGeneList} function or other method. The name of this column is required and must be exact. This function was designed for SNP markers, but if using other types of markers, you should choose the base pair location that best represents the marker (e.g., center position) and include that in this column.}

NOTE: Order of the columns in both files are not necessarily important, but correct column heading names are essential. R programming is case sensitive, so make sure it matches exactly unless otherwise noted. Other columns may be included, but will not be used by the function. Any columns included in this file will be returned with the final marker file after the \code{MapMarkers} function is completed.

3) Additional columns included in the output file of the \code{MapMarkers} function:
    \item{Distance}{The base pair distance of the marker from the closest feature identified. If the marker is located inside the feature, the distance is set to zero.}
    \item{Inside?}{The category in which the marker and feature pair fall into. This is based on the distance between the Marker and the closest feature, which is broken into 11 categories described in the next section.}

4) Categories that are included in the "Inside?" column:
    \item{Yes,_Inside_Gene}{Marker is located in the closest feature.}
    \item{Marker_is_<=_2500_bp_Before_Feature}{The closest feature is located after the marker position and is within 2,500 base pairs (bp).}
    \item{Marker_is_<=_2500_bp_After_Feature}{The closest feature is located before the marker position and is within 2,500 bp.}
    \item{Marker_is_>_2500_bp_<=5000_bp_Before_Feature}{The closest feature is located before the marker position and is between 2,500 bp and 5,000 bp from the marker.}
    \item{Marker_is_>_2500_bp_<=5000_bp_After_Feature}{The closest feature is located after the marker position and is between 2,500 bp and 5,000 bp from the marker.}
    \item{Marker_is_>_5000_bp_<=25000_bp_Before_Feature}{The closest feature is located before the marker position and is between 5,000 bp and 25,000 bp from the marker.}
    \item{Marker_is_>_5000_bp_<=25000_bp_After_Feature}{The closest feature is located after the marker position and is between 5,000 bp and 25,000 bp from the marker.}
    \item{Nearest_feature_is_>_25,000_bp_before_marker}{The closest feature is located before the marker position and is more than 25,000 bp from the marker.}
    \item{Nearest_feature_is_>_25,000_bp_after_marker}{The closest feature is located after the marker position and is more than 25,000 bp from the marker.}
    \item{Nearest_feature_is_>_1_Mb_before_marker}{The closest feature is located before the marker position and is more than 1,000,000 bp (1 Mb) from the marker.}
    \item{Nearest_feature_is_>_1_Mp_after_marker}{The closest feature is located after the marker position and is more than 1,000,000 bp (1 Mb) from the marker.}

}
\references{
Hulsman Hanna, L. L., and D. G. Riley. 2014. Mapping genomic markers to closest feature using the R package Map2NCBI. Livest. Sci. 162:59-65. doi:\href{https://doi.org/10.1016/j.livsci.2014.01.019}{10.1016/j.livsci.2014.01.019}
}
\author{
Lauren L. Hulsman Hanna and David G. Riley
}
\note{
For issues or problems with this function, please contact
Lauren Hanna at \email{Lauren.Hanna@ndsu.edu}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 Function: \code{\link{GetGeneList}}
}
\examples{
#Example 1: Step 1 includes running "GetGeneList" function.
#As this step is interactive, a dataset from Bos taurus has
#been generated and available to use in the \data folder as
#well as a subset of marker information from BTA 1. Use the
#following code to run this example:

data(GeneList_BTA1)
data(Example10MarkerFile)
Example1 = MapMarkers(GeneList_BTA1, Example10MarkerFile,
    nAut=29,other="X",savefiles = FALSE)

#Note, this example will not save the output to the working
#directory, but will return the information to "Example1"
#variable.
}
