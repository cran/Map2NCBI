\name{GetGeneList_v11}
\alias{GetGeneList_v11}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Function to Filter and Save Genomic Features from NCBI (builds prior to 2018)
}
\description{
\code{GetGeneList_v11} allows the user to access the NCBI database for the species specified using the secure ftp site, download feature information for the genome build specified (builds prior to 2018 releases), as well as filter and save feature information for future use. After this function, no other access to NCBI or the internet is required. This function is not limited to only genes, but can also be used for other genomic features like RNA, UTR, and others available for the species specified by the user. If you are needing a newer assembly, see \code{\link{GetGeneList}} function. Once the \code{GetGeneList_v11} function is complete, no other access to NCBI or the internet is required. This function requires user input to determine the assembly and features that will be retained during the filtering process.
}
\usage{
GetGeneList_v11(Species, build, featuretype = c("GENE", "PSEUDO"),
                savefiles = TRUE, destfile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Species}{
This term designates the species to be used in the function and is dependent on the scientific name.
Options: Must include in quotation marks and can either separate the genus and species by a space or an underscore (e.g., "Bos taurus" or "Bos_taurus").
}
  \item{build}{
This term designates the species' genome build to use in the analysis.
Options: Must include in quotation marks for the entire build name (e.g., "BUILD.6.1" or "ANNOTATION_RELEASE.104"). These functions may not be compatible with all versions of genome builds, and is not compatible with builds released in 2018 or later. To see if a build is available for a particular species, go to: \url{ftp://ftp.ncbi.nih.gov/genomes/MapView/}, select the species, and click on "sequence".
}
  \item{featuretype}{
This specifies which feature type(s) the user wants to use and collect. Default setting includes genes and pseudo genes, but user can specify more or less. Each feature type must be in quotation marks. Using multiple feature types must be included as an array type (e.g., c("GENE","RNA")).
Options: Can choose to include any of the following: GENE, PSEUDO, RNA, CDS, and UTR.
}
  \item{savefiles}{
Default is true. This term allows you to save the original feature list downloaded from the NCBI database as a text file as well as the filtered feature list produced from the function only if set to TRUE.
Options: Must be either TRUE or FALSE.
}
  \item{destfile}{
This is the pathway to the computer location in which files will be saved and must be specified using quotation marks (e.g., \code{"C:/Temp/"} or \code{paste(getwd(),"/",sep="")}).
}
}
\details{
In running this function, the user will be prompted to enter feedback after the file downloads to specify the primary genome build to use (if multiple builds are present) as well as the primary feature the user wants to focus on in case there is duplicate information. While waiting for the function to run, if the user presses "Enter" prematurely, this will result in the function not running correctly.

If \code{savefiles = TRUE}, then both the original file from NCBI and the filtered file the user specified will be saved in the \code{destfile} location. Once the function has run, the user can choose to either use the information at that time or call it later using the saved file. In either case, the output from the filtered file can be used with marker data to run the \code{\link{MapMarkers}} function that is also a part of this package.

The file returned contains 15 columns based on the older NCBI file structure. Those column headings and descriptions are provided below.
}
\value{
Column headings and descriptions of the file returned to the user from the \code{GetGeneList_v11} function.
  \item{tax_id}{Taxonomy id of the species in NCBI.}
  \item{chromosome}{The chromosome the feature is located on in the specified species, which can include mitochondrial DNA if applicable.}
  \item{start}{The start position of the feature on the chromosome.}
  \item{end}{The end or stop position of the feature on the chromosome.}
  \item{chr_orient}{The orientation of the feature on the chromosome (can be + or -).}
  \item{contig}{The set of overlapping DNA fragments that represent the region of DNA with the same sequence containing the feature.}
  \item{ctg_start}{The start position of the feature on the contig specified.}
  \item{ctg_stop}{The stop position of the feature on the contig specified.}
  \item{ctg_orient}{The orientation of the feature on the contig specified (can be + or -).}
  \item{feature_name}{The NCBI official abbreviation of the feature name.}
  \item{feature_id}{The feature ID on the NCBI database.}
  \item{feature_type}{The type of feature, which can be GENE, PSEUDO, RNA, CDS, and UTR.}
  \item{group_label}{The designated group label on the NCBI database.}
  \item{transcript}{The build in which the feature information is found on.}
  \item{evidence_code}{The evidence code or information, if given.}
}
\references{
Hulsman Hanna, L. L., and D. G. Riley. 2014. Mapping genomic markers to closest feature using the R package Map2NCBI. Livest. Sci. 162:59-65. doi:\href{https://doi.org/10.1016/j.livsci.2014.01.019}{10.1016/j.livsci.2014.01.019}
}
\author{
Lauren L. Hulsman Hanna and
David G. Riley
}
\note{
For issues or problems with this function, please contact Lauren Hanna at \email{Lauren.Hanna@ndsu.edu}.
}


\seealso{
  Function: \code{\link{MapMarkers}}
}
\examples{
#Example 1: Run the following example and, when prompted,
#choose [1], [n], and [1] to filter the build and feature
#information. This example is interactive and requires
#user input. Please note that pressing "Enter" prematurely
#can cause the function to not run properly.
\dontrun{
GeneList = GetGeneList_v11("Bos taurus",build="BUILD.6.1",savefiles=TRUE,destfile=path.expand("~/"))
}
}
