\name{GetGeneList}
\alias{GetGeneList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Function to Filter and Save Genomic Features from NCBI
}
\description{
GetGeneList allows the user to access the NCBI database for the species specified using the secure ftp site, download feature information for the genome build specified, filter and save feature information for future use. After this function, no other access to NCBI or the internet is required. This function is not limited to only genes, but can also be used for other genomic features like RNA, UTR, and others available for the species specified by the user.
}
\usage{
GetGeneList(Species, build, featuretype = c("GENE", "PSEUDO"), 
            savefiles = FALSE, destfile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Species}{
This term designates the species to be used in the function and is dependent on the scientific name.
Options: Must include in quotation markers and can either separate the genus and species by a space or an underscore (e.g., "Bos taurus" or "Bos_taurus").
}
  \item{build}{
This is a numerical term that designates the species' current genome build to use in the analysis.
Options: Must be in numerical form and not surrounded by quotation marks (e.g., 6.1). These functions may not be compatible with earlier versions of genome builds. To see if a current build is available for a particular species, go to: ftp://ftp.ncbi.nih.gov/genomes/MapView/
}
  \item{featuretype}{
This specifies which feature type(s) the user wants to use and collect. Default setting includes genes and pseudo genes, but user can specify more or less. Each feature type must be in quotation marks. Using multiple feature types must be included as a list (e.g., c("GENE","RNA")).
Options: Can choose to include any of the following: GENE, PSEUDO, RNA, CDS, and UTR.
}
  \item{savefiles}{
Default is false. This term allows you to save the original feature list downloaded from the NCBI database as a text file as well as the filtered feature list produced from the function only if set to TRUE.
Options: Must be either TRUE or FALSE.
}
  \item{destfile}{
This is the pathway to the folder in which files will be saved and must be specified using quotation marks (e.g., "C:/Temp/").
}
}
\details{
In running this function, the user will be prompted to enter feedback after the file downloads to specify the primary genome build to use (if multiple builds are present) as well as the primary feature the user wants to focus on in case there is duplicate information. While waiting for the function to run, if the user presses "Enter" prematurely, this will result in the function not running correctly. 

If savefiles = TRUE, then both the original file from NCBI and the filtered file the user specified will be saved in the destfile location. Once the function has run, the user can choose to either use the information at that time or call it later using the saved file. In either case, the output from the filtered file can be used with marker data to run the MapMarkers function (see separate documentation) that is also a part of this package.

The file returned contains 15 columns based on the current NCBI file structure. Those column headings and descriptions are provided below in the "Value" section.
}
\value{
Column headings and descriptions of the file returned to the user from the "GetGeneList" function.
  \item{tax_id}{Taxonomy id of the species in NCBI.}
  \item{chromosome}{The chromosome the feature is located on in the specified species, which can include mitochondrial DNA if applicable.}
  \item{chr_start}{The start position of the feature on the chromosome.}
  \item{chr_stop}{The stop position of the feature on the chromosome.}
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
Publication of journal article on the Map2NCBI package is pending.
}
\author{
Lauren L. Hulsman Hanna and
David G. Riley
}
\note{
For issues or problems with this function, please contact Lauren Hanna at Lauren.Hanna@ndsu.edu.
}


\seealso{
  Function: MapMarkers
}
\examples{
#Example 1: Run the following example and, when prompted, 
#choose [1], [n], and [1] to filter the build and feature 
#information. This example is interactive and requires 
#user input. Please note that pressing "Enter" prematurely 
#can cause the function to not run properly.
\dontrun{
GeneList = GetGeneList("Bos taurus",build=6.1,savefiles=TRUE,destfile="~/Map2NCBI/data/Test/")
}
}