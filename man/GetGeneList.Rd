\name{GetGeneList}
\alias{GetGeneList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Function to Filter and Save Genomic Features from NCBI (2018 or newer builds)
}
\description{
\code{GetGeneList} allows the user to access the NCBI database for the species specified using the secure ftp site, download feature information for the latest genome build, as well as filter and save feature information for future use. This is completed using the \pkg{\link{rentrez}} package and will only find the latest assembly for that species. If you are needing an older assembly, see \code{\link{GetGeneList_v11}} function. Once the \code{GetGeneList} function is complete, no other access to NCBI or the internet is required. This function requires user input to determine the feature and class types that will be retained during the filtering process. Note: The requirements for this function have changed slightly due to NCBI ftp site organization changes.
}
\usage{
GetGeneList(Species,savefiles = TRUE, destfile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Species}{
This term designates the species to be used in the function and is dependent on the scientific name.
Options: Must include in quotation marks, where the genus and species should be separated by a space (e.g., "Bos taurus").
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
In running this function, the user will be prompted to enter feedback after the file downloads. Items that will be requested, if multiples are present include 1) primary genome build, 2) primary feature type, and 3) primary class type to prioritize filtering the dataset on. In each case, the user can opt to keep all feature and class types. This will mean that duplicate information is available per gene ID. If filtered, all unique gene ID will be returned, where preference is given to the class feature and class types specified. Gene ID without the preferred feature and class types will be queried for their available information and added while still removing duplicates. The file returned contains 20 columns based on the current NCBI file structure. Those column headings and descriptions are provided below in the \code{Value} section.

Note: While waiting for the function to run, if the user presses "Enter" prematurely, this will result in the function not running correctly and it will have to be started over. Please read instructions carefully.

If \code{savefiles = TRUE}, then both the original file from NCBI and the filtered file the user specified will be saved in the \code{destfile} location. Once the function has run, the user can choose to either use the information at that time or call it later using the saved file. In either case, the output from the filtered file can be used with marker data to run the \code{\link{MapMarkers}} function that is also a part of this package.

}
\value{
Column headings and descriptions returned to the user from the \code{GetGeneList} function.

  \item{feature}{The type of feature based on INSDC, which can include GENE, RNA (various types), and CDS.}
  \item{class}{Gene features are subdivided into classes according to the gene biotype. ncRNA features are subdivided according to the ncRNA_class. CDS features are subdivided into with_protein and without_protein, depending on whether the CDS feature has a protein accession assigned or not. CDS features marked as without_protein include CDS features for C regions and V/D/J segments of immunoglobulin and similar genes that undergo genomic rearrangement, and pseudogenes.}
  \item{assembly}{Accession.version of the assembly.}
  \item{assembly_unit}{The name of the assembly unit, such as "Primary Assembly", "ALT_REF_LOCI_1", or "non-nuclear".}
  \item{seq_type}{The type of sequence the feature is from. Typically include chromosome, mitochondrion, plasmid, or unplaced scaffold.}
  \item{chromosome}{The chromosome the feature is located on, which can include mitochondrial DNA or unknown (blank) if applicable.}
  \item{genomic_accession}{The accession.version of that genome the feature is found on.}
  \item{start}{The start position of the feature on the chromosome.}
  \item{end}{The end position of the feature on the chromosome.}
  \item{strand}{The orientation of the feature on the chromosome (can be + or -).}
  \item{product_accession}{The accession.version of the product referenced by this feature, if it exists.}
  \item{non-redundant_refseq}{For bacteria and archaea assemblies, this column contains the non-redundant WP_ protein accession corresponding to the CDS feature. This may be the same as the previous column for RefSeq genomes annotated directly with WP_ RefSeq proteins, or may be different for genomes annotated with genome-specific protein accessions (e.g. NP_ or YP_ RefSeq proteins) that reference a WP_ RefSeq accession.}
  \item{related_accession}{For eukaryotic RefSeq annotations, this is the RefSeq protein accession corresponding to the transcript feature, or the RefSeq transcript accession corresponding to the protein feature.}
  \item{name}{For genes, this is the gene description or full name. For RNA, CDS, and some other features, this is the product name.}
  \item{symbol}{The gene symbol.}
  \item{GeneID}{The corresponding gene ID on the NCBI database the feature is located in.}
  \item{locus_tag}{No description available from NCBI. Typically a blank column.}
  \item{feature_interval_length}{This is the sum of the lengths of all intervals for the feature (i.e. the length without introns for a joined feature).}
  \item{product_length}{This is the length of the product corresponding to the accession.version in product_accession" column. Protein product lengths are in amino acid units and do not include the stop codon which is included in "feature_interval_length" column. Additionally, product_length may differ from feature_interval_length if the product contains sequence differences vs. the genome, as found for some RefSeq transcript and protein products based on mRNA sequences and also for INSDC proteins that are submitted to correct genome discrepancies.}
  \item{attributes}{A semi-colon delimited list of a controlled set of qualifiers, if available. The list currently includes: partial, pseudo, pseudogene, ribosomal_slippage, trans_splicing, anticodon=NNN (for tRNAs), old_locus_tag=XXX.}
}
\references{
Hulsman Hanna, L. L., and D. G. Riley. 2014. Mapping genomic markers to closest feature using the R package Map2NCBI. Livest. Sci. 162:59-65. doi:\href{https://doi.org/10.1016/j.livsci.2014.01.019}{10.1016/j.livsci.2014.01.019}

National Center for Biotechnology Information. 2018. Latest assembly version 'README' file, last updated 26 February 2018. Available at: \url{https://ftp.ncbi.nlm.nih.gov/genomes/refseq/README.txt} (Accessed 14 April 2019).
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
#choose [n],[1],[n],[3] to filter the build and feature
#information. This example is interactive and requires
#user input. Please note that pressing "Enter" prematurely
#can cause the function to not run properly.
\dontrun{
GeneList = GetGeneList("Bos taurus",savefiles=TRUE,destfile=path.expand("~/"))
}
}
