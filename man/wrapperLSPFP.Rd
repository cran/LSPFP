\name{wrapperLSPFP}
\alias{wrapperLSPFP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
wrapperLSPFP
}
\description{
THIS IS AN ALPHA BUILD OF 1.0.1
This function plots the positions of peptides with associated proteins from shotgun proteomics data from MaxQuant or Progenesis as input. The plots contain informations about: intensity, position, protein structure, location of protein domains, genename, protein accession, secretion score and truncation score. The plots are written to a PDF file
and a data.frame containing protein feature information is saved as a rds-file and csv-file.
}
\usage{
wrapperLSPFP(globpath, expname, sourcefiles, org, grlocationdf,
            version = "actual", species = c("HUMAN", "MOUSE", "RAT", "PIG"),
            proteomeid = c("UP000005640", "UP000000589", "UP000002494", "UP000008227"),
            taxid = c("9606", "10090", "10116", "9823"), domain = rep("Eukaryota", 4),
            forcedl = FALSE, pepstack = 2, pepque = 2, sortprint = "fcsmall",
            unipep = TRUE, localfasta = "none")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{globpath}{
Character string indicating the path of the global directory.
}
  \item{expname}{
Character string indicating the name of the directory where the files for this run are saved.
%%     ~~Describe \code{expname} here~~
}
  \item{sourcefiles}{
Character string indicating the path to the peptides-file.
%%     ~~Describe \code{sourcefiles} here~~
}
  \item{org}{
Character string specifying the organism the peptides are from.
%%     ~~Describe \code{org} here~~
}
  \item{grlocationdf}{
  Data.frame including the following columns: Expname, Location, Treatment, Sample, Group.
%%     ~~Describe \code{grlocationdf} here~~
}
  \item{version}{
  Character string indicating what version of the \file{BasicData} file should be used.
%%     ~~Describe \code{version} here~~
}
  \item{species}{
  Character strings specifying the UniProt species names of the data sets for download from UniProt to \file{BasicData}.
%%     ~~Describe \code{species} here~~
}
  \item{proteomeid}{
  Character strings specifying the UniProt proteome IDs of the organisms for download. Must be in the same order as \code{species}.
%%     ~~Describe \code{proteomeid} here~~
}
  \item{taxid}{
  Character strings containing the UniProt taxonomic IDs of the organisms for download. Must be in the same order as \code{species}.
%%     ~~Describe \code{taxid} here~~
}
  \item{domain}{
  Character strings containing the UniProt domain description. Must be in the same order as \code{species}.
%%     ~~Describe \code{domain} here~~
}
  \item{forcedl}{
Logical, TRUE: indicates that the actual \code{version} of \file{BasicData} should be downloaded again, FALSE: if not.
}
  \item{pepstack}{
  Numerical indicating the minimal number of runs per group where the same peptide should have been measured.
%%     ~~Describe \code{pepstack} here~~
}
  \item{pepque}{
  Numerical indicating the minimal number of peptides that should have been measured in each groupe.
%%     ~~Describe \code{pepque} here~~
}
  \item{sortprint}{
  Character string that indicates how the peptide plot will be sorted. "fcsmall": decreasing \code{fcsmall} values  , "lflf": decreasing \code{fclf} values, "trunc": decreasing \code{tscore} values, "acc": increasing \code{accessions}, "ntts": decreasing values of \code{NTT}*\code{NTTcov}, "ctts": decreasing values of \code{CTT}*\code{CTTcov}.
%%     ~~Describe \code{transf} here~~
}
  \item{unipep}{
  Logical that indicates if only unique peptides should be used. TRUE: only uniques, FALSE: all peptides in the file.
%%     ~~Describe \code{unipep} here~~
}
 \item{localfasta}{
  Character string indicating that a local FASTA file should be used. Default: "none", if the FASTA file should be downloaded from UniProt, "...": any valid \code{path} to a FASTA file on the system.
%%     ~~Describe \code{unipep} here~~
}

}
\details{

\bold{Worklfow}\cr
1. Create a data structure.\cr
2. Download and prepare data from UniProt.\cr
3. Calculate features.\cr
4. Print plots to PDF.\cr
5. Save features as .rds and .csv.\cr\cr

\bold{Input}\cr
---------------------------------------------------------\cr
peptide-files: the program uses the file extension to decide if the input was created by MaxQuant (.txt) or by Progenesis (.csv).\cr\cr
\bold{\code{grlocationdf}:} Should be a data.frame that contains the following columns:\cr\cr
\bold{\code{Expname}:} This column should contain all experiment names from the peptide file that should be used for the feature plotting. Put each experiment in a single row. It should be a character string. Spaces will be filled with underscores automatically.\cr

\bold{\code{Location}:} Filled with "Secretome" or "Proteome" character strings that indicates if the experiment is from a lysate or not .\cr
\bold{\code{Treatment}:} If the experiments are based on different treatments of the cell, should be marked here. The value must be a character string that matches [A-Z][A-Z].\cr
\bold{\code{Sample}:} Here should be marked from which sample each experiment is from. The values should be simple numbers of type numeric.\cr
\bold{\code{Group}:} These values are used to set experiment groups for the secretome scores. They have to be numbers of type numeric. Set at least two groups if there are to different Locations (e.g. "Secretome" = 1 and "Proteome" = 2). Choose more groups if there are dependencies between treatment, sample and location.\cr

\code{grlocationdf} is used to assign the experiments correctly during the different scoring and plotting functions.\cr\cr

\bold{Datastructure}\cr
At the beginning the following file-structure is created: \cr
\file{globalpath/BasicData}\cr
\file{globalpath/AnalysisData}\cr
If it already exists, nothing new will be created and the old one is used.
In BasicData all the different versions of UniProt download files are stored.
In AnalysisData the output-files are stored.

\bold{UniProt}\cr
For annotation and protein sequence information the organisem specific gff-files and fasta-files are downloaded from UniProt database \url{http://www.uniprot.org/}. As default the actual data-sets are retrieved. If an already downloaded dataset should be used, set \code{version} to the foldername of the existing dataset. No archived version can be downloaded from UniProt.\cr
 \cr

\bold{Output}\cr
The PDF output file will contain a plot for every unique protein a peptide was identified by MaxQuant or Progenesis. \cr

A data.frame containing values that are used for the score calculation will be at: \cr\cr \file{globpath/AnalysisData/expname/feature_table.rds} and \cr\file{globpath/AnalysisData/expname/feature_table.csv} \cr
It contains the following columns: Accession, NTT, NTTcov, CTT, CTTcov, TotalPep, ProtL, MeanSec, MeanProt, MeanSecLF, MeanProtLF, ntts, ctts, tscore, fcsmall, fclf \cr
A .csv-file of \code{grlocationdf} will be saved in \file{AnalysisData}, too.\cr
Further more \file{intenscount_table.csv} and \file{namesdf_table.csv} are stored in the \file{AnalysisData} directory. They
contain information that enables printing with \code{print_selected_peptides}.\cr
Further more there should be two tables (\file{metaIntens.csv}, \file{metaIntens.csv}) that contain different length and intensity values based on annotated extracellular and cytoplasmic regions.\cr

}
\value{
The return value is TRUE if no error occured and FALSE otherwise.
The plots and the feature table can be found in \file{globpath/AnalysisData/expname}.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%%\references{
%% ~put references to the literature/web site here ~
%%}

\author{
Rafael Dellen \cr
\email{Rafael.Dellen@uni-duesseldorf.de}
}
\note{
The download from UniProt may take a while, between 5 and 10 minutes per organism.
A detailed vignette will follow soon.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}

\examples{
#The download of gff-files and FASTA-sequences from UniProt
# might be time consuming.
# Please consider this before running the example.

\dontrun{
#please choose a path
globpath <- getwd()

expname <- "Test_Mouse"

sourcefiles <- system.file("extdata", "Mouse.txt",
                        package="LSPFP")

org <- "Mouse"


#prepare grlocationdf
expnames <- c("Lysat_PB1a","Lysat_PB2a","Lysat_PB3a","Lysat_PB4a",
  "Lysat_PB5a","Lysat_PK1a", "Lysat_PK2a","Lysat_PK3a","Lysat_PK4a",
  "Lysat_PK5a","Sekretom_PB1a","Sekretom_PB2a","Sekretom_PB3a",
  "Sekretom_PB4a","Sekretom_PB5a","Sekretom_PK1a","Sekretom_PK2a",
  "Sekretom_PK3a","Sekretom_PK4a","Sekretom_PK5a")

# Are the values from the secretome or the proteome of the cells?
explocation <- c(rep("Proteome",10),rep("Secretome",10))

# Are the cells from the same culture eg. patient?
expsample <- c(rep(1:5,4))

# Are the samples differently treated?
#(different environments, chemicals, tissue extraction technique)?
exptreatment <- c(rep("AA",5),rep("BB",5),rep("AB",5),rep("BC",5))

#Group specifies which experiments belong together
group <- c(rep(1,10),rep(2,10))
grlocationdf <- data.frame(Expname = expnames, Location = explocation,
                          Treatment = exptreatment, Sample = expsample,
                          Group = group, stringsAsFactors = FALSE)


species <- "MOUSE"
proteomeid <- "UP000000589"
taxid <- "10090"
domain <- "Eukaryota"


res <- wrapperLSPFP(globpath, expname, sourcefiles, org, grlocationdf,
                    species= species, proteomeid = proteomeid,
                    taxid = taxid, domain = domain)

}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}

