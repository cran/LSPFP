\name{printSelectedPeptides}
\alias{printSelectedPeptides}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Only print peptides
%%  ~~function to do ... ~~
}
\description{
This function allows printing of different peptide plots after running wrapperLSPFP once.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
printSelectedPeptides(path, fname, mysorteddf = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{
  Character vector giving the path to a data folder in the analysis directory.
%%     ~~Describe \code{path} here~~
}
  \item{fname}{
  Character vector containing the name of the new file.
%%     ~~Describe \code{fname} here~~
}
  \item{mysorteddf}{
  A data.frame that can be a sorted and/or shortened data.frame of the \code{feature_table} from \code{path}.
%%     ~~Describe \code{mysorteddf} here~~
}
}
\details{
This function offers the possibility to print more selectively peptide plots after wrapperLSPFP was applied.
To get a smaller plot, it is necessary to commit an existing anlaysis directory from the AnalysisData directory.
The new file is named like \code{fname} and stored in the specified directory \code{path}. If a smaller dataset should be printed, a row wise sorted or shortened \code{feature_table} is to be used as basis for printing but no columns should be removed. This new data.frame (\code{mysorteddf}) will be processed for printing.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A PDF file named like fname in the analysis directory path. Returns TRUE if the printing was successful and FALSE if not.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{
Rafael Dellen \cr
\email{Rafael.Dellen@uni-duesseldorf.de}
%%  ~~who you are~~
}
\note{
A detailed vignette will follow soon.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{

\dontrun{
#To run this example you need to run the example from wrapperLSPFP

path <- paste0(getwd(),"/AnalysisData/Test_Mouse")
test <- read.csv(paste0(path, "/feature_table.csv"))


#Rows can be deleted or sorted but columns should not be removed
#Print the first 25 rows
myprint <- test[1:25,]
#Sort by accession
myprint2 <- test[sort(test$Accession), ]

printSelectedPeptides(path, "newfeaturetable", mysorteddf = myprint)

}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
