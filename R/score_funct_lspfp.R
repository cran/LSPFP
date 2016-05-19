# Secretome Lysat Peptide Descriptor
#
# scoring functions

# calculate truncation score
# Input: data.frame featuredf
# Output: numerical vector ts
# NTT	NTTcov	CTT	CTTcov
truncscoring <- function(ftdf)
{
  ts <- rep(0, length(ftdf$Accession))
  
  ts <- (ftdf$NTT*ftdf$NTTcov)+(ftdf$CTT*ftdf$CTTcov)
  
  return(ts)
}


# calculate secretion scores
# Input: data.frame featuredf 
# Output: a list of numerical vectors 
secscoring <- function(ftdf)
{   
  
  fcsmall <- (ftdf$MeanSec+0.01)/(ftdf$MeanProt+0.01)
  fclf <- ftdf$MeanSecLF/ftdf$MeanProtLF  
   
  result <- list(fcsmall, fclf)
 
  return(result)
}