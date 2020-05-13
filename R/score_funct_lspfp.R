# Secretome Lysat Peptide Descriptor
#
# scoring functions

# calculate truncation score
# Input: data.frame featuredf
# Output: a list of numerical vectors 
# NTT	NTTcov	CTT	CTTcov
truncscoring <- function(ftdf)
{
  #ts <- rep(0, length(ftdf$Accession))
  
  #ts <- (ftdf$NTT*ftdf$NTTcov)+(ftdf$CTT*ftdf$CTTcov)
  
  ntts <- (ftdf$NTT*ftdf$NTTcov)
  
  ctts <- (ftdf$CTT*ftdf$CTTcov)
  
  ts <- ntts + ctts
  
  ts2 <- ntts + ctts + ftdf$SecCutCov
  
  return(list(ts, ntts, ctts, ts2))
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