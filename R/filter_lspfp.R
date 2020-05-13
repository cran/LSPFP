# Secretome Lysat Peptide Descriptor
#
# contains filter function

# Input: data.frame peptides, numerical vector pepstack, nuemrical vector pepque
# Output: data.frame peptides with 4 new columns:
#         numerical often,
#         numerical nrpep, 
#         logical filtered_out is true if the peptide is filtered out
#         logical delete_prot is true if there are to less peptides in this protein 

peptide_filter <- function(peptides, grlocationdf, pepstack = 2, pepque = 2)
{
  
  filtered_out <- rep(FALSE, length(peptides$Sequence))
  delete_prot <- rep(FALSE, length(peptides$Sequence))
  mypeptides <- cbind(peptides, filtered_out, delete_prot, stringsAsFactors = FALSE)
  
  # often is an new column that indicates if a peptide occures atleast in one group often enough
  # and is 0 otherwise  
  # counts the number of intensity columns for a peptide
  # how often in each group a peptide was measured
  often <- rep(0, length(mypeptides$Sequence))  
  mypeptides <- cbind(mypeptides, often, stringsAsFactors = FALSE)
  
  for (i in unique(grlocationdf$Group))
  {
    myexp <- grlocationdf$Expname[grlocationdf$Group == i]  
    myexp <- paste0("Intensity_", myexp)  
    intens <- (mypeptides[ , myexp] > 0)
    counter <- (rowSums(intens) >= pepstack)
    mypeptides$often <- mypeptides$often+counter
  }
  
  #mypeptides <- mypeptides[(mypeptides$often > 0), ]
  
  
  # counts the number of peptides for a protein for an experiment in the group
  # when there are at least x peptides in the group
  nrpep <- rep(0, length(mypeptides$Sequence))
  mypeptides <- cbind(mypeptides, nrpep, stringsAsFactors = FALSE)
  
  for (a in unique(mypeptides$Accession))
  {
    for (i in unique(grlocationdf$Group))
    {
      myexp <- grlocationdf$Expname[grlocationdf$Group == i]  
      myexp <- paste0("Intensity_", myexp)  
      mypep <- mypeptides[mypeptides$Accession == a, c("Sequence", myexp)]
      row <- rowSums(mypep[, myexp] > 0)
      if (length(mypep[(row > 0),"Sequence"]) >= pepque)
      {
        mypeptides[mypeptides$Accession == a, "nrpep"] <- mypeptides[mypeptides$Accession == a, "nrpep"]+1
      }
    }
    
  }
  #myfltdpeptides <- mypeptides[(mypeptides$nrpep > 0), ]
  #
  mypeptides$filtered_out <- !((mypeptides$nrpep > 0) | (mypeptides$often > 0))
  
  
  #test if the measurements for a protein contain enough valid peptides
  for(b in unique(mypeptides$Accession))
  {
    #if (b == "A2A432") print(bla)
    nrseq <- length(unique(mypeptides[mypeptides$Accession == b, "Sequence"]))
    
    #are there more than two peptides measured for this protein?
    # if not mark the protein/Accession for deletion
    if (nrseq < 2)
    {
      mypeptides[mypeptides$Accession == b, "delete_prot"] <- TRUE
    
    } else {
      
       intenscol <- grep("Intensity_", names(mypeptides))
       binintens <- mypeptides[mypeptides$Accession == b, intenscol]
       binintens[binintens > 0] <- 1
       binintens[binintens < 1] <- 0
       #binintens <- sapply(binintens, as.numeric)
       rs <- rowSums(binintens)
       #are there peptides that were measured in at least two experiments?
       r_nr <- sum(rs >= 2)
       
       #where less than two peptides measured twice in all experiments?        
       # if true mark the protein/Accession for deletion
       if(r_nr < 2)
       {
         mypeptides[mypeptides$Accession == b, "delete_prot"] <- TRUE
       }
         
    }
    
  }
  
  
  return(mypeptides)
  
}