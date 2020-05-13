# Secretome Lysat Peptide Descriptor

# math functions

# Prepare Data for Matrices

# create intensity-matrix for one Accession/Protein
# Input: data.frame df (Accession, Start, End,Counts, Intensities), numeric vector protlength, data.frame location
# Output: list (logical vector,  numeric matrix count , numeric matrix intensity )
create_countintensematrix<-function(df, protlength, location)
{
  rows <- length(location[,"Expname"])
  #   columnIntens <- grep("Intensity_",names(df))
  #   columnCount <- grep("Counts_",names(df))  
  columnIntens <- names(df) %in% paste0("Intensity_",location[,"Expname"])
  columnCount <- names(df) %in% paste0("Counts_",location[,"Expname"])
  countm <- matrix(c(0),nrow= rows, ncol= protlength)
  rownames(countm) <- location[,"Expname"]
  intensm <- matrix(c(0),nrow= rows, ncol= protlength)
  rownames(intensm) <- location[,"Expname"]
  #set all Intensities that are -Inf to 0
  df[which(df == -Inf, arr.ind= TRUE)] <- 0
  #   > colnames(intenscountdf)
  #   [1] "Accession"             "Sequence"              "Start_position"        "End_position"          "Counts_01_conred"     
  #   [6] "Counts_02_EPSred"      "Counts_03_conred"      "Counts_04_EPSred"      "Counts_05_conred"      "Counts_06_EPSred"     
  #   [11] "Counts_07_conblue"     "Counts_08_EPSblue"     "Counts_09_conblue"     "Counts_10_EPSblue"     "Counts_11_conblue"    
  #   [16] "Counts_12EPSblue"      "Counts_Lysat54con1"    "Counts_Lysat54con2"    "Counts_Lysat54EPS1"    "Counts_Lysat54EPS2"   
  #   [21] "Counts_Lysat70con1"    "Counts_Lysat70con2"    "Counts_Lysat70EPS1"    "Counts_Lysat70EPS2"    "Intensity_01_conred"  
  #   [26] "Intensity_02_EPSred"   "Intensity_03_conred"   "Intensity_04_EPSred"   "Intensity_05_conred"   "Intensity_06_EPSred"  
  #   [31] "Intensity_07_conblue"  "Intensity_08_EPSblue"  "Intensity_09_conblue"  "Intensity_10_EPSblue"  "Intensity_11_conblue" 
  #   [36] "Intensity_12EPSblue"   "Intensity_Lysat54con1" "Intensity_Lysat54con2" "Intensity_Lysat54EPS1" "Intensity_Lysat54EPS2"
  #   [41] "Intensity_Lysat70con1" "Intensity_Lysat70con2" "Intensity_Lysat70EPS1" "Intensity_Lysat70EPS2"
  
  for (s in seq_along(df$Sequence))
  {
    # in some cases th peptide position does not match the length of the protein from the UniProt HumanFasta File
    if (df$Start_position[s] > protlength | df$End_position[s] > protlength)
    {
      cat("For Accession: ", df$Accession[1], " one peptide is out of range of the protein!\n")
      next
    }
    else
    {
      countm[1:rows, df$Start_position[s]:df$End_position[s]] <- countm[1:rows, df$Start_position[s]:df$End_position[s]] + unlist(df[s, columnCount])
      intensm[1:rows, df$Start_position[s]:df$End_position[s]] <- intensm[1:rows, df$Start_position[s]:df$End_position[s]] + unlist(df[s, columnIntens])
    }
    
  }
  
  
  result <- TRUE
  return(list(result, countm, intensm))
}


#logistische Funktion, transformiert Werte groesser Null in das Interval [0.5,1]
# und berechnet den Mean fuer die gesamte location
#Input: numerical Vector x, data.frame loc that contains 
# "Expname"   "Location"  "Treatment" "Sample"
# Output: list of (numerical Vector ms and mp)(MeanSecLF,MeanProtLF) or FALSE 
mylogisfunc <- function(x, loc, acc, topodf)
{
  if(!is.numeric(x))
  {
    print("x is not numeric")
    return(FALSE)
  }
  if(sum(is.na(loc)) >= 1)
  {
    print("Location data frame contains NA")   
    return(FALSE)
  }
  
  # split matrix by treatment
  
  for (t in unique(loc$Treatment))
  {
    for (l in unique(loc$Location))
    {
      #
      index <- (loc$Location == l) & ( loc$Treatment == t)
      mym <- x[index, ]
      # only columns where an intensity was measured
      if (is.data.frame(mym) == TRUE)
      {
        colindex <- (colSums(mym) > 0) 
        # Vergleiche alle Spalten mit nur den Spalten, in denen eine Intensitaet gemessen wurde
        # all columns 
        
        mym <- apply(mym, 2, function(y) 1/(1+exp(-y)))                   
        x[index, ] <- mym
      }
      else
      {
        mym <- 1/(1+exp(-mym))
        x[index, ] <- mym
      }
      
    }
  }
  
  #calculate means
  mp <- colMeans(x[loc$Location == "Proteome", ])
  ms <- colMeans(x[loc$Location == "Secretome", ])
  
  #gebe die Means beider Bereiche zurueck
  result <- list(ms, mp)
  return(result)
}


