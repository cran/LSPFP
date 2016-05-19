# Secretome Lysat Peptide Descriptor
#
#Function to create a feature.table
#

#Input:character vector accession, character vector printpath, data.frame intenscountdf, logical vector printb  
#Output: logical vector result, data.frame (Accession)
#createFeatureTable <- function(acc, intenscountdf, FastaOrg, AccOrg, GffOrg, locationdf, mybw, printb = TRUE)

createFeatureTable <- function(acc, intenscountdf, FastaOrg, AccOrg, GffOrg, locationdf)
{
  
  #prepare table for feature values
  nracc <- length(acc)
  myvector <- rep(0, nracc)
    
  feature <- data.frame(Accession = acc, NTT = myvector, NTTcov = myvector,
                        CTT = myvector, CTTcov = myvector,  
                        TotalPep = myvector, ProtL = myvector,                                                 
                        MeanSec = myvector, MeanProt = myvector,                
                        MeanSecLF = myvector, MeanProtLF = myvector,                                                                                                                                                     
                        stringsAsFactors = FALSE)
  
  
  #start calculating feature values for every accession
  for (i in feature$Accession)
  {
    index <- feature$Accession == i
    #calculate start and stop positions of the peptides and the according
    # intensities and counts
    icdf <- intenscountdf[(intenscountdf$Accession == i) & (intenscountdf$filtered_out == FALSE), ]
    if(length(icdf[,1]) == 0)
    {
      cat("No peptides for Accession: ", i,"\n")
      next
    }
    
  
    #test
    if (sum(AccOrg$actProteinAcc == i) == 0)
    {
      protlength <- 0 
      cat("No protein in OrgACC for this Accession: ", i,"\n")
      next
    } else {
      
      protlength <- getLength(FastaOrg[[AccOrg[AccOrg$actProteinAcc == i, "index"]]])  
    }
    
    #topological features
    topo_res <- getTopofeatures(GffOrg, protlength, i)
    if(topo_res[[1]] == FALSE)
    {
      #cat("Could not get topological features for Accession: ",i,"\n")
      #next
    }
    topovect <- topo_res[[2]] 
    topodf <- topo_res[[3]] 
        
    #calculate matrices    
    res_matrix <- create_countintensematrix(icdf, protlength, locationdf )
    if(res_matrix[[1]] == FALSE)
    {
      #cat("Could not create matrices for Accession: ", i, "\n")
      next
    }
    countm <- res_matrix[[2]]
    intensm <- res_matrix[[3]]
    
    # set names for the Secretome and the Proteome runs
    namessec <- locationdf[locationdf$Location == "Secretome", "Expname" ]
    namesprot <-locationdf[locationdf$Location == "Proteome", "Expname" ]
            
    # a:= mean of column sum of the secretome intensities
    # b:= mean of column sum of the proteome intensities    
    a <- colMeans(intensm[namessec,])
    b <- colMeans(intensm[namesprot,])
    
    # ac:= mean of column sum of the secretome counts
    # bc:= mean of column sum of the proteome counts
    ac <- colSums(countm[namessec,])/ length(namessec)
    bc <- colSums(countm[namesprot,])/length(namesprot)    
    
    # plot counts /coverage
    #ac1 <- colSums(countm[namessec,])
    #bc1 <- colSums(countm[namesprot,])
    #covall <-ac1+bc1
        
    #---------------------------------------------------
    # feature reading
       
    #kollekt    
    feature[index, "MeanSec"] <- mean(a)
    feature[index, "MeanProt"] <- mean(b)
    
    #---------------------------------------------------------------
    #"Normalisierung" mit der logistischen Funktion   
    tmpLF <- mylogisfunc(intensm, locationdf, i, topodf)  
    
    if(is.list(tmpLF) == FALSE)
    {
      a_lf <- rep(0,protlength)
      b_lf <- rep(0,protlength)
      feature[index, "MeanSecLF"] <- 0
      feature[index, "MeanProtLF"] <- 0
      
    } else {
      
      a_lf <- tmpLF[[1]]
      b_lf <- tmpLF[[2]]
      feature[index, "MeanSecLF"] <- mean(a_lf)
      feature[index, "MeanProtLF"] <- mean(b_lf)
    }
    
    #----------------------------------------------
        
    # detect N and C terminal truncation
    #minpepl <- min(nchar(peptides$Sequence[peptides$Accession == i]))
    minpepl <- min(nchar(intenscountdf$Sequence[intenscountdf$Accession == i]))
    
    #old version with counts
    #res_tt <- detectTerminal(minpepl, ac, bc, i)
    #new version based on intensities
    ia <- intensm[namessec,]
    ib <- intensm[namesprot,]
    ia[ia > 0] <- 1
    ib[ib > 0] <- 1
    abin <- colMeans(ia[namessec,])
    bbin <- colMeans(ib[namesprot,])
    #print(bla)
    res_tt <- detectTerminal(minpepl, abin, bbin, i)
    
    feature[index,"NTT"] <- res_tt[[1]]
    feature[index,"NTTcov"] <- res_tt[[2]]
    feature[index,"CTT"] <- res_tt[[3]]    
    feature[index,"CTTcov"] <- res_tt[[4]]
        
    #Metadata
    # Total Number of Peptides for this Accession    
    feature[index, "TotalPep"] <- length(intenscountdf$Sequence[intenscountdf$Accession == i])
    # Number of Amino Acids for this Accession
    feature[index, "ProtL"] <- protlength
          
  }    
  
  return(feature)
}



# detect if there is a N or C Terminal Truncation
# Input: numerical vector minpepl, ac, bc, character vector acc
# Output: numerical vector ntt, nttcov, ctt, cttcov
# now it it should be based on a binary intensity matrix
# ac:= mean of column sum of the secretome binary intensities
# bc:= mean of column sum of the proteome binary intensities
# in the old version
# ac:= mean of column sum of the secretome counts
# bc:= mean of column sum of the proteome counts
# mimpl := length of the shortest peptidesequence for the actual accession
detectTerminal <- function(minpepl, ac, bc, acc)
{
 
  
  ntt <- 0
  ctt <- 0
  sec <- which(ac >= 0.75)
  #sec <- which(ac >= 0.5)
  if (length(sec) == 0)
  {
    return(list(0,0,0,0))
  }
  
  prot <- which(bc >= 0.75)
  #prot <- which(bc >= 0.5)
  if (length(prot) == 0)
  {
    return(list(0,0,0,0))
  }
  
  #   print(paste0("sec:",sec, "\n"))
  #   print(paste0("prot:",sec, "\n"))
  #   print(acc)
  
  # which returns the indices of a logical object
  # gibt es Signale im Proteomdatensatz die vor dem ersten Signal im
  # Sekretomdatensatz liegen? Sind diese laenger als das kuerzeste Peptid
  # fuer diese Accession?
  
  #N-Terminale Trunkierung
  # welche Position in prot ist kleiner oder gleich der ersten in sec, die einen Mittelwert ueber/= 0.75 hat
  lenprot <- length(which(prot <= sec[1]))
  nttcov <- lenprot/sec[1]
  if (lenprot >= minpepl)
  {
    #N-Terminale Trunkierung moeglich 
    cat("N-Terminale Trunkierung moeglich: ", acc, "\n")
    ntt <- 1    
  }       
  
  #C-Terminale Trunkierung
  clenprot <- length(which(prot >= sec[length(sec)]))
  cttcov <- clenprot/(length(ac)-sec[length(sec)]+1)
  if (clenprot >= minpepl)
  {
    #C-Terminale Trunkierung moeglich  
    cat("C-Terminale Trunkierung moeglich: ", acc,"\n")
    ctt <- 1
  }    
  
  return(list(ntt, nttcov, ctt, cttcov))
}



# function to get the peptide start and stop positions from the protein sequence
# Input: data.frame AccOrg, list FastaOrg, data.frame peptides
# Output: data.frame containing the start and stop position for each peptide and every protein

getPositions <- function(acc, fasta, pep)
{
  #container data.frame
  posdf <- data.frame(Accession = character(0), Sequence = character(0), Start_position = numeric(0), End_position = numeric(0) )
  
  # for every accession in the peptides df
  for (a in unique(pep$Accession))
  {
    # for every sequence for accession a 
    for (s in unique(pep[pep$Accession == a, "Sequence"]) )
    {
      # get the AA sequence of the protein
      index <- acc[acc$actProteinAcc == a, "index"]
      
      if(length(index) == 0)
      {
        print(paste0("getPositions: No index,sequence found for Accession:", a))
        next
      }
      seqprot <- getSequence(fasta[[index]], as.string = TRUE)
      startpos <- regexpr(s, seqprot)
      
      if (startpos[1] == -1)
      {
        next
      } else {
        # calculate the endposition
        # if there is more than one position for the peptide s
        for (st in startpos)
        {
          ep <- st + nchar(s)-1
          newline <- data.frame(Accession = a, Sequence = s, Start_position = st, End_position = ep)
          posdf <- rbind(posdf, newline)
        }
        
      }
      
    }
        
  }

result <- merge(pep, posdf, by = c("Accession", "Sequence"), all.x = TRUE)
return(result)
}



#--------------------------------------------------------------------------
#########################################################
## Load information from FASTA file stored in actFASTA ##
#########################################################
# Input: data.table actFASTAindex, character vector accesionNumber, list actFASTA
# actFASTAindex
#"index","actProteinAcc","actProtein" 
# actFASTA: contains entries like this: sp|P58166|INHBE_HUMAN
#
# Output: data.frame proteinnames
# "Annotation","Name","Sequence","Abrevation","Genename","AccessionFound"
# ACHTUNG actFASTA ist sehr gross, mehr als 100MB

# IDEE: uebergebe nur einen Teil von actFASTAindex 
# actFASTAindex[actProteinAcc==accessionNumber, ]

lspfp_getFastaInformation <- function(actFASTAindex, accessionNumbers, actFASTA)
{

  # create result dataframe
  filler <- rep("NO", length(unique(accessionNumbers)))
  logicalfiller <- rep(FALSE, length(unique(accessionNumbers)))
  result <- data.frame(Accession = unique(accessionNumbers), Annotation = filler,
                       Protein_names = filler, Gene_names = filler, 
                       AccessionFound = logicalfiller, stringsAsFactors = FALSE)
  
  
  # for every accession
  for (acc in unique(accessionNumbers))
  {
    actProtein <- {}
    
    if (length(actFASTAindex[actFASTAindex$actProteinAcc == acc, "actProtein"]) == 0)
    {
      
      print(" --- ")
      print("Check release date of FASTA-file and data set.")
      print(paste0("Accession ", acc, " not found in FASTA-file!"))
      print(" --- ")
      actProtein$AccessionFound <- FALSE
    
    } else {
      actProtein$Annotation <- actFASTAindex[actFASTAindex$actProteinAcc == acc, "actProtein"]
      y <- actFASTAindex[actFASTAindex$actProteinAcc == acc, "index"]
      
      actProtein$Name <- strsplit(strsplit(actProtein$Annotation, "|", 
                                           fixed=TRUE)[[1]][[3]], " OS=", fixed=TRUE)[[1]][[1]]
      
      #actProtein$Sequence <- getSequence(actFASTA[[y]], as.string = TRUE)[[1]]
      
      actProtein$Abrevation <- strsplit(actProtein$Annotation, "GN=", fixed=TRUE)[[1]][2]
      actProtein$Abrevation <- strsplit(actProtein$Abrevation, " ", fixed=TRUE)[[1]][1]
      
      actProtein$Genename <- strsplit(actProtein$Name, " ", fixed=TRUE)[[1]][1]
      actProtein$Genename <- strsplit(actProtein$Genename, "_", fixed=TRUE)[[1]][1]
      
      actProtein$Name <- paste(strsplit(actProtein$Name, " ", fixed=TRUE)[[1]][-1], collapse=" ")
      
      actProtein$AccessionFound <- TRUE
    }
    
  
    print(c(actProtein$Annotation, actProtein$Name, actProtein$Genename, actProtein$AccessionFound))
    result[result$Accession == acc, c("Annotation","Protein_names","Gene_names","AccessionFound")] <- c(actProtein$Annotation, actProtein$Name, actProtein$Genename, actProtein$AccessionFound)
    print(result[result$Accession == acc, c("Annotation","Protein_names","Gene_names","AccessionFound")])
  }
 
  
  return(result)
}
