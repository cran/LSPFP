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
                        CTT = myvector, CTTcov = myvector, SecCutCov = myvector, SecProtCov = myvector,
                        TotalPep = myvector, ProtL = myvector,
                        MeanSec = myvector, MeanProt = myvector,
                        MeanSecLF = myvector, MeanProtLF = myvector,
                        stringsAsFactors = FALSE)

  # set names for the Secretome and the Proteome runs
  namessec <- locationdf[locationdf$Location == "Secretome", "Expname" ]
  namesprot <-locationdf[locationdf$Location == "Proteome", "Expname" ]

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
    feature[index,"SecProtCov"] <- res_tt[[5]]
    feature[index,"SecCutCov"] <- res_tt[[6]]

    #Metadata
    # Total Number of Peptides for this Accession
    feature[index, "TotalPep"] <- length(intenscountdf$Sequence[intenscountdf$Accession == i])
    # Number of Amino Acids for this Accession
    feature[index, "ProtL"] <- protlength

  }

  return(feature)
}


addMetaScoreFeature <- function(acc, intenscountdf, FastaOrg, AccOrg, GffOrg, locationdf)
{
  # print(bla)
  #prepare table for feature values
  nracc <- length(acc)
  nrexp <- length(locationdf$Expname)

  # Accession, Expname, Location, IntensityEC, IntensityIP
  # Accession, Topomatch, LengthEC, LengthIP
  myrownr <- rep("", nracc * nrexp)
  metaIntens <- data.frame(Accession = myrownr, Expname = myrownr, Location = myrownr, IntensityEC = myrownr,
                           IntensityEClength = myrownr, IntensityCP = myrownr, IntensityCPlength = myrownr,
                           stringsAsFactors = FALSE)
  mynracc <- rep("", nracc)
  metaLength <- data.frame(Accession = mynracc, Topomatch = rep(FALSE, nracc), LengthEC = mynracc,
                           LengthCP = mynracc, stringsAsFactors = FALSE)

  names(metaIntens)
  names(metaLength)
  head(metaIntens)
  head(metaLength)


  names(intenscountdf)
  head(intenscountdf)
  # [1] "Accession"             "Sequence"              "Start_position"        "End_position"          "filtered_out"          "Counts_lysate_1"
  # [7] "Counts_lysate_2"       "Counts_lysate_3"       "Counts_lysate_4"       "Counts_lysate_5"       "Counts_secretome_1"    "Counts_secretome_2"
  # [13] "Counts_secretome_3"    "Counts_secretome_4"    "Counts_secretome_5"    "Intensity_lysate_1"    "Intensity_lysate_2"    "Intensity_lysate_3"
  # [19] "Intensity_lysate_4"    "Intensity_lysate_5"    "Intensity_secretome_1" "Intensity_secretome_2" "Intensity_secretome_3" "Intensity_secretome_4"
  # [25] "Intensity_secretome_5"

  # names(FastaOrg)
  # FastaOrg[[1]]
  # contains Fasta Sequences
  # [1] "SLTIIVSSVLRCQASELHVGCVGGCVYSKCGLAFELLIVGNTTIFRMNSSNPFQSLLQTLLHPLNFIVIKIVARSLAGQLHLSPRYYQLWSRLLKLDL"
  # attr(,"name")
  # [1] "tr|A0A075B5I0|A0A075B5I0_MOUSE"
  # attr(,"Annot")
  # [1] ">tr|A0A075B5I0|A0A075B5I0_MOUSE Protein Ighv1-13 OS=Mus musculus GN=Ighv1-13 PE=4 SV=1"
  # attr(,"class")
  # [1] "SeqFastaAA"

  names(AccOrg)
  # [1] "index"         "actProteinAcc" "actProtein"
  head(AccOrg)
  # index actProteinAcc                                                                                            actProtein
  # 1     1    A0A075B5I0                >tr|A0A075B5I0|A0A075B5I0_MOUSE Protein Ighv1-13 OS=Mus musculus GN=Ighv1-13 PE=4 SV=1
  # 2     2    A0A075B5I2           >tr|A0A075B5I2|A0A075B5I2_MOUSE Protein Trbv4 (Fragment) OS=Mus musculus GN=Trbv4 PE=4 SV=1
  # 3     3    A0A075B5J4           >tr|A0A075B5J4|A0A075B5J4_MOUSE Protein Trbc2 (Fragment) OS=Mus musculus GN=Trbc2 PE=4 SV=1
  # 4     4    A0A075B5J5                    >tr|A0A075B5J5|A0A075B5J5_MOUSE Protein Trbv31 OS=Mus musculus GN=Trbv31 PE=4 SV=1
  # 5     5    A0A075B5J6                  >tr|A0A075B5J6|A0A075B5J6_MOUSE Protein Gm20730 OS=Mus musculus GN=Gm20730 PE=4 SV=7
  # 6     6    A0A075B5K0 >tr|A0A075B5K0|A0A075B5K0_MOUSE Protein Igkv14-126 (Fragment) OS=Mus musculus GN=Igkv14-126 PE=4 SV=7

  names(GffOrg)
  head(GffOrg)
  # Browse[1]> names(GffOrg)
  # [1] "seqname"    "source"     "feature"    "start"      "end"        "score"      "strand"     "frame"      "attributes" "Note"       "ID"

  names(locationdf)
  # Browse[1]> names(locationdf)
  # [1] "Expname"   "Location"  "Treatment" "Sample"    "Group"


  #start calculating feature values for every accession for metaLength
  #[1] "Accession" "Topomatch" "LengthEC"  "LengthCP"
  indexML <- 1
  indexIn <- 1
  lenLoc <- length(locationdf[,1])

  for (a in acc)
  {
    # a <- acc[1]
    #a <- "Q8VDN2"
    #test
    if (sum(AccOrg$actProteinAcc == a) == 0)
    {
      protlength <- 0
      cat("No protein in OrgACC for this Accession: ", a,"\n")
      next

    } else{
      protlength <- getLength(FastaOrg[[AccOrg[AccOrg$actProteinAcc == a, "index"]]])
    }

    #topological features
    topo_res <- getTopofeatures(GffOrg, protlength, a)
    if(topo_res[[1]] == FALSE)
    {
      #cat("Could not get topological features for Accession: ",i,"\n")
      #next
    }
    topovect <- topo_res[[2]]
    topodf <- topo_res[[3]]

    if(sum(grepl("CP|EC",topodf$feature)) == 0)
    {
      print("No fitting topological feature found for MetaScore")
      metaLength[indexML, "Topomatch"] <- FALSE
      #next
    } else{
      metaLength[indexML, "Topomatch"] <- TRUE
    }

    topomatch <- grep("EC", topodf$feature)
    if (length(topomatch) > 0)
    {
      metaLength[indexML, "LengthEC"] <- sum((topodf$end[topomatch]-topodf$start[topomatch])+1)
    } else {
      metaLength[indexML, "LengthEC"] <- 0
    }

    topomatch <- grep("CP", topodf$feature)
    if (length(topomatch) > 0)
    {
      metaLength[indexML, "LengthCP"] <- sum((topodf$end[topomatch]-topodf$start[topomatch])+1)
    } else {
      metaLength[indexML, "LengthCP"] <- 0
    }

    metaLength[indexML, "Accession"] <- a

    # fill in metaIntens
    if (metaLength[indexML, "Topomatch"] == TRUE)
    {
      endpos <- (indexIn+lenLoc-1)
      metaIntens[indexIn:endpos, "Accession"] <- rep(a, lenLoc)
      metaIntens[indexIn:endpos, "Expname"] <- locationdf$Expname
      metaIntens[indexIn:endpos, "Location"] <- locationdf$Location

      #calculate matrices
      icdf <- intenscountdf[(intenscountdf$Accession == a) & (intenscountdf$filtered_out == FALSE), ]
      if(length(icdf[,1]) == 0)
      {
        cat("No peptides for Accession: ", a,"\n")
        metaIntens[indexIn:endpos, "IntensityEC"] <- 0
        metaIntens[indexIn:endpos, "IntensityCP"] <- 0
        next
      }

      #calculate matrices
      res_matrix <- create_countintensematrix(icdf, protlength, locationdf )
      if(res_matrix[[1]] == FALSE)
      {
        cat("Could not create matrices for Accession: ", a, "\n")
        metaIntens[indexIn:endpos, "IntensityEC"] <- 0
        metaIntens[indexIn:endpos, "IntensityCP"] <- 0
        next
      }
      countm <- res_matrix[[2]]
      intensm <- res_matrix[[3]]

      topodf
      intensm
      # start calculating intensities
      itopEC <- topovect == "EC"
      itopCP <- topovect == "CP"
      for(t in indexIn:endpos )
      {
        ex <- metaIntens[ t , "Expname"]
        #sum(intensm["lysate_1", itopEC])
        metaIntens[t, "IntensityEC"] <- sum(2^intensm[ex, itopEC])
        metaIntens[t, "IntensityCP"] <- sum(2^intensm[ex, itopCP])
        metaIntens[t, "IntensityEClength"] <- sum(2^intensm[ex, itopEC] > 0)
        metaIntens[t, "IntensityCPlength"] <- sum(2^intensm[ex, itopCP] > 0)
      }

      indexIn <- endpos + 1

    }else{
      next
    }

    indexML <- indexML + 1
  }

  metaLength <- metaLength[!metaLength$Accession == "", ]
  metaIntens <- metaIntens[!metaIntens$Accession == "", ]

  result <- list(metaLength, metaIntens)
  return(result)
}


# Version vor dem 19.09.2016
addMetaScoreFeature_old <- function(acc, intenscountdf, FastaOrg, AccOrg, GffOrg, locationdf)
{

  #prepare table for feature values
  nracc <- length(acc)
  #myvector <- rep(0, nracc)

  # set names for the Secretome and the Proteome runs
  namessec <- locationdf[locationdf$Location == "Secretome", "Expname" ]
  namesprot <-locationdf[locationdf$Location == "Proteome", "Expname" ]

  # set table for a metainformation based score
  #print(bla)
  tmpmeta <- matrix(data = 0, nrow = nracc, ncol = length(locationdf$Location))
  colnames(tmpmeta) <- sort(locationdf$Expname)

  # tmpmeta2 will contain all
  tmpmeta2 <- matrix(data = 0, nrow = nracc, ncol = length(locationdf$Location))
  colnames(tmpmeta2) <- sort(paste0("Meta2_", locationdf$Expname))

  metascore_int <- data.frame(Accession = acc, Topomatch = rep(FALSE, nracc), tmpmeta, tmpmeta2, stringsAsFactors = FALSE)


  #start calculating feature values for every accession
  for (i in metascore_int$Accession)
  {
    index <- metascore_int$Accession == i
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

    # start filling
    #print(bla)
    # test if something to fill
    if(sum(grepl("CP|EC",topodf$feature)) == 0)
    {
      print("No fitting topological feature found for MetaScore")
      next
    }
    metascore_int[index, "Topomatch"] <- TRUE
    # for every experiment
    myec <- topovect == "EC"
    mycp <- topovect == "CP"
    for (n in locationdf$Expname)
    {
      # EC Extracellular
      # CP Cytoplasmic
      intec <- sum(intensm[n, myec]) + 0.01
      intcp <- sum(intensm[n, mycp]) + 0.01

      metascore_int[index, n] <- intec/intcp
      # if(metascore_int[index, n] > 1000)
      # {
      #   print(bla)
      # }

    }
  }


  return(metascore_int)
}


# detect if there is a N or C Terminal Truncation
# Input: numerical vector minpepl, ac, bc, character vector acc
# Output: numerical vector ntt, nttcov, ctt, cttcov, seccov
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
  nttcov <- 0
  ctt <- 0
  cttcov <- 0
  #this is like ceiling(0.75 * numberofruns)
  #sec <- which(ac >= 0.75)
  #this is like ceiling(0.5 * numberofruns)
  sec <- which(ac >= 0.5)
  #sec <- which(ac >= 0.5)
  if (length(sec) == 0)
  {
    return(list(0,0,0,0,0,0))
  }

  #prot <- which(bc >= 0.75)
  prot <- which(bc >= 0.5)
  if (length(prot) == 0)
  {
    return(list(0,0,0,0,0,0))
  }

  #   print(paste0("sec:",sec, "\n"))
  #   print(paste0("prot:",sec, "\n"))
  #   print(acc)

  # which returns the indices of a logical object
  # gibt es Signale im Proteomdatensatz die vor dem ersten Signal im
  # Sekretomdatensatz liegen? Sind diese laenger als das kuerzeste Peptid
  # fuer diese Accession?

  #N-Terminale Trunkierung
  # welche Position in prot ist kleiner oder gleich der ersten in sec, die einen Mittelwert ueber/= 0.5 hat
  lenprot <- length(which(prot <= sec[1]))

  if (lenprot >= minpepl)
  {
    #N-Terminale Trunkierung moeglich
    cat("N-Terminale Trunkierung moeglich: ", acc, "\n")
    ntt <- 1
    nttcov <- lenprot/sec[1]
  }

  #C-Terminale Trunkierung
  clenprot <- length(which(prot >= sec[length(sec)]))

  if (clenprot >= minpepl)
  {
    #C-Terminale Trunkierung moeglich
    cat("C-Terminale Trunkierung moeglich: ", acc,"\n")
    ctt <- 1
    cttcov <- clenprot/(length(ac)-sec[length(sec)]+1)
  }

  #SecProtCov
  secprotcov <- length(sec)/length(ac)
  seccutcov <- 0

  if(secprotcov >= 0.2)
  {
    # seccutcov
    # N and C terminal truncation
    if ((ntt+ctt) == 2 )
    {
      seccutcov <- length(sec)/(sec[length(sec)]-sec[1]+1)

    } else if (ntt == 1)
      #only N terminal truncation
    {
      seccutcov <- length(sec)/(length(ac)-sec[1]+1)
    } else if (ctt == 1)
      #only C terminal truncation
    {
      seccutcov <- length(sec)/sec[length(sec)]
    } else
    {
      seccutcov <- 0
    }
  }


  # if (seccov < 0.2)
  # {
  #   seccov <- 0
  # }

  return(list(ntt, nttcov, ctt, cttcov, secprotcov, seccutcov))
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
      #print(bla)
      actProtein$Annotation <- actFASTAindex[actFASTAindex$actProteinAcc == acc, "actProtein"]
      y <- actFASTAindex[actFASTAindex$actProteinAcc == acc, "index"]

      actProtein$Name <- strsplit(strsplit(actProtein$Annotation, "|",
                                           fixed=TRUE)[[1]][[3]], " OS=", fixed=TRUE)[[1]][[1]]

      #actProtein$Sequence <- getSequence(actFASTA[[y]], as.string = TRUE)[[1]]

      # actProtein$Abbreviation <- strsplit(actProtein$Annotation, "GN=", fixed=TRUE)[[1]][2]
      # actProtein$Abbreviation <- strsplit(actProtein$Abbreviation, " ", fixed=TRUE)[[1]][1]

      actProtein$Genename <- strsplit(actProtein$Annotation, "GN=", fixed=TRUE)[[1]][2]
      actProtein$Genename <- strsplit(actProtein$Genename, " ", fixed=TRUE)[[1]][1]

      actProtein$Abbreviation <- strsplit(actProtein$Name, " ", fixed=TRUE)[[1]][1]
      actProtein$Abbreviation <- strsplit(actProtein$Abbreviation, "_", fixed=TRUE)[[1]][1]

      # actProtein$Genename <- strsplit(actProtein$Name, " ", fixed=TRUE)[[1]][1]
      # actProtein$Genename <- strsplit(actProtein$Genename, "_", fixed=TRUE)[[1]][1]

      actProtein$Name <- paste(strsplit(actProtein$Name, " ", fixed=TRUE)[[1]][-1], collapse=" ")

      actProtein$AccessionFound <- TRUE
    }


    print(c(actProtein$Annotation, actProtein$Name, actProtein$Genename, actProtein$AccessionFound))
    result[result$Accession == acc, c("Annotation","Protein_names","Gene_names","AccessionFound")] <- c(actProtein$Annotation, actProtein$Name, actProtein$Genename, actProtein$AccessionFound)
    print(result[result$Accession == acc, c("Annotation","Protein_names","Gene_names","AccessionFound")])
  }


  return(result)
}
