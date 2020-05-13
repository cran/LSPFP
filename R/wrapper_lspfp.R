# Secretome and Lysat Peptide Feature Plotter
#
# Wrapper functions
#
# encapsulates workflows

#---------------------------------------------------------------------------------
# function wrapperLSPFP wrapps the whole workflow from
# downloading the Uniprot files, processing them, doing the filtering and scoring
# and finaly the plots
#
# Input:
# PLM part
# character vector globpath with the path to the global directory of BasicData and AnalysisData
# character vector expname with the experiment name
# character vector sourcefiles with one or two sourcefiles pathes
# character vector version indicates what version out of the BasicData file should be used
# character vector species giving the uniprot species names that should be downloaded to BasicData
# character vector proteomeid giving the uniprot proteome id of the organisem, must have the same order as species
# character vector taxid containing the uniprot taxonomic IDs, must have the same order as species
# character vector domain conatining the uniprot domain description, must have the same order as species
# logical vector forcedl indicates if the download should run again
#
# SecScore part
# PLM part
# logical vector TRUE if everything is ok
# list experiment
# logical vector unipep indicates if only unique peptides will be used, FALSE if not

#
# Output:
# PLM part


wrapperLSPFP <- function(globpath, expname, sourcefiles, org, grlocationdf,
                        version = "actual",
                        species = c("HUMAN", "MOUSE", "RAT", "PIG"),
                        proteomeid = c("UP000005640","UP000000589","UP000002494","UP000008227"),
                        taxid = c("9606","10090","10116","9823"),
                        domain = rep("Eukaryota",4),
                        forcedl = FALSE,
                        pepstack = 2, pepque = 2, sortprint = "fcsmall",
                        unipep = TRUE, localfasta = "none")
{
  #test input
  #test for missing values
  if (missing(globpath))
  {
    cat("wrapperLSPFP: ","No global path was specified!","\n")
    return(FALSE)
  }
  if(class(globpath)!= "character")
  {
    cat("wrapperLSPFP: ","globpath should be character!","\n")
    return(FALSE)
  }

  if (missing(expname))
  {
    cat("wrapperLSPFP: ","No experiment name was specified!","\n")
    return(FALSE)
  }
  if(class(expname)!= "character")
  {
    cat("wrapperLSPFP: ","expname should be character!","\n")
    return(FALSE)
  }

  if (missing(sourcefiles))
  {
    cat("wrapperLSPFP: ","No sourcefile was specified!","\n")
    return(FALSE)
  }
  if(class(sourcefiles)!= "character")
  {
    cat("wrapperLSPFP: ","sourcefiles should be character!","\n")
    return(FALSE)
  }

  if (missing(org))
  {
    cat("wrapperLSPFP: ","No organisem was specified!","\n")
    return(FALSE)
  }
  if(class(org)!= "character")
  {
    cat("wrapperLSPFP: ","org should be character!","\n")
    return(FALSE)
  }

  if (missing(grlocationdf))
  {
    cat("wrapperLSPFP: ","No grlocationdf was specified!","\n")
    return(FALSE)
  }
  if(class(grlocationdf)!= "data.frame")
  {
    cat("wrapperLSPFP: ","grlocationdf should be a data.frame!","\n")
    return(FALSE)
  }
  if(class(grlocationdf$Expname)!= "character")
  {
    cat("wrapperLSPFP: ","grlocationdf$expname should be character","\n")
    return(FALSE)
  }
  if(class(grlocationdf$Location)!= "character")
  {
    cat("wrapperLSPFP: ","grlocationdf$location should be character","\n")
    return(FALSE)
  }
  if(sum(grepl("Secretome|Proteome",grlocationdf$Location)) != length(grlocationdf$Location))
  {
    cat("wrapperLSPFP: ","grlocationdf$location should be Secretome or Proteome!","\n")
    return(FALSE)
  }
  if(class(grlocationdf$Treatment)!= "character")
  {
    cat("wrapperLSPFP: ","grlocationdf$Treatment should be character","\n")
    return(FALSE)
  }
  if(sum(grepl("[A-Z][A-Z]",grlocationdf$Treatment)) != length(grlocationdf$Treatment))
  {
    cat("wrapperLSPFP: ","grlocationdf$Treatment should be [A-Z][A-Z]!","\n")
    return(FALSE)
  }
  if(!is.numeric(grlocationdf$Sample))
  {
    cat("wrapperLSPFP: ","grlocationdf$Sample should be numeric","\n")
    return(FALSE)
  }
  if(!is.numeric(grlocationdf$Group))
  {
    cat("wrapperLSPFP: ","grlocationdf$Group should be numeric","\n")
    return(FALSE)
  }

  if(class(version)!= "character")
  {
    cat("wrapperLSPFP: ","version should be character!","\n")
    return(FALSE)
  }
  if(class(species)!= "character")
  {
    cat("wrapperLSPFP: ","species should be character!","\n")
    return(FALSE)
  }

  if(class(proteomeid)!= "character")
  {
    cat("wrapperLSPFP: ","proteomeid should be character!","\n")
    return(FALSE)
  }

  if(class(taxid)!= "character")
  {
    cat("wrapperLSPFP: ","taxid should be character!","\n")
    return(FALSE)
  }

  if(class(domain)!= "character")
  {
    cat("wrapperLSPFP: ","domain should be character!","\n")
    return(FALSE)
  }

  if(class(forcedl)!= "logical")
  {
    cat("wrapperLSPFP: ","forcedl should be logical!","\n")
    return(FALSE)
  }

  if(class(sortprint)!= "character")
  {
    cat("wrapperLSPFP: ","sortprint should be character!","\n")
    return(FALSE)
  }
  if( !grepl("^fcsmall$|^fclf$|^trunc$|^acc$|^ntts$|^ctts$|^trunc2$",sortprint))
  {
    print("wrapperLSPFP: ERROR sortprint is not correct")
    return(FALSE)
  }

  if(class(unipep)!= "logical")
  {
    cat("wrapperLSPFP: ","unipep should be logical!","\n")
    return(FALSE)
  }


  if(!is.numeric(pepstack))
  {
    cat("wrapperLSPFP: ","pepstack should be numeric!","\n")
    return(FALSE)
  }

  if(!is.numeric(pepque))
  {
    cat("wrapperLSPFP: ","pepque should be numeric!","\n")
    return(FALSE)
  }

  if (class(localfasta) != "character"  )
  {
    cat("wrapperLSPFP: ","localfasta should be character!","\n")
    return(FALSE)
  }


  #print(bla)
  print("wrapperLSPFP: Start wrapperPLM")
  test <- wrapperPLM(globpath, expname, sourcefiles, org,
             version, species, proteomeid, taxid,
             domain, forcedl, localfasta)
  #return(test)

  if (test[[1]] == FALSE)
  {
    print("wrapperLSPFP: wrapperPLM stopped with an error")
    return(FALSE)
  }
  print("wrapperLSPFP: End wrapperPLM")


  print("wrapperLSPFP: Start wrapperSecScore")
  exp <- test[[2]]


  #replace spaces in expname
  grlocationdf$Expname <- gsub(" ","_",grlocationdf$Expname)
  test2 <- wrapperSecScore(org, exp, grlocationdf, globpath, sortprint = sortprint, unipep = unipep )

  if(test2 == FALSE)
  {
    print("wrapperLSPFP: wrapperSecScore stopped with an error")
    return(FALSE)
  }

  print("wrapperLSPFP: End wrapperSecScore")
  return(TRUE)
}


#---------------------------------------------------------------------------------
# function wrapper PLM wrapps the rest of the functionality of PLM
# Input and Output see above


wrapperPLM <- function(globpath, expname, sourcefiles, org,
                       version = "actual",
                       species = c("HUMAN", "MOUSE", "RAT", "PIG"),
                       proteomeid = c("UP000005640","UP000000589","UP000002494","UP000008227"),
                       taxid = c("9606","10090","10116","9823"),
                       domain = rep("Eukaryota",4),
                       forcedl = FALSE, localfasta = "none")
{
  #test for missing values
  if (missing(globpath))
  {
    cat("wrapperPLM: ","No global path was specified!","\n")
    return(FALSE)
  }

  if (missing(expname))
  {
    cat("wrapperPLM: ","No experiment name was specified!","\n")
    return(FALSE)
  }
  if (missing(sourcefiles))
  {
    cat("wrapperPLM: ","No sourcefile was specified!","\n")
    return(FALSE)
  }

  cat("\n","Start: Creating global directory structure.","\n")
  #create main directory with subfolder BasicData and AnalysisData
  setDirStruc(globpath)


  #print(getwd())
  #change wd globpath
  userwd <- getwd()
  setwd(globpath)
  #print(getwd())


  #load needed data sets from UniProt or uses archived ones (local)
  #Datensaetze von UniProt herunter laden oder zuvor gespeicherte verwenden
  # ACHTUNG: es koennen keine archivierten Datensaetze herunter geladen werden,
  #         sollte also eine aeltere Version eine Org gebraucht werden,
  #         muss dieser manuell eingefuegt werden.

  #Test if the user wants to use his own FASTA file
  if (localfasta == "none")
  {
    cat("wrapperPLM:","Start: Prepare fasta and gff files from UniProt.","\n")
    if (version=="actual")
    {
      t <- fasta2rds(species,proteomeid,taxid,domain,forcedl)
      t2 <- gff2rds(species, forcedl)

    }else if (dir.exists(paste0(globpath,"/BasicData/", version))) {

      t <- paste0(globpath,"/BasicData/", version)
      t2 <- paste0(globpath,"/BasicData/", version)

    }else {

      cat("wrapperPLM STOP: your BasicData version does not exist!","\n")
      #change wd back to user
      setwd(userwd)
      return(FALSE)
    }

  # user wants his own FASTA file
  }else {
    cat("wrapperPLM ","Start: Prepare local FASTA file.","\n")
#     if(!version == "actual" )
#     {
#       cat("warpperPLM STOP: just the actual version is allowed if you want to use your own FASTA file!","\n")
#       #change wd back to user
#       setwd(userwd)
#       return(FALSE)
#
#     } else

# test the used variables
    if ( !(length(species) = 1))
    {

      cat("wrapperPLM STOP: just one species is allowed if you want to use your own FASTA file!","\n")
      #change wd back to user
      setwd(userwd)
      return(FALSE)

    } else if (!file.exists(localfasta)){

      cat("wrapperPLM STOP: your localfasta FASTA file does not exist!","\n")
      #change wd back to user
      setwd(userwd)
      return(FALSE)

    } else {
      #now decide which version the user wants
      #if (version == "actual" )

        t <- localfasta2rds(localfasta, version, org)
        if(t == "ERROR")
        {
          cat("wrapperPLM STOP: something went wrong during localfasta2rds!","\n")
          #change wd back to user
          setwd(userwd)
          return(FALSE)
        }
        #cat("t:",t,"\n")
        t2 <- localgff2rds(t, version, org)
        #cat("t2:",t2,"\n")
        if(t2 == "ERROR")
        {
          cat("wrapperPLM STOP: something went wrong during localgff2rds!","\n")
          #change wd back to user
          setwd(userwd)
          return(FALSE)
        }
#         t <- paste0("BasicData/",t)
#         t2 <- paste0("BasicData/",t2)

    }

  }




  cat("\n","Start: Create new experiment.","\n")
  #create new experiment structure /Folder Sources and expinfo.txt
  tmp <- createNewExp(expname, sourcefiles)
  if (tmp[[1]] == FALSE)
  {
    cat("STOP: something went wrong during createNewExp","\n")
    #change wd back to user
    setwd(userwd)
    return(FALSE)
  }


  cat("\n","Load Existing experiment.","\n")
  #load existing expinfo.txt
  tmp2 <- LoadExp(expname)
  if (tmp2[[1]] == FALSE )
  {
    cat("STOP 2: something went wrong during LoadExp","\n")
    #change wd back to user
    setwd(userwd)
    return(FALSE)
  }
  experiment <- tmp2[[2]]

  cat("\n","Start: Load data from fasta and gff files","\n")
  #load Data form gff and fasta?
  listDataExperiment <- loadData(experiment)
  if (listDataExperiment[[1]] == FALSE )
  {
    cat("STOP 2b: something went wrong during loadData","\n")
    #change wd back to user
    setwd(userwd)
    return(FALSE)
  }

  data <- listDataExperiment[[2]]
  experiment <- listDataExperiment[[3]]


  #Reste von userInterfaceSettings
  # Irgendwie unnoetig
  # openProgress(experiment, "PDF")
  experiment$OutputType <- "PDF"
  experiment$RunFrom <- 1
  experiment$RunTo <- length(experiment$AccessionStack)


  # setze Pfade zu den Datensaetzen von UniProt
  experiment$UniprotrelFASTA <- t
  experiment$UniprotrelGFF <- t2

  # save the lastest version of experimentinfo to the experiment folder
  save(experiment, file=experiment$InfoFile, ascii=TRUE)


  #table output

  # proteinsAnalysed data.table
  # actFASTAindex data.table
  # actFASTA list
  # actGFF data.table
  # seqData list
  # data ?
  #result <- list(TRUE, res_plot[[2]], res_sD[[3]], res_sD[[4]], res_sD[[5]], res_sD[[6]], data )
  result <- list(TRUE,experiment)

  #change wd back to user
  setwd(userwd)
  #print(getwd())

  return(result)
}


#--------------------------------------------------------------------
# wrapper for the secretome score functions

wrapperSecScore <- function(org, exp, grlocationdf, globpath, pepstack = 2, pepque = 2, sortprint = "fcsmall", unipep = TRUE)
{
  #print(bla)
  #test for missing values
  if (missing(org))
  {
    cat("\n","No organisem was specified!","\n")
    return(FALSE)
  }

  #change wd globpath
  userwd <- getwd()
  setwd(globpath)

  #-------------
  #load BasicData
  print("wrapperSecScore: load data from files")
  # list of SeqFastaAA objects
  #
  # > str(myFastaHuman[[1]])
  # Class 'SeqFastaAA'  atomic [1:1] MGAPLLSPGWGAGAAGRRWWMLLAPLLPALLLVRPAGALVEGLYCGTRDCYEVLGVSRSAGKAEIARAYRQLARRYHPDRYRPQPGDEGPGRTPQSAEEAFLLVATAYETLKVSQAAAELQQYCMQNACKDALLVGVPAGSNPFREPRSCALL
  # ..- attr(*, "name")= chr "tr|A0A024R161|A0A024R161_HUMAN"
  # ..- attr(*, "Annot")= chr ">tr|A0A024R161|A0A024R161_HUMAN Guanine nucleotide-binding protein subunit gamma OS=Homo sapiens GN=DNAJC25-GNG10 PE=3 SV=1"
  #
  Fasta <- readRDS( paste0(exp$UniprotrelFASTA, "/FASTA", tolower(org), ".rds"))

  #Uniprot GFF annotations from BasicData folder
  #
  #> names(myGffHuman)
  # [1] "seqname"    "source"     "feature"    "start"      "end"        "score"      "strand"     "frame"      "attributes"
  # [10] "Note"       "ID"
  #
  Gff <- as.data.frame(readRDS(paste0(exp$UniprotrelGFF,"/GFF",tolower(org),".rds")))

  #Indexfile: indicates at whitch position an Accession is stored in the list FASTAhuman.rds
  #
  # > names(myACCHuman)
  # [1] "index"         "actProteinAcc" "actProtein"
  #
  Acc <- as.data.frame(readRDS(paste0(exp$UniprotrelFASTA, "/", toupper(org), ".ACC.RDS")))


  #-------
  # load experiment data

  if (exp$FileType[[1]] == ".csv")
  {
    # load Progenesis Peptides
    peptides <- as.data.frame(ss_openCSV(paste0(exp$Path,"/",exp$FileName)))

    #fuege hier nun start und stopposition der Daten hinzu
    peptides <- getPositions(Acc, Fasta, peptides)

  } else {
    # load MaxQuant Peptides
    peptides <- as.data.frame(ss_openTXT(paste0(exp$Path,"/",exp$FileName)))
  }

  #filter out non unique peptides if wanted
  if (unipep)
  {
    peptides <- peptides[peptides$Unique == TRUE, ]
  }

  #if peptides are empty
  if(length(peptides[,1]) < 1)
  {
    print("wrapperSecScore_ERROR: peptides are empty")
    return(FALSE)
  }


  #test if peptides do belong to the given org
 # print("wrapperSecScore: test organism")
 # acclen <- nchar(peptides$Accession)
 # accindex <- which(acclen == 6)
 # url <- paste0("https://www.uniprot.org/uniprot/?query=",peptides[accindex[1],"Accession"],"&columns=organism&format=tab")

 # bdown(url,paste0(exp$Path,"/orgtest"))
 # orgtest <- read.delim(paste0(exp$Path,"/orgtest"))
 # morg <- grepl(org, orgtest[,1], ignore.case = TRUE)
 # if(sum(morg) < 1)
  if(5 < 1)
  {
    print("wrapperSecScore_ERROR: org does not match with organism of accession in peptides-file")
    return(FALSE)
  }


  #print(bla)
  #peptide_filter <- function(peptides, pepstack = 2, pepque = 2)
  print("wrapperSecScore: start filtering")
  fpeptides <- peptide_filter(peptides, grlocationdf, pepstack, pepque)
  fpeptides <- fpeptides[fpeptides$delete_prot == FALSE, ]

  print("wrapperSecScore: start combine peptides")
  #   #log scale auf Intensitaeten
  #   t1 <- fpeptides[ , grep("Intensity_",names(fpeptides))]
  #   # add 1 to all intensity values to avoid that values between [0,1] are going to be negative log(1) = 0 log(0)= -Inf
  #   t1 <- t1+1
  #   t1 <- ceiling(log2(t1))
  #
  #   #set counts that are NA to 0
  #   t2 <- fpeptides[ , grep("Counts_", names(fpeptides))]
  #   t2[is.na(t2)] <- 0
  intenscountdf <- data.frame(fpeptides[ , c("Accession", "Sequence","Start_position","End_position","filtered_out")],
                              fpeptides[ , grep("Counts_", names(fpeptides))],
                              fpeptides[ , grep("Intensity_",names(fpeptides))] , stringsAsFactors = FALSE)

  #in case there is more than one peptide at the same position for the same protein
  # filtered_out is TRUE if the peptide is filtered out
  intenscountdf_test <- intenscountdf
  #key <- paste0(intenscountdf[ ,  c("Accession", "Sequence","Start_position","End_position")], collapse = "")
  key <- rep("nix",length(intenscountdf[ , 1]))
  for(kt in seq_along(key))
  {
    key[kt] <- paste0(intenscountdf[kt,  c("Accession", "Sequence","Start_position","End_position")], collapse = "_")
  }

  index <- seq_along(intenscountdf[ , 1])
  keyindex <- data.frame(key, index, stringsAsFactors = FALSE)
  #tmpIntens <- intenscountdf[1:length(index), ]
  tmpIntens <- intenscountdf[1:length(unique(key)), ]
  posIntens <- grep("Intensity_",names(tmpIntens))
  posCounts <- grep("Counts_",names(tmpIntens))
  i <- 1
  for(k in unique(keyindex$key))
  {
    t <- intenscountdf[keyindex[keyindex$key == k, "index"], ]
    # if there is just one peptide
    if(length(t[ ,1]) == 1)
    {
      tmpIntens[i, ] <- t

    } else if(length(t[,1]) > 1) {
      tmpIntens[i, c("Accession","Sequence","Start_position","End_position")] <- t[1,
                          c("Accession","Sequence","Start_position","End_position")]
      tmpIntens[i, posIntens] <- ceiling(colSums(t[ , posIntens], na.rm = TRUE)/length(t[,1]))
      tmpIntens[i, posCounts] <- ceiling(colSums(t[ , posCounts], na.rm = TRUE)/length(t[,1]))
      tmpIntens[i, "filtered_out"] <- (sum(t[ , "filtered_out"], na.rm = TRUE) >= 1)

    } else {
      next
    }
    i <- i+1
  }
  #print(bla)
  intenscountdf <- tmpIntens

  print("wrapperSecScore: start preparing Intensities and Counts")
  #log scale auf Intensitaeten
  t1 <- grep("Intensity_",names(intenscountdf))
  # add 1 to all intensity values to avoid that values between [0,1] are going to be negative log(1) = 0 log(0)= -Inf
  intenscountdf[ , t1] <- intenscountdf[ , t1]+1
  # new version, to avoid loss of detail
  intenscountdf[ , t1] <- log2(intenscountdf[ , t1])
  # old version with ceiling
  # intenscountdf[ , t1] <- ceiling(log2(intenscountdf[ , t1]))

  #set counts that are NA to 0
  t2 <- grep("Counts_", names(intenscountdf))
  for (ts in t2)
  {
    intenscountdf[is.na(intenscountdf[ , ts]), ts] <- 0
  }


  #just the names for annotation
  # erstelle namesdf ueber die annotations, da im Progenesis diese nicht vorkommen
  #namesdf <- unique(fpeptides[ , c("Accession","Gene_names","Protein_names")])
  acc_flt <- unique(fpeptides$Accession)
  namesdf <- lspfp_getFastaInformation(Acc, acc_flt, Fasta)

  #print(bla)

  print("wrapperSecScore: calculate features")
  myfeatures <- createFeatureTable(acc_flt, intenscountdf, Fasta, Acc, Gff, grlocationdf)

  #----------------------------------------------------------------------------------------
  mymeta <- addMetaScoreFeature(acc_flt, intenscountdf, Fasta, Acc, Gff, grlocationdf)
  write.csv(mymeta[[1]], file=paste0(exp$Path,"/","metaLength.csv"), row.names = FALSE)
  write.csv(mymeta[[2]], file=paste0(exp$Path,"/","metaIntens.csv"), row.names = FALSE)
  #print(bla)

  # mymeta[mymeta$Topomatch == TRUE,][1,]

  #ttest
  # var(mymeta[mymeta$Topomatch == TRUE,][1,3:7])
  # sone <-mymeta[mymeta$Topomatch == TRUE,]
  # sone <- as.vector(t(sone[2,3:7]))
  # stwo <-mymeta[mymeta$Topomatch == TRUE,]
  # stwo <- as.vector(t(stwo[2,8:12]))
  # myt <- t.test(sone,stwo)
  #
  #sam
  # samdata <- mymeta[mymeta$Topomatch == TRUE, 3:12]
  # gruppe <- c(rep("proteome",5), rep("secretome",5))
  # library(siggenes)
  # samOutCat <- sam(samdata, gruppe)
  # samdatares <- cbind(Accession = mymeta[mymeta$Topomatch == TRUE, "Accession"], samdata, samOutCat@d, samOutCat@p.value, stringsAsFactors = FALSE)
  #
  # samMeta <- merge(samdatares, mymeta, by = "Accession", all.y = TRUE)
  # myfeatures <- merge(myfeatures, samMeta, by = "Accession", all.x = TRUE)
  #
  # mymeta[mymeta$Topomatch == TRUE,]
  # summary(mymeta[mymeta$Topomatch == TRUE,])
  #------------------------------------------------------------------------------------

  #return(list(ts, ntts, ctts))
  tmptscore <- truncscoring(myfeatures)
  tscore <- tmptscore[[1]]
  ntts <- tmptscore[[2]]
  ctts <- tmptscore[[3]]
  tscore2 <- tmptscore[[4]]

  tmpsscore <- secscoring(myfeatures)
  fcsmall <- tmpsscore[[1]]
  fclf <- tmpsscore[[2]]

  myfeatures2 <- cbind(myfeatures, ntts, ctts, tscore, tscore2, fcsmall, fclf, stringsAsFactors =  FALSE)

  myfeatures2<- merge(myfeatures2, namesdf[, c("Accession","Gene_names","Protein_names")], by = "Accession", all.x = TRUE)

  # plot the log or log lf sorting
  if(sortprint == "fcsmall")
  {
    sorteddf <- myfeatures2[order(myfeatures2[ ,"fcsmall"], decreasing = TRUE ), ]
    printp <- paste0(exp$Path,"/","peptide_plot_sorted_fcsmall")

  } else if(sortprint == "fclf") {
    sorteddf <- myfeatures2[order(myfeatures2[ ,"fclf"], decreasing = TRUE ), ]
    printp <- paste0(exp$Path,"/","peptide_plot_sorted_fclf")

  } else if(sortprint == "acc") {
    sorteddf <- myfeatures2[order(myfeatures2[ ,"Accession"] ), ]
    printp <- paste0(exp$Path,"/","peptide_plot_sorted_accession")

  } else if(sortprint == "trunc") {
    sorteddf <- myfeatures2[order(myfeatures2[ ,"tscore"], decreasing = TRUE ), ]
    printp <- paste0(exp$Path,"/","peptide_plot_sorted_trunc")

  } else if(sortprint == "trunc2") {
    sorteddf <- myfeatures2[order(myfeatures2[ ,"tscore2"], decreasing = TRUE ), ]
    printp <- paste0(exp$Path,"/","peptide_plot_sorted_trunc2")

  } else if(sortprint == "ntts") {
    sorteddf <- myfeatures2[order(myfeatures2[ ,"ntts"], decreasing = TRUE ), ]
    printp <- paste0(exp$Path,"/","peptide_plot_sorted_ntts")

  } else if(sortprint == "ctts") {
    sorteddf <- myfeatures2[order(myfeatures2[ ,"ctts"], decreasing = TRUE ), ]
    printp <- paste0(exp$Path,"/","peptide_plot_sorted_ctts")

  } else {
    print("wrapperSecScore_ERROR: sortprint not given correct")
    return(FALSE)
  }


  print("wrapperSecScore: save data")
  saveRDS(sorteddf, paste0(exp$Path,"/","feature_table.rds"))
  write.csv(sorteddf, file=paste0(exp$Path,"/","feature_table.csv"), row.names = FALSE)
  write.csv(intenscountdf, file=paste0(exp$Path,"/","intenscount_table.csv"), row.names = FALSE)
  write.csv(grlocationdf, file=paste0(exp$Path,"/","grlocationdf_table.csv"), row.names = FALSE)
  write.csv(namesdf, file=paste0(exp$Path,"/","namesdf_table.csv"), row.names = FALSE)
  # save the lastest version of experimentinfo to the experiment folder
  exp$org <- org
  experiment <- exp
  save(experiment, file=exp$InfoFile, ascii=TRUE)

  print("wrapperSecScore: plot peptides to pdf")

  print_peptides(sorteddf, intenscountdf, Fasta, Acc, Gff, grlocationdf, printp, namesdf)
  dev.off()

  #set back wd
  setwd(userwd)

  return(TRUE)
}
