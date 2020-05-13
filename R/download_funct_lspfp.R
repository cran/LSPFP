# Secretome Lysat Peptide Descriptor
#
# Collection of Uniprot download functions

# download a proeteom fasta file from Uniprot-FTP-Server and save it to BasicData
# This function first checks the actual release version of Uniprot FTP Server files.
# If there is a new version, a new Release Folder will be created in BasicData.
# Than the FASTA files of the required Proteom and its additions are downloaded.
# Both are unzipped and loaded again with read.fasta from "seqinr".
# These results are stored as rds files in the BasicData folder


# Input: Domain (vector),Proteom_ID (vector),Tax_ID (vector),speciesname (vector),forcedl (logical)
# Output: char vector result contains path to directory with FASTA-rds-files
#Proteom_ID,Tax_ID,Species
#UP000005640 9606    HUMAN
#UP000000589 10090   MOUSE
#UP000002494 10116   RAT
#UP000008227 9823    PIG
#Domain: Archaea,Bacteria,Eukaryota,Viruses#
#forcedl: if old databasic/existing files should be replaced

#fasta2rds <- function(species=c("HUMAN", "MOUSE", "RAT", "PIG")){
fasta2rds <- function(species=c("HUMAN", "MOUSE", "RAT", "PIG"),
                      proteomeid=c("UP000005640","UP000000589","UP000002494","UP000008227"),
                      taxid=c("9606","10090","10116","9823"),
                      domain=rep("Eukaryota",4), forcedl = FALSE)
{
  #wd <- paste0("bin/")
  wd <- c("BasicData/")
  FASTAEndingGZ <- ".fasta.gz"
  FASTAEnding <- ".fasta"
  #FASTAPath <- "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/proteomes/"
  FASTAPath <- "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/"


  #erstelle eine aktuelle version der Datenbank
  #ACHTUNG: hier auch ueberpruefen ob Datei schon vorhanden
  #save Uniprot Readme
  #check version/ release of Uniprot Data on uniprot and in BasicData folder
  cat("Check version \n")
  bdown("ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/relnotes.txt",
        paste0(wd, "relnotes.txt"))
  version <- gsub(" ", "", readLines(paste0(wd, "relnotes.txt"), 1))
  cat("\n")

  if (!dir.exists(paste0(wd, version)))
  {
    dir.create(paste0(wd, version))
  }

  if (!file.exists(paste0(wd, version, "/relnotes.txt")) | forcedl)
  {
    file.copy(paste0(wd, "relnotes.txt"), paste0(wd, version, "/relnotes.txt"), overwrite = FALSE)
  }

  file.remove(paste0(wd, "relnotes.txt"))
  result <- paste0(wd, version)

  cat("Starting FASTA file download, extraction and RDS creation...\n")
  uniprot <- data.frame(species, proteomeid, taxid, domain, stringsAsFactors = FALSE)


  for (i in seq_along(uniprot$species))
  {
    s <- uniprot[i,"species"]
    pid <- uniprot[i,"proteomeid"]
    tid <- uniprot[i,"taxid"]
    d <- uniprot[i,"domain"]
    ################################################
    ## Download FASTA files
    ################################################
    FASTAsource <- paste0(FASTAPath, d, "/", pid, "_", tid, FASTAEndingGZ)
    FASTAdest <- paste0(wd, version,"/", s,"_",pid, "_", tid, FASTAEndingGZ)
    FASTAfile <- paste0(wd, version,"/", s,"_",pid, "_", tid, FASTAEnding)
    FASTAsourceAdd <- paste0(FASTAPath, d, "/", pid, "_", tid, "_additional", FASTAEndingGZ)
    FASTAdestAdd <- paste0(wd, version,"/", s,"_",pid, "_", tid, "_additional", FASTAEndingGZ)
    FASTAfileAdd <- paste0(wd, version,"/", s,"_",pid, "_", tid, "_additional", FASTAEnding)

    FASTArds <- paste0(wd, version, "/", "FASTA", tolower(s), ".rds")

    #print(bla)

    cat(paste("Start download check:", s, "\n"))

    if(!file.exists(FASTAfile) | forcedl )
    {
      ret <- bdown(FASTAsource, FASTAdest)
      cat("\n")
      cat("Extracting:", s, "\n")
      gunzip(FASTAdest, remove=TRUE, overwrite=TRUE)
    }
    else
    {
      cat(FASTAfile, " is up to date. \n")
    }

    if(!file.exists(FASTAfileAdd) | forcedl )
    {
      retadd <- bdown(FASTAsourceAdd, FASTAdestAdd)
      cat("\n")
      cat("Extracting:", s, "\n")
      gunzip(FASTAdestAdd, remove=TRUE, overwrite=TRUE)
    }
    else
    {
      cat(FASTAfileAdd, " is up to date. \n")
    }

    if (!file.exists(FASTArds) | forcedl )
    {
      cat("Reading FASTA file for", s, "\n")

      ###and save as RDS files ##
      FASTAopen <- read.fasta(FASTAfile, as.string=TRUE, seqtype="AA", )
      FASTAopenAdd <- read.fasta(FASTAfileAdd, as.string=TRUE, seqtype="AA", )
      #zusammenfuegen
      FASTAopenTg <- c(FASTAopen,FASTAopenAdd)
      cat("Saving RDS file for", s, "\n")
      saveRDS(FASTAopenTg, file=FASTArds)
    }
    else
    {
      cat(FASTArds, " is up to date. \n")
    }


    ######################################################
    ## Create and save index, accessions, protein names ##
    ######################################################
    if (!file.exists(paste0(wd,version,"/", s,".ACC.RDS"))| forcedl)
    {
      actProteinAcc <- {}
      actProtein <- {}
      actFASTAindex <- {}
      actFASTA <- FASTAopenTg
      index <- seq(1:length(actFASTA))
      lenFASTA <- length(actFASTA)
      cat("Creating Accession table for ", s, "   #", length(actFASTA), "\n")
      barWidth <- (getOption("width")-20)
      pb <- txtProgressBar(style=3, , width=barWidth)
      m <- 0
      actProtein <- character(lenFASTA)
      actProteinAcc <- character(lenFASTA)
      for (y in 1:lenFASTA){
        actProtein[y] <- getAnnot(actFASTA[[y]], as.string = TRUE)
        actProtein[y] <- as.character(actProtein[y])
        actProteinAcc[y] <- strsplit(actProtein[y], "|", fixed=TRUE)[[1]][2]
        progress <- (y / lenFASTA)
        setTxtProgressBar(pb, progress)
      }
      cat("\nSaving result...")
      actFASTAindex <- data.table(index, actProteinAcc, actProtein)
      saveRDS(actFASTAindex, paste0(wd,version,"/", s,".ACC.RDS"))
      close(pb)
      cat(paste0("Finished Accession table for ", s, "\n"))
      cat("______________________________________________\n\n")
    }
    else
    {
      cat(paste0(wd,version,"/", s,".ACC.RDS"), " is up to date. \n")
    }
  }
  cat("Finished process for\n")
  cat(species, "\n")
  cat("______________________________________________\n\n")

  return(result)
}


#localfasta2rds allows to use a local fasta file.
# Input: character vector pathlf (path to local fasta file)
#         character vector version
# Output: char vector result contains path to directory with FASTA-rds-files
localfasta2rds <- function(pathlf, version, org)
{
  if ( !file.exists(pathlf) )
  {
    cat("localfasta2rds ERROR: There is no file called:\n",pathlf,"\n")
    result <- "ERROR"
    return(result)
  }

  pathsplit <- unlist(strsplit(pathlf,"/"))
  fastaname <- pathsplit[length(pathsplit)]

  if(version=="actual")
  {
    cat("Check version \n")
    bdown("ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/relnotes.txt",
          paste0("BasicData/", "relnotes.txt"))
    actversion <- gsub(" ", "", readLines(paste0("BasicData/", "relnotes.txt"), 1))
    cat("\n")
    newdir <- paste0("BasicData/", actversion, "_", fastaname, "_", org)

    if (dir.exists(newdir))
    {
      cat("localfasta2rds ERROR: there is already a directory:\n", newdir,"\n")
      result <- "ERROR"
      return(result)
    }

    dir.create(newdir)
    file.copy(pathlf, paste0(newdir, "/"))


  } else if (dir.exists(paste0("BasicData/", version))) {

    newdir <- paste0("BasicData/",version,"_",fastaname,"_",org)
    if (dir.exists(newdir))
    {
      cat("localfasta2rds ERROR: there is already a directory:\n", newdir,"\n")
      result <- "ERROR"
      return(result)
    }

    dir.create(newdir)
    file.copy(pathlf, paste0( newdir, "/"))


    #version was not correct
  } else {
    cat("localfasta2rds ERROR: no directory matches to the given version \n")
    result <- "ERROR"
    return(result)
  }

  # save fasta in RDS
  FASTAopen <- read.fasta(pathlf, as.string = TRUE, seqtype = "AA", )
  cat("localfasta2rds: Saving FASTA-RDS file for",fastaname , "\n")
  FASTArds <- paste0(newdir, "/FASTA", tolower(org), ".rds")
  saveRDS(FASTAopen, file = FASTArds)

  #create index
  if (!file.exists(paste0(newdir, "/", toupper(org),".ACC.RDS")))
  {
    actProteinAcc <- {}
    actProtein <- {}
    actFASTAindex <- {}
    actFASTA <- FASTAopen
    index <- seq(1:length(actFASTA))
    lenFASTA <- length(actFASTA)
    cat("Creating Accession table for ", fastaname, "   #", length(actFASTA), "\n")
    barWidth <- (getOption("width")-20)
    pb <- txtProgressBar(style=3, , width=barWidth)
    m <- 0
    actProtein <- character(lenFASTA)
    actProteinAcc <- character(lenFASTA)
    for (y in 1:lenFASTA){
      actProtein[y] <- getAnnot(actFASTA[[y]], as.string = TRUE)
      actProtein[y] <- as.character(actProtein[y])
      actProteinAcc[y] <- strsplit(actProtein[y], "|", fixed=TRUE)[[1]][2]
      progress <- (y / lenFASTA)
      setTxtProgressBar(pb, progress)
    }
    cat("\nSaving result...")
    actFASTAindex <- data.table(index, actProteinAcc, actProtein)
    saveRDS(actFASTAindex, paste0(newdir,"/", toupper(org),".ACC.RDS"))
    close(pb)
    cat(paste0("Finished Accession table for ", fastaname, "\n"))
    cat("______________________________________________\n\n")
  }
  else
  {
    cat(paste0(newdir,"/", toupper(org),".ACC.RDS"), " is up to date. \n")
  }


  result <- newdir
  return(result)
}

#########################################
## Download function with progress bar ##
#########################################
bdown <- function(url, file)
{

  h <- basicTextGatherer()
  curlPerform(url=url, customrequest='HEAD', header=0L, nobody=1L, headerfunction=h$update)

  if(grepl('Content-Length: ', h$value()))
  {
    #print(paste0("bdown Ergebnis Abfrgae Contenlength", grepl("Content-Length:",h$value())))
    size <- as.numeric(strsplit(strsplit(h$value(),'Content-Length: ')[[1]][2],"\r")[[1]][1])
  }
  else
  {
    size <- 1
  }

  #size <-2000
  f <- CFILE(file, mode="wb")
  #print(paste0("bdown size: ", size))
  if (size > 1000){
    barWidth <- (getOption("width")-20)
    bar <- txtProgressBar(0, size, style=3, , width=barWidth)
    a <- curlPerform(url=url, noprogress=0L,
                     writedata = f@ref,
                     progressfunction=function(down,up)
                       setTxtProgressBar(bar, down[2]))
  }else{
    a <- curlPerform(url=url, noprogress=0L,
                     writedata = f@ref,
                     progressfunction=function(down,up)
                       cat("Downloading ", format(down[2], big.mark=",")," bytes\r"))
  }

  close(f)

  return(a)
}


#-----------------------------------------------------------------------------------------------------------
# gff download
#------------------

## Saving imported GFF-files to RDS
# this function downloads species specific gff files from uniprot and saves them as rds in
# the BasicData folder for this release.
# Input: char vector speciesName, logical vector forcedl
#       forcedl indicates if the files should be downloaded even if there is already a file with this release
# Output: char vector result contains path to directory with gff-rds-files

gff2rds <- function(speciesName = c("HUMAN", "MOUSE", "RAT", "PIG"), forcedl)
{



  GFFEnding <- ".gff"
  wd <- "BasicData/"
  GFFEndingGZ <- ".gff.gz"
  GFFPath <- "https://www.uniprot.org/uniprot/?format=gff&compress=yes&query=organism:"


  #erstelle eine aktuelle Version der Datenbank
  #ACHTUNG: hier auch ueberpruefen ob Datei schon vorhanden
  #save Uniprot Readme
  # Idee speichere die gffs passend zu den Aktualisierungen der Proteome

  # check version/ release of Uniprot Data on uniprot and in BasicData folder
  cat("Check version \n")
  bdown("ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/relnotes.txt",
        paste0(wd, "relnotes.txt"))
  version <- gsub(" ", "", readLines(paste0(wd, "relnotes.txt"), 1))
  cat("\n")

  if (!dir.exists(paste0(wd, version)))
  {
    dir.create(paste0(wd, version))
  }

  if (!file.exists(paste0(wd, version, "/relnotes.txt")) | forcedl)
  {
    file.copy(paste0(wd, "relnotes.txt"), paste0(wd, version, "/relnotes.txt"), overwrite = FALSE)
  }

  file.remove(paste0(wd, "relnotes.txt"))
  result <- paste0(wd, version)

  cat("Starting GFF file download, extraction and RDS creation...\n")
  for (s in speciesName){
    ##############################################
    ## Download GFF files and save as RDS files ##
    ##############################################
    GFFsource <- paste0(GFFPath, s)
    GFFdest <- paste0(wd, version, "/", s, GFFEndingGZ)
    GFFfile <- paste0(wd, version, "/", s, GFFEnding)
    GFFrds <- paste0(wd, version, "/GFF", tolower(s), ".rds")

    if(!file.exists(GFFfile) | forcedl )
    {
      cat("Downloading:", s, "\n")
      ret <- bdown(GFFsource, GFFdest)
      cat("\n")
      cat("- Finished downloading ", s, "\n")
      cat("Extracting:", s, "\n")
      gunzip(GFFdest, remove=TRUE, overwrite=TRUE)

    }
    else
    {
      cat(GFFfile, " already exists. \n")
    }


    if(!file.exists(GFFrds) | forcedl )
    {
      cat("- Using GFF-file for species: ", s, "\n")
      gfffile <- paste0(wd, version, "/", s,".gff")
      gff <- gffRead(gfffile)
      #names(gff)
      #[1] "seqname"    "source"     "feature"    "start"      "end"        "score"      "strand"     "frame"      "attributes" NA


      cat("- Extracting attributes from GFF-file (this will take a while):\n")
      cat("   Notes...\n")
      gff$Note <- getAttributeField(gff$attributes, "Note")
      cat("   IDs...\n")
      gff$ID <- getAttributeField(gff$attributes, "ID")
      gff <- (as.data.table(gff))
      cat("- GFF-parsing finished", "\n")
      cat("- Saving GFF data table as RDS file for", s, "\n")
      #GFFrds <- paste0(wd, version, "/GFF", tolower(s), ".rds")
      saveRDS(gff, file=GFFrds)
      cat("- Finished GFF.RDS generation for ", s, "\n")
      cat("___________________________________________\n")

    }
    else
    {
      cat(GFFrds, " already exists. \n")
    }


  }
  cat("- Finished GFF.RDS generation for ", "\n")
  cat(speciesName, "\n")
  cat("___________________________________________\n")


  return(result)
}

#-----------------------------------------------------------
## Saving imported GFF-files to RDS
# this function downloads species specific gff files from uniprot and saves them as rds in
# the BasicData folder for this release.
# Input: char vector pathdir
#         char vector version
#         char vector org
#
# Output: char vector result contains path to directory with gff-rds-files

localgff2rds <- function(pathdir, version, org)
{
  # check uniprot version/ release of Uniprot Data on uniprot and in BasicData folder
  cat("localgff2rds: Check version \n")
  bdown("ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/relnotes.txt",
        paste0("BasicData/", "relnotes.txt"))
  actversion <- gsub(" ", "", readLines(paste0("BasicData/", "relnotes.txt"), 1))
  cat("\n")


  if (version == "actual")
  {
    # there is already a directory
    if (dir.exists(paste0("BasicData/", actversion)))
    {
      file.copy(paste0("BasicData/", actversion,"/GFF",tolower(org),".rds"), paste0(pathdir,"/"))
      file.copy(paste0("BasicData/", actversion,"/",toupper(org),".gff"), paste0(pathdir,"/"))

      #there is no directory but there should be a newly made one from localfasta2rds
    } else if (dir.exists(pathdir)){

          # run download
          ##############################################
          ## Download GFF files and save as RDS files ##
          ##############################################
          GFFEnding <- ".gff"
          GFFEndingGZ <- ".gff.gz"
          GFFPath <- "http://www.uniprot.org/uniprot/?format=gff&compress=yes&query=organism:"

          GFFsource <- paste0(GFFPath, toupper(org))
          GFFdest <- paste0(pathdir, "/", org, GFFEndingGZ)
          GFFfile <- paste0(pathdir, "/", toupper(org), GFFEnding)
          GFFrds <- paste0(pathdir, "/GFF", tolower(org), ".rds")

          if(!file.exists(GFFfile) )
          {
            cat("localgff2rds: Downloading:", org, "\n")
            ret <- bdown(GFFsource, GFFdest)
            cat("\n")
            cat("localgff2rds: Finished downloading ", org, "\n")
            cat("localgff2rds: Extracting:", org, "\n")
            gunzip(GFFdest, remove = TRUE, overwrite = TRUE)

          }
          else
          {
            cat("localgff2rds: ",GFFfile, " already exists. \n")
          }


          if(!file.exists(GFFrds))
          {
            cat("localgff2rds: Using GFF-file for species: ", org, "\n")
            #gfffile <- paste0(wd, version, "/", s,".gff")
            #gff <- gffRead(gfffile)
            gff <- gffRead(GFFfile)
            #names(gff)
            #[1] "seqname"    "source"     "feature"    "start"      "end"        "score"      "strand"     "frame"      "attributes" NA


            cat("- Extracting attributes from GFF-file (this will take a while):\n")
            cat("   Notes...\n")
            gff$Note <- getAttributeField(gff$attributes, "Note")
            cat("   IDs...\n")
            gff$ID <- getAttributeField(gff$attributes, "ID")
            gff <- (as.data.table(gff))
            cat("- GFF-parsing finished", "\n")
            cat("- Saving GFF data table as RDS file for", org, "\n")
            #GFFrds <- paste0(wd, version, "/GFF", tolower(s), ".rds")
            saveRDS(gff, file=GFFrds)
            cat("- Finished GFF.RDS generation for ", org, "\n")
            cat("___________________________________________\n")

          }
          else
          {
            cat(GFFrds, " already exists. \n")
          }


    } else {
      #none of the given directories exist
      cat("localgff2rds ERROR: there is no pathdir with the given name\n")
      return("ERROR")
    }

    # there is an archive file for this version
  } else if (dir.exists(paste0("BasicData/", version))){

    file.copy(paste0("BasicData/", version,"/GFF",tolower(org),".rds"), paste0(pathdir,"/"))
    file.copy(paste0("BasicData/", version,"/",toupper(org),".gff"), paste0(pathdir,"/"))


  } else {

    cat("localgff2rds ERROR: there is no archive with the given name\n")
    return("ERROR")

  }



  cat("- Finished GFF.RDS generation for ", org, "\n")

  result <- pathdir
  return(result)

}



############################################################
# get values from a choosen attribute from gff

getAttributeField <- function (x, field, attrsep = ";")
{
  s = strsplit(x, split = attrsep, fixed = TRUE)
  sapply(s, function(atts)
  {
    a = strsplit(atts, split = "=", fixed = TRUE)
    m = match(field, sapply(a, "[", 1))
    if (!is.na(m))
    {
      rv = a[[m]][2]
    }
    else
    {
      rv = as.character(NA)
    }
    return(rv)
  }
  )
}

###############################################
gffRead <- function(gffFile, nrows = -1)
{
  cat("Reading ", gffFile, ": ", sep="")
  gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
                   header=FALSE, comment.char="#", nrows = nrows,
                   colClasses=c("character", "character", "character", "integer",
                                "integer",
                                "character", "character", "character", "character"))
  # es entsteht eine leere Spalte am Ende der Tabelle, da im gff ein \t am Ende steht
  # falls diese vorhanden ist loesche diese

  if (sum(gff[, length(gff)] %in% "") == length(gff[, length(gff)]))
  {
    gff[, length(gff)] <- NULL
  }

  colnames(gff) = c("seqname", "source", "feature", "start", "end",
                    "score", "strand", "frame", "attributes")
  cat("found", nrow(gff), "rows with classes:",
      paste(sapply(gff, class), collapse=", "), "\n")
  stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))

  return(gff)
}
