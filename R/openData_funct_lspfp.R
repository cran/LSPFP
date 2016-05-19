# Secretome Lysat Peptide Descriptor
#
# sourcefile functions


#--------------------------------------------------------------------------------------
# open txt file of an experiment fitted to the newest MaxQuant output
# Input: list experiment, numeric vector fileNumber indicating the fileName stored in experiment$FileName
# Output: list of data.table data, list experiment

openTXT <- function (experiment, fileNumber) 
{
  txt <- "- Starting TXT-read..."
  #printTcl(txt)
  print(txt)
  expNamesTmp <- ""
  data <- {}
  fileName <- paste0(experiment$Path, experiment$FileName[fileNumber])
  characterCols <- c("Charges", "Oxidation (M) site IDs", "MS/MS IDs", "Evidence IDs", "Best MS/MS", 
                     "Mod. peptide IDs", "Protein group IDs", "Reverse", "Potential contaminant", "id",
                     "Carbamidomethyl (C) site IDs","NEM site IDs")
  header <- readLines(fileName,n=1)
  header2 <- unlist(strsplit(header,"\t"))
  
  data <- suppressWarnings(data.table::fread(fileName, 
                            colClasses=list(character=characterCols[characterCols %in% header2]),
                            integer64 = "double"))
  
  
  # replace all spaces
  setnames(data, gsub(" ", "_", colnames(data)))
  setnames(data, "Leading_razor_protein", "Accession")
  setnames(data, "Unique_(Groups)", "Unique")
  setnames(data, gsub("Experiment", "Counts", colnames(data)))
  
  # Delete rows with Reverse or Contaminat = "+"
  #data <- data[Reverse != "+"]
  #data <- data[Potential_contaminant != "+"]  
  # try to avoid warning
  data <- data["Reverse" != "+"]
  data <- data["Potential_contaminant" != "+"]  
  
  # Include index in table
  #data[,dataIndex:=(1:nrow(data))] 
  #try to avoid warning
  data[,"dataIndex":=(1:nrow(data))] 
  
  # Dimensions of data.frame 
  dataRows <- nrow(data)
  dataCols <- ncol(data)
  
  # Change column name for the Accession numbers
  #   setnames(data, "Leading razor protein", "Accession")
  #   setnames(data, "Unique (Groups)", "Unique")
  
  # Setting Unique as logical. FALSE + 0 = 0, TRUE + 0 = 1. 
  #   set(data,which(data[,Unique]=="yes"),"Unique","TRUE")
  #   set(data,which(data[,Unique]=="no"),"Unique","FALSE")
  #try to avoid warning
  #print(bla)
  Unique <- NULL
  set(data,which(data[,Unique]=="yes"),"Unique","TRUE")
  set(data,which(data[,Unique]=="no"),"Unique","FALSE")
  data$Unique <- as.logical(data$Unique) 
  #  setnames(data, gsub(" ", "_", colnames(data)))
  #  setnames(data, gsub("Experiment", "Counts", colnames(data)))
  
  # Extracting the number of different experimental groups. 
  # The 1st data row it is checked for "Normalized abundance" and "Spectral counts". 
  # The 2nd row is checked for the number and the names of experimental groups.
  
  expColsInt <- which(substring((colnames(data)),1, 10) == "Intensity_")
  experiment$StartColIntensities[fileNumber] <- head(expColsInt, 1)
  
  expColsCnt <- which(substring((colnames(data)),1, 7) == "Counts_")
  experiment$StartColCounts[fileNumber] <- head(expColsCnt, 1)
  
  experiment$NumExps[fileNumber] <- length(expColsInt)
  if (length(length(expColsInt)) != length(length(expColsCnt))) 
  {
    txt <- "- Warning: Experiment and Intensity numbers do not match ! ! !"
    #tkinsert(txtWidget,"0.0",paste(txt, "\n"))
    print(txt)
  }
  
  experiment$ExpNames[fileNumber] <- list(substring(colnames(data)[expColsInt], 11))
  
  # suppressWarnings(
  #   for (j in names(data[,startColCounts[fileNumber]:(startColCounts[fileNumber]+(numExps[fileNumber]*2)-1), with=FALSE])){
  #     set(data,which(is.na(data[[j]])),j,0)
  #   }
  # )
  
  setnames(data, paste0(colnames(data), "_File", fileNumber))
  experiment$ColumnNamesInt[fileNumber] <- list(colnames(data)[expColsInt])
  experiment$ColumnNamesCnt[fileNumber] <- list(colnames(data)[expColsCnt])
  
  return(list(data, experiment))
}


#--------------------------------------------
# open txt file of an experiment fittet to the newest MaxQuant output
# Input: path to peptide.txt
# Output: data.table data

ss_openTXT <- function (path) 
{
  txt <- "- Starting TXT-read..."  
  print(txt)
  expNamesTmp <- ""
  data <- {}
  #fileName <- paste0(experiment$Path, experiment$FileName[fileNumber])
  #characterCols <- c("Charges", "Oxidation (M) site IDs", "MS/MS IDs", "Evidence IDs", "Best MS/MS", 
   #                  "Mod. peptide IDs", "Protein group IDs", "Reverse", "Potential contaminant", "id") 
  characterCols <- c("Charges", "Oxidation (M) site IDs", "MS/MS IDs", "Evidence IDs", "Best MS/MS", 
                      "Mod. peptide IDs", "Protein group IDs", "Reverse", "Potential contaminant", "id",
                      "Carbamidomethyl (C) site IDs","NEM site IDs")
  header <- readLines(path,n=1)
  header2 <- unlist(strsplit(header,"\t"))
  
  data <- suppressWarnings(data.table::fread(path, 
                            colClasses=list(character=characterCols[characterCols %in% header2]),
                            integer64 = "double"))
  
#   data <- data.table::fread(path, 
#                             colClasses=list(character=characterCols),
#                             integer64 = "double",verbose=TRUE)
  
  
  # replace all spaces
  setnames(data, gsub(" ", "_", colnames(data)))
  #nnames <- gsub(" ","_", colnames(data))
  #   nnames(grep("Unique_\\(Groups\\)",nnames)) <- Unique
  #   nnames(grep("Leading_razor_protein",nnames)) <- Unique
  setnames(data, "Leading_razor_protein", "Accession")
  setnames(data, "Unique_(Groups)", "Unique")
  setnames(data, gsub("Experiment", "Counts", colnames(data)))
  
  # Delete rows with Reverse or Contaminat = "+" 
#   data <- data[Reverse != "+"]
#   data <- data[Potential_contaminant != "+"]
  #try to avoid warning
  data <- data["Reverse" != "+"]
  data <- data["Potential_contaminant" != "+"]
  
  # Include index in table
  #data[,dataIndex:=(1:nrow(data))] 
  # try to avoid warnings
  data[,"dataIndex":=(1:nrow(data))] 
  
  # Dimensions of data.frame 
  dataRows <- nrow(data)
  dataCols <- ncol(data)
  
  # Change column name for the Accession numbers
  #   setnames(data, "Leading razor protein", "Accession")
  #   setnames(data, "Unique (Groups)", "Unique")
  
  # Setting Unique as logical. FALSE + 0 = 0, TRUE + 0 = 1. 
  #set(data,which(data[,Unique]=="yes"),"Unique","TRUE")
  #set(data,which(data[,Unique]=="no"),"Unique","FALSE")
  #try to avoid warnings
  #print(bla)
  Unique <- NULL
  set(data,which(data[, Unique]=="yes"),"Unique","TRUE")
  set(data,which(data[, Unique]=="no"),"Unique","FALSE")
  data$Unique <- as.logical(data$Unique)
  #setattr(data, "Unique", logical)

  #data[, "Unique":=as.logical("Unique")]
  #  setnames(data, gsub(" ", "_", colnames(data)))
  #  setnames(data, gsub("Experiment", "Counts", colnames(data)))
  
  # Extracting the number of different experimental groups. 
  # The 1st data row it is checked for "Normalized abundance" and "Spectral counts". 
  # The 2nd row is checked for the number and the names of experimental groups.
  
  expColsInt <- which(substring((colnames(data)),1, 10) == "Intensity_")
  #experiment$StartColIntensities[fileNumber] <- head(expColsInt, 1)
  
  expColsCnt <- which(substring((colnames(data)),1, 7) == "Counts_")
  #experiment$StartColCounts[fileNumber] <- head(expColsCnt, 1)
  
  #experiment$NumExps[fileNumber] <- length(expColsInt)
  if (length(length(expColsInt)) != length(length(expColsCnt))) 
  {
    txt <- "- Warning: Experiment and Intensity numbers do not match ! ! !"    
    print(txt)
  }
  
  
  return(data)
}


# open csv file of an experiment
# Input: list experiment, numeric vector fileNumber indicating the fileName stored in experiment$FileName
# Output: list of logical vector result, data.table data, list experiment

openCSV <- function (experiment, fileNumber) 
{
  data <- {}
  result <- FALSE
  
  txt <- "- Starting CSV-read..."
  txt <- paste(txt, experiment$FileName[fileNumber])  
  print(txt)

  fileName <- paste0(experiment$Path, experiment$FileName[fileNumber])
  headers <- read.table(file = fileName, check.names = FALSE, sep = ',', header = FALSE, stringsAsFactors = FALSE, nrows = 3)
  
  # warum wird hier die letzte Spalte abgeschnitten? da stehen Daten drin!
  # fuehrt weiter unten auch zu Folgefehlern
  #   headers <- headers[-ncolHeaders]
  #   ncolHeaders <- ncol(headers)
  experiment$StartColIntensities[fileNumber] <- which(headers[1,]=="Normalized abundance")
  
  
  #test if there are Raw abundance values
  until <- which(headers[1,]=="Spectral counts")-1
  raw <- 0
  if(sum(headers[1,]=="Raw abundance", na.rm = TRUE) >= 1)
  {
    raw <- which(headers[1,]=="Raw abundance")
    #test order of the columns
    if ( raw < experiment$StartColIntensities[fileNumber])
    {
      cat("\n", "Unexspected order of the columns Normalized abundance and Raw abundance!", "\n")
      cat("Normalized abundance should be directly before Raw abundance!", "\n")
      result <- FALSE
      return(list(result,data,experiment))
    } 
    else if (raw >= until )
    {
      cat("\n", "Unexspected order of the columns Raw abundance Spectral counts!", "\n")
      cat("Raw abundance should be directly before Spectral counts!", "\n")
      result <- FALSE
      return(list(result,data,experiment))
    }
    else
    {
      #delete Raw Abundance from data      
      headers[ ,(raw:until)] <- list(NULL)      
    }   
  }
  
  ncolHeaders <- ncol(headers)
  experiment$StartColCounts[fileNumber] <- which(headers[1,]=="Spectral counts")
  distCntInt <- experiment$StartColCounts[fileNumber] - experiment$StartColIntensities[fileNumber]
  sequenceCol <- which(headers[3,]=="Sequence")
  accessionCol <- which(headers[3,]=="Accession")
  numExpGroups <- length(which(headers[2,]!=""))
  numExpGroups <- numExpGroups / 2
  expGroups <- headers[2,][headers[2,]!=""][1:numExpGroups]
  
  expNamesTmp <-  headers[3,experiment$StartColIntensities[fileNumber]:
                            (experiment$StartColCounts[fileNumber]-1)]
  
  experiment$ColumnNamesInt[fileNumber] <- list(paste0("Intensity_", 
                                                       headers[3,experiment$StartColIntensities[fileNumber]:
                                                                 (experiment$StartColCounts[fileNumber]-1)]))
  
  experiment$ColumnNamesCnt[fileNumber] <- list(paste0("Counts_", 
                                                       headers[3,experiment$StartColCounts[fileNumber]:ncolHeaders]))
  
  headers[3,experiment$StartColIntensities[fileNumber]:
            (experiment$StartColCounts[fileNumber]-1)] <- 
    unlist(experiment$ColumnNamesInt[fileNumber])
  
  headers[3,experiment$StartColCounts[fileNumber]:ncolHeaders] <- 
    unlist(experiment$ColumnNamesCnt[fileNumber])
  
  headers[3,]  <- gsub(" ", "_", headers[3,])
  
  # Dimensions of data.frame 
  headersRows <- nrow(headers)
  headersCols <- ncol(headers)
  
  experiment$ExpNames[fileNumber] <- list(expNamesTmp)
  data <- read.table(file=fileName, check.names=TRUE, sep=',', 
                     header=FALSE, stringsAsFactors=FALSE, skip=3)
  
  #if there are raw abundances in the data  
  if (raw > 0)
  {
    #delete Raw Abundance from data      
    data[ ,(raw:until)] <- list(NULL)
  }
  #print(bla)
  
  colnames(data) <- headers[3,]
  
  # Dimensions of data.frame 
  dataRows <- nrow(data)
  dataCols <- ncol(data)
  
  intensities <- data[,(experiment$StartColIntensities[fileNumber]):
                        (experiment$StartColCounts[fileNumber]-1)]
  
  counts <- data[,experiment$StartColCounts[fileNumber]:dataCols-1]
  
  #hier wird auch wieder etwas abgeschnitten
  #data <- data[,1:dataCols-1]
  data <- data[,1:dataCols]
  
  data <- data.table(data)
  
  experiment$NumExps[fileNumber] <- distCntInt
  
  setnames(data, "Use_in_quantitation", "Unique")
  # Setting Unique as logical. FALSE + 0 = 0, TRUE + 0 = 1. 
  #data[, Unique:=as.logical(Unique)]
  #try to avoid no visible binding
  #print(bla)
  #setattr(data, "Unique", logical)
  #data[, "Unique":=as.logical("Unique")]
  data$Unique <- as.logical(data$Unique)
  setkey(data, "Accession", "Sequence")
  ## Get Counts and Intensities
  DTintensities <- data[, lapply(.SD, sum), by=c("Accession", "Sequence"), 
                        .SDcols=(experiment$StartColIntensities[fileNumber]):
                          (experiment$NumExps[fileNumber]+
                             experiment$StartColIntensities[fileNumber]-1)]
  
  #ACHTUNG der Bereich von .SDcols ist groesser als die data.table data
  #   Browse[1]> (experiment$StartColCounts[fileNumber]):
  #     +     (experiment$NumExps[fileNumber]+
  #              +          experiment$StartColCounts[fileNumber]-1)
  #   [1] 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50
  # Browse[1]> length(colnames(data))
  # [1] 41
  
  DTcounts <- data[, lapply(.SD, sum), by=c("Accession", "Sequence"), 
                   .SDcols=(experiment$StartColCounts[fileNumber]):
                     (experiment$NumExps[fileNumber]+
                        experiment$StartColCounts[fileNumber]-1)] 
  
  ## Get column IDs from columns with annotations
  annotationCols <- 1:(experiment$StartColIntensities[fileNumber]-1)
  ## Exclude Unique column
  annotationCols <- annotationCols[-(which(colnames(data)=="Unique"))]
  DTannotations <- data[, lapply(.SD, list), 
                        by=c("Accession", "Sequence"), .SDcols=annotationCols]
  
  DTannotations$Accession <- 
    #substring(data[, list(Accession=paste(Accession, collapse="")), by=c("Accession", "Sequence")]$Accession, 1,6)
    #try to avoid binding variable warning
    substring(data[, list("Accession"=paste("Accession", collapse="")), 
                   by=c("Accession", "Sequence")]$Accession, 1,6)
  
  DTannotations <- DTannotations[,lapply(.SD, as.character), by=c("Accession", "Sequence")]
  
  #print(bla)
  
  DTunique <- data[, lapply(.SD, all), 
                   by=c("Accession", "Sequence"), .SDcols="Unique"]
  
  data <- cbind(DTannotations, DTunique, DTintensities, DTcounts)
  ## Remove additional Sequence columns
  removeCols <- which(colnames(data)=="Sequence")
  data <- data[,!removeCols[-1] , with=FALSE]
  ## Remove additional Accession columns
  removeCols <- which(colnames(data)=="Accession")
  data <- data[,!removeCols[-1] , with=FALSE]
  
  #data[,dataIndex:=(1:nrow(data))]
  #try to avoid warning
  data[,"dataIndex":=(1:nrow(data))]
  setnames(data, paste0(colnames(data), "_File", fileNumber))
  experiment$ColumnNamesInt[fileNumber] <- 
    list(paste0(experiment$ColumnNamesInt[[fileNumber]],  "_File", fileNumber))
  
  experiment$ColumnNamesCnt[fileNumber] <- 
    list(paste0(experiment$ColumnNamesCnt[[fileNumber]],  "_File", fileNumber))
  
  result <- TRUE
  return(list(result, data, experiment))
}


# open csv file of an experiment for the secscoreWrapper
# Input: character vector with the path of the csv-file (Progenesis)
# Output:  data.table data

ss_openCSV <- function (path) 
{
  data <- {}
  result <- FALSE
  
  txt <- "- Starting CSV-read..." 
  print(txt)
  
  headers <- read.table(file = path, check.names = FALSE, sep = ',', header = FALSE, stringsAsFactors = FALSE, nrows = 3)
  # warum wird hier die letzte Spalte abgeschnitten? da stehen Daten drin!
  # fuehrt weiter unten auch zu Folgefehlern
  #   headers <- headers[-ncolHeaders]
  #   ncolHeaders <- ncol(headers)
  
  #remove columns taht contain NA
  naheader <- grep("NA",headers[1,])
  headers[ , naheader] <- NULL
  
  StartColIntensities <- which(headers[1,] == "Normalized abundance")
  
  #test if there are Raw abundance values
  until <- which(headers[1,]=="Spectral counts")-1
  raw <- 0
  if(sum(headers[1,]=="Raw abundance",na.rm = TRUE) >= 1)
  {
    raw <- which(headers[1,]=="Raw abundance")
    #test order of the columns
    if ( raw < StartColIntensities)
    {
      cat("\n", "Unexspected order of the columns Normalized abundance and Raw abundance!", "\n")
      cat("Normalized abundance should be directly before Raw abundance!", "\n")
      result <- FALSE
      return(result)
    } 
    else if (raw >= until )
    {
      cat("\n", "Unexspected order of the columns Raw abundance Spectral counts!", "\n")
      cat("Raw abundance should be directly before Spectral counts!", "\n")
      result <- FALSE
      return(result)
    }
    else
    {
      #delete Raw Abundance from data      
      headers[ ,(raw:until)] <- list(NULL)      
    }   
  }
  
  ncolHeaders <- ncol(headers)
  StartColCounts <- which(headers[1,]=="Spectral counts")
  distCntInt <- StartColCounts - StartColIntensities
  sequenceCol <- which(headers[3,]=="Sequence")
  accessionCol <- which(headers[3,]=="Accession")
  numExpGroups <- length(which(headers[2,]!=""))
  numExpGroups <- numExpGroups / 2
  expGroups <- headers[2,][headers[2,]!=""][1:numExpGroups]
  
  expNamesTmp <-  headers[3,StartColIntensities:
                            (StartColCounts-1)]
  
  ColumnNamesInt <- list(paste0("Intensity_", headers[3, StartColIntensities : (StartColCounts-1)]))
  
  ColumnNamesCnt <- list(paste0("Counts_", headers[3, StartColCounts:ncolHeaders]))
  
  headers[3, StartColIntensities : (StartColCounts-1)] <- unlist(ColumnNamesInt)
  
  headers[3,StartColCounts : ncolHeaders] <- unlist(ColumnNamesCnt)
  
  headers[3,]  <- gsub(" ", "_", headers[3,])
  
  # Dimensions of data.frame 
  headersRows <- nrow(headers)
  headersCols <- ncol(headers)
  
  ExpNames <- list(expNamesTmp)
  
  data <- read.table(file=path, check.names=TRUE, sep=',', 
                     header=FALSE, stringsAsFactors=FALSE, skip=3)
  #print(bla)
  
  #remove NA columns (see above)
  data[ ,naheader] <- NULL
  
  #if there are raw abundances in the data  
  if (raw > 0)
  {
    #delete Raw Abundance from data      
    data[ ,(raw:until)] <- list(NULL)
  }
  
  colnames(data) <- headers[3,]
  
  # Dimensions of data.frame 
  dataRows <- nrow(data)
  dataCols <- ncol(data)
  
  intensities <- data[ ,(StartColIntensities) : (StartColCounts-1)]
  
  counts <- data[ ,StartColCounts:dataCols-1]
  
  #hier wird auch wieder etwas abgeschnitten
  #data <- data[,1:dataCols-1]
  data <- data[,1:dataCols]
  
  colnames(data)[colnames(data) == "Use_in_quantitation"] <- "Unique"
  data$Unique <- as.logical(data$Unique)
  
  
  
  return(data)
}