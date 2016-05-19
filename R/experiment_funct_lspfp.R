# Secretome Lysat Peptide Descriptor
# 
# Experiment functions
# here you will find all functions, that
# are involved with the experiment creation 
# and administration

#-----------------------------------------------------------------------
#Function to create the experiment basic file- and data-structure
# Input:
# char vector expname, 
# char vector sourcefiles (max two sourcefiles)
#                         one for Proteom and one for Secretome
# Output: list of logical vector result and list experiment
# list experiment
# result == FALSE if something went wrong 

createNewExp <- function(expname, sourcefiles) 
{
  experiment <- {} 
  result <- TRUE
  #print("start test")
  # Test input
  if (grepl("[^A-Za-z0-9_]", expname) == TRUE )
  { 
    print("Experiment name is invalid!")
    print("Allowed pattern for the experimentname: [A-Za-z0-9_]")
    result <- FALSE
    return(list(result,experiment))
  }  
  
  if(length(sourcefiles) == 0)
  {
    print("No sourcefilepath is given")
    result <- FALSE
    return(list(result,experiment))
  }  
  
  fesource <- file.exists(sourcefiles)
  
  if(sum(fesource) < length(sourcefiles))
  {
    print("Following sourcefilepath is wrong")
    print(sourcefiles[fesource %in% FALSE])
    result <- FALSE
    return(list(result,experiment))
  }
  
  path <- paste0("AnalysisData/",expname)
  
  if(dir.exists(path))
  {
    print("Experiment already exists!")
    print(path)
    print("Please choose an other name for your new experiment,")
    print("or delete the existing one!")
    result <- FALSE
    return(list(result,experiment))
  }
  
  #start creating experiment  
  print("start creating experiment folder")
  
  dir.create(path)
  experiment$Name <- expname  
  experiment$Path <- paste0(path,"/")
  experiment$InfoFile <- paste0(experiment$Path, "expinfo.txt")
  
  tmp <- getExperimentSourceFiles(experiment, sourcefiles)
  
  if (tmp[[1]] == FALSE)
  {
    print("createNewExp: Something went wrong with the sourcefiles")
    print("createNewExp: No experiment was created")
    result <- FALSE
    file.remove(path, recursive = TRUE)
    
    return(list(result,experiment))
  }
  
  experiment <- tmp[[2]]
  save(experiment, file=experiment$InfoFile, ascii=TRUE)      
  
  
  return(list(result,experiment))
}


## Get name and path of experiment data files ##
# Input: list experiment, char vector sourcefiles (max two sourcefiles)
#                         one for Proteom and one for Secretome
# Output: list of logical vector result and list experiment
# 

getExperimentSourceFiles <- function(experiment, sourcefiles)
{
  result <- TRUE
  
  #Filters <- matrix(c("CSV or Text", "*.csv;*.txt"), 1, 2, byrow = TRUE)  
  #experiment$Path <- paste0("data/", experiment$Name, "/")
  #experiment$InfoFile <- paste0(experiment$Path, "expinfo.txt")
  
  # choose.file is only available on Windows, Value: character vector giving zero or more file paths
  #   experiment$SourceFile <- choose.files(default=getwd(), caption = "Select files", 
  #                                         multi = TRUE, filters = Filters, index = nrow(Filters))
  
  #experiment$SourceFile <- sourcefiles
  
  #if (length(experiment$SourceFile) > 2 | length(experiment$SourceFile) < 1) {
  lnsrc <- length(sourcefiles)
  if (lnsrc > 2 | lnsrc < 1)
  {
    #     ReturnVal <- tkmessageBox(title = "Error in file selection",
    #                               message = "Please select one or two files.", icon = "info", type = "ok")
    print("getExpSF: Error in file selection, Please select one or two files.")
    result <- FALSE
    return(list(result, experiment))
  }
  
  experiment$SourceFile <- sourcefiles
  
  # get the file extension.
  experiment$FileName <- basename(experiment$SourceFile)
  experiment$NumberOfFiles <- length(experiment$FileName)
  
  experiment$FileType <- sub('.*(?=.{4}$)', '', experiment$FileName, perl=TRUE)
  
  if (length(experiment$FileName) == 2 & 
        experiment$FileType[1] != experiment$FileType[2]) {
    #     ReturnVal <- tkmessageBox(title = "Error in datatypes",
    #                               message = "Files do not have the same file extension.", icon = "info", type = "ok")
    print("getExpSF: Error in datatypes, Files do not have the same file extension.")
    result <- FALSE
    return(list(result, experiment))
  }
  
  if (!(experiment$FileType[[1]] == ".csv" | experiment$FileType[[1]] == ".txt"))
  {
    stop("getExpSF: Data file in wrong format")  
    result <- FALSE
    return(list(result, experiment))
  }
  #####################################
  ## Copy data files into /data/ dir ##
  #####################################
  #wozu ist das gut?
  test <- file.copy(experiment$SourceFile, paste0(experiment$Path, experiment$FileName))
  
  #assign("experiment", experiment, envir=globalenv())
  
  return(list(result, experiment))
}


## Load existing experiment desrciption from expinfo.txt in the experiment folder
# Input: char vector expName
# Output: list of logical vector result (something went wrong = FALSE), list experiment
LoadExp <- function(expName){
  
  result <- TRUE
  experiment <-{}
  
  if(length(expName) > 1 )
  {
    print("Just one experiment name allowed")
    result <- FALSE
    return(list(result,experiment))
  }
  
  path <- paste0("AnalysisData/", expName, "/", "expinfo.txt")
  
  if (file.exists(path))
  {
    # loaded file was named experiment before   
    # compare createNewExp()
    #experiment <- {}
    load(file=path)
    #result <- experiment
  }
  else
  {
    print("Could not find Info-File!")
    result <- FALSE
    #return(list(result,experiment))
  }
  
  
  return(list(result,experiment))
}


## Load data into "data" data.table##
#Auch hier wird experiment wieder uebergeben
# Input: list experiment
# Output: list of logical vector result, data.table data and list experiment

loadData <- function (experiment)
{
  result <- FALSE
  tmpData <- {}
  data <- {}
  
  #fileNumber point to the index of the FileName in experiment
  # you now know the sourcefile, just coding  the sourcefile
  for (fileNumber in 1:length(experiment$FileName))
  {
    if (experiment$FileType[[1]] == ".csv")
    {
      mytmp <- openCSV(experiment, fileNumber)
      if(mytmp[[1]] == FALSE) 
      {
        cat("\n", "Error during openCSV!","\n")
        result <- FALSE
        return(list(result, data, experiment))
      }
      #print(bla)
      listDataExperiment <- list(mytmp [[2]], mytmp[[3]])
      
    }
    else
    {
      listDataExperiment <- openTXT(experiment, fileNumber)
    }
    
    data <- listDataExperiment[[1]]
    experiment <- listDataExperiment[[2]]
    tmpData[fileNumber] <- list(data)
    colnames(tmpData[[fileNumber]])
    colNameAccession <- paste0("Accession_File", fileNumber)
    experiment$AccessionNumbers[fileNumber] <- unique(tmpData[[fileNumber]][,colNameAccession, with=FALSE])
    colNameSequence <- paste0("Sequence_File", fileNumber)
    experiment$Sequences[fileNumber] <- unique(tmpData[[fileNumber]][,colNameSequence, with=FALSE])
    setnames(tmpData[[fileNumber]], colNameSequence, "Sequence")
    tmpData[[fileNumber]]$Accession <- tmpData[[fileNumber]][, colNameAccession, with=FALSE]    
  }
  
  #add tables together from to or more files
  tmpDataAll <- as.data.table(tmpData[[1]])
  setkey(tmpDataAll,  "Accession", "Sequence")
  
  if (experiment$NumberOfFiles > 1)
  {
    for (fileNumber in 2:experiment$NumberOfFiles)
    {
      setkey(tmpData[[fileNumber]], "Accession", "Sequence")
      tmpDataCur <- as.data.table(tmpData[[fileNumber]])
      tmpDataAll <- merge(tmpDataAll, tmpDataCur, all=TRUE)
    }
    
    rm(tmpDataCur)
  }
  
  data <- tmpDataAll
  rm(tmpDataAll)
  nrowData <- nrow(data)
  data$DataIndexAll <- (1:nrowData)
  
  #############################################
  ## Combine Accession columns to one column ##
  #############################################
  accessionCols <- which(substring(colnames(data), 1, 9)=="Accession")
  
  if (experiment$NumberOfFiles > 1)
  {
    for (fileNumber in 2:experiment$NumberOfFiles)
    {      
      data[is.na(data$Accession), "Accession"] <- 
        data[is.na(data$Accession),accessionCols[fileNumber], with=FALSE]
    }
  }
  
  ##########################################
  ## Combine Unique columns to one column ##
  ##########################################
  uniqueCols <- which(substring(colnames(data), 1, 11)=="Unique_File")
  
  if (length(uniqueCols) != experiment$NumberOfFiles)
  {
    print("Number of files not equal number of Unique columns")
  }
  
  uniques <- (data[,uniqueCols, with=FALSE])
  uniques[is.na(uniques)] <- TRUE
  data$Unique <- Reduce("&", uniques)
  
  ## experiment$AccessionStack holds the accessions of all the proteins which should be process in the given order
  ## For now, if two files are loaded, this is the union of the two accession lists.
  ## Sorted: intersect, Diff1, Diff2
  #   experiment <- combiningData(experiment)
  experiment$AccessionStack <- unique(data[, "Accession", with=FALSE])
  
  result<- TRUE
  return(list(result, data, experiment))
}


#was macht die Funktion genau?
# Input: list experiment, char vector output mit Outputtyp
# Output:
openProgress <- function(experiment, output = "PDF"){ 
  
  factorExists <- experiment$Factor
  factorLevels <- experiment$Levels
  experiment$OutputType <- output
  
  for (fileNumber in 1:experiment$NumberOfFiles){
    accessions <- data[, paste0("Accession_File", fileNumber), with = FALSE]
    ## If NOT Unique in one file, the peptide will be accounted for as NOT unique in all files
    ## Thus, Unique column is used and not the file specific column
    accessions <- cbind(accessions, data[, "Unique", with = FALSE])
    nrowNonUnique <- length(unique(na.exclude(accessions$Accession)))
    nrowUnique <- length(unique(na.exclude(accessions[accessions$Unique == TRUE]$Accession)))
    
    txt <- paste0("Dataset ", fileNumber, ": ", nrowNonUnique, " (", nrowUnique, ")")
    
    print(txt)
  }
  
}