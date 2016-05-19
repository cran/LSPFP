# Secretome Lysat Peptide Descriptor
#
# printing functions

#printSelectedPeptides will be used to print the features after all of them are calculated
# it should be possible to choose between different sorting types and the amount of Data that will be shown
# Input: character vector path to the data folder in the analysis directory
#       character vector fname is the name of the new file
#       data.frame mysorteddf can be a how ever sortedet and shortened data.frame of the feature table
# Output: pdf_plot-file in the folder of choise

printSelectedPeptides <- function(path, fname, mysorteddf = NULL)
{
  result <- FALSE  
  if (missing(path))
  {
    cat("printSelectedPeptides ERROR: ","No path was specified!","\n")
    return(FALSE)
  }
  if(class(path)!= "character")
  {
    cat("printSelectedPeptides ERROR: ","path should be character!","\n")
    return(FALSE)
  }
  
  if (missing(fname))
  {
    cat("printSelectedPeptides ERROR: ","No fname was specified!","\n")
    return(FALSE)
  }
  if(class(fname)!= "character")
  {
    cat("printSelectedPeptides ERROR: ","fname should be character!","\n")
    return(FALSE)
  }
  
  
  #test if all needed files can be found
  print("printSelectedPeptides: check needed data")
  if(!dir.exists(path))
  {
    cat("printSelectedPeptides ERROR: there is no directory at path: ","\n",path,"\n")
    return(result)
  }   
  
  
  if(!file.exists(paste0(path,"/expinfo.txt")))
  {
    cat("printSelectedPeptides ERROR: there is no expinfo.txt in the given folder.","\n")
    return(result)
  }
  
  tmp <- unlist(strsplit(path,"/"))
  expname <- tmp[length(tmp)]
  globalpath <- paste0(tmp[1:(length(tmp)-2)], collapse = "/")
  experiment <- {}
  #org <- {}
  load(file = paste0(path,"/expinfo.txt"))    
  print("printSelectedPeptides: load data from files")
  org <- experiment$org
  
  #print(bla)
  
  fastapath <- paste0(globalpath,"/" ,experiment$UniprotrelFASTA, "/FASTA", tolower(org), ".rds")
  if(!file.exists(fastapath))
  {
    cat("printSelectedPeptides ERROR: there is no matching fasta file in the BasicData directory.","\n")
    return(result)
  }  
  # list of SeqFastaAA objects
  #
  # > str(myFastaHuman[[1]])
  # Class 'SeqFastaAA'  atomic [1:1] MGAPLLSPGWGAGAAGRRWWMLLAPLLPALLLVRPAGALVEGLYCGTRDCYEVLGVSRSAGKAEIARAYRQLARRYHPDRYRPQPGDEGPGRTPQSAEEAFLLVATAYETLKVSQAAAELQQYCMQNACKDALLVGVPAGSNPFREPRSCALL
  # ..- attr(*, "name")= chr "tr|A0A024R161|A0A024R161_HUMAN"
  # ..- attr(*, "Annot")= chr ">tr|A0A024R161|A0A024R161_HUMAN Guanine nucleotide-binding protein subunit gamma OS=Homo sapiens GN=DNAJC25-GNG10 PE=3 SV=1"
  #
  Fasta <- readRDS(fastapath)
  
  
  gffpath <- paste0(globalpath,"/", experiment$UniprotrelGFF,"/GFF",tolower(org),".rds")
  if(!file.exists(gffpath))
  {
    cat("printSelectedPeptides ERROR: there is no matching gff file in the BasicData directory.","\n")
    return(result)
  }  
  #Uniprot GFF annotations from BasicData folder
  #
  #> names(myGffHuman)
  # [1] "seqname"    "source"     "feature"    "start"      "end"        "score"      "strand"     "frame"      "attributes"
  # [10] "Note"       "ID"
  #
  Gff <- as.data.frame(readRDS(gffpath))
  
  
  accpath <- paste0(globalpath, "/", experiment$UniprotrelFASTA, "/", toupper(org), ".ACC.RDS")
  if(!file.exists(accpath))
  {
    cat("printSelectedPeptides ERROR: there is no matching acc file in the BasicData directory.","\n")
    return(result)
  }  
  #Indexfile: indicates at which position an Accession is stored in the list FASTAhuman.rds
  #
  # > names(myACCHuman)
  # [1] "index"         "actProteinAcc" "actProtein" 
  #
  Acc <- as.data.frame(readRDS(accpath))
  
  
  sortpath <- paste0(path, "/","feature_table.csv")
  if(!file.exists(sortpath))
  {
    cat("printSelectedPeptides ERROR: there is no feature_table.csv in the given folder.","\n")
    return(result)
  }
  sorteddf <- read.csv(file = sortpath, stringsAsFactors = FALSE)
  
  if (!is.null(mysorteddf))
  {
    if(sum(names(sorteddf) %in% names(mysorteddf)) == length(names(sorteddf)))
    {
      sorteddf <- mysorteddf  
    } else {
      cat("printSelectedPeptides ERROR: a column is missing in mysorteddf.","\n",
          "please support the same columns as in the feature_table.csv","\n")
      return(result)
      
    }
    
  }
  
  icpath <- paste0(path ,"/","intenscount_table.csv")
  if(!file.exists(icpath))
  {
    cat("printSelectedPeptides ERROR: there is no intenscount_table.csv in the given folder.","\n")
    return(result)
  }   
  intenscountdf <- read.csv(file = icpath , stringsAsFactors = FALSE)
  
  
  grlocpath <- paste0(path,"/","grlocationdf_table.csv")
  if(!file.exists(grlocpath))
  {
    cat("printSelectedPeptides ERROR: there is no grlocationdf_table.csv in the given folder.","\n")
    return(result)
  }   
  grlocationdf <- read.csv(file = grlocpath, stringsAsFactors = FALSE)
  
  
  namepath <- paste0(path,"/","namesdf_table.csv")
  if(!file.exists(namepath))
  {
    cat("printSelectedPeptides ERROR: there is no namesdf_table.csv in the given folder.","\n")
    return(result)
  }   
  namesdf <- read.csv(file = namepath, stringsAsFactors = FALSE)
  
  
  #print(bla)
  
  #prepare data for printing
  
  
  #run print
  printpath <- paste0(path,"/",fname)
  if(file.exists(paste0(printpath,".pdf")))
  {
    cat("printSelectedPeptides ERROR: a pdf with the given name:",fname ,"\n","already exists","\n")
    return(result)
  }   
  
  print_peptides(sorteddf,intenscountdf,Fasta,Acc,Gff,grlocationdf,printpath,namesdf)
  
  dev.off()
  
  #finish  
  result <- TRUE
  return(result)
}



# plots a given sorted data.frame in peptide line style
# Input:
# Output: pdf file
print_peptides <- function(sorteddf,intenscountdf,FastaOrg,AccOrg,GffOrg,locationdf,printp,namesdf)
{
  intenspos <- grep("Intensity_", names(intenscountdf))
  tpos <- intenscountdf[ , intenspos] <= 0
  intenscountdf[,intenspos][tpos] <- 0
  #prepare PDF as dev for plotoutputs
  generalpardefault <- par(no.readonly = TRUE)
  #cat("print_peptides generaldefault1:",generalpardefault[["new"]],"\n")
  #cat("print_peptides parnew1:",par("new"),"\n")
  
  pdfFile <- paste0(printp,".pdf")  
  pdf(paper = "a4r", file = pdfFile, width = 12, height = 10) 
  pdfpardefault <- par(no.readonly = TRUE)
  #cat("print_peptides parnew2:",par("new"),"\n")
  
  for (acc in sorteddf$Accession)
  {
    if (sum(AccOrg$actProteinAcc == acc)== 0)
    {
      protlength <- 0 
      cat("No protein in OrgACC for this Accession: ", acc,"\n")
      next
    }
    else
    {
      protlength <- getLength(FastaOrg[[AccOrg[AccOrg$actProteinAcc == acc, "index"]]])  
      protseq <- getSequence(FastaOrg[[AccOrg[AccOrg$actProteinAcc == acc, "index"]]], as.string = TRUE)
    }
        
    #-----------------------------------------------------------------------
    # calculate total number of lanes for the actual accession
    nrlanes <- 1
    for (run in locationdf$Expname)
    {
      
      #lane container
      lanes <- c(-1)
      
      stend <- intenscountdf[ (intenscountdf$Accession == acc) & (intenscountdf[ ,paste0("Intensity_",run)] >= 0)
                              & (intenscountdf[ ,paste0("Counts_",run)] > 0), 
                              c("Start_position","End_position",paste0("Intensity_",run))]
      
      stend <- stend[ order(stend$Start_position, stend$End_position), ]
      
      # for each peptide
      for (p in seq_along(stend[ , 1]))
      {
        x <- stend[p , c("Start_position","End_position")]
        
        # calculate best nr of lanes axe
        lanepos <- which(lanes < (stend$Start_position[p]-1))
        if(length(lanepos)== 0)
        {
          lanes <- c(lanes, stend$End_position[p])
          #cat("ACC: ", acc," Run: ", run, " Added one field to lane.\n")
          next
        }
        lanes[lanepos[1]] <- stend$End_position[p]
        
      }
      
      if (nrlanes < length(lanes))
      {
        nrlanes <- length(lanes)
      }
    }
    
    #cat("Accession: ", acc, " nrlanes: ", nrlanes,"\n")
    nrruns <- length(locationdf[,1])
    lyaxe <- (nrruns*nrlanes*2)+1
    
    #-------------------------------------------------------
    #calculate own colortable
    #grep("Intensity_", names(intenscountdf))
    valtab <- as.matrix(intenscountdf[intenscountdf$Accession == acc ,grep("Intensity_", names(intenscountdf))])
    valtab[valtab == -Inf] <- 0
    valtab <- table(valtab)
    
    lbiggernull <- sum(as.numeric(names(valtab)) > 0)
    
    # if all values are zero
    if (lbiggernull == 0 )
    {
      #mywhite <- rgb(255,255,255,0, maxColorValue = 255)
      #coldf <- data.frame(Intval = 0, Color = mywhite, stringsAsFactors = FALSE)
      mygrey <- rgb(204, 204, 204, alpha = 250, maxColorValue = 255)
      coldf <- data.frame(Intval = 0, Color = mygrey, stringsAsFactors = FALSE)
    }
    else
    {
      #nrcol is the stepwidth for the colours
      nrcol<-floor(255/lbiggernull)
      
      red<-seq(0, 255, nrcol)
      len<-length(red)
      #green<-seq(0,255,nrcol)
      #green<-c(seq(0,255,2*nrcol),seq(255,0,-2*nrcol))[1:len]
      green <- rep(0, nrcol)
      blue<- seq(255, 0, -nrcol)
      alpha <- rep(250, len)
      mycol<-rgb(red, green, blue, alpha, maxColorValue = 255)
      #mygrey <- rgb(204, 204, 204, alpha = 250, maxColorValue = 255)
      #intenscol <- c(mygrey, mycol[1:lbiggernull])
      intenscol <- mycol[1:lbiggernull]
      indexcol <- as.numeric(names(valtab))
      indexcol <- indexcol[indexcol > 0]
      #indexcol <- c(0,indexcol)
      indexcol <- indexcol
      coldf <- data.frame(Intval = indexcol, Color = intenscol, stringsAsFactors = FALSE)
      # barplot(1:length(intenscol),col=intenscol)
      #actcol[biggernull+1] <- mycol[1:lbiggernull]   
    }
    
    
    #plot basic
    par(oma=c(1,1,1,4))
    # GScoreT  GScoreS  tscore	sscore	sandtscore
    pos <- sorteddf$Accession == acc
    sub <- paste0("FCLF: ",  round(sorteddf$fclf[pos],8),
                  " FCsmall: ", round(sorteddf$fcsmall[pos],8),
                  " truncated: ",round(sorteddf$tscore[pos],8))
   
    nameindex <- namesdf$Accession %in% acc
    
    # ein name reicht
    if(sum(nameindex) > 1)
    {
      ind <- grep(TRUE,nameindex)
      nameindex[ind[2:length(ind)]] <- FALSE
    }
        
    mymain <- paste("Accession:", acc,"Genename: ",namesdf$Gene_names[nameindex],
                    "\n",namesdf$Protein_names[nameindex],"\n","Number of Peptides:", sorteddf$TotalPep[pos])
    plot(range(1:protlength), range(1:lyaxe), type = "n", main = mymain,
         xlab = "Amino acid position", ylab = "",  sub = sub,  yaxt = "n", xaxt = "n")
    
    #add numebers at the beginning and at the end
    #myat <- seq(from = 1, to = protlength, by = 100)
    lenIntervall <- (if (protlength < 51) 10 
                     else if (protlength < 101) 50 
                     else if (protlength < 501) 100 
                     else if (protlength < 2001) 200 
                     else if (protlength < 5001) 500 
                     else if (protlength >= 5001) 1000)
    
    #axis(1, at = c(1,protlength), labels = c(1,protlength), lwd = 1, cex.axis = 0.7)
    
    axis(1, at = 1, cex.axis = 1.6)
    axis(1, at = protlength, cex.axis = 1.6)
    axis(1, at = seq(lenIntervall, (protlength-1), by = lenIntervall), cex.axis = 0.6, col.axis = "grey")
    
    #plot runnames
    mylabels <- locationdf$Expname
    myat <-  seq(from = 1, to = lyaxe - (nrlanes*2), by = nrlanes*2)
    axis(2, at = myat, labels = mylabels, las = 2, lwd = 1, cex.axis = 0.7)
    
    #plot topological areas
    topo_res <- getTopofeatures(GffOrg, protlength, acc)
    if(topo_res[[1]] == FALSE)
    {
      #cat("Could not get topological features for Accession: ", acc, "\n")
      #next
    }
    topovect <- topo_res[[2]] 
    topodf <- topo_res[[3]] 
    print_topoareas(topodf, lyaxe+1)
    
    
    #--------------------------------------------
    #print peptide lines
    # for each run
    runc <- 0
    
    
    for (run in locationdf$Expname)
    {
      
      #lane container
      lanes <- rep(-1, nrlanes)
           
      stend <- intenscountdf[ (intenscountdf$Accession == acc) & (intenscountdf[ ,paste0("Intensity_",run)] >= 0)
                              & (intenscountdf[ ,paste0("Counts_",run)] > 0), 
                              c("Start_position","End_position",paste0("Intensity_",run),"Sequence")]
      
      stend <- stend[ order(stend$Start_position, stend$End_position), ]
      
      # for each peptide
      for (p in seq_along(stend[ , 1]))
      {
        x <- stend[p , c("Start_position","End_position")]
        
        # calculate y axe
        lanepos <- which(lanes < (stend$Start_position[p]-1))        
        
        if(length(lanepos)== 0)
        {
          #cat("No free lane left! Run: ", run, " Acc: ", acc,"\n")          
          next
        }
        
        lanes[lanepos[1]] <- stend$End_position[p]
        y <- rep((runc*nrlanes*2)+lanepos[1], 2)         
       
        mycol <- coldf[coldf$Intval ==  stend[p, paste0("Intensity_",run)], "Color"]           
                
        lines(x, y, col = mycol, lwd = 1)
        
        #-------------------------------------
        # print semi tryptic cleavage site
        # see above
        # protlength <- getLength(FastaOrg[[AccOrg[AccOrg$actProteinAcc == acc, "index"]]]) 
        # protseq <- getSequence(FastaOrg[[AccOrg[AccOrg$actProteinAcc == acc, "index"]]])
        # N-terminal
              
        if (!grepl("R$|K$", substr(protseq, x[1]-1, x[1]-1)))
        {
          lines(c(x[1], x[1]), c(y[1]-0.5, y[1]+0.5), col = "black", lwd = 1)
        }      
        
        # C-terminal
        pepseq <- stend[p, "Sequence"]
        if(!grepl("R$|K$", pepseq) & (x[2] < protlength))
        {
          lines(c(x[2], x[2]), c(y[2]-0.5, y[2]+0.5), col = "black", lwd = 1)
        }
        
      }
      
      runc <- runc+1  
    }
    
    #---------
    #start here drawing of structural Info
    my_y <- runc*nrlanes*2
    
    #print helix rectangles    
    helixdf <- GffOrg[(GffOrg$seqname == acc & GffOrg$feature == "Helix"), c("start","end")]
    
    for (h in seq_along(helixdf[ , 1]))
    {
      rect(helixdf[h, "start"], my_y, helixdf[h, "end"], my_y+3, col = "olivedrab", border = NA)      
    }
    
    # print beta strand
    betadf <- GffOrg[(GffOrg$seqname == acc & GffOrg$feature == "Beta strand"), c("start","end")]
    
    for (h in seq_along(helixdf[ , 1]))
    {
      rect(betadf[h, "start"], my_y, betadf[h, "end"], my_y+3, col = "orange", border = NA)      
    }
    
    # print turn
    turndf <- GffOrg[(GffOrg$seqname == acc & GffOrg$feature == "Turn"), c("start","end")]
    
    for (h in seq_along(turndf[ , 1]))
    {
      rect(turndf[h, "start"], my_y, turndf[h, "end"], my_y+3, col = "black", border = NA)      
    }
    
    #---------
    
    par(mar= c(0,0,0,0), fig=c(0,1,0,1), oma= c(0,0,0,0), new = TRUE)
    #coldf <- data.frame(Intval = indexcol, Color = intenscol)
    
    #legend for the structure    
    legend("topleft",  legend = c("a-Helix","b-Strand","Turn") , col = c("olivedrab","orange","black"), xjust = 1, 
           ncol = 1, cex = 0.5, lwd = 4, text.font = 1, text.col = "black",
           title="Structure", title.col = "black", xpd = TRUE)
    
    #legend for Intensities
    legend("topright",  legend = coldf$Intval , col = coldf$Color, xjust = 1, 
           ncol = 1, cex = 0.5, lwd = 4, text.font = 1, text.col = "black",
           title="log2 Intensities", title.col = "black", xpd = TRUE)
    #legend for TM
    rectAlpha <- 0.1
    rectRed <- rgb(1,0,0,rectAlpha)
    rectBlue <- rgb(0,0,1,rectAlpha)
    rectGreen <- rgb(0,1,0,rectAlpha)   
    rectYellow <- rgb(1,1,0,rectAlpha)
    mylegd <- c("ExtraCellular", "SignalPeptide", "CytoPlasm", "TransMembrane")
    mycol <- c(rectRed, rectBlue, rectGreen, rectYellow)    
    legend("bottomright",  legend = mylegd , col = mycol, xjust = 1, 
           ncol = 1, cex = 0.5, lwd = 6, text.font = 1, text.col = "black",
           title="Locations", title.col = "black", xpd = TRUE)
    
    par(pdfpardefault)
  }  
  
  
  #set back output parameters
  #cat("print_peptides parnew3:",par("new"),"\n")
  #cat("print_peptides generaldefault2:",generalpardefault[["new"]],"\n")
  dev.off()
  #cat("print_peptides parnew3:",par("new"),"\n")
  #cat("print_peptides generaldefault2:",generalpardefault[["new"]],"\n")
  par(generalpardefault)
  #cat("print_peptides generaldefault3:",generalpardefault[["new"]],"\n")
  #cat("print_peptides parnew4:",par("new"),"\n")
  result <- TRUE
  return(result)
  
}



#print fucntion for topological domains
# Input: data.frame topodf, numerical vector maxhight
# Output: non
# device should be on

print_topoareas <- function(topodf,maxh)
{
  # set basic rectangle parameters
  rectBorder <- FALSE
  rectAlpha <- 0.1
  rectRed <- rgb(1,0,0,rectAlpha)
  rectBlue <- rgb(0,0,1,rectAlpha)
  rectGreen <- rgb(0,1,0,rectAlpha)
  #rectYellow <- rgb(255,215,0,rectAlpha, maxColorValue = 255)
  rectYellow <- rgb(1,1,0,rectAlpha)
  
  for (i in seq_along(topodf[,1]))
  {
        
    #print signal peptide in blue
    if (topodf$feature[i] == "SP")
    {
      #rect(topodf[i,"start"],0,topodf[i,"end"],maxh, col = rectBlue, lwd = 0, border = rectBorder)  
      rect(topodf[i,"start"],0,topodf[i,"end"],maxh, col = rectBlue, lwd = 0, border = rectBlue)  
    }  
    
    
    #print extracellular in red
    if (topodf$feature[i] == "EC")
    {
      #rect(topodf[i,"start"],0,topodf[i,"end"],maxh, col = rectRed, lwd = 0, border = rectBorder)  
      rect(topodf[i,"start"],0,topodf[i,"end"],maxh, col = rectRed, lwd = 0, border = rectRed)
    }  
    
    
    #print cytoplasma in green
    if (topodf$feature[i] == "CP")
    {
      #rect(topodf[i,"start"],0,topodf[i,"end"],maxh, col = rectGreen, lwd = 0, border = rectBorder)  
      rect(topodf[i,"start"],0,topodf[i,"end"],maxh, col = rectGreen, lwd = 0, border = rectGreen)  
    }  
    
    
    #print transmembrane in yellow
    if (topodf$feature[i] == "TM")
    {
      #rect(topodf[i,"start"],0,topodf[i,"end"],maxh, col = rectYellow, lwd = 0, border = rectBorder)  
      rect(topodf[i,"start"],0,topodf[i,"end"],maxh, col = rectYellow, lwd = 0, border = rectYellow)  
    }  
    
    
  }
  
}
