# Secretome Lysat Peptide Descriptor
#
# functions to get infos from the gff-data.frame

# getTopofeatures extract a vector that shows the belonging to the topological domains for
# every aminoacid of the protein
# Input: data.frame GffHuman, numerical vector proteinlength, character vector accession
# Output: logical vector result, character vector topo, data.frame topodf

getTopofeatures <- function(gff,protl, acc)
{
  topo <- rep("no",protl)
  topodf <- data.frame(feature=character(), start= numeric(), end=numeric())
  mygff <- gff[gff$seqname == acc, ]
  result <- FALSE
      
  signal <- mygff[mygff$feature == "Signal peptide", c("feature", "start", "end")]  
  if(nrow(signal) >0)
  {
    for (i in seq_along(signal$start))
    {
      topo[signal$start[1]:signal$end[i]] <- "SP"
    }
    signal$feature <- "SP"
    topodf <- rbind(topodf,signal) 
    result <- TRUE
  }  
  
  
  transm <- mygff[mygff$feature == "Transmembrane", c("feature","start", "end")]  
  if(nrow(transm) >0)
  {
    for (i in seq_along(transm$start))
    {
      topo[transm$start[1]:transm$end[i]] <- "TM"
    }
    transm$feature <- "TM"
    topodf <- rbind(topodf,transm)
    result <- TRUE
  }  
  
  
  extra <- mygff[mygff$feature == "Topological domain" & mygff$Note == "Extracellular", c("feature","start", "end")]  
  if(nrow(extra) >0)
  {
    for (i in seq_along(extra$start))
    {
      topo[extra$start[1]:extra$end[i]] <- "EC"
    }
    extra$feature <- "EC"
    topodf <- rbind(topodf,extra)
    result <- TRUE
  }  
  
  
  cyto <- mygff[mygff$feature == "Topological domain" & mygff$Note == "Cytoplasmic", c("feature","start", "end")]  
  if(nrow(cyto) >0)
  {
    for (i in seq_along(cyto$start))
    {
      topo[cyto$start[1]:cyto$end[i]] <- "CP"
    }
    cyto$feature <- "CP"
    topodf <- rbind(topodf,cyto)
    result <- TRUE
  }  
  
  
  
  return <- list(result,topo,topodf)
}