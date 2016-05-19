# Secretome Lysat Peptide Descriptor
#
# general controll functions
# for
# file structure and so on


# Input: path (character)
# creates filestructure for all the runs
# path as main dir
# |-BasicData
#     |-Uniprotversion folders   
# |-AnalysisData
#     |-Experiment  (name choosen by user)
#       

setDirStruc <- function (path)
{
  
  basicd <- paste0(path,"/BasicData")
  anad <- paste0(path,"/AnalysisData")
  
  if (dir.exists(path))
  {
    print(paste0("Directory already exists: ",path))
    if (!dir.exists(basicd))
    {
      print(paste0("Create directory: ", basicd))
      dir.create(basicd)
    }
    else
    {
      print(paste0("Directory already exists: ", basicd))
    }
    
    if (!dir.exists(anad))
    {
      print(paste0("Create directory: ", anad))
      dir.create(anad)
    }
    else
    {
      print(paste0("Directory already exists: ", anad))
    }
    
  }
  else
  {
    print("Create directories:")
    print(path)
    print(basicd)
    print(anad)
    
    dir.create(path)  
    dir.create(basicd)
    dir.create(anad)
  }
}
