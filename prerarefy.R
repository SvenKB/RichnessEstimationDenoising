.libPaths("./R/library/")

args = commandArgs(TRUE)
input1 = args[1]
input2 = args[2]
input3 = args[3]
input4 = as.numeric(args[4])


print(input1)
print(input2)
print(input3)
print(input4)


library(dada2)
library(ShortRead)
library(parallel)

## Prepare function to rarefy fastq files
   generateRarefiedFastq <- function(path,out.path,depth=1000,seed=1406) {
    
    f <- FastqSampler(path, depth,ordered = T)
    set.seed(seed)
    r <- yield(f)    # sample of size n=50
    close(f)
    name <- strsplit(basename(path),".fastq")[[1]][[1]]
    
    dir.create(paste0(out.path,"/Rarefied/"), showWarnings = T)
    dir.create(paste0(out.path,"/Rarefied/rf",depth), showWarnings = T)
    
    curr_path <- paste0(out.path,"/Rarefied/rf",depth,"/",name,depth,".fastq")
    
    writeFastq(object = r,file = curr_path,mode = "w")
  }
  
  ##########################
  #### A - PRJNA679485  ####
  ##########################
  
  ### Set path
  
  path <- input1
  out.path <- path
  
  
  pattern_fs <- input2
  pattern_rs <- input3
  depths <- input4
  
  # Forward and reverse fastq filenames have format: SAMPLENAME_1.fastq.gz and SAMPLENAME_2.fastq.gz
  
  fnFs <- sort(list.files(path, pattern=paste0(pattern_fs,".fastq.gz"), full.names = TRUE))
  print(fnFs)
  fnRs <- sort(list.files(path, pattern=paste0(pattern_rs,".fastq.gz"), full.names = TRUE))
  print(fnRs)
  # Rarefy using prepared function
  
  nc <- 36
  print(nc)
  cl <- makeForkCluster(nc)

  parLapply(cl=cl,fnFs, generateRarefiedFastq,depth=depths,out.path=out.path)
  parLapply(cl=cl,fnRs, generateRarefiedFastq,depth=depths,out.path=out.path)
  stopCluster(cl)
  