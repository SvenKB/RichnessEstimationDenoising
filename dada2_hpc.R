#### Load packages ####

.libPaths("~/R/library/")
library(dada2)
library(ShortRead)
library(phyloseq)
library(tidyverse)


args = commandArgs(TRUE)
curr_path <- args[1]
curr_rf <- args[2]
curr_pool <- as.logical(args[3])
pattern1 <- args[4]
pattern2 <- args[5]
cores <- as.numeric(args[6])
curr_prior <- as.logical(args[7])
curr_priorset <- args[8]

print(curr_path)
print(curr_rf)
print(curr_pool)
print(pattern1)
print(pattern2)
print(cores)
print(curr_prior)
print(curr_priorset)

pooled <- ifelse(curr_pool==T,"pooled","unpooled")

## Prepare pipeline


curr_priorset <- readRDS(curr_priorset)


print(head(curr_priorset))

singletons <- function(x) {sum(x==1)}

myPipeline <- function(path,
                       rarefaction_level = "",
                       pool=F,
                       phyl_out = F,
                       prior = character(0),
                       return_priors=F,
                       self_contained_prior=F,
                       pattern_fs = "",
                       pattern_rs = "") {
  
  
  #path <- paste0(path,rarefaction_level)
  
  # Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
  fnFs <- sort(list.files(path, pattern=paste0(pattern_fs,rarefaction_level,".fastq"), full.names = TRUE))
  print(fnFs)
  fnRs <- sort(list.files(path, pattern=paste0(pattern_rs,rarefaction_level,".fastq"), full.names = TRUE))
  print(fnRs)
  # Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
  sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
  print(sample.names)
  
  # Place filtered files in filtered/ subdirectory
  filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
  names(filtFs) <- sample.names
  names(filtRs) <- sample.names
  
  print(filtFs)
  print(filtRs)
  
  print(length(filtFs))
  print(length(filtRs))
  
  print(length(unique(filtFs)))
  print(length(unique(filtRs)))
  
  
  out <- filterAndTrim(fwd=fnFs, filt=filtFs, rev=fnRs, filt.rev=filtRs, truncLen=c(250,250),
                       maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                       compress=TRUE, multithread=cores,trimLeft = c(17,21))
  
  
  # Learn error rates
  print("Start learning forward error rates:")
  errF <- learnErrors(filtFs, multithread=cores)
  print("Start learning beackward error rates:")
  errR <- learnErrors(filtRs, multithread=cores)
  
  # Dereplicate
  print("Start dereplication forward")
  derepFs <- derepFastq(filtFs, verbose=TRUE)
  print("Start dereplication backward")
  derepRs <- derepFastq(filtRs, verbose=TRUE)
  
  
  # Measure number of singletons per sample
  singletonsFs <- lapply(derepFs, function(x) sum(x$uniques == 1))
  singletonsRs <- lapply(derepRs, function(x) sum(x$uniques == 1))
  
  # Sample inference
  
  print("Sample inference forward reads")
  dadaFs <- dada(derepFs, err=errF, multithread=cores,pool = pool,priors = prior)
  print("Sample inference backward reads")
  dadaRs <- dada(derepRs, err=errR, multithread=cores,pool = pool,priors = prior)
  
  
  # Get prior reads for rarefied data
  if(return_priors == T) {
    FsReads <- unique(unlist(lapply(dadaFs,function(x) names(x[[1]]))))
    RsReads <- unique(unlist(lapply(dadaRs,function(x) names(x[[1]]))))
  }
  
  if(self_contained_prior == T) {
    prior <- unique(c(FsReads,RsReads))
    dadaFs <- dada(derepFs, err=errF, multithread=cores,pool = pool,priors = prior)
    dadaRs <- dada(derepRs, err=errR, multithread=cores,pool = pool,priors = prior)
  }
  
  print("Merge reads")
  # Merge forward and reverse reads
  mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
  
  print("Make sequence table")
  # Construct sequence table
  seqtab <- makeSequenceTable(mergers)
  
  print("Remove bimeras")
  # Remove chimeric reads
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=cores, verbose=TRUE)
  
  print("Assign taxonomy")
  # Assign taxonomy
  taxa <- assignTaxonomy(seqtab.nochim, "/scratch/tmp/kleineba/PreRarefaction/Data/tax/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
  
  print("Create sample data and augment diversities")
  # Create sample df
  samples.out <- rownames(seqtab.nochim)
  subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
  #gender <- substr(subject,1,1)
  #subject <- substr(subject,2,999)
  #day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
  samdf <- data.frame(Subject=subject)
  #samdf$When <- "Early"
  #samdf$When[samdf$Day>100] <- "Late"
  rownames(samdf) <- samples.out
  samdf$FsS <- unlist(singletonsFs)
  samdf$RsS <- unlist(singletonsRs)
  samdf$rf_level <- rarefaction_level
  samdf$pooled <- pool
  samdf$singleton <- apply(seqtab.nochim,1,singletons)
  
  
  # Add diversity 
  comm <- seqtab.nochim
  
  samdf$Richness <- hillR::hill_taxa(comm=comm,MARGIN=1,q=0)
  samdf$Shannon <-  hillR::hill_taxa(comm=comm,MARGIN=1,q=1)
  samdf$Simpson <-  hillR::hill_taxa(comm=comm,MARGIN=1,q=2)
  samdf$ace <-  apply(comm,1,fossil::ACE)
  samdf$chao1 <- apply(comm,1,fossil::chao1)
  samdf$goods <- QsRutils::goods(comm)[[3]]
  samdf$margalef <- apply(comm,1,abdiv::margalef)
  samdf$menhinick <- apply(comm,1,abdiv::menhinick)
  
  if (phyl_out == T) {
    ps <- phyloseq(otu_table(t(seqtab.nochim), taxa_are_rows=T), 
                   sample_data(samdf), 
                   tax_table(taxa))
    #ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
    
    out <- ps
  } else {
    out <- samdf
  }
  
  if(return_priors == T) {
    out <- list("out" = out,
                "FsPriors" = FsReads,
                "RsPriors" = RsReads)
  }
  
  return(out)
  
}



out <- myPipeline(path = curr_path,
                  rarefaction_level = curr_rf,
                  phyl_out = T,
                  prior = curr_priorset,
                  return_priors = curr_prior,
                  pool=curr_pool,
                  pattern_fs = pattern1,
                  pattern_rs = pattern2)

dir.create(paste0(curr_path,"/results"))

saveRDS(out,paste0(curr_path,"/results/results_",curr_rf,"_",pooled,".RDS"))