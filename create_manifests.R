files <- list.files("/media/sven/Data/Projects/PreRarefaction/A/filtered")[1:100]

sample_id <- unique(sapply(strsplit(files, "_"), `[`, 1)[1:100])
direction <- rep(c("forward","reverse"),50)

manifest <- list("sample-id" = sample_id,
                       "forward-absolute-filepath" = paste0("/media/sven/Data/Projects/PreRarefaction/A/filtered/",sample_id,"_F_filt.fastq.gz"),
                       "reverse-absolute-filepath" = paste0("/media/sven/Data/Projects/PreRarefaction/A/filtered/",sample_id,"_R_filt.fastq.gz"))

write.table(manifest,"/media/sven/Data/Projects/PreRarefaction/A/filtered/PairedEndManifest",quote=F,row.names=F,sep="\t")



files <- list.files("/media/sven/Data/Projects/PreRarefaction/A/Rarefied/")

sample_id <- unique(sapply(strsplit(files, "_"), `[`, 1)[1:100])

manifest <- list("sample-id" = sample_id,
                 "forward-absolute-filepath" = paste0("/media/sven/Data/Projects/PreRarefaction/A/Rarefied/",sample_id,"_F_filt.fastq.gz"),
                 "reverse-absolute-filepath" = paste0("/media/sven/Data/Projects/PreRarefaction/A/Rarefied/",sample_id,"_R_filt.fastq.gz"))

write.table(manifest,"/media/sven/Data/Projects/PreRarefaction/A/Rarefied/Manifest_A_rarefied",
            quote=F,row.names=F,sep="\t",
            col.names = c("sample-id",
                          "forward-absolute-filepath",
                          "reverse-absolute-filepath"))


files <- list.files("/media/sven/Data/Projects/PreRarefaction/D/filtered/")

sample_id <- unique(sapply(strsplit(files, "_"), `[`, 1))
sample_id <- subset(sample_id,substr(sample_id,nchar(sample_id)-1,nchar(sample_id)) %in% c("-N"))

manifest <- list("sample-id" = sample_id,
                 "forward-absolute-filepath" = paste0("/media/sven/Data/Projects/PreRarefaction/D/filtered/",sample_id,"_F_filt.fastq.gz"),
                 "reverse-absolute-filepath" = paste0("/media/sven/Data/Projects/PreRarefaction/D/filtered/",sample_id,"_R_filt.fastq.gz"))

write.table(manifest,"/media/sven/Data/Projects/PreRarefaction/D/filtered/Manifest_D_full",
            quote=F,row.names=F,sep="\t",
            col.names = c("sample-id",
                          "forward-absolute-filepath",
                          "reverse-absolute-filepath"))


files <- list.files("/media/sven/Data/Projects/PreRarefaction/D/Rarefied/filtered/")

sample_id <- unique(sapply(strsplit(files, "_"), `[`, 1))
sample_id <- subset(sample_id,substr(sample_id,nchar(sample_id)-1,nchar(sample_id)) %in% c("-N"))

manifest <- list("sample-id" = sample_id,
                 "forward-absolute-filepath" = paste0("/media/sven/Data/Projects/PreRarefaction/D/Rarefied/filtered/",sample_id,"_F_filt.fastq.gz"),
                 "reverse-absolute-filepath" = paste0("/media/sven/Data/Projects/PreRarefaction/D/Rarefied/filtered/",sample_id,"_R_filt.fastq.gz"))

write.table(manifest,"/media/sven/Data/Projects/PreRarefaction/D/Rarefied/filtered/Manifest_D_rarefied",
            quote=F,row.names=F,sep="\t",
            col.names = c("sample-id",
                          "forward-absolute-filepath",
                          "reverse-absolute-filepath"))


files <- list.files("/media/sven/Data/Projects/PreRarefaction/B/filtered/")

sample_id <- unique(sapply(strsplit(files, "_"), `[`, 1))

manifest <- list("sample-id" = sample_id,
                 "forward-absolute-filepath" = paste0("/media/sven/Data/Projects/PreRarefaction/B/filtered/",sample_id,"_F_filt.fastq.gz"),
                 "reverse-absolute-filepath" = paste0("/media/sven/Data/Projects/PreRarefaction/B/filtered/",sample_id,"_R_filt.fastq.gz"))

write.table(manifest,"/media/sven/Data/Projects/PreRarefaction/B/filtered/Manifest_B",
            quote=F,row.names=F,sep="\t",
            col.names = c("sample-id",
                          "forward-absolute-filepath",
                          "reverse-absolute-filepath"))

files <- list.files("/media/sven/Data/Projects/PreRarefaction/B/Rarefied/filtered/")

sample_id <- unique(sapply(strsplit(files, "_"), `[`, 1))

manifest <- list("sample-id" = sample_id,
                 "forward-absolute-filepath" = paste0("/media/sven/Data/Projects/PreRarefaction/B/Rarefied/filtered/",sample_id,"_F_filt.fastq.gz"),
                 "reverse-absolute-filepath" = paste0("/media/sven/Data/Projects/PreRarefaction/B/Rarefied/filtered/",sample_id,"_R_filt.fastq.gz"))

write.table(manifest,"/media/sven/Data/Projects/PreRarefaction/B/Rarefied/filtered/Manifest_B_rarefied",
            quote=F,row.names=F,sep="\t",
            col.names = c("sample-id",
                          "forward-absolute-filepath",
                          "reverse-absolute-filepath"))


