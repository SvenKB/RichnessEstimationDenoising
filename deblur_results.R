
library(tidyverse)
library(qiime2R)
library(phyloseq)
library(ggpubr)

`%notin%` <- Negate(`%in%`)

#### Prepare functions

addDiverstiy <- function(phyl) {
  
  comm <- t(phyloseq::otu_table(phyl))
  samdf <- phyloseq::sample_data(phyl)
  
  samdf$Richness <- hillR::hill_taxa(comm=comm,MARGIN=1,q=0)
  samdf$Shannon <-  hillR::hill_taxa(comm=comm,MARGIN=1,q=1)
  samdf$Simpson <-  hillR::hill_taxa(comm=comm,MARGIN=1,q=2)
  samdf$ace <-  apply(comm,1,fossil::ACE)
  samdf$chao1 <- apply(comm,1,fossil::chao1)
  samdf$goods <- QsRutils::goods(comm)[[3]]
  samdf$margalef <- apply(comm,1,abdiv::margalef)
  samdf$menhinick <- apply(comm,1,abdiv::menhinick)
  samdf$singleton <- apply(comm,1,function(x) sum(x==1))
  
  phyloseq::sample_data(phyl) <- samdf
  
  return(phyl)
  
}


prepareDeblur <- function(phyl,deblur_otu,rf = "Original",id = "A",single = "No") { 
  
  tax_table(phyl) <- NULL
  otu_table(phyl) <- otu_table(deblur_otu,taxa_are_rows=T)
  phyl <- addDiverstiy(phyl)
  
  if(rf == "Original") { 
  sample_data(phyl)$orig_depth <- sample_sums(phyl)
  } else { sample_data(phyl)$orig_depth <- NA}
  
  
  df <- sample_data(phyl) %>% as_tibble() %>% mutate(rf_type = rf,
                                                     depth = sample_sums(phyl),
                                                     study_id = id,
                                                     with_singletons = single)
  
  return(list(df,
              phyl))
} 


deblurRarefy <- function(phyl,depth=11150,single = "No",rf="Post",id = "A") {
  
  phyl_rf <- rarefy_even_depth(phyl,sample.size = depth) %>% addDiverstiy()
  df <- sample_data(phyl_rf) %>% as.tibble() %>% mutate(rf_type = rf,
                                                        depth = sample_sums(phyl_rf),
                                                        with_singletons = single,
                                                        study_id = id)
  return(df)
}


prepareAll <- function(phyl_path,
                       otu_path,
                       otu_rf_path,
                       id = "A",
                       depth = 11150,
                       single = "No") {
  
  orig <- readRDS(phyl_path)[[1]]
  
  
  otu <- qiime2R::read_qza(otu_path)$data %>% as.data.frame()
  otu_rf <- qiime2R::read_qza(otu_rf_path)$data %>% as.data.frame()
  
  deb <- prepareDeblur(orig,otu,id=id,single = single)[[1]]
  phyl <- prepareDeblur(orig,otu,id=id,single = single)[[2]]
  deb_pre <- prepareDeblur(orig,otu_rf,rf = "Pre",id=id,single = single)[[1]] 
  deb_post <- deblurRarefy(phyl,depth=depth,single=single,id=id)
  
  
  out <- reduce(list(deb,
                     deb_pre,
                     deb_post
  ), bind_rows)
  
  return(out)
  
}


myRFplot <- function(out) {
  
  ggiNEXT(out, type=1) + theme_minimal() + xlab("Sequencing depth") -> p
  
  plot_df <- p$data %>% filter(method=="interpolated")
  
  meta_df <-  plot_df %>% group_by(site) %>% summarize(x = max(x),
                                                       y=max(y))
  
  ggplot(aes(x=x,y=y,col=site),data=plot_df) +
    geom_line(size=1.5) +
    #geom_label(aes(label=site),data=meta_df,hjust=0,color="black") +
    geom_ribbon(aes(ymin=y.lwr,ymax=y.upr,fill=site),color=0,alpha=.2) +
    theme_minimal() +
    scale_x_continuous(labels = function(x) format(x, scientific = F)) +
    theme(legend.position = "none",
          axis.title.y = element_blank()) +
    xlab("Sequencing depth") -> pout
  
  return(pout)
}

###################
#### Dataset A ####
###################


A <- prepareAll(phyl_path = "A/results/unpooled_noprior_A.RDS",
           otu_path = "Data/Deblur/table-deblur_A.qza",
           otu_rf_path = "Data/Deblur/table-deblur_A_rf.qza",
           id = "A",
           depth = 11150,
           single = "No")

A_singleton <- prepareAll(phyl_path = "A/results/unpooled_noprior_A.RDS",
                          otu_path = "Data/Deblur/table-deblur_A_singleton.qza",
                          otu_rf_path = "Data/Deblur/table-deblur_A_rf_singleton.qza",
                          id = "A",
                          depth = 11150,
                          single = "Yes")


###################
#### Dataset B ####
###################

B <- prepareAll(phyl_path = "B/results/unpooled_noprior_B.RDS",
                otu_path = "Data/Deblur/table-deblur_B.qza",
                otu_rf_path = "Data/Deblur/table-deblur_B_rf.qza",
                id = "B",
                depth = 21092,
                single = "No")

B_singleton <- prepareAll(phyl_path = "B/results/unpooled_noprior_B.RDS",
                          otu_path = "Data/Deblur/table-deblur_B_singleton.qza",
                          otu_rf_path = "Data/Deblur/table-deblur_B_rf_singleton.qza",
                          id = "B",
                          depth = 21092,
                          single = "Yes")

###################
#### Dataset D ####
###################

D <- prepareAll(phyl_path = "D/results/unpooled_noprior_D.RDS",
                otu_path = "Data/Deblur/table-deblur_D.qza",
                otu_rf_path = "Data/Deblur/table-deblur_D_rf.qza",
                id = "D",
                depth = 21092,
                single = "No")

D_singleton <- prepareAll(phyl_path = "D/results/unpooled_noprior_D.RDS",
                          otu_path = "Data/Deblur/table-deblur_D_singleton.qza",
                          otu_rf_path = "Data/Deblur/table-deblur_D_rf_singleton.qza",
                          id = "D",
                          depth = 21092,
                          single = "Yes")



#### Cmbine data

full_df <- reduce(list(A,
                       A_singleton,
                       B,
                       B_singleton,
                       D,
                       D_singleton),bind_rows) %>% mutate(study_name = case_when(study_id == "A" ~ "Almeida-Santos et al., 2021",
                                                                                 study_id == "B" ~"Manor et al., 2018",
                                                                                 study_id == "D" ~ "Aho et al., 2019"))




#### Update original (not sub-sampled) sequencing depth across scenarios
full_df$orig_depth[full_df$rf_type=="Pre"] <- NA
full_df %>% group_by(Subject) %>% fill(orig_depth) -> full_df

rm(list=setdiff(ls(),"full_df"))


############################################
#### Preprare correlation table & plots ####
############################################


#### Table 
full_df %>% 
  group_by(study_name,with_singletons,rf_type) %>% 
  summarize(rho = cor(Richness,orig_depth,method = "spearman"),
            Mean = round(mean(Richness),2),
            SD = round(sqrt(var(Richness)),2),
            Singletons = round(mean(singleton),2)) %>%
  mutate(Mean = paste0(Mean," (",SD,")"),.after = rho,.keep = "unused") %>% xtable::xtable()



#### Plots
ggplot(aes(x=orig_depth,y=Richness,group=rf_type,color=rf_type),data=full_df %>% filter(study_id== "A" & orig_depth < 250000)) +
  geom_point(alpha=.3) +
  geom_smooth(method="lm") +
  facet_wrap(~with_singletons,scales="free") +
  theme_minimal() +
  scale_colour_brewer(palette = "Set1", name= "Sub-sampling strategy", labels=c("None","Post-ASV","Pre-ASV")) +
  theme(axis.title.x = element_blank()) +
  ylab("Observed richness") -> p_A

ggplot(aes(x=orig_depth,y=Richness,group=rf_type,color=rf_type),data=full_df %>% filter(study_id== "B" & orig_depth < 400000)) +
  geom_point(alpha=.3) +
  geom_smooth(method="lm") +
  facet_wrap(~with_singletons,scales="free") +
  theme_minimal() +
  scale_colour_brewer(palette = "Set1", name= "Sub-sampling strategy", labels=c("None","Post-ASV","Pre-ASV")) +
  theme(axis.title.x = element_blank()) +
  ylab("Observed richness") -> p_B

ggplot(aes(x=orig_depth,y=Richness,group=rf_type,color=rf_type),data=full_df %>% filter(study_id== "D" & orig_depth < 400000)) +
  geom_point(alpha=.3) +
  geom_smooth(method="lm") +
  facet_wrap(~with_singletons,scales="free") +
  theme_minimal() +
  scale_colour_brewer(palette = "Set1", name= "Sub-sampling strategy", labels=c("None","Post-ASV","Pre-ASV")) +
  theme(axis.title.x = element_blank()) +
  ylab("Observed richness") -> p_D


ggarrange(p_A,p_B,p_D,ncol=1,common.legend = T,labels = c("A","B","C")) -> p

ggsave("Plots/corr_allstudies_deblur.PDF",plot = p,device = cairo_pdf,units = "in",width = 8.3,height=11.7)

##########################
#### Prepare boxplots ####
##########################

indices <- c("Richness","Shannon","Simpson","ace","chao1","goods","margalef","menhinick")
names <- c("Observed richness","Shannon","Simpson","ACE","Chao1","Good's coverage","Margalef","Menhinick")

full_df %>%
  filter(study_id=="A") %>%
  rename_at(all_of(indices),~names) %>%
  dplyr::select(study_id,rf_type,with_singletons,orig_depth,all_of(names)) %>%
  pivot_longer(all_of(names),names_to="Index") %>%
  mutate(Index = factor(Index,levels = names)) %>%
  ggplot(aes(x=with_singletons,y=value,fill=rf_type)) +
  geom_point(position=position_jitterdodge(dodge.width=.8,jitter.width = .2), size=1, alpha=0.5, 
             aes(colour=rf_type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_log10() +
  facet_wrap(~Index,nrow=1,scales="free") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1", name= "Sub-sampling strategy", labels=c("None","Post-ASV","Pre-ASV")) +
  scale_colour_brewer(palette = "Set1", name= "Sub-sampling strategy", labels=c("None","Post-ASV","Pre-ASV")) +
  theme(legend.position = "bottom",
        axis.text.x =  element_blank(),
        axis.title.y = element_blank())  -> p_A

full_df %>%
  filter(study_id=="B") %>%
  rename_at(all_of(indices),~names) %>%
  dplyr::select(study_id,rf_type,with_singletons,orig_depth,all_of(names)) %>%
  pivot_longer(all_of(names),names_to="Index") %>%
  mutate(Index = factor(Index,levels = names)) %>%
  ggplot(aes(x=with_singletons,y=value,fill=rf_type)) +
  geom_point(position=position_jitterdodge(dodge.width=.8,jitter.width = .2), size=1, alpha=0.5, 
             aes(colour=rf_type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_log10() +
  facet_wrap(~Index,nrow=1,scales="free") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1", name= "Sub-sampling strategy", labels=c("None","Post-ASV","Pre-ASV")) +
  scale_colour_brewer(palette = "Set1", name= "Sub-sampling strategy", labels=c("None","Post-ASV","Pre-ASV")) +
  theme(legend.position = "bottom",
        axis.text.x =  element_blank(),
        axis.title.y = element_blank()) -> p_B

full_df %>%
  filter(study_id=="D") %>%
  rename_at(all_of(indices),~names) %>%
  dplyr::select(study_id,rf_type,with_singletons,orig_depth,all_of(names)) %>%
  pivot_longer(all_of(names),names_to="Index") %>%
  mutate(Index = factor(Index,levels = names)) %>%
  ggplot(aes(x=with_singletons,y=value,fill=rf_type)) +
  geom_point(position=position_jitterdodge(dodge.width=.8,jitter.width = .2), size=1, alpha=0.5, 
             aes(colour=rf_type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_log10() +
  facet_wrap(~Index,nrow=1,scales="free") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1", name= "Sub-sampling strategy", labels=c("None","Post-ASV","Pre-ASV")) +
  scale_colour_brewer(palette = "Set1", name= "Sub-sampling strategy", labels=c("None","Post-ASV","Pre-ASV")) +
  theme(legend.position = "bottom",
        axis.text.x =  element_blank(),
        axis.title.y = element_blank()) -> p_D


ggarrange(p_A + rremove("xlab"),p_B + rremove("xlab"),p_D,common.legend = T,ncol=1,labels = c("A","B","C")) -> p 

ggsave("Plots/indices_boxplot_deblur.PDF",plot = p,device = cairo_pdf,units = "in",width = 11.7,height=8.3)

###########################
#### Rarefaction plots ####
###########################
library(iNEXT)

#### A ####

A <- prepareDeblur(phyl = readRDS("A/results/unpooled_noprior_A.RDS")[[1]],
                   deblur_otu = qiime2R::read_qza("Data/Deblur/table-deblur_A.qza")$data %>% as.data.frame(),
                   rf = "Original",
                   id = "A",
                   single = "No")[[2]]

A_sig <- prepareDeblur(phyl = readRDS("A/results/unpooled_noprior_A.RDS")[[1]],
                       deblur_otu = qiime2R::read_qza("Data/Deblur/table-deblur_A_singleton.qza")$data %>% as.data.frame(),
                       rf = "Original",
                       id = "A",
                       single = "Yes")[[2]]




A <- subset_samples(A ,sample_sums(A )>10000) %>% otu_table %>% as.data.frame
A_sig <- subset_samples(A_sig,sample_sums(A_sig)>10000) %>% otu_table %>% as.data.frame

out_A <- iNEXT(A, q=c(0), datatype="abundance",endpoint = 40000)
out_Asig <- iNEXT(A_sig, q=c(0), datatype="abundance",endpoint = 40000)


#### B ####

B <- prepareDeblur(phyl = readRDS("B/results/unpooled_noprior_B.RDS")[[1]],
                   deblur_otu = qiime2R::read_qza("Data/Deblur/table-deblur_B.qza")$data %>% as.data.frame(),
                   rf = "Original",
                   id = "B",
                   single = "No")[[2]]

B_sig <- prepareDeblur(phyl = readRDS("B/results/unpooled_noprior_B.RDS")[[1]],
                       deblur_otu = qiime2R::read_qza("Data/Deblur/table-deblur_B_singleton.qza")$data %>% as.data.frame(),
                       rf = "Original",
                       id = "B",
                       single = "Yes")[[2]]




B <- subset_samples(B ,sample_sums(B)>10000) %>% otu_table %>% as.data.frame
B_sig <- subset_samples(B_sig,sample_sums(B_sig)>10000) %>% otu_table %>% as.data.frame

out_B <- iNEXT(B, q=c(0), datatype="abundance",endpoint = 40000)
out_Bsig <- iNEXT(B_sig, q=c(0), datatype="abundance",endpoint = 40000)




#### D ####
D <- prepareDeblur(phyl = readRDS("D/results/unpooled_noprior_D.RDS")[[1]],
                   deblur_otu = qiime2R::read_qza("Data/Deblur/table-deblur_D.qza")$data %>% as.data.frame(),
                   rf = "Original",
                   id = "D",
                   single = "No")[[2]]

D_sig <- prepareDeblur(phyl = readRDS("D/results/unpooled_noprior_D.RDS")[[1]],
                       deblur_otu = qiime2R::read_qza("Data/Deblur/table-deblur_D_singleton.qza")$data %>% as.data.frame(),
                       rf = "Original",
                       id = "D",
                       single = "Yes")[[2]]




D <- subset_samples(D ,sample_sums(D )>10000) %>% otu_table %>% as.data.frame
D_sig <- subset_samples(D_sig,sample_sums(D_sig)>10000) %>% otu_table %>% as.data.frame

out_D <- iNEXT(D, q=c(0), datatype="abundance",endpoint = 40000)
out_Dsig <- iNEXT(D_sig, q=c(0), datatype="abundance",endpoint = 40000)


myRFplot(out_A) -> A1
myRFplot(out_Asig) -> A2
myRFplot(out_B) -> B1
myRFplot(out_Bsig) -> B2
myRFplot(out_D) -> D1
myRFplot(out_Dsig) -> D2

library(patchwork)


p_A <- A1 / A2
p_B <- B1 / B2
p_D <- D1 / D2

p_A | p_B | p_D -> p

plot[[1]] <- plot[[1]] + plot_layout(tag_level = "new")
plot[[2]] <- plot[[2]] + plot_layout(tag_level = "new")
plot[[3]] <- plot[[3]] + plot_layout(tag_level = "new")


plot + plot_annotation(tag_levels = c("A","1")) -> plot

ggsave("Plots/rarefy_curve_all_deblur.PDF",plot = plot,device = cairo_pdf,units = "in",width = 8.3,height=6)
