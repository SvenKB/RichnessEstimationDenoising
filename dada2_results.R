#### Prepare priors obtained by full (pooled) pipeline and data for pre-rearefaction 

#### Load unpooled data
A <- readRDS("Data/A/results/unpooled_noprior_A.RDS")
B <- readRDS("Data/B/results/unpooled_noprior_B.RDS")
D <- readRDS("Data/D/results/unpooled_noprior_D.RDS")


#### Obtain minimal sequencing depth to be used as sub-sampling depth 
## Threshold set to at leat 10.000 reads to ensure good quality samples
min(sample_sums(A[[1]])[sample_sums(A[[1]]) > 10000]) # 11151
min(sample_sums(B[[1]])[sample_sums(B[[1]]) > 10000]) # 21092
min(sample_sums(D[[1]])[sample_sums(D[[1]]) > 10000]) # 12341

#### Prepare priors
unique(unlist(A[2:3])) -> A
unique(unlist(B[2:3])) -> B
unique(unlist(D[2:3])) -> D


#### Save priors
saveRDS(A,"Data/A/results/unpooled_priors_A.RDS")
saveRDS(B,"Data/B/results/unpooled_priors_B.RDS")
saveRDS(D,"Data/D/results/unpooled_priors_D.RDS")


#### Load pooled data
A <- readRDS("Data/A/results/pooled_noprior_A.RDS")
B <- readRDS("Data/B/results/pooled_noprior_B.RDS")
D <- readRDS("Data/D/results/pooled_noprior_D.RDS")


#### prepare priors
unique(unlist(A[2:3])) -> A
unique(unlist(B[2:3])) -> B
unique(unlist(D[2:3])) -> D

#### Save priors
saveRDS(A,"Data/A/results/pooled_priors_A.RDS")
saveRDS(B,"Data/B/results/pooled_priors_B.RDS")
saveRDS(D,"Data/D/results/pooled_priors_D.RDS")

####################################
#### Prepare post-rarefied data ####
####################################

##################
#### No Prior ####
##################

#### Load data
A_unpooled <- readRDS("Data/A/results/unpooled_noprior_A.RDS")[[1]]
B_unpooled <- readRDS("Data/B/results/unpooled_noprior_B.RDS")[[1]]
D_unpooled <- readRDS("Data/D/results/unpooled_noprior_D.RDS")[[1]]

A_pooled <- readRDS("Data/A/results/pooled_noprior_A.RDS")[[1]]
B_pooled <- readRDS("Data/B/results/pooled_noprior_B.RDS")[[1]]
D_pooled <- readRDS("Data/D/results/pooled_noprior_D.RDS")[[1]]

#### Sub-sample data to minimal seuqnecing depth
A_unpooled_rf <- phyloseq::rarefy_even_depth(A_unpooled, sample.size = 11150)
B_unpooled_rf <- phyloseq::rarefy_even_depth(B_unpooled, sample.size = 21092)
D_unpooled_rf <- phyloseq::rarefy_even_depth(D_unpooled, sample.size = 12341)

A_pooled_rf <- phyloseq::rarefy_even_depth(A_pooled, sample.size = 11150)
B_pooled_rf <- phyloseq::rarefy_even_depth(B_pooled, sample.size = 21092)
D_pooled_rf <- phyloseq::rarefy_even_depth(D_pooled, sample.size = 12341)


#### Add diversity indices

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

A_unpooled_rf <- addDiverstiy(A_unpooled_rf)
B_unpooled_rf <- addDiverstiy(B_unpooled_rf)
D_unpooled_rf <- addDiverstiy(D_unpooled_rf)
A_pooled_rf <- addDiverstiy(A_pooled_rf)
B_pooled_rf <- addDiverstiy(B_pooled_rf)
D_pooled_rf <- addDiverstiy(D_pooled_rf)
  
#### save Sub-sampled data
saveRDS(A_unpooled_rf,"Data/A/results/unpooled_noprior_A_postrf.RDS")
saveRDS(B_unpooled_rf,"Data/B/results/unpooled_noprior_B_postrf.RDS")
saveRDS(D_unpooled_rf,"Data/D/results/unpooled_noprior_D_postrf.RDS")

saveRDS(A_pooled_rf,"Data/A/results/pooled_noprior_A_postrf.RDS")
saveRDS(B_pooled_rf,"Data/B/results/pooled_noprior_B_postrf.RDS")
saveRDS(D_pooled_rf,"Data/D/results/pooled_noprior_D_postrf.RDS")


###############
#### Prior ####
###############

#### Load data
A_unpooled <- readRDS("Data/A/results/unpooled_prior_A.RDS")
B_unpooled <- readRDS("Data/B/results/unpooled_prior_B.RDS")
D_unpooled <- readRDS("Data/D/results/unpooled_prior_D.RDS")

A_pooled <- readRDS("Data/A/results/pooled_prior_A.RDS")
B_pooled <- readRDS("Data/B/results/pooled_prior_B.RDS")
D_pooled <- readRDS("Data/D/results/pooled_prior_D.RDS")

#### Subsample data to minimal seuqnecing depth
A_unpooled_rf <- phyloseq::rarefy_even_depth(A_unpooled, sample.size = 11150)
B_unpooled_rf <- phyloseq::rarefy_even_depth(B_unpooled, sample.size = 21092)
D_unpooled_rf <- phyloseq::rarefy_even_depth(D_unpooled, sample.size = 12341)

A_pooled_rf <- phyloseq::rarefy_even_depth(A_pooled, sample.size = 11150)
B_pooled_rf <- phyloseq::rarefy_even_depth(B_pooled, sample.size = 21092)
D_pooled_rf <- phyloseq::rarefy_even_depth(D_pooled, sample.size = 12341)

#### Add diversities

A_unpooled_rf <- addDiverstiy(A_unpooled_rf)
B_unpooled_rf <- addDiverstiy(B_unpooled_rf)
D_unpooled_rf <- addDiverstiy(D_unpooled_rf)

A_pooled_rf <- addDiverstiy(A_pooled_rf)
B_pooled_rf <- addDiverstiy(B_pooled_rf)
D_pooled_rf <- addDiverstiy(D_pooled_rf)


#### save Sub-sampled data
saveRDS(A_unpooled_rf,"Data/A/results/unpooled_prior_A_postrf.RDS")
saveRDS(B_unpooled_rf,"Data/B/results/unpooled_prior_B_postrf.RDS")
saveRDS(D_unpooled_rf,"Data/D/results/unpooled_prior_D_postrf.RDS")

saveRDS(A_pooled_rf,"Data/A/results/pooled_prior_A_postrf.RDS")
saveRDS(B_pooled_rf,"Data/B/results/pooled_prior_B_postrf.RDS")
saveRDS(D_pooled_rf,"Data/D/results/pooled_prior_D_postrf.RDS")


##################################
#### Prepare combined dataset ####
##################################

library(tidyverse)

####################
#### Unrarefied ####
####################



prepareData <- function(path,rf,study_name,study_id,prior) {
  
  phyl <- readRDS(path)
  
  if(is.list(phyl)) {
  phyl <- readRDS(path)[[1]]
  }
  
  samdf <- phyloseq::sample_data(phyl)
  samdf$depth <- phyloseq::sample_sums(phyl)
  
  samdf %>%
    as_tibble() %>%
    mutate(rf_level = rf,
           study_name = study_name,
           study_id = study_id,
           prior = prior,.before="Subject") -> samdf
  
  return(samdf)
  
}


#### Unrarefied / no prior / unpooled

A_unpooled_noprior <- prepareData(path="Data/A/results/unpooled_noprior_A.RDS",
                                  rf = "Original",
                                  study_name = "Almeida-Santos et al., 2021",
                                  study_id = "A",
                                  prior = "No prior")
                      
                      
B_unpooled_noprior <- prepareData("Data/B/results/unpooled_noprior_B.RDS",
                                  rf = "Original",
                                  study_name = "Manor et al., 2018",
                                  study_id = "B",
                                  prior = "No prior")
                      
                      
D_unpooled_noprior <- prepareData("Data/D/results/unpooled_noprior_D.RDS",
                                  rf = "Original",
                                  study_name = "Aho et al., 2019",
                                  study_id = "D",
                                  prior = "No prior")

                      
                      

#### Unrarefied / no prior / pooled
A_pooled_noprior <- prepareData("Data/A/results/pooled_noprior_A.RDS",
                                rf = "Original",
                                study_name = "Almeida-Santos et al., 2021",
                                study_id = "A",
                                prior = "No prior")
                    
                    

B_pooled_noprior <- prepareData("Data/B/results/pooled_noprior_B.RDS",
                                rf = "Original",
                                study_name = "Manor et al., 2018",
                                study_id = "B",
                                prior = "No prior")
                    
                    

D_pooled_noprior <- prepareData("Data/D/results/pooled_noprior_D.RDS",
                                rf = "Original",
                                study_name = "Aho et al., 2019",
                                study_id = "D",
                                prior = "No prior")
                    
                    

#### Unrarefied / unpooled / prior

A_unpooled_prior <- prepareData(path="Data/A/results/unpooled_prior_A.RDS",
                                  rf = "Original",
                                  study_name = "Almeida-Santos et al., 2021",
                                  study_id = "A",
                                  prior = "Prior")


B_unpooled_prior <- prepareData("Data/B/results/unpooled_prior_B.RDS",
                                  rf = "Original",
                                  study_name = "Manor et al., 2018",
                                  study_id = "B",
                                  prior = "Prior")


D_unpooled_prior <- prepareData("Data/D/results/unpooled_prior_D.RDS",
                                  rf = "Original",
                                  study_name = "Aho et al., 2019",
                                  study_id = "D",
                                  prior = "Prior")



#### Unrarefied / prior / pooled

A_pooled_prior <- prepareData(path="Data/A/results/pooled_prior_A.RDS",
                                rf = "Original",
                                study_name = "Almeida-Santos et al., 2021",
                                study_id = "A",
                                prior = "Prior")


B_pooled_prior <- prepareData("Data/B/results/pooled_prior_B.RDS",
                                rf = "Original",
                                study_name = "Manor et al., 2018",
                                study_id = "B",
                                prior = "Prior")


D_pooled_prior <- prepareData("Data/D/results/pooled_prior_D.RDS",
                                rf = "Original",
                                study_name = "Aho et al., 2019",
                                study_id = "D",
                                prior = "Prior")



#### Combine unrarefied data

unrarefied <- reduce(list(A_unpooled_noprior,
                          B_unpooled_noprior,
                          D_unpooled_noprior,
                          A_pooled_noprior,
                          B_pooled_noprior,
                          D_pooled_noprior,
                          A_unpooled_prior,
                          B_unpooled_prior,
                          D_unpooled_prior,
                          A_pooled_prior,
                          B_pooled_prior,
                          D_pooled_prior),bind_rows) %>%
  mutate(orig_depth = depth)


rm(A_unpooled_noprior,
        B_unpooled_noprior,
        D_unpooled_noprior,
        A_pooled_noprior,
        B_pooled_noprior,
        D_pooled_noprior,
        A_unpooled_prior,
        B_unpooled_prior,
        D_unpooled_prior,
        A_pooled_prior,
        B_pooled_prior,
        D_pooled_prior)

######################
#### Pre-rarefied ####
######################

#### Pre-Rarefied / no prior / unpooled

A_unpooled_noprior <- prepareData(path="Data/A/results/unpooled_noprior_A_prerf.RDS",
                                  rf = "Pre",
                                  study_name = "Almeida-Santos et al., 2021",
                                  study_id = "A",
                                  prior = "No prior")


B_unpooled_noprior <- prepareData("Data/B/results/unpooled_noprior_B_prerf.RDS",
                                  rf = "Pre",
                                  study_name = "Manor et al., 2018",
                                  study_id = "B",
                                  prior = "No prior")


D_unpooled_noprior <- prepareData("Data/D/results/unpooled_noprior_D_prerf.RDS",
                                  rf = "Pre",
                                  study_name = "Aho et al., 2019",
                                  study_id = "D",
                                  prior = "No prior")

#### Pre-rarefied / prior / unpooled

A_unpooled_prior <- prepareData(path="Data/A/results/unpooled_prior_A_prerf.RDS",
                                rf = "Pre",
                                study_name = "Almeida-Santos et al., 2021",
                                study_id = "A",
                                prior = "Prior")


B_unpooled_prior <- prepareData("Data/B/results/unpooled_prior_B_prerf.RDS",
                                rf = "Pre",
                                study_name = "Manor et al., 2018",
                                study_id = "B",
                                prior = "Prior")


D_unpooled_prior <- prepareData("Data/D/results/unpooled_prior_D_prerf.RDS",
                                rf = "Pre",
                                study_name = "Aho et al., 2019",
                                study_id = "D",
                                prior = "Prior")


#### Pre-Rarefied / no prior / pooled

A_pooled_noprior <- prepareData("Data/A/results/pooled_noprior_A_prerf.RDS",
                                rf = "Pre",
                                study_name = "Almeida-Santos et al., 2021",
                                study_id = "A",
                                prior = "No prior")



B_pooled_noprior <- prepareData("Data/B/results/pooled_noprior_B_prerf.RDS",
                                rf = "Pre",
                                study_name = "Manor et al., 2018",
                                study_id = "B",
                                prior = "No prior")



D_pooled_noprior <- prepareData("Data/D/results/pooled_noprior_D_prerf.RDS",
                                rf = "Pre",
                                study_name = "Aho et al., 2019",
                                study_id = "D",
                                prior = "No prior")

#### Pre-rarefied / prior / pooled

A_pooled_prior <- prepareData(path="Data/A/results/pooled_prior_A_prerf.RDS",
                              rf = "Pre",
                              study_name = "Almeida-Santos et al., 2021",
                              study_id = "A",
                              prior = "Prior")


B_pooled_prior <- prepareData("Data/B/results/pooled_prior_B_prerf.RDS",
                              rf = "Pre",
                              study_name = "Manor et al., 2018",
                              study_id = "B",
                              prior = "Prior")


D_pooled_prior <- prepareData("Data/D/results/pooled_prior_D_prerf.RDS",
                              rf = "Pre",
                              study_name = "Aho et al., 2019",
                              study_id = "D",
                              prior = "Prior")


#### Combine pre-rarefied data

prerarefied <- reduce(list(A_unpooled_noprior,
                          B_unpooled_noprior,
                          D_unpooled_noprior,
                          A_pooled_noprior,
                          B_pooled_noprior,
                          D_pooled_noprior,
                          A_unpooled_prior,
                          B_unpooled_prior,
                          D_unpooled_prior,
                          A_pooled_prior,
                          B_pooled_prior,
                          D_pooled_prior),bind_rows)

rm(A_unpooled_noprior,
   B_unpooled_noprior,
   D_unpooled_noprior,
   A_pooled_noprior,
   B_pooled_noprior,
   D_pooled_noprior,
   A_unpooled_prior,
   B_unpooled_prior,
   D_unpooled_prior,
   A_pooled_prior,
   B_pooled_prior,
   D_pooled_prior)

#######################
#### Post-rarefied ####
#######################

#### Post-rarefied / no prior / unpooled

A_unpooled_noprior <- prepareData(path="Data/A/results/unpooled_noprior_A_postrf.RDS",
                                  rf = "Post",
                                  study_name = "Almeida-Santos et al., 2021",
                                  study_id = "A",
                                  prior = "No prior")

B_unpooled_noprior <- prepareData(path="Data/B/results/unpooled_noprior_B_postrf.RDS",
                                  rf = "Post",
                                  study_name = "Manor et al., 2018",
                                  study_id = "B",
                                  prior = "No prior")


D_unpooled_noprior <- prepareData(path="Data/D/results/unpooled_noprior_D_postrf.RDS",
                                  rf = "Post",
                                  study_name = "Aho et al., 2019",
                                  study_id = "D",
                                  prior = "No prior")


#### Post-rarefied / prior / unpooled

A_unpooled_prior <- prepareData(path="Data/A/results/unpooled_prior_A_postrf.RDS",
                                rf = "Post",
                                study_name = "Almeida-Santos et al., 2021",
                                study_id = "A",
                                prior = "Prior")


B_unpooled_prior <- prepareData("Data/B/results/unpooled_prior_B_postrf.RDS",
                                rf = "Post",
                                study_name = "Manor et al., 2018",
                                study_id = "B",
                                prior = "Prior")


D_unpooled_prior <- prepareData("Data/D/results/unpooled_prior_D_postrf.RDS",
                                rf = "Post",
                                study_name = "Aho et al., 2019",
                                study_id = "D",
                                prior = "Prior")



#### Post-rarefied / no prior / pooled

A_pooled_noprior <- prepareData("Data/A/results/pooled_noprior_A_postrf.RDS",
                                rf = "Post",
                                study_name = "Almeida-Santos et al., 2021",
                                study_id = "A",
                                prior = "No prior")



B_pooled_noprior <- prepareData("Data/B/results/pooled_noprior_B_postrf.RDS",
                                rf = "Post",
                                study_name = "Manor et al., 2018",
                                study_id = "B",
                                prior = "No prior")



D_pooled_noprior <- prepareData("Data/D/results/pooled_noprior_D_postrf.RDS",
                                rf = "Post",
                                study_name = "Aho et al., 2019",
                                study_id = "D",
                                prior = "No prior")

#### Post-rarefied / prior / pooled

A_pooled_prior <- prepareData(path="Data/A/results/pooled_prior_A_postrf.RDS",
                              rf = "Post",
                              study_name = "Almeida-Santos et al., 2021",
                              study_id = "A",
                              prior = "Prior")


B_pooled_prior <- prepareData("Data/B/results/pooled_prior_B_postrf.RDS",
                              rf = "Post",
                              study_name = "Manor et al., 2018",
                              study_id = "B",
                              prior = "Prior")


D_pooled_prior <- prepareData("Data/D/results/pooled_prior_D_postrf.RDS",
                              rf = "Post",
                              study_name = "Aho et al., 2019",
                              study_id = "D",
                              prior = "Prior")

#### Combine post-rarefied datasets


postrarefied <- reduce(list(A_unpooled_noprior,
                           B_unpooled_noprior,
                           D_unpooled_noprior,
                           A_pooled_noprior,
                           B_pooled_noprior,
                           D_pooled_noprior,
                           A_unpooled_prior,
                           B_unpooled_prior,
                           D_unpooled_prior,
                           A_pooled_prior,
                           B_pooled_prior,
                           D_pooled_prior),bind_rows)

rm(A_unpooled_noprior,
   B_unpooled_noprior,
   D_unpooled_noprior,
   A_pooled_noprior,
   B_pooled_noprior,
   D_pooled_noprior,
   A_unpooled_prior,
   B_unpooled_prior,
   D_unpooled_prior,
   A_pooled_prior,
   B_pooled_prior,
   D_pooled_prior)

#######################
#### Final dataset ####
#######################

df <- reduce(list(unrarefied,
                  postrarefied,
                  prerarefied),bind_rows)


#### Mutate variables

df %>%
  arrange(Subject) %>%
  group_by(pooled) %>%
  fill(orig_depth) %>% ungroup -> df


df %>%
  mutate(pooled = as.factor(ifelse(pooled==TRUE,"Pooled","Unpooled")),
         pooled = factor(pooled, levels =c("Unpooled","Pooled"))) -> df

#### Quality checks and filtering



df <- df %>%
  filter(orig_depth > 10000  & orig_depth < 500000) %>%
  filter(depth > 2000) %>%
  filter(!(study_id=="D" & grepl("O",Subject)))

table(df$Subject) -> filt

filt <- names(filt[filt==12])

df %>% filter(Subject %in% filt) -> df


######################
#### Prepare plots ###
######################


### Boxplots

indices <- c("Richness","Shannon","Simpson","ace","chao1","goods","margalef","menhinick")
names <- c("Observed richness","Shannon","Simpson","ACE","Chao1","Good's coverage","Margalef","Menhinick")



df %>%
  filter(prior=="No prior") %>%
  filter(study_id=="A") %>%
  rename_at(all_of(indices),~names) %>%
  dplyr::select(study_name,prior,rf_level,pooled,singleton,orig_depth,all_of(names)) %>%
  pivot_longer(all_of(names),names_to="Index") %>%
  mutate(Index = factor(Index,levels = names)) %>%
  ggplot(aes(x=pooled,y=value,fill=rf_level)) +
  geom_point(position=position_jitterdodge(dodge.width=.8,jitter.width = .2), size=1, alpha=0.5, 
             aes(colour=rf_level)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_log10() +
  facet_wrap(~Index,ncol=2,scales="free") +
  theme_clean() +
  scale_fill_brewer(palette = "Set1") +
  scale_colour_brewer(palette = "Set1") +
  theme(legend.position = "bottom") -> p_A

ggsave("Plots/indices_boxplot_A.PDF",plot = p_A,device = cairo_pdf,units = "in",width = 8.3,height=11.7)



df %>%
  filter(prior=="No prior") %>%
  filter(study_id=="B") %>%
  rename_at(all_of(indices),~names) %>%
  dplyr::select(study_name,prior,rf_level,pooled,singleton,orig_depth,all_of(names)) %>%
  pivot_longer(all_of(names),names_to="Index") %>%
  mutate(Index = factor(Index,levels = names)) %>%
  ggplot(aes(x=pooled,y=value,fill=rf_level)) +
  geom_point(position=position_jitterdodge(dodge.width=.8,jitter.width = .2), size=1, alpha=0.5, 
             aes(colour=rf_level)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_log10() +
  facet_wrap(~Index,ncol=2,scales="free") +
  theme_clean() +
  scale_fill_brewer(palette = "Set1") +
  scale_colour_brewer(palette = "Set1") +
  theme(legend.position = "bottom") -> p_B

ggsave("Plots/indices_boxplot_B.PDF",plot = p_B,device = cairo_pdf,units = "in",width = 8.3,height=11.7)



df %>%
  filter(prior=="No prior") %>%
  filter(study_id=="D") %>%
  rename_at(all_of(indices),~names) %>%
  dplyr::select(study_name,prior,rf_level,pooled,singleton,orig_depth,all_of(names)) %>%
  pivot_longer(all_of(names),names_to="Index") %>%
  mutate(Index = factor(Index,levels = names)) %>%
  ggplot(aes(x=pooled,y=value,fill=rf_level)) +
  geom_point(position=position_jitterdodge(dodge.width=.8,jitter.width = .2), size=1, alpha=0.5, 
             aes(colour=rf_level)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_log10() +
  facet_wrap(~Index,ncol=2,scales="free") +
  theme_clean() +
  scale_fill_brewer(palette = "Set1") +
  scale_colour_brewer(palette = "Set1") +
  theme(legend.position = "bottom") -> p_D

ggsave("Plots/indices_boxplot_D.PDF",plot = p_D,device = cairo_pdf,units = "in",width = 8.3,height=11.7)


#### Combined version

df %>%
  filter(prior=="No prior") %>%
  filter(study_id=="A") %>%
  rename_at(all_of(indices),~names) %>%
  dplyr::select(study_name,prior,rf_level,pooled,singleton,orig_depth,all_of(names)) %>%
  pivot_longer(all_of(names),names_to="Index") %>%
  mutate(Index = factor(Index,levels = names)) %>%
  ggplot(aes(x=pooled,y=value,fill=rf_level)) +
  geom_point(position=position_jitterdodge(dodge.width=.8,jitter.width = .2), size=1, alpha=0.5, 
             aes(colour=rf_level)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_log10() +
  facet_wrap(~Index,nrow=1,scales="free") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1",name="Sub-sampling strategy",labels=c("None","Post-ASV","Pre-ASV")) +
  scale_colour_brewer(palette = "Set1",name="Sub-sampling strategy",labels=c("None","Post-ASV","Pre-ASV")) +
  theme(legend.position = "bottom",
        axis.text.x =  element_blank(),
        axis.title.y = element_blank()) -> p_A




df %>%
  filter(prior=="No prior") %>%
  filter(study_id=="B") %>%
  rename_at(all_of(indices),~names) %>%
  dplyr::select(study_name,prior,rf_level,pooled,singleton,orig_depth,all_of(names)) %>%
  pivot_longer(all_of(names),names_to="Index") %>%
  mutate(Index = factor(Index,levels = names)) %>%
  ggplot(aes(x=pooled,y=value,fill=rf_level)) +
  geom_point(position=position_jitterdodge(dodge.width=.8,jitter.width = .2), size=1, alpha=0.5, 
             aes(colour=rf_level)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_log10() +
  facet_wrap(~Index,nrow=1,scales="free") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1",name="Sub-sampling strategy",labels=c("None","Post-ASV","Pre-ASV")) +
  scale_colour_brewer(palette = "Set1",name="Sub-sampling strategy",labels=c("None","Post-ASV","Pre-ASV")) +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.title.y = element_blank()) -> p_B



df %>%
  filter(prior=="No prior") %>%
  filter(study_id=="D") %>%
  rename_at(all_of(indices),~names) %>%
  dplyr::select(study_name,prior,rf_level,pooled,singleton,orig_depth,all_of(names)) %>%
  pivot_longer(all_of(names),names_to="Index") %>%
  mutate(Index = factor(Index,levels = names)) %>%
  ggplot(aes(x=pooled,y=value,fill=rf_level)) +
  geom_point(position=position_jitterdodge(dodge.width=.8,jitter.width = .2), size=1, alpha=0.5, 
             aes(colour=rf_level)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_log10() +
  facet_wrap(~Index,nrow=1,scales="free") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1",name="Sub-sampling strategy",labels=c("None","Post-ASV","Pre-ASV")) +
  scale_colour_brewer(palette = "Set1",name="Sub-sampling strategy",labels=c("None","Post-ASV","Pre-ASV")) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = -40,hjust=0),
        axis.title.y = element_blank()) -> p_D


ggarrange(p_A + rremove("xlab"),p_B + rremove("xlab"),p_D,common.legend = T,ncol=1,labels = c("A","B","C")) -> p 


ggsave("Plots/indices_boxplot_combined.PDF",plot = p,device = cairo_pdf,units = "in",width = 11.7,height=8.3)



df %>%
  filter(prior=="Prior") %>%
  filter(study_id=="A") %>%
  rename_at(all_of(indices),~names) %>%
  dplyr::select(study_name,prior,rf_level,pooled,singleton,orig_depth,all_of(names)) %>%
  pivot_longer(all_of(names),names_to="Index") %>%
  mutate(Index = factor(Index,levels = names)) %>%
  ggplot(aes(x=pooled,y=value,fill=rf_level)) +
  geom_point(position=position_jitterdodge(dodge.width=.8,jitter.width = .2), size=1, alpha=0.5, 
             aes(colour=rf_level)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_log10() +
  facet_wrap(~Index,nrow=1,scales="free") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1",name="Sub-sampling strategy",labels=c("None","Post-ASV","Pre-ASV")) +
  scale_colour_brewer(palette = "Set1",name="Sub-sampling strategy",labels=c("None","Post-ASV","Pre-ASV")) +
  theme(legend.position = "bottom",
        axis.text.x =  element_blank(),
        axis.title.y = element_blank())-> p_A




df %>%
  filter(prior=="Prior") %>%
  filter(study_id=="B") %>%
  rename_at(all_of(indices),~names) %>%
  dplyr::select(study_name,prior,rf_level,pooled,singleton,orig_depth,all_of(names)) %>%
  pivot_longer(all_of(names),names_to="Index") %>%
  mutate(Index = factor(Index,levels = names)) %>%
  ggplot(aes(x=pooled,y=value,fill=rf_level)) +
  geom_point(position=position_jitterdodge(dodge.width=.8,jitter.width = .2), size=1, alpha=0.5, 
             aes(colour=rf_level)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_log10() +
  facet_wrap(~Index,nrow=1,scales="free") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1",name="Sub-sampling strategy",labels=c("None","Post-ASV","Pre-ASV")) +
  scale_colour_brewer(palette = "Set1",name="Sub-sampling strategy",labels=c("None","Post-ASV","Pre-ASV")) +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.title.y = element_blank())-> p_B



df %>%
  filter(prior=="Prior") %>%
  filter(study_id=="D") %>%
  rename_at(all_of(indices),~names) %>%
  dplyr::select(study_name,prior,rf_level,pooled,singleton,orig_depth,all_of(names)) %>%
  pivot_longer(all_of(names),names_to="Index") %>%
  mutate(Index = factor(Index,levels = names)) %>%
  ggplot(aes(x=pooled,y=value,fill=rf_level)) +
  geom_point(position=position_jitterdodge(dodge.width=.8,jitter.width = .2), size=1, alpha=0.5, 
             aes(colour=rf_level)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_log10() +
  facet_wrap(~Index,nrow=1,scales="free") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1",name="Sub-sampling strategy",labels=c("None","Post-ASV","Pre-ASV")) +
  scale_colour_brewer(palette = "Set1",name="Sub-sampling strategy",labels=c("None","Post-ASV","Pre-ASV")) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = -40,hjust=0),
        axis.title.y = element_blank()) -> p_D


ggarrange(p_A + rremove("xlab"),p_B + rremove("xlab"),p_D,common.legend = T,ncol=1,labels = c("A","B","C")) -> p 


ggsave("Plots/indices_boxplot_combined_prior.PDF",plot = p,device = cairo_pdf,units = "in",width = 11.7,height=8.3)

#############################
#### Prepare tables etc. ####
#############################


### Richness and correlation with original seq.depth


df %>%
  filter(study_id == "A") %>%
  filter(prior == "No prior") %>%
  group_by(pooled,rf_level) %>%
  summarize(rho = cor(Richness,orig_depth,method = "spearman"),
            Richness = round(mean(Richness),2),
            Singletons = round(mean(singleton),2)) %>% xtable::xtable()


df %>%
  filter(study_id == "B") %>%
  filter(prior == "No prior") %>%
  group_by(pooled,rf_level) %>%
  summarize(rho = cor(Richness,orig_depth,method = "spearman"),
            Richness = round(mean(Richness),2),
            Singletons = round(mean(singleton),2))  %>% xtable::xtable()



df %>%
  filter(study_id == "D") %>%
  filter(prior == "No prior") %>%
  group_by(pooled,rf_level) %>%
  summarize(rho = cor(Richness,orig_depth,method = "spearman"),
            Richness = round(mean(Richness),2),
            Singletons = round(mean(singleton),2))  %>% xtable::xtable()

df %>%
  filter(prior == "No prior") %>%
  group_by(study_id,pooled,rf_level) %>%
  summarize(rho = cor(Richness,orig_depth,method = "spearman"),
            Mean = round(mean(Richness),2),
            SD = round(sqrt(var(Richness)),3),
            Singletons = round(mean(singleton),2))  



df %>%
  filter(prior == "Prior") %>%
  group_by(study_id,pooled,rf_level) %>%
  summarize(rho = cor(Richness,orig_depth,method = "spearman"),
            Mean = round(mean(Richness),2),
            SD = round(sqrt(var(Richness)),3),
            Singletons = round(mean(singleton),2)) |> 
  mutate(Mean = paste0(Mean," (",SD,")"),.after = rho,.keep = "unused") |> xtable::xtable()



ggplot(aes(x=orig_depth,y=Richness,group=rf_level,color=rf_level),data=df %>% filter(study_id== "A") %>% filter(prior == "No prior")) +
  geom_point(alpha=.3) +
  geom_smooth(method="lm") +
  facet_wrap(~pooled,scales="free") +
  theme_minimal() +
  scale_colour_brewer(palette = "Set1",name="Sub-sampling strategy",labels=c("None","Post-ASV","Pre-ASV")) +
  theme(axis.title.x = element_blank()) +
  ylab("Observed richness") -> p_A



ggplot(aes(x=orig_depth,y=Richness,group=rf_level,color=rf_level),data=df %>% filter(study_id== "B") %>% filter(prior == "No prior")) +
  geom_point(alpha=.3) +
  geom_smooth(method="lm") +
  facet_wrap(~pooled,scales="free") +
  theme_minimal() +
  scale_colour_brewer(palette = "Set1",name="Sub-sampling strategy",labels=c("None","Post-ASV","Pre-ASV")) +
  theme(axis.title.x = element_blank()) +
  ylab("Observed richness") -> p_B


ggplot(aes(x=orig_depth,y=Richness,group=rf_level,color=rf_level),data=df %>% filter(study_id== "D") %>% filter(prior == "No prior")) +
  geom_point(alpha=.3) +
  geom_smooth(method="lm") +
  facet_wrap(~pooled,scales="free") +
  theme_minimal() +
  scale_colour_brewer(palette = "Set1",name="Sub-sampling strategy",labels=c("None","Post-ASV","Pre-ASV")) +
  xlab("Original sequencing depth") +
  ylab("Observed richness") -> p_D



ggarrange(p_A,p_B,p_D,ncol=1,common.legend = T,labels = c("A","B","C")) -> p

ggsave("Plots/corr_allstudies.PDF",plot = p,device = cairo_pdf,units = "in",width = 8.3,height=11.7)




ggplot(aes(x=orig_depth,y=Richness,group=rf_level,color=rf_level),data=df %>% filter(study_id== "A") %>% filter(prior == "Prior")) +
  geom_point(alpha=.3) +
  geom_smooth(method="lm") +
  facet_wrap(~pooled,scales="free") +
  theme_minimal() +
  scale_colour_brewer(palette = "Set1",name="Sub-sampling strategy",labels=c("None","Post-ASV","Pre-ASV")) +
  theme(axis.title.x = element_blank()) +
  ylab("Observed richness") -> p_A



ggplot(aes(x=orig_depth,y=Richness,group=rf_level,color=rf_level),data=df %>% filter(study_id== "B") %>% filter(prior == "Prior")) +
  geom_point(alpha=.3) +
  geom_smooth(method="lm") +
  facet_wrap(~pooled,scales="free") +
  theme_minimal() +
  scale_colour_brewer(palette = "Set1",name="Sub-sampling strategy",labels=c("None","Post-ASV","Pre-ASV")) +
  theme(axis.title.x = element_blank()) +
  ylab("Observed richness") -> p_B


ggplot(aes(x=orig_depth,y=Richness,group=rf_level,color=rf_level),data=df %>% filter(study_id== "D") %>% filter(prior == "Prior")) +
  geom_point(alpha=.3) +
  geom_smooth(method="lm") +
  geom_smooth(method="lm") +
  facet_wrap(~pooled,scales="free") +
  theme_minimal() +
  scale_colour_brewer(palette = "Set1",name="Sub-sampling strategy",labels=c("None","Post-ASV","Pre-ASV")) +
  xlab("Original sequencing depth") +
  ylab("Observed richness") -> p_D



ggarrange(p_A,p_B,p_D,ncol=1,common.legend = T,labels = c("A","B","C")) -> p

ggsave("Plots/corr_allstudies_prior.PDF",plot = p,device = cairo_pdf,units = "in",width = 8.3,height=11.7)



### Rarefaction curves

library(iNEXT)
library(phyloseq)

#### Load data
A_unpooled <- readRDS("Data/A/results/unpooled_noprior_A.RDS")[[1]] 
B_unpooled <- readRDS("Data/B/results/unpooled_noprior_B.RDS")[[1]] 
D_unpooled <- readRDS("Data/D/results/unpooled_noprior_D.RDS")[[1]] 

A_unpooled <- subset_samples(A_unpooled,sample_sums(A_unpooled)>10000) %>% otu_table %>% as.data.frame
B_unpooled <- subset_samples(B_unpooled,sample_sums(B_unpooled)>10000) %>% otu_table %>% as.data.frame
D_unpooled <- subset_samples(D_unpooled,sample_sums(D_unpooled)>10000) %>% otu_table %>% as.data.frame



A_pooled <- readRDS("Data/A/results/pooled_noprior_A.RDS")[[1]] 
B_pooled <- readRDS("Data/B/results/pooled_noprior_B.RDS")[[1]] 
D_pooled <- readRDS("Data/D/results/pooled_noprior_D.RDS")[[1]] 

A_pooled <- subset_samples(A_pooled,sample_sums(A_pooled)>10000) %>% otu_table %>% as.data.frame
B_pooled <- subset_samples(B_pooled,sample_sums(B_pooled)>10000) %>% otu_table %>% as.data.frame
D_pooled <- subset_samples(D_pooled,sample_sums(D_pooled)>10000) %>% otu_table %>% as.data.frame



A_out_unpooled <- iNEXT(A_unpooled, q=c(0), datatype="abundance",endpoint = 20000)
A_out_pooled <- iNEXT(A_pooled, q=c(0), datatype="abundance",endpoint = 20000)

B_out_unpooled <- iNEXT(B_unpooled, q=c(0), datatype="abundance")
B_out_pooled <- iNEXT(B_pooled, q=c(0), datatype="abundance")

D_out_unpooled <- iNEXT(D_unpooled, q=c(0), datatype="abundance")
D_out_pooled <- iNEXT(D_pooled, q=c(0), datatype="abundance")



ggiNEXT(A_out_unpooled, type=1) + theme_minimal() + xlab("Sequencing depth") -> p_A_unpooled
ggiNEXT(A_out_pooled, type=1) + theme_minimal() + xlab("Sequencing depth") -> p_A_pooled


ggiNEXT(B_out_unpooled, type=1) + theme_minimal() + xlab("Sequencing depth") -> p_B_unpooled
ggiNEXT(B_out_pooled, type=1) + theme_minimal() + xlab("Sequencing depth") -> p_B_pooled

ggiNEXT(D_out_unpooled, type=1) + theme_minimal() + xlab("Sequencing depth") -> p_D_unpooled
ggiNEXT(D_out_pooled, type=1) + theme_minimal() + xlab("Sequencing depth") -> p_D_pooled


myRFplot <- function(p) {
  
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



myRFplot(p_A_unpooled) + theme(axis.title.x = element_blank(),axis.text.x = element_text(angle=-35)) -> p_A1
myRFplot(p_A_pooled) + theme(axis.text.x = element_text(angle=-35)) -> p_A2

myRFplot(p_B_unpooled) + theme(axis.title.x = element_blank(),axis.text.x = element_text(angle=-35)) -> p_A1 -> p_B1
myRFplot(p_B_pooled) + theme(axis.text.x = element_text(angle=-35)) -> p_B2

myRFplot(p_D_unpooled) + theme(axis.title.x = element_blank(),axis.text.x = element_text(angle=-35)) -> p_A1 -> p_D1
myRFplot(p_D_pooled) + theme(axis.text.x = element_text(angle=-35)) -> p_D2


p_A <- p_A1 / p_A2
p_B <- p_B1 / p_B2
p_D <- p_D1 / p_D2

p_A | p_B | p_D -> plot



plot[[1]] <- plot[[1]] + plot_layout(tag_level = "new")
plot[[2]] <- plot[[2]] + plot_layout(tag_level = "new")
plot[[3]] <- plot[[3]] + plot_layout(tag_level = "new")


plot + plot_annotation(tag_levels = c("A","1")) -> plot

ggsave("Plots/rarefy_curve_all.PDF",plot = plot,device = cairo_pdf,units = "in",width = 8.3,height=6)


############################
#### Graphical abstract ####
############################
df1 <- df

df1$pooled <- factor(df1$pooled, levels = c("Unpooled","Pooled"), 
                     labels = c("Unpooled processing DADA2","Pooled processing in DADA2"))

ggplot(aes(x=orig_depth,y=Richness,group=rf_level,color=rf_level),data=df1 %>% filter(study_id== "B") %>% filter(prior == "No prior")) +
  geom_point(alpha=.3) +
  geom_smooth(method="lm") +
  facet_wrap(~pooled,scales="free") +
  theme_minimal() +
  scale_colour_brewer(palette = "Set1",name="Sub-sampling strategy",labels=c("No sub-sampling","Sub-sampling after deriving ASVs","Sub-sampling before deriving ASVs")) +
  theme(axis.title.x = element_blank()) +
  ylab("Observed richness") -> p_abstract

ggsave("Plots/graphical_abstract.PDF",plot = p_abstract, device = cairo_pdf, units = "in",width=12,height = 6)
