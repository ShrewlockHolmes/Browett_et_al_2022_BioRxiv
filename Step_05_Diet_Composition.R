#### Libraries ####

library(stringr)
library(ggplot2)
library(phyloseq)
#library(FSA)
#library(microbiome)
library(vegan)
library(gridExtra)
library(knitr)
library(kableExtra)
#library(VennDiagram)
#library(spaa)
library(dplyr)
library(beepr)
library(gridExtra)

#### Load Data ####

# Browett Method 1
load("./browett_et_al_2021_assignment_restrictions/filter_method_1/taxonomically_assigned_filtered.RData")
dir.create("./browett_et_al_2021_assignment_restrictions/filter_method_1/betadiversity_output")
dir.path <- "./browett_et_al_2021_assignment_restrictions/filter_method_1/betadiversity_output/"
# Browett Method 2
load("./browett_et_al_2021_assignment_restrictions/filter_method_2/taxonomically_assigned_filtered.RData")
dir.create("./browett_et_al_2021_assignment_restrictions/filter_method_2/betadiversity_output")
dir.path <- "./browett_et_al_2021_assignment_restrictions/filter_method_2/betadiversity_output/"

# Alberdi Method 1
load("./alberdi_taxa_assignment_restrictions/filter_method_1/taxonomically_assigned_filtered.RData")
dir.create("./alberdi_taxa_assignment_restrictions/filter_method_2/betadiversity_output")
dir.path <- "./alberdi_taxa_assignment_restrictions/filter_method_2/betadiversity_output/"

# Alberdi Method 2
load("./alberdi_taxa_assignment_restrictions/filter_method_2/taxonomically_assigned_filtered.RData")
dir.create("./alberdi_taxa_assignment_restrictions/filter_method_2/betadiversity_output")
dir.path <- "./alberdi_taxa_assignment_restrictions/filter_method_2/betadiversity_output/"

#### Compare metric ####

dir.create("./comparing_metric")

# This is code taken from Deagle et al.

paths <- c(b.meth1 = "./browett_et_al_2021_assignment_restrictions/filter_method_1/taxonomically_assigned_filtered.RData",
           b.meth2 = "./browett_et_al_2021_assignment_restrictions/filter_method_2/taxonomically_assigned_filtered.RData",
           a.meth1 = "./alberdi_taxa_assignment_restrictions/filter_method_1/taxonomically_assigned_filtered.RData",
           a.meth2 = "./alberdi_taxa_assignment_restrictions/filter_method_2/taxonomically_assigned_filtered.RData")

for(i in 1:length(paths)){
  
  mt.nam <- names(paths[i])
  load(paths[i])
  
  
  mrg.ra <- subset_samples(final.diet, shrew_country == "Pygmy_Ireland")
  #mrg.ra <- final.diet
  mrg.ra <- prune_taxa(taxa_sums(mrg.ra) > 0, mrg.ra)
  mrg.ra <- tax_glom(mrg.ra, taxrank = "order")
  
  otudf <- as.data.frame(as.matrix(otu_table(mrg.ra)))
  otudf <- as.matrix(otu_table(mrg.ra))
  rownames(otudf) <- tax_table(mrg.ra)[,"order"]
  colnames(otudf) <- sample_data(mrg.ra)$shrew_country
  
  BBA_Dprop=prop.table(otudf,2)
  dim(BBA_Dprop)
  apply(BBA_Dprop,2,sum)
  
  PA_BBA_D= BBA_Dprop
  PA_BBA_D[PA_BBA_D>.0099999999]=1
  PA_BBA_D[PA_BBA_D<.01]=0
  apply(PA_BBA_D,2,sum)
  
  FOO_All_BBA_=apply(PA_BBA_D,1,sum)/length(PA_BBA_D[1,])
  Ord_BBA_=order(-FOO_All_BBA_) #
  
  P_FOO_BBA_=  apply(PA_BBA_D,1,sum)/sum(apply(PA_BBA_D,1,sum))
  P_FOO_BBA_
  sum(P_FOO_BBA_)
  
  P_wFOO_BBA_=apply(prop.table(PA_BBA_D,2),1,mean,na.rm=T)
  P_wFOO_BBA_
  sum(P_wFOO_BBA_)
  
  P_RRA_BBA_=apply(prop.table(as.matrix(otudf),2),1,mean)
  P_RRA_BBA_
  sum(P_RRA_BBA_)
  
  par(mfrow = c(2,2))
  #jpeg(paste("./comparing_metric/metric_comparison_", mt.nam[i], ".jpeg", sep = ""), width = 600, height = 450)
  plot(P_FOO_BBA_[Ord_BBA_][1:25],pch=15,col=2,type="b",ylim= c(0,0.3),ylab="", xlab="Prey Taxa",xlim=c(0,26))
  points(26, sum(P_FOO_BBA_[Ord_BBA_][26:length(P_FOO_BBA_)]),pch=15,col=2,cex=1.2 )
  
  points(P_wFOO_BBA_[Ord_BBA_][1:25],pch=16,col=3,type="b")
  points(26, sum(P_wFOO_BBA_[Ord_BBA_][26:length(P_wFOO_BBA_)]),pch=16,col=3,cex=1.2 )
  sum(P_wFOO_BBA_[Ord_BBA_][1:15])
  
  points(P_RRA_BBA_[Ord_BBA_][1:25],pch=17,type="b")
  points(26, sum(P_RRA_BBA_[Ord_BBA_][26:length(P_RRA_BBA_)]),pch=17,cex=1.2 )
  text(4, P_RRA_BBA_[Ord_BBA_][4], labels = names(P_RRA_BBA_[Ord_BBA_][4]), pos = 3)
  text(6, P_RRA_BBA_[Ord_BBA_][6], labels = names(P_RRA_BBA_[Ord_BBA_][6]), pos = 4)
  text(7, P_RRA_BBA_[Ord_BBA_][7], labels = names(P_RRA_BBA_[Ord_BBA_][7]), pos = 2)
  
  #legend('top',c('POO','wPOO','RRA'),pch=c(15,16,17),col=c(2,3,1))
  #legend('topright', "S. minutus - Ireland")
  
  
  # ---
  ### France
  
  mrg.ra <- subset_samples(final.diet, shrew_country == "GWTS_Ireland")
  #mrg.ra <- final.diet
  mrg.ra <- prune_taxa(taxa_sums(mrg.ra) > 0, mrg.ra)
  mrg.ra <- tax_glom(mrg.ra, taxrank = "order")
  
  otudf <- as.data.frame(as.matrix(otu_table(mrg.ra)))
  otudf <- as.matrix(otu_table(mrg.ra))
  rownames(otudf) <- tax_table(mrg.ra)[,"order"]
  colnames(otudf) <- sample_data(mrg.ra)$shrew_season
  
  BBA_Dprop=prop.table(otudf,2)
  dim(BBA_Dprop)
  apply(BBA_Dprop,2,sum)
  
  PA_BBA_D= BBA_Dprop
  PA_BBA_D[PA_BBA_D>.0099999999]=1
  PA_BBA_D[PA_BBA_D<.01]=0
  apply(PA_BBA_D,2,sum)
  
  FOO_All_BBA_=apply(PA_BBA_D,1,sum)/length(PA_BBA_D[1,])
  Ord_BBA_=order(-FOO_All_BBA_) #
  
  P_FOO_BBA_=  apply(PA_BBA_D,1,sum)/sum(apply(PA_BBA_D,1,sum))
  P_FOO_BBA_
  sum(P_FOO_BBA_)
  
  P_wFOO_BBA_=apply(prop.table(PA_BBA_D,2),1,mean,na.rm=T)
  P_wFOO_BBA_
  sum(P_wFOO_BBA_)
  
  P_RRA_BBA_=apply(prop.table(as.matrix(otudf),2),1,mean)
  P_RRA_BBA_
  sum(P_RRA_BBA_)
  
  plot(P_FOO_BBA_[Ord_BBA_][1:25],pch=15,col=2,type="b",ylim= c(0,0.5),ylab="", xlab="Prey Taxa",xlim=c(0,26))
  points(26, sum(P_FOO_BBA_[Ord_BBA_][26:length(P_FOO_BBA_)]),pch=15,col=2,cex=1.2 )
  
  points(P_wFOO_BBA_[Ord_BBA_][1:25],pch=16,col=3,type="b")
  points(26, sum(P_wFOO_BBA_[Ord_BBA_][26:length(P_wFOO_BBA_)]),pch=16,col=3,cex=1.2 )
  sum(P_wFOO_BBA_[Ord_BBA_][1:15])
  
  points(P_RRA_BBA_[Ord_BBA_][1:25],pch=17,type="b")
  points(26, sum(P_RRA_BBA_[Ord_BBA_][26:length(P_RRA_BBA_)]),pch=17,cex=1.2 )
  text(1, P_RRA_BBA_[Ord_BBA_][1], labels = names(P_RRA_BBA_[Ord_BBA_][1]), pos = 4)
  
  #legend('top',c('POO','wPOO','RRA'),pch=c(15,16,17),col=c(2,3,1))
  #legend('topright', "C. russula - Ireland")
  
  
  # ---
  ### Pygmies
  
  mrg.ra <- subset_samples(final.diet, shrew_country == "Pygmy_France")
  #mrg.ra <- final.diet
  mrg.ra <- prune_taxa(taxa_sums(mrg.ra) > 0, mrg.ra)
  mrg.ra <- tax_glom(mrg.ra, taxrank = "order")
  
  otudf <- as.data.frame(as.matrix(otu_table(mrg.ra)))
  otudf <- as.matrix(otu_table(mrg.ra))
  rownames(otudf) <- tax_table(mrg.ra)[,"order"]
  colnames(otudf) <- sample_data(mrg.ra)$shrew_season
  
  BBA_Dprop=prop.table(otudf,2)
  dim(BBA_Dprop)
  apply(BBA_Dprop,2,sum)
  
  PA_BBA_D= BBA_Dprop
  PA_BBA_D[PA_BBA_D>.0099999999]=1
  PA_BBA_D[PA_BBA_D<.01]=0
  apply(PA_BBA_D,2,sum)
  
  FOO_All_BBA_=apply(PA_BBA_D,1,sum)/length(PA_BBA_D[1,])
  Ord_BBA_=order(-FOO_All_BBA_) #
  
  P_FOO_BBA_=  apply(PA_BBA_D,1,sum)/sum(apply(PA_BBA_D,1,sum))
  P_FOO_BBA_
  sum(P_FOO_BBA_)
  
  P_wFOO_BBA_=apply(prop.table(PA_BBA_D,2),1,mean,na.rm=T)
  P_wFOO_BBA_
  sum(P_wFOO_BBA_)
  
  P_RRA_BBA_=apply(prop.table(as.matrix(otudf),2),1,mean)
  P_RRA_BBA_
  sum(P_RRA_BBA_)
  
  plot(P_FOO_BBA_[Ord_BBA_][1:25],pch=15,col=2,type="b",ylim= c(0,0.3),ylab="", xlab="Prey Taxa",xlim=c(0,26))
  points(26, sum(P_FOO_BBA_[Ord_BBA_][26:length(P_FOO_BBA_)]),pch=15,col=2,cex=1.2 )
  
  points(P_wFOO_BBA_[Ord_BBA_][1:25],pch=16,col=3,type="b")
  points(26, sum(P_wFOO_BBA_[Ord_BBA_][26:length(P_wFOO_BBA_)]),pch=16,col=3,cex=1.2 )
  sum(P_wFOO_BBA_[Ord_BBA_][1:15])
  
  points(P_RRA_BBA_[Ord_BBA_][1:25],pch=17,type="b")
  points(26, sum(P_RRA_BBA_[Ord_BBA_][26:length(P_RRA_BBA_)]),pch=17,cex=1.2 )
  text(1, P_RRA_BBA_[Ord_BBA_][1], labels = names(P_RRA_BBA_[Ord_BBA_][1]), pos = 4)
  #text(2, P_RRA_BBA_[Ord_BBA_][2], labels = names(P_RRA_BBA_[Ord_BBA_][2]), pos = 1)
  text(3, P_RRA_BBA_[Ord_BBA_][3], labels = names(P_RRA_BBA_[Ord_BBA_][3]), pos = 4)
  
  #legend('top',c('POO','wPOO','RRA'),pch=c(15,16,17),col=c(2,3,1))
  #legend('topright', "S.minutus - Belle Ile")
  
  # ---
  ### GWTS
  
  mrg.ra <- subset_samples(final.diet, shrew_country == "GWTS_France")
  #mrg.ra <- final.diet
  mrg.ra <- prune_taxa(taxa_sums(mrg.ra) > 0, mrg.ra)
  mrg.ra <- tax_glom(mrg.ra, taxrank = "order")
  
  otudf <- as.data.frame(as.matrix(otu_table(mrg.ra)))
  otudf <- as.matrix(otu_table(mrg.ra))
  rownames(otudf) <- tax_table(mrg.ra)[,"order"]
  colnames(otudf) <- sample_data(mrg.ra)$shrew_season
  
  BBA_Dprop=prop.table(otudf,2)
  dim(BBA_Dprop)
  apply(BBA_Dprop,2,sum)
  
  PA_BBA_D= BBA_Dprop
  PA_BBA_D[PA_BBA_D>.0099999999]=1
  PA_BBA_D[PA_BBA_D<.01]=0
  apply(PA_BBA_D,2,sum)
  
  FOO_All_BBA_=apply(PA_BBA_D,1,sum)/length(PA_BBA_D[1,])
  Ord_BBA_=order(-FOO_All_BBA_) #
  
  P_FOO_BBA_=  apply(PA_BBA_D,1,sum)/sum(apply(PA_BBA_D,1,sum))
  P_FOO_BBA_
  sum(P_FOO_BBA_)
  
  P_wFOO_BBA_=apply(prop.table(PA_BBA_D,2),1,mean,na.rm=T)
  P_wFOO_BBA_
  sum(P_wFOO_BBA_)
  
  P_RRA_BBA_=apply(prop.table(as.matrix(otudf),2),1,mean)
  P_RRA_BBA_
  sum(P_RRA_BBA_)
  
  plot(P_FOO_BBA_[Ord_BBA_][1:25],pch=15,col=2,type="b",ylim= c(0,0.5),ylab="", xlab="Prey Taxa",xlim=c(0,26))
  points(26, sum(P_FOO_BBA_[Ord_BBA_][26:length(P_FOO_BBA_)]),pch=15,col=2,cex=1.2 )
  #text(P_FOO_BBA_[Ord_BBA_][1], labels = names(P_FOO_BBA_[Ord_BBA_][1]))
  
  points(P_wFOO_BBA_[Ord_BBA_][1:25],pch=16,col=3,type="b")
  points(26, sum(P_wFOO_BBA_[Ord_BBA_][26:length(P_wFOO_BBA_)]),pch=16,col=3,cex=1.2 )
  sum(P_wFOO_BBA_[Ord_BBA_][1:15])
  
  points(P_RRA_BBA_[Ord_BBA_][1:25],pch=17,type="b")
  points(26, sum(P_RRA_BBA_[Ord_BBA_][26:length(P_RRA_BBA_)]),pch=17,cex=1.2 )
  text(1, P_RRA_BBA_[Ord_BBA_][1], labels = names(P_RRA_BBA_[Ord_BBA_][1]), pos = 4)
  text(3, P_RRA_BBA_[Ord_BBA_][3], labels = names(P_RRA_BBA_[Ord_BBA_][3]), pos = 4)
  text(4, P_RRA_BBA_[Ord_BBA_][4], labels = names(P_RRA_BBA_[Ord_BBA_][4]), pos = 2)
  
  #legend('top',c('POO','wPOO','RRA'),pch=c(15,16,17),col=c(2,3,1))
  #legend('topright', "C.russula - Belle Ile")
  #dev.off()
  
}

beep("coin")


#### colour palette ####

pal.o = c("#f0a3ff", # Araneae  # Spare palette from https://graphicdesign.stackexchange.com/questions/3682/where-can-i-find-a-large-palette-set-of-contrasting-colors-for-coloring-many-d
          "#0075dc", # Arch * # Palette is from study to create most contrasting colours
          "#993f00", # Coleoptera
          "#4c005c", # Diptera *
          "#191919", # Ento *
          "#005c31", # Geo *
          "#2bce48", # Glom
          "#ffcc99", # Hap *
          "#808080", # Hemip
          "#94ffb5", # Hymen
          "#8f7c00", # Isopoda
          "#9dcc00", # Julida
          "#c20088", # Lepidop
          "blue", # Lithobiomorpha *
          "#ffa405", # Mesostig
          "#ffa8bb", # Neuropter *
          "#426600", # Oppiliones *
          "#ff0010", # Orbitida *
          "#5ef1f2",# Orthoptera *
          "#00998f", # Poly *
          "#e0ff66", # psocoptera
          "indianred",     # Sarc
          "#003380",    # stylo
          "green", # trom *
          "khaki4",
          "darkred",
          "coral4",
          "violetred2",
          "#0075dc",
          "#993f00")


#### Composition (POO metric) ####
## Paths to files

paths <- c(b.meth1 = "./browett_et_al_2021_assignment_restrictions/filter_method_1/taxonomically_assigned_filtered.RData",
           b.meth2 = "./browett_et_al_2021_assignment_restrictions/filter_method_2/taxonomically_assigned_filtered.RData",
           a.meth1 = "./alberdi_taxa_assignment_restrictions/filter_method_1/taxonomically_assigned_filtered.RData",
           a.meth2 = "./alberdi_taxa_assignment_restrictions/filter_method_2/taxonomically_assigned_filtered.RData")

plot_list1 <- list()
#plot_list2 <- list()


for(i in 1:length(paths)){
  
  mt.nam <- names(paths[i])
  load(paths[i])
  
  # Agglomerate to order level
  
  mrg.ra <- tax_glom(final.diet, taxrank = "order")
  
  otudf <- as.data.frame(as.matrix(otu_table(mrg.ra)))
  otudf <- as.matrix(otu_table(mrg.ra))
  rownames(otudf) <- tax_table(mrg.ra)[,"order"]
  colnames(otudf) <- sample_data(mrg.ra)$shrew_zone
  
  BBA_Dprop=prop.table(otudf,2)
  dim(BBA_Dprop)
  apply(BBA_Dprop,2,sum)
  
  PA_BBA_D= BBA_Dprop
  PA_BBA_D[PA_BBA_D>.0099999999]=1
  PA_BBA_D[PA_BBA_D<.01]=0
  apply(PA_BBA_D,2,sum)
  
  FOO_All_BBA_=apply(PA_BBA_D,1,sum)/length(PA_BBA_D[1,]) * 100
  Ord_BBA_=order(-FOO_All_BBA_) #
  POO_All_BBA_pfr=apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "Pygmy_France")],1,sum)/sum(apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "Pygmy_France")],1,sum)) * 100
  POO_All_BBA_cfr=apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "GWTS_France")],1,sum)/sum(apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "GWTS_France")],1,sum)) * 100
  POO_All_BBA_ped=apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "Pygmy_Edge")],1,sum)/sum(apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "Pygmy_Edge")],1,sum)) * 100
  POO_All_BBA_ced=apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "GWTS_Edge")],1,sum)/sum(apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "GWTS_Edge")],1,sum)) * 100
  POO_All_BBA_pout=apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "Pygmy_Outside")],1,sum)/sum(apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "Pygmy_Outside")],1,sum)) * 100
  POO_All_BBA_cin=apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "GWTS_Inside")],1,sum)/sum(apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "GWTS_Inside")],1,sum)) * 100
  
  chck2 <- data.frame(order = rep(c(names(POO_All_BBA_pfr[1:19]), "Other"), 6),
                      prev = c(POO_All_BBA_pfr[1:19], sum(POO_All_BBA_pfr[20:ntaxa(mrg.ra)]),
                               POO_All_BBA_cfr[1:19], sum(POO_All_BBA_cfr[20:ntaxa(mrg.ra)]),
                               POO_All_BBA_ped[1:19], sum(POO_All_BBA_ped[20:ntaxa(mrg.ra)]),
                               POO_All_BBA_pout[1:19], sum(POO_All_BBA_pout[20:ntaxa(mrg.ra)]),
                               POO_All_BBA_ced[1:19], sum(POO_All_BBA_ced[20:ntaxa(mrg.ra)]),
                               POO_All_BBA_cin[1:19], sum(POO_All_BBA_cin[20:ntaxa(mrg.ra)])),
                      species = c(rep("Pygmy", 20), rep("GWTS", 20),rep("Pygmy", 40), rep("GWTS", 40)),
                      country = c(rep("Belle Ile", 40), rep("Ireland", 80)),
                      zone = c(rep("Pygmy", 20), rep("GWTS", 20), rep("Pygmy Edge", 20), rep("Pygmy Outside", 20),
                               rep("GWTS Edge", 20), rep("GWTS Inside", 20)))
  
  chck2$zone <- factor(chck2$zone, levels = c("GWTS", "Pygmy", "GWTS Inside", "GWTS Edge", "Pygmy Edge", "Pygmy Outside"))
  
  p.1 <- ggplot(chck2, aes(x= zone, y=prev, fill = order)) +
    geom_bar(stat="identity", color="black") +
    scale_fill_manual(values = pal.o) +
    facet_wrap(~ country, scale = "free_x") +
    theme_classic() + # appearance composition
    theme(strip.text.x = element_text(size = 15, face = "bold")) + # facet text
    scale_x_discrete(labels=c("GWTS" = "C.russula",
                              "Pygmy" = "S.minutus",
                              "GWTS Inside" = "C.russula\nInside",
                              "GWTS Edge" = "C.russula\nEdge",
                              "Pygmy Edge" = "S.minutus\nEdge",
                              "Pygmy Outside" = "S.minutus\nOutside")) +
    theme(legend.position = "right",
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          legend.key.size = unit(0.8, "cm")) +
    labs(y = "Percentage of Occurence") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 20,
                                      face = "bold")) +
    theme(axis.text.x = element_text(size = 10 ,
                                     face = "bold",
                                     color = "black",
                                     angle = 0,
                                     hjust = 0.5),
          axis.text.y = element_text(size = 12,
                                     face = "bold"))
  
  p.1
  
  
  plot_list1[[mt.nam]] <- p.1
  
}

beep("coin")

# Sample size in each group
table(sample_data(final.diet)$shrew_zone)

plot_list1[[1]]
plot_list1[[2]]
plot_list1[[3]]
plot_list1[[4]]


#### Composition (RRA metric) ####

paths <- c(b.meth1 = "./browett_et_al_2021_assignment_restrictions/filter_method_1/taxonomically_assigned_filtered.RData",
           b.meth2 = "./browett_et_al_2021_assignment_restrictions/filter_method_2/taxonomically_assigned_filtered.RData",
           a.meth1 = "./alberdi_taxa_assignment_restrictions/filter_method_1/taxonomically_assigned_filtered.RData",
           a.meth2 = "./alberdi_taxa_assignment_restrictions/filter_method_2/taxonomically_assigned_filtered.RData")

plot_list1 <- list() # Individuals
plot_list2 <- list() # Shrew_Country
plot_list3 <- list() # Shrew_Country_Season
plot_list4 <- list() # Shrew_Zone

for(i in 1:length(paths)){
  
  mt.nam <- names(paths[i])
  load(paths[i])
  
  
  top25ord <- final.diet
  
  # Section to make it match POO #
  # isolate taxa that we want to group as 'other'
  othr.ord <- subset_taxa(final.diet, !(order %in% c("Araneae", "Class_Arachnida", "Class_Insecta", "Class_Malacostraca", 
                                                     "Coleoptera", "Diptera", "Entomobryomorpha", "Geophilomorpha", "Glomerida", "Haplotaxida",
                                                     "Hemiptera", "Hymenoptera", "Isopoda", "Julida", "Lepidoptera", "Opiliones", "Polydesmida", 
                                                     "Sarcoptiformes", "Stylommatophora")))
  
  n <- taxa_names(othr.ord)
  tax_table(top25ord)[n, 1:4] <- "Other"
  
  top25ord <- tax_glom(top25ord, taxrank = "order")
  
  # Section to use top taxa according to RRA #
  
  sampl.filt <- final.diet
  
  order.glom <- tax_glom(sampl.filt, taxrank= "order")
  ntx <- as.numeric(ntaxa(order.glom))
  top25ord <- prune_taxa(names(sort(taxa_sums(order.glom),TRUE)[1:22]), order.glom)
  othr.ord <- prune_taxa(names(sort(taxa_sums(order.glom),TRUE)[23:ntx]), order.glom)
  
  n <- taxa_names(othr.ord)
  tax_table(order.glom)[n, 1:4] <- "Other"
  n <- grep("Class_", tax_table(order.glom)[,"order"])
  tax_table(order.glom)[n, 1:4] <- "Other"
  top25ord <- tax_glom(order.glom, taxrank= "order")
  
  #- Create Barplot -#
  
  ra.samples = transform_sample_counts(top25ord, function(x) 100 * x/sum(x))
  ra.samples.bar <- phyloseq::plot_bar(ra.samples) # extracts information needed for barplots
  ra.samples.bar.data <- ra.samples.bar$data
  
  ra.samples.bar.data$Country <- str_replace(ra.samples.bar.data$Country, "France", "Belle Ile")
  ra.samples.bar.data$Zone <- str_replace(ra.samples.bar.data$Zone, "France", "Belle Ile")
  ra.samples.bar.data$Shrew <- str_replace(ra.samples.bar.data$Shrew, "Pygmy", "S. minutus")
  ra.samples.bar.data$Shrew <- str_replace(ra.samples.bar.data$Shrew, "GWTS", "C. russula")
  ra.samples.bar.data$Zone <- factor(ra.samples.bar.data$Zone, 
                                     levels = c("Belle Ile", "Inside", "Outside", "Edge"))
  
  p1 <- ggplot(ra.samples.bar.data, aes(x= Sample, y=Abundance, fill = order)) +
    geom_bar(stat="identity", color="black") +
    scale_fill_manual(values = pal.o) +
    facet_wrap(~ Shrew + Zone, scale = "free_x") +
    theme_classic() + # appearance composition
    theme(strip.text.x = element_text(size = 15, face = "bold")) + # facet text
    theme(legend.position = "right",
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.key.size = unit(1, "cm")) +
    labs(y = "Relative Abundance (%)") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 20, face = "bold")) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(color = "black")) 
  
  
  #leg <- cowplot::get_legend(p1)
  #grid.newpage()
  #grid.draw(leg)
  
  #--
  
  merge_samples_mean <- function(physeq, group){
    group_sums <- as.matrix(table(sample_data(physeq)[ ,group]))[,1]
    merged <- merge_samples(physeq, group)
    x <- as.matrix(otu_table(merged))
    if(taxa_are_rows(merged)){ x<-t(x) }
    out <- t(x/group_sums)
    out <- otu_table(out, taxa_are_rows = TRUE)
    otu_table(merged) <- out
    return(merged)
  }
  
  #--
  
  top25ord.ra = transform_sample_counts(top25ord, function(x) 100 * x/sum(x))
  ra.samples = merge_samples_mean(top25ord.ra, "shrew_country")
  sample_data(ra.samples)$shrew_country <- rownames(sample_data(ra.samples))
  sample_data(ra.samples)$Shrew <- c("GWTS", "GWTS", "Pygmy", "Pygmy")
  sample_data(ra.samples)$Country <- rep(c("France", "Ireland"), 2)
  
  ra.samples.bar <- phyloseq::plot_bar(ra.samples) # extracts information needed for barplots
  ra.samples.bar.data <- ra.samples.bar$data
  p2 <- ggplot(ra.samples.bar.data, aes(x= Sample, y=Abundance, fill = order)) +
    geom_bar(stat="identity", color="black") +
    scale_fill_manual(values = pal.o) +
    facet_wrap(~ Country, scale = "free_x") +
    theme_classic() + # appearance composition
    theme(strip.text.x = element_text(size = 15, face = "bold")) + # facet text
    scale_x_discrete(labels=c("GWTS_France" = "C.russula",
                              "Pygmy_France" = "S.minutus",
                              "GWTS_Ireland" = "C.russula",
                              "Pygmy_Ireland" = "S.minutus")) +
    theme(legend.position = "right",
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          legend.key.size = unit(0.8, "cm")) +
    labs(y = "Relative Abundance (%)") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 20,
                                      face = "bold")) +
    theme(axis.text.x = element_text(size = 12 ,
                                     face = "bold",
                                     angle = 0,
                                     hjust = 0.5),
          axis.text.y = element_text(size = 12,
                                     face = "bold"))
  
  
  
  ra.samples <- merge_samples_mean(top25ord.ra, "shrew_country_season")
  sample_data(ra.samples)$shrew_country_season <- levels(sample_data(top25ord.ra)$shrew_country_season)
  sample_data(ra.samples)$Shrew <- c(rep("GWTS", 4), rep("Pygmy", 4))
  sample_data(ra.samples)$Country <- rep(c("France", "France", "Ireland", "Ireland"), 2)
  sample_data(ra.samples)$Season <- rep(c("Summer", "Winter"), 4)
  sample_data(ra.samples)$country_season <- paste(sample_data(ra.samples)$Country, 
                                                  sample_data(ra.samples)$Season,
                                                  sep = " ")
  
  
  ra.samples.bar <- phyloseq::plot_bar(ra.samples) # extracts information needed for barplots
  ra.samples.bar.data <- ra.samples.bar$data
  p3 <- ggplot(ra.samples.bar.data, aes(x= Sample, y=Abundance, fill = order)) +
    geom_bar(stat="identity", color="black") +
    scale_fill_manual(values = pal.o) +
    facet_wrap(~ Country + Shrew, scale = "free_x", nrow = 1) +
    theme_classic() + # appearance composition
    theme(strip.text.x = element_text(size = 15, face = "bold")) + # facet text
    scale_x_discrete(labels=c("GWTS_France_Summer" = "C.russula Sum",
                              "Pygmy_France_Summer" = "S.minutus Sum",
                              "GWTS_France_Winter" = "C.russla Win",
                              "Pygmy_France_Winter" = "S.minutus Win",
                              "GWTS_Ireland_Summer" = "C.russula Sum",
                              "Pygmy_Ireland_Summer" = "S.minutus Sum",
                              "GWTS_Ireland_Winter" = "C.russla Win",
                              "Pygmy_Ireland_Winter" = "S.minutus Win")) +
    theme(legend.position = "right",
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          legend.key.size = unit(0.8, "cm")) +
    labs(y = "Relative Abundance (%)") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 20,
                                      face = "bold")) +
    theme(axis.text.x = element_text(size = 12 ,
                                     face = "bold",
                                     angle = 45,
                                     hjust = 1),
          axis.text.y = element_text(size = 12,
                                     face = "bold"))
  
  
  ra.samples <- merge_samples_mean(top25ord.ra, "shrew_zone")
  sample_data(ra.samples)$shrew_zone <- levels(sample_data(top25ord.ra)$shrew_zone)
  sample_data(ra.samples)$Shrew <- c(rep("GWTS", 3), rep("Pygmy", 3))
  sample_data(ra.samples)$Country <- rep(c("Ireland", "France", "Ireland"), 2)
  sample_data(ra.samples)$Zone <- c("Edge", "France", "Inside", "Edge", "France", "Outside")
  
  
  ra.samples.bar <- phyloseq::plot_bar(ra.samples) # extracts information needed for barplots
  ra.samples.bar.data <- ra.samples.bar$data
  ra.samples.bar.data$Sample <- factor(ra.samples.bar.data$Sample, levels = c("GWTS_France", "Pygmy_France", "GWTS_Inside", "GWTS_Edge", "Pygmy_Edge", "Pygmy_Outside"))
  p4 <- ggplot(ra.samples.bar.data, aes(x= Sample, y=Abundance, fill = order)) +
    geom_bar(stat="identity", color="black") +
    scale_fill_manual(values = pal.o) +
    facet_wrap(~ Country, scale = "free_x", nrow = 1) +
    theme_classic() + # appearance composition
    theme(strip.text.x = element_text(size = 15, face = "bold")) + # facet text
    scale_x_discrete(labels=c("GWTS_France_Summer" = "C.russula Sum",
                              "Pygmy_France_Summer" = "S.minutus Sum",
                              "GWTS_France_Winter" = "C.russla Win",
                              "Pygmy_France_Winter" = "S.minutus Win",
                              "GWTS_Ireland_Summer" = "C.russula Sum",
                              "Pygmy_Ireland_Summer" = "S.minutus Sum",
                              "GWTS_Ireland_Winter" = "C.russla Win",
                              "Pygmy_Ireland_Winter" = "S.minutus Win")) +
    theme(legend.position = "right",
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          legend.key.size = unit(0.8, "cm")) +
    labs(y = "Relative Abundance (%)") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 20,
                                      face = "bold")) +
    theme(axis.text.x = element_text(size = 12 ,
                                     face = "bold",
                                     angle = 45,
                                     hjust = 1),
          axis.text.y = element_text(size = 12,
                                     face = "bold"))
  
  
  plot_list1[[mt.nam]] <- p1
  plot_list2[[mt.nam]] <- p2
  plot_list3[[mt.nam]] <- p3
  plot_list4[[mt.nam]] <- p4
  
  
}

beep("coin")

plot_list1[[4]]
plot_list2[[1]]
plot_list3[[1]]
plot_list4[[4]]

#### Taxa tables for break down of different orders ####

merge_samples_mean <- function(physeq, group){
  group_sums <- as.matrix(table(sample_data(physeq)[ ,group]))[,1]
  merged <- merge_samples(physeq, group)
  x <- as.matrix(otu_table(merged))
  if(taxa_are_rows(merged)){ x<-t(x) }
  out <- t(x/group_sums)
  out <- otu_table(out, taxa_are_rows = TRUE)
  otu_table(merged) <- out
  return(merged)
}

#---------------#

final.diet.ra = transform_sample_counts(final.diet, function(x) 100 * x/sum(x))
zone.samp <- merge_samples_mean(final.diet.ra, "shrew_zone")
sample_data(zone.samp)$shrew_zone <- levels(sample_data(final.diet.ra)$shrew_zone)
sample_data(zone.samp)$Shrew <- c(rep("GWTS", 3), rep("Pygmy", 3))
sample_data(zone.samp)$Country <- rep(c("Ireland", "France", "Ireland"), 2)
sample_data(zone.samp)$Zone <- c("Edge", "France", "Inside", "Edge", "France", "Outside")



# Coleoptera

sbphy <- subset_taxa(zone.samp, (order %in% c("Coleoptera")))

tx.col <- as.data.frame(as.matrix(tax_table(sbphy)[,1:7]))
tx.col <- cbind(tx.col, as.matrix(otu_table(sbphy)))
write.csv(tx.col, file = "coleoptera.csv")

# Araneae

sbphy <- subset_taxa(zone.samp, (order %in% c("Araneae")))

tx.col <- as.data.frame(as.matrix(tax_table(sbphy)[,1:7]))
tx.col <- cbind(tx.col, as.matrix(otu_table(sbphy)))
write.csv(tx.col, file = "araneae.csv")


# Isopoda

sbphy <- subset_taxa(zone.samp, (order %in% c("Isopoda")))

tx.col <- as.data.frame(as.matrix(tax_table(sbphy)[,1:7]))
tx.col <- cbind(tx.col, as.matrix(otu_table(sbphy)))
write.csv(tx.col, file = "isopoda.csv")


# Lepidoptera

sbphy <- subset_taxa(zone.samp, (order %in% c("Lepidoptera")))

tx.col <- as.data.frame(as.matrix(tax_table(sbphy)[,1:7]))
tx.col <- cbind(tx.col, as.matrix(otu_table(sbphy)))
write.csv(tx.col, file = "lepidoptera.csv")

# Stylommatophora

sbphy <- subset_taxa(zone.samp, (order %in% c("Stylommatophora")))

tx.col <- as.data.frame(as.matrix(tax_table(sbphy)[,1:7]))
tx.col <- cbind(tx.col, as.matrix(otu_table(sbphy)))
write.csv(tx.col, file = "Stylommatophora.csv")


# Entomobryomorpha

sbphy <- subset_taxa(zone.samp, (order %in% c("Entomobryomorpha")))

tx.col <- as.data.frame(as.matrix(tax_table(sbphy)[,1:7]))
tx.col <- cbind(tx.col, as.matrix(otu_table(sbphy)))
write.csv(tx.col, file = "Entomobryomorpha.csv")



#### Ireland - Seasons #####

mrg.ra <- subset_samples(final.diet, Country == "Ireland")
#mrg.ra <- final.diet
mrg.ra <- prune_taxa(taxa_sums(mrg.ra) > 0, mrg.ra)
mrg.ra <- tax_glom(mrg.ra, taxrank = "order")

otudf <- as.data.frame(as.matrix(otu_table(mrg.ra)))
otudf <- as.matrix(otu_table(mrg.ra))
rownames(otudf) <- tax_table(mrg.ra)[,"order"]
colnames(otudf) <- sample_data(mrg.ra)$shrew_season

BBA_Dprop=prop.table(otudf,2)
dim(BBA_Dprop)
apply(BBA_Dprop,2,sum)

PA_BBA_D= BBA_Dprop
PA_BBA_D[PA_BBA_D>.0099999999]=1
PA_BBA_D[PA_BBA_D<.01]=0
apply(PA_BBA_D,2,sum)

FOO_All_BBA_=apply(PA_BBA_D,1,sum)/length(PA_BBA_D[1,]) * 100
Ord_BBA_=order(-FOO_All_BBA_) #
#FOO_All_BBA_psum=apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "Pygmy_Summer")],1,sum)/length(PA_BBA_D[1,which(colnames(PA_BBA_D) == "Pygmy_Summer")]) * 100
POO_All_BBA_psum=apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "Pygmy_Summer")],1,sum)/sum(apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "Pygmy_Summer")],1,sum)) * 100
#FOO_All_BBA_csum=apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "GWTS_Summer")],1,sum)/length(PA_BBA_D[1,which(colnames(PA_BBA_D) == "GWTS_Summer")]) * 100
POO_All_BBA_csum=apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "GWTS_Summer")],1,sum)/sum(apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "GWTS_Summer")],1,sum)) * 100
#FOO_All_BBA_pwin=apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "Pygmy_Winter")],1,sum)/length(PA_BBA_D[1,which(colnames(PA_BBA_D) == "Pygmy_Winter")]) * 100
POO_All_BBA_pwin=apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "Pygmy_Winter")],1,sum)/sum(apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "Pygmy_Winter")],1,sum)) * 100
#FOO_All_BBA_cwin=apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "GWTS_Winter")],1,sum)/length(PA_BBA_D[1,which(colnames(PA_BBA_D) == "GWTS_Winter")]) * 100
POO_All_BBA_cwin=apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "GWTS_Winter")],1,sum)/sum(apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "GWTS_Winter")],1,sum)) * 100

#chck <- data.frame(order = rep(names(FOO_All_BBA_psum[1:20]), 4),
#                   prev = c(FOO_All_BBA_psum[1:20], FOO_All_BBA_csum[1:20],
#                            FOO_All_BBA_pwin[1:20], FOO_All_BBA_cwin[1:20]),
#                   country = c(rep("PygmyS", 20), rep("GWTSS", 20),rep("PygmyW", 20), rep("GWTSW", 20)))

chck <- data.frame(order = rep(names(POO_All_BBA_psum[1:20]), 4),
                   prev = c(POO_All_BBA_psum[1:20], POO_All_BBA_pwin[1:20],
                            POO_All_BBA_csum[1:20], POO_All_BBA_cwin[1:20]),
                   country = c(rep("PygmyS", 20), rep("PygmyW", 20),rep("GWTSS", 20), rep("GWTSW", 20)))


#chck$order <- factor(chck$order[1:20], levels = chck$order[order(-chck$prev[1:20])])
chck$country <- factor(chck$country, levels = c("PygmyS", "PygmyW", "GWTSS", "GWTSW"))
chck$order <- factor(chck$order, levels = c("Diptera", "Entomobryomorpha", "Coleoptera",
                                            "Araneae", "Hemiptera", "Opiliones",
                                            "Class_Malacostraca", "Isopoda", "Class_Insecta",
                                            "Stylommatophora", "Lepidoptera", "Hymenoptera",
                                            "Class_Arachnida", "Archaeognatha", "Geophilomorpha",
                                            "Lithobiomorpha", "Mesostigmata", "Haplotaxida",
                                            "Sarcoptiformes", "Glomerida"))

p11 <- ggplot(chck, aes(x = order, y= prev)) +
  geom_bar(aes(fill = country), colour = "black", stat = "identity", position = "dodge") +
  theme_classic() +
  scale_fill_manual(values = c("lightblue", "blue", "pink", "red"),
                    #labels = c("S.minutus - Summer", "C.russula - Summer", "S.minutus - Winter", "C.russula - Winter"),
                    labels = c("S.minutus - Summer", "S.minutus - Winter", "C.russula - Summer", "C.russula - Winter"),
                    name = expression(paste(bold("A) - "), bold("Ireland")))) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0,25)) +
  theme(axis.text.x = element_text(size = 9, angle = 45, hjust = 1,
                                   face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold")) +
  labs(y = "Percentage of Occurence") +
  #theme(legend.position = c(0.15,0.8)) +
  theme(legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 10)) +
  #geom_text(x = 16, y =90, label = "France", size = 8) +
  theme(legend.position = c(0.75, 0.75)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10, face = "bold"))

p11

#### France - Seasons #####

mrg.ra <- subset_samples(final.diet, Country == "France")
#mrg.ra <- final.diet
mrg.ra <- prune_taxa(taxa_sums(mrg.ra) > 0, mrg.ra)
mrg.ra <- tax_glom(mrg.ra, taxrank = "order")

otudf <- as.data.frame(as.matrix(otu_table(mrg.ra)))
otudf <- as.matrix(otu_table(mrg.ra))
rownames(otudf) <- tax_table(mrg.ra)[,"order"]
colnames(otudf) <- sample_data(mrg.ra)$shrew_season

BBA_Dprop=prop.table(otudf,2)
dim(BBA_Dprop)
apply(BBA_Dprop,2,sum)

PA_BBA_D= BBA_Dprop
PA_BBA_D[PA_BBA_D>.0099999999]=1
PA_BBA_D[PA_BBA_D<.01]=0
apply(PA_BBA_D,2,sum)

FOO_All_BBA_=apply(PA_BBA_D,1,sum)/length(PA_BBA_D[1,]) * 100
Ord_BBA_=order(-FOO_All_BBA_) #
#FOO_All_BBA_psum=apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "Pygmy_Summer")],1,sum)/length(PA_BBA_D[1,which(colnames(PA_BBA_D) == "Pygmy_Summer")]) * 100
POO_All_BBA_psum=apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "Pygmy_Summer")],1,sum)/sum(apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "Pygmy_Summer")],1,sum)) * 100
#FOO_All_BBA_csum=apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "GWTS_Summer")],1,sum)/length(PA_BBA_D[1,which(colnames(PA_BBA_D) == "GWTS_Summer")]) * 100
POO_All_BBA_csum=apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "GWTS_Summer")],1,sum)/sum(apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "GWTS_Summer")],1,sum)) * 100
#FOO_All_BBA_pwin=apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "Pygmy_Winter")],1,sum)/length(PA_BBA_D[1,which(colnames(PA_BBA_D) == "Pygmy_Winter")]) * 100
POO_All_BBA_pwin=apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "Pygmy_Winter")],1,sum)/sum(apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "Pygmy_Winter")],1,sum)) * 100
#FOO_All_BBA_cwin=apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "GWTS_Winter")],1,sum)/length(PA_BBA_D[1,which(colnames(PA_BBA_D) == "GWTS_Winter")]) * 100
POO_All_BBA_cwin=apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "GWTS_Winter")],1,sum)/sum(apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == "GWTS_Winter")],1,sum)) * 100

#chck <- data.frame(order = rep(names(FOO_All_BBA_psum[1:20]), 4),
#                   prev = c(FOO_All_BBA_psum[1:20], FOO_All_BBA_csum[1:20],
#                            FOO_All_BBA_pwin[1:20], FOO_All_BBA_cwin[1:20]),
#                   country = c(rep("PygmyS", 20), rep("GWTSS", 20),rep("PygmyW", 20), rep("GWTSW", 20)))

chck <- data.frame(order = rep(names(POO_All_BBA_psum[1:20]), 4),
                   prev = c(POO_All_BBA_psum[1:20], POO_All_BBA_pwin[1:20],
                            POO_All_BBA_csum[1:20], POO_All_BBA_cwin[1:20]),
                   country = c(rep("PygmyS", 20), rep("PygmyW", 20),rep("GWTSS", 20), rep("GWTSW", 20)))


#chck$order <- factor(chck$order[1:20], levels = chck$order[order(-chck$prev[1:20])])
chck$country <- factor(chck$country, levels = c("PygmyS", "PygmyW", "GWTSS", "GWTSW"))
chck$order <- factor(chck$order, levels = c("Diptera", "Entomobryomorpha", "Coleoptera",
                                            "Araneae", "Hemiptera", "Opiliones",
                                            "Class_Malacostraca", "Isopoda", "Class_Insecta",
                                            "Stylommatophora", "Lepidoptera", "Hymenoptera",
                                            "Class_Arachnida", "Julida", "Geophilomorpha",
                                            "Trombidiformes", "Polydesmida", "Haplotaxida",
                                            "Sarcoptiformes", "Class_Diplopoda"))
p12 <- ggplot(chck, aes(x = order, y= prev)) +
  geom_bar(aes(fill = country), colour = "black", stat = "identity", position = "dodge") +
  theme_classic() +
  scale_fill_manual(values = c("lightblue", "blue", "pink", "red"),
                    #labels = c("S.minutus - Summer", "C.russula - Summer", "S.minutus - Winter", "C.russula - Winter"),
                    labels = c("S.minutus - Summer", "S.minutus - Winter", "C.russula - Summer", "C.russula - Winter"),
                    name = expression(paste(bold("B) - "), bold("Belle Ile")))) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0,25)) +
  theme(axis.text.x = element_text(size = 9, angle = 45, hjust = 1,
                                   face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold")) +
  labs(y = "Percentage of Occurence") +
  #theme(legend.position = c(0.15,0.8)) +
  theme(legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_blank(),
        legend.key.size = unit(0, "mm")) +
  #geom_text(x = 16, y =90, label = "France", size = 8) +
  theme(legend.position = c(0.74, 0.85)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10, face = "bold"))

p12

#### Grid Figures ####

grid.arrange(p11, p12, nrow = 2)

#### NMDS ####

lv <- c("species", "genus", "family", "order")

for(i in 1:length(paths)){
  
  mt.nam <- names(paths[i])
  load(paths[i])
  
  # RRA Bray Curtis - MOTU level
  
  dir.create("./NMDS")
  path1 <- paste("./NMDS/", mt.nam, "_NMDS_glom_", "MOTU", "_Bray_Curtis", sep = "")
  dir.create(path1)
  
  mrg.ra = transform_sample_counts(final.diet, function(x) 100 * x/sum(x))
  set.seed(738)
  ord_diet <- phyloseq::ordinate(mrg.ra,
                                 method = "NMDS",
                                 distance = "bray",
                                 k = 6)
  
  par.tab <- data.frame(plot = c("All_Bray", "Ireland_Bray"),
                        converged = NA,
                        dimensions = as.numeric(1:2),
                        stress = as.numeric(1:2))
  
  par.tab[1,2] <- ord_diet$converged # Did iterations converge?
  par.tab[1,3] <- ord_diet$ndim # number of dimensions
  par.tab[1,4] <- ord_diet$stress # stress
  
  scrs <- scores(ord_diet, display='sites')
  smpl.df <- as.data.frame(as.matrix(sample_data(mrg.ra)))
  scrs <- cbind(as.data.frame(scrs), smpl.df)
  
  # First we have to identify the 'outter most' samples to define the polygon shape.
  
  grp.a <- scrs[scrs$shrew_country == "GWTS_France", ][chull(scrs[scrs$shrew_country == "GWTS_France", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
  grp.b <- scrs[scrs$shrew_country == "GWTS_Ireland", ][chull(scrs[scrs$shrew_country == "GWTS_Ireland", c("NMDS1", "NMDS2")]), ]
  grp.c <- scrs[scrs$shrew_country == "Pygmy_France", ][chull(scrs[scrs$shrew_country == "Pygmy_France", c("NMDS1", "NMDS2")]), ]
  grp.d <- scrs[scrs$shrew_country == "Pygmy_Ireland", ][chull(scrs[scrs$shrew_country == "Pygmy_Ireland", c("NMDS1", "NMDS2")]), ]
  
  hull.data <- rbind(grp.a, grp.b, grp.c, grp.d)  #combine grp.a and grp.b
  hull.data
  
  p.1.2 <- ggplot(scrs, aes(x = NMDS1, y = NMDS2)) +
    geom_polygon(data = hull.data, aes(x=NMDS1, y=NMDS2, fill = shrew_country, group = shrew_country), alpha=0.30) +
    scale_fill_manual(values= c("#E7298A", "#D95F02", "#7570B3", "#1B9E77")) +
    geom_point(data = scrs, aes(colour = shrew_country), size=3, alpha=0.9) +
    scale_colour_manual(values= c("GWTS_France" = "#E7298A",
                                  "GWTS_Ireland" = "#D95F02",
                                  "Pygmy_France" = "#7570B3",
                                  "Pygmy_Ireland" = "#1B9E77"),
                        labels = c("GWTS_France" = "C.russula\nBelle Ile\n",
                                   "GWTS_Ireland" = "C.russula\nIreland\n",
                                   "Pygmy_France" = "S.minutus\nBelle Ile\n",
                                   "Pygmy_Ireland" = "S.minutus\nIreland"),
                        name = element_blank()) +
    theme_bw() +
    theme(axis.text = element_text(size = 14, face= "bold", colour = "black"),
          axis.title=element_text(size=18),
          legend.title= element_text(size= 18),
          legend.text = element_text(size = 14),
          legend.key = element_rect(fill = "white", colour = NA),
          legend.background = element_rect(fill = "white"),
          panel.grid.major = element_line(colour = "grey95"),
          panel.grid.minor = element_line(colour = "grey95")
    ) +
    guides(colour = guide_legend(override.aes = list(size=5)),
           colour = guide_legend("shrew_country"), fill = FALSE) +
    theme(legend.position = "right") 
  
  #p.1.2 + geom_label(label = "K = 7\nStress < 0.1", 
  #                   x = -0.9, y=0.9,
  #                   label.size = 0.35)
  
  ggsave(paste(path1, "/All_samples_NMDS_1_2_", "MOTU", "_level.jpeg", sep = ""), 
         plot = p.1.2, device = "jpeg")
  
  grp.a <- scrs[scrs$shrew_country == "GWTS_France", ][chull(scrs[scrs$shrew_country == "GWTS_France", c("NMDS1", "NMDS3")]), ]  # hull values for grp A
  grp.b <- scrs[scrs$shrew_country == "GWTS_Ireland", ][chull(scrs[scrs$shrew_country == "GWTS_Ireland", c("NMDS1", "NMDS3")]), ]
  grp.c <- scrs[scrs$shrew_country == "Pygmy_France", ][chull(scrs[scrs$shrew_country == "Pygmy_France", c("NMDS1", "NMDS3")]), ]
  grp.d <- scrs[scrs$shrew_country == "Pygmy_Ireland", ][chull(scrs[scrs$shrew_country == "Pygmy_Ireland", c("NMDS1", "NMDS3")]), ]
  
  hull.data <- rbind(grp.a, grp.b, grp.c, grp.d)  #combine grp.a and grp.b
  hull.data
  
  
  p.1.3 <- ggplot(scrs, aes(x = NMDS1, y = NMDS3)) +
    geom_polygon(data = hull.data, aes(x=NMDS1, y=NMDS3, fill = shrew_country, group = shrew_country), alpha=0.30) +
    scale_fill_manual(values= c("#E7298A", "#D95F02", "#7570B3", "#1B9E77")) +
    geom_point(data = scrs, aes(colour = shrew_country), size=3, alpha=0.9) +
    scale_colour_manual(values= c("GWTS_France" = "#E7298A",
                                  "GWTS_Ireland" = "#D95F02",
                                  "Pygmy_France" = "#7570B3",
                                  "Pygmy_Ireland" = "#1B9E77"),
                        labels = c("GWTS_France" = "C.russula\nBelle Ile",
                                   "GWTS_Ireland" = "C.russula\nIreland",
                                   "Pygmy_France" = "S.minutus\nBelle Ile",
                                   "Pygmy_Ireland" = "S.minutus\nIreland"),
                        name = element_blank()) +
    theme_bw() +
    theme(axis.text = element_text(size = 14, face= "bold", colour = "black"),
          axis.title=element_text(size=18),
          legend.title= element_text(size= 18),
          legend.text = element_text(size = 14),
          legend.key = element_rect(fill = "white", colour = NA),
          legend.background = element_rect(fill = "white"),
          panel.grid.major = element_line(colour = "grey95"),
          panel.grid.minor = element_line(colour = "grey95")
    ) +
    guides(colour = guide_legend(override.aes = list(size=5)),
           colour = guide_legend("shrew_country"), fill = FALSE) +
    theme(legend.position = "right")
  
  p.1.3
  
  ggsave(paste(path1, "/All_samples_NMDS_1_3_", "MOTU", "_level.jpeg", sep = ""), 
         plot = p.1.3, device = "jpeg")
  
  ## Ireland
  
  final.ire <- subset_samples(final.diet, Country == "Ireland")
  #final.ire <- subset_samples(final.ire, alt.sample != "")
  mrg.ra = transform_sample_counts(final.ire, function(x) 100 * x/sum(x))
  
  set.seed(738)
  ord_diet <- phyloseq::ordinate(mrg.ra,
                                 method = "NMDS",
                                 distance = "bray",
                                 k = 5)
  
  
  par.tab[2,2] <- ord_diet$converged # Did iterations converge?
  par.tab[2,3] <- ord_diet$ndim # number of dimensions
  par.tab[2,4] <- ord_diet$stress # stress
  
  
  
  scrs <- scores(ord_diet, display='sites')
  smpl.df <- as.data.frame(as.matrix(sample_data(mrg.ra)))
  scrsire <- cbind(as.data.frame(scrs), smpl.df)
  
  grp.a <- scrsire[scrsire$shrew_zone == "GWTS_Edge", ][chull(scrsire[scrsire$shrew_zone == "GWTS_Edge", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
  grp.b <- scrsire[scrsire$shrew_zone == "GWTS_Inside", ][chull(scrsire[scrsire$shrew_zone == "GWTS_Inside", c("NMDS1", "NMDS2")]), ]
  grp.c <- scrsire[scrsire$shrew_zone == "Pygmy_Edge", ][chull(scrsire[scrsire$shrew_zone == "Pygmy_Edge", c("NMDS1", "NMDS2")]), ]
  grp.d <- scrsire[scrsire$shrew_zone == "Pygmy_Outside", ][chull(scrsire[scrsire$shrew_zone == "Pygmy_Outside", c("NMDS1", "NMDS2")]), ]
  
  hull.data <- rbind(grp.a, grp.b, grp.c, grp.d)  #combine grp.a and grp.b
  hull.data
  
  
  p.1.4 <- ggplot(scrsire, aes(x = NMDS1, y = NMDS2)) +
    geom_polygon(data = hull.data, aes(x=NMDS1, y=NMDS2, fill = shrew_zone, group = shrew_zone), alpha=0.30) +
    scale_fill_manual(values= c("#66A61E", "#E6AB02", "#A6761D", "#666666")) +
    geom_point(data = scrsire, aes(colour = shrew_zone), size=3, alpha=0.9) +
    scale_colour_manual(values= c("GWTS_Edge" = "#66A61E",
                                  "GWTS_Inside" = "#E6AB02",
                                  "Pygmy_Edge" = "#A6761D",
                                  "Pygmy_Outside" = "#666666"),
                        labels = c("GWTS_Inside" = "C.russula\nInside\n",
                                   "GWTS_Edge" = "C.russula\nEdge\n",
                                   "Pygmy_Edge" = "S.minutus\nEdge\n",
                                   "Pygmy_Outside" = "S.minutus\nOutside"),
                        name = element_blank()) +
    theme_bw() +
    theme(axis.text = element_text(size = 14, face= "bold", colour = "black"),
          axis.title=element_text(size=18),
          legend.title= element_text(size= 18),
          legend.text = element_text(size = 14),
          legend.key = element_rect(fill = "white", colour = NA),
          legend.background = element_rect(fill = "white"),
          panel.grid.major = element_line(colour = "grey95"),
          panel.grid.minor = element_line(colour = "grey95")
    ) +
    guides(colour = guide_legend(override.aes = list(size=5)),
           colour = guide_legend("shrew_zone"), fill = FALSE) +
    theme(legend.position = "right")
  
  #p.1.4 + geom_label(label = "K = 7\nStress < 0.09", 
  #                   x = -0.9, y = -0.7,
  #                   label.size = 0.35)
  
  ggsave(paste(path1, "/Ireland_samples_NMDS_1_2_", "MOTU", "_level.jpeg", sep = ""), 
         plot = p.1.4, device = "jpeg")
  
  ## Ireland NMDS 1 and NMDS 3
  
  grp.a <- scrsire[scrsire$shrew_zone == "GWTS_Edge", ][chull(scrsire[scrsire$shrew_zone == "GWTS_Edge", c("NMDS1", "NMDS3")]), ]  # hull values for grp A
  grp.b <- scrsire[scrsire$shrew_zone == "GWTS_Inside", ][chull(scrsire[scrsire$shrew_zone == "GWTS_Inside", c("NMDS1", "NMDS3")]), ]
  grp.c <- scrsire[scrsire$shrew_zone == "Pygmy_Edge", ][chull(scrsire[scrsire$shrew_zone == "Pygmy_Edge", c("NMDS1", "NMDS3")]), ]
  grp.d <- scrsire[scrsire$shrew_zone == "Pygmy_Outside", ][chull(scrsire[scrsire$shrew_zone == "Pygmy_Outside", c("NMDS1", "NMDS3")]), ]
  
  hull.data <- rbind(grp.a, grp.b, grp.c, grp.d)  #combine grp.a and grp.b
  hull.data
  
  
  p.1.5 <- ggplot(scrsire, aes(x = NMDS1, y = NMDS3)) +
    geom_polygon(data = hull.data, aes(x=NMDS1, y=NMDS3, fill = shrew_zone, group = shrew_zone), alpha=0.30) +
    scale_fill_manual(values= c("#66A61E", "#E6AB02", "#A6761D", "#666666")) +
    geom_point(data = scrsire, aes(colour = shrew_zone), size=3, alpha=0.9) +
    scale_colour_manual(values= c("GWTS_Edge" = "#66A61E",
                                  "GWTS_Inside" = "#E6AB02",
                                  "Pygmy_Edge" = "#A6761D",
                                  "Pygmy_Outside" = "#666666"),
                        labels = c("GWTS_Inside" = "C.russula\nInside",
                                   "GWTS_Edge" = "C.russula\nEdge",
                                   "Pygmy_Edge" = "S.minutus\nEdge",
                                   "Pygmy_Outside" = "S.minutus\nOutside"),
                        name = element_blank()) +
    theme_bw() +
    theme(axis.text = element_text(size = 14, face= "bold", colour = "black"),
          axis.title=element_text(size=18),
          legend.title= element_text(size= 18),
          legend.text = element_text(size = 14),
          legend.key = element_rect(fill = "white", colour = NA),
          legend.background = element_rect(fill = "white"),
          panel.grid.major = element_line(colour = "grey95"),
          panel.grid.minor = element_line(colour = "grey95")
    ) +
    guides(colour = guide_legend(override.aes = list(size=5)),
           colour = guide_legend("shrew_zone"), fill = FALSE) +
    theme(legend.position = "right")
  
  #p.1.5
  
  ggsave(paste(path1, "/Ireland_samples_NMDS_1_3_", "MOTU", "_level.jpeg", sep = ""), 
         plot = p.1.5, device = "jpeg")
  
  write.csv(par.tab, file = paste(path1, "/parameter_table_bray_c_", "MOTU", "_level.csv", sep = ""))
  
  #### PA Jacard ####
  
  path2 <- paste("./NMDS/", mt.nam, "_NMDS_glom_", "MOTU", "_Jaccard", sep = "")
  dir.create(path2)
  
  mrg.ra = transform_sample_counts(final.diet, function(x) 100 * x/sum(x))
  o.tab <- otu_table(mrg.ra)
  o.tab[1:10,1:10]
  o.tab[o.tab>.0099999999]=1
  o.tab[o.tab<.01]=0
  o.tab[1:10,1:10]
  otu_table(mrg.ra) <- o.tab
  otu_table(mrg.ra)[1:10,1:10]
  
  set.seed(738)
  ord_diet <- phyloseq::ordinate(mrg.ra,
                                 method = "NMDS",
                                 distance = "jaccard",
                                 k = 5)
  
  par.tab <- data.frame(plot = c("All_Jaccard", "Ireland_Jaccard"),
                        converged = NA,
                        dimensions = as.numeric(1:2),
                        stress = as.numeric(1:2))
  
  par.tab[1,2] <- ord_diet$converged # Did iterations converge?
  par.tab[1,3] <- ord_diet$ndim # number of dimensions
  par.tab[1,4] <- ord_diet$stress # stress
  
  scrs <- scores(ord_diet, display='sites')
  smpl.df <- as.data.frame(as.matrix(sample_data(mrg.ra)))
  scrs <- cbind(as.data.frame(scrs), smpl.df)
  
  # First we have to identify the 'outter most' samples to define the polygon shape.
  
  grp.a <- scrs[scrs$shrew_country == "GWTS_France", ][chull(scrs[scrs$shrew_country == "GWTS_France", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
  grp.b <- scrs[scrs$shrew_country == "GWTS_Ireland", ][chull(scrs[scrs$shrew_country == "GWTS_Ireland", c("NMDS1", "NMDS2")]), ]
  grp.c <- scrs[scrs$shrew_country == "Pygmy_France", ][chull(scrs[scrs$shrew_country == "Pygmy_France", c("NMDS1", "NMDS2")]), ]
  grp.d <- scrs[scrs$shrew_country == "Pygmy_Ireland", ][chull(scrs[scrs$shrew_country == "Pygmy_Ireland", c("NMDS1", "NMDS2")]), ]
  
  hull.data <- rbind(grp.a, grp.b, grp.c, grp.d)  #combine grp.a and grp.b
  hull.data
  
  
  p.2.1 <- ggplot(scrs, aes(x = NMDS1, y = NMDS2)) +
    geom_polygon(data = hull.data, aes(x=NMDS1, y=NMDS2, fill = shrew_country, group = shrew_country), alpha=0.30) +
    scale_fill_manual(values= c("#E7298A", "#D95F02", "#7570B3", "#1B9E77")) +
    geom_point(data = scrs, aes(colour = shrew_country), size=3, alpha=0.9) +
    scale_colour_manual(values= c("GWTS_France" = "#E7298A",
                                  "GWTS_Ireland" = "#D95F02",
                                  "Pygmy_France" = "#7570B3",
                                  "Pygmy_Ireland" = "#1B9E77"),
                        labels = c("GWTS_France" = "C.russula\nBelle Ile\n",
                                   "GWTS_Ireland" = "C.russula\nIreland\n",
                                   "Pygmy_France" = "S.minutus\nBelle Ile\n",
                                   "Pygmy_Ireland" = "S.minutus\nIreland"),
                        name = element_blank()) +
    theme_bw() +
    theme(axis.text = element_text(size = 14, face= "bold", colour = "black"),
          axis.title=element_text(size=18),
          legend.title= element_text(size= 18),
          legend.text = element_text(size = 14),
          legend.key = element_rect(fill = "white", colour = NA),
          legend.background = element_rect(fill = "white"),
          panel.grid.major = element_line(colour = "grey95"),
          panel.grid.minor = element_line(colour = "grey95")
    ) +
    guides(colour = guide_legend(override.aes = list(size=5)),
           colour = guide_legend("shrew_country"), fill = FALSE) +
    theme(legend.position = "right") 
  
  #p.2.1 + geom_label(label = "K = 7\nStress < 0.1", 
  #                   x = -0.9, y=0.9,
  #                   label.size = 0.35)
  
  ggsave(paste(path2, "/All_samples_NMDS_1_2_", "MOTU", "_level.jpeg", sep = ""), 
         plot = p.2.1, device = "jpeg")
  
  grp.a <- scrs[scrs$shrew_country == "GWTS_France", ][chull(scrs[scrs$shrew_country == "GWTS_France", c("NMDS1", "NMDS3")]), ]  # hull values for grp A
  grp.b <- scrs[scrs$shrew_country == "GWTS_Ireland", ][chull(scrs[scrs$shrew_country == "GWTS_Ireland", c("NMDS1", "NMDS3")]), ]
  grp.c <- scrs[scrs$shrew_country == "Pygmy_France", ][chull(scrs[scrs$shrew_country == "Pygmy_France", c("NMDS1", "NMDS3")]), ]
  grp.d <- scrs[scrs$shrew_country == "Pygmy_Ireland", ][chull(scrs[scrs$shrew_country == "Pygmy_Ireland", c("NMDS1", "NMDS3")]), ]
  
  hull.data <- rbind(grp.a, grp.b, grp.c, grp.d)  #combine grp.a and grp.b
  hull.data
  
  
  p.2.2 <- ggplot(scrs, aes(x = NMDS1, y = NMDS3)) +
    geom_polygon(data = hull.data, aes(x=NMDS1, y=NMDS3, fill = shrew_country, group = shrew_country), alpha=0.30) +
    scale_fill_manual(values= c("#E7298A", "#D95F02", "#7570B3", "#1B9E77")) +
    geom_point(data = scrs, aes(colour = shrew_country), size=3, alpha=0.9) +
    scale_colour_manual(values= c("GWTS_France" = "#E7298A",
                                  "GWTS_Ireland" = "#D95F02",
                                  "Pygmy_France" = "#7570B3",
                                  "Pygmy_Ireland" = "#1B9E77"),
                        labels = c("GWTS_France" = "C.russula\nBelle Ile",
                                   "GWTS_Ireland" = "C.russula\nIreland",
                                   "Pygmy_France" = "S.minutus\nBelle Ile",
                                   "Pygmy_Ireland" = "S.minutus\nIreland"),
                        name = element_blank()) +
    theme_bw() +
    theme(axis.text = element_text(size = 14, face= "bold", colour = "black"),
          axis.title=element_text(size=18),
          legend.title= element_text(size= 18),
          legend.text = element_text(size = 14),
          legend.key = element_rect(fill = "white", colour = NA),
          legend.background = element_rect(fill = "white"),
          panel.grid.major = element_line(colour = "grey95"),
          panel.grid.minor = element_line(colour = "grey95")
    ) +
    guides(colour = guide_legend(override.aes = list(size=5)),
           colour = guide_legend("shrew_country"), fill = FALSE) +
    theme(legend.position = "right")
  
  #p.2.2
  
  ggsave(paste(path2, "/All_samples_NMDS_1_3_", "MOTU", "_level.jpeg", sep = ""), 
         plot = p.2.2, device = "jpeg")
  
  ## Ireland
  
  final.ire <- subset_samples(final.diet, Country == "Ireland")
  #final.ire <- subset_samples(final.ire, alt.sample != "")
  mrg.ra = transform_sample_counts(final.ire, function(x) 100 * x/sum(x))
  o.tab <- otu_table(mrg.ra)
  o.tab[1:10,1:10]
  o.tab[o.tab>.0099999999]=1
  o.tab[o.tab<.01]=0
  o.tab[1:10,1:10]
  otu_table(mrg.ra) <- o.tab
  otu_table(mrg.ra)[1:10,1:10]
  
  set.seed(738)
  ord_diet <- phyloseq::ordinate(mrg.ra,
                                 method = "NMDS",
                                 distance = "jaccard",
                                 k = 5)
  
  par.tab[2,2] <- ord_diet$converged # Did iterations converge?
  par.tab[2,3] <- ord_diet$ndim # number of dimensions
  par.tab[2,4] <- ord_diet$stress # stress
  
  scrs <- scores(ord_diet, display='sites')
  smpl.df <- as.data.frame(as.matrix(sample_data(mrg.ra)))
  scrsire <- cbind(as.data.frame(scrs), smpl.df)
  
  grp.a <- scrsire[scrsire$shrew_zone == "GWTS_Edge", ][chull(scrsire[scrsire$shrew_zone == "GWTS_Edge", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
  grp.b <- scrsire[scrsire$shrew_zone == "GWTS_Inside", ][chull(scrsire[scrsire$shrew_zone == "GWTS_Inside", c("NMDS1", "NMDS2")]), ]
  grp.c <- scrsire[scrsire$shrew_zone == "Pygmy_Edge", ][chull(scrsire[scrsire$shrew_zone == "Pygmy_Edge", c("NMDS1", "NMDS2")]), ]
  grp.d <- scrsire[scrsire$shrew_zone == "Pygmy_Outside", ][chull(scrsire[scrsire$shrew_zone == "Pygmy_Outside", c("NMDS1", "NMDS2")]), ]
  
  hull.data <- rbind(grp.a, grp.b, grp.c, grp.d)  #combine grp.a and grp.b
  hull.data
  
  
  p.2.3 <- ggplot(scrsire, aes(x = NMDS1, y = NMDS2)) +
    geom_polygon(data = hull.data, aes(x=NMDS1, y=NMDS2, fill = shrew_zone, group = shrew_zone), alpha=0.30) +
    scale_fill_manual(values= c("#66A61E", "#E6AB02", "#A6761D", "#666666")) +
    geom_point(data = scrsire, aes(colour = shrew_zone), size=3, alpha=0.9) +
    scale_colour_manual(values= c("GWTS_Edge" = "#66A61E",
                                  "GWTS_Inside" = "#E6AB02",
                                  "Pygmy_Edge" = "#A6761D",
                                  "Pygmy_Outside" = "#666666"),
                        labels = c("GWTS_Inside" = "C.russula\nInside\n",
                                   "GWTS_Edge" = "C.russula\nEdge\n",
                                   "Pygmy_Edge" = "S.minutus\nEdge\n",
                                   "Pygmy_Outside" = "S.minutus\nOutside"),
                        name = element_blank()) +
    theme_bw() +
    theme(axis.text = element_text(size = 14, face= "bold", colour = "black"),
          axis.title=element_text(size=18),
          legend.title= element_text(size= 18),
          legend.text = element_text(size = 14),
          legend.key = element_rect(fill = "white", colour = NA),
          legend.background = element_rect(fill = "white"),
          panel.grid.major = element_line(colour = "grey95"),
          panel.grid.minor = element_line(colour = "grey95")
    ) +
    guides(colour = guide_legend(override.aes = list(size=5)),
           colour = guide_legend("shrew_zone"), fill = FALSE) +
    theme(legend.position = "right")
  
  #p.2.3 + geom_label(label = "K = 7\nStress < 0.09", 
  #                   x = -0.9, y = -0.7,
  #                   label.size = 0.35)
  
  ggsave(paste(path2, "/Ireland_samples_NMDS_1_2_", "MOTU", "_level.jpeg", sep = ""), 
         plot = p.2.3, device = "jpeg")
  
  
  ## Ireland NMDS 1 and NMDS 3
  
  grp.a <- scrsire[scrsire$shrew_zone == "GWTS_Edge", ][chull(scrsire[scrsire$shrew_zone == "GWTS_Edge", c("NMDS1", "NMDS3")]), ]  # hull values for grp A
  grp.b <- scrsire[scrsire$shrew_zone == "GWTS_Inside", ][chull(scrsire[scrsire$shrew_zone == "GWTS_Inside", c("NMDS1", "NMDS3")]), ]
  grp.c <- scrsire[scrsire$shrew_zone == "Pygmy_Edge", ][chull(scrsire[scrsire$shrew_zone == "Pygmy_Edge", c("NMDS1", "NMDS3")]), ]
  grp.d <- scrsire[scrsire$shrew_zone == "Pygmy_Outside", ][chull(scrsire[scrsire$shrew_zone == "Pygmy_Outside", c("NMDS1", "NMDS3")]), ]
  
  hull.data <- rbind(grp.a, grp.b, grp.c, grp.d)  #combine grp.a and grp.b
  hull.data
  
  
  p.2.4 <- ggplot(scrsire, aes(x = NMDS1, y = NMDS3)) +
    geom_polygon(data = hull.data, aes(x=NMDS1, y=NMDS3, fill = shrew_zone, group = shrew_zone), alpha=0.30) +
    scale_fill_manual(values= c("#66A61E", "#E6AB02", "#A6761D", "#666666")) +
    geom_point(data = scrsire, aes(colour = shrew_zone), size=3, alpha=0.9) +
    scale_colour_manual(values= c("GWTS_Edge" = "#66A61E",
                                  "GWTS_Inside" = "#E6AB02",
                                  "Pygmy_Edge" = "#A6761D",
                                  "Pygmy_Outside" = "#666666"),
                        labels = c("GWTS_Inside" = "C.russula\nInside",
                                   "GWTS_Edge" = "C.russula\nEdge",
                                   "Pygmy_Edge" = "S.minutus\nEdge",
                                   "Pygmy_Outside" = "S.minutus\nOutside"),
                        name = element_blank()) +
    theme_bw() +
    theme(axis.text = element_text(size = 14, face= "bold", colour = "black"),
          axis.title=element_text(size=18),
          legend.title= element_text(size= 18),
          legend.text = element_text(size = 14),
          legend.key = element_rect(fill = "white", colour = NA),
          legend.background = element_rect(fill = "white"),
          panel.grid.major = element_line(colour = "grey95"),
          panel.grid.minor = element_line(colour = "grey95")
    ) +
    guides(colour = guide_legend(override.aes = list(size=5)),
           colour = guide_legend("shrew_zone"), fill = FALSE) +
    theme(legend.position = "right")
  
  #p.2.4
  
  ggsave(paste(path2, "/Ireland_samples_NMDS_1_3_", "MOTU", "_level.jpeg", sep = ""), 
         plot = p.2.4, device = "jpeg")
  
  write.csv(par.tab, file = paste(path2, "/parameter_table_jaccard_", "MOTU", "_level.csv", sep = ""))
  
  # --- ---#
  
  for(k in 1:length(lv)){
    
    glom.diet <- tax_glom(final.diet, taxrank = lv[k])
    
    # RRA Bray Curtis
    
    path3 <- paste("./NMDS/", mt.nam, "_NMDS_glom_", lv[k], "_Bray_Curtis", sep = "")
    dir.create(path3)
    
    mrg.ra = transform_sample_counts(glom.diet, function(x) 100 * x/sum(x))
    set.seed(738)
    ord_diet <- phyloseq::ordinate(mrg.ra,
                                   method = "NMDS",
                                   distance = "bray",
                                   k = 6)
    
    par.tab <- data.frame(plot = c("All_Bray", "Ireland_Bray"),
                          converged = NA,
                          dimensions = as.numeric(1:2),
                          stress = as.numeric(1:2))
    
    par.tab[1,2] <- ord_diet$converged # Did iterations converge?
    par.tab[1,3] <- ord_diet$ndim # number of dimensions
    par.tab[1,4] <- ord_diet$stress # stress
    
    scrs <- scores(ord_diet, display='sites')
    smpl.df <- as.data.frame(as.matrix(sample_data(mrg.ra)))
    scrs <- cbind(as.data.frame(scrs), smpl.df)
    
    # First we have to identify the 'outter most' samples to define the polygon shape.
    
    grp.a <- scrs[scrs$shrew_country == "GWTS_France", ][chull(scrs[scrs$shrew_country == "GWTS_France", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
    grp.b <- scrs[scrs$shrew_country == "GWTS_Ireland", ][chull(scrs[scrs$shrew_country == "GWTS_Ireland", c("NMDS1", "NMDS2")]), ]
    grp.c <- scrs[scrs$shrew_country == "Pygmy_France", ][chull(scrs[scrs$shrew_country == "Pygmy_France", c("NMDS1", "NMDS2")]), ]
    grp.d <- scrs[scrs$shrew_country == "Pygmy_Ireland", ][chull(scrs[scrs$shrew_country == "Pygmy_Ireland", c("NMDS1", "NMDS2")]), ]
    
    hull.data <- rbind(grp.a, grp.b, grp.c, grp.d)  #combine grp.a and grp.b
    hull.data
    
    
    p.1.2 <- ggplot(scrs, aes(x = NMDS1, y = NMDS2)) +
      geom_polygon(data = hull.data, aes(x=NMDS1, y=NMDS2, fill = shrew_country, group = shrew_country), alpha=0.30) +
      scale_fill_manual(values= c("#E7298A", "#D95F02", "#7570B3", "#1B9E77")) +
      geom_point(data = scrs, aes(colour = shrew_country), size=3, alpha=0.9) +
      scale_colour_manual(values= c("GWTS_France" = "#E7298A",
                                    "GWTS_Ireland" = "#D95F02",
                                    "Pygmy_France" = "#7570B3",
                                    "Pygmy_Ireland" = "#1B9E77"),
                          labels = c("GWTS_France" = "C.russula\nBelle Ile\n",
                                     "GWTS_Ireland" = "C.russula\nIreland\n",
                                     "Pygmy_France" = "S.minutus\nBelle Ile\n",
                                     "Pygmy_Ireland" = "S.minutus\nIreland"),
                          name = element_blank()) +
      theme_bw() +
      theme(axis.text = element_text(size = 14, face= "bold", colour = "black"),
            axis.title=element_text(size=18),
            legend.title= element_text(size= 18),
            legend.text = element_text(size = 14),
            legend.key = element_rect(fill = "white", colour = NA),
            legend.background = element_rect(fill = "white"),
            panel.grid.major = element_line(colour = "grey95"),
            panel.grid.minor = element_line(colour = "grey95")
      ) +
      guides(colour = guide_legend(override.aes = list(size=5)),
             colour = guide_legend("shrew_country"), fill = FALSE) +
      theme(legend.position = "right") 
    
    #p.1.2 + geom_label(label = "K = 7\nStress < 0.1", 
    #                   x = -0.9, y=0.9,
    #                   label.size = 0.35)
    
    ggsave(paste(path3, "/All_samples_NMDS_1_2_", lv[k], "_level.jpeg", sep = ""), 
           plot = p.1.2, device = "jpeg")
    
    grp.a <- scrs[scrs$shrew_country == "GWTS_France", ][chull(scrs[scrs$shrew_country == "GWTS_France", c("NMDS1", "NMDS3")]), ]  # hull values for grp A
    grp.b <- scrs[scrs$shrew_country == "GWTS_Ireland", ][chull(scrs[scrs$shrew_country == "GWTS_Ireland", c("NMDS1", "NMDS3")]), ]
    grp.c <- scrs[scrs$shrew_country == "Pygmy_France", ][chull(scrs[scrs$shrew_country == "Pygmy_France", c("NMDS1", "NMDS3")]), ]
    grp.d <- scrs[scrs$shrew_country == "Pygmy_Ireland", ][chull(scrs[scrs$shrew_country == "Pygmy_Ireland", c("NMDS1", "NMDS3")]), ]
    
    hull.data <- rbind(grp.a, grp.b, grp.c, grp.d)  #combine grp.a and grp.b
    hull.data
    
    
    p.1.3 <- ggplot(scrs, aes(x = NMDS1, y = NMDS3)) +
      geom_polygon(data = hull.data, aes(x=NMDS1, y=NMDS3, fill = shrew_country, group = shrew_country), alpha=0.30) +
      scale_fill_manual(values= c("#E7298A", "#D95F02", "#7570B3", "#1B9E77")) +
      geom_point(data = scrs, aes(colour = shrew_country), size=3, alpha=0.9) +
      scale_colour_manual(values= c("GWTS_France" = "#E7298A",
                                    "GWTS_Ireland" = "#D95F02",
                                    "Pygmy_France" = "#7570B3",
                                    "Pygmy_Ireland" = "#1B9E77"),
                          labels = c("GWTS_France" = "C.russula\nBelle Ile",
                                     "GWTS_Ireland" = "C.russula\nIreland",
                                     "Pygmy_France" = "S.minutus\nBelle Ile",
                                     "Pygmy_Ireland" = "S.minutus\nIreland"),
                          name = element_blank()) +
      theme_bw() +
      theme(axis.text = element_text(size = 14, face= "bold", colour = "black"),
            axis.title=element_text(size=18),
            legend.title= element_text(size= 18),
            legend.text = element_text(size = 14),
            legend.key = element_rect(fill = "white", colour = NA),
            legend.background = element_rect(fill = "white"),
            panel.grid.major = element_line(colour = "grey95"),
            panel.grid.minor = element_line(colour = "grey95")
      ) +
      guides(colour = guide_legend(override.aes = list(size=5)),
             colour = guide_legend("shrew_country"), fill = FALSE) +
      theme(legend.position = "right")
    
    #p.1.3
    
    ggsave(paste(path3, "/All_samples_NMDS_1_3_", lv[k], "_level.jpeg", sep = ""), 
           plot = p.1.3, device = "jpeg")
    
    
    ## Ireland
    
    final.ire <- subset_samples(glom.diet, Country == "Ireland")
    #final.ire <- subset_samples(final.ire, alt.sample != "")
    mrg.ra = transform_sample_counts(final.ire, function(x) 100 * x/sum(x))
    
    set.seed(738)
    ord_diet <- phyloseq::ordinate(mrg.ra,
                                   method = "NMDS",
                                   distance = "bray",
                                   k = 5)
    
    par.tab[2,2] <- ord_diet$converged # Did iterations converge?
    par.tab[2,3] <- ord_diet$ndim # number of dimensions
    par.tab[2,4] <- ord_diet$stress # stress
    
    scrs <- scores(ord_diet, display='sites')
    smpl.df <- as.data.frame(as.matrix(sample_data(mrg.ra)))
    scrsire <- cbind(as.data.frame(scrs), smpl.df)
    
    grp.a <- scrsire[scrsire$shrew_zone == "GWTS_Edge", ][chull(scrsire[scrsire$shrew_zone == "GWTS_Edge", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
    grp.b <- scrsire[scrsire$shrew_zone == "GWTS_Inside", ][chull(scrsire[scrsire$shrew_zone == "GWTS_Inside", c("NMDS1", "NMDS2")]), ]
    grp.c <- scrsire[scrsire$shrew_zone == "Pygmy_Edge", ][chull(scrsire[scrsire$shrew_zone == "Pygmy_Edge", c("NMDS1", "NMDS2")]), ]
    grp.d <- scrsire[scrsire$shrew_zone == "Pygmy_Outside", ][chull(scrsire[scrsire$shrew_zone == "Pygmy_Outside", c("NMDS1", "NMDS2")]), ]
    
    hull.data <- rbind(grp.a, grp.b, grp.c, grp.d)  #combine grp.a and grp.b
    hull.data
    
    
    p.1.4 <- ggplot(scrsire, aes(x = NMDS1, y = NMDS2)) +
      geom_polygon(data = hull.data, aes(x=NMDS1, y=NMDS2, fill = shrew_zone, group = shrew_zone), alpha=0.30) +
      scale_fill_manual(values= c("#66A61E", "#E6AB02", "#A6761D", "#666666")) +
      geom_point(data = scrsire, aes(colour = shrew_zone), size=3, alpha=0.9) +
      scale_colour_manual(values= c("GWTS_Edge" = "#66A61E",
                                    "GWTS_Inside" = "#E6AB02",
                                    "Pygmy_Edge" = "#A6761D",
                                    "Pygmy_Outside" = "#666666"),
                          labels = c("GWTS_Inside" = "C.russula\nInside\n",
                                     "GWTS_Edge" = "C.russula\nEdge\n",
                                     "Pygmy_Edge" = "S.minutus\nEdge\n",
                                     "Pygmy_Outside" = "S.minutus\nOutside"),
                          name = element_blank()) +
      theme_bw() +
      theme(axis.text = element_text(size = 14, face= "bold", colour = "black"),
            axis.title=element_text(size=18),
            legend.title= element_text(size= 18),
            legend.text = element_text(size = 14),
            legend.key = element_rect(fill = "white", colour = NA),
            legend.background = element_rect(fill = "white"),
            panel.grid.major = element_line(colour = "grey95"),
            panel.grid.minor = element_line(colour = "grey95")
      ) +
      guides(colour = guide_legend(override.aes = list(size=5)),
             colour = guide_legend("shrew_zone"), fill = FALSE) +
      theme(legend.position = "right")
    
    #p.1.4 + geom_label(label = "K = 7\nStress < 0.09", 
    #                   x = -0.9, y = -0.7,
    #                   label.size = 0.35)
    
    ggsave(paste(path3, "/Ireland_samples_NMDS_1_2_", lv[k], "_level.jpeg", sep = ""), 
           plot = p.1.4, device = "jpeg")
    
    
    ## Ireland
    
    grp.a <- scrsire[scrsire$shrew_zone == "GWTS_Edge", ][chull(scrsire[scrsire$shrew_zone == "GWTS_Edge", c("NMDS1", "NMDS3")]), ]  # hull values for grp A
    grp.b <- scrsire[scrsire$shrew_zone == "GWTS_Inside", ][chull(scrsire[scrsire$shrew_zone == "GWTS_Inside", c("NMDS1", "NMDS3")]), ]
    grp.c <- scrsire[scrsire$shrew_zone == "Pygmy_Edge", ][chull(scrsire[scrsire$shrew_zone == "Pygmy_Edge", c("NMDS1", "NMDS3")]), ]
    grp.d <- scrsire[scrsire$shrew_zone == "Pygmy_Outside", ][chull(scrsire[scrsire$shrew_zone == "Pygmy_Outside", c("NMDS1", "NMDS3")]), ]
    
    hull.data <- rbind(grp.a, grp.b, grp.c, grp.d)  #combine grp.a and grp.b
    hull.data
    
    
    p.1.5 <- ggplot(scrsire, aes(x = NMDS1, y = NMDS3)) +
      geom_polygon(data = hull.data, aes(x=NMDS1, y=NMDS3, fill = shrew_zone, group = shrew_zone), alpha=0.30) +
      scale_fill_manual(values= c("#66A61E", "#E6AB02", "#A6761D", "#666666")) +
      geom_point(data = scrsire, aes(colour = shrew_zone), size=3, alpha=0.9) +
      scale_colour_manual(values= c("GWTS_Edge" = "#66A61E",
                                    "GWTS_Inside" = "#E6AB02",
                                    "Pygmy_Edge" = "#A6761D",
                                    "Pygmy_Outside" = "#666666"),
                          labels = c("GWTS_Inside" = "C.russula\nInside",
                                     "GWTS_Edge" = "C.russula\nEdge",
                                     "Pygmy_Edge" = "S.minutus\nEdge",
                                     "Pygmy_Outside" = "S.minutus\nOutside"),
                          name = element_blank()) +
      theme_bw() +
      theme(axis.text = element_text(size = 14, face= "bold", colour = "black"),
            axis.title=element_text(size=18),
            legend.title= element_text(size= 18),
            legend.text = element_text(size = 14),
            legend.key = element_rect(fill = "white", colour = NA),
            legend.background = element_rect(fill = "white"),
            panel.grid.major = element_line(colour = "grey95"),
            panel.grid.minor = element_line(colour = "grey95")
      ) +
      guides(colour = guide_legend(override.aes = list(size=5)),
             colour = guide_legend("shrew_zone"), fill = FALSE) +
      theme(legend.position = "right")
    
    p.1.5
    
    ggsave(paste(path3, "/Ireland_samples_NMDS_1_3_", lv[k], "_level.jpeg", sep = ""), 
           plot = p.1.5, device = "jpeg")
    
    write.csv(par.tab, file = paste(path3, "/parameter_table_bray_curtis_", lv[k], "_level.csv", sep = ""))
    
    #### PA Jacard ####
    
    path4 <- paste("./NMDS/", mt.nam, "_NMDS_glom_", lv[k], "_Jaccard", sep = "")
    dir.create(path4)
    
    mrg.ra = transform_sample_counts(glom.diet, function(x) 100 * x/sum(x))
    o.tab <- otu_table(mrg.ra)
    o.tab[1:10,1:10]
    o.tab[o.tab>.0099999999]=1
    o.tab[o.tab<.01]=0
    o.tab[1:10,1:10]
    otu_table(mrg.ra) <- o.tab
    otu_table(mrg.ra)[1:10,1:10]
    
    set.seed(738)
    ord_diet <- phyloseq::ordinate(mrg.ra,
                                   method = "NMDS",
                                   distance = "jaccard",
                                   k = 5)
    
    par.tab <- data.frame(plot = c("All_Jaccard", "Ireland_Jaccard"),
                          converged = NA,
                          dimensions = as.numeric(1:2),
                          stress = as.numeric(1:2))
    
    par.tab[1,2] <- ord_diet$converged # Did iterations converge?
    par.tab[1,3] <- ord_diet$ndim # number of dimensions
    par.tab[1,4] <- ord_diet$stress # stress
    
    scrs <- scores(ord_diet, display='sites')
    smpl.df <- as.data.frame(as.matrix(sample_data(mrg.ra)))
    scrs <- cbind(as.data.frame(scrs), smpl.df)
    
    # First we have to identify the 'outter most' samples to define the polygon shape.
    
    grp.a <- scrs[scrs$shrew_country == "GWTS_France", ][chull(scrs[scrs$shrew_country == "GWTS_France", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
    grp.b <- scrs[scrs$shrew_country == "GWTS_Ireland", ][chull(scrs[scrs$shrew_country == "GWTS_Ireland", c("NMDS1", "NMDS2")]), ]
    grp.c <- scrs[scrs$shrew_country == "Pygmy_France", ][chull(scrs[scrs$shrew_country == "Pygmy_France", c("NMDS1", "NMDS2")]), ]
    grp.d <- scrs[scrs$shrew_country == "Pygmy_Ireland", ][chull(scrs[scrs$shrew_country == "Pygmy_Ireland", c("NMDS1", "NMDS2")]), ]
    
    hull.data <- rbind(grp.a, grp.b, grp.c, grp.d)  #combine grp.a and grp.b
    hull.data
    
    p.2.2 <- ggplot(scrs, aes(x = NMDS1, y = NMDS2)) +
      geom_polygon(data = hull.data, aes(x=NMDS1, y=NMDS2, fill = shrew_country, group = shrew_country), alpha=0.30) +
      scale_fill_manual(values= c("#E7298A", "#D95F02", "#7570B3", "#1B9E77")) +
      geom_point(data = scrs, aes(colour = shrew_country), size=3, alpha=0.9) +
      scale_colour_manual(values= c("GWTS_France" = "#E7298A",
                                    "GWTS_Ireland" = "#D95F02",
                                    "Pygmy_France" = "#7570B3",
                                    "Pygmy_Ireland" = "#1B9E77"),
                          labels = c("GWTS_France" = "C.russula\nBelle Ile\n",
                                     "GWTS_Ireland" = "C.russula\nIreland\n",
                                     "Pygmy_France" = "S.minutus\nBelle Ile\n",
                                     "Pygmy_Ireland" = "S.minutus\nIreland"),
                          name = element_blank()) +
      theme_bw() +
      theme(axis.text = element_text(size = 14, face= "bold", colour = "black"),
            axis.title=element_text(size=18),
            legend.title= element_text(size= 18),
            legend.text = element_text(size = 14),
            legend.key = element_rect(fill = "white", colour = NA),
            legend.background = element_rect(fill = "white"),
            panel.grid.major = element_line(colour = "grey95"),
            panel.grid.minor = element_line(colour = "grey95")
      ) +
      guides(colour = guide_legend(override.aes = list(size=5)),
             colour = guide_legend("shrew_country"), fill = FALSE) +
      theme(legend.position = "right") 
    
    #p.2.2 + geom_label(label = "K = 7\nStress < 0.1", 
    #                   x = -0.9, y=0.9,
    #                   label.size = 0.35)
    
    ggsave(paste(path4, "/All_samples_NMDS_1_2_", lv[k], "_level.jpeg", sep = ""), 
           plot = p.2.2, device = "jpeg")
    
    
    grp.a <- scrs[scrs$shrew_country == "GWTS_France", ][chull(scrs[scrs$shrew_country == "GWTS_France", c("NMDS1", "NMDS3")]), ]  # hull values for grp A
    grp.b <- scrs[scrs$shrew_country == "GWTS_Ireland", ][chull(scrs[scrs$shrew_country == "GWTS_Ireland", c("NMDS1", "NMDS3")]), ]
    grp.c <- scrs[scrs$shrew_country == "Pygmy_France", ][chull(scrs[scrs$shrew_country == "Pygmy_France", c("NMDS1", "NMDS3")]), ]
    grp.d <- scrs[scrs$shrew_country == "Pygmy_Ireland", ][chull(scrs[scrs$shrew_country == "Pygmy_Ireland", c("NMDS1", "NMDS3")]), ]
    
    hull.data <- rbind(grp.a, grp.b, grp.c, grp.d)  #combine grp.a and grp.b
    hull.data
    
    
    p.2.3 <- ggplot(scrs, aes(x = NMDS1, y = NMDS3)) +
      geom_polygon(data = hull.data, aes(x=NMDS1, y=NMDS3, fill = shrew_country, group = shrew_country), alpha=0.30) +
      scale_fill_manual(values= c("#E7298A", "#D95F02", "#7570B3", "#1B9E77")) +
      geom_point(data = scrs, aes(colour = shrew_country), size=3, alpha=0.9) +
      scale_colour_manual(values= c("GWTS_France" = "#E7298A",
                                    "GWTS_Ireland" = "#D95F02",
                                    "Pygmy_France" = "#7570B3",
                                    "Pygmy_Ireland" = "#1B9E77"),
                          labels = c("GWTS_France" = "C.russula\nBelle Ile",
                                     "GWTS_Ireland" = "C.russula\nIreland",
                                     "Pygmy_France" = "S.minutus\nBelle Ile",
                                     "Pygmy_Ireland" = "S.minutus\nIreland"),
                          name = element_blank()) +
      theme_bw() +
      theme(axis.text = element_text(size = 14, face= "bold", colour = "black"),
            axis.title=element_text(size=18),
            legend.title= element_text(size= 18),
            legend.text = element_text(size = 14),
            legend.key = element_rect(fill = "white", colour = NA),
            legend.background = element_rect(fill = "white"),
            panel.grid.major = element_line(colour = "grey95"),
            panel.grid.minor = element_line(colour = "grey95")
      ) +
      guides(colour = guide_legend(override.aes = list(size=5)),
             colour = guide_legend("shrew_country"), fill = FALSE) +
      theme(legend.position = "right")
    
    #p.2.3
    
    ggsave(paste(path4, "/All_samples_NMDS_1_3_", lv[k], "_level.jpeg", sep = ""), 
           plot = p.2.3, device = "jpeg")
    
    
    ## Ireland
    
    final.ire <- subset_samples(glom.diet, Country == "Ireland")
    #final.ire <- subset_samples(final.ire, alt.sample != "") # option to remove 'trouble makers'
    mrg.ra = transform_sample_counts(final.ire, function(x) 100 * x/sum(x))
    o.tab <- otu_table(mrg.ra)
    o.tab[1:10,1:10]
    o.tab[o.tab>.0099999999]=1
    o.tab[o.tab<.01]=0
    o.tab[1:10,1:10]
    otu_table(mrg.ra) <- o.tab
    otu_table(mrg.ra)[1:10,1:10]
    
    set.seed(738)
    ord_diet <- phyloseq::ordinate(mrg.ra,
                                   method = "NMDS",
                                   distance = "jaccard",
                                   k = 5)
    
    par.tab[2,2] <- ord_diet$converged # Did iterations converge?
    par.tab[2,3] <- ord_diet$ndim # number of dimensions
    par.tab[2,4] <- ord_diet$stress # stress
    
    scrs <- scores(ord_diet, display='sites')
    smpl.df <- as.data.frame(as.matrix(sample_data(mrg.ra)))
    scrsire <- cbind(as.data.frame(scrs), smpl.df)
    
    grp.a <- scrsire[scrsire$shrew_zone == "GWTS_Edge", ][chull(scrsire[scrsire$shrew_zone == "GWTS_Edge", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
    grp.b <- scrsire[scrsire$shrew_zone == "GWTS_Inside", ][chull(scrsire[scrsire$shrew_zone == "GWTS_Inside", c("NMDS1", "NMDS2")]), ]
    grp.c <- scrsire[scrsire$shrew_zone == "Pygmy_Edge", ][chull(scrsire[scrsire$shrew_zone == "Pygmy_Edge", c("NMDS1", "NMDS2")]), ]
    grp.d <- scrsire[scrsire$shrew_zone == "Pygmy_Outside", ][chull(scrsire[scrsire$shrew_zone == "Pygmy_Outside", c("NMDS1", "NMDS2")]), ]
    
    hull.data <- rbind(grp.a, grp.b, grp.c, grp.d)  #combine grp.a and grp.b
    hull.data
    
    
    p.2.4 <- ggplot(scrsire, aes(x = NMDS1, y = NMDS2)) +
      geom_polygon(data = hull.data, aes(x=NMDS1, y=NMDS2, fill = shrew_zone, group = shrew_zone), alpha=0.30) +
      scale_fill_manual(values= c("#66A61E", "#E6AB02", "#A6761D", "#666666")) +
      geom_point(data = scrsire, aes(colour = shrew_zone), size=3, alpha=0.9) +
      scale_colour_manual(values= c("GWTS_Edge" = "#66A61E",
                                    "GWTS_Inside" = "#E6AB02",
                                    "Pygmy_Edge" = "#A6761D",
                                    "Pygmy_Outside" = "#666666"),
                          labels = c("GWTS_Inside" = "C.russula\nInside\n",
                                     "GWTS_Edge" = "C.russula\nEdge\n",
                                     "Pygmy_Edge" = "S.minutus\nEdge\n",
                                     "Pygmy_Outside" = "S.minutus\nOutside"),
                          name = element_blank()) +
      theme_bw() +
      theme(axis.text = element_text(size = 14, face= "bold", colour = "black"),
            axis.title=element_text(size=18),
            legend.title= element_text(size= 18),
            legend.text = element_text(size = 14),
            legend.key = element_rect(fill = "white", colour = NA),
            legend.background = element_rect(fill = "white"),
            panel.grid.major = element_line(colour = "grey95"),
            panel.grid.minor = element_line(colour = "grey95")
      ) +
      guides(colour = guide_legend(override.aes = list(size=5)),
             colour = guide_legend("shrew_zone"), fill = FALSE) +
      theme(legend.position = "right")
    
    #p.2.4 + geom_label(label = "K = 7\nStress < 0.09", 
    #                   x = -0.9, y = -0.7,
    #                   label.size = 0.35)
    
    ggsave(paste(path4, "/Ireland_samples_NMDS_1_2_", lv[k], "_level.jpeg", sep = ""), 
           plot = p.2.4, device = "jpeg")
    
    
    ## Ireland NMDS 1 and NMDS 3
    
    grp.a <- scrsire[scrsire$shrew_zone == "GWTS_Edge", ][chull(scrsire[scrsire$shrew_zone == "GWTS_Edge", c("NMDS1", "NMDS3")]), ]  # hull values for grp A
    grp.b <- scrsire[scrsire$shrew_zone == "GWTS_Inside", ][chull(scrsire[scrsire$shrew_zone == "GWTS_Inside", c("NMDS1", "NMDS3")]), ]
    grp.c <- scrsire[scrsire$shrew_zone == "Pygmy_Edge", ][chull(scrsire[scrsire$shrew_zone == "Pygmy_Edge", c("NMDS1", "NMDS3")]), ]
    grp.d <- scrsire[scrsire$shrew_zone == "Pygmy_Outside", ][chull(scrsire[scrsire$shrew_zone == "Pygmy_Outside", c("NMDS1", "NMDS3")]), ]
    
    hull.data <- rbind(grp.a, grp.b, grp.c, grp.d)  #combine grp.a and grp.b
    hull.data
    
    
    p.2.5 <- ggplot(scrsire, aes(x = NMDS1, y = NMDS3)) +
      geom_polygon(data = hull.data, aes(x=NMDS1, y=NMDS3, fill = shrew_zone, group = shrew_zone), alpha=0.30) +
      scale_fill_manual(values= c("#66A61E", "#E6AB02", "#A6761D", "#666666")) +
      geom_point(data = scrsire, aes(colour = shrew_zone), size=3, alpha=0.9) +
      scale_colour_manual(values= c("GWTS_Edge" = "#66A61E",
                                    "GWTS_Inside" = "#E6AB02",
                                    "Pygmy_Edge" = "#A6761D",
                                    "Pygmy_Outside" = "#666666"),
                          labels = c("GWTS_Inside" = "C.russula\nInside",
                                     "GWTS_Edge" = "C.russula\nEdge",
                                     "Pygmy_Edge" = "S.minutus\nEdge",
                                     "Pygmy_Outside" = "S.minutus\nOutside"),
                          name = element_blank()) +
      theme_bw() +
      theme(axis.text = element_text(size = 14, face= "bold", colour = "black"),
            axis.title=element_text(size=18),
            legend.title= element_text(size= 18),
            legend.text = element_text(size = 14),
            legend.key = element_rect(fill = "white", colour = NA),
            legend.background = element_rect(fill = "white"),
            panel.grid.major = element_line(colour = "grey95"),
            panel.grid.minor = element_line(colour = "grey95")
      ) +
      guides(colour = guide_legend(override.aes = list(size=5)),
             colour = guide_legend("shrew_zone"), fill = FALSE) +
      theme(legend.position = "right")
    
    #p.2.5
    
    ggsave(paste(path4, "/Ireland_samples_NMDS_1_3_", lv[k], "_level.jpeg", sep = ""), 
           plot = p.2.5, device = "jpeg")
    
    write.csv(par.tab, file = paste(path4, "/parameter_table_Jaccard_", lv[k], "_level.csv", sep = ""))
    
  }
}

beep("coin")


#### Homogeneity ####

paths <- c(b.meth1 = "./browett_et_al_2021_assignment_restrictions/filter_method_1/taxonomically_assigned_filtered.RData",
           b.meth2 = "./browett_et_al_2021_assignment_restrictions/filter_method_2/taxonomically_assigned_filtered.RData",
           a.meth1 = "./alberdi_taxa_assignment_restrictions/filter_method_1/taxonomically_assigned_filtered.RData",
           a.meth2 = "./alberdi_taxa_assignment_restrictions/filter_method_2/taxonomically_assigned_filtered.RData")


lv <- c("MOTU", "species", "genus", "family", "order")

## Bray Curtis
for(i in 1:length(paths)){
  
  mt.nam <- names(paths[i])
  load(paths[i])
  
  for(k in 1:length(lv)){
    
    if(lv[k] == "MOTU"){glom.diet <- final.diet}
    else {glom.diet <- tax_glom(final.diet, taxrank = lv[k])}
    
    print(paste("Running ", mt.nam, " at the ", lv[k], " level"))
    # RRA Bray Curtis - MOTU level
    
    dir.create("./Homogeneity_Bray_Curtis")
    dir.create(paste("./Homogeneity_Bray_Curtis/", mt.nam, sep = ""))
    path1 <- paste("./Homogeneity_Bray_Curtis/", mt.nam, "/Homogeneity_glom_", lv[k], "_Bray_Curtis", sep = "")
    dir.create(path1)
    
    mrg.ra = transform_sample_counts(glom.diet, function(x) 100 * x/sum(x))
    options(scipen = 999)
    otu <- abundances(mrg.ra)
    meta <- meta(mrg.ra)
    dist_matrix <- phyloseq::distance(mrg.ra,
                                      method = "bray")
    
    # ---
    ### Groups: Shrew - Country 
    mod <- betadisper(dist_matrix, paste(meta$Shrew,
                                         meta$Country))
    
    dis.stat <- data.frame(group = mod$group, dist = mod$distances)
    
    G.1 <- ggplot(dis.stat, aes(group, dist)) +
      geom_boxplot(fill= c("#440154ff", '#21908dff', "indianred", "#fde725ff"), alpha = 0.35) +
      geom_jitter(aes(colour = group), position = position_jitter(width = 0.1), size=3) +
      stat_summary(fun.y = mean, geom="point", shape=20, size = 10) +
      scale_colour_manual(values= c("#440154ff", '#21908dff', "indianred", "#fde725ff")) + # Manually sets colours according to object "palette"
      scale_x_discrete(labels= c("C. russla\nFrance",
                                 "C. russla\nIreland",
                                 "S. minutus\nFrance",
                                 "S. minutus\nIreland")) + # Manually labels groups using set described in "key" object
      # scale_y_continuous(limits = c(0, 0.55), expand = c(0, 0)) + # Sets Graph Value Limits
      theme_classic() + # Chooses Theme
      theme(axis.text = element_text(size = 10, face= "bold"), # Size of Axis measurements/Values Text
            axis.title=element_text(size=14), # Sets size of Axis Titles. face= "bold" will change the font to bold
            panel.grid.major = element_line(colour = "grey95"), # Sets colour of major grid lines
            panel.grid.minor = element_line(colour = "grey95") # Sets colour of minor grid lines
      ) +
      labs(y= "Level of Dispersion", x = "") + # Labels Axis'
      #facet_wrap(~Shrew, scale = "free_x") +
      theme(legend.position="none")
    
    G.1
    
    ggsave(paste(path1, "/Homogeneity_Shrew_Country_", lv[k], "_level.jpeg", sep = ""), 
           plot = G.1, device = "jpeg")
    
    a1 <- anova(mod)
    a1
    capture.output(a1, file = paste(path1, "_ANOVA_Shrew_Country_Result.txt"))
    tk1 <- TukeyHSD(mod)
    tk1
    capture.output(tk1, file = paste(path1, "_Tukey_Shrew_Country_Result.txt"))
    #plot(TukeyHSD(mod))
    
    set.seed(29)
    pm <- permutest(mod, permutations = 9999)
    pm
    capture.output(pm, file = paste(path1, "_Permutest_Shrew_Country_Result.txt"))
    # Permutest pairwise post-hoc
    ptst <- permutest(mod, pairwise = TRUE, permutations = 9999)
    ptst
    capture.output(ptst, file = paste(path1, "_PostHoc_Pairwise_Permutest_Shrew_Country_Result.txt"))
    # List the significant values (i.e. p < 0.05)
    ptst$pairwise$observed[which(ptst$pairwise$observed < 0.05)]
    lst.mn <- aggregate(mod$distances, list(mod$group), mean)
    lst.mn
    capture.output(lst.mn, file = paste(path1, "_Mean_Values_Shrew_Country_Result.txt"))
    lst.sd <- aggregate(mod$distances, list(mod$group), sd)
    lst.sd
    capture.output(lst.sd, file = paste(path1, "_Standard_Deviation_Values_Shrew_Country_Result.txt"))
    
    
    
    # ---
    ### Groups: Shrew - Country - Season
    mod <- betadisper(dist_matrix, paste(meta$Shrew,
                                         meta$Country,
                                         meta$Season))
    
    dis.stat <- data.frame(group = mod$group, dist = mod$distances)
    
    G.2 <- ggplot(dis.stat, aes(group, dist)) +
      geom_boxplot(aes(fill= group), alpha = 0.35) +
      geom_jitter(aes(colour = group), position = position_jitter(width = 0.1), size=3) +
      stat_summary(fun.y = mean, geom="point", shape=20, size = 10) +
      theme_classic() + # Chooses Theme
      theme(axis.text = element_text(size = 10, face= "bold"), # Size of Axis measurements/Values Text
            axis.title=element_text(size=14), # Sets size of Axis Titles. face= "bold" will change the font to bold
            panel.grid.major = element_line(colour = "grey95"), # Sets colour of major grid lines
            panel.grid.minor = element_line(colour = "grey95")) + # Sets colour of minor grid lines
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(y= "Level of Dispersion", x = "") + # Labels Axis'
      theme(legend.position="none")
    
    G.2
    
    ggsave(paste(path1, "/Homogeneity_Shrew_Country_Season_", lv[k], "_level.jpeg", sep = ""), 
           plot = G.2, device = "jpeg")
    
    a1 <- anova(mod)
    a1
    capture.output(a1, file = paste(path1, "_ANOVA_Shrew_Country_Season_Result.txt"))
    tk1 <- TukeyHSD(mod)
    tk1
    capture.output(tk1, file = paste(path1, "_Tukey_Shrew_Country_Season_Result.txt"))
    #plot(TukeyHSD(mod))
    
    set.seed(29)
    pm <- permutest(mod, permutations = 9999)
    pm
    capture.output(pm, file = paste(path1, "_Permutest_Shrew_Country_Season_Result.txt"))
    # Permutest pairwise post-hoc
    ptst <- permutest(mod, pairwise = TRUE, permutations = 9999)
    ptst
    capture.output(ptst, file = paste(path1, "_PostHoc_Pairwise_Permutest_Shrew_Country_Season_Result.txt"))
    # List the significant values (i.e. p < 0.05)
    ptst$pairwise$observed[which(ptst$pairwise$observed < 0.05)]
    lst.mn <- aggregate(mod$distances, list(mod$group), mean)
    lst.mn
    capture.output(lst.mn, file = paste(path1, "_Mean_Values_Shrew_Country_Season_Result.txt"))
    lst.sd <- aggregate(mod$distances, list(mod$group), sd)
    lst.sd
    capture.output(lst.sd, file = paste(path1, "_Standard_Deviation_Values_Shrew_Country_Season_Result.txt"))
    
    
    
    # ---
    ### Groups: Shrew - Country - Season - Zone
    
    mod <- betadisper(dist_matrix, paste(meta$Shrew,
                                         meta$Country,
                                         meta$Season,
                                         meta$Zone))
    
    dis.stat <- data.frame(group = mod$group, dist = mod$distances)
    
    G.3 <- ggplot(dis.stat, aes(group, dist)) +
      geom_boxplot(aes(fill= group), alpha = 0.35) +
      geom_jitter(aes(colour = group), position = position_jitter(width = 0.1), size=3) +
      stat_summary(fun.y = mean, geom="point", shape=20, size = 10) +
      theme_classic() + # Chooses Theme
      theme(axis.text = element_text(size = 10, face= "bold"), # Size of Axis measurements/Values Text
            axis.title=element_text(size=14), # Sets size of Axis Titles. face= "bold" will change the font to bold
            panel.grid.major = element_line(colour = "grey95"), # Sets colour of major grid lines
            panel.grid.minor = element_line(colour = "grey95")) + # Sets colour of minor grid lines
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(y= "Level of Dispersion", x = "") + # Labels Axis'
      theme(legend.position="none")
    
    G.3
    
    ggsave(paste(path1, "/Homogeneity_Shrew_Country_Season_Zone_", lv[k], "_level.jpeg", sep = ""), 
           plot = G.3, device = "jpeg")
    
    a1 <- anova(mod)
    a1
    capture.output(a1, file = paste(path1, "_ANOVA_Shrew_Country_Season_Zone_Result.txt"))
    tk1 <- TukeyHSD(mod)
    tk1
    capture.output(tk1, file = paste(path1, "_Tukey_Shrew_Country_Season_Zone_Result.txt"))
    #plot(TukeyHSD(mod))
    
    set.seed(29)
    pm <- permutest(mod, permutations = 9999)
    pm
    capture.output(pm, file = paste(path1, "_Permutest_Shrew_Country_Season_Zone_Result.txt"))
    # Permutest pairwise post-hoc
    ptst <- permutest(mod, pairwise = TRUE, permutations = 9999)
    ptst
    capture.output(ptst, file = paste(path1, "_PostHoc_Pairwise_Permutest_Shrew_Country_Season_Zone_Result.txt"))
    # List the significant values (i.e. p < 0.05)
    ptst$pairwise$observed[which(ptst$pairwise$observed < 0.05)]
    lst.mn <- aggregate(mod$distances, list(mod$group), mean)
    lst.mn
    capture.output(lst.mn, file = paste(path1, "_Mean_Values_Shrew_Country_Season_Zone_Result.txt"))
    lst.sd <- aggregate(mod$distances, list(mod$group), sd)
    lst.sd
    capture.output(lst.sd, file = paste(path1, "_Standard_Deviation_Values_Shrew_Country_Season_Zone_Result.txt"))
    
    # ---
    ### Groups: Shrew - Country - Season - Transect
    
    mod <- betadisper(dist_matrix, paste(meta$Shrew,
                                         meta$Country,
                                         meta$Season,
                                         meta$Transect))
    
    dis.stat <- data.frame(group = mod$group, dist = mod$distances)
    
    G.4 <- ggplot(dis.stat, aes(group, dist)) +
      geom_boxplot(aes(fill= group), alpha = 0.35) +
      geom_jitter(aes(colour = group), position = position_jitter(width = 0.1), size=3) +
      stat_summary(fun.y = mean, geom="point", shape=20, size = 10) +
      theme_classic() + # Chooses Theme
      theme(axis.text = element_text(size = 10, face= "bold"), # Size of Axis measurements/Values Text
            axis.title=element_text(size=14), # Sets size of Axis Titles. face= "bold" will change the font to bold
            panel.grid.major = element_line(colour = "grey95"), # Sets colour of major grid lines
            panel.grid.minor = element_line(colour = "grey95")) + # Sets colour of minor grid lines
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(y= "Level of Dispersion", x = "") + # Labels Axis'
      theme(legend.position="none")
    
    G.4
    
    ggsave(paste(path1, "/Homogeneity_Shrew_Country_Season_Zone_", lv[k], "_level.jpeg", sep = ""), 
           plot = G.4, device = "jpeg")
    
    a1 <- anova(mod)
    a1
    capture.output(a1, file = paste(path1, "_ANOVA_Shrew_Country_Season_Transect_Result.txt"))
    tk1 <- TukeyHSD(mod)
    tk1
    capture.output(tk1, file = paste(path1, "_Tukey_Shrew_Country_Season_Transect_Result.txt"))
    #plot(TukeyHSD(mod))
    
    set.seed(29)
    pm <- permutest(mod, permutations = 9999)
    pm
    capture.output(pm, file = paste(path1, "_Permutest_Shrew_Country_Season_Transect_Result.txt"))
    # Permutest pairwise post-hoc
    ptst <- permutest(mod, pairwise = TRUE, permutations = 9999)
    ptst
    capture.output(ptst, file = paste(path1, "_PostHoc_Pairwise_Permutest_Shrew_Country_Season_Transect_Result.txt"))
    # List the significant values (i.e. p < 0.05)
    ptst$pairwise$observed[which(ptst$pairwise$observed < 0.05)]
    lst.mn <- aggregate(mod$distances, list(mod$group), mean)
    lst.mn
    capture.output(lst.mn, file = paste(path1, "_Mean_Values_Shrew_Country_Season_Transect_Result.txt"))
    lst.sd <- aggregate(mod$distances, list(mod$group), sd)
    lst.sd
    capture.output(lst.sd, file = paste(path1, "_Standard_Deviation_Values_Shrew_Country_Season_Transect_Result.txt"))
    
  }
  
  
}


## Jaccard
for(i in 1:length(paths)){
  
  mt.nam <- names(paths[i])
  load(paths[i])
  
  for(k in 1:length(lv)){
    
    if(lv[k] == "MOTU"){glom.diet <- final.diet}
    else {glom.diet <- tax_glom(final.diet, taxrank = lv[k])}
    
    print(paste("Running ", mt.nam, " at the ", lv[k], " level"))
    # RRA Bray Curtis - MOTU level
    
    dir.create("./Homogeneity_Jaccard")
    dir.create(paste("./Homogeneity_Jaccard/", mt.nam, sep = ""))
    path1 <- paste("./Homogeneity_Jaccard/", mt.nam, "/Homogeneity_glom_", lv[k], "_Jaccard", sep = "")
    dir.create(path1)
    
    mrg.ra = transform_sample_counts(glom.diet, function(x) 100 * x/sum(x))
    options(scipen = 999)
    otu <- abundances(mrg.ra)
    otu[otu>.0099999999]=1
    otu[otu<.01]=0
    meta <- meta(mrg.ra)
    dist_matrix <- phyloseq::distance(mrg.ra,
                                      method = "jaccard")
    
    
    # ---
    ### Groups: Shrew - Country 
    mod <- betadisper(dist_matrix, paste(meta$Shrew,
                                         meta$Country))
    
    dis.stat <- data.frame(group = mod$group, dist = mod$distances)
    
    G.1 <- ggplot(dis.stat, aes(group, dist)) +
      geom_boxplot(fill= c("#440154ff", '#21908dff', "indianred", "#fde725ff"), alpha = 0.35) +
      geom_jitter(aes(colour = group), position = position_jitter(width = 0.1), size=3) +
      stat_summary(fun.y = mean, geom="point", shape=20, size = 10) +
      scale_colour_manual(values= c("#440154ff", '#21908dff', "indianred", "#fde725ff")) + # Manually sets colours according to object "palette"
      scale_x_discrete(labels= c("C. russla\nFrance",
                                 "C. russla\nIreland",
                                 "S. minutus\nFrance",
                                 "S. minutus\nIreland")) + # Manually labels groups using set described in "key" object
      # scale_y_continuous(limits = c(0, 0.55), expand = c(0, 0)) + # Sets Graph Value Limits
      theme_classic() + # Chooses Theme
      theme(axis.text = element_text(size = 10, face= "bold"), # Size of Axis measurements/Values Text
            axis.title=element_text(size=14), # Sets size of Axis Titles. face= "bold" will change the font to bold
            panel.grid.major = element_line(colour = "grey95"), # Sets colour of major grid lines
            panel.grid.minor = element_line(colour = "grey95") # Sets colour of minor grid lines
      ) +
      labs(y= "Level of Dispersion", x = "") + # Labels Axis'
      #facet_wrap(~Shrew, scale = "free_x") +
      theme(legend.position="none")
    
    G.1
    
    ggsave(paste(path1, "/Homogeneity_Shrew_Country_", lv[k], "_level.jpeg", sep = ""), 
           plot = G.1, device = "jpeg")
    
    a1 <- anova(mod)
    a1
    capture.output(a1, file = paste(path1, "_ANOVA_Shrew_Country_Result.txt"))
    tk1 <- TukeyHSD(mod)
    tk1
    capture.output(tk1, file = paste(path1, "_Tukey_Shrew_Country_Result.txt"))
    #plot(TukeyHSD(mod))
    
    set.seed(29)
    pm <- permutest(mod, permutations = 9999)
    pm
    capture.output(pm, file = paste(path1, "_Permutest_Shrew_Country_Result.txt"))
    # Permutest pairwise post-hoc
    ptst <- permutest(mod, pairwise = TRUE, permutations = 9999)
    ptst
    capture.output(ptst, file = paste(path1, "_PostHoc_Pairwise_Permutest_Shrew_Country_Result.txt"))
    # List the significant values (i.e. p < 0.05)
    ptst$pairwise$observed[which(ptst$pairwise$observed < 0.05)]
    lst.mn <- aggregate(mod$distances, list(mod$group), mean)
    lst.mn
    capture.output(lst.mn, file = paste(path1, "_Mean_Values_Shrew_Country_Result.txt"))
    lst.sd <- aggregate(mod$distances, list(mod$group), sd)
    lst.sd
    capture.output(lst.sd, file = paste(path1, "_Standard_Deviation_Values_Shrew_Country_Result.txt"))
    
    
    
    # ---
    ### Groups: Shrew - Country - Season
    mod <- betadisper(dist_matrix, paste(meta$Shrew,
                                         meta$Country,
                                         meta$Season))
    
    dis.stat <- data.frame(group = mod$group, dist = mod$distances)
    
    G.2 <- ggplot(dis.stat, aes(group, dist)) +
      geom_boxplot(aes(fill= group), alpha = 0.35) +
      geom_jitter(aes(colour = group), position = position_jitter(width = 0.1), size=3) +
      stat_summary(fun.y = mean, geom="point", shape=20, size = 10) +
      theme_classic() + # Chooses Theme
      theme(axis.text = element_text(size = 10, face= "bold"), # Size of Axis measurements/Values Text
            axis.title=element_text(size=14), # Sets size of Axis Titles. face= "bold" will change the font to bold
            panel.grid.major = element_line(colour = "grey95"), # Sets colour of major grid lines
            panel.grid.minor = element_line(colour = "grey95")) + # Sets colour of minor grid lines
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(y= "Level of Dispersion", x = "") + # Labels Axis'
      theme(legend.position="none")
    
    G.2
    
    ggsave(paste(path1, "/Homogeneity_Shrew_Country_Season_", lv[k], "_level.jpeg", sep = ""), 
           plot = G.2, device = "jpeg")
    
    a1 <- anova(mod)
    a1
    capture.output(a1, file = paste(path1, "_ANOVA_Shrew_Country_Season_Result.txt"))
    tk1 <- TukeyHSD(mod)
    tk1
    capture.output(tk1, file = paste(path1, "_Tukey_Shrew_Country_Season_Result.txt"))
    #plot(TukeyHSD(mod))
    
    set.seed(29)
    pm <- permutest(mod, permutations = 9999)
    pm
    capture.output(pm, file = paste(path1, "_Permutest_Shrew_Country_Season_Result.txt"))
    # Permutest pairwise post-hoc
    ptst <- permutest(mod, pairwise = TRUE, permutations = 9999)
    ptst
    capture.output(ptst, file = paste(path1, "_PostHoc_Pairwise_Permutest_Shrew_Country_Season_Result.txt"))
    # List the significant values (i.e. p < 0.05)
    ptst$pairwise$observed[which(ptst$pairwise$observed < 0.05)]
    lst.mn <- aggregate(mod$distances, list(mod$group), mean)
    lst.mn
    capture.output(lst.mn, file = paste(path1, "_Mean_Values_Shrew_Country_Season_Result.txt"))
    lst.sd <- aggregate(mod$distances, list(mod$group), sd)
    lst.sd
    capture.output(lst.sd, file = paste(path1, "_Standard_Deviation_Values_Shrew_Country_Season_Result.txt"))
    
    
    
    # ---
    ### Groups: Shrew - Country - Season - Zone
    
    mod <- betadisper(dist_matrix, paste(meta$Shrew,
                                         meta$Country,
                                         meta$Season,
                                         meta$Zone))
    
    dis.stat <- data.frame(group = mod$group, dist = mod$distances)
    
    G.3 <- ggplot(dis.stat, aes(group, dist)) +
      geom_boxplot(aes(fill= group), alpha = 0.35) +
      geom_jitter(aes(colour = group), position = position_jitter(width = 0.1), size=3) +
      stat_summary(fun.y = mean, geom="point", shape=20, size = 10) +
      theme_classic() + # Chooses Theme
      theme(axis.text = element_text(size = 10, face= "bold"), # Size of Axis measurements/Values Text
            axis.title=element_text(size=14), # Sets size of Axis Titles. face= "bold" will change the font to bold
            panel.grid.major = element_line(colour = "grey95"), # Sets colour of major grid lines
            panel.grid.minor = element_line(colour = "grey95")) + # Sets colour of minor grid lines
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(y= "Level of Dispersion", x = "") + # Labels Axis'
      theme(legend.position="none")
    
    G.3
    
    ggsave(paste(path1, "/Homogeneity_Shrew_Country_Season_Zone_", lv[k], "_level.jpeg", sep = ""), 
           plot = G.3, device = "jpeg")
    
    a1 <- anova(mod)
    a1
    capture.output(a1, file = paste(path1, "_ANOVA_Shrew_Country_Season_Zone_Result.txt"))
    tk1 <- TukeyHSD(mod)
    tk1
    capture.output(tk1, file = paste(path1, "_Tukey_Shrew_Country_Season_Zone_Result.txt"))
    #plot(TukeyHSD(mod))
    
    set.seed(29)
    pm <- permutest(mod, permutations = 9999)
    pm
    capture.output(pm, file = paste(path1, "_Permutest_Shrew_Country_Season_Zone_Result.txt"))
    # Permutest pairwise post-hoc
    ptst <- permutest(mod, pairwise = TRUE, permutations = 9999)
    ptst
    capture.output(ptst, file = paste(path1, "_PostHoc_Pairwise_Permutest_Shrew_Country_Season_Zone_Result.txt"))
    # List the significant values (i.e. p < 0.05)
    ptst$pairwise$observed[which(ptst$pairwise$observed < 0.05)]
    lst.mn <- aggregate(mod$distances, list(mod$group), mean)
    lst.mn
    capture.output(lst.mn, file = paste(path1, "_Mean_Values_Shrew_Country_Season_Zone_Result.txt"))
    lst.sd <- aggregate(mod$distances, list(mod$group), sd)
    lst.sd
    capture.output(lst.sd, file = paste(path1, "_Standard_Deviation_Values_Shrew_Country_Season_Zone_Result.txt"))
    
    # ---
    ### Groups: Shrew - Country - Season - Transect
    
    mod <- betadisper(dist_matrix, paste(meta$Shrew,
                                         meta$Country,
                                         meta$Season,
                                         meta$Transect))
    
    dis.stat <- data.frame(group = mod$group, dist = mod$distances)
    
    G.4 <- ggplot(dis.stat, aes(group, dist)) +
      geom_boxplot(aes(fill= group), alpha = 0.35) +
      geom_jitter(aes(colour = group), position = position_jitter(width = 0.1), size=3) +
      stat_summary(fun.y = mean, geom="point", shape=20, size = 10) +
      theme_classic() + # Chooses Theme
      theme(axis.text = element_text(size = 10, face= "bold"), # Size of Axis measurements/Values Text
            axis.title=element_text(size=14), # Sets size of Axis Titles. face= "bold" will change the font to bold
            panel.grid.major = element_line(colour = "grey95"), # Sets colour of major grid lines
            panel.grid.minor = element_line(colour = "grey95")) + # Sets colour of minor grid lines
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(y= "Level of Dispersion", x = "") + # Labels Axis'
      theme(legend.position="none")
    
    G.4
    
    ggsave(paste(path1, "/Homogeneity_Shrew_Country_Season_Zone_", lv[k], "_level.jpeg", sep = ""), 
           plot = G.4, device = "jpeg")
    
    a1 <- anova(mod)
    a1
    capture.output(a1, file = paste(path1, "_ANOVA_Shrew_Country_Season_Transect_Result.txt"))
    tk1 <- TukeyHSD(mod)
    tk1
    capture.output(tk1, file = paste(path1, "_Tukey_Shrew_Country_Season_Transect_Result.txt"))
    #plot(TukeyHSD(mod))
    
    set.seed(29)
    pm <- permutest(mod, permutations = 9999)
    pm
    capture.output(pm, file = paste(path1, "_Permutest_Shrew_Country_Season_Transect_Result.txt"))
    # Permutest pairwise post-hoc
    ptst <- permutest(mod, pairwise = TRUE, permutations = 9999)
    ptst
    capture.output(ptst, file = paste(path1, "_PostHoc_Pairwise_Permutest_Shrew_Country_Season_Transect_Result.txt"))
    # List the significant values (i.e. p < 0.05)
    ptst$pairwise$observed[which(ptst$pairwise$observed < 0.05)]
    lst.mn <- aggregate(mod$distances, list(mod$group), mean)
    lst.mn
    capture.output(lst.mn, file = paste(path1, "_Mean_Values_Shrew_Country_Season_Transect_Result.txt"))
    lst.sd <- aggregate(mod$distances, list(mod$group), sd)
    lst.sd
    capture.output(lst.sd, file = paste(path1, "_Standard_Deviation_Values_Shrew_Country_Season_Transect_Result.txt"))
    
  }
  
  
}


#### PERMANOVA ####

final.diet.g <- tax_glom(final.diet, taxrank = "order")

mrg.ra = transform_sample_counts(final.diet.g, function(x) 100 * x/sum(x))
options(scipen = 999)
otu <- abundances(mrg.ra)
meta <- meta(mrg.ra)
dist_matrix <- phyloseq::distance(mrg.ra,
                                  method = "jaccard")

set.seed(59009)
permanova <- adonis(t(otu) ~ Country * Shrew,
                    data = meta,
                    permutations=999,
                    method = "jaccard")
#method = "jaccard")

permanova


## PERMANOVAS smaller groups

mrg.ra = transform_sample_counts(final.diet.g, function(x) 100 * x/sum(x))
mrg.ra <- subset_samples(mrg.ra, shrew_country == "GWTS_Ireland")
otu <- abundances(mrg.ra)
meta <- meta(mrg.ra)

set.seed(1234)
permanova <- adonis(t(otu) ~ Season + Transect + Zone + Trap_Site,# strata = meta$Trap_Site,
                    data = meta,
                    permutations=9999,
                    method = "jaccard")
#method = "jaccard")
permanova

#pairwise.adonis(t(otu), meta$Transect)
#pairwise.adonis(t(otu), meta$Trap_Site)

permanova <- adonis(t(otu) ~ Season + Zone + Transect + Trap_Site,# strata = meta$Trap_Site,
                    data = meta,
                    permutations=9999,
                    method = "jaccard")
#method = "jaccard")
permanova

mrg.ra = transform_sample_counts(final.diet.g, function(x) 100 * x/sum(x))
mrg.ra <- subset_samples(mrg.ra, shrew_country == "GWTS_France")
otu <- abundances(mrg.ra)
meta <- meta(mrg.ra)

set.seed(4526)
permanova <- adonis(t(otu) ~ Season + Trap_Site,# strata = meta$Trap_Site,
                    data = meta,
                    permutations=9999,
                    method = "jaccard")
#method = "jaccard")
permanova


mrg.ra = transform_sample_counts(final.diet.g, function(x) 100 * x/sum(x))
mrg.ra <- subset_samples(mrg.ra, shrew_country == "Pygmy_Ireland")
otu <- abundances(mrg.ra)
meta <- meta(mrg.ra)

set.seed(4526)
permanova <- adonis(t(otu) ~ Season + Transect + Zone + Trap_Site,# strata = meta$Trap_Site,
                    data = meta,
                    permutations=9999,
                    method = "jaccard")
#method = "jaccard")
permanova


permanova <- adonis(t(otu) ~ Season + Zone + Transect + Trap_Site, #strata = meta$Transect,
                    data = meta,
                    permutations=9999,
                    method = "jaccard")
#method = "jaccard")
permanova

