
#### Load Libraries ####

library(phyloseq)
library(ggplot2)
library(dplyr)
library(vegan)
library(decontam)
library(stringr)
library(knitr)
library(hilldiv)
library(kableExtra)
library(reshape2)
library(gridExtra)

#### Custom Pallette for future plots ####

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

# to stop numbersbeing shown as a power to etc
options(scipen = 999)

#### Load Dataset ####

path <- "browett_et_al_2021_assignment_restrictions/taxonomically_assigned_not_filter.RData" # our filter set
path <- "alberdi_taxa_assignment_restrictions/taxonomically_assigned_not_filter.RData" # alberdi inspired

load(path)

gillet.phylo

#### Check Sample Sizes and Read Depth ####

df <- as.data.frame(sample_data(gillet.phylo)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(gillet.phylo)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Shrew)) + 
  geom_point() 

samples.phylo <- subset_samples(gillet.phylo, Control != "Yes")
median(sample_sums(samples.phylo)) # Median number of reads per individual
range(sample_sums(samples.phylo)) # Range of read count per individual
mean(sample_sums(samples.phylo)) # Average read count (in total) for each sample


#### Use Decontam to identify potential contaminants ####

sample_data(gillet.phylo)$is.neg <- sample_data(gillet.phylo)$Control == "Yes"
contamdf.prev <- isContaminant(gillet.phylo, 
                               method="prevalence", 
                               neg="is.neg",
                               threshold = 0.5)
#table(contamdf.prev$contaminant)

tax_table(gillet.phylo) <- cbind(tax_table(gillet.phylo), contaminant = contamdf.prev$contaminant) # adds contaminant TRUE or FALSE column to taxa_table

#df <- as.data.frame(tax_table(phylo))

cont <- which(contamdf.prev$contaminant)
con_table <- as.data.frame(tax_table(gillet.phylo)[cont,2:7])
for (i in 1:length(cont)){
  con_table$prev[i] <- sum(otu_table(gillet.phylo)[cont[i],] > 0)
  con_table$count[i] <- taxa_sums(gillet.phylo)[cont[i]]
}


con_table
#kable(con_table, caption = "List of MOTUs identified by 'decontam' as contaminants")


#### Exploring the blanks ####

cntrl <- subset_samples(gillet.phylo, Control == "Yes")
cntrl <- prune_taxa(taxa_sums(cntrl) > 0, cntrl)
ntaxa(cntrl)
median(sample_sums(cntrl)) # Median number of reads per individual
range(sample_sums(cntrl)) # Range of read count per individual
mean(sample_sums(cntrl)) # Average read count (in total) for each sample
sum(sample_sums(cntrl)) # Total number of reads in blanks

cntrl.glom <- tax_glom(cntrl, taxrank = "class")
cntrl.bar <- phyloseq::plot_bar(cntrl.glom) # extracts information needed for barplots
data <- cntrl.bar$data
#data <- psmelt(cntrl.bar)
p1 <- ggplot(data, aes(x= Sample, y=Abundance, fill=class))

p1 <- p1 + geom_bar(stat="identity", color="black") +
  scale_fill_manual(values = pal.o) +
  facet_wrap(~Shrew, scale = "free_x") +
  ggtitle("Read abundance in negative controls") +
  theme_gray() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ra.cntrl = transform_sample_counts(cntrl.glom, function(x) 100 * x/sum(x))
ra.cntrl.bar <- phyloseq::plot_bar(ra.cntrl) # extracts information needed for barplots
data <- ra.cntrl.bar$data
p2 <- ggplot(data, aes(x= Sample, y=Abundance, fill=class))

p2 <- p2 + geom_bar(stat="identity", color="black") +
  scale_fill_manual(values = pal.o) +
  facet_wrap(~Shrew, scale = "free_x") +
  ggtitle("Relative read abundance in negative controls") +
  theme_gray() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

grid.arrange(p1, p2)


blnk.tx <- as.data.frame(as.matrix(tax_table(cntrl)[,4:7]))
blnk.ot <- as.data.frame(as.matrix(otu_table(cntrl)[,]))
blnk.all <- cbind(blnk.tx, blnk.ot)
blnk.all$all.blnk <- taxa_sums(cntrl)
blnk.all$total <- taxa_sums(gillet.phylo)[taxa_names(cntrl)[]]
for (i in 1:nrow(blnk.all)) {blnk.all$perc[i] <- (blnk.all$all.blnk[i] / blnk.all$total[i]) * 100 }

blnk.all
kable(blnk.all, caption = "List of MOTUs identified in blanks")

#### Prominence of GWTS reads in Pygmy samples and vice versa ####

host.phylo <- subset_taxa(gillet.phylo, (order %in% c("Eulipotyphla")))
host.phylo <- tax_glom(host.phylo, taxrank = "genus")
tax_table(host.phylo)
taxa_sums(host.phylo)

gwts.phylo <- subset_samples(host.phylo, Shrew == "GWTS")
ps.phylo <- subset_samples(host.phylo, Shrew == "Pygmy")

otu_table(gwts.phylo)
otu_table(ps.phylo)

(mean(otu_table(gwts.phylo)[2,]) / as.numeric(taxa_sums(host.phylo)[2])) * 100
(mean(otu_table(gwts.phylo)[4,]) / as.numeric(taxa_sums(host.phylo)[4])) * 100

(mean(otu_table(ps.phylo)[1,]) / as.numeric(taxa_sums(host.phylo)[1])) * 100
(mean(otu_table(ps.phylo)[3,]) / as.numeric(taxa_sums(host.phylo)[3])) * 100


#### Method 1: Removing MOTUs that are prominant in Blanks ####

### This will make a table with all the MOTUs found in the negative controls, 
## of which over 2% of the total reads were found in the negative controls

table(blnk.all[,"perc"] > 2) # 2 is for 2%
table(blnk.all[,"perc"] > 1)
table(blnk.all[,"perc"] > 0.1)
table(blnk.all[,"perc"] > 0.01)

tab.nam <- blnk.all[,"perc"] > 0.1 
tab.df <- blnk.all[tab.nam,]

# don't include MOTU_70 and MOTU_81, which correspond to rows 2 and 3
tab.df <- tab.df[-c(2,3),]
removeTaxa <- rownames(tab.df) # Lists the MOTUs to remove 
phy.obj <- subset_taxa(gillet.phylo, !(taxa_names(gillet.phylo) %in% removeTaxa))
phy.obj

## Filtering Low Read MOTUs ##

samples.phylo <- subset_samples(phy.obj, Control != "Yes")
diet.prey <- subset_taxa(samples.phylo, !(class %in% c("Mammalia", 
                                                       "none",
                                                       "Actinopteri",
                                                       "Bdelloidea",
                                                       "Udeonychophora", # velvet worms
                                                       "Merostomata", # horse shoe crabs
                                                       "Gammaproteobacteria", # bacteria
                                                       "Magnoliopsida", # plants
                                                       "Monogononta", # rotifers
                                                       "Dothideomycetes", # fungi
                                                       "Trebouxiophyceae", # green algae
                                                       "Chondrichthyes", # Cartilaginous fish
                                                       "Mucoromycetes", # fungi
                                                       "Phylum_Endomyxa", # micro things
                                                       "Eutardigrada", # tartigrades!!
                                                       "Elardia", # Amoebas
                                                       "Cephalopoda", # Cephalopods
                                                       "Amphibia", # Amphibians
                                                       "Aves", # Birds
                                                       "Chromadorea", # roundworms
                                                       "Hexanauplia",  # parasitic crustaceans
                                                       "Kingdom_Metazoa",
                                                       "Kingdom_",
                                                       "Phylum_Discosea", # amoebas
                                                       "Branchiopoda", # marine crustaceans
                                                       "Phylum_Nematoda")))

#df <- data.frame(sort(sample_sums(diet.prey), decreasing = TRUE))
#df <- data.frame(sample_sums(diet.prey))
#write.csv(df, file = "./check_sample_sizes.csv")
sampl.filt <- prune_samples(sample_sums(diet.prey) > 1000, diet.prey)
otu.tab <- as.data.frame(otu_table(sampl.filt))
new.otu.tab <- copy_filt(otu.tab, 0.001)
new.otu.tab <- as.matrix(new.otu.tab)
otu_table(sampl.filt) <- otu_table(new.otu.tab, taxa_are_rows = TRUE)
sampl.filt
final.diet <- prune_taxa(taxa_sums(sampl.filt) > 0, sampl.filt)
final.diet
sort(sample_sums(final.diet))

#dir.create("./browett_et_al_2021_assignment_restrictions/filter_method_1")
#dir.create("./alberdi_taxa_assignment_restrictions/filter_method_1")

#save(final.diet, file = "./browett_et_al_2021_assignment_restrictions/filter_method_1/taxonomically_assigned_filtered.RData")
save(final.diet, file = "./alberdi_taxa_assignment_restrictions/filter_method_1/taxonomically_assigned_filtered.RData")


#### Method 2: Removing the max number of reads in each contaminant from dataset ####

phy.obj <- gillet.phylo

samples.phylo <- subset_samples(phy.obj, Control != "Yes")
diet.prey <- subset_taxa(samples.phylo, !(class %in% c("Mammalia", 
                                                       "none",
                                                       "Actinopteri",
                                                       "Bdelloidea",
                                                       "Udeonychophora", # velvet worms
                                                       "Merostomata", # horse shoe crabs
                                                       "Gammaproteobacteria", # bacteria
                                                       "Magnoliopsida", # plants
                                                       "Monogononta", # rotifers
                                                       "Dothideomycetes", # fungi
                                                       "Trebouxiophyceae", # green algae
                                                       "Chondrichthyes", # Cartilaginous fish
                                                       "Mucoromycetes", # fungi
                                                       "Phylum_Endomyxa", # micro things
                                                       "Eutardigrada", # tartigrades!!
                                                       "Elardia", # Amoebas
                                                       "Cephalopoda", # Cephalopods
                                                       "Amphibia", # Amphibians
                                                       "Aves", # Birds
                                                       "Chromadorea", # roundworms
                                                       "Hexanauplia",  # parasitic crustaceans
                                                       "Kingdom_Metazoa",
                                                       "Kingdom_",
                                                       "Phylum_Discosea", # amoebas
                                                       "Branchiopoda", # marine crustaceans
                                                       "Phylum_Nematoda")))

#df <- data.frame(sort(sample_sums(diet.prey), decreasing = TRUE))
#df <- data.frame(sample_sums(diet.prey))
#write.csv(df, file = "./check_sample_sizes.csv")
sampl.filt <- prune_samples(sample_sums(diet.prey) > 1000, diet.prey)
otu.tab <- as.data.frame(otu_table(sampl.filt))

for(i in 1:nrow(blnk.all)){
  motu.nm <- rownames(blnk.all[i,])
  if (is.na(otu.tab[motu.nm,1]) == FALSE){
    otu.tab[motu.nm,] <- otu.tab[motu.nm,] - max(blnk.all[i,5:29])
    otu.tab[motu.nm,][otu.tab[motu.nm,] < 0] <- 0 
  }
}

new.otu.tab <- copy_filt(otu.tab, 0.001)
new.otu.tab <- as.matrix(new.otu.tab)
otu_table(sampl.filt) <- otu_table(new.otu.tab, taxa_are_rows = TRUE)
sampl.filt
final.diet <- prune_taxa(taxa_sums(sampl.filt) > 0, sampl.filt)
final.diet
sort(sample_sums(final.diet))
# Sample FR62 seemed problematic and decreased to a read depth of 70 when MOTUs were filtered
# The next two lines are to remove FR62 and it's MOTUs
final.diet <- subset_samples(final.diet, sample_names(final.diet) != "FR62")
#final.diet <- prune_samples(sample_sums(final.diet) > 950, final.diet)
final.diet <- prune_taxa(taxa_sums(final.diet) > 0, final.diet)
final.diet

## Merge MOTUs belonging to the same species (result of over clustering)
n <- which(as.numeric(tax_table(final.diet)[,"pident"]) > 98)
m <- taxa_names(final.diet)[n]
p <- taxa_names(final.diet)[-n]

chck <- prune_taxa(m, final.diet)
chck2 <- tax_glom(chck, taxrank="species")
chck2
chck3 <- prune_taxa(p, final.diet)
chck3
chck4 <- merge_phyloseq(chck2, chck3)
chck4

final.diet <- chck4


#dir.create("./browett_et_al_2021_assignment_restrictions/filter_method_2")
#dir.create("./alberdi_taxa_assignment_restrictions/filter_method_2")

#save(final.diet, file = "./browett_et_al_2021_assignment_restrictions/filter_method_2/taxonomically_assigned_filtered.RData")
save(final.diet, file = "./alberdi_taxa_assignment_restrictions/filter_method_2/taxonomically_assigned_filtered.RData")






