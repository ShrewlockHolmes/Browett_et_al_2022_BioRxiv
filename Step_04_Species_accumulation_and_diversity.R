#### libraries ####

library(vegan)
library(phyloseq)
library(spaa)
library(ggplot2)
library(gridExtra)
library(hilldiv)

#### Load Dataset ####


# Browett Method 1
load("./browett_et_al_2021_assignment_restrictions/filter_method_1/taxonomically_assigned_filtered.RData")

# Browett Method 2
load("./browett_et_al_2021_assignment_restrictions/filter_method_2/taxonomically_assigned_filtered.RData")

# Alberdi Method 1
load("./alberdi_taxa_assignment_restrictions/filter_method_1/taxonomically_assigned_filtered.RData")

# Alberdi Method 2
load("./alberdi_taxa_assignment_restrictions/filter_method_2/taxonomically_assigned_filtered.RData")

final.diet

#### Seq Depth ####

new.otu.tab <- otu_table(final.diet)

# hill_div packages assessment of read depth per sample, 
# according to a species richness equivilent
depth_cov(new.otu.tab,
          qvalue = 0)

# hill_div packages assessment of read depth per sample, 
# according to a shannon diversity equivilent
depth_cov(new.otu.tab,
          qvalue = 1)

#### Accumulation Curves ####

dir.create("./accumulation_curves")

#### Pygmy Ireland ####
par(mar = c(4.1, 4.1, 2.1, 2.1))
sbset <- subset_samples(final.diet, shrew_country == "Pygmy_Ireland")
sbset <- prune_taxa(taxa_sums(sbset) > 0, sbset)

new <- t(as.data.frame(as.matrix(otu_table(sbset))))
specpool(new)

#build the species accumulation curve & rarefaction curve (expected)
aurora.specaccum <- specaccum(new[,],method = "rarefaction")
#build a expected curve (randomization for boxplot comparison)
aurora.specaccum.rand <- specaccum(new[,], "random")
#plot both curves ("observed" vs "randomized")
jpeg("./accumulation_curves/Pygmy_Ireland_accumulation_MOTUs.jpeg", width = 600, height = 450)
plot(aurora.specaccum.rand,ci.type="poly", 
     col="blue", lwd=2, ci.lty=0, ci.col="lightblue",
     xlab = "Number of Shrews", ylab = "Number of Taxa")
boxplot(aurora.specaccum.rand, col="yellow", add=TRUE, pch="+")

text(20, 590, "Chao = 1670.299 +/- 169.39", cex = 1)
legend('bottomright', "S. minutus - Ireland\nMOTUs\n")
dev.off()


# Species
sbset.gn <- tax_glom(sbset, taxrank = "species")
new <- t(as.data.frame(as.matrix(otu_table(sbset.gn))))
specpool(new)

#build the species accumulation curve & rarefaction curve (expected)
aurora.specaccum <- specaccum(new[,],method = "rarefaction")
#build a expected curve (randomization for boxplot comparison)
aurora.specaccum.rand <- specaccum(new[,], "random")
#plot both curves ("observed" vs "randomized")
jpeg("./accumulation_curves/Pygmy_Ireland_accumulation_species.jpeg", width = 600, height = 450)
plot(aurora.specaccum.rand,ci.type="poly", 
     col="blue", lwd=2, ci.lty=0, ci.col="lightblue",
     xlab = "Number of Shrews", ylab = "Number of Taxa")
boxplot(aurora.specaccum.rand, col="yellow", add=TRUE, pch="+")

text(20, 300, "Chao = 712.589 +/- 76.11", cex = 1)
legend('bottomright', "S. minutus - Ireland\nSpecies\n")
dev.off()

sbset.gn <- tax_glom(sbset, taxrank = "genus")
new <- t(as.data.frame(as.matrix(otu_table(sbset.gn))))
specpool(new)

#build the species accumulation curve & rarefaction curve (expected)
aurora.specaccum <- specaccum(new[,],method = "rarefaction")
#build a expected curve (randomization for boxplot comparison)
aurora.specaccum.rand <- specaccum(new[,], "random")
#plot both curves ("observed" vs "randomized")
jpeg("./accumulation_curves/Pygmy_Ireland_accumulation_genus.jpeg", width = 600, height = 450)
plot(aurora.specaccum.rand,ci.type="poly", 
     col="blue", lwd=2, ci.lty=0, ci.col="lightblue",
     xlab = "Number of Shrews", ylab = "Number of Taxa")
boxplot(aurora.specaccum.rand, col="yellow", add=TRUE, pch="+")

text(20, 250, "Chao = 470.55 +/- 44.89", cex = 1)
legend('bottomright', "S. minutus - Ireland\nGenus\n")
dev.off()

# Family
sbset.fm <- tax_glom(sbset, taxrank = "family")
new <- t(as.data.frame(as.matrix(otu_table(sbset.fm))))
specpool(new)

#build the species accumulation curve & rarefaction curve (expected)
aurora.specaccum <- specaccum(new[,],method = "rarefaction")
#build a expected curve (randomization for boxplot comparison)
aurora.specaccum.rand <- specaccum(new[,], "random")
#plot both curves ("observed" vs "randomized")
jpeg("./accumulation_curves/Pygmy_Ireland_accumulation_family.jpeg", width = 600, height = 450)
plot(aurora.specaccum.rand,ci.type="poly", 
     col="blue", lwd=2, ci.lty=0, ci.col="lightblue",
     xlab = "Number of Shrews", ylab = "Number of Taxa")
boxplot(aurora.specaccum.rand, col="yellow", add=TRUE, pch="+")

text(20, 120, "Chao = 178.4 +/- 14.29", cex = 1)
legend('bottomright', "S. minutus - Ireland\nFamily\n")
dev.off()

# Order
sbset.or <- tax_glom(sbset, taxrank = "order")
new <- t(as.data.frame(as.matrix(otu_table(sbset.or))))
specpool(new)

#build the species accumulation curve & rarefaction curve (expected)
aurora.specaccum <- specaccum(new[,],method = "rarefaction")
#build a expected curve (randomization for boxplot comparison)
aurora.specaccum.rand <- specaccum(new[,], "random")
#plot both curves ("observed" vs "randomized")
jpeg("./accumulation_curves/Pygmy_Ireland_accumulation_order.jpeg", width = 600, height = 450)
plot(aurora.specaccum.rand,ci.type="poly", 
     col="blue", lwd=2, ci.lty=0, ci.col="lightblue",
     xlab = "Number of Shrews", ylab = "Number of Taxa")
boxplot(aurora.specaccum.rand, col="yellow", add=TRUE, pch="+")

text(18, 30, "Chao = 42.90 +/- 10.06", cex = 1)
legend('bottomright', "S. minutus - Ireland\nOrder\n")
dev.off()

#### Pygmy France #####

sbset <- subset_samples(final.diet, shrew_country == "Pygmy_France")
sbset <- prune_taxa(taxa_sums(sbset) > 0, sbset)

new <- t(as.data.frame(as.matrix(otu_table(sbset))))
specpool(new)

#build the species accumulation curve & rarefaction curve (expected)
aurora.specaccum <- specaccum(new[,],method = "rarefaction")
#build a expected curve (randomization for boxplot comparison)
aurora.specaccum.rand <- specaccum(new[,], "random")
#plot both curves ("observed" vs "randomized")
jpeg("./accumulation_curves/Pygmy_France_accumulation_MOTU.jpeg", width = 600, height = 450)
plot(aurora.specaccum.rand,ci.type="poly", 
     col="blue", lwd=2, ci.lty=0, ci.col="lightblue",
     xlab = "Number of Shrews", ylab = "Number of Taxa")
boxplot(aurora.specaccum.rand, col="yellow", add=TRUE, pch="+")

text(6, 210, "Chao = 576.08 +/- 90.03", cex = 1)
legend('bottomright', "S. minutus - Belle Ile\nMOTUs\n")
dev.off()

# Species
sbset.gn <- tax_glom(sbset, taxrank = "species")
new <- t(as.data.frame(as.matrix(otu_table(sbset.gn))))
specpool(new)

#build the species accumulation curve & rarefaction curve (expected)
aurora.specaccum <- specaccum(new[,],method = "rarefaction")
#build a expected curve (randomization for boxplot comparison)
aurora.specaccum.rand <- specaccum(new[,], "random")
#plot both curves ("observed" vs "randomized")
jpeg("./accumulation_curves/Pygmy_France_accumulation_species.jpeg", width = 600, height = 450)
plot(aurora.specaccum.rand,ci.type="poly", 
     col="blue", lwd=2, ci.lty=0, ci.col="lightblue",
     xlab = "Number of Shrews", ylab = "Number of Taxa")
boxplot(aurora.specaccum.rand, col="yellow", add=TRUE, pch="+")

text(6, 140, "Chao = 263.69 +/- 38.48", cex = 1)
legend('bottomright', "S. minutus - Belle Ile\nSpecies\n")
dev.off()

# Genus
sbset.gn <- tax_glom(sbset, taxrank = "genus")
new <- t(as.data.frame(as.matrix(otu_table(sbset.gn))))
specpool(new)

#build the species accumulation curve & rarefaction curve (expected)
aurora.specaccum <- specaccum(new[,],method = "rarefaction")
#build a expected curve (randomization for boxplot comparison)
aurora.specaccum.rand <- specaccum(new[,], "random")
#plot both curves ("observed" vs "randomized")
jpeg("./accumulation_curves/Pygmy_France_accumulation_genus.jpeg", width = 600, height = 450)
plot(aurora.specaccum.rand,ci.type="poly", 
     col="blue", lwd=2, ci.lty=0, ci.col="lightblue",
     xlab = "Number of Shrews", ylab = "Number of Taxa")
boxplot(aurora.specaccum.rand, col="yellow", add=TRUE, pch="+")

text(6, 140, "Chao = 263.69 +/- 38.48", cex = 1)
legend('bottomright', "S. minutus - Belle Ile\nGenus\n")
dev.off()

# Family
sbset.fm <- tax_glom(sbset, taxrank = "family")
new <- t(as.data.frame(as.matrix(otu_table(sbset.fm))))
specpool(new)

#build the species accumulation curve & rarefaction curve (expected)
aurora.specaccum <- specaccum(new[,],method = "rarefaction")
#build a expected curve (randomization for boxplot comparison)
aurora.specaccum.rand <- specaccum(new[,], "random")
#plot both curves ("observed" vs "randomized")
jpeg("./accumulation_curves/Pygmy_France_accumulation_family.jpeg", width = 600, height = 450)
plot(aurora.specaccum.rand,ci.type="poly", 
     col="blue", lwd=2, ci.lty=0, ci.col="lightblue",
     xlab = "Number of Shrews", ylab = "Number of Taxa")
boxplot(aurora.specaccum.rand, col="yellow", add=TRUE, pch="+")

text(6, 95, "Chao = 148.44 +/- 21.68", cex = 1)
legend('bottomright', "S. minutus - Belle Ile\nFamily\n")
dev.off()

# Order
sbset.or <- tax_glom(sbset, taxrank = "order")
new <- t(as.data.frame(as.matrix(otu_table(sbset.or))))
specpool(new)

#build the species accumulation curve & rarefaction curve (expected)
aurora.specaccum <- specaccum(new[,],method = "rarefaction")
#build a expected curve (randomization for boxplot comparison)
aurora.specaccum.rand <- specaccum(new[,], "random")
#plot both curves ("observed" vs "randomized")
jpeg("./accumulation_curves/Pygmy_France_accumulation_Order.jpeg", width = 600, height = 450)
plot(aurora.specaccum.rand,ci.type="poly", 
     col="blue", lwd=2, ci.lty=0, ci.col="lightblue",
     xlab = "Number of Shrews", ylab = "Number of Taxa")
boxplot(aurora.specaccum.rand, col="yellow", add=TRUE, pch="+")

text(6, 31, "Chao = 39.64 +/- 8.09", cex = 1)
legend('bottomright', "S. minutus - Belle Ile\nOrder\n")
dev.off()

#### GWTS Ireland #####

sbset <- subset_samples(final.diet, shrew_country == "GWTS_Ireland")
sbset <- prune_taxa(taxa_sums(sbset) > 0, sbset)

new <- t(as.data.frame(as.matrix(otu_table(sbset))))
specpool(new)

#build the species accumulation curve & rarefaction curve (expected)
aurora.specaccum <- specaccum(new[,],method = "rarefaction")
#build a expected curve (randomization for boxplot comparison)
aurora.specaccum.rand <- specaccum(new[,], "random")
#plot both curves ("observed" vs "randomized")
jpeg("./accumulation_curves/GWTS_Ireland_accumulation_MOTUs.jpeg", width = 600, height = 450)
plot(aurora.specaccum.rand,ci.type="poly", 
     col="blue", lwd=2, ci.lty=0, ci.col="lightblue",
     xlab = "Number of Shrews", ylab = "Number of Taxa")
boxplot(aurora.specaccum.rand, col="yellow", add=TRUE, pch="+")

text(7, 280, "Chao = 773.60 +/- 107.88", cex = 1)
legend('bottomright', "C. russula - Ireland\nMOTUs\n")
dev.off()

# Species
sbset.gn <- tax_glom(sbset, taxrank = "species")
new <- t(as.data.frame(as.matrix(otu_table(sbset.gn))))
specpool(new)

#build the species accumulation curve & rarefaction curve (expected)
aurora.specaccum <- specaccum(new[,],method = "rarefaction")
#build a expected curve (randomization for boxplot comparison)
aurora.specaccum.rand <- specaccum(new[,], "random")
#plot both curves ("observed" vs "randomized")
jpeg("./accumulation_curves/GWTS_Ireland_accumulation_species.jpeg", width = 600, height = 450)
plot(aurora.specaccum.rand,ci.type="poly", 
     col="blue", lwd=2, ci.lty=0, ci.col="lightblue",
     xlab = "Number of Shrews", ylab = "Number of Taxa")
boxplot(aurora.specaccum.rand, col="yellow", add=TRUE, pch="+")

text(7, 130, "Chao = 223.66 +/- 31.23", cex = 1)
legend('bottomright', "C. russula - Ireland\nSpecies\n")
dev.off()


# Genus
sbset.gn <- tax_glom(sbset, taxrank = "genus")
new <- t(as.data.frame(as.matrix(otu_table(sbset.gn))))
specpool(new)

#build the species accumulation curve & rarefaction curve (expected)
aurora.specaccum <- specaccum(new[,],method = "rarefaction")
#build a expected curve (randomization for boxplot comparison)
aurora.specaccum.rand <- specaccum(new[,], "random")
#plot both curves ("observed" vs "randomized")
jpeg("./accumulation_curves/GWTS_Ireland_accumulation_genus.jpeg", width = 600, height = 450)
plot(aurora.specaccum.rand,ci.type="poly", 
     col="blue", lwd=2, ci.lty=0, ci.col="lightblue",
     xlab = "Number of Shrews", ylab = "Number of Taxa")
boxplot(aurora.specaccum.rand, col="yellow", add=TRUE, pch="+")

text(7, 130, "Chao = 223.66 +/- 31.23", cex = 1)
legend('bottomright', "C. russula - Ireland\nGenus\n")
dev.off()

# Family
sbset.fm <- tax_glom(sbset, taxrank = "family")
new <- t(as.data.frame(as.matrix(otu_table(sbset.fm))))
specpool(new)

#build the species accumulation curve & rarefaction curve (expected)
aurora.specaccum <- specaccum(new[,],method = "rarefaction")
#build a expected curve (randomization for boxplot comparison)
aurora.specaccum.rand <- specaccum(new[,], "random")
#plot both curves ("observed" vs "randomized")
jpeg("./accumulation_curves/GWTS_Ireland_accumulation_family.jpeg", width = 600, height = 450)
plot(aurora.specaccum.rand,ci.type="poly", 
     col="blue", lwd=2, ci.lty=0, ci.col="lightblue",
     xlab = "Number of Shrews", ylab = "Number of Taxa")
boxplot(aurora.specaccum.rand, col="yellow", add=TRUE, pch="+")

text(7, 90, "Chao = 120.22 +/- 14.14", cex = 1)
legend('bottomright', "C. russula - Ireland\nFamily\n")
dev.off()

# Order
sbset.or <- tax_glom(sbset, taxrank = "order")
new <- t(as.data.frame(as.matrix(otu_table(sbset.or))))
specpool(new)

#build the species accumulation curve & rarefaction curve (expected)
aurora.specaccum <- specaccum(new[,],method = "rarefaction")
#build a expected curve (randomization for boxplot comparison)
aurora.specaccum.rand <- specaccum(new[,], "random")
#plot both curves ("observed" vs "randomized")
jpeg("./accumulation_curves/GWTS_Ireland_accumulation_order.jpeg", width = 600, height = 450)
plot(aurora.specaccum.rand,ci.type="poly", 
     col="blue", lwd=2, ci.lty=0, ci.col="lightblue",
     xlab = "Number of Shrews", ylab = "Number of Taxa")
boxplot(aurora.specaccum.rand, col="yellow", add=TRUE, pch="+")

text(7, 30, "Chao = 30.94 +/- 2.58", cex = 1)
legend('bottomright', "C. russula - Ireland\nOrder\n")
dev.off()

#### GWTS France ####

sbset <- subset_samples(final.diet, shrew_country == "GWTS_France")
sbset <- prune_taxa(taxa_sums(sbset) > 0, sbset)

new <- t(as.data.frame(as.matrix(otu_table(sbset))))
specpool(new)

#build the species accumulation curve & rarefaction curve (expected)
aurora.specaccum <- specaccum(new[,],method = "rarefaction")
#build a expected curve (randomization for boxplot comparison)
aurora.specaccum.rand <- specaccum(new[,], "random")
#plot both curves ("observed" vs "randomized")
jpeg("./accumulation_curves/GWTS_France_accumulation_MOTUs.jpeg", width = 600, height = 450)
plot(aurora.specaccum.rand,ci.type="poly", 
     col="blue", lwd=2, ci.lty=0, ci.col="lightblue",
     xlab = "Number of Shrews", ylab = "Number of Taxa")
boxplot(aurora.specaccum.rand, col="yellow", add=TRUE, pch="+")

text(5, 215, "Chao = 625.14 +/- 103.97", cex = 1)
legend('bottomright', "C. russula - Belle Ile\nMOTUs\n")
dev.off()

# Species
sbset.gn <- tax_glom(sbset, taxrank = "species")
new <- t(as.data.frame(as.matrix(otu_table(sbset.gn))))
specpool(new)

#build the species accumulation curve & rarefaction curve (expected)
aurora.specaccum <- specaccum(new[,],method = "rarefaction")
#build a expected curve (randomization for boxplot comparison)
aurora.specaccum.rand <- specaccum(new[,], "random")
#plot both curves ("observed" vs "randomized")
jpeg("./accumulation_curves/GWTS_France_accumulation_species.jpeg", width = 600, height = 450)
plot(aurora.specaccum.rand,ci.type="poly", 
     col="blue", lwd=2, ci.lty=0, ci.col="lightblue",
     xlab = "Number of Shrews", ylab = "Number of Taxa")
boxplot(aurora.specaccum.rand, col="yellow", add=TRUE, pch="+")

text(5, 105, "Chao = 164.79 +/- 20.61", cex = 1)
legend('bottomright', "C. russula - Belle Ile\nSpecies\n")
dev.off()


# Genus
sbset.gn <- tax_glom(sbset, taxrank = "genus")
new <- t(as.data.frame(as.matrix(otu_table(sbset.gn))))
specpool(new)

#build the species accumulation curve & rarefaction curve (expected)
aurora.specaccum <- specaccum(new[,],method = "rarefaction")
#build a expected curve (randomization for boxplot comparison)
aurora.specaccum.rand <- specaccum(new[,], "random")
#plot both curves ("observed" vs "randomized")
jpeg("./accumulation_curves/GWTS_France_accumulation_genus.jpeg", width = 600, height = 450)
plot(aurora.specaccum.rand,ci.type="poly", 
     col="blue", lwd=2, ci.lty=0, ci.col="lightblue",
     xlab = "Number of Shrews", ylab = "Number of Taxa")
boxplot(aurora.specaccum.rand, col="yellow", add=TRUE, pch="+")

text(5, 105, "Chao = 164.79 +/- 20.61", cex = 1)
legend('bottomright', "C. russula - Belle Ile\nGenus\n")
dev.off()

# Family
sbset.fm <- tax_glom(sbset, taxrank = "family")
new <- t(as.data.frame(as.matrix(otu_table(sbset.fm))))
specpool(new)

#build the species accumulation curve & rarefaction curve (expected)
aurora.specaccum <- specaccum(new[,],method = "rarefaction")
#build a expected curve (randomization for boxplot comparison)
aurora.specaccum.rand <- specaccum(new[,], "random")
#plot both curves ("observed" vs "randomized")
jpeg("./accumulation_curves/GWTS_France_accumulation_family.jpeg", width = 600, height = 450)
plot(aurora.specaccum.rand,ci.type="poly", 
     col="blue", lwd=2, ci.lty=0, ci.col="lightblue",
     xlab = "Number of Shrews", ylab = "Number of Taxa")
boxplot(aurora.specaccum.rand, col="yellow", add=TRUE, pch="+")

text(5, 85, "Chao = 129.10 +/- 19.19", cex = 1)
legend('bottomright', "C. russula - Belle Ile\nFamily\n")
dev.off()

# Order
sbset.or <- tax_glom(sbset, taxrank = "order")
new <- t(as.data.frame(as.matrix(otu_table(sbset.or))))
specpool(new)

#build the species accumulation curve & rarefaction curve (expected)
aurora.specaccum <- specaccum(new[,],method = "rarefaction")
#build a expected curve (randomization for boxplot comparison)
aurora.specaccum.rand <- specaccum(new[,], "random")
#plot both curves ("observed" vs "randomized")
jpeg("./accumulation_curves/GWTS_France_accumulation_order.jpeg", width = 600, height = 450)
plot(aurora.specaccum.rand,ci.type="poly", 
     col="blue", lwd=2, ci.lty=0, ci.col="lightblue",
     xlab = "Number of Shrews", ylab = "Number of Taxa")
boxplot(aurora.specaccum.rand, col="yellow", add=TRUE, pch="+")

text(5, 30, "Chao = 52.48 +/- 29.83", cex = 1)
legend('bottomright', "C. russula - Belle Ile\nOrder\n")
dev.off()

par(c(5.1, 4.1, 4.1, 2.1))


#### Alpha Diversity: Individuals ####

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


paths <- c(b.meth1 = "./browett_et_al_2021_assignment_restrictions/filter_method_1/taxonomically_assigned_filtered.RData",
           b.meth2 = "./browett_et_al_2021_assignment_restrictions/filter_method_2/taxonomically_assigned_filtered.RData",
           a.meth1 = "./alberdi_taxa_assignment_restrictions/filter_method_1/taxonomically_assigned_filtered.RData",
           a.meth2 = "./alberdi_taxa_assignment_restrictions/filter_method_2/taxonomically_assigned_filtered.RData")

plot_list1 <- list()
plot_list2 <- list()

ave.df <- data.frame(group = character(), method = character(), 
                     mean.richness = character(), mean.rich.sd = character(),
                     mean.eveness = character(), mean.even.sd = character(), stringsAsFactors = FALSE)
stat_list <- list()
posthoc_test <- list()
ovrlp.df <- data.frame(group = character(), method = character(), richness = character(),
                       shannon = character(), levins = character(), stnd.levins = character(), stringsAsFactors = FALSE)


length(paths)

stat_list_rich <- list()
stat_list_shan <- list()
posthoc_test_rich <- list()
posthoc_test_shan <- list()

for(i in 1:length(paths)){
  
  print(paste("running alpha diversity for group", i, "of", length(paths), sep = ""))
  
  mt.nam <- names(paths[i])
  load(paths[i])
  
  sep.gn <- final.diet
  
  min_lib <- min(sample_sums(sep.gn)) # establish smallest sampling depth
  nsamp = nsamples(sep.gn)
  trials = 100 # how many trials to run
  richness <- matrix(nrow = nsamp, ncol = trials)
  row.names(richness) <- sample_names(sep.gn)
  evenness <- matrix(nrow = nsamp, ncol = trials)
  row.names(evenness) <- sample_names(sep.gn)
  set.seed(3)
  
  # loop to rerun rarifying and diversity estimates
  for (j in 1:trials) {
    #Subsample
    r <- rarefy_even_depth(sep.gn, sample.size = min_lib, verbose = FALSE, replace = TRUE)
    
    #Calculate richness
    rich <- as.numeric(as.matrix(estimate_richness(r, measures = "Observed")))
    richness[ ,j] <- rich
    
    #Calculate evenness
    even <- as.numeric(as.matrix(estimate_richness(r, measures = "Shannon")))
    evenness[ ,j] <- even
  }
  
  # Combine diversity estimates with Sample Data
  data <- as.data.frame(as.matrix(sample_data(sep.gn)))
  Sample <- row.names(richness)
  rich.mean <- apply(richness, 1, mean)
  rich.sd <- apply(richness, 1, sd)
  #measure <- rep("Richness", nsamp)
  even.mean <- apply(evenness, 1, mean)
  even.sd <- apply(evenness, 1, sd)
  #rich_stats <- data.frame(Sample, mean, sd, measure)
  rich_stats <- cbind(data, rich.mean, rich.sd, even.mean, even.sd)
  
  
  G.1 <- ggplot(rich_stats, aes(Zone, rich.mean)) +
    geom_boxplot(fill= c("#440154ff", '#21908dff', "indianred", "#fde725ff", "violetred2", "#8f7c00"), alpha = 0.25) +
    geom_jitter(aes(colour = Shrew), position = position_jitter(width = 0.1), size=3) +
    stat_summary(fun.y = mean, geom="point", shape=20, size = 10) +
    # scale_colour_manual(values= pal3) + # Manually sets colours according to object "palette"
    #scale_x_discrete(labels= c("Inside", "Middle", "Edge")) + # Manually labels groups using set described in "key" object
    # scale_y_continuous(limits = c(0, 0.55), expand = c(0, 0)) + # Sets Graph Value Limits
    theme_classic() + # Chooses Theme
    theme(axis.text = element_text(size = 10, face= "bold"), # Size of Axis measurements/Values Text
          axis.title=element_text(size=14), # Sets size of Axis Titles. face= "bold" will change the font to bold
          panel.grid.major = element_line(colour = "grey95"), # Sets colour of major grid lines
          panel.grid.minor = element_line(colour = "grey95") # Sets colour of minor grid lines
    ) +
    labs(y= "Species\nrichness", x = "") + # Labels Axis'
    facet_wrap(~ Shrew, scale = "free_x") +
    theme(strip.text.x = element_text(size = 12, face = "bold.italic")) +
    theme(legend.position="none")
  
  
  G.2 <- ggplot(rich_stats, aes(Zone, even.mean)) +
    geom_boxplot(fill= c("#440154ff", '#21908dff', "indianred", "#fde725ff", "violetred2", "#8f7c00"), alpha = 0.25) +
    geom_jitter(aes(colour = Shrew), position = position_jitter(width = 0.1), size=3) +
    stat_summary(fun.y = mean, geom="point", shape=20, size = 10) +
    # scale_colour_manual(values= pal3) + # Manually sets colours according to object "palette"
    #scale_x_discrete(labels= c("Inside", "Middle", "Edge")) + # Manually labels groups using set described in "key" object
    # scale_y_continuous(limits = c(0, 0.55), expand = c(0, 0)) + # Sets Graph Value Limits
    theme_classic() + # Chooses Theme
    theme(axis.text = element_text(size = 10, face= "bold"), # Size of Axis measurements/Values Text
          axis.title=element_text(size=14), # Sets size of Axis Titles. face= "bold" will change the font to bold
          panel.grid.major = element_line(colour = "grey95"), # Sets colour of major grid lines
          panel.grid.minor = element_line(colour = "grey95") # Sets colour of minor grid lines
    ) +
    labs(y= "Shannon\nDiversity Index", x = "") + # Labels Axis'
    facet_wrap(~ Shrew, scale = "free_x") +
    theme(strip.text.x = element_text(size = 12, face = "bold.italic")) +
    theme(legend.position="none")
  
  plot_list1[[mt.nam]] <- G.1
  plot_list2[[mt.nam]] <- G.2
  
  ave <- matrix(nrow = 4,
                ncol = 6)
  ave[,1] <- levels(rich_stats$shrew_country)
  ave[,2] <- names(paths)[i]
  ave[,3] <- aggregate(rich_stats$rich.mean, list(rich_stats$shrew_country), mean)[,2]
  ave[,4] <- aggregate(rich_stats$rich.mean, list(rich_stats$shrew_country), sd)[,2]
  ave[,5] <- aggregate(rich_stats$even.mean, list(rich_stats$shrew_country), mean)[,2]
  ave[,6] <- aggregate(rich_stats$even.mean, list(rich_stats$shrew_country), sd)[,2]
  
  colnames(ave) <- colnames(ave.df)
  ave.df <- rbind(ave.df, as.data.frame(ave))
  
  stat_list_rich[[mt.nam]] <- kruskal.test(rich.mean ~ shrew_zone, data = rich_stats)  
  stat_list_shan[[mt.nam]] <- kruskal.test(even.mean ~ shrew_zone, data = rich_stats)  
  
  posthoc_test_rich[[mt.nam]] <- pairwise.wilcox.test(rich_stats$rich.mean,
                                                      rich_stats$shrew_zone,
                                                      p.adjust.method = "BH")
  posthoc_test_shan[[mt.nam]] <- pairwise.wilcox.test(rich_stats$even.mean,
                                                      rich_stats$shrew_zone,
                                                      p.adjust.method = "BH")
  
  
}

beep("coin")

grid.arrange(plot_list1[["b.meth1"]],
             plot_list2[["b.meth1"]])

grid.arrange(plot_list1[["b.meth2"]],
             plot_list2[["b.meth2"]])

grid.arrange(plot_list1[["a.meth1"]],
             plot_list2[["a.meth1"]])

grid.arrange(plot_list1[["a.meth2"]],
             plot_list2[["a.meth2"]])

stat_list_rich[[4]]
stat_list_shan[[4]]
posthoc_test_rich[[4]]
posthoc_test_shan[[4]]

save(plot_list1, plot_list2, stat_list_rich, stat_list_shan, posthoc_test_rich, posthoc_test_shan,
     file = "./alpha_diversity_individual_plots_stats.RData")

#### Alpha Diversity: Group: RRA metric ####

vr <- "shrew_zone"

final.diet.ra = transform_sample_counts(final.diet, function(x) x/sum(x))
otudf <- as.matrix(otu_table(final.diet.ra))
colnames(otudf) <- as.character(as.matrix(sample_data(final.diet.ra)[,vr]))
grps <- unique(colnames(otudf))

rra.tab <- t(otudf)
rra.tab[1:5,1:5]
#rra.tab <- apply(array, margin, ...)
rra.tab <- aggregate(rra.tab, list(rownames(rra.tab)), mean)
rra.tab[,1:5]
rownames(rra.tab) <- rra.tab$Group.1
rra.tab <- rra.tab[,-1]
rra.tab[,1:5]
grps <- sort(unique(rownames(rra.tab)))

pootab <- t(rra.tab)
ovrlp <- matrix(nrow = 5,
                ncol = length(grps))
ovrlp[1,] <- table(sample_data(final.diet)[,vr])
ovrlp[2,] <- specnumber(pootab, MARGIN = 2)
ovrlp[3,] <- as.numeric(niche.width(pootab, method = "shannon"))  
ovrlp[4,] <- as.numeric(niche.width(pootab, method = "levins"))
for (i in 1:ncol(pootab)){
  B <- ovrlp[4,i]
  n <- ovrlp[2,i]
  Ba <- (B - 1) / (n - 1)
  ovrlp[5,i] <- Ba
}


rownames(ovrlp) <- c("sample size", "species richness", "shannon", "levins", "stnd.levins")  
colnames(ovrlp) <- grps
ovrlp.fin <- t(as.data.frame(ovrlp))
ovrlp.fin <- cbind(glom = "MOTU", ovrlp.fin)

# - glom different levels - #

lv <- c("species", "genus", "family", "order")

for(m in 1:length(lv)){
  
  glom.diet <- tax_glom(final.diet, taxrank = lv[m])
  
  final.diet.ra = transform_sample_counts(glom.diet, function(x) x/sum(x))
  otudf <- as.matrix(otu_table(final.diet.ra))
  colnames(otudf) <- as.character(as.matrix(sample_data(final.diet.ra)[,vr]))
  grps <- unique(colnames(otudf))
  
  rra.tab <- t(otudf)
  rra.tab[1:5,1:5]
  #rra.tab <- apply(array, margin, ...)
  rra.tab <- aggregate(rra.tab, list(rownames(rra.tab)), mean)
  rra.tab[,1:5]
  rownames(rra.tab) <- rra.tab$Group.1
  rra.tab <- rra.tab[,-1]
  rra.tab[,1:5]
  grps <- sort(unique(rownames(rra.tab)))
  
  pootab <- t(rra.tab)
  ovrlp <- matrix(nrow = 5,
                  ncol = length(grps))
  ovrlp[1,] <- table(sample_data(final.diet)[,vr])
  ovrlp[2,] <- specnumber(pootab, MARGIN = 2)
  ovrlp[3,] <- as.numeric(niche.width(pootab, method = "shannon"))  
  ovrlp[4,] <- as.numeric(niche.width(pootab, method = "levins"))
  for (i in 1:ncol(pootab)){
    B <- ovrlp[4,i]
    n <- ovrlp[2,i]
    Ba <- (B - 1) / (n - 1)
    ovrlp[5,i] <- Ba
  }
  
  
  rownames(ovrlp) <- c("sample size", "species richness", "shannon", "levins", "stnd.levins")  
  colnames(ovrlp) <- grps
  ovrlp <- t(as.data.frame(ovrlp))
  ovrlp <- cbind(glom = lv[m], ovrlp)
  ovrlp.fin <- rbind(ovrlp.fin, ovrlp)
  
}


write.csv(ovrlp.fin, file = paste("./RRA_metric_", vr, "_group_alpha_diversity.csv", sep = ""))


#### Alpha Diversity: Group: POO metric ####

vr <- "shrew_country"

otudf <- as.data.frame(as.matrix(otu_table(final.diet)))
otudf <- as.matrix(otu_table(final.diet))
#colnames(otudf) <- sample_data(final.diet)$shrew_zone
colnames(otudf) <- as.character(as.matrix(sample_data(final.diet)[,vr]))
#colnames(otudf) <- as.character(unlist(sample_data(final.diet)[,vr])) # same as previous line
grps <- sort(unique(colnames(otudf)))

BBA_Dprop=prop.table(otudf,2)
dim(BBA_Dprop)
apply(BBA_Dprop,2,sum)

# Using a threshold of 0.1% to declare presence/absence 
PA_BBA_D= BBA_Dprop
PA_BBA_D[PA_BBA_D>.00099999999]=1 # 0.01 = 1%
PA_BBA_D[PA_BBA_D<.001]=0 # remember to change this too
apply(PA_BBA_D,2,sum)

FOO_All_BBA_=apply(PA_BBA_D,1,sum)/length(PA_BBA_D[1,]) * 100
Ord_BBA_=order(-FOO_All_BBA_) #

grp.list <- matrix(nrow = length(grps),
                   ncol = length(Ord_BBA_))

grp.list <- as.data.frame(grp.list)
colnames(grp.list) <- rownames(PA_BBA_D[Ord_BBA_,1:2])
rownames(grp.list) <- grps


for(i in 1:length(grps)){
  
  POO_= apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == grps[i])],1,sum)/sum(apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == grps[i])],1,sum)) * 100
  grp.list[i,] <- POO_
  
}

grp.list[,1:5]

pootab <- t(grp.list)
ovrlp <- matrix(nrow = 5,
                ncol = length(grps))
ovrlp[1,] <- table(sample_data(final.diet)[,vr])
ovrlp[2,] <- specnumber(pootab, MARGIN = 2)
ovrlp[3,] <- as.numeric(niche.width(pootab, method = "shannon"))  
ovrlp[4,] <- as.numeric(niche.width(pootab, method = "levins"))
for (i in 1:ncol(pootab)){
  B <- ovrlp[4,i]
  n <- ovrlp[2,i]
  Ba <- (B - 1) / (n - 1)
  ovrlp[5,i] <- Ba
}


rownames(ovrlp) <- c("sample size", "species richness", "shannon", "levins", "stnd.levins")  
colnames(ovrlp) <- grps
ovrlp.fin <- t(as.data.frame(ovrlp))
ovrlp.fin <- cbind(glom = "MOTU", ovrlp.fin)

# - glom different levels - #

lv <- c("species", "genus", "family", "order")

for(m in 1:length(lv)){
  glom.diet <- tax_glom(final.diet, taxrank = lv[m])
  
  
  #otudf <- as.data.frame(as.matrix(otu_table(glom.diet)))
  otudf <- as.matrix(otu_table(glom.diet))
  #colnames(otudf) <- sample_data(final.diet)$shrew_zone
  colnames(otudf) <- as.character(as.matrix(sample_data(glom.diet)[,vr]))
  #colnames(otudf) <- as.character(unlist(sample_data(final.diet)[,vr])) # same as previous line
  grps <- sort(unique(colnames(otudf)))
  
  BBA_Dprop=prop.table(otudf,2)
  dim(BBA_Dprop)
  apply(BBA_Dprop,2,sum)
  
  # Using a threshold of 1% to declare presence/absence 
  PA_BBA_D= BBA_Dprop
  PA_BBA_D[PA_BBA_D>.00099999999]=1 # 0.01 = 1%
  PA_BBA_D[PA_BBA_D<.001]=0 # remember to change this too
  apply(PA_BBA_D,2,sum)
  
  FOO_All_BBA_=apply(PA_BBA_D,1,sum)/length(PA_BBA_D[1,]) * 100
  Ord_BBA_=order(-FOO_All_BBA_) #
  
  grp.list <- matrix(nrow = length(grps),
                     ncol = length(Ord_BBA_))
  
  grp.list <- as.data.frame(grp.list)
  colnames(grp.list) <- rownames(PA_BBA_D[Ord_BBA_,1:2])
  rownames(grp.list) <- grps
  
  
  for(i in 1:length(grps)){
    
    POO_= apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == grps[i])],1,sum)/sum(apply(PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == grps[i])],1,sum)) * 100
    grp.list[i,] <- POO_
    
  }
  
  grp.list[,1:5]
  
  pootab <- t(grp.list)
  ovrlp <- matrix(nrow = 5,
                  ncol = length(grps))
  ovrlp[1,] <- table(sample_data(final.diet)[,vr])
  ovrlp[2,] <- specnumber(pootab, MARGIN = 2)
  ovrlp[3,] <- as.numeric(niche.width(pootab, method = "shannon"))  
  ovrlp[4,] <- as.numeric(niche.width(pootab, method = "levins"))
  for (i in 1:ncol(pootab)){
    B <- ovrlp[4,i]
    n <- ovrlp[2,i]
    Ba <- (B - 1) / (n - 1)
    ovrlp[5,i] <- Ba
  }
  
  rownames(ovrlp) <- c("sample size", "species richness", "shannon", "levins", "stnd.levins")  
  colnames(ovrlp) <- grps
  ovrlp <- t(as.data.frame(ovrlp))
  ovrlp <- cbind(glom = lv[m], ovrlp)
  ovrlp.fin <- rbind(ovrlp.fin, ovrlp)
  
}

ovrlp.fin

#write.csv(ovrlp.fin, file = paste("./POO_metric_", vr, "_group_alpha_diversity.csv", sep = ""))
write.csv(ovrlp.fin, file = paste("./POO_metric_", vr, "_group_alpha_diversity_01_percent_presence_rate.csv", sep = ""))

#### Alpha Diversity: Group: Subsampled: POO metric ####

#final.diet <- subset_samples(final.diet, Shrew == "Pygmy")

vr <- "shrew_country"
#n <- length(table(sample_data(final.diet)[,vr]))
min.sz <- min(table(sample_data(final.diet)[,vr]))
#grps <- sort(unique(sample_data(final.diet)[,"shrew_zone"]))
grps <- sort(unique(as.matrix(sample_data(final.diet)[,vr])))
lv <- c("species", "genus", "family", "order")
trials = 50

rich <- matrix(nrow = length(grps) * (length(lv) + 1), ncol = trials)
shan <- matrix(nrow = length(grps) * (length(lv) + 1), ncol = trials)
lev <- matrix(nrow = length(grps) * (length(lv) + 1), ncol = trials)
st.lev <- matrix(nrow = length(grps) * (length(lv) + 1), ncol = trials)

set.seed(1)

for(k in 1:trials){
  
  print(paste("running subsample ", k, " of ", trials))
  otudf <- as.data.frame(as.matrix(otu_table(final.diet)))
  otudf <- as.matrix(otu_table(final.diet))
  #colnames(otudf) <- sample_data(final.diet)$shrew_zone
  colnames(otudf) <- as.character(as.matrix(sample_data(final.diet)[,vr]))
  #colnames(otudf) <- as.character(unlist(sample_data(final.diet)[,vr])) # same as previous line
  grps <- sort(unique(colnames(otudf)))
  
  BBA_Dprop=prop.table(otudf,2)
  dim(BBA_Dprop)
  apply(BBA_Dprop,2,sum)
  
  # Using a threshold of 1% to declare presence/absence 
  PA_BBA_D= BBA_Dprop
  PA_BBA_D[PA_BBA_D>.00099999999]=1 # 0.01 = 1%
  PA_BBA_D[PA_BBA_D<.001]=0 # Change here too
  apply(PA_BBA_D,2,sum)
  
  FOO_All_BBA_=apply(PA_BBA_D,1,sum)/length(PA_BBA_D[1,]) * 100
  Ord_BBA_=order(-FOO_All_BBA_) #
  
  grp.list <- matrix(nrow = length(grps),
                     ncol = length(Ord_BBA_))
  
  grp.list <- as.data.frame(grp.list)
  colnames(grp.list) <- rownames(PA_BBA_D[Ord_BBA_,1:2])
  rownames(grp.list) <- grps
  
  
  for(i in 1:length(grps)){
    
    POO_n <- PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == grps[i])]
    colnames(POO_n) <- 1:ncol(POO_n)
    POO_n <- POO_n[,sample(colnames(POO_n), min.sz, replace = FALSE)]
    POO_ <- apply(POO_n,1,sum)/sum(apply(POO_n,1,sum)) * 100
    grp.list[i,] <- POO_
    
    
  }
  
  grp.list[,1:5]
  
  pootab <- t(grp.list)
  ovrlp <- matrix(nrow = 5,
                  ncol = length(grps))
  ovrlp[1,] <- min.sz
  ovrlp[2,] <- specnumber(pootab, MARGIN = 2)
  ovrlp[3,] <- as.numeric(niche.width(pootab, method = "shannon"))  
  ovrlp[4,] <- as.numeric(niche.width(pootab, method = "levins"))
  for (i in 1:ncol(pootab)){
    B <- ovrlp[4,i]
    n <- ovrlp[2,i]
    Ba <- (B - 1) / (n - 1)
    ovrlp[5,i] <- Ba
  }
  
  
  rownames(ovrlp) <- c("sample size", "species richness", "shannon", "levins", "stnd.levins")  
  colnames(ovrlp) <- grps
  ovrlp.fin <- t(as.data.frame(ovrlp))
  ovrlp.fin <- cbind(glom = "MOTU", ovrlp.fin)
  
  # - glom different levels - #
  
  lv <- c("species", "genus", "family", "order")
  
  for(m in 1:length(lv)){
    glom.diet <- tax_glom(final.diet, taxrank = lv[m])
    
    
    #otudf <- as.data.frame(as.matrix(otu_table(glom.diet)))
    otudf <- as.matrix(otu_table(glom.diet))
    #colnames(otudf) <- sample_data(final.diet)$shrew_zone
    colnames(otudf) <- as.character(as.matrix(sample_data(glom.diet)[,vr]))
    #colnames(otudf) <- as.character(unlist(sample_data(final.diet)[,vr])) # same as previous line
    grps <- sort(unique(colnames(otudf)))
    
    BBA_Dprop=prop.table(otudf,2)
    dim(BBA_Dprop)
    apply(BBA_Dprop,2,sum)
    
    # Using a threshold of 0.01% to declare presence/absence 
    PA_BBA_D= BBA_Dprop
    PA_BBA_D[PA_BBA_D>.00099999999]=1 # 0.01 is 1%
    PA_BBA_D[PA_BBA_D<.001]=0 # Change here too
    apply(PA_BBA_D,2,sum)
    
    FOO_All_BBA_=apply(PA_BBA_D,1,sum)/length(PA_BBA_D[1,]) * 100
    Ord_BBA_=order(-FOO_All_BBA_) #
    
    grp.list <- matrix(nrow = length(grps),
                       ncol = length(Ord_BBA_))
    
    grp.list <- as.data.frame(grp.list)
    colnames(grp.list) <- rownames(PA_BBA_D[Ord_BBA_,1:2])
    rownames(grp.list) <- grps
    
    
    for(i in 1:length(grps)){
      
      POO_n <- PA_BBA_D[Ord_BBA_,which(colnames(PA_BBA_D) == grps[i])]
      colnames(POO_n) <- 1:ncol(POO_n)
      POO_n <- POO_n[,sample(colnames(POO_n), min.sz, replace = FALSE)]
      POO_ <- apply(POO_n,1,sum)/sum(apply(POO_n,1,sum)) * 100
      grp.list[i,] <- POO_
      
    }
    
    grp.list[,1:5]
    
    pootab <- t(grp.list)
    ovrlp <- matrix(nrow = 5,
                    ncol = length(grps))
    ovrlp[1,] <- min.sz
    ovrlp[2,] <- specnumber(pootab, MARGIN = 2)
    ovrlp[3,] <- as.numeric(niche.width(pootab, method = "shannon"))  
    ovrlp[4,] <- as.numeric(niche.width(pootab, method = "levins"))
    for (i in 1:ncol(pootab)){
      B <- ovrlp[4,i]
      n <- ovrlp[2,i]
      Ba <- (B - 1) / (n - 1)
      ovrlp[5,i] <- Ba
    }
    
    rownames(ovrlp) <- c("sample size", "species richness", "shannon", "levins", "stnd.levins")  
    colnames(ovrlp) <- grps
    ovrlp <- t(as.data.frame(ovrlp))
    ovrlp <- cbind(glom = lv[m], ovrlp)
    ovrlp.fin <- rbind(ovrlp.fin, ovrlp)
    
  }
  
  #ovrlp.fin
  
  rich[,k] <- as.numeric(ovrlp.fin[,"species richness"])
  shan[,k] <- as.numeric(ovrlp.fin[,"shannon"])
  lev[,k] <- as.numeric(ovrlp.fin[,"levins"])
  st.lev[,k] <- as.numeric(ovrlp.fin[,"stnd.levins"])
  
}

ovrlp.fin.ave <- ovrlp.fin[,1:2]
rich.mean <- apply(rich, 1, mean)
rich.sd <- apply(rich, 1, sd)
shan.mean <- apply(shan, 1, mean)
shan.sd <- apply(shan, 1, sd)
lev.mean <- apply(lev, 1, mean)
lev.sd <- apply(lev, 1, sd)
st.lev.mean <- apply(st.lev, 1, mean)
st.lev.sd <- apply(st.lev, 1, sd)

#rich_stats <- data.frame(Sample, mean, sd, measure)
#rich_stats <- cbind(data, rich.mean, rich.sd, even.mean, even.sd)
ovrlp.fin.ave <- cbind(ovrlp.fin.ave, rich.mean, rich.sd,
                       shan.mean, shan.sd,
                       lev.mean, lev.sd,
                       st.lev.mean, st.lev.sd)

write.csv(ovrlp.fin.ave, file = paste("./POO_metric_", vr, "subsampled_group_alpha_diversity_01_percent_presence_rate.csv", sep = ""))
#write.csv(ovrlp.fin.ave, file = paste("./POO_metric_", vr, "subsampled_group_alpha_diversity_01_percent_presence_rate_Sorex_only.csv", sep = ""))
