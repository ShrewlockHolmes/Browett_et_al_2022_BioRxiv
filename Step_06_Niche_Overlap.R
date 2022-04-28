#library(EcoSimR)
library(gtools)

dir.create("./ecosimr_output")

source("./ecosimr/EcoSimR - Main Source.R")

my.params <- Param.List
my.params$Plot.Output <- "file"
my.params$Print.Output <- "file"
my.params$N.Reps <- 10000

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

# Browett Method 1
load("./browett_et_al_2021_assignment_restrictions/filter_method_1/taxonomically_assigned_filtered.RData")
dir.create("./browett_et_al_2021_assignment_restrictions/filter_method_1/ecosimr_output")
dir.path <- "./browett_et_al_2021_assignment_restrictions/filter_method_1/ecosimr_output/"

# Browett Method 2
load("./browett_et_al_2021_assignment_restrictions/filter_method_2/taxonomically_assigned_filtered.RData")
dir.create("./browett_et_al_2021_assignment_restrictions/filter_method_2/ecosimr_output")
dir.path <- "./browett_et_al_2021_assignment_restrictions/filter_method_2/ecosimr_output/"

# Alberdi Method 1
load("./alberdi_taxa_assignment_restrictions/filter_method_1/taxonomically_assigned_filtered.RData")

# Alberdi Method 2
load("./alberdi_taxa_assignment_restrictions/filter_method_2/taxonomically_assigned_filtered.RData")


dir.create("./ecosimr_output")

source("./ecosimr/EcoSimR - Main Source.R")

my.params <- Param.List

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

# Browett Method 1
load("./browett_et_al_2021_assignment_restrictions/filter_method_1/taxonomically_assigned_filtered.RData")

paths <- c(b.meth1 = "./browett_et_al_2021_assignment_restrictions/filter_method_1/taxonomically_assigned_filtered.RData",
           b.meth2 = "./browett_et_al_2021_assignment_restrictions/filter_method_2/taxonomically_assigned_filtered.RData",
           a.meth1 = "./alberdi_taxa_assignment_restrictions/filter_method_1/taxonomically_assigned_filtered.RData",
           a.meth2 = "./alberdi_taxa_assignment_restrictions/filter_method_2/taxonomically_assigned_filtered.RData")

load(paths[4])



#### Automate ecosimr: RRA Metric ####

vr <- "shrew_country"

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
grps <- unique(rownames(rra.tab))

library(gtools)
library(stringr)
chck <- combinations(n=length(grps),r=2,v=grps,repeats.allowed=FALSE)

chck[1,1]
chck[1,2]

inpt <- nrow(as.data.frame(chck[,]))

df <- data.frame(group_1 = 1:inpt,
                 group_2 = NA,
                 Seed = NA,
                 Num_of_reps = NA,
                 Elapsed_time = NA,
                 Metric = NA,
                 Algorithm = NA,
                 Observed = NA,
                 Mean_sim_index = NA,
                 Var_sim_index = NA,
                 low_95_1tail = NA,
                 upper_95_1tail = NA,
                 low_95_2tail = NA,
                 upper_95_2tail = NA,
                 P_obs_less_Null = NA,
                 P_obs_great_Null = NA,
                 P_obs_equal_Null = NA)

for (i in 1:inpt){
  
  chck1 <- chck[i,1]
  chck2 <- chck[i,2]
  pooeded <- rbind(rra.tab[chck1,],
                   rra.tab[chck2,])
  pooeded <- pooeded[,which(apply(pooeded,2,sum) > 0)]
  write.csv(pooeded, file = "file_for_ecosimr.csv")
  my.params$Data.File <- "file_for_ecosimr.csv"
  
  Output.Results(my.params, Null.Model.Engine(my.params))
  
  x1 <- read.delim(file = "./Niche Overlap Output.txt", header = FALSE, sep = "\t")
  x1$V1 <- str_replace(x1$V1, " =  ", ":")
  x1$V1 <- str_replace(x1$V1, " > ", ":")
  x1$V1 <- str_replace(x1$V1, " < ", ":")
  x2 <- as.matrix(str_split(x1$V1, ":", simplify = TRUE))
  x2 <- as.matrix(x2[2:18,1:2])
  df1 <- as.data.frame(t(x2))
  #colnames(df1) <- colnames(df)
  #df <- rbind(df, df1[2,])
  df[i,] <- as.character(unlist(df1[2,]))
  df[i,1] <- chck1
  df[i,2] <- chck2
  
}

beep("treasure")

df <- cbind(glom = "MOTU", df)
df$glom <- as.character(df$glom)

#write.csv(df, file = paste("./", vr, "_ecosimr_results.csv", sep = ""))

# - # Agglom

# MOTU level: All MOTUs

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
  grps <- unique(rownames(rra.tab))
  
  #library(gtools)
  chck <- combinations(n=length(grps),r=2,v=grps,repeats.allowed=FALSE)
  
  chck[1,1]
  chck[1,2]
  
  inpt <- nrow(as.data.frame(chck[,]))
  
  for (i in 1:inpt){
    
    chck1 <- chck[i,1]
    chck2 <- chck[i,2]
    pooeded <- rbind(rra.tab[chck1,],
                     rra.tab[chck2,])
    pooeded <- pooeded[,which(apply(pooeded,2,sum) > 0)]
    write.csv(pooeded, file = "file_for_ecosimr.csv")
    my.params$Data.File <- "file_for_ecosimr.csv"
    
    Output.Results(my.params, Null.Model.Engine(my.params))
    
    x1 <- read.delim(file = "./Niche Overlap Output.txt", header = FALSE, sep = "\t")
    x1$V1 <- str_replace(x1$V1, " =  ", ":")
    x1$V1 <- str_replace(x1$V1, " > ", ":")
    x1$V1 <- str_replace(x1$V1, " < ", ":")
    x2 <- as.matrix(str_split(x1$V1, ":", simplify = TRUE))
    x2 <- as.matrix(x2[2:18,1:2])
    df1 <- as.data.frame(t(x2))
    #colnames(df1) <- colnames(df)
    #df <- rbind(df, df1[2,])
    n <- as.character(c(lv[m], chck1, chck2, as.character(unlist(df1[2,3:17]))))
    df <- rbind(df, n)
    #df.g[i,] <- as.character(unlist(df1[2,]))
    #df.g[i,1] <- chck1
    #df.g[i,2] <- chck2
    
  }
  
}

beep("treasure")

write.csv(df, file = paste("./", vr, "_ecosimr_results_RRA.csv", sep = ""))



#### Automate ecosimr: RRA Metric: Core ####

library(microbiome)

# Core Filter no. 1
filt.1 <- subset_samples(final.diet, shrew_country == "Pygmy_Ireland")
core.members.1 <- core_members(filt.1, detection = 0, prevalence = 2/nsamples(filt.1))
filt.2 <- subset_samples(final.diet, shrew_country == "Pygmy_France")
core.members.2 <- core_members(filt.2, detection = 0, prevalence = 2/nsamples(filt.2))
filt.3 <- subset_samples(final.diet, shrew_country == "GWTS_Ireland")
core.members.3 <- core_members(filt.3, detection = 0, prevalence = 2/nsamples(filt.3))
filt.4 <- subset_samples(final.diet, shrew_country == "GWTS_France")
core.members.4 <- core_members(filt.4, detection = 0, prevalence = 2/nsamples(filt.4))

core.members <- unique(c(core.members.1,core.members.2,core.members.3,core.members.4))

final.diet.cr <- prune_taxa(core.members, final.diet)

# Core Filter no. 2
core.members <- core_members(final.diet, detection = 0, prevalence = 2/nsamples(final.diet))

final.diet.cr <- prune_taxa(core.members, final.diet)


vr <- "shrew_country"

final.diet.ra = transform_sample_counts(final.diet.cr, function(x) x/sum(x))
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
grps <- unique(rownames(rra.tab))

library(gtools)
chck <- combinations(n=length(grps),r=2,v=grps,repeats.allowed=FALSE)

chck[1,1]
chck[1,2]

inpt <- nrow(as.data.frame(chck[,]))

df <- data.frame(group_1 = 1:inpt,
                 group_2 = NA,
                 Seed = NA,
                 Num_of_reps = NA,
                 Elapsed_time = NA,
                 Metric = NA,
                 Algorithm = NA,
                 Observed = NA,
                 Mean_sim_index = NA,
                 Var_sim_index = NA,
                 low_95_1tail = NA,
                 upper_95_1tail = NA,
                 low_95_2tail = NA,
                 upper_95_2tail = NA,
                 P_obs_less_Null = NA,
                 P_obs_great_Null = NA,
                 P_obs_equal_Null = NA)

for (i in 1:inpt){
  
  chck1 <- chck[i,1]
  chck2 <- chck[i,2]
  pooeded <- rbind(rra.tab[chck1,],
                   rra.tab[chck2,])
  pooeded <- pooeded[,which(apply(pooeded,2,sum) > 0)]
  write.csv(pooeded, file = "file_for_ecosimr.csv")
  my.params$Data.File <- "file_for_ecosimr.csv"
  
  Output.Results(my.params, Null.Model.Engine(my.params))
  
  x1 <- read.delim(file = "./Niche Overlap Output.txt", header = FALSE, sep = "\t")
  x1$V1 <- str_replace(x1$V1, " =  ", ":")
  x1$V1 <- str_replace(x1$V1, " > ", ":")
  x1$V1 <- str_replace(x1$V1, " < ", ":")
  x2 <- as.matrix(str_split(x1$V1, ":", simplify = TRUE))
  x2 <- as.matrix(x2[2:18,1:2])
  df1 <- as.data.frame(t(x2))
  #colnames(df1) <- colnames(df)
  #df <- rbind(df, df1[2,])
  df[i,] <- as.character(unlist(df1[2,]))
  df[i,1] <- chck1
  df[i,2] <- chck2
  
}

beep("treasure")

df <- cbind(glom = "MOTU", df)
df$glom <- as.character(df$glom)

#write.csv(df, file = paste("./", vr, "_ecosimr_results.csv", sep = ""))

# - # Agglom


lv <- c("species", "genus", "family", "order")

for(m in 1:length(lv)){
  
  glom.diet <- tax_glom(final.diet.cr, taxrank = lv[m])
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
  grps <- unique(rownames(rra.tab))
  
  #library(gtools)
  chck <- combinations(n=length(grps),r=2,v=grps,repeats.allowed=FALSE)
  
  chck[1,1]
  chck[1,2]
  
  inpt <- nrow(as.data.frame(chck[,]))
  
  for (i in 1:inpt){
    
    chck1 <- chck[i,1]
    chck2 <- chck[i,2]
    pooeded <- rbind(rra.tab[chck1,],
                     rra.tab[chck2,])
    pooeded <- pooeded[,which(apply(pooeded,2,sum) > 0)]
    write.csv(pooeded, file = "file_for_ecosimr.csv")
    my.params$Data.File <- "file_for_ecosimr.csv"
    
    Output.Results(my.params, Null.Model.Engine(my.params))
    
    x1 <- read.delim(file = "./Niche Overlap Output.txt", header = FALSE, sep = "\t")
    x1$V1 <- str_replace(x1$V1, " =  ", ":")
    x1$V1 <- str_replace(x1$V1, " > ", ":")
    x1$V1 <- str_replace(x1$V1, " < ", ":")
    x2 <- as.matrix(str_split(x1$V1, ":", simplify = TRUE))
    x2 <- as.matrix(x2[2:18,1:2])
    df1 <- as.data.frame(t(x2))
    #colnames(df1) <- colnames(df)
    #df <- rbind(df, df1[2,])
    n <- as.character(c(lv[m], chck1, chck2, as.character(unlist(df1[2,3:17]))))
    df <- rbind(df, n)
    #df.g[i,] <- as.character(unlist(df1[2,]))
    #df.g[i,1] <- chck1
    #df.g[i,2] <- chck2
    
  }
  
}

beep("treasure")

write.csv(df, file = paste("./", vr, "_ecosimr_results_RRA_Core.csv", sep = ""))





#### Automate ecosimr: RRA Metric: Interspecific ####

vr <- "shrew_country"

final.diet.ra = transform_sample_counts(final.diet, function(x) x/sum(x))
otudf <- as.matrix(otu_table(final.diet.ra))
colnames(otudf) <- as.character(as.matrix(sample_data(final.diet.ra)[,vr]))
grps <- unique(colnames(otudf))

rra.tab <- t(otudf)
rra.tab[1:5,1:5]

df <- data.frame(group_1 = 1:length(grps),
                 group_2 = "Interspecific",
                 Seed = NA,
                 Num_of_reps = NA,
                 Elapsed_time = NA,
                 Metric = NA,
                 Algorithm = NA,
                 Observed = NA,
                 Mean_sim_index = NA,
                 Var_sim_index = NA,
                 low_95_1tail = NA,
                 upper_95_1tail = NA,
                 low_95_2tail = NA,
                 upper_95_2tail = NA,
                 P_obs_less_Null = NA,
                 P_obs_great_Null = NA,
                 P_obs_equal_Null = NA)

for (i in 1:length(grps)){
  
  pooeded <- rra.tab[which(rownames(rra.tab) == grps[i])]
  pooeded <- pooeded[,which(apply(pooeded,2,sum) > 0)]
  rownames(pooeded) <- 1:nrow(pooeded) # for interspecific, we can't have duplicate rownames
  write.csv(pooeded, file = "file_for_ecosimr.csv") ###############
  my.params$Data.File <- "file_for_ecosimr.csv"
  
  Output.Results(my.params, Null.Model.Engine(my.params))
  
  x1 <- read.delim(file = "./Niche Overlap Output.txt", header = FALSE, sep = "\t")
  x1$V1 <- str_replace(x1$V1, " =  ", ":")
  x1$V1 <- str_replace(x1$V1, " > ", ":")
  x1$V1 <- str_replace(x1$V1, " < ", ":")
  x2 <- as.matrix(str_split(x1$V1, ":", simplify = TRUE))
  x2 <- as.matrix(x2[2:18,1:2])
  df1 <- as.data.frame(t(x2))
  #colnames(df1) <- colnames(df)
  #df <- rbind(df, df1[2,])
  df[i,3:17] <- as.character(unlist(df1[2,3:17]))
  df[i,1] <- grps[i]
  
}

beep("treasure")

df <- cbind(glom = "MOTU", df)
df$glom <- as.character(df$glom)

#write.csv(df, file = paste("./", vr, "_ecosimr_results.csv", sep = ""))

# - # Agglom

# MOTU level: All MOTUs

lv <- c("species", "genus", "family", "order")

for(m in 1:length(lv)){
  
  glom.diet <- tax_glom(final.diet, taxrank = lv[m])
  final.diet.ra = transform_sample_counts(glom.diet, function(x) x/sum(x))
  otudf <- as.matrix(otu_table(final.diet.ra))
  colnames(otudf) <- as.character(as.matrix(sample_data(final.diet.ra)[,vr]))
  grps <- unique(colnames(otudf))
  
  rra.tab <- t(otudf)
  rra.tab[1:5,1:5]
  
  for (i in 1:length(grps)){
    
    print(paste("Running ", grps[i], " at ", lv[m], " level", sep = ""))
    
    pooeded <- rra.tab[which(rownames(rra.tab) == grps[i])]
    pooeded <- pooeded[,which(apply(pooeded,2,sum) > 0)]
    rownames(pooeded) <- 1:nrow(pooeded) # for interspecific, we can't have duplicate rownames
    write.csv(pooeded, file = "file_for_ecosimr.csv") ###############
    my.params$Data.File <- "file_for_ecosimr.csv"
    
    Output.Results(my.params, Null.Model.Engine(my.params))
    
    x1 <- read.delim(file = "./Niche Overlap Output.txt", header = FALSE, sep = "\t")
    x1$V1 <- str_replace(x1$V1, " =  ", ":")
    x1$V1 <- str_replace(x1$V1, " > ", ":")
    x1$V1 <- str_replace(x1$V1, " < ", ":")
    x2 <- as.matrix(str_split(x1$V1, ":", simplify = TRUE))
    x2 <- as.matrix(x2[2:18,1:2])
    df1 <- as.data.frame(t(x2))
    #colnames(df1) <- colnames(df)
    #df <- rbind(df, df1[2,])
    n <- as.character(c(lv[m], grps[i], "Interspecific", as.character(unlist(df1[2,3:17]))))
    df <- rbind(df, n)
    #df.g[i,] <- as.character(unlist(df1[2,]))
    #df.g[i,1] <- chck1
    #df.g[i,2] <- chck2
    
  }
  
}

beep("treasure")

write.csv(df, file = paste("./", vr, "_ecosimr_results_RRA_Interspecific.csv", sep = ""))


#### Automate ecosimr: POO Metric ####

source("./ecosimr/EcoSimR - Main Source.R")

my.params <- Param.List

## Run Analysis
#my.params$Data.File <- paste("./ecosimr_output/", "all.motus.france.csv", sep = "")
my.params$Plot.Output <- "file"
my.params$Print.Output <- "file"
my.params$N.Reps <- 10000

# MOTU level: All MOTUs

otudf <- as.data.frame(as.matrix(otu_table(final.diet)))
otudf <- as.matrix(otu_table(final.diet))
#colnames(otudf) <- sample_data(final.diet)$shrew_zone
colnames(otudf) <- as.character(as.matrix(sample_data(final.diet)[,vr]))
#colnames(otudf) <- as.character(unlist(sample_data(final.diet)[,vr])) # same as previous line
grps <- unique(colnames(otudf))

BBA_Dprop=prop.table(otudf,2)
dim(BBA_Dprop)
apply(BBA_Dprop,2,sum)

# Using a threshold of 1% to declare presence/absence 
PA_BBA_D= BBA_Dprop
PA_BBA_D[PA_BBA_D>.0099999999]=1
PA_BBA_D[PA_BBA_D<.01]=0
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

library(gtools)

#permutations(n=length(grps),r=2,v=grps,repeats.allowed=FALSE)
chck <- combinations(n=length(grps),r=2,v=grps,repeats.allowed=FALSE)

chck[1,1]
chck[1,2]

inpt <- nrow(as.data.frame(chck[,]))

df <- data.frame(group_1 = 1:inpt,
                 group_2 = NA,
                 Seed = NA,
                 Num_of_reps = NA,
                 Elapsed_time = NA,
                 Metric = NA,
                 Algorithm = NA,
                 Observed = NA,
                 Mean_sim_index = NA,
                 Var_sim_index = NA,
                 low_95_1tail = NA,
                 upper_95_1tail = NA,
                 low_95_2tail = NA,
                 upper_95_2tail = NA,
                 P_obs_less_Null = NA,
                 P_obs_great_Null = NA,
                 P_obs_equal_Null = NA)

for (i in 1:inpt){
  
  chck1 <- chck[i,1]
  chck2 <- chck[i,2]
  pooeded <- rbind(grp.list[chck1,],
                   grp.list[chck2,])
  pooeded <- pooeded[,which(apply(pooeded,2,sum) > 0)]
  write.csv(pooeded, file = "file_for_ecosimr.csv")
  my.params$Data.File <- "file_for_ecosimr.csv"
  
  Output.Results(my.params, Null.Model.Engine(my.params))
  
  x1 <- read.delim(file = "./Niche Overlap Output.txt", header = FALSE, sep = "\t")
  x1$V1 <- str_replace(x1$V1, " =  ", ":")
  x1$V1 <- str_replace(x1$V1, " > ", ":")
  x1$V1 <- str_replace(x1$V1, " < ", ":")
  x2 <- as.matrix(str_split(x1$V1, ":", simplify = TRUE))
  x2 <- as.matrix(x2[2:18,1:2])
  df1 <- as.data.frame(t(x2))
  #colnames(df1) <- colnames(df)
  #df <- rbind(df, df1[2,])
  df[i,] <- as.character(unlist(df1[2,]))
  df[i,1] <- chck1
  df[i,2] <- chck2
  
}

beep("treasure")

df <- cbind(glom = "MOTU", df)
df$glom <- as.character(df$glom)

#write.csv(df, file = paste("./", vr, "_ecosimr_results.csv", sep = ""))


# - # Agglom

# MOTU level: All MOTUs

lv <- c("species", "genus", "family", "order")

for(m in 1:length(lv)){
  
  glom.diet <- tax_glom(final.diet, taxrank = lv[m])
  
  otudf <- as.data.frame(as.matrix(otu_table(glom.diet)))
  otudf <- as.matrix(otu_table(glom.diet))
  #colnames(otudf) <- sample_data(final.diet)$shrew_zone
  colnames(otudf) <- as.character(as.matrix(sample_data(glom.diet)[,vr]))
  #colnames(otudf) <- as.character(unlist(sample_data(final.diet)[,vr])) # same as previous line
  grps <- unique(colnames(otudf))
  
  BBA_Dprop=prop.table(otudf,2)
  dim(BBA_Dprop)
  apply(BBA_Dprop,2,sum)
  
  # Using a threshold of 1% to declare presence/absence 
  PA_BBA_D= BBA_Dprop
  PA_BBA_D[PA_BBA_D>.0099999999]=1
  PA_BBA_D[PA_BBA_D<.01]=0
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
  
  
  
  #permutations(n=length(grps),r=2,v=grps,repeats.allowed=FALSE)
  chck <- combinations(n=length(grps),r=2,v=grps,repeats.allowed=FALSE)
  
  #chck[1,1]
  #chck[1,2]
  
  inpt <- nrow(as.data.frame(chck[,]))
  
  for (i in 1:inpt){
    
    chck1 <- chck[i,1]
    chck2 <- chck[i,2]
    pooeded <- rbind(grp.list[chck1,],
                     grp.list[chck2,])
    pooeded <- pooeded[,which(apply(pooeded,2,sum) > 0)]
    write.csv(pooeded, file = "file_for_ecosimr.csv")
    my.params$Data.File <- "file_for_ecosimr.csv"
    
    Output.Results(my.params, Null.Model.Engine(my.params))
    
    x1 <- read.delim(file = "./Niche Overlap Output.txt", header = FALSE, sep = "\t")
    x1$V1 <- str_replace(x1$V1, " =  ", ":")
    x1$V1 <- str_replace(x1$V1, " > ", ":")
    x1$V1 <- str_replace(x1$V1, " < ", ":")
    x2 <- as.matrix(str_split(x1$V1, ":", simplify = TRUE))
    x2 <- as.matrix(x2[2:18,1:2])
    df1 <- as.data.frame(t(x2))
    #colnames(df1) <- colnames(df)
    #df <- rbind(df, df1[2,])
    n <- as.character(c(lv[m], chck1, chck2, as.character(unlist(df1[2,3:17]))))
    df <- rbind(df, n)
    #df.g[i,] <- as.character(unlist(df1[2,]))
    #df.g[i,1] <- chck1
    #df.g[i,2] <- chck2
    
  }
  
  
}

beep("treasure")

write.csv(df, file = paste("./", vr, "_ecosimr_results.csv", sep = ""))


#### Automate ecosimr: POO Metric: Core Diet ####

library(microbiome)

# Core Filter no. 1
filt.1 <- subset_samples(final.diet, shrew_country == "Pygmy_Ireland")
core.members.1 <- core_members(filt.1, detection = 0, prevalence = 2/nsamples(filt.1))
filt.2 <- subset_samples(final.diet, shrew_country == "Pygmy_France")
core.members.2 <- core_members(filt.2, detection = 0, prevalence = 2/nsamples(filt.2))
filt.3 <- subset_samples(final.diet, shrew_country == "GWTS_Ireland")
core.members.3 <- core_members(filt.3, detection = 0, prevalence = 2/nsamples(filt.3))
filt.4 <- subset_samples(final.diet, shrew_country == "GWTS_France")
core.members.4 <- core_members(filt.4, detection = 0, prevalence = 2/nsamples(filt.4))

core.members <- unique(c(core.members.1,core.members.2,core.members.3,core.members.4))

final.diet.cr <- prune_taxa(core.members, final.diet)

# Core Filter no. 2
core.members <- core_members(final.diet, detection = 0, prevalence = 2/nsamples(final.diet))

final.diet.cr <- prune_taxa(core.members, final.diet)



source("./ecosimr/EcoSimR - Main Source.R")

my.params <- Param.List

## Run Analysis
#my.params$Data.File <- paste("./ecosimr_output/", "all.motus.france.csv", sep = "")
my.params$Plot.Output <- "file"
my.params$Print.Output <- "file"
my.params$N.Reps <- 10000

# MOTU level: All MOTUs

vr <- "shrew_country"

otudf <- as.data.frame(as.matrix(otu_table(final.diet.cr)))
otudf <- as.matrix(otu_table(final.diet.cr))
#colnames(otudf) <- sample_data(final.diet)$shrew_zone
colnames(otudf) <- as.character(as.matrix(sample_data(final.diet.cr)[,vr]))
#colnames(otudf) <- as.character(unlist(sample_data(final.diet)[,vr])) # same as previous line
grps <- unique(colnames(otudf))

BBA_Dprop=prop.table(otudf,2)
dim(BBA_Dprop)
apply(BBA_Dprop,2,sum)

# Using a threshold of 1% to declare presence/absence 
PA_BBA_D= BBA_Dprop
PA_BBA_D[PA_BBA_D>.0099999999]=1
PA_BBA_D[PA_BBA_D<.01]=0
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

library(gtools)

#permutations(n=length(grps),r=2,v=grps,repeats.allowed=FALSE)
chck <- combinations(n=length(grps),r=2,v=grps,repeats.allowed=FALSE)

chck[1,1]
chck[1,2]

inpt <- nrow(as.data.frame(chck[,]))

df <- data.frame(group_1 = 1:inpt,
                 group_2 = NA,
                 Seed = NA,
                 Num_of_reps = NA,
                 Elapsed_time = NA,
                 Metric = NA,
                 Algorithm = NA,
                 Observed = NA,
                 Mean_sim_index = NA,
                 Var_sim_index = NA,
                 low_95_1tail = NA,
                 upper_95_1tail = NA,
                 low_95_2tail = NA,
                 upper_95_2tail = NA,
                 P_obs_less_Null = NA,
                 P_obs_great_Null = NA,
                 P_obs_equal_Null = NA)

for (i in 1:inpt){
  
  chck1 <- chck[i,1]
  chck2 <- chck[i,2]
  pooeded <- rbind(grp.list[chck1,],
                   grp.list[chck2,])
  pooeded <- pooeded[,which(apply(pooeded,2,sum) > 0)]
  write.csv(pooeded, file = "file_for_ecosimr.csv")
  my.params$Data.File <- "file_for_ecosimr.csv"
  
  Output.Results(my.params, Null.Model.Engine(my.params))
  
  x1 <- read.delim(file = "./Niche Overlap Output.txt", header = FALSE, sep = "\t")
  x1$V1 <- str_replace(x1$V1, " =  ", ":")
  x1$V1 <- str_replace(x1$V1, " > ", ":")
  x1$V1 <- str_replace(x1$V1, " < ", ":")
  x2 <- as.matrix(str_split(x1$V1, ":", simplify = TRUE))
  x2 <- as.matrix(x2[2:18,1:2])
  df1 <- as.data.frame(t(x2))
  #colnames(df1) <- colnames(df)
  #df <- rbind(df, df1[2,])
  df[i,] <- as.character(unlist(df1[2,]))
  df[i,1] <- chck1
  df[i,2] <- chck2
  
}

beep("treasure")

df <- cbind(glom = "MOTU", df)
df$glom <- as.character(df$glom)

#write.csv(df, file = paste("./", vr, "_ecosimr_results.csv", sep = ""))


# - # Agglom

# MOTU level: All MOTUs

lv <- c("species", "genus", "family", "order")

for(m in 1:length(lv)){
  
  glom.diet <- tax_glom(final.diet.cr, taxrank = lv[m])
  
  otudf <- as.data.frame(as.matrix(otu_table(glom.diet)))
  otudf <- as.matrix(otu_table(glom.diet))
  #colnames(otudf) <- sample_data(final.diet)$shrew_zone
  colnames(otudf) <- as.character(as.matrix(sample_data(glom.diet)[,vr]))
  #colnames(otudf) <- as.character(unlist(sample_data(final.diet)[,vr])) # same as previous line
  grps <- unique(colnames(otudf))
  
  BBA_Dprop=prop.table(otudf,2)
  dim(BBA_Dprop)
  apply(BBA_Dprop,2,sum)
  
  # Using a threshold of 1% to declare presence/absence 
  PA_BBA_D= BBA_Dprop
  PA_BBA_D[PA_BBA_D>.0099999999]=1
  PA_BBA_D[PA_BBA_D<.01]=0
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
  
  
  
  #permutations(n=length(grps),r=2,v=grps,repeats.allowed=FALSE)
  chck <- combinations(n=length(grps),r=2,v=grps,repeats.allowed=FALSE)
  
  #chck[1,1]
  #chck[1,2]
  
  inpt <- nrow(as.data.frame(chck[,]))
  
  for (i in 1:inpt){
    
    chck1 <- chck[i,1]
    chck2 <- chck[i,2]
    pooeded <- rbind(grp.list[chck1,],
                     grp.list[chck2,])
    pooeded <- pooeded[,which(apply(pooeded,2,sum) > 0)]
    write.csv(pooeded, file = "file_for_ecosimr.csv")
    my.params$Data.File <- "file_for_ecosimr.csv"
    
    Output.Results(my.params, Null.Model.Engine(my.params))
    
    x1 <- read.delim(file = "./Niche Overlap Output.txt", header = FALSE, sep = "\t")
    x1$V1 <- str_replace(x1$V1, " =  ", ":")
    x1$V1 <- str_replace(x1$V1, " > ", ":")
    x1$V1 <- str_replace(x1$V1, " < ", ":")
    x2 <- as.matrix(str_split(x1$V1, ":", simplify = TRUE))
    x2 <- as.matrix(x2[2:18,1:2])
    df1 <- as.data.frame(t(x2))
    #colnames(df1) <- colnames(df)
    #df <- rbind(df, df1[2,])
    n <- as.character(c(lv[m], chck1, chck2, as.character(unlist(df1[2,3:17]))))
    df <- rbind(df, n)
    #df.g[i,] <- as.character(unlist(df1[2,]))
    #df.g[i,1] <- chck1
    #df.g[i,2] <- chck2
    
  }
  
  
}

beep("treasure")

write.csv(df, file = paste("./", vr, "_ecosimr_results_core_diet.csv", sep = ""))

