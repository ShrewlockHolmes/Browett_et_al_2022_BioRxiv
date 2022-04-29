#### Shrew Weight ####

df <- read.csv(file = "Body_Measurements.csv", sep=",")

df.wgt <- df[which(df[,"Bodypart"] == "Weight"),]

#df.pyg$loc_bod <- paste(df.pyg$Location, df.pyg$Bodypart, sep = "_")
#df.wgt$Zone <- stringr::str_replace(df.wgt$Zone, "France", "Belle Ile")
df.wgt$shr_zne <- paste(df.wgt$Species, df.wgt$Zone, sep = "_")
#df.wgt$shr_zne <- factor(df.wgt$shr_zne, levels = c("GWTS_France", "GWTS_Inside", "GWTS_Edge", "Pygmy_France", "Pygmy_Outside", "Pygmy_Edge"))
df.wgt$shr_zne <- factor(df.wgt$shr_zne, levels = c("GWTS_France", "GWTS_Edge", "GWTS_Inside", "Pygmy_Outside", "Pygmy_Edge", "Pygmy_France"))
#df.wgt$shr_zne <- factor(df.wgt$shr_zne, levels = c("GWTS_France", "GWTS_Edge", "GWTS_Inside", "Pygmy_Edge", "Pygmy_Outside", "Pygmy_France"))
df.wgt$Zone <- factor(df.wgt$Zone, levels = c("France", "Inside", "Edge", "Outside"))
df.wgt$shr_cntry <- paste(df.wgt$Species, df.wgt$Location, sep = "_")
df.wgt$shr_cntry <- factor(df.wgt$shr_cntry, levels = c("GWTS_France", "GWTS_Ireland", "Pygmy_France", "Pygmy_Ireland"))


G.1 <- ggplot(df.wgt, aes(Species, Measure, group = shr_zne)) +
  geom_boxplot(aes(fill= Zone), alpha = 0.75) +
  geom_jitter(aes(colour = Zone), position = position_jitterdodge(), size=3) +
  stat_summary(fun = mean, geom="point", position=position_dodge(.75), shape=20, size = 7) +
  #scale_colour_brewer(palette = "Dark2") + # Manually sets colours according to object "palette"
  #scale_colour_manual(values = (rep("black",4))) +
  scale_colour_manual(values= c("France" = "#1B9E77",
                                "Inside" = "#D95F02",
                                "Edge" = "#7570B3",
                                "Outside" = "#E7298A"),
                      labels = c("France" = "Belle Ile",
                                 "Inside" = "Inside",
                                 "Edge" = "Edge",
                                 "Outside" = "Outside"),
                      name = element_blank()) +
  #scale_fill_brewer(palette = "Dark2") +
  scale_fill_manual(values= c("France" = "#1B9E77",
                              "Inside" = "#D95F02",
                              "Edge" = "#7570B3",
                              "Outside" = "#E7298A"),
                    labels = c("France" = "Belle Ile",
                               "Inside" = "Inside",
                               "Edge" = "Edge",
                               "Outside" = "Outside"),
                    name = element_blank()) +
  scale_x_discrete(labels= c("C. russula", "S. minutus")) + # Manually labels groups using set described in "key" object
  # scale_y_continuous(limits = c(0, 0.55), expand = c(0, 0)) + # Sets Graph Value Limits
  theme_classic() + # Chooses Theme
  theme(axis.text = element_text(size = 10, face= "bold"), # Size of Axis measurements/Values Text
        axis.title=element_text(size=14), # Sets size of Axis Titles. face= "bold" will change the font to bold
        panel.grid.major = element_line(colour = "grey95"), # Sets colour of major grid lines
        panel.grid.minor = element_line(colour = "grey95") # Sets colour of minor grid lines
  ) +
  labs(y= "Weight (grams)", x = "") + # Labels Axis'
  #theme(strip.text.x = element_text(size = 12, face = "bold.italic")) +
  #theme(legend.position = "right", # c(0.8,0.8),
  #      legend.title = element_blank(),
  #     legend.text = element_text(size = 15),
  #    legend.key.size = unit(1, "cm")) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20,
                                    face = "bold")) +
  theme(axis.text.x = element_text(size = 12 ,
                                   face = "bold",
                                   angle = 45,
                                   color = "black",
                                   hjust = 1),
        axis.text.y = element_text(size = 12,
                                   face = "bold"))

G.1

legd <- cowplot::get_legend(G.1 + theme(legend.position = "right",
                                        legend.title = element_blank(),
                                        legend.text = element_text(size = 15),
                                        legend.key.size = unit(1, "cm")))
grid.newpage()
grid.draw(legd)

aggregate(df.wgt$Measure, list(df.wgt$shr_zne), mean)
aggregate(df.wgt$Measure, list(df.wgt$shr_zne), sd)

options(scipen = 999)
# ANOVA using the species richness
a1 <- aov(df.wgt$Measure ~ df.wgt$shr_zne)
TukeyHSD(x=a1, 'df.wgt$shr_zne', conf.level=0.95)
plot(TukeyHSD(x=a1, 'df.wgt$shr_zne', conf.level=0.95))

aggregate(df.wgt$Measure, list(df.wgt$shr_cntry), mean)
aggregate(df.wgt$Measure, list(df.wgt$shr_cntry), sd)

options(scipen = 999)
# ANOVA using the species richness
a1 <- aov(df.wgt$Measure ~ df.wgt$shr_cntry)
TukeyHSD(x=a1, 'df.wgt$shr_cntry', conf.level=0.95)
plot(TukeyHSD(x=a1, 'df.wgt$shr_cntry', conf.level=0.95))


#### Alternative Size Figure ####

G.1 <- ggplot(df.wgt, aes(shr_zne, Measure, group = shr_zne)) +
  geom_boxplot(aes(fill= Species), alpha = 0.75) +
  geom_jitter(aes(colour = Species), position = position_jitterdodge(), size=2) +
  stat_summary(fun = mean, geom="point", position=position_dodge(.75), shape=20, size = 7) +
  #scale_colour_brewer(palette = "Dark2") + # Manually sets colours according to object "palette"
  #scale_colour_manual(values = (rep("black",4))) +
  facet_wrap(~ Species, scale = "free_x") +
  scale_colour_manual(values= c("GWTS" = "#1B9E77",
                                "Pygmy" = "#D95F02"),
                      #"Edge" = "#7570B3",
                      #"Outside" = "#E7298A"),
                      labels = c("GWTS" = "C. russula",
                                 "Pygmy" = "S. minutus"),
                      #"Edge" = "Edge",
                      #"Outside" = "Outside"),
                      name = element_blank()) +
  #scale_fill_brewer(palette = "Dark2") +
  scale_fill_manual(values= c("GWTS" = "#1B9E77",
                              "Pygmy" = "#D95F02"),
                    #"Edge" = "#7570B3",
                    #"Outside" = "#E7298A"),
                    labels = c("GWTS" = "C. russula",
                               "Pygmy" = "S. minutus"),
                    #"Edge" = "Edge",
                    #"Outside" = "Outside"),
                    name = element_blank()) +
  scale_x_discrete(labels=c("GWTS_France" = "C.russula\nBelle Ile",
                            "Pygmy_France" = "S.minutus\nBelle Ile",
                            "GWTS_Inside" = "C.russula\nInside",
                            "GWTS_Edge" = "C.russula\nEdge",
                            "Pygmy_Edge" = "S.minutus\nEdge",
                            "Pygmy_Outside" = "S.minutus\nOutside")) +
  # scale_y_continuous(limits = c(0, 0.55), expand = c(0, 0)) + # Sets Graph Value Limits
  theme_classic() + # Chooses Theme
  theme(axis.text = element_text(size = 10, face= "bold"), # Size of Axis measurements/Values Text
        axis.title=element_text(size=14), # Sets size of Axis Titles. face= "bold" will change the font to bold
        panel.grid.major = element_line(colour = "grey95"), # Sets colour of major grid lines
        panel.grid.minor = element_line(colour = "grey95") # Sets colour of minor grid lines
  ) +
  labs(y= "Weight (grams)", x = "") + # Labels Axis'
  #theme(strip.text.x = element_text(size = 12, face = "bold.italic")) +
  #theme(legend.position = "right", # c(0.8,0.8),
  #      legend.title = element_blank(),
  #     legend.text = element_text(size = 15),
  #    legend.key.size = unit(1, "cm")) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20,
                                    face = "plain")) +
  theme(axis.text.x = element_text(size = 12 ,
                                   face = "bold",
                                   #angle = 45,
                                   color = "black"),
        #hjust = 1),
        axis.text.y = element_text(size = 12,
                                   face = "bold"))

G.1

legd <- cowplot::get_legend(G.1 + theme(legend.position = "right",
                                        legend.title = element_blank(),
                                        legend.text = element_text(size = 15),
                                        legend.key.size = unit(1, "cm")))
grid.newpage()
grid.draw(legd)



#### Body Length ####

df.wgt <- df[which(df[,"Bodypart"] == "Total"),]

#df.pyg$loc_bod <- paste(df.pyg$Location, df.pyg$Bodypart, sep = "_")
#df.wgt$Zone <- stringr::str_replace(df.wgt$Zone, "France", "Belle Ile")
df.wgt$shr_zne <- paste(df.wgt$Species, df.wgt$Zone, sep = "_")
#df.wgt$shr_zne <- factor(df.wgt$shr_zne, levels = c("GWTS_France", "GWTS_Edge", "GWTS_Inside", "Pygmy_France", "Pygmy_Edge", "Pygmy_Outside"))
df.wgt$shr_zne <- factor(df.wgt$shr_zne, levels = c("GWTS_France", "GWTS_Edge", "GWTS_Inside", "Pygmy_Outside", "Pygmy_Edge", "Pygmy_France"))
df.wgt$Zone <- factor(df.wgt$Zone, levels = c("France", "Inside", "Edge", "Outside"))
df.wgt$shr_cntry <- paste(df.wgt$Species, df.wgt$Location, sep = "_")
df.wgt$shr_cntry <- factor(df.wgt$shr_cntry, levels = c("GWTS_France", "GWTS_Ireland", "Pygmy_France", "Pygmy_Ireland"))


G.2 <- ggplot(df.wgt, aes(Species, Measure, group = shr_zne)) +
  geom_boxplot(aes(fill= Zone), alpha = 0.75) +
  geom_jitter(aes(colour = Zone), position = position_jitterdodge(), size=3) +
  stat_summary(fun = mean, geom="point", position=position_dodge(.75), shape=20, size = 7) +
  #scale_colour_brewer(palette = "Dark2") + # Manually sets colours according to object "palette"
  #scale_colour_manual(values = (rep("black",4))) +
  scale_colour_manual(values= c("France" = "#1B9E77",
                                "Inside" = "#D95F02",
                                "Edge" = "#7570B3",
                                "Outside" = "#E7298A"),
                      labels = c("France" = "Belle Ile",
                                 "Inside" = "Inside",
                                 "Edge" = "Edge",
                                 "Outside" = "Outside"),
                      name = element_blank()) +
  #scale_fill_brewer(palette = "Dark2") +
  scale_fill_manual(values= c("France" = "#1B9E77",
                              "Inside" = "#D95F02",
                              "Edge" = "#7570B3",
                              "Outside" = "#E7298A"),
                    labels = c("France" = "Belle Ile",
                               "Inside" = "Inside",
                               "Edge" = "Edge",
                               "Outside" = "Outside"),
                    name = element_blank()) +
  scale_x_discrete(labels= c("C. russula", "S. minutus")) + # Manually labels groups using set described in "key" object
  # scale_y_continuous(limits = c(0, 0.55), expand = c(0, 0)) + # Sets Graph Value Limits
  theme_classic() + # Chooses Theme
  theme(axis.text = element_text(size = 10, face= "bold"), # Size of Axis measurements/Values Text
        axis.title=element_text(size=14), # Sets size of Axis Titles. face= "bold" will change the font to bold
        panel.grid.major = element_line(colour = "grey95"), # Sets colour of major grid lines
        panel.grid.minor = element_line(colour = "grey95") # Sets colour of minor grid lines
  ) +
  labs(y= "Weight (grams)", x = "") + # Labels Axis'
  #theme(strip.text.x = element_text(size = 12, face = "bold.italic")) +
  #theme(legend.position = "right", # c(0.8,0.8),
  #      legend.title = element_blank(),
  #     legend.text = element_text(size = 15),
  #    legend.key.size = unit(1, "cm")) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20,
                                    face = "bold")) +
  theme(axis.text.x = element_text(size = 12 ,
                                   face = "bold",
                                   angle = 45,
                                   color = "black",
                                   hjust = 1),
        axis.text.y = element_text(size = 12,
                                   face = "bold"))

G.2

legd <- cowplot::get_legend(G.1 + theme(legend.position = "right",
                                        legend.title = element_blank(),
                                        legend.text = element_text(size = 15),
                                        legend.key.size = unit(1, "cm")))
grid.newpage()
grid.draw(legd)

# Facet Labels
supp.labels <- c("C. russula", "S. minutus")
names(supp.labels) <- c("GWTS", "Pygmy")

G.3 <- ggplot(df.wgt, aes(shr_zne, Measure, group = shr_zne)) +
  geom_boxplot(aes(fill= Species), alpha = 0.75) +
  geom_jitter(aes(colour = Species), position = position_jitterdodge(), size=2) +
  stat_summary(fun = mean, geom="point", position=position_dodge(.75), shape=20, size = 7) +
  #scale_colour_brewer(palette = "Dark2") + # Manually sets colours according to object "palette"
  #scale_colour_manual(values = (rep("black",4))) +
  facet_wrap(~ Species, 
             labeller = labeller(Species = supp.labels),
             scale = "free_x") +
  scale_colour_manual(values= c("GWTS" = "#1B9E77",
                                "Pygmy" = "#D95F02"),
                      #"Edge" = "#7570B3",
                      #"Outside" = "#E7298A"),
                      labels = c("GWTS" = "C. russula",
                                 "Pygmy" = "S. minutus"),
                      #"Edge" = "Edge",
                      #"Outside" = "Outside"),
                      name = element_blank()) +
  #scale_fill_brewer(palette = "Dark2") +
  scale_fill_manual(values= c("GWTS" = "#1B9E77",
                              "Pygmy" = "#D95F02"),
                    #"Edge" = "#7570B3",
                    #"Outside" = "#E7298A"),
                    labels = c("GWTS" = "C. russula",
                               "Pygmy" = "S. minutus"),
                    #"Edge" = "Edge",
                    #"Outside" = "Outside"),
                    name = element_blank()) +
  scale_x_discrete(labels=c("GWTS_France" = "Belle Ile",
                            "Pygmy_France" = "Belle Ile",
                            "GWTS_Edge" = "'During'",
                            "GWTS_Inside" = "'After'",
                            "Pygmy_Edge" = "'During'",
                            "Pygmy_Outside" = "'Before'")) +
  # scale_y_continuous(limits = c(0, 0.55), expand = c(0, 0)) + # Sets Graph Value Limits
  theme_classic() + # Chooses Theme
  theme(axis.text = element_text(size = 10, face= "bold"), # Size of Axis measurements/Values Text
        axis.title=element_text(size=10), # Sets size of Axis Titles. face= "bold" will change the font to bold
        panel.grid.major = element_line(colour = "grey95"), # Sets colour of major grid lines
        panel.grid.minor = element_line(colour = "grey95") # Sets colour of minor grid lines
  ) +
  labs(y= "Total Body Length (mm)", x = "") + # Labels Axis'
  theme(strip.text.x = element_text(size = 12, face = "bold.italic")) +
  #theme(legend.position = "right", # c(0.8,0.8),
  #      legend.title = element_blank(),
  #     legend.text = element_text(size = 15),
  #    legend.key.size = unit(1, "cm")) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14,
                                    face = "plain")) +
  theme(axis.text.x = element_text(size = 12 ,
                                   face = "bold",
                                   #angle = 45,
                                   color = "black"),
        #hjust = 1),
        axis.text.y = element_text(size = 12,
                                   face = "bold"))

G.3


aggregate(df.wgt$Measure, list(df.wgt$shr_zne), mean)
aggregate(df.wgt$Measure, list(df.wgt$shr_zne), sd)

sts <- data.frame(table(df.wgt$shr_zne))
sts$mean <- as.matrix(aggregate(df.wgt$Measure, list(df.wgt$shr_zne), mean))[,2]
sts$mean <- as.numeric(sts$mean)
sts$stnd_dev <- as.matrix(aggregate(df.wgt$Measure, list(df.wgt$shr_zne), sd))[,2]
sts$stnd_dev <- as.numeric(sts$stnd_dev)
str(sts)
sts$SE <- (sts$stnd_dev / sqrt(sts$Freq))
sts$SE

options(scipen = 999)
# ANOVA using the species richness
a1 <- aov(df.wgt$Measure ~ df.wgt$shr_zne)
TukeyHSD(x=a1, 'df.wgt$shr_zne', conf.level=0.95)
plot(TukeyHSD(x=a1, 'df.wgt$shr_zne', conf.level=0.95))

aggregate(df.wgt$Measure, list(df.wgt$shr_cntry), mean)
aggregate(df.wgt$Measure, list(df.wgt$shr_cntry), sd)

sts <- data.frame(table(df.wgt$shr_cntry))
sts$mean <- as.matrix(aggregate(df.wgt$Measure, list(df.wgt$shr_cntry), mean))[,2]
sts$mean <- as.numeric(sts$mean)
sts$stnd_dev <- as.matrix(aggregate(df.wgt$Measure, list(df.wgt$shr_cntry), sd))[,2]
sts$stnd_dev <- as.numeric(sts$stnd_dev)
str(sts)
sts$SE <- (sts$stnd_dev / sqrt(sts$Freq))
sts$SE


options(scipen = 999)
# ANOVA using the species richness
a1 <- aov(df.wgt$Measure ~ df.wgt$shr_cntry)
TukeyHSD(x=a1, 'df.wgt$shr_cntry', conf.level=0.95)
plot(TukeyHSD(x=a1, 'df.wgt$shr_cntry', conf.level=0.95))

#### Weight by Season ####

library(reshape)

df.wd <- read.csv(file = "./Full_Data_Wide_Format_Season.csv", sep = ',')

df.sn <- melt(data = df.wd,
              id = c("Species",
                     "Zone",
                     "Season"))

df.wgt.s <- df.sn[which(df.sn[,"variable"] == "Weight"),]

df.wgt.s$zone_season <- paste(df.wgt.s$Zone,
                              df.wgt.s$Season,
                              sep = '_')

df.wgt.s$species_zone <- paste(df.wgt.s$Species,
                              df.wgt.s$Zone,
                              sep = '_')
df.wgt.s$species_zone <- factor(df.wgt.s$species_zone, 
                                       levels = c("C.russula_Belle Ile",
                                                  "C.russula_During",
                                                  "C.russula_After",
                                                  "S.minutus_Before",
                                                  "S.minutus_During",
                                                  "S.minutus_Belle Ile"))


df.wgt.s$species_zone_season <- paste(df.wgt.s$Species,
                                      df.wgt.s$Zone,
                                      df.wgt.s$Season,
                                      sep = '_')

#df.wgt.s$shr_zne <- paste(df.wgt.s$Species, df.wgt$Zone, sep = "_")

df.wgt.s$Zone <- factor(df.wgt.s$Zone, levels = c("Belle Ile", "During", "After", "Before"))
df.wgt.s$Season <- factor(df.wgt.s$Season, levels = c("Summer", "Spring"))
df.wgt.s$species_zone_season <- factor(df.wgt.s$species_zone_season, 
                                       levels = c("C.russula_Belle Ile_Summer",
                                                  "C.russula_Belle Ile_Spring",
                                                  "C.russula_During_Summer",
                                                  "C.russula_During_Spring",
                                                  "C.russula_After_Summer",
                                                  "C.russula_After_Spring",
                                                  "S.minutus_Before_Summer",
                                                  "S.minutus_Before_Spring",
                                                  "S.minutus_During_Summer",
                                                  "S.minutus_During_Spring",
                                                  "S.minutus_Belle Ile_Summer",
                                                  "S.minutus_Belle Ile_Spring"))

#facet species
#group zone
# colour season

G.1 <- ggplot(df.wgt.s, aes(species_zone, value, group = species_zone_season)) +
  geom_boxplot(aes(fill= Season), alpha = 0.75) +
  geom_jitter(aes(colour = Season), position = position_jitterdodge(), size=3) +
  stat_summary(fun = mean, geom="point", position=position_dodge(.75), shape=20, size = 7) +
  facet_wrap(~ Species, scales = 'free_x') +
  #scale_colour_brewer(palette = "Dark2") +
  #scale_fill_brewer(palette = "Dark2") +
  scale_colour_manual(values = c("#7570B3", "#E7298A")) +
  scale_fill_manual(values = c("#7570B3", "#E7298A")) +
  theme_classic() + # Chooses Theme
  theme(axis.text = element_text(size = 10, face= "bold"), # Size of Axis measurements/Values Text
        axis.title=element_text(size=14), # Sets size of Axis Titles. face= "bold" will change the font to bold
        panel.grid.major = element_line(colour = "grey95"), # Sets colour of major grid lines
        panel.grid.minor = element_line(colour = "grey95") # Sets colour of minor grid lines
  ) +
  labs(y= "Weight (grams)", x = "") + # Labels Axis'
  theme(strip.text.x = element_text(size = 12, face = "bold.italic")) +
  #theme(legend.position = "right", # c(0.8,0.8),
  #      legend.title = element_blank(),
  #     legend.text = element_text(size = 15),
  #    legend.key.size = unit(1, "cm")) +
  theme(legend.position = c(0.90,0.85)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20,
                                    face = "bold")) +
  theme(axis.text.x = element_text(size = 12 ,
                                   face = "bold",
                                   angle = 45,
                                   color = "black",
                                   hjust = 1),
        axis.text.y = element_text(size = 12,
                                   face = "bold"))


G.1

aggregate(df.wgt.s$value, list(df.wgt.s$species_zone_season), mean)
aggregate(df.wgt.s$value, list(df.wgt.s$species_zone_season), sd)

sts <- data.frame(table(df.wgt.s$species_zone_season))
sts$mean <- as.matrix(aggregate(df.wgt.s$value, list(df.wgt.s$species_zone_season), mean))[,2]
sts$mean <- as.numeric(sts$mean)
sts$stnd_dev <- as.matrix(aggregate(df.wgt.s$value, list(df.wgt.s$species_zone_season), sd))[,2]
sts$stnd_dev <- as.numeric(sts$stnd_dev)
str(sts)
sts$SE <- (sts$stnd_dev / sqrt(sts$Freq))
sts$SE

#### Length by Season ####

df.wgt.s <- df.sn[which(df.sn[,"variable"] == "Total_Length"),]

df.wgt.s$zone_season <- paste(df.wgt.s$Zone,
                              df.wgt.s$Season,
                              sep = '_')

df.wgt.s$species_zone <- paste(df.wgt.s$Species,
                               df.wgt.s$Zone,
                               sep = '_')
df.wgt.s$species_zone <- factor(df.wgt.s$species_zone, 
                                levels = c("C.russula_Belle Ile",
                                           "C.russula_During",
                                           "C.russula_After",
                                           "S.minutus_Before",
                                           "S.minutus_During",
                                           "S.minutus_Belle Ile"))

df.wgt.s$species_zone_season <- paste(df.wgt.s$Species,
                                      df.wgt.s$Zone,
                                      df.wgt.s$Season,
                                      sep = '_')

#df.wgt.s$shr_zne <- paste(df.wgt.s$Species, df.wgt$Zone, sep = "_")

df.wgt.s$Zone <- factor(df.wgt.s$Zone, levels = c("Belle Ile", "During", "After", "Before"))
df.wgt.s$Season <- factor(df.wgt.s$Season, levels = c("Summer", "Spring"))
df.wgt.s$species_zone_season <- factor(df.wgt.s$species_zone_season, 
                                       levels = c("C.russula_Belle Ile_Summer",
                                                  "C.russula_Belle Ile_Spring",
                                                  "C.russula_During_Summer",
                                                  "C.russula_During_Spring",
                                                  "C.russula_After_Summer",
                                                  "C.russula_After_Spring",
                                                  "S.minutus_Before_Summer",
                                                  "S.minutus_Before_Spring",
                                                  "S.minutus_During_Summer",
                                                  "S.minutus_During_Spring",
                                                  "S.minutus_Belle Ile_Summer",
                                                  "S.minutus_Belle Ile_Spring"))
#facet species
#group zone
# colour season

G.2 <- ggplot(df.wgt.s, aes(species_zone, value, group = species_zone_season)) +
  geom_boxplot(aes(fill= Season), alpha = 0.75) +
  geom_jitter(aes(colour = Season), position = position_jitterdodge(), size=3) +
  stat_summary(fun = mean, geom="point", position=position_dodge(.75), shape=20, size = 7) +
  facet_wrap(~ Species, scales = 'free_x') +
  #scale_colour_brewer(palette = "Dark2") +
  #scale_fill_brewer(palette = "Dark2") +
  scale_colour_manual(values = c("#7570B3", "#E7298A")) +
  scale_fill_manual(values = c("#7570B3", "#E7298A")) +
  theme_classic() + # Chooses Theme
  theme(axis.text = element_text(size = 10, face= "bold"), # Size of Axis measurements/Values Text
        axis.title=element_text(size=14), # Sets size of Axis Titles. face= "bold" will change the font to bold
        panel.grid.major = element_line(colour = "grey95"), # Sets colour of major grid lines
        panel.grid.minor = element_line(colour = "grey95") # Sets colour of minor grid lines
  ) +
  labs(y= "Length (mm)", x = "") + # Labels Axis'
  theme(strip.text.x = element_text(size = 12, face = "bold.italic")) +
  #theme(legend.position = "right", # c(0.8,0.8),
  #      legend.title = element_blank(),
  #     legend.text = element_text(size = 15),
  #    legend.key.size = unit(1, "cm")) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20,
                                    face = "bold")) +
  theme(axis.text.x = element_text(size = 12 ,
                                   face = "bold",
                                   angle = 45,
                                   color = "black",
                                   hjust = 1),
        axis.text.y = element_text(size = 12,
                                   face = "bold"))


G.2

aggregate(df.wgt.s$value, list(df.wgt.s$species_zone_season), mean)
aggregate(df.wgt.s$value, list(df.wgt.s$species_zone_season), sd)

sts <- data.frame(table(df.wgt.s$species_zone_season))
sts$mean <- as.matrix(aggregate(df.wgt.s$value, list(df.wgt.s$species_zone_season), mean))[,2]
sts$mean <- as.numeric(sts$mean)
sts$stnd_dev <- as.matrix(aggregate(df.wgt.s$value, list(df.wgt.s$species_zone_season), sd))[,2]
sts$stnd_dev <- as.numeric(sts$stnd_dev)
str(sts)
sts$SE <- (sts$stnd_dev / sqrt(sts$Freq))
sts$SE


grid.arrange(G.1,
             G.2,
             nrow = 2)

grid.arrange(G.1 + theme(axis.text.x = element_blank()),
             G.2 +
               theme(strip.background.x = element_blank(),
                     strip.text.x = element_blank()) +
               theme(axis.text.x = element_blank()),
             nrow = 2)
