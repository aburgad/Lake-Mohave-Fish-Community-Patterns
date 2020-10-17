#####################################
#####################################
# Lake Mohave Fish Community Dynamics
#####################################
#####################################

# Install required packages if not already
package_list <- c("tidyverse","tidyr","vegan",
                  "dplyr","reshape","clustsig",
                  "cowplot","ggrepel","RVAideMemoire")

new_packages <- package_list[!(package_list %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load required packages
#======================================
library(tidyverse) # Plots and dplyr
library(tidyr) # Gather function
library(vegan) # NMDS
library(dplyr) # Wrangling data
library(reshape) # Melting data
library(clustsig) # SIMPROF cluster
library(cowplot) # Combine plots 
library(ggrepel) # Repel text in plots 
library(RVAideMemoire) # Shapiro tests
#======================================

# Organized in order of results presented
# 1: Results section and abundance of taxa
# 2: Long-term abundance and richness trends
# 3: Trajectory (NMDS), alternative states (SIMPROF), and turnover (DCA)

#==========================================================================================================
# 1 - Results section and abundance of taxa
#===========================================

# Read in CPUE data frame 
df1 <- read.csv("spring_fall_combined_data_cpue.csv", header = T, row.names = 1)

# Check structure and summary
str(df1)
summary(df1)

# Transform data from long to short 
abundance_df <- gather(df1, key = "variable", value = "CPUE", -season)

# Plot
ggplot(abundance_df, aes(x = reorder(variable, -CPUE), CPUE)) +
  geom_boxplot()
  
# Subset 10 species with highest CPUE
df2 <- df1 %>%
  select(CYCA,XYTE,ONMY,
         MISA,MOSA,ICPU,
         MIDO,LEMA,AMNA,DOCE)

# Melt data to long format
cpue_melt <- melt(df2)

# Boxplot of abundance
p1 <- ggplot(cpue_melt, aes(x=reorder(variable,-value,FUN = mean),y=value)) + # Reorder from highest to lowest CPUE
  geom_boxplot(color = "black", size = 0.5,
               outlier.colour = "black", 
               fill = "gray",
               outlier.size = 1,
               outlier.shape = 1) +
  theme_classic() + 
  theme(axis.line.x = element_line(colour = 'black', size = 0.1, linetype = 'solid'), 
        axis.line.y = element_line(colour = 'black', size = 0.1, linetype = 'solid'),
        axis.text.x = element_text(colour = "black", size = 15, angle = 40, hjust = 1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.title.y = element_text(colour = "black", size = 16)) +
  labs(y = "Mean CPUE")

p1

#-------------------------------------------------
# Stacked plot of non-native and native abundance
#-------------------------------------------------

# Read in csv 
df3 <- read.csv("native_non_abundance.csv", header = T)

# Calculate relative abundance
comm_abund <- sweep(df3[,2:3], 1, rowSums(df3), '/')*100

# Bind datasets 
comm_abund <- cbind(comm_abund,df3$year)
colnames(comm_abund)[3] <- "year"

# Melt data
abun_melt <- melt(comm_abund, id.vars=c("year"))

# Stacked line graph of native and nonnative abundance
p2 <- ggplot(abun_melt, aes(x=year,y=value,group=variable,fill=variable)) + 
  geom_area(position="fill") +
  scale_fill_brewer(palette="Greys", breaks=rev(levels(abun_melt$variable))) +
  theme_classic() + 
  theme(axis.line.x = element_line(colour = 'black', size = 0.1, linetype = 'solid'), 
        axis.line.y = element_line(colour = 'black', size = 0.1, linetype = 'solid'),
        axis.title = element_text(size = 15, colour = "black"),
        legend.text = element_text(size = 15, colour = "black"),
        legend.title = element_blank(),
        axis.text.x = element_text(colour = "black", size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 15)) + 
  scale_x_continuous(breaks = c(1980,1985,1990,1995,2000,2005,2010,2015,2020)) +
  labs(y = "Proportion of catch")

p2

# Combine plots
catch <- plot_grid(p2, p1, ncol = 1, labels = c("a","b"))

catch

# Save image
# ggsave("Figure X.tiff", width = 4, height = 7, units = 'in', res = 500)

#==========================================================================================================
# 2 - Long-term abundance and richness trends
#=============================================

# Check distributions
lapply(df1[ ,2:19], shapiro.test)

# Long-term abundance trends df
df4 <- cbind(df3$year, df1)
colnames(df4)[1] <- "year"

# Correlation matrix of species
cor_df <- as.data.frame(cor(df4[,3:20]))

# Create function for plots 
abundance_plot <- function(species, title){

  plot <- ggplot(df4, aes(x = year, y = species)) +
    geom_point(colour = "black", size = 2) +
    theme_classic() + 
    theme(axis.line.x = element_line(colour = 'black', size = 0.1, linetype = 'solid'), 
          axis.line.y = element_line(colour = 'black', size = 0.1, linetype = 'solid'),
          axis.title = element_text(size = 12, colour = "black"),
          legend.title = element_blank(),
          axis.text.x = element_text(colour = "black", size = 12, angle = 30, hjust = 1),
          axis.text.y = element_text(colour = "black", size = 12)) +
    labs(y = title, x = "Year")
  
  return(plot)
}

# Correlations (abundance x time) 
for (i in 3:length(df4)){
  a <- cor.test(df4$year, df4[,i], method = "spearman")
  print(paste(colnames(df4)[i], " rho", a$estimate, " p-value:", a$p.value))
}

###############################
# Species abundance trend plots
###############################

# Dorosoma cepedianum
print(DOCE_plot <- abundance_plot(df4$DOCE, expression(paste(italic("Dorosoma cepedianum"), " CPUE"))))

# Dorosoma petenense
print(DOPE_plot <- abundance_plot(df4$DOPE, expression(paste(italic("Dorosoma petenense"), " CPUE"))))

# Oncorhynchus clarkii henshaw
print(ONCLHE_plot <- abundance_plot(df4$ONCLHE, expression(paste(italic("Oncorhynchus clarkii henshawi"), " CPUE"))))

# Oncorhynchus mykiss
print(ONMY_plot <- abundance_plot(df4$ONMY, expression(paste(italic("Oncorhynchus mykiss"), " CPUE"))))

# Salvelinus fontinalis
print(SAFO_plot <- abundance_plot(df4$SAFO, expression(paste(italic("Salvelinus fontinalis"), " CPUE"))))

# Gila elegans
print(GIEL_plot <- abundance_plot(df4$GIEL, expression(paste(italic("Gila elegans"), " CPUE"))))

# Cyprinus carpio
print(CYCA_plot <- abundance_plot(df4$CYCA, expression(paste(italic("Cyprinus carpio"), " CPUE"))))

# Xyrauchen texanus
print(XYTE_plot <- abundance_plot(df4$XYTE, expression(paste(italic("Xyrauchen texanus"), " CPUE"))))

# Ameiurus melas
print(AMME_plot <- abundance_plot(df4$AMME, expression(paste(italic("Ameiurus melas"), " CPUE"))))

# Ameiurus natalis
print(AMNA_plot <- abundance_plot(df4$AMNA, expression(paste(italic("Ameiurus natalis"), " CPUE"))))

# Ictalurus punctatus
print(ICPU_plot <- abundance_plot(df4$ICPU, expression(paste(italic("Ictalurus punctatus"), " CPUE"))))

# Lepomis cyanellus
print(LECY_plot <- abundance_plot(df4$LECY, expression(paste(italic("Lepomis cyanellus"), " CPUE"))))

# Lepomis macrochirus
print(LEMA_plot <- abundance_plot(df4$LEMA, expression(paste(italic("Lepomis macrochirus"), " CPUE"))))

# Micropterus dolomieu
print(MIDO_plot <- abundance_plot(df4$MIDO, expression(paste(italic("Micropterus dolomieu"), " CPUE"))))

# Micropterus salmoides
print(MISA_plot <- abundance_plot(df4$MISA, expression(paste(italic("Micropterus salmoides"), " CPUE"))))

# Pomoxis nigromaculatus
print(PONI_plot <- abundance_plot(df4$PONI, expression(paste(italic("Pomoxis nigromaculatus"), " CPUE"))))

# Morone saxatilis
print(MOSA_plot <- abundance_plot(df4$MOSA, expression(paste(italic("Morone saxatilis"), " CPUE"))))

#########################
# Total CPUE and richness
#########################

# Read in csv
df5 <- read.csv("cpue_richness.csv", header = T)

# Check for normal distribution
lapply(df5[,2:3], shapiro.test)

# CPUE plot
cpue <- ggplot(df5, aes(x = year, y = cpue)) +
  geom_point(colour = "black", size = 2) +
  theme_classic() + 
  theme(axis.line.x = element_line(colour = 'black', size = 1.0, linetype = 'solid'), 
        axis.line.y = element_line(colour = 'black', size = 1.0, linetype = 'solid'),
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),
        legend.text = element_text(size = 15, colour = "black"),
        legend.title = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  labs(y = "Mean CPUE", x = "Year")

cpue

# Richness plot
richness <- ggplot(df5, aes(x = year, y = richness)) +
  geom_point(colour = "black", size = 2) +
  theme_classic() + 
  theme(axis.line.x = element_line(colour = 'black', size = 1.0, linetype = 'solid'), 
        axis.line.y = element_line(colour = 'black', size = 1.0, linetype = 'solid'),
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 15, angle = 30, hjust = 1, colour = "black"),
        axis.title = element_text(size = 15, colour = "black"),
        legend.text = element_text(size = 15, colour = "black"),
        legend.title = element_blank(),
        legend.position = "bottom") +
  labs(y = "Total species richness", x = "Year") +
  scale_x_continuous(breaks = c(1980,1985,1990,1995,2000,2005,2010,2015,2020))

richness

# Plot together
plot_grid(cpue,richness,nrow = 2, labels = c("a","b"))

# Correlations
cor.test(df5$year, df5$richness, method = "spearman")

cor.test(df5$year, df5$cpue, method = "spearman")

#==========================================================================================================
# 3 - Fish community trajectory, alternative states, and turnover 
#================================================================

# Create new dataframe and remove cutthroat spp.
df6 <- df1 %>%
  select(-Cutthroat.spp.)

# Remove species that occur in less than 5% of samples 
df7 <- df6[, colSums(df6[,2:18] != 0) > 1]

# Fourth-root transformation
df7_tran <- sqrt(sqrt(df7[,2:17]))

# Establish a seed so results are consistent
set.seed(101)

# NMDS ordination
ord <- metaMDS(df7_tran,distance = "bray", k = 2, try = 20, trymax = 1000, autotransform = FALSE, wascores = TRUE, plot = TRUE)

# Start from previous best solution
ord_best <- metaMDS(df7_tran,distance = "bray", k = 2, try = 20, trymax = 1000, autotransform = FALSE, wascores = TRUE, plot = TRUE, previous.best = ord)

print(ord_best)

stressplot(ord_best)

# Goodness of fit
gof <- goodness(ord_best)
plot(ord_best, display = "sites", type = "n")
points(ord_best, display = "sites", cex = 2*gof/mean(gof))

#######################
# NMDS trajectory plot 
#######################

# Extract NMDS scores 
fish_scores <- as.data.frame(scores(ord_best))

# Arrows
fish_scores$d <- NA
fish_scores$d[40] <- fish_scores$NMDS1[40]
fish_scores$e <- NA
fish_scores$e[40] <- fish_scores$NMDS2[40]


p3 <- ggplot(data = fish_scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point() +
  theme_classic() +
  geom_text_repel(label=rownames(fish_scores), size = 6) + # Repel year labels if overlapping
  geom_segment(aes(xend=c(tail(NMDS1,n=-1),NA), yend=c(tail(NMDS2,n=-1),NA)),) + # Add trajectory arrow
  geom_segment(aes(xend=c(tail(d,n=-1),NA), yend=c(tail(e,n=-1),NA)),
               arrow = arrow(length=unit(0.6,"cm"), type = "closed")) +
  theme(axis.text = element_text(size = 22, colour = "black"),
        axis.title = element_text(size = 22, colour = "black"),
        legend.text = element_text(size = 22),
        legend.title = element_blank()) 
p3

# Save image
# ggsave("Figure X.tiff", width = 5, height = 5, units = 'in', res = 500)

################################################
# SIMPROF cluster analysis - alternative states
################################################

# Simprof cluster
sim <- simprof(df7_tran, num.expected=2000, num.simulated=999,
              method.cluster="average", method.distance="braycurtis",
              method.transform="identity", alpha=0.01,
              sample.orientation="row", const=0,
              silent=TRUE, increment=100,
              undef.zero=TRUE)

simprof.plot(sim, leafcolors=NA, plot=TRUE, fill=TRUE,
             leaflab="perpendicular", siglinetype=1)

#############################################################################
# NMDS plot with convex hulls based on SIMPROF (alternative states)
#############################################################################

# Create list of years for convex hulls 
year <- as.data.frame(unlist(sim$significantclusters))
colnames(year) <- "year"

# Group as factor
group <- as.data.frame(factor(c(rep("one", 10), rep("two", 15), rep("three",15))))
colnames(group) <- "group"

# Create dataframe
sim_df <- cbind(fish_scores[,1:2], year, group)

# Create convex hulls 
hull1 <- fish_scores[sim_df$group == "one",]

hull2 <- fish_scores[sim_df$group == "two",]

hull3 = fish_scores[sim_df$group == "three",]

# Plot
p4 <- ggplot(data = sim_df, aes(x = NMDS1, y = NMDS2)) +
  geom_point() +
  theme_classic() +
  geom_text_repel(label=rownames(sim_df), size = 3.5) +
  geom_polygon(inherit.aes=F, data=hull1[chull(hull1),], aes(x=NMDS1, y=NMDS2), 
               alpha=0.3, fill="gray", linetype=1, color="black") + 
  geom_polygon(inherit.aes=F, data=hull2[chull(hull2),], aes(x=NMDS1, y=NMDS2), 
               alpha=0.3, fill="gray", linetype=1, color="black") +
  geom_polygon(inherit.aes=F, data=hull3[chull(hull3),], aes(x=NMDS1, y=NMDS2), 
               alpha=0.3, fill="gray", linetype=1, color="black") +
  theme(axis.line.x = element_line(colour = 'black', size = 1.0, linetype = 'solid'), 
        axis.line.y = element_line(colour = 'black', size = 1.0, linetype = 'solid'),
        axis.text = element_text(size = 15, colour = "black"),
        axis.title = element_text(size = 15, colour = "black")) +
  geom_vline(xintercept = 0.0, linetype = 2,colour = "black") +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "black") 

p4

###################################
# NMDS with species fit as vectors 
###################################

# Extract species vectors 
ef <- envfit(ord_best$points, df7_tran, perm = 999)

ef

# Create dataframe
ef.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r))

# Make species as rownames
ef.df$species=rownames(ef.df)

# Only significant pvalues
# Shortcutting ef$vectors
A <- as.list(ef$vectors)

# Create dataframe
pvals <- as.data.frame(A$pvals)
arrows <- as.data.frame(A$arrows*sqrt(A$r))
C <- cbind(arrows, pvals)

# Subset
Cred <- subset(C,pvals<0.05)
Cred <- cbind(Cred, Species = rownames(Cred))

# Multiplier for arrow length 
mult <- 0.5

p5 <- ggplot(data = sim_df, aes(x = NMDS1, y = NMDS2)) +
  theme_classic() +
  geom_polygon(inherit.aes=F, data=hull1[chull(hull1),], aes(x=NMDS1, y=NMDS2), 
               alpha=0.3, fill=NA, linetype=1, color="black") + 
  geom_polygon(inherit.aes=F, data=hull2[chull(hull2),], aes(x=NMDS1, y=NMDS2), 
               alpha=0.3, fill=NA, linetype=1, color="black") +
  geom_polygon(inherit.aes=F, data=hull3[chull(hull3),], aes(x=NMDS1, y=NMDS2), 
               alpha=0.3, fill=NA, linetype=1, color="black") +
  geom_segment(data = Cred, aes(x = 0, xend = mult*MDS1, y = 0, yend = mult*MDS2), 
               arrow = arrow(length = unit(0.2, "cm")), size = 0.7,colour = "black") +
  geom_text_repel(data = Cred, aes(x=mult*MDS1, y=mult*MDS2, label = Species), size = 4, colour = "black") +
  theme(axis.line.x = element_line(colour = 'black', size = 1.0, linetype = 'solid'), 
        axis.line.y = element_line(colour = 'black', size = 1.0, linetype = 'solid'),
        axis.text = element_text(size = 15, colour = "black"),
        axis.title = element_text(size = 15, colour = "black")) 

p5

# Plot together 
nmds_grid <- plot_grid(p4,p5, ncol = 1, labels = c("a","b"))

nmds_grid

# Save image
# ggsave("Figure X.tiff", width = 5, height = 7, units = 'in', res = 500)

# DCA - turnover
dca <- decorana(df7_tran, iweigh=0, ira=0, before=NULL, after=NULL)
summary(dca)
plot(dca, display = c("sites"))



