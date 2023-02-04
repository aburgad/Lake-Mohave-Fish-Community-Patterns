#####################################
#####################################
# Lake Mohave Fish Community Dynamics
#####################################
#####################################

# Install required packages if not already
package_list <- c("tidyverse","vegan",
                  "reshape","clustsig",
                  "cowplot","ggrepel",
                  "RVAideMemoire")

new_packages <- package_list[!(package_list %in% installed.packages()[,"Package"])]
if(length(new_packages)) 
  install.packages(new_packages)

# Load required packages
#======================================
library(tidyverse) # Plots and dplyr
library(vegan) # NMDS
library(reshape) # Melting data
library(clustsig) # SIMPROF cluster
library(cowplot) # Combine plots 
library(ggrepel) # Repel text in plots 
library(RVAideMemoire) # Shapiro tests
#======================================

# Organized in order of results
# 1: Results section and abundance of taxa
# 2: Long-term abundance and richness trends
# 3: Trajectory (NMDS), alternative states (SIMPROF), and turnover (DCA)

#===========================================
# 1 - Results section and abundance of taxa
#===========================================

df <- read.csv(file = "mohave_roundup_catch_data.csv")

tot_sp_catch <- df %>%
  select(-year, -season, 
         -date_begin, -date_end, 
         -net_units) %>%
  colSums(na.rm = FALSE) %>%
  data.frame() 

colnames(tot_sp_catch) <- "n"

tot_sp_catch$relative_abundance <- round((tot_sp_catch$n/sum(tot_sp_catch$n))*100, 3)

# Summary of total individuals 
sum(tot_sp_catch$n)

# Calculate CPUE for each species 
for (i in 1:ncol(df[,6:22])){
 df[,6:22][,i] <- df[,6:22][,i] / df$net_units
}

total_cpue_sp <- df %>%
  select(-year, -season, 
         -date_begin, -date_end, 
         -net_units) %>%
  colSums(na.rm = FALSE) %>%
  data.frame() 

mean_cpue_sp <- df %>%
  select(-year, -season, 
          -date_begin, -date_end, 
          -net_units) %>%
  colMeans(na.rm = FALSE) %>%
  data.frame() 

# Sum across years 
total_cpue <- df %>%
  select(-season,-date_begin,-date_end) %>%
  group_by(year) %>%
  summarise_all(sum)

num_samp <- df %>%
  group_by(year) %>%
  tally()

total_cpue <- merge(num_samp, total_cpue, by = "year")

# Function to calculate average CPUE
average_cpue <- function(){  
  if(total_cpue$n[i] > 1)
for (i in 1:ncol(total_cpue[,4:20])) {
  total_cpue[,4:20][i] <- round(total_cpue[,4:20][i] / total_cpue$n, 3)
  }
  return(total_cpue)
}

df1 <- average_cpue()

# Check structure and summary
str(df1)
summary(df1)

# Transform data from long format to short 
abundance_df <- gather(df1, key = "variable", value = "CPUE", -year,-n,-net_units)

# Plot
ggplot(abundance_df, aes(x = reorder(variable, -CPUE), CPUE)) +
  geom_boxplot()
  
# Subset 10 species with highest CPUE (core species)
core_df <- df1 %>%
  select(CYCA,XYTE,ONMY,
         MISA,MOSA,ICPU,
         MIDO,LEMA,AMNA,DOCE)

# Melt data to long format
cpue_melt <- melt(core_df)

# Boxplot of abundance
p1 <- ggplot(cpue_melt, aes(x=reorder(variable,-value,FUN = mean),y=value)) + # Reorder from highest to lowest CPUE
  geom_boxplot(color = "black", size = 0.5,
               outlier.colour = "black", 
               fill = "gray",
               outlier.size = 1,
               outlier.shape = 1) +
  theme_classic() + 
  theme(axis.line.x = element_line(colour = 'black', size = 1.0, linetype = 'solid'), 
        axis.line.y = element_line(colour = 'black', size = 1.0, linetype = 'solid'),
        axis.text.x = element_text(colour = "black", size = 15, angle = 40, hjust = 1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.title.y = element_text(colour = "black", size = 16)) +
  labs(y = "Mean CPUE")

print(p1)

# Stacked plot of non-native and native abundance
#-------------------------------------------------

comm_abund <- df %>%
  select(-season,-date_begin,-date_end,-net_units) %>%
  gather(key = "species", value = "n", -year) %>%
  mutate(native = ifelse(species == "XYTE" | species == "GIEL", "Native", "Nonnative")) %>%
  group_by(native,year) %>%
  summarise(total = sum(n)) %>%
  spread(native,total)

# Calculate relative abundance for native and nonnative
rel_abund <- comm_abund %>%
  mutate(Total = Native + Nonnative,
         Native = Native/Total*100,
         Nonnative = Nonnative/Total*100) %>%
  gather(key = "variable", value = "value", -year,-Total)

# Stacked line graph of native and nonnative abundance
p2 <- ggplot(rel_abund, aes(x=year,y=value,group=variable,fill=variable)) + 
  geom_area(position="fill") +
  theme_classic() + 
  scale_fill_manual(values = c("gray85","gray40")) + 
  theme(axis.line.x = element_line(colour = 'black', size = 1.0, linetype = 'solid'), 
        axis.line.y = element_line(colour = 'black', size = 1.0, linetype = 'solid'),
        axis.title = element_text(size = 15, colour = "black"),
        legend.text = element_text(size = 15, colour = "black"),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text.x = element_text(colour = "black", size = 15, angle = 30, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 15)) + 
  scale_x_continuous(breaks = c(1980,1985,1990,1995,2000,2005,2010,2015,2020)) +
  labs(y = "Proportion of catch")

# Combine plots
catch <- plot_grid(p2, p1, ncol = 1, labels = c("a","b"))

print(catch)

# Save image
ggsave("Figure 2.tiff", device = 'tiff',width = 4, height = 7, units = 'in', dpi = 500)

#=============================================
# 2 - Long-term abundance and richness trends
#=============================================

# Check distributions
lapply(df1[ ,4:20], shapiro.test)

# Correlation matrix of species
cor_df <- as.data.frame(cor(df1[,4:20]), method = 'spearman')

# Plot of species correlations
cor_df$species = row.names(cor_df)
cor_df2 = melt(cor_df, id = "species")
cor_df2$species <- factor(cor_df2$species,levels=unique(cor_df2$species))

ggplot(cor_df2, aes(x = variable, y = species)) + 
  geom_tile(aes(fill = value), colour = "grey45") + 
  coord_equal() + 
  scale_fill_gradient(low = "gray75", high = "black") + 
  theme(axis.text.y = element_text(face = "bold", colour = "grey25"), 
        legend.title = element_text(size = 10, face = "bold"),legend.position = "right", 
        axis.text.x = element_text(angle = 90, face = "bold",colour = "grey25", vjust = 0.5, hjust = 0), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = NA), 
        axis.ticks = element_blank()) + 
  labs(x= "", y = "", fill = "Spearman's Correlation") + 
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev(levels(cor_df2$species))) 

# Spearman's correlation output (abundance x time) 
for (i in 3:length(df1)){
  a <- cor.test(df1$year, df1[,i], method = "spearman")
  print(paste(colnames(df1)[i], " rho", a$estimate, " p-value:", a$p.value))
}


# Species abundance individual trend plots
#------------------------------------------

# Create function for plots 
abundance_plot <- function(species, title){
  
  plot <- ggplot(df1, aes(x = year, y = species)) +
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

# Dorosoma cepedianum
print(abundance_plot(df1$DOCE, expression(paste(italic("Dorosoma cepedianum"), " CPUE"))))

# Dorosoma petenense
print(abundance_plot(df1$DOPE, expression(paste(italic("Dorosoma petenense"), " CPUE"))))

# Oncorhynchus clarkii henshaw
print(abundance_plot(df1$ONCLHE, expression(paste(italic("Oncorhynchus clarkii henshawi"), " CPUE"))))

# Oncorhynchus mykiss
print(abundance_plot(df1$ONMY, expression(paste(italic("Oncorhynchus mykiss"), " CPUE"))))

# Salvelinus fontinalis
print(abundance_plot(df1$SAFO, expression(paste(italic("Salvelinus fontinalis"), " CPUE"))))

# Gila elegans
print(abundance_plot(df1$GIEL, expression(paste(italic("Gila elegans"), " CPUE"))))

# Cyprinus carpio
print(abundance_plot(df1$CYCA, expression(paste(italic("Cyprinus carpio"), " CPUE"))))

# Xyrauchen texanus
print(abundance_plot(df1$XYTE, expression(paste(italic("Xyrauchen texanus"), " CPUE"))))

# Ameiurus melas
print(abundance_plot(df1$AMME, expression(paste(italic("Ameiurus melas"), " CPUE"))))

# Ameiurus natalis
print(abundance_plot(df1$AMNA, expression(paste(italic("Ameiurus natalis"), " CPUE"))))

# Ictalurus punctatus
print(abundance_plot(df1$ICPU, expression(paste(italic("Ictalurus punctatus"), " CPUE"))))

# Lepomis cyanellus
print(abundance_plot(df1$LECY, expression(paste(italic("Lepomis cyanellus"), " CPUE"))))

# Lepomis macrochirus
print(abundance_plot(df1$LEMA, expression(paste(italic("Lepomis macrochirus"), " CPUE"))))

# Micropterus dolomieu
print(abundance_plot(df1$MIDO, expression(paste(italic("Micropterus dolomieu"), " CPUE"))))

# Micropterus salmoides
print(abundance_plot(df1$MISA, expression(paste(italic("Micropterus salmoides"), " CPUE"))))

# Pomoxis nigromaculatus
print(abundance_plot(df1$PONI, expression(paste(italic("Pomoxis nigromaculatus"), " CPUE"))))

# Morone saxatilis
print(abundance_plot(df1$MOSA, expression(paste(italic("Morone saxatilis"), " CPUE"))))


# Species trend plot
#-------------------------

# Setting independent y-axis scales due to large variations in CPUE 

trend_plot <- df1 %>%
  gather(key = "species", value = "CPUE", -year, -n, -net_units) %>%
  ggplot(., aes(x = year, y = CPUE)) +
  geom_point(size = 1) +
  facet_wrap(~species,scales="free") +
  #geom_smooth(colour = 'black',se=FALSE) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 75, hjust=1,colour = "black"),
        axis.text.y = element_text(colour = "black")) +
  labs(x = 'Year', y = 'CPUE')

print(trend_plot)

# Total CPUE and richness
-------------------------

# Sum CPUE (appears to be a bug that requires highlighting code to run)
cpue_df <- df1 %>%
  select(-n,-net_units) %>%
  rowwise() %>%
  mutate(cpue = sum(c_across(DOCE:MOSA))) %>%
  select(year,cpue)

# Calculate richness
rich_df <- df1 %>%
  select(-year,-n,-net_units) %>%
  mutate_if(is.numeric, ~1 * (. > 0)) %>%
  rowwise() %>%
  mutate(richness = sum(c_across(DOCE:MOSA))) %>%
  select(richness)

# Combine columns
cpue_rich_df <- cbind(cpue_df,rich_df)

# Check for normal distribution
lapply(cpue_rich_df[,2:3], shapiro.test)

# CPUE plot
cpue <- ggplot(cpue_rich_df, aes(x = year, y = cpue)) +
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
  labs(y = "Total CPUE", x = "Year")

# Richness plot
richness <- ggplot(cpue_rich_df, aes(x = year, y = richness)) +
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

# Plot together
cpue_rich_p <- plot_grid(cpue,richness,nrow = 2, labels = c("a","b"))

print(cpue_rich_p)

# Save image
ggsave("Figure 3.tiff", device = 'tiff',width = 6, height = 7, units = 'in', dpi = 500)

# Mean and SD 
mean(cpue_rich_df$richness)
sd(cpue_rich_df$richness)

mean(cpue_rich_df$cpue)
sd(cpue_rich_df$cpue)

# Correlations
cor.test(cpue_rich_df$year, cpue_rich_df$richness, method = "spearman")

cor.test(cpue_rich_df$year, cpue_rich_df$cpue, method = "spearman")

#================================================================
# 3 - Fish community trajectory, alternative states, and turnover 
#================================================================

# Remove species that occur in less than 5% of samples for multivariate analysis
comm_df <- df1[,4:20][, colSums(df1[,4:20] != 0) > 1]

# Fourth-root transformation
comm_df_tran <- sqrt(sqrt(comm_df))

row.names(comm_df_tran) <- df1$year

# Establish a seed so results are consistent
set.seed(101)

# NMDS ordination
ord <- metaMDS(comm_df_tran,distance = "bray", k = 2, 
               try = 20, trymax = 1000, 
               autotransform = FALSE, 
               wascores = TRUE, plot = TRUE)

# Start from previous best solution
ord_best <- metaMDS(comm_df_tran,distance = "bray", k = 2, 
                    try = 20, trymax = 1000, 
                    autotransform = FALSE, wascores = TRUE, 
                    plot = TRUE, previous.best = ord)

print(ord_best)

stressplot(ord_best)

# Goodness of fit
gof <- goodness(ord_best)
plot(ord_best, display = "sites", type = "n")
points(ord_best, display = "sites", cex = 2*gof/mean(gof))

# NMDS trajectory plot 
#----------------------

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
  scale_x_continuous(limits = c(-0.8, 0.6)) +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title = element_text(size = 18, colour = "black"),
        legend.text = element_text(size = 18),
        legend.title = element_blank()) 

print(p3)

# Save image
ggsave("Figure 4.tiff", width = 8, height = 5, units = 'in', dpi = 500)


# SIMPROF cluster analysis - alternative states
#-----------------------------------------------

# Simprof cluster
sim <- simprof(comm_df_tran, num.expected=2000, num.simulated=999,
              method.cluster="average", method.distance="braycurtis",
              method.transform="identity", alpha=0.01,
              sample.orientation="row", const=0,
              silent=FALSE, increment=100,
              undef.zero=TRUE)

simprof.plot(sim, leafcolors=NA, plot=TRUE, fill=TRUE,
             leaflab="perpendicular", siglinetype=1)

# NMDS plot with convex hulls based on SIMPROF (alternative states)
#-------------------------------------------------------------------

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
  scale_x_continuous(limits = c(-0.8, 0.6)) +
  theme(axis.line.x = element_line(colour = 'black', size = 1.0, linetype = 'solid'), 
        axis.line.y = element_line(colour = 'black', size = 1.0, linetype = 'solid'),
        axis.text = element_text(size = 15, colour = "black"),
        axis.title = element_text(size = 15, colour = "black")) +
  geom_vline(xintercept = 0.0, linetype = 2,colour = "black") +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "black") 

print(p4)


# NMDS with species fit as vectors 
#----------------------------------

# Extract species vectors 
ef <- envfit(ord_best$points, comm_df_tran, perm = 999)

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
  scale_x_continuous(limits = c(-0.8, 0.6)) +
  theme(axis.line.x = element_line(colour = 'black', size = 1.0, linetype = 'solid'), 
        axis.line.y = element_line(colour = 'black', size = 1.0, linetype = 'solid'),
        axis.text = element_text(size = 15, colour = "black"),
        axis.title = element_text(size = 15, colour = "black")) 

# Plot together 
nmds_grid <- plot_grid(p4,p5, ncol = 1, labels = c("a","b"))

print(nmds_grid)

# Save image
ggsave("Figure 5.tiff", width = 5, height = 7, units = 'in', dpi = 500)


# DCA - turnover
#---------------

dca <- decorana(comm_df_tran, iweigh=0, ira=0, before=NULL, after=NULL)
summary(dca)
plot(dca, display = c("sites"))



