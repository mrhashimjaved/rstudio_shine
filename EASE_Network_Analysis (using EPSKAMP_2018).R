
#removing old variables if any existing
rm(list=ls())

#✅ 1. Prepare Your Data
#Set working directory
setwd("D:/GIHD DATA/Analysis/Network Analysis")

#read datafile into R
library(foreign)
#ease_df <- read.spss("ease.sav", to.data.frame = TRUE)
#ease_df <- read.spss("2. EASE (n=5856).sav", to.data.frame = TRUE)

#selecting PSC variables for analysis
#psc_df1 <- ease_df[,c(12,19:53,104:105)]
#psc_df1 <- ease_df[,c(9,14:48,93:94)]
psc_df1 <- read.spss("2. EASE (n=5233) (PSC Nvalid=35).sav", to.data.frame = TRUE)

library("tidyverse")
#psc_df2 <- psc_df1 %>% filter(psc_total_cat=='At Distress')
psc_df2 <- psc_df1[,1:35]
#psc_df2 <- psc_df1[,2:36]


#removing cases with missing values in data
psc_df2 <- na.omit(psc_df2)
#renaming psc_1 to psc_01.....
names(psc_df2)[names(psc_df2) == "psc_1"] <- "psc_01"
names(psc_df2)[names(psc_df2) == "psc_2"] <- "psc_02"
names(psc_df2)[names(psc_df2) == "psc_3"] <- "psc_03"
names(psc_df2)[names(psc_df2) == "psc_4"] <- "psc_04"
names(psc_df2)[names(psc_df2) == "psc_5"] <- "psc_05"
names(psc_df2)[names(psc_df2) == "psc_6"] <- "psc_06"
names(psc_df2)[names(psc_df2) == "psc_7"] <- "psc_07"
names(psc_df2)[names(psc_df2) == "psc_8"] <- "psc_08"
names(psc_df2)[names(psc_df2) == "psc_9"] <- "psc_09"

#Creating Legend by defining groups of variables
psc_gps <- scan("PSC_Groups.txt", what="character", sep="\n")

#Creating Item descriptions
psc_items <- scan("PSC_Items.txt", what="character", sep="\n")

############################################################
#GPT Guided code form (EPSKAMP et al, 2018)           ######
############################################################

psc_data <- psc_df2

#✅ 2. Estimate Network (with EBICglasso)
library(bootnet)
library(qgraph)

# 🔹 Step 1: Estimate the network using regularized partial correlations
network <- estimateNetwork(psc_data, default = "EBICglasso", corMethod = "cor_auto")

#🔹 Step 2: Plot and Save the qgraph Object
    # Plot and save the qgraph object
network_plot <- qgraph(network$graph, layout = "spring", labels = colnames(psc_data), theme = "colorblind")

#🔹 Step 3: Compute Expected Influence
ei <- expectedInfluence(network_plot)
  #If you want expected influence 2-step (i.e., including indirect connections), do:
ei_2step <- expectedInfluence(network_plot, order = 2)

#🔹 Step 4: View and Sort the Results
# View EI
print(ei)

# Top 5 most central items
sort(ei, decreasing = TRUE)[1:5]



#✅ 3. Visualize the Network
# Plot with labels and better layout
qgraph(network$graph, layout = "spring", labels = colnames(psc_data), theme = "colorblind")

#✅ 4. Analyze Centrality
centralityPlot(network)
#🔹 Option 1: Get Full Centrality Metrics Manually (Separately)
centrality <- centrality(network)
centrality$Betweenness
centrality$Closeness
centrality$Strength  # Already plotted
#Option 1: Get Full Centrality Metrics Manually (together)
centralityPlot(network, include = c("Strength", "Closeness", "Betweenness"))

#🔹 Option 2: Use Expected Influence (Recommended!)
EI <- centrality_auto(network, include = "expectedInfluence")
EI$expectedInfluence
#Plot it:
centralityPlot(network, include = "ExpectedInfluence")


#✅ 5. Test Stability with Bootstrapping
boot_results <- bootnet(network, nBoots = 1000, type = "nonparametric")

# Plot confidence intervals around edge weights
plot(boot_results, "edge")

# Centrality stability
plot(boot_results, "strength")



