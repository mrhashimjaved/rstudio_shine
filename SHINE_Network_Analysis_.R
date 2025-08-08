
#removing old variables if any existing
rm(list=ls())

#Set working directory
setwd("D:/GIHD DATA/Analysis/Network Analysis")

#read datafile into R
library(foreign)
shine_sdq <- read.spss("1. SHINE_TSDQ (n=9685).sav", to.data.frame = TRUE)

#selecting PSC variables for analysis
sdq_df1 <- shine_sdq[,c(31:55,56:62)]

library("tidyverse")
sdq_df2 <- sdq_df1[,1:25]


#removing cases with missing values in data
sdq_df2 <- na.omit(sdq_df2)
#renaming psc_1 to psc_01.....
names(sdq_df2)[names(sdq_df2)== "sdq_1"] <- "sdq_01"
names(sdq_df2)[names(sdq_df2)== "sdq_2"] <- "sdq_02"
names(sdq_df2)[names(sdq_df2)== "sdq_3"] <- "sdq_03"
names(sdq_df2)[names(sdq_df2)== "sdq_4"] <- "sdq_04"
names(sdq_df2)[names(sdq_df2)== "sdq_5"] <- "sdq_05"
names(sdq_df2)[names(sdq_df2)== "sdq_6"] <- "sdq_06"
names(sdq_df2)[names(sdq_df2)== "sdq_7"] <- "sdq_07"
names(sdq_df2)[names(sdq_df2)== "sdq_8"] <- "sdq_08"
names(sdq_df2)[names(sdq_df2)== "sdq_9"] <- "sdq_09"
#Creating Legend by defining groups of variables
sdq_gps <- scan("SDQ_Groups.txt", what="character", sep="\n")

#Creating Item descriptions
SDQ_items <- scan("SDQ_Items.txt", what="character", sep="\n")

library(psych)
poly_values = polychoric(sdq_df2)
polyMat <- poly_values$rho

library(qgraph)
#partial correlation matrix
sdq_corMat <- cor_auto(sdq_df2)

#Unregularized partial correlation netwrok
sdq_pcor <- qgraph(cor(sdq_df2), layout="spring", tuning=0.25, groups=sdq_gps,
                   legend.cex=0.50, vsize=4, esize=15)

#Regularized partial correlation Network
psc_Glasso <- qgraph(sdq_corMat, graph="glasso", layout="spring", tuning=0.25, 
                     sampleSize=nrow(sdq_df2), groups=sdq_gps, nodeNames=sdq_items, 
                     legend.cex=0.20, vsize=4, esize=15, minimum=0.075)

#checking for redundant nodes
library(networktools)
goldbricker(sdq_df2,p=0.05)

#Checking for highly influential nodes in network (Important Symptoms)
centralityPlot(psc_Glasso, include="all", decreasing=TRUE)

#checking for Expected Influence using networktools package
sdq_Glasso <- qgraph(sdq_corMat, graph="glasso", layout="spring", tuning=0.25, sampleSize=nrow(sdq_df2))
psc_expInf <- expectedInf(psc_Glasso)
plot(sdq_expInf)

#Testing for Network accuracy and Centrality stability
library(bootnet)
#i.Network estimation
sdq_net <- estimateNetwork(sdq_df2, default="EBICglasso")
plot(psc_net, layout="spring", labels=TRUE)
centralityPlot(sdq_net, include = c("Betweenness", "Closeness", "Strength"), 
               decreasing=TRUE)

#ii. Network accuracy (Using bootstrapping)
sdq_boot1 <- bootnet(sdq_net, nBoots=2500, nCores=8)
plot(sdq_boot1, labels=FALSE, order="sample")

#iii. Network stability (Using bootstrapping)
sdq_boot2 <- bootnet(sdq_net, nBoots=2500, type="case", nCores=8)
plot(sdq_boot2)
corStability(sdq_boot2)

#iv. Testing for significance differences
differenceTest(sdq_boot1, 22,33,"strength")

#testing between all pairs 
#Warning!!
#Requires a lot of RAM !!

plot(sdq_boot1, "edge", plot="difference", onlyNonZero=TRUE, order="sample")


plot(sdq_boot1)

#Simulation Studies


#Comparing Plots
Layout <- averageLayout(psc_pcor,psc_Glasso)
layout(t(1:2))
sdq_pcor <- qgraph(sdq_corMat, graph="sdq_pcor", layout="spring")
sdq_Glasso <- qgraph(sdq_corMat, graph="glasso", layout="spring", tuning=0.25, sampleSize=nrow(sdq_df2))

#stop split layout
layout(1)

#EGAnet :: Explorartory Graph Analysis (To analyze communities)

library(EGAnet)
EGA(psc_Glasso)

# Eigenvalue Decomposition
plot(eigen(sdq_corMat)$values, type = "b")
abline(h=1, col="red", lty=3)

# SPINGLASS ALGORITHM
library(igraph)
g = as.igraph(sdq_Glasso, attributes=TRUE)
sgc <- spinglass.community(g)
sgc$membership

# Exploratory Graph Analysis
# library("devtools")
# devtools::install_github('hfgolino/EGA')
library("EGA")
ega<-EGA(sdq_df2, plot.EGA = TRUE)
ega$membership

group.spinglass <- list(c(6,15,20,23,27,29,30,31,35), c(1,2,3,4,7,8,11,12,21,22,26)
                        , c(13,19,24), c(9,10,14,17,18), c(5,16,25,28,32,33,34))

group.ega <- list(c(1,2,3,4,8,9,10,11,12,13,21,22,26), c(6,15,17,29,31,35), c(5,14,
                                                                              18,20,23,27,30), c(7,19,24,25,28), c(16,32,33,34))

group.psc <- list(c(1,2,3,5,6,9,10,12,15,17,18,21,23,24,25,26,28,30), c(4,7,8,9,14)
                  , c(11,13,19,22,27), c(16,29,31,32,33,34,35))

par(mfrow = c(1,2))
set.seed(5)
graphEGA <- qgraph(psc_corMat, graph="glasso", layout="spring", tuning=0.25, 
                   sampleSize=nrow(psc_df2), groups=group.ega,  
                   legend.cex=0.20, vsize=4, esize=15, minimum=0.075, title="ega walktrap")
#graph4 <- qgraph(psc_corMat, graph="glasso", layout="spring", sampleSize=nrow(psc_df2), 
#                 groups=group.ega, vsize=7, cut=0, maximum=.45, border.width=1.5,
#                 color=c("red", "green", "blue", "orange", "white"), title="ega walktrap",
#                 layout.par=list(init=matrix(rnorm(nNode*2),nNode,2)))
set.seed(5)
graphSG <- qgraph(psc_corMat, graph="glasso", layout="spring", tuning=0.25, 
                  sampleSize=nrow(psc_df2), groups=group.spinglass,  
                  legend.cex=0.20, vsize=4, esize=15, minimum=0.075, title="igraph spinglass")
#graph4 <- qgraph(psc_corMat, graph="glasso", layout="spring", sampleSize=nrow(psc_df2), 
#                 groups=group.spinglass, vsize=7, cut=0, maximum=.45, border.width=1.5,
#                 color=c("red", "green", "blue", "orange", "white"), title="ega spinglass",
#                 layout.par=list(init=matrix(rnorm(nNode*2),nNode,2)))
set.seed(5)
graphc <- qgraph(psc_corMat, graph="glasso", layout="spring", tuning=0.25, 
                 sampleSize=nrow(psc_df2), groups=group.psc,  
                 legend.cex=0.20, vsize=4, esize=15, minimum=0.075, title="PSC Clusters")

#********************************************************************************************************#
#  Gender Wise  |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#  ===========  |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#  PSC - Male   |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#********************************************************************************************************#

#removing old variables if any existing
rm(list=ls())

#read datafile into R
library(foreign)
ease_df <- read.spss("ease.sav", to.data.frame = TRUE)

#selecting PSC variables for analysis
psc_df1 <- ease_df[,c(12,19:53,104:105)]
library("tidyverse")
psc_df_m <- psc_df1[,1:36] %>% filter(gender=='Male')
psc_df_m <- psc_df_m[,2:36]



#removing cases with missing values in data
psc_df_m <- na.omit(psc_df_m)
#renaming psc_1 to psc_01.....
names(psc_df_m)[names(psc_df_m) == "psc_1"] <- "psc_01"
names(psc_df_m)[names(psc_df_m) == "psc_2"] <- "psc_02"
names(psc_df_m)[names(psc_df_m) == "psc_3"] <- "psc_03"
names(psc_df_m)[names(psc_df_m) == "psc_4"] <- "psc_04"
names(psc_df_m)[names(psc_df_m) == "psc_5"] <- "psc_05"
names(psc_df_m)[names(psc_df_m) == "psc_6"] <- "psc_06"
names(psc_df_m)[names(psc_df_m) == "psc_7"] <- "psc_07"
names(psc_df_m)[names(psc_df_m) == "psc_8"] <- "psc_08"
names(psc_df_m)[names(psc_df_m) == "psc_9"] <- "psc_09"

#Creating Legend by defining groups of variables
psc_gps <- scan("PSC_Groups.txt", what="character", sep="\n")

#Creating Item descriptions
psc_items <- scan("PSC_Items.txt", what="character", sep="\n")

library(qgraph)
#partial correlation matrix
psc_corMat_m <- cor_auto(psc_df_m)

#Regularized partial correlation Network
psc_Glasso_m <- qgraph(psc_corMat_m, graph="glasso", layout="spring", tuning=0.25, 
                       sampleSize=nrow(psc_df_m), groups=psc_gps, nodeNames=psc_items, 
                       legend.cex=0.20, vsize=5, esize=15, minimum=0.075)

#checking for redundant nodes
library(networktools)
goldbricker(psc_df_m,p=0.05)

#Checking for highly influential nodes in network (Important Symptoms)
centralityPlot(psc_Glasso_m, include="all", decreasing=TRUE)

#Testing for Network accuracy and Centrality stability
library(bootnet)
#i.Network estimation
psc_net_m <- estimateNetwork(psc_df_m, default="EBICglasso")
plot(psc_net_m, layout="spring", labels=TRUE)
centralityPlot(psc_net, include = c("Betweenness", "Closeness", "Strength"), 
               decreasing=TRUE)

#ii. Network accuracy (Using bootstrapping)
psc_boot1_m <- bootnet(psc_net_m, nBoots=2500, nCores=8)
plot(psc_boot1_m, labels=FALSE, order="sample")

#iii. Network stability (Using bootstrapping)
psc_boot2_m <- bootnet(psc_net_m, nBoots=2500, type="case", nCores=8)
plot(psc_boot2_m)
corStability(psc_boot2_m)

#iv. Testing for significance differences
differenceTest(psc_boot1_m, 22,33,"strength")



#********************************************************************************************************#
#  Gender Wise  |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#  ===========  |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#  PSC - Female |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#********************************************************************************************************#

psc_df_f <- psc_df1[,1:36] %>% filter(gender=='Female')
psc_df_f <- psc_df_f[,2:36]


#removing cases with missing values in data
psc_df_f <- na.omit(psc_df_f)
#renaming psc_1 to psc_01.....
names(psc_df_f)[names(psc_df_f) == "psc_1"] <- "psc_01"
names(psc_df_f)[names(psc_df_f) == "psc_2"] <- "psc_02"
names(psc_df_f)[names(psc_df_f) == "psc_3"] <- "psc_03"
names(psc_df_f)[names(psc_df_f) == "psc_4"] <- "psc_04"
names(psc_df_f)[names(psc_df_f) == "psc_5"] <- "psc_05"
names(psc_df_f)[names(psc_df_f) == "psc_6"] <- "psc_06"
names(psc_df_f)[names(psc_df_f) == "psc_7"] <- "psc_07"
names(psc_df_f)[names(psc_df_f) == "psc_8"] <- "psc_08"
names(psc_df_f)[names(psc_df_f) == "psc_9"] <- "psc_09"


library(qgraph)
#partial correlation matrix
psc_corMat_f <- cor_auto(psc_df_f)

#Regularized partial correlation Network
psc_Glasso_f <- qgraph(psc_corMat_f, graph="glasso", layout="spring", tuning=0.25, 
                       sampleSize=nrow(psc_df_f), groups=psc_gps, nodeNames=psc_items, 
                       legend.cex=0.20, vsize=5, esize=15, minimum=0.075)

#checking for redundant nodes
library(networktools)
goldbricker(psc_df_f,p=0.05)

#Checking for highly influential nodes in network (Important Symptoms)
centralityPlot(psc_Glasso_f, include="all", decreasing=TRUE)

#Testing for Network accuracy and Centrality stability
library(bootnet)
#i.Network estimation
psc_net_f <- estimateNetwork(psc_df_f, default="EBICglasso")
plot(psc_net_f, layout="spring", labels=TRUE)
centralityPlot(psc_net_f, include = c("Betweenness", "Closeness", "Strength"), 
               decreasing=TRUE)

#ii. Network accuracy (Using bootstrapping)
psc_boot1_f <- bootnet(psc_net_f, nBoots=2500, nCores=8)
plot(psc_boot1_f, labels=FALSE, order="sample")

#iii. Network stability (Using bootstrapping)
psc_boot2_f <- bootnet(psc_net_f, nBoots=2500, type="case", nCores=8)
plot(psc_boot2_f)
corStability(psc_boot2_f)

#iv. Testing for significance differences
differenceTest(psc_boot1_f, 22,33,"strength")




#*******************************************************#
# PSC - Comparing Plots    |||||||||||||||||||||||||||||#
#*******************************************************#

Layout <- averageLayout(psc_Glasso_m,psc_Glasso_f)
layout(t(1:2))
psc_Glasso_m <- qgraph(psc_corMat_m, graph="glasso", layout="spring", tuning=0.25, 
                       sampleSize=nrow(psc_df_m), groups=psc_gps, nodeNames=psc_items, 
                       legend.cex=0.20, vsize=5, esize=15, borders=FALSE,
                       minimum=0.075, legend=FALSE)
psc_Glasso_f <- qgraph(psc_corMat_f, graph="glasso", layout="spring", tuning=0.25, 
                       sampleSize=nrow(psc_df_f), groups=psc_gps, nodeNames=psc_items, 
                       legend.cex=0.20, vsize=5, esize=15, borders=FALSE, 
                       minimum=0.075, legend=FALSE)

layout(1) #stop split layout




#********************************************************************************************************#
#        ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#  SDQ   ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#        ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#********************************************************************************************************#


#selecting SDQ variables for analysis
sdq_df1 <- shine_sdq[,c(31:55,56:62)]
sdq_df2 <- sdq_df1[,1:25]
#library(tidyverse)
#sdq_df_m <- sdq_df1 %>% filter(gender=="Male") #Male
#sdq_df_f <- sdq_df1 %>% filter(gender=="Female") #FeMale


#removing cases with missing values in data
sdq_df2 <- na.omit(sdq_df2)
#renaming sdq_1 to sdq_01.....
names(sdq_df2)[names(sdq_df2) == "sdq_1"] <- "sdq_01"
names(sdq_df2)[names(sdq_df2) == "sdq_2"] <- "sdq_02"
names(sdq_df2)[names(sdq_df2) == "sdq_3"] <- "sdq_03"
names(sdq_df2)[names(sdq_df2) == "sdq_4"] <- "sdq_04"
names(sdq_df2)[names(sdq_df2) == "sdq_5"] <- "sdq_05"
names(sdq_df2)[names(sdq_df2) == "sdq_6"] <- "sdq_06"
names(sdq_df2)[names(sdq_df2) == "sdq_7"] <- "sdq_07"
names(sdq_df2)[names(sdq_df2) == "sdq_8"] <- "sdq_08"
names(sdq_df2)[names(sdq_df2) == "sdq_9"] <- "sdq_09"

#Creating Legend by defining groups of variables
sdq_gps <- scan("D:/GIHD DATA/Analysis/Network Analysis/SDQ_Groups.txt", what="character", sep="\n")

#Creating Item descriptions
sdq_items <- scan("D:/GIHD DATA/Analysis/Network Analysis/sdq_Items.txt", what="character", sep="\n")

library(psych)
poly_values = polychoric(sdq_df2)
polyMat <- poly_values$rho

library(qgraph)
#partial correlation matrix
sdq_corMat <- cor_auto(sdq_df2)

#Plotting Graph with item description
sdq_Glasso <- qgraph(sdq_corMat, graph="glasso", layout="spring", tuning=0.25, 
                     sampleSize=nrow(sdq_df2), groups=sdq_gps, nodeNames=sdq_items, 
                     legend.cex=0.20, vsize=5, esize=15, minimum=0.075)

#checking for redundant nodes
library(networktools)
goldbricker(sdq_df2,p=0.05)

#Checking for highly influential nodes in network (Important Symptoms)
centralityPlot(sdq_Glasso, include="all", decreasing=TRUE)

#checking for Expected Influence using networktools package
sdq_Glasso <- qgraph(sdq_corMat, graph="glasso", layout="spring", tuning=0.25, sampleSize=nrow(sdq_df2))
sdq_expInf <- expectedInf(sdq_Glasso)
plot(sdq_expInf)

#Testing for Network accuracy and Centrality stability
library(bootnet)
#i.Network estimation
sdq_net <- estimateNetwork(sdq_df2, default="EBICglasso")
plot(sdq_net, layout="spring", labels=TRUE)
centralityPlot(sdq_net, include = c("Betweenness", "Closeness", "Strength"), 
               decreasing=TRUE)

#ii. Network accuracy (Using bootstrapping)
sdq_boot1 <- bootnet(sdq_net, nBoots=500, nCores=8)
plot(sdq_boot1, labels=FALSE, order="sample")

#iii. Network stability (Using bootstrapping)
sdq_boot2 <- bootnet(sdq_net, nBoots=2500, type="case", nCores=8)
plot(sdq_boot2)
corStability(sdq_boot2)

#iv. Testing for significance differences
differenceTest(sdq_boot1, 10,11,"strength")

#testing between all pairs 
#Warning!!
#Requires a lot of RAM !!

plot(sdq_boot1, "edge", plot="difference", onlyNonZero=TRUE, order="sample")


plot(sdq_boot1)

#Simulation Studies


#Comparing Plots
Layout <- averageLayout(sdq_pcor,sdq_Glasso)
layout(t(1:2))
sdq_pcor <- qgraph(sdq_corMat, graph="sdq_pcor", layout="spring")
sdq_Glasso <- qgraph(sdq_corMat, graph="glasso", layout="spring", tuning=0.25, sampleSize=nrow(sdq_df2))

#stop split layout
layout(1)





#********************************************************************************************************#
#  Gender Wise  |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#  ===========  |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#  SDQ - Male   |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#********************************************************************************************************#

sdq_df1 <- shine_tsdq[,c(24,31,55,56:62)]
library(tidyverse)
sdq_df_m <- sdq_df1 %>% filter(demo_ch_gender=="Male") #Male
sdq_df_m <- sdq_df_m[,2:26]

#removing cases with missing values in data
sdq_df_m <- na.omit(sdq_df_m)
#renaming sdq_1 to sdq_01.....
names(sdq_df_m)[names(sdq_df_m) == "sdq_1"] <- "sdq_01"
names(sdq_df_m)[names(sdq_df_m) == "sdq_2"] <- "sdq_02"
names(sdq_df_m)[names(sdq_df_m) == "sdq_3"] <- "sdq_03"
names(sdq_df_m)[names(sdq_df_m) == "sdq_4"] <- "sdq_04"
names(sdq_df_m)[names(sdq_df_m) == "sdq_5"] <- "sdq_05"
names(sdq_df_m)[names(sdq_df_m) == "sdq_6"] <- "sdq_06"
names(sdq_df_m)[names(sdq_df_m) == "sdq_7"] <- "sdq_07"
names(sdq_df_m)[names(sdq_df_m) == "sdq_8"] <- "sdq_08"
names(sdq_df_m)[names(sdq_df_m) == "sdq_9"] <- "sdq_09"

#Creating Legend by defining groups of variables
sdq_gps <- scan("D:/Network Analysis/sdq_Groups.txt", what="character", sep="\n")

#Creating Item descriptions
sdq_items <- scan("D:/Network Analysis/sdq_Items.txt", what="character", sep="\n")

library(qgraph)
#partial correlation matrix
sdqm_corMat <- cor_auto(sdq_df_m)

#Plotting Graph with item description
sdqm_Glasso <- qgraph(sdqm_corMat, graph="glasso", layout="spring", tuning=0.25, 
                      sampleSize=nrow(sdq_df_m), groups=sdq_gps, nodeNames=sdq_items, 
                      legend.cex=0.20, vsize=5, esize=15, minimum=0.075)

#checking for redundant nodes
library(networktools)
goldbricker(sdq_df_m,p=0.05)

#Checking for highly influential nodes in network (Important Symptoms)
centralityPlot(sdqm_Glasso, include="all", decreasing=TRUE)

#Testing for Network accuracy and Centrality stability
library(bootnet)
#i.Network estimation
sdqm_net <- estimateNetwork(sdq_df_m, default="EBICglasso")
plot(sdqm_net, layout="spring", labels=TRUE)

#ii. Network accuracy (Using bootstrapping)
sdqm_boot1 <- bootnet(sdqm_net, nBoots=2500, nCores=8)
plot(sdqm_boot1, labels=FALSE, order="sample")

#iii. Network stability (Using bootstrapping)
sdqm_boot2 <- bootnet(sdqm_net, nBoots=2500, type="case", nCores=8)
plot(sdqm_boot2)
corStability(sdqm_boot2)


#********************************************************************************************************#
#  SDQ - Female   |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#********************************************************************************************************#

library(tidyverse)
sdq_df_f <- sdq_df1 %>% filter(demo_ch_gender=="Female") #FeMale
sdq_df_f <- sdq_df_f[,2:26]

#removing cases with missing values in data
sdq_df_f <- na.omit(sdq_df_f)
#renaming sdq_1 to sdq_01.....
names(sdq_df_f)[names(sdq_df_f) == "sdq_1"] <- "sdq_01"
names(sdq_df_f)[names(sdq_df_f) == "sdq_2"] <- "sdq_02"
names(sdq_df_f)[names(sdq_df_f) == "sdq_3"] <- "sdq_03"
names(sdq_df_f)[names(sdq_df_f) == "sdq_4"] <- "sdq_04"
names(sdq_df_f)[names(sdq_df_f) == "sdq_5"] <- "sdq_05"
names(sdq_df_f)[names(sdq_df_f) == "sdq_6"] <- "sdq_06"
names(sdq_df_f)[names(sdq_df_f) == "sdq_7"] <- "sdq_07"
names(sdq_df_f)[names(sdq_df_f) == "sdq_8"] <- "sdq_08"
names(sdq_df_f)[names(sdq_df_f) == "sdq_9"] <- "sdq_09"

#Creating Legend by defining groups of variables
sdq_gps <- scan("D:/Network Analysis/sdq_Groups.txt", what="character", sep="\n")

#Creating Item descriptions
sdq_items <- scan("D:/Network Analysis/sdq_Items.txt", what="character", sep="\n")

library(qgraph)
#partial correlation matrix
sdqf_corMat <- cor_auto(sdq_df_f)

#Plotting Graph with item description
sdqf_Glasso <- qgraph(sdqf_corMat, graph="glasso", layout="spring", tuning=0.25, 
                      sampleSize=nrow(sdq_df_f), groups=sdq_gps, nodeNames=sdq_items, 
                      legend.cex=0.20, vsize=5, esize=15, minimum=0.075)

#checking for redundant nodes
library(networktools)
goldbricker(sdq_df_f,p=0.05)

#Checking for highly influential nodes in network (Important Symptoms)
centralityPlot(sdqf_Glasso, include="all", decreasing=TRUE)

#Testing for Network accuracy and Centrality stability
library(bootnet)
#i.Network estimation
sdqf_net <- estimateNetwork(sdq_df_f, default="EBICglasso")
plot(sdqf_net, layout="spring", labels=TRUE)
centralityPlot(sdqf_net, include = c("Betweenness", "Closeness", "Strength"), 
               decreasing=TRUE)

#ii. Network accuracy (Using bootstrapping)
sdqf_boot1 <- bootnet(sdqf_net, nBoots=2500, nCores=8)
plot(sdqf_boot1, labels=FALSE, order="sample")

#iii. Network stability (Using bootstrapping)
sdqf_boot2 <- bootnet(sdqf_net, nBoots=2500, type="case", nCores=8)
plot(sdqf_boot2)
corStability(sdqf_boot2)

#iv. Testing for significance differences
differenceTest(sdqf_boot1, 22,33,"strength")


#*******************************************************#
#Comparing Plots    ||||||||||||||||||||||||||||||||||||#
#*******************************************************#
Layout <- averageLayout(sdqm_Glasso,sdqf_Glasso)
layout(t(1:2))
sdqm_Glasso <- qgraph(sdqm_corMat, graph="glasso", layout="spring", tuning=0.25,
                      sampleSize=nrow(sdq_df_m), groups=sdq_gps, nodeNames=sdq_items, 
                      legend.cex=0.20, vsize=5, esize=15, minimum=0.075, legend=FALSE)
sdqf_Glasso <- qgraph(sdqf_corMat, graph="glasso", layout="spring", tuning=0.25,
                      sampleSize=nrow(sdq_df_f), groups=sdq_gps, nodeNames=sdq_items, 
                      legend.cex=0.20, vsize=5, esize=15, minimum=0.075, legend=FALSE)

#stop split layout
layout(1)