
###*****************************************Microbiome analysis***********************************************************************

##Library loading
library(tidyverse)
library(rstatix)
library(microbiome)
library("readxl")       # necessary to import the data from Excel file
library(vegan)
library("phyloseq")
library("devtools")
#install_github("joey711")
library("ggplot2")
library("ggpubr")
library(rtk)
library(dplyr)
library(reshape2)
library(ape)
library(gridExtra)
library(plyr)


#To load extra-packages in r
scripts <- c("graphical_methods.R",
             "tree_methods.R",
             "plot_merged_trees.R",
             "specificity_methods.R",
             "ternary_plot.R",
             "richness.R",
             "edgePCA.R",
             "copy_number_correction.R",
             "import_frogs.R",
             "prevalence.R",
             "compute_niche.R")
urls <- paste0("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/R/", scripts)

for (url in urls) {
  source(url)
}

###Move to the suitable directory***********
setwd("~/Documents/marilyne/DataVisualization")

##data loading
TABLEOTU <- read_csv("combined_samples.OTU.filtered_Mahe.csv")#OTU table combined_samples.OTU.filtered_Mahe.csv with transformation of sample names: from OR.RT-1.. to RT1 for instance

NEW_TABLEOTU <- read.csv("MyOTUFilteredTab.csv", sep = "")#First column modification: To OTU 1...n to OTU OTU1...OTUn with command line sed -e 's/^/OTU/' Oracle_soils_16S_341F_785R_87_samples.OTU.filtered.cleaved.nosubstringOTUs.lulu.uchime.csv > MyOTUFilteredTab.csv

#The metadata with environnemental variables informations
ORACLE_METADATA<-read_csv("PratiquesCuturales.csv",na = c("", "NA"), quoted_na = TRUE,comment = "")#my metadata  

# create OTU_mat for phyloseq
otu_mat <- cbind(NEW_TABLEOTU[1],TABLEOTU[-1])

# create variable taxmat for phyloseq
new_cols <- c("Domain", "Phylum", "Class", "Ordre", "Family", "Genus", "Species")# New column headers
tax_mat <- NEW_TABLEOTU%>%
  separate(taxonomy, new_cols) %>%  # Split taxonomy column
  select(1,12:18)
tax_mat <-  tax_mat %>% drop_na("Domain","Phylum", "Class", "Ordre", "Family", "Genus", "Species")

#create metatadata for phyloseq
samples_df <- as.data.frame(ORACLE_METADATA)

# define the rownames from the OTU table
row.names(otu_mat) <- otu_mat$OTU

# remove the column otu since it is now used as a row name
otu_mat <- otu_mat %>% select (-OTU) 

#colnames(otu_mat)
# Idem for the two other dataframe

row.names(tax_mat) <- tax_mat$OTU
tax_mat <- tax_mat %>% select (-OTU) 

row.names(samples_df) <- samples_df$sample
samples_df <- samples_df %>% select (-sample) 

# Transformation into matrixes

otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

#Transform phyloseq objects
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df) 

# Define the variable coracle for phyloseq
oracle <- phyloseq(OTU, TAX, samples)
oracle

#data visualization
sample_names(oracle)
rank_names(oracle)
sample_variables(oracle)

#Data rarefaction
oracle<-rarefy_even_depth(oracle, rngseed = 1121983)
## check with sample_sums 
sample_sums(oracle)[1:5] 

#OTU Plylum barplot
plot_bar(oracle, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

##Microbome composition
#Select Bacteria(at Domain level) and aggregate by Phylum according to "Historiques_3ans"
p <- plot_composition(oracle, "Domain", "Bacteria", "Phylum", numberOfTaxa = 38, fill = "Phylum")
p <- p + facet_wrap(~Historique_3ans, scales = "free_x", nrow = 1)
plot(p)

#Select Bacteria(at Domain level) and aggregate by Phylum according to "Pratiques_culturales".
p <- plot_composition(oracle, "Domain", "Bacteria", "Phylum", numberOfTaxa = 38, fill = "Phylum")
p <- p + facet_wrap(~Pratiques_culturales, scales = "free_x", nrow = 1)
plot(p)

#Select Bacteria(at Domain level) and aggregate by Phylum according to "Type_sol_prod".
p <- plot_composition(oracle, "Domain", "Bacteria", "Phylum", numberOfTaxa = 10, fill = "Phylum")
p <- p + facet_wrap(~Type_sol_prod, scales = "free_x", nrow = 1)
plot(p)

#Select Bacteria(at Domain level) and aggregate by Phylum according to "Typo_Scientifique_sol".
p <- plot_composition(oracle, "Domain", "Bacteria", "Phylum", numberOfTaxa = 10, fill = "Phylum")
p <- p + facet_wrap(~Typo_Scientifique_sol, scales = "free_x", nrow = 1)
plot(p)


### Estimate richness: alpha diversity*****************************************************************************
alpha.diversity <-  estimate_richness(oracle, measures = c("Observed", "Shannon","Simpson", "InvSimpson"))
head(alpha.diversity)

## plot_richness : Historique_3ans
p <- plot_richness(oracle , color = "Historique_3ans", x = "Historique_3ans", measures = c("Observed", "Shannon", "Simpson", "InvSimpson"))
## plot as boxplot
p <- p + geom_boxplot(aes(fill = Historique_3ans), alpha=0.2)
plot(p)

## plot_richness : Pratiques_culturales
p <- plot_richness(oracle , color = "Pratiques_culturales", x = "Pratiques_culturales", measures = c("Observed", "Shannon", "Simpson", "InvSimpson"))
## plot as boxplot
p <- p + geom_boxplot(aes(fill = Pratiques_culturales), alpha=0.2)
plot(p)

## plot_richness : Type_sol_prod
p <- plot_richness(oracle , color = "Type_sol_prod", x = "Type_sol_prod", measures = c("Observed", "Shannon", "Simpson", "InvSimpson"))
## plot as boxplot
p <- p + geom_boxplot(aes(fill = Type_sol_prod), alpha=0.2)
plot(p)

## plot_richness : Typo_Scientifique_sol
p <- plot_richness(oracle , color = "Typo_Scientifique_sol", x = "Typo_Scientifique_sol", measures = c("Observed", "Shannon", "Simpson", "InvSimpson"))
## plot as boxplot
p <- p + geom_boxplot(aes(fill = Typo_Scientifique_sol), alpha=0.2)
plot(p)

#Plot taxa phylum according to "Historique_3ans" ******************************************************************
plot_bar(oracle, x="Phylum", fill = "Phylum")+
  #, facet_grid = "Historique_3ans") +
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack",inherit.aes = TRUE)+
  facet_wrap(~ Historique_3ans)+
  scale_y_continuous(trans = 'log10')

#Plot taxa phylum according to "Type_sol_prod" *******************************************************************
plot_bar(oracle, x="Phylum", fill = "Phylum")+
  #, facet_grid = "Type_sol_prod") 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack",inherit.aes = TRUE)+
  #scale_y_continuous(trans = 'log10')+
  facet_wrap(~ Type_sol_prod, scales = "free")

#Plot taxa phylum according to "Pratiques_culturales" ******************************************  
plot_bar(oracle, x="Phylum", fill = "Phylum") +
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack",inherit.aes = TRUE)+
  facet_wrap(~ Pratiques_culturales)+
  scale_y_continuous(trans = 'log10')


### Estimate beta diversity*****************************************************************************

#stressplot with vegan package
otu.oracle <- otu_table(oracle)
matrix <- otu.oracle
matrix_use<-as.matrix(matrix[,colSums(matrix)>=1])
library(vegan)
#NMDS <- metaMDS(matrix_use, distance = "bray", trymax = 100)

ord <- ordinate(oracle, method = "NMDS", distance = "bray")

#Plotting beta diversity on "Pratiques_culturales"
p <- plot_ordination(oracle, ord, color = "Pratiques_culturales")+
  geom_point(size= 4)
p <- p + theme_bw() + ggtitle("stress = 0.164                       NMDS + Bray-Curtis")
plot(p)

#Plotting beta diversity on "Historiques_3ans"
P <- plot_ordination(oracle, ord, color = "Historique_3ans")+
  geom_point(size=4)
P <- P + theme_bw() + ggtitle("stress = 0.164           NMDS + Bray-Curtis")
plot(P) 

#Plotting beta diversity on "Type_sol_prod"
P <- plot_ordination(oracle, ord, color = "Type_sol_prod")+
  geom_point(size=4)
P <- P + theme_bw() + ggtitle("stress = 0.164                     NMDS + Bray-Curtis") 
plot(P) 

#Plotting beta diversity on  "Typo_Scientifique_sol"
P <- plot_ordination(oracle, ord, color =  "Typo_Scientifique_sol")+
  geom_point(size=4)
P <- P + theme_bw() + ggtitle("stress = 0.164                     NMDS + Bray-Curtis") 
plot(P) 

##Phyla Barchart according to factor

#Plot taxa phylum according to "Historique_3ans" *********************************************************************************
plot_bar(oracle, x="Phylum", fill = "Phylum")+
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack",inherit.aes = TRUE)+
  facet_wrap(~ Historique_3ans,scale = "free")

#Plot taxa phylum according to "Type_sol_prod" *********************************************************************************
plot_bar(oracle, x="Phylum", fill = "Phylum")+
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack",inherit.aes = TRUE)+
  facet_wrap(~ Type_sol_prod,scale = "free")

#Plot taxa phylum according to "Typo_Scientifique_sol" ****************************************************************************
plot_bar(oracle, x="Phylum", fill = "Phylum")+
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack",inherit.aes = TRUE)+
  facet_wrap(~ Typo_Scientifique_sol,scale = "free")

#Plot taxa phylum according to "Pratiques_culturales" ******************************************  
plot_bar(oracle, x="Phylum", fill = "Phylum") +
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack",inherit.aes = TRUE)+
  facet_wrap(~ Pratiques_culturales,scales = "free")


#Corrélation entre variables quantitatives et qualitatives***********
## corrélation entre Historiques_3ans des précédents culturaux et les paramètres physico-chimiques
#Plotting pHeau
ggboxplot(ORACLE_METADATA, x = "Historique_3ans", y = "pHeau", 
          color = "Historique_3ans", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("Ass_continuous", "Rot_Cer", "Rot_Leg"),
          ylab = "pH", xlab = "Historiques_3ans")
ggsave("CorrelationPH&PrecCult.png", width = 16, height = 16, units = "cm")

#Plotting pHKCl
ggboxplot(ORACLE_METADATA, x = "Historique_3ans", y = "pHKCl", 
          color = "Historique_3ans", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("Ass_continuous", "Rot_Cer", "Rot_Leg"),
          ylab = "pHKCl", xlab = "Historiques_3ans")
ggsave("CorrelationPHKCl&PrecCult.png", width = 16, height = 16, units = "cm")

#Plotting C(g/kg)
ggboxplot(ORACLE_METADATA, x = "Historique_3ans", y = "C", 
          color = "Historique_3ans", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("Ass_continuous", "Rot_Cer", "Rot_Leg"),
          ylab = "C(g/kg)", xlab = "Historiques_3ans")
ggsave("CorrelationC&PrecCult.png", width = 16, height = 16, units = "cm")

#Plotting M.O(%)
ggboxplot(ORACLE_METADATA, x = "Historique_3ans", y = "MO", 
          color = "Historique_3ans", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("Ass_continuous", "Rot_Cer", "Rot_Leg"),
          ylab = "M.O(%)", xlab = "Historiques_3ans")
ggsave("CorrelationMO&PrecCult.png", width = 16, height = 16, units = "cm")

#Plotting Nt(%)
ggboxplot(ORACLE_METADATA, x = "Historique_3ans", y = "Nt", 
          color = "Historique_3ans", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("Ass_continuous", "Rot_Cer", "Rot_Leg"),
          ylab = "Nt(g/kg)", xlab = "Historiques_3ans")
ggsave("CorrelationNt&PrecCult.png", width = 16, height = 16, units = "cm")

#Plotting Kt(mg/kg)
ggboxplot(ORACLE_METADATA, x = "Historique_3ans", y = "Kt", 
          color = "Historique_3ans", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("Ass_continuous", "Rot_Cer", "Rot_Leg"),
          ylab = "Kt(mg/kg)", xlab = "Historiques_3ans")
ggsave("CorrelationKt&PrecCult.png", width = 16, height = 16, units = "cm")

#Plotting C/N
ggboxplot(ORACLE_METADATA, x = "Historique_3ans", y = "C/Npt", 
          color = "Historique_3ans", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("Ass_continuous", "Rot_Cer", "Rot_Leg"),
          ylab = "C/N", xlab = "Historiques_3ans")
ggsave("CorrelationCN&PrecCult.png", width = 16, height = 16, units = "cm")

#Plotting Pt
ggboxplot(ORACLE_METADATA, x = "Historique_3ans", y = "Pt", 
          color = "Historique_3ans", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("Ass_continuous", "Rot_Cer", "Rot_Leg"),
          ylab = "Pt(mg/kg)", xlab = "Historiques_3ans")
ggsave("CorrelationPt&PrecCult.png", width = 16, height = 16, units = "cm")

#Plotting CEC
ggboxplot(ORACLE_METADATA, x = "Historique_3ans", y = "CEC", 
          color = "Historique_3ans", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("Ass_continuous", "Rot_Cer", "Rot_Leg"),
          ylab = "CEC(cmolc/kg)", xlab = "Historiques_3ans")
ggsave("CorrelationCEC&PrecCult.png", width = 16, height = 16, units = "cm")

#Plotting Ca
ggboxplot(ORACLE_METADATA, x = "Historique_3ans", y = "Ca", 
          color = "Historique_3ans", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("Ass_continuous", "Rot_Cer", "Rot_Leg"),
          ylab = "Ca(cmol/kg)", xlab = "Historiques_3ans")
ggsave("CorrelationCa&PrecCult.png", width = 16, height = 16, units = "cm")


#Plotting Mg
ggboxplot(ORACLE_METADATA, x = "Historique_3ans", y = "Mg", 
          color = "Historique_3ans", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("Ass_continuous", "Rot_Cer", "Rot_Leg"),
          ylab = "Mg(cmolc/kg)", xlab = "Historiques_3ans")
ggsave("CorrelationMg&PrecCult.png", width = 16, height = 16, units = "cm")

#Plotting K
ggboxplot(ORACLE_METADATA, x = "Historique_3ans", y = "K", 
          color = "Historique_3ans", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("Ass_continuous", "Rot_Cer", "Rot_Leg"),
          ylab = "K(cmolc/kg)", xlab = "Historiques_3ans")
ggsave("CorrelationK&PrecCult.png", width = 16, height = 16, units = "cm")

#Plotting SBE
ggboxplot(ORACLE_METADATA, x = "Historique_3ans", y = "SBE", 
          color = "Historique_3ans", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("Ass_continuous", "Rot_Cer", "Rot_Leg"),
          ylab = "SBE(cmolc/kg)", xlab = "Historiques_3ans")
ggsave("CorrelationSBE&PrecCult.png", width = 16, height = 16, units = "cm")

#Plotting TS
ggboxplot(ORACLE_METADATA, x = "Historique_3ans", y = "TS", 
          color = "Historique_3ans", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("Ass_continuous", "Rot_Cer", "Rot_Leg"),
          ylab = "TS(%)", xlab = "Historiques_3ans")
ggsave("CorrelationTS&PrecCult.png", width = 16, height = 16, units = "cm")

#Plotting A%
ggboxplot(ORACLE_METADATA, x = "Historique_3ans", y = "A", 
          color = "Historique_3ans", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("Ass_continuous", "Rot_Cer", "Rot_Leg"),
          ylab = "A%", xlab = "Historiques_3ans")
ggsave("CorrelationA&PrecCult.png", width = 16, height = 16, units = "cm")

#Plotting L%
ggboxplot(ORACLE_METADATA, x = "Historique_3ans", y = "L", 
          color = "Historique_3ans", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("Ass_continuous", "Rot_Cer", "Rot_Leg"),
          ylab = "L(%)", xlab = "Historiques_3ans")
ggsave("CorrelationL&PrecCult.png", width = 16, height = 16, units = "cm")

#Plotting S%
ggboxplot(ORACLE_METADATA, x = "Historique_3ans", y = "S", 
          color = "Historique_3ans", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("Ass_continuous", "Rot_Cer", "Rot_Leg"),
          ylab = "S(%)", xlab = "Historiques_3ans")
ggsave("CorrelationS&PrecCult.png", width = 16, height = 16, units = "cm")

## corrélation entre Type_sol_prod des précédents culturaux et les paramètres physico-chimiques**********************
#"#B2182B" "#D6604D" "#F4A582" "#FDDBC7" "#D1E5F0" "#92C5DE" "#4393C3" "#2166AC"
#Plotting pHeau
ggboxplot(ORACLE_METADATA, x = "Type_sol_prod", y = "pHeau", 
          color = "Type_sol_prod", palette = c("#B2182B","#D6604D","#F4A582","#FDDBC7","#D1E5F0","#92C5DE","#4393C3"),
          order = c("Argileux", "Argilo-sableux", "Caillouteux","Gravillonnaire","Sableux","Sablo-argileux","Sablo-limoneux"),las =5,
          ylab = "pH", xlab = "Type_sol_prod")
ggsave("CorrelationPH&Type_sol_prod.png", width = 16, height = 16, units = "cm")

#Plotting pHKCl
ggboxplot(ORACLE_METADATA, x = "Type_sol_prod", y = "pHKCl", 
          color = "Type_sol_prod", palette = c("#B2182B","#D6604D","#F4A582","#FDDBC7","#D1E5F0","#92C5DE","#4393C3"),
          order = c("Argileux", "Argilo-sableux", "Caillouteux","Gravillonnaire","Sableux","Sablo-argileux","Sablo-limoneux"),
          ylab = "pHKCl", xlab = "Type_sol_prod")
ggsave("CorrelationPHKCl&Type_sol_prod.png", width = 16, height = 16, units = "cm")

#Plotting C(g/kg)
ggboxplot(ORACLE_METADATA, x = "Type_sol_prod", y = "C", 
          color = "Type_sol_prod", palette = c("#B2182B","#D6604D","#F4A582","#FDDBC7","#D1E5F0","#92C5DE","#4393C3"),
          order = c("Argileux", "Argilo-sableux", "Caillouteux","Gravillonnaire","Sableux","Sablo-argileux","Sablo-limoneux"),
          ylab = "C(g/kg)", xlab = "Type_sol_prod")
ggsave("CorrelationC&Type_sol_prod.png", width = 16, height = 16, units = "cm")

#Plotting M.O(%)
ggboxplot(ORACLE_METADATA, x = "Type_sol_prod", y = "MO", 
          color = "Type_sol_prod", palette = c("#B2182B","#D6604D","#F4A582","#FDDBC7","#D1E5F0","#92C5DE","#4393C3"),
          order = c("Argileux", "Argilo-sableux", "Caillouteux","Gravillonnaire","Sableux","Sablo-argileux","Sablo-limoneux"),
          ylab = "M.O(%)", xlab = "Type_sol_prod")
ggsave("CorrelationMO&Type_sol_prod.png", width = 16, height = 16, units = "cm")

#Plotting Nt(%)
ggboxplot(ORACLE_METADATA, x = "Type_sol_prod", y = "Nt", 
          color = "Type_sol_prod", palette = c("#B2182B","#D6604D","#F4A582","#FDDBC7","#D1E5F0","#92C5DE","#4393C3"),
          order = c("Argileux", "Argilo-sableux", "Caillouteux","Gravillonnaire","Sableux","Sablo-argileux","Sablo-limoneux"),
          ylab = "Nt(g/kg)", xlab = "Type_sol_prod")
ggsave("CorrelationNt&Type_sol_prod.png", width = 16, height = 16, units = "cm")

#Plotting Kt(mg/kg)
ggboxplot(ORACLE_METADATA, x = "Type_sol_prod", y = "Kt", 
          color = "Type_sol_prod", palette = c("#B2182B","#D6604D","#F4A582","#FDDBC7","#D1E5F0","#92C5DE","#4393C3"),
          order = c("Argileux", "Argilo-sableux", "Caillouteux","Gravillonnaire","Sableux","Sablo-argileux","Sablo-limoneux"),
          ylab = "Kt(mg/kg)", xlab = "Type_sol_prod")
ggsave("CorrelationKt&Type_sol_prod.png", width = 16, height = 16, units = "cm")

#Plotting C/N
ggboxplot(ORACLE_METADATA, x = "Type_sol_prod", y = "C/N", 
          color = "Type_sol_prod", palette = c("#B2182B","#D6604D","#F4A582","#FDDBC7","#D1E5F0","#92C5DE","#4393C3"),
          order = c("Argileux", "Argilo-sableux", "Caillouteux","Gravillonnaire","Sableux","Sablo-argileux","Sablo-limoneux"),
          ylab = "C/N", xlab = "Type_sol_prod")
ggsave("CorrelationCN&Type_sol_prod.png", width = 16, height = 16, units = "cm")

#Plotting Pt
ggboxplot(ORACLE_METADATA, x = "Type_sol_prod", y = "Pt", 
          color = "Type_sol_prod", palette = c("#B2182B","#D6604D","#F4A582","#FDDBC7","#D1E5F0","#92C5DE","#4393C3"),
          order = c("Argileux", "Argilo-sableux", "Caillouteux","Gravillonnaire","Sableux","Sablo-argileux","Sablo-limoneux"),
          ylab = "Pt(mg/kg)", xlab = "Type_sol_prod")
ggsave("CorrelationPt&Type_sol_prod.png", width = 16, height = 16, units = "cm")

#Plotting CEC
ggboxplot(ORACLE_METADATA, x = "Type_sol_prod", y = "CEC", 
          color = "Type_sol_prod", palette = c("#B2182B","#D6604D","#F4A582","#FDDBC7","#D1E5F0","#92C5DE","#4393C3"),
          order = c("Argileux", "Argilo-sableux", "Caillouteux","Gravillonnaire","Sableux","Sablo-argileux","Sablo-limoneux"),
          ylab = "CEC(cmolc/kg)", xlab = "Type_sol_prod")
ggsave("CorrelationCEC&Type_sol_prod.png", width = 16, height = 16, units = "cm")

#Plotting Ca
ggboxplot(ORACLE_METADATA, x = "Type_sol_prod", y = "Ca", 
          color = "Type_sol_prod", palette = c("#B2182B","#D6604D","#F4A582","#FDDBC7","#D1E5F0","#92C5DE","#4393C3"),
          order = c("Argileux", "Argilo-sableux", "Caillouteux","Gravillonnaire","Sableux","Sablo-argileux","Sablo-limoneux"),
          ylab = "Ca(cmol/kg)", xlab = "Type_sol_prod")
ggsave("CorrelationCa&Type_sol_prod.png", width = 16, height = 16, units = "cm")


#Plotting Mg
ggboxplot(ORACLE_METADATA, x = "Type_sol_prod", y = "Mg", 
          color = "Type_sol_prod", palette = c("#B2182B","#D6604D","#F4A582","#FDDBC7","#D1E5F0","#92C5DE","#4393C3"),
          order = c("Argileux", "Argilo-sableux", "Caillouteux","Gravillonnaire","Sableux","Sablo-argileux","Sablo-limoneux"),
          ylab = "Mg(cmolc/kg)", xlab = "Type_sol_prod")
ggsave("CorrelationMg&Type_sol_prod.png", width = 16, height = 16, units = "cm")

#Plotting K
ggboxplot(ORACLE_METADATA, x = "Type_sol_prod", y = "K", 
          color = "Type_sol_prod", palette = c("#B2182B","#D6604D","#F4A582","#FDDBC7","#D1E5F0","#92C5DE","#4393C3"),
          order = c("Argileux", "Argilo-sableux", "Caillouteux","Gravillonnaire","Sableux","Sablo-argileux","Sablo-limoneux"),
          ylab = "K(cmolc/kg)", xlab = "Type_sol_prod")
ggsave("CorrelationK&Type_sol_prod.png", width = 16, height = 16, units = "cm")

#Plotting SBE
ggboxplot(ORACLE_METADATA, x = "Type_sol_prod", y = "SBE", 
          color = "Type_sol_prod", palette = c("#B2182B","#D6604D","#F4A582","#FDDBC7","#D1E5F0","#92C5DE","#4393C3"),
          order = c("Argileux", "Argilo-sableux", "Caillouteux","Gravillonnaire","Sableux","Sablo-argileux","Sablo-limoneux"),
          ylab = "SBE(cmolc/kg)", xlab = "Type_sol_prod")
ggsave("CorrelationSBE&Type_sol_prod.png", width = 16, height = 16, units = "cm")

#Plotting TS
ggboxplot(ORACLE_METADATA, x = "Type_sol_prod", y = "TS", 
          color = "Type_sol_prod", palette = c("#B2182B","#D6604D","#F4A582","#FDDBC7","#D1E5F0","#92C5DE","#4393C3"),
          order = c("Argileux", "Argilo-sableux", "Caillouteux","Gravillonnaire","Sableux","Sablo-argileux","Sablo-limoneux"),
          ylab = "TS(%)", xlab = "Type_sol_prod")
ggsave("CorrelationTS&Type_sol_prod.png", width = 16, height = 16, units = "cm")

#Plotting A%
ggboxplot(ORACLE_METADATA, x = "Type_sol_prod", y = "A", 
          color = "Type_sol_prod", palette = c("#B2182B","#D6604D","#F4A582","#FDDBC7","#D1E5F0","#92C5DE","#4393C3"),
          order = c("Argileux", "Argilo-sableux", "Caillouteux","Gravillonnaire","Sableux","Sablo-argileux","Sablo-limoneux"),
          ylab = "A%", xlab = "Type_sol_prod")
ggsave("CorrelationA&Type_sol_prod.png", width = 16, height = 16, units = "cm")

#Plotting L%
ggboxplot(ORACLE_METADATA, x = "Type_sol_prod", y = "L", 
          color = "Type_sol_prod", palette = c("#B2182B","#D6604D","#F4A582","#FDDBC7","#D1E5F0","#92C5DE","#4393C3"),
          order = c("Argileux", "Argilo-sableux", "Caillouteux","Gravillonnaire","Sableux","Sablo-argileux","Sablo-limoneux"),
          ylab = "L(%)", xlab = "Type_sol_prod")
ggsave("CorrelationL&Type_sol_prod.png", width = 16, height = 16, units = "cm")

#Plotting S%
ggboxplot(ORACLE_METADATA, x = "Type_sol_prod", y = "S", 
          color = "Type_sol_prod", palette = c("#B2182B","#D6604D","#F4A582","#FDDBC7","#D1E5F0","#92C5DE","#4393C3"),
          order = c("Argileux", "Argilo-sableux", "Caillouteux","Gravillonnaire","Sableux","Sablo-argileux","Sablo-limoneux"),
          ylab = "S(%)", xlab = "Type_sol_prod")
ggsave("CorrelationS&Type_sol_prod.png", width = 16, height = 16, units = "cm")

## corrélation entre Pratiques_culturales et les paramètres physico-chimiques*********
#"#B2182B" "#D6604D" "#F4A582" "#FDDBC7" "#D1E5F0" "#92C5DE" "#4393C3" "#2166AC"
#Plotting pHeau
ggboxplot(ORACLE_METADATA, x = "Pratiques_culturales", y = "pHeau", 
          color = "Pratiques_culturales", palette = c("#009999", "#0000FF","#F4A582"),
          names.arg = c("CA","TM","Zaï"),
          #order = c("CA","TM","Zaï"),
          ylab = "pH", xlab = "Pratiques_culturales")
ggsave("CorrelationPH&Pratiques_culturales.png", width = 16, height = 16, units = "cm")

#Plotting pHKCl
ggboxplot(ORACLE_METADATA, x = "Pratiques_culturales", y = "pHKCl", 
          color = "Pratiques_culturales", palette = c("#009999", "#0000FF","#F4A582","#FDDBC7","#D1E5F0"),
          order = c("CA", "CA+Zaï", "TM","Zaï","Zaï+CA"),
          ylab = "pHKCl", xlab = "Pratiques_culturales")
ggsave("CorrelationPHKCl&Pratiques_culturales.png", width = 16, height = 16, units = "cm")

#Plotting C(g/kg)
ggboxplot(ORACLE_METADATA, x = "Pratiques_culturales", y = "C", 
          color = "Pratiques_culturales", palette = c("#009999", "#0000FF","#F4A582","#FDDBC7","#D1E5F0"),
          order = c("CA", "CA+Zaï", "TM","Zaï","Zaï+CA"),
          ylab = "C(g/kg)", xlab = "Pratiques_culturales")
ggsave("CorrelationC&Pratiques_culturales.png", width = 16, height = 16, units = "cm")

#Plotting M.O(%)
ggboxplot(ORACLE_METADATA, x = "Pratiques_culturales", y = "MO", 
          color = "Pratiques_culturales", palette = c("#009999", "#0000FF","#F4A582","#FDDBC7","#D1E5F0"),
          order = c("CA", "CA+Zaï", "TM","Zaï","Zaï+CA"),
          ylab = "M.O(%)", xlab = "Pratiques_culturales")
ggsave("CorrelationMO&Pratiques_culturales.png", width = 16, height = 16, units = "cm")

#Plotting Nt(%)
ggboxplot(ORACLE_METADATA, x = "Pratiques_culturales", y = "Nt", 
          color = "Pratiques_culturales", palette = c("#009999", "#0000FF","#F4A582","#FDDBC7","#D1E5F0"),
          order = c("CA", "CA+Zaï", "TM","Zaï","Zaï+CA"),
          ylab = "Nt(g/kg)", xlab = "Pratiques_culturales")
ggsave("CorrelationNt&Pratiques_culturales.png", width = 16, height = 16, units = "cm")

#Plotting Kt(mg/kg)
ggboxplot(ORACLE_METADATA, x = "Pratiques_culturales", y = "Kt", 
          color = "Pratiques_culturales", palette = c("#009999", "#0000FF","#F4A582","#FDDBC7","#D1E5F0"),
          order = c("CA", "CA+Zaï", "TM","Zaï","Zaï+CA"),
          ylab = "Kt(mg/kg)", xlab = "Pratiques_culturales")
ggsave("CorrelationKt&Pratiques_culturales.png", width = 16, height = 16, units = "cm")

#Plotting C/N
ggboxplot(ORACLE_METADATA, x = "Pratiques_culturales", y = "C/N", 
          color = "Pratiques_culturales", palette = c("#009999", "#0000FF","#F4A582","#FDDBC7","#D1E5F0"),
          order = c("CA", "CA+Zaï", "TM","Zaï","Zaï+CA"),
          ylab = "C/N", xlab = "Pratiques_culturales")
ggsave("CorrelationCN&Pratiques_culturales.png", width = 16, height = 16, units = "cm")

#Plotting Pt
ggboxplot(ORACLE_METADATA, x = "Pratiques_culturales", y = "Pt", 
          color = "Pratiques_culturales", palette = c("#009999", "#0000FF","#F4A582","#FDDBC7","#D1E5F0"),
          order = c("CA", "CA+Zaï", "TM","Zaï","Zaï+CA"),
          ylab = "Pt(mg/kg)", xlab = "Pratiques_culturales")
ggsave("CorrelationPt&Pratiques_culturales.png", width = 16, height = 16, units = "cm")

#Plotting CEC
ggboxplot(ORACLE_METADATA, x = "Pratiques_culturales", y = "CEC", 
          color = "Pratiques_culturales", palette = c("#009999", "#0000FF","#F4A582","#FDDBC7","#D1E5F0"),
          order = c("CA", "CA+Zaï", "TM","Zaï","Zaï+CA"),
          ylab = "CEC(cmolc/kg)", xlab = "Pratiques_culturales")
ggsave("CorrelationCEC&Pratiques_culturales.png", width = 16, height = 16, units = "cm")

#Plotting Ca
ggboxplot(ORACLE_METADATA, x = "Pratiques_culturales", y = "Ca", 
          color = "Pratiques_culturales", palette = c("#009999", "#0000FF","#F4A582","#FDDBC7","#D1E5F0"),
          order = c("CA", "CA+Zaï", "TM","Zaï","Zaï+CA"),
          ylab = "Ca(cmol/kg)", xlab = "Pratiques_culturales")
ggsave("CorrelationCa&Pratiques_culturales.png", width = 16, height = 16, units = "cm")


#Plotting Mg
ggboxplot(ORACLE_METADATA, x = "Pratiques_culturales", y = "Mg", 
          color = "Pratiques_culturales", palette = c("#009999", "#0000FF","#F4A582","#FDDBC7","#D1E5F0"),
          order = c("CA", "CA+Zaï", "TM","Zaï","Zaï+CA"),
          ylab = "Mg(cmolc/kg)", xlab = "Pratiques_culturales")
ggsave("CorrelationMg&Pratiques_culturales.png", width = 16, height = 16, units = "cm")

#Plotting K
ggboxplot(ORACLE_METADATA, x = "Pratiques_culturales", y = "K", 
          color = "Pratiques_culturales", palette = c("#009999", "#0000FF","#F4A582","#FDDBC7","#D1E5F0"),
          order = c("CA", "CA+Zaï", "TM","Zaï","Zaï+CA"),
          ylab = "K(cmolc/kg)", xlab = "Pratiques_culturales")
ggsave("CorrelationK&Pratiques_culturales.png", width = 16, height = 16, units = "cm")

#Plotting SBE
ggboxplot(ORACLE_METADATA, x = "Pratiques_culturales", y = "SBE", 
          color = "Pratiques_culturales", palette = c("#009999", "#0000FF","#F4A582","#FDDBC7","#D1E5F0"),
          order = c("CA", "CA+Zaï", "TM","Zaï","Zaï+CA"),
          ylab = "SBE(cmolc/kg)", xlab = "Pratiques_culturales")
ggsave("CorrelationSBE&Pratiques_culturales.png", width = 16, height = 16, units = "cm")

#Plotting TS
ggboxplot(ORACLE_METADATA, x = "Pratiques_culturales", y = "TS", 
          color = "Pratiques_culturales", palette = c("#009999", "#0000FF","#F4A582","#FDDBC7","#D1E5F0"),
          order = c("CA", "CA+Zaï", "TM","Zaï","Zaï+CA"),
          ylab = "TS(cmolc/kg)", xlab = "Pratiques_culturales")
ggsave("CorrelationTS&Pratiques_culturales.png", width = 16, height = 16, units = "cm")

#Plotting A%
ggboxplot(ORACLE_METADATA, x = "Pratiques_culturales", y = "A", 
          color = "Pratiques_culturales", palette = c("#009999", "#0000FF","#F4A582","#FDDBC7","#D1E5F0"),
          order = c("CA", "CA+Zaï", "TM","Zaï","Zaï+CA"),
          ylab = "A%", xlab = "Pratiques_culturales")
ggsave("CorrelationA&Pratiques_culturales.png", width = 16, height = 16, units = "cm")

#Plotting L%
ggboxplot(ORACLE_METADATA, x = "Pratiques_culturales", y = "L", 
          color = "Pratiques_culturales", palette = c("#009999", "#0000FF","#F4A582","#FDDBC7","#D1E5F0"),
          order = c("CA", "CA+Zaï", "TM","Zaï","Zaï+CA"),
          ylab = "L(%)", xlab = "Pratiques_culturales")
ggsave("CorrelationL&Pratiques_culturales.png", width = 16, height = 16, units = "cm")

#Plotting S%
ggboxplot(ORACLE_METADATA, x = "Pratiques_culturales", y = "S", 
          color = "Pratiques_culturales", palette = c("#009999", "#0000FF","#F4A582","#FDDBC7","#D1E5F0"),
          order = c("CA", "CA+Zaï", "TM","Zaï","Zaï+CA"),
          ylab = "S(%)", xlab = "Pratiques_culturales")
ggsave("CorrelationS&Pratiques_culturales.png", width = 16, height = 16, units = "cm")

