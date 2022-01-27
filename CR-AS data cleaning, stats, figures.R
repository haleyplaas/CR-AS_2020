#CR-AS data cleaning, stats, and figures for manuscript
rm(list=ls()) #clear environment
setwd("/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs") #set working directory 
#install all necessary packages
library(dplyr); library(tidyr); library(ggplot2); library(gam); library(gamm4); library(cowplot); library(mgcv); library(reshape2); library(MuMIn); library(stringr); library(ISLR); library(voxel); library(gridExtra); library(purrr); library(data.table); library(phyloseq); library(decontam); library(RAM); library(tidyverse); library(DESeq2); library(microbiome); library(vegan); library(viridis); library(patchwork); library(gapminder); library(tidyverse); library(ape); library(RColorBrewer)

#Reading in the dada2 datasets 
aerosol <- read.csv("Aerosol_16S_counts.csv", header = TRUE, na.strings = c(""," ", ".", "NA"))
names(aerosol)<-sapply(str_remove_all(colnames(aerosol),"AERO_"),"[")
water <- read.csv("Water_16S_counts.csv", header = TRUE, na.strings = c(""," ", ".", "NA"))  
names(water)<-sapply(str_remove_all(colnames(water),"W_"),"[")
ASV.w.taxonomy <- read.csv("/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/PhyloseqReadyFiles/Taxonomy_Table.csv", header=TRUE)
?lm()
#Removing the top ASVs found in the field blanks entirely from PM analysis:
# AEROSOL
all.ASVs <- aerosol %>% select(ASV)
FB.taxa.removed.from.matrix <- aerosol %>% filter(BLANK == 0)
#PM.with.all <- aerosol %>% anti_join(removal, by = "ASV") %>% mutate_if(is.numeric, ~replace_na(., 0)) %>% dplyr::select(-BLANK_1, -BLANK_2) #using average sequence depth
clean.ASVs <- FB.taxa.removed.from.matrix %>% select(ASV) 
clean.PM <- clean.ASVs %>% left_join(aerosol, by = "ASV") %>% mutate_if(is.numeric, ~replace_na(., 0)) %>% dplyr::select(-BLANK) #removing all altogether
PM.with.all <- all.ASVs %>% left_join(clean.PM, by = "ASV") %>% mutate_if(is.numeric, ~replace_na(., 0))
# PM.with.all <- all.ASVs %>% left_join(FB.removed, by = "ASV") %>% mutate_if(is.numeric, ~replace_na(., 0)) %>% dplyr::select(-BLANK_1, -BLANK_2) #seeing taxonomy
colnames(PM.with.all) <- c('ASV', 'QFF_Blank','S1_A', 'S2_A', 'S2_B', 'S3_A', 'S3_B', 'S4_A', 'S4_B', 'S5_A', 'S5_B', 'S6_B', 'S7_A', 'S7_B', 'S8_A', 'S8_B')
taxa.on.field.blank <- ASV.w.taxonomy %>% anti_join(FB.taxa.removed.from.matrix, by = "ASV")
write.csv(taxa.on.field.blank, "/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/files.for.manipulation/fieldblankQFF.contaminants.csv")

# PM phyloseq object 
# SEQUENCES
PM.with.all.mat <- PM.with.all %>% column_to_rownames("ASV") %>% data.matrix(PM.with.all)
sapply(PM.with.all.mat, as.numeric) # data matrix with sequence reads
seq.tab.PM <- otu_table(PM.with.all.mat, taxa_are_rows=TRUE)
# TAXONOMY
taxa_matrix <- read.csv("/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/PhyloseqReadyFiles/Taxonomy_Table.csv", row.names=1, header=TRUE)
taxa_matrix.1 <- as.matrix(taxa_matrix) # matrix assigning taxonomy to ASVs
taxa.tab <- tax_table(taxa_matrix.1) 
# METADATA
metadata <- read.csv("/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/PhyloseqReadyFiles/both.sample.types_metadata.csv", row.names = 1, header = TRUE)
metadata.PM <- metadata %>% filter(Sample_Type == "PM")
meta.tab.PM <- sample_data(metadata.PM)
# pre-decontam PM phyloseq object
PM.ps <- phyloseq(seq.tab.PM, meta.tab.PM, taxa.tab)
phylo.tree.PM <- rtree(ntaxa(PM.ps), rooted=TRUE, tip.label=taxa_names(PM.ps))
phyloseq.ob.PM <- phyloseq(seq.tab.PM, meta.tab.PM, taxa.tab, phylo.tree.PM)
# Removal of contaminants in PM extraction blank (QF/F) using R package "decontam" 
# visualization
decontam.df.PM <- as.data.frame(sample_data(phyloseq.ob.PM)) # Put sample_data into a ggplot-friendly data.frame
decontam.df.PM$LibrarySize <- sample_sums(phyloseq.ob.PM)
decontam.df.PM <- decontam.df.PM[order(decontam.df.PM$LibrarySize),]
decontam.df.PM$Index <- seq(nrow(decontam.df.PM))
ggplot(data=decontam.df.PM, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()
#FREQUENCY METHOD
contamdf.freq.PM <- isContaminant(phyloseq.ob.PM, method="frequency", conc="quant_reading")
head(contamdf.freq.PM)
table(contamdf.freq.PM$contaminant)
head(which(contamdf.freq.PM$contaminant))
set.seed(100)
plot_frequency(phyloseq.ob.PM, taxa_names(phyloseq.ob.PM)[sample(which(contamdf.freq.PM$contaminant),3)], conc="quant_reading") + xlab("DNA Concentration (ng/L)")
phyloseq.ob.no.contam.PM <- prune_taxa(!contamdf.freq.PM$contaminant, phyloseq.ob.PM)
phyloseq.ob.no.contam.PM
contaminants.on.QFF <- filter(contamdf.freq.PM, contaminant == TRUE) %>% rownames_to_column("ASV")
decontam.qff.contaminants <- contaminants.on.QFF %>% left_join(ASV.w.taxonomy, by = "ASV")
QFF.ASVs.to.remove <- contaminants.on.QFF %>% select("ASV")
write.csv(decontam.qff.contaminants, "/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/files.for.manipulation/extraction.blank.QFF.contaminants.csv")

# WATER 
H2O.with.all <- water %>% mutate_if(is.numeric, ~replace_na(., 0))
# Water phyloseq object 
# SEQUENCES
H2O.with.all.mat <- H2O.with.all %>% column_to_rownames("ASV") %>% data.matrix(H2O.with.all)
sapply(H2O.with.all.mat, as.numeric) # data matrix with sequence reads
seq.tab.H2O <- otu_table(H2O.with.all.mat, taxa_are_rows=TRUE)
# METADATA
metadata.H2O <- metadata %>% filter(Sample_Type == "Water")
meta.tab.H2O <- sample_data(metadata.H2O)
# pre-decontam H2O phyloseq object
H2O.ps <- phyloseq(seq.tab.H2O, meta.tab.H2O, taxa.tab)
phylo.tree.H2O <- rtree(ntaxa(H2O.ps), rooted=TRUE, tip.label=taxa_names(H2O.ps))
phyloseq.ob.H2O <- phyloseq(seq.tab.H2O, meta.tab.H2O, taxa.tab, phylo.tree.H2O)
# Removal of contaminants in H2O extraction blank (Supor) using R package "decontam" 
# visualization
decontam.df.H2O <- as.data.frame(sample_data(phyloseq.ob.H2O)) # Put sample_data into a ggplot-friendly data.frame
decontam.df.H2O$LibrarySize <- sample_sums(phyloseq.ob.H2O)
decontam.df.H2O <- decontam.df.H2O[order(decontam.df.H2O$LibrarySize),]
decontam.df.H2O$Index <- seq(nrow(decontam.df.H2O))
ggplot(data=decontam.df.H2O, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()
#FREQUENCY METHOD
contamdf.freq.H2O <- isContaminant(phyloseq.ob.H2O, method="frequency", conc="quant_reading")
head(contamdf.freq.H2O)
table(contamdf.freq.H2O$contaminant)
head(which(contamdf.freq.H2O$contaminant))
set.seed(100)
plot_frequency(phyloseq.ob.H2O, taxa_names(phyloseq.ob.H2O)[sample(which(contamdf.freq.H2O$contaminant),3)], conc="quant_reading") + xlab("DNA Concentration (ng/L)")
phyloseq.ob.no.contam.H2O <- prune_taxa(!contamdf.freq.H2O$contaminant, phyloseq.ob.H2O)
phyloseq.ob.no.contam.H2O
contaminants.on.Supor <- filter(contamdf.freq.H2O, contaminant == TRUE) %>% rownames_to_column("ASV") 
decontam.Supor <- contaminants.on.Supor %>% left_join(ASV.w.taxonomy, by = "ASV")
Supor.ASVs.to.remove <- decontam.Supor %>% select("ASV")
write.csv(decontam.Supor, "/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/files.for.manipulation/extractionSupor.contaminants.csv")

# COMBINED SAMPLE TYPES 
H2O.with.all.1 <- H2O.with.all %>% select(ASV:`10.01_B`, -Supor_Blank)
H2O.with.all.2 <- H2O.with.all.1 %>% anti_join(Supor.ASVs.to.remove, by = "ASV")
PM.with.all.1 <- PM.with.all %>% select(ASV:`S8_B`, -QFF_Blank)
PM.with.all.2 <- PM.with.all.1 %>% anti_join(QFF.ASVs.to.remove, by = "ASV")
both.sample.types <- H2O.with.all.2 %>% left_join(PM.with.all.2, by = c("ASV"), keep = FALSE)
seq.taxa.combined <- both.sample.types %>% left_join(ASV.w.taxonomy, by = "ASV", keep = FALSE) %>% column_to_rownames("ASV") 

# ASV data from PM and H2O into combined phyloseq object 
# SEQUENCES
both.sample.types.1 <- both.sample.types %>% dplyr::relocate('ASV','S1_A', 'S2_A', 'S2_B', 'S3_A', 'S3_B', 'S4_A', 'S4_B', 'S5_A', 'S5_B', 'S6_B', 'S7_A', 'S7_B', 'S8_A', 'S8_B', "06.11_A", "06.11_B", "06.23_A", "06.23_B", "07.07_A", "07.07_B", "07.21_A", "07.21_B", "08.12_A", "08.12_B", "08.18_A", "08.18_B", "09.01_A", "09.01_B","09.15_A", "09.15_B", "10.01_A", "10.01_B")
both.sample.types.2 <- all.ASVs %>% left_join(both.sample.types.1, by = "ASV") %>% mutate_if(is.numeric, ~replace_na(., 0)) %>% column_to_rownames("ASV")
both.sample.types.2 <- data.matrix(both.sample.types.2) 
sapply(both.sample.types.2, as.numeric) # data matrix with sequence reads
seq.tab <- otu_table(both.sample.types.2, taxa_are_rows=TRUE) 
taxa_names(seq.tab)
# METADATA
metadata.both <- metadata %>% filter(Sample_or_Control == "SAMPLE")
meta.tab <- sample_data(metadata.both)
# Final Phyloseq object
both.sample.types.ps <- phyloseq(seq.tab, meta.tab, taxa.tab)
phylo.tree <- rtree(ntaxa(both.sample.types.ps), rooted=TRUE, tip.label=taxa_names(both.sample.types.ps))
phyloseq.ob <- phyloseq(seq.tab, meta.tab, taxa.tab, phylo.tree)

#Showing sequencing depths
sequence.depths <- data.table(as(sample_data(phyloseq.ob), "data.frame"),
                      TotalReads = sample_sums(phyloseq.ob), keep.rownames = TRUE)
setnames(sequence.depths, "rn", "SampleID")
ggplot(sequence.depths, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth") + facet_wrap("Sample_Type")

#Cumulative Sum Scaling (to normalize based on variable sequencing depths between PM and Water samples)
CSS <- phyloseq::phyloseq_to_metagenomeSeq(phyloseq.ob)
CSS.1 <- metagenomeSeq::MRcounts(CSS, norm = TRUE, log = TRUE)
CSS.2 <- as.data.frame(CSS.1)
#removing ASVs with no reads from df
CSS.3 <- CSS.2 %>% rownames_to_column("ASV")
RA.df <- CSS.3 %>% left_join(ASV.w.taxonomy, by = "ASV")

RA.df <- seq.taxa.combined.1 %>% dplyr::relocate('S1_A', 'S2_A', 'S2_B', 'S3_A', 'S3_B', 'S4_A', 'S4_B', 'S5_A', 'S5_B', 'S6_B', 'S7_A', 'S7_B', 'S8_A', 'S8_B', "06.11_A", "06.11_B", "06.23_A", "06.23_B", "07.07_A", "07.07_B", "07.21_A", "07.21_B", "08.12_A", "08.12_B", "08.18_A", "08.18_B", "09.01_A", "09.01_B","09.15_A", "09.15_B", "10.01_A", "10.01_B") %>% mutate_if(is.numeric, ~replace_na(., 0)) 

#Estimating Water values during aerosol sampling periods rather than at interval points
write.csv(RA.df, "/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/files.for.manipulation/to.calculate.water.during.sampling.periods.csv")
water.sampling.periods <- read.csv("/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/water.sampling.period.csv", header=TRUE)
RA.df.0 <- RA.df %>% left_join(water.sampling.periods, by = "ASV", keep = FALSE) 
RA.df.SPs <- RA.df.0 %>% dplyr::select("ASV", 'S1_A', 'S2_A', 'S2_B', 'S3_A', 'S3_B', 'S4_A', 'S4_B', 'S5_A', 'S5_B', 'S6_B', 'S7_A', 'S7_B', 'S8_A', 'S8_B', "W1_A", "W2_A" , "W2_B" , "W3_A" , "W3_B" , "W4_A" , "W4_B" , "W5_A" , "W5_B",  "W6_B",  "W7_A"  ,"W7_B"  ,"W8_A",  "W8_B", "Kingdom", "Phylum"  ,"Class" ,  "Order" ,  "Family" , "Genus"  , "Species", "NA." )

#Plot Color Schemes
brewer.pal(n = 10, name = 'Paired')

class.colors <- c("Acidobacteria" = "#A6CEE3",
           "Actinobacteria" = "#1F78B4",
           "Bacteroidetes"   = "#B2DF8A",
           "Cyanobacteria"   = "#33A02C",
           "Firmicutes"  = "#FB9A99",
           "Planctomycetes" = "#E31A1C",
           "Proteobacteria" = "#FDBF6F",
           "z.Other" = "Gray")

genus.colors <- c("Caenarcaniphilales_X" = "#A6CEE3",
                    "Chroococcidiopsaceae" = "#1F78B4",
                    "Cyanobiaceae"   = "#B2DF8A",
                    "Microcystaceae"  = "#33A02C",
                    "Nostocaceae" = "#FB9A99",
                    "Obscuribacterales_X" = "#E31A1C",
                    "Phormidiaceae" = "#FDBF6F",
                    "Synechococcaceae" = "#FF7F00",
                    "z.Not.Assigned" = "#CAB2D6",
                    "z.Other" = "Gray")

species.colors <- c("Aphanizomenon_NIES81" = "indianred1", 
                    "Caenarcaniphilales_XX" = "#A6CEE3",
                    "Chroococcidiopsis_SAG_2023" = "#1F78B4",
                    "Cyanobium_PCC6307" = "#B2DF8A",
                    "Dolichospermum_NIES41" = "#FB9A99",
                    "Microcystis_PCC7914" = "#33A02C",
                    "Obscuribacterales_XX" = "#E31A1C",
                    "Synechococcus_PCC7942" = "#FF7F00",
                    "Tolypothrix_PCC7601" = "#6A3D9A",
                    "z.Not.Assigned" = "#CAB2D6",
                    "z.Other" = "Gray")
  
#plotting relative abundance by CSS, with water values averaged for aerosol sampling periods
RA.df.1 <- RA.df.SPs %>% group_by(Class) %>% summarise_if(is.numeric, ~sum(na.exclude(.)))
RA.df.2 <- RA.df.1 %>% dplyr::select(S1_A:`W8_B`)
RA.df.3 <- RA.df.2 %>% summarise(across(S1_A:`W8_B`, sum))
RA.df.4 <- mapply(`/`, RA.df.2, RA.df.3)
Class.26 <- RA.df.1 %>% dplyr::select("Class") %>% unique()
RA.df.5 <- cbind(Class.26, RA.df.4)
RA.df.6 <- RA.df.5 %>% mutate_at(vars(Class), ~replace_na(., "z.Not.Assigned")) 
write.csv(RA.df.6, "/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/files.for.manipulation/Relative.Abundance.Class.csv")
RA.df.7 <- RA.df.6 %>% mutate(Class = case_when(
  Class == "Acidobacteria" ~ "Acidobacteria",
  Class == "Actinobacteria" ~ "Actinobacteria",
  Class ==  "Bacteroidetes" ~ "Bacteroidetes", 
  Class == "Cyanobacteria" ~ "Cyanobacteria",
  Class ==  "Firmicutes" ~ "Firmicutes", 
  Class ==  "Planctomycetes" ~ "Planctomycetes", 
  Class ==  "Proteobacteria" ~	"Proteobacteria", 
  Class ==  "Apicomplexa:apic" ~ "z.Other",
  Class ==  "Armatimonadetes" ~ "z.Other",
  Class ==  "Chlamydiae"     ~ "z.Other",
  Class ==   "Chloroflexi"    ~ "z.Other",
  Class ==   "Chlorophyta:plas"  ~ "z.Other",
  Class ==   "Cryptophyta:plas"  ~ "z.Other",
  Class ==   "DeinococcusThermus" ~ "z.Other",
  Class ==    "Dependentiae"       ~ "z.Other",
  Class ==    "Dinoflagellata:plas" ~ "z.Other",
  Class ==  "Elusimicrobia" ~ "z.Other",
  Class ==   "Entotheonellaeota"  ~ "z.Other",
  Class ==   "Epsilonbacteraeota" ~ "z.Other",
  Class ==   "Fibrobacteres"   ~ "z.Other",
  Class ==   "Fusobacteria"   ~ "z.Other",
  Class ==  "Gemmatimonadetes"  ~ "z.Other",
  Class ==  "Hydrogenedentes"   ~ "z.Other",
  Class ==   "Kiritimatiellaeota"  ~ "z.Other",
  Class ==   "Nitrospirae"   ~ "z.Other",
  Class ==   "Ochrophyta:plas" ~ "z.Other",
  Class ==   "Patescibacteria"  ~ "z.Other",
  Class ==   "Spirochaetes"   ~ "z.Other",
  Class ==   "Streptophyta:mito" ~ "z.Other",
  Class ==  "Streptophyta:plas" ~ "z.Other",
  Class ==   "Synergistetes"    ~ "z.Other",
  Class == "Verrucomicrobia"  ~ "z.Other",
  Class == "z.Not.Assigned"  ~ "z.Other")) 
RA.df.8 <- RA.df.7 %>% group_by(Class) %>% summarise(across(S1_A:`W8_B`, sum))
RA.df.9 <- pivot_longer(RA.df.8, cols = c('S1_A', 'S2_A', 'S2_B', 'S3_A', 'S3_B', 'S4_A', 'S4_B', 'S5_A', 'S5_B', 'S6_B', 'S7_A', 'S7_B', 'S8_A', 'S8_B', "W1_A", "W2_A" , "W2_B" , "W3_A" , "W3_B" , "W4_A" , "W4_B" , "W5_A" , "W5_B",  "W6_B",  "W7_A"  ,"W7_B"  ,"W8_A",  "W8_B"), names_to = "Sample") %>% mutate(Sample_Type = case_when(str_detect(Sample, "S") ~ "PM", str_detect(Sample, "W") ~ "Water")) %>% mutate(Site = case_when(str_detect(Sample, "_A") ~ "A", str_detect(Sample, "_B") ~ "B")) 
RA.df.PM <- RA.df.9 %>% dplyr::filter(Sample_Type == "PM")
RA.df.Water <- RA.df.9 %>% dplyr::filter(Sample_Type == "Water")
RA.plot.PM <- ggplot(RA.df.PM, aes(fill = Class, y=value, x=Sample)) + 
  geom_bar(position="stack", stat="identity", colour = "black") + 
  theme_half_open() + 
  theme(text = element_text(size=10), 
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(3, "mm"),
        legend.position = "none") + 
  scale_fill_manual(values=class.colors) + 
  labs(title = expression(paste(PM[2.5]), element_text(face = "bold")), y = "Relative Abundance", x = "Sample") + 
  guides(fill = guide_legend(ncol = 1)) 

RA.plot.Water <- ggplot(RA.df.Water, aes(fill = Class, y=value, x=Sample)) + 
  geom_bar(position="stack", stat="identity", colour = "black") + 
  theme_half_open() + 
  theme(text = element_text(size=10), 
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(3, "mm")) + 
  scale_fill_manual(values=class.colors) + 
  labs(title = "Water", y = "Relative Abundance", x = "Sample") + 
  guides(fill = guide_legend(ncol = 1))

RA.plot.PM + RA.plot.Water

# Resolved to Genus within the Class of Cyanobacteria
genus.RA.df.1 <- RA.df.SPs %>% filter(Class == "Cyanobacteria") %>% group_by(Genus) %>% summarise_if(is.numeric, ~sum(na.exclude(.)))
genus.RA.df.2 <- genus.RA.df.1 %>% dplyr::select(S1_A:`W8_B`)
genus.RA.df.3 <- genus.RA.df.2 %>% summarise(across(S1_A:`W8_B`, sum))
genus.RA.df.4 <- mapply(`/`, genus.RA.df.2, genus.RA.df.3)
Genus <- genus.RA.df.1 %>% dplyr::select("Genus") %>% unique()
genus.RA.df.5 <- cbind(Genus, genus.RA.df.4)
genus.RA.df.6 <- genus.RA.df.5 %>% mutate_at(vars(Genus), ~replace_na(., "z.Not.Assigned")) 
write.csv(genus.RA.df.6, "/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/files.for.manipulation/Relative.Abundance.Genus.csv")
genus.RA.df.7 <- genus.RA.df.6 %>% mutate(Genus = case_when(
  Genus == "Caenarcaniphilales_X" ~ "Caenarcaniphilales_X",
  Genus ==   "Chroococcales_X"    ~ "z.Other",
  Genus ==   "Coleofasciculaceae"  ~ "z.Other",
  Genus ==  "Chroococcidiopsaceae" ~ "Chroococcidiopsaceae", 
  Genus == "Cyanobiaceae" ~ "Cyanobiaceae",
  Genus ==   "Cyanobacteriaceae"  ~ "z.Other",
  Genus ==   "Gastranaerophilales_X" ~ "z.Other",
  Genus ==  "Microcystaceae" ~ "Microcystaceae", 
  Genus ==  "Nostocaceae" ~ "Nostocaceae", 
  Genus ==  "Obscuribacterales_X" ~	"Obscuribacterales_X", 
  Genus ==   "Oscillatoriales_X" ~ "z.Other",
  Genus ==  "Phormidiaceae" ~ "Phormidiaceae",
  Genus ==  "Pleurocapsales_X" ~ "z.Other",
  Genus ==  "Pseudanabaenaceae" ~ "z.Other",
  Genus ==  "Synechococcaceae" ~ "Synechococcaceae",
  Genus ==  "z.Not.Assigned" ~ "z.Not.Assigned")) %>% group_by(Genus) %>% summarise(across(S1_A:`W8_B`, sum))
genus.RA.df.8 <- pivot_longer(genus.RA.df.7, cols = c('S1_A', 'S2_A', 'S2_B', 'S3_A', 'S3_B', 'S4_A', 'S4_B', 'S5_A', 'S5_B', 'S6_B', 'S7_A', 'S7_B', 'S8_A', 'S8_B', "W1_A", "W2_A" , "W2_B" , "W3_A" , "W3_B" , "W4_A" , "W4_B" , "W5_A" , "W5_B",  "W6_B",  "W7_A"  ,"W7_B"  ,"W8_A",  "W8_B"), names_to = "Sample") %>% mutate(Sample_Type = case_when(str_detect(Sample, "S") ~ "PM", str_detect(Sample, "W") ~ "Water")) %>% mutate(Site = case_when(str_detect(Sample, "_A") ~ "A", str_detect(Sample, "_B") ~ "B")) 
genus.RA.df.PM <- genus.RA.df.8 %>% dplyr::filter(Sample_Type == "PM")
genus.RA.df.Water <- genus.RA.df.8 %>% dplyr::filter(Sample_Type == "Water")
genus.RA.plot.PM <- ggplot(genus.RA.df.PM, aes(fill = Genus, y=value, x=Sample)) + 
  geom_bar(position="stack", stat="identity", colour = "black") + 
  theme_half_open() + 
  theme(text = element_text(size=10), 
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(3, "mm"),
        legend.position = "none") + 
  scale_fill_manual(values=genus.colors) + 
  labs(title = expression(paste(PM[2.5]), element_text(face = "bold")), y = "Relative Abundance", x = "Sample") + 
  guides(fill = guide_legend(ncol = 1)) 

genus.RA.plot.Water <- ggplot(genus.RA.df.Water, aes(fill = Genus, y=value, x=Sample)) + 
  geom_bar(position="stack", stat="identity", colour = "black") + 
  theme_half_open() + 
  theme(text = element_text(size=10), 
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(3, "mm")) + 
  scale_fill_manual(values=genus.colors) + 
  labs(title = "Water", y = "Relative Abundance", x = "Sample") + 
  guides(fill = guide_legend(ncol = 1))

genus.RA.plot.PM + genus.RA.plot.Water

# Resolved to Species within the Class of Cyanobacteria
species.RA.df.1 <- RA.df.SPs %>% filter(Class == "Cyanobacteria") %>% group_by(Species) %>% summarise_if(is.numeric, ~sum(na.exclude(.)))
species.RA.df.2 <- species.RA.df.1 %>% dplyr::select(S1_A:`W8_B`)
species.RA.df.3 <- species.RA.df.2 %>% summarise(across(S1_A:`W8_B`, sum))
species.RA.df.4 <- mapply(`/`, species.RA.df.2, species.RA.df.3)
Species <- species.RA.df.1 %>% dplyr::select("Species") %>% unique()
species.RA.df.5 <- cbind(Species, species.RA.df.4)
species.RA.df.6 <- species.RA.df.5 %>% mutate_at(vars(Species), ~replace_na(., "z.Not.Assigned")) 
write.csv(species.RA.df.6, "/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/files.for.manipulation/Relative.Abundance.species.csv")
species.RA.df.7 <- species.RA.df.6 %>% mutate(Species = case_when(
  Species ==  "Anabaena_PCC7108" ~ "z.Other",
  Species == "Aphanizomenon_NIES81" ~ "Aphanizomenon_NIES81",
  Species == "Caenarcaniphilales_XX" ~ "Caenarcaniphilales_XX",
  Species == "Caldora_VP642b"~ "z.Other",
  Species == "Calothrix_PCC6303"~ "z.Other",
  Species == "Cephalothrix_SAG_75.79"~ "z.Other",    
  Species == "Chroococcidiopsis_PCC_7203"~ "z.Other",
  Species == "Chroococcidiopsis_SAG_2023"~ "Chroococcidiopsis_SAG_2023",
  Species == "Cyanobium_PCC6307" ~ "Cyanobium_PCC6307",
  Species == "Cyanothece"~ "z.Other",
  Species == "Dolichospermum_NIES41" ~ "Dolichospermum_NIES41",
  Species == "Gastranaerophilales_XX"~ "z.Other",
  Species == "Geminocystis_PCC6308"~ "z.Other",
  Species == "Gloeocapsa"~ "z.Other",
  Species == "Microcystis_PCC7914" ~ "Microcystis_PCC7914",
  Species == "Nodularia_PCC9350"~ "z.Other",
  Species == "Nostoc_PCC73102"~ "z.Other",
  Species == "Nostoc_PCC7524"~ "z.Other",
  Species == "Obscuribacterales_XX" ~ "Obscuribacterales_XX",
  Species == "Phormidium"~ "z.Other",
  Species == "Phormidium_IAM_M71"~ "z.Other",        
  Species == "Pleurocapsa"~ "z.Other",
  Species == "Pseudanabaena_PCC7429"~ "z.Other",
  Species == "Synechococcus_CC9902" ~ "z.Other",
  Species == "Synechococcus_PCC7942" ~ "Synechococcus_PCC7942",
  Species == "Tolypothrix_IAM_M259" ~ "z.Other",
  Species == "Tolypothrix_PCC7601" ~ "Tolypothrix_PCC7601",
  Species == "z.Not.Assigned" ~ "z.Not.Assigned")) %>% group_by(Species) %>% summarise(across(S1_A:`W8_B`, sum))
species.RA.df.8 <- pivot_longer(species.RA.df.7, cols = c('S1_A', 'S2_A', 'S2_B', 'S3_A', 'S3_B', 'S4_A', 'S4_B', 'S5_A', 'S5_B', 'S6_B', 'S7_A', 'S7_B', 'S8_A', 'S8_B', "W1_A", "W2_A" , "W2_B" , "W3_A" , "W3_B" , "W4_A" , "W4_B" , "W5_A" , "W5_B",  "W6_B",  "W7_A"  ,"W7_B"  ,"W8_A",  "W8_B"), names_to = "Sample") %>% mutate(Sample_Type = case_when(str_detect(Sample, "S") ~ "PM", str_detect(Sample, "W") ~ "Water")) %>% mutate(Site = case_when(str_detect(Sample, "_A") ~ "A", str_detect(Sample, "_B") ~ "B")) %>% mutate_if(is.numeric, ~replace_na(., 0))
species.RA.df.PM <- species.RA.df.8 %>% dplyr::filter(Sample_Type == "PM")
species.RA.df.Water <- species.RA.df.8 %>% dplyr::filter(Sample_Type == "Water")
species.RA.plot.PM <- ggplot(species.RA.df.PM, aes(fill = Species, y=value, x=Sample)) + 
  geom_bar(position="stack", stat="identity", colour = "black") + 
  theme_half_open() + 
  theme(text = element_text(size=10), 
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(3, "mm"),
        legend.position = "none") + 
  scale_fill_manual(values=species.colors) + 
  labs(title = expression(paste(PM[2.5]), element_text(face = "bold")), y = "Relative Abundance", x = "Sample") + 
  guides(fill = guide_legend(ncol = 1)) 

species.RA.plot.Water <- ggplot(species.RA.df.Water, aes(fill = Species, y=value, x=Sample)) + 
  geom_bar(position="stack", stat="identity", colour = "black") + 
  theme_half_open() + 
  theme(text = element_text(size=10), 
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(3, "mm")) + 
  scale_fill_manual(values=species.colors) + 
  labs(title = "Water", y = "Relative Abundance", x = "Sample") + 
  guides(fill = guide_legend(ncol = 1))

species.RA.plot.PM + species.RA.plot.Water

#examining Alpha and Beta Diversity 
alpha.diversity <- estimate_richness(phyloseq.ob)
alpha.diversity 

Observed <- plot_richness(phyloseq.ob,x = "Sample_Type" , color = "Sample_Type", measures = c("Observed")) + geom_boxplot() + theme_bw() + theme(text = element_text(size=7), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7, angle = 45, hjust = 1), legend.text = element_text(size = 6), legend.key.size = unit(2, "mm"), legend.position = "none") + labs(title = "Alpha Diversity")
Shannon <- plot_richness(phyloseq.ob,x = "Sample_Type" , color = "Sample_Type", measures = c("Shannon")) + geom_boxplot()+ theme_bw() + theme(text = element_text(size=7), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7, angle = 45, hjust = 1), legend.text = element_text(size = 6), legend.key.size = unit(2, "mm"), legend.position = "none") 
InvSimpson <- plot_richness(phyloseq.ob,x = "Sample_Type" , color = "Sample_Type", measures = c("InvSimpson")) + geom_boxplot()+ theme_bw() + theme(text = element_text(size=7), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7, angle = 45, hjust = 1), legend.text = element_text(size = 6), legend.key.size = unit(2, "mm"), legend.position = "none")

Observed + Shannon + InvSimpson

#T-test of Diversity Index results
alpha.diversity.1 <- rownames_to_column(alpha.diversity)
anova.run.1 <- pivot_longer(alpha.diversity.1, cols = c("Observed", "Chao1", "se.chao1", "ACE", "se.ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"), names_to = "Index") 
anova.run.2 <-anova.run.1 %>% mutate(Sample_Type = case_when(
  str_detect(anova.run.1$rowname, "X") ~ "Water", 
  str_detect(anova.run.1$rowname, "S") ~ "PM"))
Observed.test <- anova.run.2 %>% dplyr::filter(Index == "Observed")
t.Observed <- t.test(value ~ Sample_Type, data = Observed.test)
t.Observed
Shannon.test <- anova.run.2 %>% dplyr::filter(Index == "Shannon")
t.Shannon <- t.test(value ~ Sample_Type, data = Shannon.test)
t.Shannon
Invsimpson.test <- anova.run.2 %>% dplyr::filter(Index == "InvSimpson")
t.Invsimpson <- t.test(value ~ Sample_Type, data = Invsimpson.test)
t.Invsimpson

#beta diversity 
wunifrac_dist <- phyloseq::distance(phyloseq.ob, method = "unifrac", weighted = F)
ordination <- ordinate(phyloseq.ob, method="PCoA", distance=wunifrac_dist)
plot_ordination(phyloseq.ob, ordination, color="Sample_Type") + theme(aspect.ratio=1)

# Calculating aerosolization factors (AF) for all classes of bacteria, and genera and species of cyanobacteria
RA.df.AF <- pivot_longer(RA.df.6, cols = c('S1_A', 'S2_A', 'S2_B', 'S3_A', 'S3_B', 'S4_A', 'S4_B', 'S5_A', 'S5_B', 'S6_B', 'S7_A', 'S7_B', 'S8_A', 'S8_B', "W1_A", "W2_A" , "W2_B" , "W3_A" , "W3_B" , "W4_A" , "W4_B" , "W5_A" , "W5_B",  "W6_B",  "W7_A"  ,"W7_B"  ,"W8_A",  "W8_B"), names_to = "Sample") %>% mutate(Sample_Type = case_when(str_detect(Sample, "S") ~ "PM", str_detect(Sample, "W") ~ "Water")) %>% mutate(Site = case_when(str_detect(Sample, "_A") ~ "A", str_detect(Sample, "_B") ~ "B")) %>% mutate_if(is.numeric, ~replace_na(., 0)) %>% mutate(Sampling_Period = case_when(str_detect(Sample, "1") ~ "1", str_detect(Sample, "2") ~ "2", str_detect(Sample, "3") ~ "3", str_detect(Sample, "4") ~ "4", str_detect(Sample, "5") ~ "5", str_detect(Sample, "6") ~ "6", str_detect(Sample, "7") ~ "7", str_detect(Sample, "8") ~ "8")) %>% dplyr::select(-Sample) 
class.RAs <- RA.df.AF %>% group_by(Class, Sampling_Period, Site) %>% pivot_wider(names_from = Sample_Type, values_from = value)  
class.AFs <- class.RAs %>% mutate(`r-AF` = PM/Water) %>% mutate(`log.r-AF` = log(`r-AF`)) %>% mutate(`r-AF`= na_if(`r-AF`, "Inf")) %>% mutate(`r-AF`= na_if(`r-AF`, "NaN")) %>% unite(Sampling_Period, Sampling_Period, Site, remove= TRUE) %>% mutate(`log.r-AF` = log(`r-AF`)) %>% mutate_all(~replace(., is.na(.), 0)) %>% mutate_all(~replace(., is.infinite(.), 0)) 

genus.RA.df.AF <- pivot_longer(genus.RA.df.6, cols = c('S1_A', 'S2_A', 'S2_B', 'S3_A', 'S3_B', 'S4_A', 'S4_B', 'S5_A', 'S5_B', 'S6_B', 'S7_A', 'S7_B', 'S8_A', 'S8_B', "W1_A", "W2_A" , "W2_B" , "W3_A" , "W3_B" , "W4_A" , "W4_B" , "W5_A" , "W5_B",  "W6_B",  "W7_A"  ,"W7_B"  ,"W8_A",  "W8_B"), names_to = "Sample") %>% mutate(Sample_Type = case_when(str_detect(Sample, "S") ~ "PM", str_detect(Sample, "W") ~ "Water")) %>% mutate(Site = case_when(str_detect(Sample, "_A") ~ "A", str_detect(Sample, "_B") ~ "B")) %>% mutate_if(is.numeric, ~replace_na(., 0)) %>% mutate(Sampling_Period = case_when(str_detect(Sample, "1") ~ "1", str_detect(Sample, "2") ~ "2", str_detect(Sample, "3") ~ "3", str_detect(Sample, "4") ~ "4", str_detect(Sample, "5") ~ "5", str_detect(Sample, "6") ~ "6", str_detect(Sample, "7") ~ "7", str_detect(Sample, "8") ~ "8")) %>% dplyr::select(-Sample)
genus.RAs <- genus.RA.df.AF %>% group_by(Genus, Sampling_Period, Site) %>% pivot_wider(names_from = Sample_Type, values_from = value)  
genus.AFs <- genus.RAs %>% mutate(`r-AF` = PM/Water) %>%  mutate(`log.r-AF` = log(`r-AF`)) %>% mutate(`r-AF`= na_if(`r-AF`, "Inf")) %>% mutate(`r-AF`= na_if(`r-AF`, "NaN")) %>% unite(Sampling_Period, Sampling_Period, Site, remove= TRUE) %>% mutate(`log.r-AF` = log(`r-AF`)) %>% mutate_all(~replace(., is.na(.), 0)) %>% mutate_all(~replace(., is.infinite(.), 0)) 

species.RA.df.AF <- pivot_longer(species.RA.df.6, cols = c('S1_A', 'S2_A', 'S2_B', 'S3_A', 'S3_B', 'S4_A', 'S4_B', 'S5_A', 'S5_B', 'S6_B', 'S7_A', 'S7_B', 'S8_A', 'S8_B', "W1_A", "W2_A" , "W2_B" , "W3_A" , "W3_B" , "W4_A" , "W4_B" , "W5_A" , "W5_B",  "W6_B",  "W7_A"  ,"W7_B"  ,"W8_A",  "W8_B"), names_to = "Sample") %>% mutate(Sample_Type = case_when(str_detect(Sample, "S") ~ "PM", str_detect(Sample, "W") ~ "Water")) %>% mutate(Site = case_when(str_detect(Sample, "_A") ~ "A", str_detect(Sample, "_B") ~ "B")) %>% mutate_if(is.numeric, ~replace_na(., 0)) %>% mutate(Sampling_Period = case_when(str_detect(Sample, "1") ~ "1", str_detect(Sample, "2") ~ "2", str_detect(Sample, "3") ~ "3", str_detect(Sample, "4") ~ "4", str_detect(Sample, "5") ~ "5", str_detect(Sample, "6") ~ "6", str_detect(Sample, "7") ~ "7", str_detect(Sample, "8") ~ "8")) %>% dplyr::select(-Sample)
species.RAs <- species.RA.df.AF %>% group_by(Species, Sampling_Period, Site) %>% pivot_wider(names_from = Sample_Type, values_from = value)
species.AFs <- species.RAs %>% mutate(`r-AF` = PM/Water) %>% mutate(`log.r-AF` = log(`r-AF`)) %>% mutate(`r-AF`= na_if(`r-AF`, "Inf")) %>% mutate(`r-AF`= na_if(`r-AF`, "NaN")) %>% unite(Sampling_Period, Sampling_Period, Site, remove= TRUE) %>% mutate_all(~replace(., is.na(.), 0)) %>% mutate_all(~replace(., is.infinite(.), 0)) 

#heat map to show enrichment in aerosol
classes.for.heat.viz <- class.AFs %>% filter(Class == "Acidobacteria" | Class == "Actinobacteria" | Class == "Bacteroidetes" | Class == "Cyanobacteria" | Class == "Firmicutes" | Class == "Planctomycetes" | Class == "Proteobacteria")
ggplot(classes.for.heat.viz, aes(Sampling_Period, Class, fill= `log.r-AF`)) + geom_tile()

genera.for.heat.viz <- genus.AFs %>% filter(Genus == "Caenarcaniphilales_X" | Genus == "Chroococcidiopsaceae" | Genus == "Cyanobiaceae" | Genus == "Microcystaceae" | Genus == "Nostocaceae" | Genus == "Obscuribacterales_X" | Genus == "Phormidiaceae"| Genus == "Synechococcaceae")
ggplot(genera.for.heat.viz, aes(Sampling_Period, Genus, fill= `log.r-AF`)) + geom_tile()

species.for.heat.viz <- species.AFs %>% filter(Species == "Aphanizomenon_NIES81" | Species == "Caenarcaniphilales_XX" | Species == "Chroococcidiopsis_SAG_2023" | Species == "Cyanobium_PCC6307" | Species == "Dolichospermum_NIES41" | Species == "Microcystis_PCC7914" | Species == "Obscuribacterales_XX"| Species == "Synechococcus_PCC7942"| Species  == "Tolypothrix_PCC7601")
ggplot(species.for.heat.viz, aes(Sampling_Period, Species, fill= `log.r-AF`)) + geom_tile()  

# Univariate Regressions
# loading in the water quality metadata
water_metadata <- read.csv("/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/master_metadata.csv", header = TRUE, na.strings = c(""," ", ".", "NA")) %>% dplyr::select(-SAMPLE.ID) 
Sonde.Data <- read.csv("/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/Sonde_Data.csv", header = TRUE, na.strings = c(""," ", ".", "NA")) %>% mutate(Site = recode(Site, `1`="A", `2`="B"))
all.water.metadata <- cbind(water_metadata, Sonde.Data) %>% select(-c(22,23)) %>% unite("Site_Date", c(Site,Date), sep = "_") 
#Triplicates into averages with standard deviations 
avg.metadata <- all.water.metadata %>% group_by(Site_Date) %>% dplyr::summarise_all(~mean(.), na.rm= FALSE)
colnames(avg.metadata) <- paste(colnames(avg.metadata),"avg",sep="_") 
avg.metadata.1 <- avg.metadata %>% dplyr::rename("Site_Date" = "Site_Date_avg") %>% tidyr::separate(col = Site_Date, into = c("Site", "Date"), sep = "_", convert = TRUE) 
Sd.metadata <- all.water.metadata %>% group_by(Site_Date) %>% dplyr::summarise_all(~sd(.), na.rm=TRUE)
colnames(Sd.metadata) <- paste(colnames(Sd.metadata),"sd",sep="_")
Sd.metadata.1 <- Sd.metadata %>% dplyr::rename("Site_Date" = "Site_Date_sd") %>% tidyr::separate(col = Site_Date, into = c("Site", "Date"), sep = "_", convert = TRUE) 
metadata.with.sd <- avg.metadata.1 %>% dplyr::left_join(Sd.metadata.1, by = c("Site", "Date"), keep = FALSE)
write.csv(metadata.with.sd, "/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/all.water.metadata.csv") # this data frame was exported to manipulate externally in Excel. To approximate the water conditions during the course of each aerosol sampling, water metadata collected at the start and end of each aerosol sampling period was averaged. 
Sampling_Period_water_metadata.with.sd <- read.csv("/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/all.water.metadata.2.csv", header = TRUE, na.strings = c(""," ", ".", "NA","#DIV/0!"))
Sampling_Period_water_metadata <- Sampling_Period_water_metadata.with.sd %>% select(1:29)

# loading in the atmosphere metadata
met.station.data <- read.csv("/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/met.station.data.csv", header = TRUE, na.strings = c(""," ", ".", "NA", "NaN")) 
met.station.data.1 <- met.station.data  %>% mutate(Date = as.Date(Date, "%Y-%m-%d")) 
met.station.data.2 <- met.station.data.1 %>% group_by(Sampling_Period, Site) %>% dplyr::summarise_all((~mean(na.omit(.)))) %>% select(-Date)
met.station.data.2[is.na(met.station.data.2)] <- NA

# combining metadata
all.metadata <- Sampling_Period_water_metadata %>% dplyr::left_join(met.station.data.2, by = c("Site", "Sampling_Period"), keep = FALSE) 

# adding relative abundances to the metadata 
class.RAs$Sampling_Period <- sub("^", "S", class.RAs$Sampling_Period) #don't run twice
class.RAs.2 <- class.RAs %>% filter(Class == "Cyanobacteria") %>% pivot_wider(names_from = Class, values_from = c(PM, Water))
species.RAs$Sampling_Period <- sub("^", "S", species.RAs$Sampling_Period) #don't run twice
species.RAs.2 <- species.RAs %>% pivot_wider(names_from = Species, values_from = c(PM, Water)) %>% mutate(PM_toxin.producers = PM_Microcystis_PCC7914 + PM_Aphanizomenon_NIES81, Water_toxin.producers = Water_Microcystis_PCC7914 + Water_Aphanizomenon_NIES81)
relative.abundances <- class.RAs.2 %>% left_join(species.RAs.2, by = c("Site", "Sampling_Period"))
regression.df.0 <- all.metadata %>% dplyr::left_join(relative.abundances, by = c("Site", "Sampling_Period"))

# the linear regression function 
regression.df <- regression.df.0 %>% dplyr::select(-Site, -Sampling_Period) %>% mutate_if(is.integer, as.numeric)
colnames(regression.df) #need this for copying and pasting strings (column names) 
# Plot functions
linear.regression.plot <- function(outcome.var, predictor.var) {
  lm(outcome.var ~ predictor.var, data = regression.df)
  ggplot(regression.df, aes(x =  predictor.var, y = outcome.var)) + 
    geom_point() +
    stat_smooth(method = "lm", col = "black") + theme_bw() 
} 
# Plot with Site specifications
linear.regression.plot.by.site <- function(outcome.var, predictor.var)  { 
  ggplot(regression.df.0, aes(x =  predictor.var, y = outcome.var, color = Site)) + geom_point(aes(color= Site)) + stat_smooth(method = "lm", fill = NA) + theme_bw() + scale_color_manual(values = c("Black", "Gray"))
}
# specifically looking at how Wind Direction impacts each outcome variable by Site (important because onshore wind direction is variable between Site A and B)
wind.direction.plot <- function(outcome.var)  { 
  ggplot(regression.df.0, aes(x =  Wind.Direction, y = outcome.var, color = Site)) + geom_point(aes(color= Site)) + stat_smooth(method = "lm", fill = NA) + theme_bw() + scale_color_manual(values = c("Black", "Gray"))
}

#Outcome Variable: Cyanobacteria in PM
regression.df  <- regression.df %>% select(PM_Cyanobacteria, everything()) #reordering df for for loop
columns <- colnames(regression.df)
lm.test <- list()
for(i in 2:ncol(regression.df)){
  lm.test[[i]] <- lm(reformulate(columns[i], "PM_Cyanobacteria"), na.action=na.exclude, data = regression.df)
}
summaries.1 <- lapply(lm.test, broom::glance)
df.1 <- as.data.frame(do.call(rbind,summaries.1))
summaries.3 <- lapply(lm.test, summary)
df.3 <- as.data.frame(do.call(rbind,summaries.3))
as.data.frame(df.3)
df.4 <- df.3 %>% select("terms")
df.5 <- apply(df.4, 2, as.character) 
model <- df.5[-1,]
model <- data.frame(model) 
significant.predictors.cyanos <- cbind(model, df.1) %>% dplyr::filter(p.value < .1)
all.predictors.cyanos <- cbind(model, df.1)
#visualizations
cyano.1 <- linear.regression.plot(regression.df$PM_Cyanobacteria, regression.df$Diatom.CHLA_avg) + labs(x = "Diatom CHLA", y = "RA Cyanobacteria in PM")
cyano.2 <- linear.regression.plot(regression.df$PM_Cyanobacteria, regression.df$Cryptophyte.CHLA_avg) + labs(x = "Cryptophyte CHLA", y = "RA Cyanobacteria in PM")
cyano.3 <- linear.regression.plot(regression.df$PM_Cyanobacteria, regression.df$PM_Cyanobium_PCC6307) + labs(x = "RA Cyanobium in PM", y = "RA Cyanobacteria in PM")
cyano.4 <- linear.regression.plot(regression.df$PM_Cyanobacteria, regression.df$Water_Cyanobacteria) + labs(x = "RA Cyanobacteria in Water", y = "RA Cyanobacteria in PM") #see labels for any I need to print 
cyano.1 + cyano.2 + cyano.3 + cyano.4

#Outcome Variable: Cyanobacteria in water
regression.df  <- regression.df %>% select(Water_Cyanobacteria, everything()) #reordering df for for loop
columns <- colnames(regression.df)
lm.test <- list()
for(i in 2:ncol(regression.df)){
  lm.test[[i]] <- lm(reformulate(columns[i], "Water_Cyanobacteria"), na.action=na.exclude, data = regression.df)
}
summaries.1 <- lapply(lm.test, broom::glance)
df.1 <- as.data.frame(do.call(rbind,summaries.1))
summaries.3 <- lapply(lm.test, summary)
df.3 <- as.data.frame(do.call(rbind,summaries.3))
as.data.frame(df.3)
df.4 <- df.3 %>% select("terms")
df.5 <- apply(df.4, 2, as.character) 
model <- df.5[-1,]
model <- data.frame(model) 
significant.predictors.cyanos.water <- cbind(model, df.1) %>% dplyr::filter(p.value < .1)
all.predictors.cyanos.water <- cbind(model, df.1)
#visualizations
w.cyano.1 <- linear.regression.plot(regression.df$Water_Cyanobacteria, regression.df$NOx_avg) + labs(x = "NOx conc.", y = "RA Cyanobacteria in Water")
w.cyano.2 <- linear.regression.plot(regression.df$Water_Cyanobacteria, regression.df$CHLA_avg) + labs(x = "Diatom CHLA", y = "RA Cyanobacteria in Water")
w.cyano.3 <- linear.regression.plot(regression.df$Water_Cyanobacteria, regression.df$Salinity_avg) + labs(x = "Salinity", y = "RA Cyanobacteria in Water")
w.cyano.4 <- linear.regression.plot(regression.df$Water_Cyanobacteria, regression.df$Salinity_avg) + labs(x = "Salinity", y = "RA Cyanobacteria in Water")
w.cyano.5 <- linear.regression.plot(regression.df$Water_Cyanobacteria, regression.df$Temp_avg) + labs(x = "Water Temp.", y = "RA Cyanobacteria in Water")
w.cyano.6 <- linear.regression.plot(regression.df$Water_Cyanobacteria, regression.df$Water_Gloeocapsa) + labs(x = "RA Gloecapsa in H2O", y = "RA Cyanobacteria in Water")
w.cyano.1 + w.cyano.2 + w.cyano.3 + w.cyano.4 + w.cyano.5 + w.cyano.6

#Outcome Variable: Microcystis in PM
regression.df  <- regression.df %>% select(PM_Microcystis_PCC7914, everything()) #reordering df for for loop
columns <- colnames(regression.df)
lm.test <- list()
for(i in 2:ncol(regression.df)){
  lm.test[[i]] <- lm(reformulate(columns[i], "PM_Microcystis_PCC7914"), na.action=na.exclude, data = regression.df)
}
summaries.1 <- lapply(lm.test, broom::glance)
df.1 <- as.data.frame(do.call(rbind,summaries.1))
summaries.3 <- lapply(lm.test, summary)
df.3 <- as.data.frame(do.call(rbind,summaries.3))
as.data.frame(df.3)
df.4 <- df.3 %>% select("terms")
df.5 <- apply(df.4, 2, as.character) 
model <- df.5[-1,]
model <- data.frame(model) 
significant.predictors.PM_Microcystis_PCC7914 <- cbind(model, df.1) %>% dplyr::filter(p.value < .1)
all.predictors.PM_Microcystis_PCC7914 <- cbind(model, df.1)
#visualizations
micro.1 <- linear.regression.plot(regression.df$PM_Microcystis_PCC7914, regression.df$PM_Dolichospermum_NIES41) + labs(x = "RA Dolichospermum in PM", y = "RA Microcystis in PM")
micro.2 <- linear.regression.plot(regression.df$PM_Microcystis_PCC7914, regression.df$PM_Cyanobium_PCC6307) + labs(x = "RA Cyanobium in PM", y = "RA Microcystis in PM")
micro.3 <- linear.regression.plot(regression.df$PM_Microcystis_PCC7914, regression.df$PM_z.Not.Assigned)+ labs(x = "RA N.A. cyanobacterial species in PM", y = "RA Microcystis in PM")
micro.4 <- linear.regression.plot(regression.df$PM_Microcystis_PCC7914, regression.df$Water_Gloeocapsa) + labs(x = "RA Gloeocapsa in H2O", y = "RA Microcystis in PM")
micro.5 <- linear.regression.plot(regression.df$PM_Microcystis_PCC7914, regression.df$NH4_avg) + labs(x = "NH4 concentration", y = "RA Microcystis in PM")
micro.6 <- linear.regression.plot(regression.df$PM_Microcystis_PCC7914, regression.df$Water_Microcystis_PCC7914) + labs(x = "RA Microcystis in H2O", y = "RA Microcystis in PM")
micro.1 + micro.2 + micro.3 + micro.4 + micro.5 + micro.6

#Outcome Variable: Microcystis in water
regression.df  <- regression.df %>% select(Water_Microcystis_PCC7914, everything()) #reordering df for for loop
columns <- colnames(regression.df)
lm.test <- list()
for(i in 2:ncol(regression.df)){
  lm.test[[i]] <- lm(reformulate(columns[i], "Water_Microcystis_PCC7914"), na.action=na.exclude, data = regression.df)
}
summaries.1 <- lapply(lm.test, broom::glance)
df.1 <- as.data.frame(do.call(rbind,summaries.1))
summaries.3 <- lapply(lm.test, summary)
df.3 <- as.data.frame(do.call(rbind,summaries.3))
as.data.frame(df.3)
df.4 <- df.3 %>% select("terms")
df.5 <- apply(df.4, 2, as.character) 
model <- df.5[-1,]
model <- data.frame(model) 
significant.predictors.Water_Microcystis_PCC7914 <- cbind(model, df.1) %>% dplyr::filter(p.value < .1)
all.predictors.Water_Microcystis_PCC7914 <- cbind(model, df.1)
#visualizations
w.micro.1 <- linear.regression.plot(regression.df$Water_Microcystis_PCC7914, regression.df$NH4_avg) + labs(x = "NH4 conc.", y = "RA Microcystis in H2O")
w.micro.2 <- linear.regression.plot(regression.df$Water_Microcystis_PCC7914, regression.df$TDN_avg) + labs(x = "TDN conc.", y = "RA Microcystis in H2O")
w.micro.3 <- linear.regression.plot(regression.df$Water_Microcystis_PCC7914, regression.df$PM_Cyanobium_PCC6307) + labs(x = "RA Cyanobium in PM", y = "RA Microcystis in H2O")
w.micro.4 <- linear.regression.plot(regression.df$Water_Microcystis_PCC7914, regression.df$Water_Aphanizomenon_NIES81) + labs(x = "RA Aphanizomenon in H2O", y = "RA Microcystis in H2O")
w.micro.5 <- linear.regression.plot(regression.df$Water_Microcystis_PCC7914, regression.df$PM_Microcystis_PCC7914) + labs(x = "RA Microcystis in PM", y = "RA Microcystis in H2O")
w.micro.1 + w.micro.2 + w.micro.3 + w.micro.4 + w.micro.5 

#Outcome Variable: Dolichospermum_NIES41 in PM
regression.df  <- regression.df %>% select(PM_Dolichospermum_NIES41 , everything()) #reordering df for for loop
columns <- colnames(regression.df)
lm.test <- list()
for(i in 2:ncol(regression.df)){
  lm.test[[i]] <- lm(reformulate(columns[i], "PM_Dolichospermum_NIES41"), na.action=na.exclude, data = regression.df)
}
summaries.1 <- lapply(lm.test, broom::glance)
df.1 <- as.data.frame(do.call(rbind,summaries.1))
summaries.3 <- lapply(lm.test, summary)
df.3 <- as.data.frame(do.call(rbind,summaries.3))
as.data.frame(df.3)
df.4 <- df.3 %>% select("terms")
df.5 <- apply(df.4, 2, as.character) 
model <- df.5[-1,]
model <- data.frame(model) 
significant.predictors.PM_Dolichospermum_NIES41 <- cbind(model, df.1) %>% dplyr::filter(p.value < .10)
all.predictors.PM_Dolichospermum_NIES41 <- cbind(model, df.1)
#visualizations
dolicho.1 <- linear.regression.plot(regression.df$PM_Dolichospermum_NIES41, regression.df$PM_Microcystis_PCC7914) + labs(x = "RA Microcystis in PM", y = "RA Dolichospermum in PM")
dolicho.2 <- linear.regression.plot(regression.df$PM_Dolichospermum_NIES41, regression.df$C.to.N_avg) + labs(x = "C:N Molar Ratio", y = "RA Dolichospermum in PM")
dolicho.3 <- linear.regression.plot(regression.df$PM_Dolichospermum_NIES41, regression.df$Relative.Humidity) + labs(x = "Relative Humidity", y = "RA Dolichospermum in PM")
dolicho.4 <- linear.regression.plot(regression.df$PM_Dolichospermum_NIES41, regression.df$PM_Aphanizomenon_NIES81) + labs(x = "RA Aphanizomenon in PM", y = "RA Dolichospermum in PM")
dolicho.5 <- linear.regression.plot(regression.df$PM_Dolichospermum_NIES41, regression.df$Water_Dolichospermum_NIES41) + labs(x = "RA Dolichospermum in H2O", y = "RA Dolichospermum in PM")
dolicho.6 <- linear.regression.plot(regression.df$PM_Dolichospermum_NIES41, regression.df$Solar.Rad) + labs(x = "Solar Irradiance", y = "RA Dolichospermum in PM")
dolicho.7 <- linear.regression.plot(regression.df$PM_Dolichospermum_NIES41, regression.df$NH4_avg) + labs(x = "NH4", y = "RA Dolichospermum in PM")
dolicho.8 <- linear.regression.plot(regression.df$PM_Dolichospermum_NIES41, regression.df$PM_Cyanobium_PCC6307) + labs(x = "PM_Cyanobium", y = "RA Dolichospermum in PM")
dolicho.9 <- linear.regression.plot(regression.df$PM_Dolichospermum_NIES41, regression.df$PM_Chroococcidiopsis_SAG_2023) + labs(x = "RA Chroococcidiopsis_SAG_2023 in PM", y = "RA Dolichospermum in PM")
dolicho.10 <- linear.regression.plot(regression.df$PM_Dolichospermum_NIES41, regression.df$PM_Chroococcidiopsis_PCC_7203) + labs(x = "RA Chroococcidiopsis_PCC_7203 in PM", y = "RA Dolichospermum in PM")

dolicho.1 + dolicho.2 + dolicho.3 + dolicho.4 + dolicho.5 + dolicho.6 + dolicho.7 + dolicho.8 + dolicho.9 + dolicho.10 

#Looking at each function dolicho relative abundance in PM was highest when wind direction fell between 75-100 degrees (from the east) -- need to think up some sort of best visualization for this for manuscript, but might need to look at this by Site
dolicho.wind.direction <- linear.regression.plot(regression.df.0$PM_Dolichospermum_NIES41, regression.df$Wind.Direction)
wind.direction.plot(regression.df.0$PM_toxin.producers)

#Outcome Variable: Dolichospermum_NIES41 in Water
regression.df  <- regression.df %>% select(Water_Dolichospermum_NIES41 , everything()) #reordering df for for loop
columns <- colnames(regression.df)
lm.test <- list()
for(i in 2:ncol(regression.df)){
  lm.test[[i]] <- lm(reformulate(columns[i], "Water_Dolichospermum_NIES41"), na.action=na.exclude, data = regression.df)
}
summaries.1 <- lapply(lm.test, broom::glance)
df.1 <- as.data.frame(do.call(rbind,summaries.1))
summaries.3 <- lapply(lm.test, summary)
df.3 <- as.data.frame(do.call(rbind,summaries.3))
as.data.frame(df.3)
df.4 <- df.3 %>% select("terms")
df.5 <- apply(df.4, 2, as.character) 
model <- df.5[-1,]
model <- data.frame(model) 
significant.predictors.Water_Dolichospermum_NIES41 <- cbind(model, df.1) %>% dplyr::filter(p.value < .1)
all.predictors.Water_Dolichospermum_NIES41 <- cbind(model, df.1)
#visualizations
w.dolicho.1 <- linear.regression.plot(regression.df$Water_Dolichospermum_NIES41, regression.df$PM_Microcystis_PCC7914) + labs(x = "RA Microcystis in PM", y = "RA Dolichospermum in H2O")
w.dolicho.2 <- linear.regression.plot(regression.df$Water_Dolichospermum_NIES41, regression.df$Solar.Rad) + labs(x = "Solar Irradiance", y = "RA Dolichospermum in H2O")
w.dolicho.3 <- linear.regression.plot(regression.df$Water_Dolichospermum_NIES41, regression.df$Water_Cyanobium_PCC6307) + labs(x = "RA Cyanobium in H2O", y = "RA Dolichospermum in H2O")
w.dolicho.4 <- linear.regression.plot(regression.df$Water_Dolichospermum_NIES41, regression.df$Water_Caenarcaniphilales_XX) + labs(x = "RA Caenarcaphilales in H2O", y = "RA Dolichospermum in H2O")
w.dolicho.5 <- linear.regression.plot(regression.df$Water_Dolichospermum_NIES41, regression.df$NH4_avg) + labs(x = "NH4 conc.", y = "RA Dolichospermum in H2O")
w.dolicho.6 <- linear.regression.plot(regression.df$Water_Dolichospermum_NIES41, regression.df$C.to.N_avg) + labs(x = "C:N Molar Ratio", y = "RA Dolichospermum in H2O")
w.dolicho.7 <- linear.regression.plot(regression.df$Water_Dolichospermum_NIES41, regression.df$PO4_avg) + labs(x = "PO4 conc.", y = "RA Dolichospermum in H2O")
w.dolicho.8 <- linear.regression.plot(regression.df$Water_Dolichospermum_NIES41, regression.df$CHLA_avg) + labs(x = "CHLA conc.", y = "RA Dolichospermum in H2O")
w.dolicho.9 <- linear.regression.plot(regression.df$Water_Dolichospermum_NIES41, regression.df$PM.avg) + labs(x = "ambient PM conc.", y = "RA Dolichospermum in H2O")
w.dolicho.10 <- linear.regression.plot(regression.df$Water_Dolichospermum_NIES41, regression.df$Relative.Humidity) + labs(x = "RH", y = "RA Dolichospermum in H2O")
w.dolicho.11 <- linear.regression.plot(regression.df$Water_Dolichospermum_NIES41, regression.df$TDN_avg) + labs(x = "TDN conc.", y = "RA Dolichospermum in H2O")
w.dolicho.12 <- linear.regression.plot(regression.df$Water_Dolichospermum_NIES41, regression.df$PM_Cyanobium_PCC6307) + labs(x = "RA Cyanobium in PM", y = "RA Dolichospermum in H2O")
linear.regression.plot(regression.df$Water_Dolichospermum_NIES41, regression.df$Water_Aphanizomenon_NIES81) + labs(x = "aphani", y = "RA Dolichospermum in H2O")

w.dolicho.1 + w.dolicho.2 + w.dolicho.3 + w.dolicho.4 + w.dolicho.5 + w.dolicho.6 + w.dolicho.7 + w.dolicho.8 + w.dolicho.9 + w.dolicho.10 + w.dolicho.11 + w.dolicho.12 

#Outcome Variable: Aphanizomenon in PM
regression.df  <- regression.df %>% select(PM_Aphanizomenon_NIES81 , everything()) #reordering df for for loop
columns <- colnames(regression.df)
lm.test <- list()
for(i in 2:ncol(regression.df)){
  lm.test[[i]] <- lm(reformulate(columns[i], "PM_Aphanizomenon_NIES81"), na.action=na.exclude, data = regression.df)
}
summaries.1 <- lapply(lm.test, broom::glance)
df.1 <- as.data.frame(do.call(rbind,summaries.1))
summaries.3 <- lapply(lm.test, summary)
df.3 <- as.data.frame(do.call(rbind,summaries.3))
as.data.frame(df.3)
df.4 <- df.3 %>% select("terms")
df.5 <- apply(df.4, 2, as.character) 
model <- df.5[-1,]
model <- data.frame(model) 
significant.predictors.PM_Aphanizomenon_NIES81 <- cbind(model, df.1) %>% dplyr::filter(p.value < .1)
all.predictors.PM_Aphanizomenon_NIES81 <- cbind(model, df.1)
#visualizations
Aphan.1 <- linear.regression.plot(regression.df$PM_Aphanizomenon_NIES81, regression.df$PM_Dolichospermum_NIES41) + labs(x = "RA Dolichospermum in PM", y = "RA Aphanizomenon in PM")
Aphan.2 <- linear.regression.plot(regression.df$PM_Aphanizomenon_NIES81, regression.df$Solar.Rad) + labs(x = "Solar Irradiance", y = "RA Aphanizomenon in PM")
Aphan.3 <- linear.regression.plot(regression.df$PM_Aphanizomenon_NIES81, regression.df$Water_Dolichospermum_NIES41) + labs(x = "RA Dolichospermum in H2O", y = "RA Aphanizomenon in PM")
Aphan.4 <- linear.regression.plot(regression.df$PM_Aphanizomenon_NIES81, regression.df$Water_Caenarcaniphilales_XX) + labs(x = "RA Caenarcaniphilales in H2O", y = "RA Aphanizomenon in PM")
Aphan.5 <- linear.regression.plot(regression.df$PM_Aphanizomenon_NIES81, regression.df$NH4_avg) + labs(x = "NH4", y = "RA Aphanizomenon in PM")
Aphan.6 <- linear.regression.plot(regression.df$PM_Aphanizomenon_NIES81, regression.df$Relative.Humidity) + labs(x = "Relative Humidity", y = "RA Aphanizomenon in PM")
Aphan.7 <- linear.regression.plot(regression.df$PM_Aphanizomenon_NIES81, regression.df$TDN_avg) + labs(x = "TDN", y = "RA Aphanizomenon in PM")
Aphan.8 <- linear.regression.plot(regression.df$PM_Aphanizomenon_NIES81, regression.df$PM_Cyanobium_PCC6307) + labs(x = "RA Cyanobium in PM", y = "RA Aphanizomenon in PM")
Aphan.9 <- linear.regression.plot(regression.df$PM_Aphanizomenon_NIES81, regression.df$Water_Aphanizomenon_NIES81) + labs(x = "RA Aphanizomenon in H2O", y = "RA Aphanizomenon in PM")

Aphan.1 + Aphan.2 + Aphan.3 + Aphan.4 + Aphan.5 + Aphan.6 + Aphan.7 + Aphan.8 + Aphan.9

#Outcome Variable: Aphanizomenon in Water
regression.df  <- regression.df %>% select(Water_Aphanizomenon_NIES81 , everything()) #reordering df for for loop
columns <- colnames(regression.df)
lm.test <- list()
for(i in 2:ncol(regression.df)){
  lm.test[[i]] <- lm(reformulate(columns[i], "Water_Aphanizomenon_NIES81"), na.action=na.exclude, data = regression.df)
}
summaries.1 <- lapply(lm.test, broom::glance)
df.1 <- as.data.frame(do.call(rbind,summaries.1))
summaries.3 <- lapply(lm.test, summary)
df.3 <- as.data.frame(do.call(rbind,summaries.3))
as.data.frame(df.3)
df.4 <- df.3 %>% select("terms")
df.5 <- apply(df.4, 2, as.character) 
model <- df.5[-1,]
model <- data.frame(model) 
significant.predictors.Water_Aphanizomenon_NIES81 <- cbind(model, df.1) %>% dplyr::filter(p.value < .1)
all.predictors.Water_Aphanizomenon_NIES81 <- cbind(model, df.1)
#visualizations
w.Aphan.1 <- linear.regression.plot(regression.df$Water_Aphanizomenon_NIES81, regression.df$Water_Microcystis_PCC7914) + labs(x = "RA Microcystis in H2O", y = "RA Aphanizomenon in H2O")
w.Aphan.2 <- linear.regression.plot(regression.df$Water_Aphanizomenon_NIES81, regression.df$Turbidity_avg) + labs(x = "Turbidity", y = "RA Aphanizomenon in H2O")
w.Aphan.3 <- linear.regression.plot(regression.df$Water_Aphanizomenon_NIES81, regression.df$PO4_avg)+ labs(x = "PO4", y = "RA Aphanizomenon in H2O")
w.Aphan.4 <- linear.regression.plot(regression.df$Water_Aphanizomenon_NIES81, regression.df$NH4_avg)+ labs(x = "NH4", y = "RA Aphanizomenon in H2O")
w.Aphan.5 <- linear.regression.plot(regression.df$Water_Aphanizomenon_NIES81, regression.df$Water_Dolichospermum_NIES41)+ labs(x = "RA Dolichospermum in H2O", y = "RA Aphanizomenon in H2O")
w.Aphan.6 <- linear.regression.plot(regression.df$Water_Aphanizomenon_NIES81, regression.df$CHLA_avg)+ labs(x = "CHLA", y = "RA Aphanizomenon in H2O")
w.Aphan.7 <- linear.regression.plot(regression.df$Water_Aphanizomenon_NIES81, regression.df$MC_avg)+ labs(x = "Microcystin conc.", y = "RA Aphanizomenon in H2O")
w.Aphan.8 <- linear.regression.plot(regression.df$Water_Aphanizomenon_NIES81, regression.df$TDN_avg)+ labs(x = "TDN", y = "RA Aphanizomenon in H2O")
w.Aphan.9 <- linear.regression.plot(regression.df$Water_Aphanizomenon_NIES81, regression.df$pH_avg)+ labs(x = "pH", y = "RA Aphanizomenon in H2O")
w.Aphan.10 <- linear.regression.plot(regression.df$Water_Aphanizomenon_NIES81, regression.df$PM_Cyanobium_PCC6307)+ labs(x = "RA Cyanobium in PM", y = "RA Aphanizomenon in H2O")

w.Aphan.1 + w.Aphan.2 + w.Aphan.3 + w.Aphan.4 + w.Aphan.5 + w.Aphan.6 + w.Aphan.7 + w.Aphan.8 + w.Aphan.9 + w.Aphan.10

#Outcome Variable: toxin producers in PM
colnames(regression.df)
regression.df  <- regression.df %>% select(PM_toxin.producers , everything()) #reordering df for for loop
columns <- colnames(regression.df)
lm.test <- list()
for(i in 2:ncol(regression.df)){
  lm.test[[i]] <- lm(reformulate(columns[i], "PM_toxin.producers"), na.action=na.exclude, data = regression.df)
}
summaries.1 <- lapply(lm.test, broom::glance)
df.1 <- as.data.frame(do.call(rbind,summaries.1))
summaries.3 <- lapply(lm.test, summary)
df.3 <- as.data.frame(do.call(rbind,summaries.3))
as.data.frame(df.3)
df.4 <- df.3 %>% select("terms")
df.5 <- apply(df.4, 2, as.character) 
model <- df.5[-1,]
model <- data.frame(model) 
significant.predictors.PM_toxin.producers <- cbind(model, df.1) %>% dplyr::filter(p.value < .1)
all.predictors.PM_toxin.producers <- cbind(model, df.1)
#visualizations
toxin.prod.1 <- linear.regression.plot(regression.df$PM_toxin.producers, regression.df$PM_Cyanobium_PCC6307)
toxin.prod.2 <- linear.regression.plot(regression.df$PM_toxin.producers, regression.df$NH4_avg)
toxin.prod.3 <- linear.regression.plot(regression.df$PM_toxin.producers, regression.df$TDN_avg)
toxin.prod.4 <- linear.regression.plot(regression.df$PM_toxin.producers, regression.df$EUK.S_avg)
toxin.prod.5 <- linear.regression.plot(regression.df$PM_toxin.producers, regression.df$NH4_avg)
toxin.prod.6 <- linear.regression.plot(regression.df$PM_toxin.producers, regression.df$PC.L_avg)
toxin.prod.7 <- linear.regression.plot(regression.df$PM_toxin.producers, regression.df$PC.S_avg)
toxin.prod.8 <- linear.regression.plot(regression.df$PM_toxin.producers, regression.df$C.to.N_avg)
toxin.prod.9 <- linear.regression.plot(regression.df$PM_toxin.producers, regression.df$Water_Dolichospermum_NIES41)

toxin.prod.1 + toxin.prod.2 + toxin.prod.3 + toxin.prod.4 + toxin.prod.5 + toxin.prod.6 + toxin.prod.7 + toxin.prod.8 + toxin.prod.9

#Outcome Variable: toxin producers in water
colnames(regression.df)
regression.df  <- regression.df %>% select(Water_toxin.producers , everything()) #reordering df for for loop
columns <- colnames(regression.df)
lm.test <- list()
for(i in 2:ncol(regression.df)){
  lm.test[[i]] <- lm(reformulate(columns[i], "Water_toxin.producers"), na.action=na.exclude, data = regression.df)
}
summaries.1 <- lapply(lm.test, broom::glance)
df.1 <- as.data.frame(do.call(rbind,summaries.1))
summaries.3 <- lapply(lm.test, summary)
df.3 <- as.data.frame(do.call(rbind,summaries.3))
as.data.frame(df.3)
df.4 <- df.3 %>% select("terms")
df.5 <- apply(df.4, 2, as.character) 
model <- df.5[-1,]
model <- data.frame(model) 
significant.predictors.Water_toxin.producers <- cbind(model, df.1) %>% dplyr::filter(p.value < .1)
all.predictors.Water_toxin.producers <- cbind(model, df.1)
#visualizations
w.toxin.prod.1 <- linear.regression.plot(regression.df$Water_toxin.producers, regression.df$NH4_avg)
w.toxin.prod.2 <- linear.regression.plot(regression.df$Water_toxin.producers, regression.df$TDN_avg)
w.toxin.prod.4 <- linear.regression.plot(regression.df$Water_toxin.producers, regression.df$PM_Cyanobium_PCC6307)
w.toxin.prod.5 <- linear.regression.plot(regression.df$Water_toxin.producers, regression.df$Water_Dolichospermum_NIES41)
w.toxin.prod.6 <- linear.regression.plot(regression.df$Water_toxin.producers, regression.df$MC_avg)
w.toxin.prod.7 <- linear.regression.plot(regression.df$Water_toxin.producers, regression.df$PO4_avg)
w.toxin.prod.8 <- linear.regression.plot(regression.df$Water_toxin.producers, regression.df$PM_toxin.producers)

w.toxin.prod.1 + w.toxin.prod.2 + w.toxin.prod.3 + w.toxin.prod.4 + w.toxin.prod.5 + w.toxin.prod.6 + w.toxin.prod.7

#Outcome Variable: PM.avg in aerosol
regression.df  <- regression.df %>% select(PM.avg , everything()) #reordering df for for loop
columns <- colnames(regression.df)
lm.test <- list()
for(i in 2:ncol(regression.df)){
  lm.test[[i]] <- lm(reformulate(columns[i], "PM.avg"), na.action=na.exclude, data = regression.df)
}
summaries.1 <- lapply(lm.test, broom::glance)
df.1 <- as.data.frame(do.call(rbind,summaries.1))
summaries.3 <- lapply(lm.test, summary)
df.3 <- as.data.frame(do.call(rbind,summaries.3))
as.data.frame(df.3)
df.4 <- df.3 %>% select("terms")
df.5 <- apply(df.4, 2, as.character) 
model <- df.5[-1,]
model <- data.frame(model) 
significant.predictors.PM.avg <- cbind(model, df.1) %>% dplyr::filter(p.value < .1)
all.predictors.PM.avg <- cbind(model, df.1)
#visualizations
PM.1 <- linear.regression.plot(regression.df$PM.avg, regression.df$Solar.Rad ) + labs(x = "Solar Irradiance", y = "Ambient PM conc.")
PM.2 <- linear.regression.plot(regression.df$PM.avg, regression.df$Air.Temp ) + labs(x = "Air Temp.", y = "Ambient PM conc.")
PM.3 <- linear.regression.plot(regression.df$PM.avg, regression.df$Salinity_avg ) + labs(x = "Salinity", y = "Ambient PM conc.")
PM.4 <- linear.regression.plot(regression.df$PM.avg, regression.df$Relative.Humidity ) + labs(x = "Relative Humidity", y = "Ambient PM conc.")
PM.5 <- linear.regression.plot(regression.df$PM.avg, regression.df$Water_Pseudanabaena_PCC7429 ) + labs(x = "RA Pseudanabaena in H2O", y = "Ambient PM conc.")
PM.6 <- linear.regression.plot(regression.df$PM.avg, regression.df$PM_z.Not.Assigned ) + labs(x = "RA Unassigned cyanobacteria in PM", y = "Ambient PM conc.")
PM.7 <- linear.regression.plot(regression.df$PM.avg, regression.df$Water_Dolichospermum_NIES41 ) + labs(x = "RA Dolicospermum in H2O", y = "Ambient PM conc.")
PM.8 <- linear.regression.plot(regression.df$PM.avg, regression.df$Water_Caenarcaniphilales_XX ) + labs(x = "RA Caenarcaniphilales in H2O", y = "Ambient PM conc.")
PM.9 <- linear.regression.plot(regression.df$PM.avg, regression.df$CHLA_avg ) + labs(x = "Chlorophyll a conc.", y = "Ambient PM conc.")
PM.10 <- linear.regression.plot(regression.df$PM.avg, regression.df$Wind.Speed) + labs(x = "Wind Speed", y = "Ambient PM conc.")
PM.11 <- linear.regression.plot(regression.df$PM.avg, regression.df$Water_Cyanobacteria ) + labs(x = "RA Cyanobacteria in H2O", y = "Ambient PM conc.")

PM.1 + PM.2 + PM.3 + PM.4 + PM.5 + PM.6 + PM.7 + PM.8 + PM.9 + PM.10 + PM.11

#Outcome Variable: Salinity
regression.df  <- regression.df %>% select(Salinity_avg , everything()) #reordering df for for loop
columns <- colnames(regression.df)
lm.test <- list()
for(i in 2:ncol(regression.df)){
  lm.test[[i]] <- lm(reformulate(columns[i], "Salinity_avg"), na.action=na.exclude, data = regression.df)
}
summaries.1 <- lapply(lm.test, broom::glance)
df.1 <- as.data.frame(do.call(rbind,summaries.1))
summaries.3 <- lapply(lm.test, summary)
df.3 <- as.data.frame(do.call(rbind,summaries.3))
as.data.frame(df.3)
df.4 <- df.3 %>% select("terms")
df.5 <- apply(df.4, 2, as.character) 
model <- df.5[-1,]
model <- data.frame(model) 
significant.predictors.Salinity_avg <- cbind(model, df.1) %>% dplyr::filter(p.value < .1)
all.predictors.Salinity_avg <- cbind(model, df.1)
#visualizations
salt.1 <- linear.regression.plot(regression.df$Salinity_avg, regression.df$PM.avg)
salt.2 <- linear.regression.plot(regression.df$Salinity_avg, regression.df$Water_Cyanobacteria)
salt.3 <- linear.regression.plot(regression.df$Salinity_avg, regression.df$Water_Pseudanabaena_PCC7429)
salt.4 <- linear.regression.plot(regression.df$Salinity_avg, regression.df$Bact_avg)
salt.5 <- linear.regression.plot(regression.df$Salinity_avg, regression.df$PM_z.Not.Assigned)
salt.6 <- linear.regression.plot(regression.df$Salinity_avg, regression.df$Diatom.CHLA_avg)

salt.1 + salt.2 + salt.3 + salt.4 + salt.5 + salt.6

#Outcome Variable: Microcystin 
colnames(regression.df)
regression.df  <- regression.df %>% select(MC_avg , everything()) #reordering df for for loop
columns <- colnames(regression.df)
lm.test <- list()
for(i in 2:ncol(regression.df)){
  lm.test[[i]] <- lm(reformulate(columns[i], "MC_avg"), na.action=na.exclude, data = regression.df)
}
summaries.1 <- lapply(lm.test, broom::glance)
df.1 <- as.data.frame(do.call(rbind,summaries.1))
summaries.3 <- lapply(lm.test, summary)
df.3 <- as.data.frame(do.call(rbind,summaries.3))
as.data.frame(df.3)
df.4 <- df.3 %>% select("terms")
df.5 <- apply(df.4, 2, as.character) 
model <- df.5[-1,]
model <- data.frame(model) 
significant.predictors.MC_avg <- cbind(model, df.1) %>% dplyr::filter(p.value < .1)
all.predictors.MC_avg <- cbind(model, df.1)
#visualizations
MC.1 <- linear.regression.plot(regression.df$MC_avg, regression.df$PO4_avg) + labs(x = "PO4", y = "Microcystin conc.")
MC.2 <- linear.regression.plot(regression.df$MC_avg, regression.df$CHLA_avg) + labs(x = "CHLA", y = "Microcystin conc.")
MC.3 <- linear.regression.plot(regression.df$MC_avg, regression.df$Water_Aphanizomenon_NIES81) + labs(x = "RA Aphanizomenon in H2O", y = "Microcystin conc.")
MC.4 <- linear.regression.plot(regression.df$MC_avg, regression.df$particulate.N_avg) + labs(x = "particulate N", y = "Microcystin conc.")
MC.5 <- linear.regression.plot(regression.df$MC_avg, regression.df$NH4_avg) + labs(x = "NH4", y = "Microcystin conc.")
MC.6 <- linear.regression.plot(regression.df$MC_avg, regression.df$TDN_avg) + labs(x = "TDN", y = "Microcystin conc.")

MC.1 + MC.2 + MC.3 + MC.4 + MC.5  + MC.6

# Visualizations by Site for robustness checks (ensuring trends aren't different between Sites)

# Cyanobacteria in PM
cyano.1.s <- linear.regression.plot.by.site(regression.df$PM_Cyanobacteria, regression.df$Diatom.CHLA_avg)
cyano.2.s <- linear.regression.plot.by.site(regression.df$PM_Cyanobacteria, regression.df$Cryptophyte.CHLA_avg)
cyano.3.s <- linear.regression.plot.by.site(regression.df$PM_Cyanobacteria, regression.df$PM_Cyanobium_PCC6307)
cyano.1.s + cyano.2.s + cyano.3.s 

# ambient PM
PM.1.s <- linear.regression.plot.by.site(regression.df$PM.avg, regression.df$Solar.Rad )
PM.2.s <- linear.regression.plot.by.site(regression.df$PM.avg, regression.df$Air.Temp )
PM.3.s <- linear.regression.plot.by.site(regression.df$PM.avg, regression.df$Salinity_avg )
PM.4.s <- linear.regression.plot.by.site(regression.df$PM.avg, regression.df$Relative.Humidity )
PM.5.s <- linear.regression.plot.by.site(regression.df$PM.avg, regression.df$Water_Pseudanabaena_PCC7429 )
PM.6.s <- linear.regression.plot.by.site(regression.df$PM.avg, regression.df$PM_z.Not.Assigned )
PM.7.s <- linear.regression.plot.by.site(regression.df$PM.avg, regression.df$Water_Dolichospermum_NIES41 )
PM.8.s <- linear.regression.plot.by.site(regression.df$PM.avg, regression.df$Water_Caenarcaniphilales_XX )
PM.9.s <- linear.regression.plot.by.site(regression.df$PM.avg, regression.df$CHLA_avg )
PM.10.s <- linear.regression.plot.by.site(regression.df$PM.avg, regression.df$Wind.Speed)
PM.11.s <- linear.regression.plot.by.site(regression.df$PM.avg, regression.df$Water_Cyanobacteria )
PM.1.s + PM.2.s + PM.3.s + PM.4.s + PM.5.s + PM.6.s + PM.7.s + PM.8.s + PM.9.s

# Examining impacts of Wind Direction



