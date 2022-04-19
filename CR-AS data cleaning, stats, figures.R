#CR-AS data cleaning, stats, and figures for manuscript
# *** NOTE WHEN RE-RUNNING CODE NEED TO TAKE # OFF OF LINES 492 & 494 -> PUT IN PLACE SO THEY ARE NOT RE-RUN MID SCRIPT ***
rm(list=ls()) #clear environment
setwd("/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs") #set working directory 
#install all necessary packages
library(dplyr); library(tidyr); library(ggplot2); library(gam); library(gamm4); library(cowplot); library(mgcv); library(reshape2); library(MuMIn); library(stringr); library(ISLR); library(voxel); library(gridExtra); library(purrr); library(data.table); library(phyloseq); library(decontam); library(RAM); library(tidyverse); library(DESeq2); library(microbiome); library(vegan); library(viridis); library(patchwork); library(gapminder); library(tidyverse); library(ape); library(RColorBrewer); library(ggpubr);library(rstatix)

#Reading in the dada2 datasets 
aerosol <- read.csv("Aerosol_16S_counts.csv", header = TRUE, na.strings = c(""," ", ".", "NA"))
names(aerosol)<-sapply(str_remove_all(colnames(aerosol),"AERO_"),"[")
water <- read.csv("Water_16S_counts.csv", header = TRUE, na.strings = c(""," ", ".", "NA"))  
names(water)<-sapply(str_remove_all(colnames(water),"W_"),"[")
ASV.w.taxonomy <- read.csv("/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/PhyloseqReadyFiles/Taxonomy_Table.csv", header=TRUE)
#inputting ASVs of interest species assignments from BLASTn search 
ASV.w.taxonomy$Species[0196]="Leptolyngbya"
ASV.w.taxonomy$Species[0389]="Anabaena"
ASV.w.taxonomy$Species[0108]="Haloleptolyngbya"
ASV.w.taxonomy$Species[1999]="Anabaena"

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
RA.df <- CSS.3 %>% left_join(ASV.w.taxonomy, by = "ASV") %>% mutate_if(is.numeric, ~replace_na(., 0)) 

#counting overlapping ASVs between PM and water samples with taxonomic descriptions
PM.samples <- RA.df %>% dplyr::select('ASV', 'S1_A', 'S2_A', 'S2_B', 'S3_A', 'S3_B', 'S4_A', 'S4_B', 'S5_A', 'S5_B', 'S6_B', 'S7_A', 'S7_B', 'S8_A', 'S8_B') %>% column_to_rownames("ASV")
PM.samples.1 <- PM.samples[rowSums(PM.samples[])>0,]
PM.samples.2 <- rownames_to_column(PM.samples.1) 
water.samples <- RA.df %>% dplyr::select('ASV', "06.11_A", "06.11_B", "06.23_A", "06.23_B", "07.07_A", "07.07_B", "07.21_A", "07.21_B", "08.12_A", "08.12_B", "08.18_A", "08.18_B", "09.01_A", "09.01_B","09.15_A", "09.15_B", "10.01_A", "10.01_B") %>% column_to_rownames("ASV")
water.samples.1 <- water.samples[rowSums(water.samples[])>0,]
water.samples.2 <- rownames_to_column(water.samples.1) 
overlapping.ASVs <- left_join(water.samples.2, PM.samples.2, by = "rowname") %>% na.omit() %>% mutate(ASV = rowname) %>% left_join(ASV.w.taxonomy) %>% select('ASV', 'S1_A', 'S2_A', 'S2_B', 'S3_A', 'S3_B', 'S4_A', 'S4_B', 'S5_A', 'S5_B', 'S6_B', 'S7_A', 'S7_B', 'S8_A', 'S8_B', "06.11_A", "06.11_B", "06.23_A", "06.23_B", "07.07_A", "07.07_B", "07.21_A", "07.21_B", "08.12_A", "08.12_B", "08.18_A", "08.18_B", "09.01_A", "09.01_B","09.15_A", "09.15_B", "10.01_A", "10.01_B", "Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species", "NA.")
overlapping.cyanos <- overlapping.ASVs %>% dplyr::filter(Class == "Cyanobacteria")

#Estimating water values during aerosol sampling periods rather than at interval points
write.csv(RA.df, "/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/files.for.manipulation/to.calculate.water.during.sampling.periods.csv")
water.sampling.periods <- read.csv("/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/water.sampling.period.csv", header=TRUE)
RA.df.0 <- RA.df %>% left_join(water.sampling.periods, by = "ASV", keep = FALSE) 
RA.df.SPs <- RA.df.0 %>% dplyr::select("ASV", 'S1_A', 'S2_A', 'S2_B', 'S3_A', 'S3_B', 'S4_A', 'S4_B', 'S5_A', 'S5_B', 'S6_B', 'S7_A', 'S7_B', 'S8_A', 'S8_B', "W1_A", "W2_A" , "W2_B" , "W3_A" , "W3_B" , "W4_A" , "W4_B" , "W5_A" , "W5_B",  "W6_B",  "W7_A"  ,"W7_B"  ,"W8_A",  "W8_B", "Kingdom", "Phylum"  ,"Class" ,  "Order" ,  "Family" , "Genus"  , "Species", "NA." )
overlapping.cyanos.ASVs <- overlapping.cyanos %>% dplyr::select(ASV)
overlapping.cyanos.ASVs <- overlapping.cyanos.ASVs %>% left_join(RA.df.SPs, by = "ASV") %>% dplyr::select("ASV", "W1_A", "W2_A" , "W2_B" , "W3_A" , "W3_B" , "W4_A" , "W4_B" , "W5_A" , "W5_B",  "W6_B",  "W7_A"  ,"W7_B"  ,"W8_A",  "W8_B", "Kingdom", "Phylum"  ,"Class" ,  "Order" ,  "Family" , "Genus"  , "Species", "NA.")
overlapping.cyanos.1 <- overlapping.cyanos %>% left_join(overlapping.cyanos.ASVs) %>% select("ASV", 'S1_A', 'S2_A', 'S2_B', 'S3_A', 'S3_B', 'S4_A', 'S4_B', 'S5_A', 'S5_B', 'S6_B', 'S7_A', 'S7_B', 'S8_A', 'S8_B', "W1_A", "W2_A" , "W2_B" , "W3_A" , "W3_B" , "W4_A" , "W4_B" , "W5_A" , "W5_B",  "W6_B",  "W7_A"  ,"W7_B"  ,"W8_A",  "W8_B", "Kingdom", "Phylum"  ,"Class" ,  "Order" ,  "Family" , "Genus"  , "Species", "NA.")

## SUGGESTION FROM COMMITEE -- RE-RUN ALL ANALYSES ONLY THE START TIME WATER SAMPLES AND THE END TIME WATER SAMPLES RATHER THAN AN AVERAGE TO SEE IF TRENDS / FINDINGS ARE MAJORLY DIFFERENT 
#using interval start and stop time samples rather than averages for water conditions
overlapping.cyanos.start.times <- RA.df.0 %>% dplyr::select("ASV", 'S1_A', 'S2_A', 'S2_B', 'S3_A', 'S3_B', 'S4_A', 'S4_B', 'S5_A', 'S5_B', 'S6_B', 'S7_A', 'S7_B', 'S8_A', 'S8_B',"06.11_A", "06.11_B", "06.23_A", "06.23_B", "07.07_A", "07.07_B", "07.21_A", "07.21_B", "08.12_A", "08.12_B", "08.18_A", "08.18_B", "09.01_A", "09.01_B", "09.15_A", "09.15_B", "10.01_A", "10.01_B", "Kingdom", "Phylum"  ,"Class" ,  "Order" ,  "Family" , "Genus"  , "Species", "NA.")
overlapping.cyanos.start.times.ASVs <- overlapping.cyanos.ASVs %>% left_join(overlapping.cyanos.start.times, by = "ASV") %>% dplyr::select("ASV", "06.11_A", "06.11_B", "06.23_A", "06.23_B", "07.07_A", "07.07_B", "07.21_A", "07.21_B", "08.12_A", "08.12_B", "08.18_A", "08.18_B", "09.01_A", "09.01_B", "09.15_A", "09.15_B", "10.01_A", "10.01_B", "Kingdom", "Phylum"  ,"Class" ,  "Order" ,  "Family" , "Genus"  , "Species", "NA.")
overlapping.cyanos.start.times.ASVs.1 <- overlapping.cyanos %>% left_join(overlapping.cyanos.ASVs) %>% select("ASV", 'S1_A', 'S2_A', 'S2_B', 'S3_A', 'S3_B', 'S4_A', 'S4_B', 'S5_A', 'S5_B', 'S6_B', 'S7_A', 'S7_B', 'S8_A', 'S8_B', "06.11_A", "06.11_B", "06.23_A", "06.23_B", "07.07_A", "07.07_B", "07.21_A", "07.21_B", "08.12_A", "08.12_B", "08.18_A", "08.18_B", "09.01_A", "09.01_B", "09.15_A", "09.15_B", "10.01_A", "10.01_B", "Kingdom", "Phylum"  ,"Class" ,  "Order" ,  "Family" , "Genus"  , "Species", "NA.")

#Plot Color Schemes
brewer.pal(n = 12, name = 'Paired')
class.colors <- c("Acidobacteria" = "#4E79A7",
           "Actinobacteria" = "#F28E2B",
           "Bacteroidetes"   = "#E15759",
           "Cyanobacteria"   = "#76B7B2",
           "Firmicutes"  = "#59A14F",
           "Planctomycetes" = "#EDC948",
           "Proteobacteria" = "#FF9DA7",
           "z.Other" = "#BAB0AC")

species.colors <- c("Anabaena" = "#4E79A7", 
                    "Aphanizomenon" = "#F28E2B",
                    "Caenarcaniphilales" = "#E15759",
                    "Chroococcidiopsis" = "#009E73",
                    "Cyanobium"  = "#8CD17D",
                    "Dolichospermum" = "#EDC948",
                    "Haloleptolyngbya" = "#B07AA1",
                    "Leptolyngbya" = "#FF9DA7",  
                    "Microcystis" = "#9C755F",
                    "Pseudanabaena" = "#A0CBE8",
                    "Tolypothrix" = "#FFBE7D",
                    "z.Not.Assigned" = "#D4A6C8",
                    "z.Other" = "#BAB0AC")

# plotting relative abundance of bacterial communities by CSS, with water values shown at intervals
RA.df.1 <- RA.df %>% group_by(Class) %>% summarise_if(is.numeric, ~sum(na.exclude(.)))
RA.df.2 <- RA.df.1 %>% dplyr::select(S1_A:`10.01_B`)
RA.df.3 <- RA.df.2 %>% summarise(across(S1_A:`10.01_B`, sum))
RA.df.4 <- mapply(`/`, RA.df.2, RA.df.3)
Class.26 <- RA.df.1 %>% dplyr::select("Class") %>% unique()
RA.df.5 <- cbind(Class.26, RA.df.4)
RA.df.6 <- RA.df.5 %>% mutate_at(vars(Class), ~replace_na(., "z.Not.Assigned"))
write.csv(RA.df.6, "/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/files.for.manipulation/Relative.Abundance.Class.interval.csv")
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
RA.df.8 <- RA.df.7 %>% group_by(Class) %>% summarise(across(S1_A:`10.01_B`, sum))
RA.df.9 <- pivot_longer(RA.df.8, cols = c('S1_A', 'S2_A', 'S2_B', 'S3_A', 'S3_B', 'S4_A', 'S4_B', 'S5_A', 'S5_B', 'S6_B', 'S7_A', 'S7_B', 'S8_A', 'S8_B',"06.11_A", "06.11_B", "06.23_A","06.23_B", "07.07_A", "07.07_B", "07.21_A", "07.21_B", "08.12_A", "08.12_B", "08.18_A", "08.18_B","09.01_A", "09.01_B", "09.15_A", "09.15_B", "10.01_A", "10.01_B"), names_to = "Sample") %>% mutate(Sample_Type = case_when(str_detect(Sample, "S") ~ "PM", str_detect(Sample, "0") ~ "Water")) %>% mutate(Site = case_when(str_detect(Sample, "_A") ~ "A", str_detect(Sample, "_B") ~ "B")) 
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

RA.plot.Water.1 <- ggplot(RA.df.Water, aes(fill = Class, y=value, x=Sample)) + 
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

RA.plot.PM + RA.plot.Water.1
  
#plotting relative abundance of bacterial communities by CSS, with water values averaged for aerosol sampling periods
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

## Plotting relative abundance of cyanobacterial genera by CSS, with water values shown at intervals
species.RA.df.1 <- RA.df %>% filter(Class == "Cyanobacteria") %>% group_by(Species) %>% summarise_if(is.numeric, ~sum(na.exclude(.)))
species.RA.df.2 <- species.RA.df.1 %>% dplyr::select(S1_A:`10.01_B`)
species.RA.df.3 <- species.RA.df.2 %>% summarise(across(S1_A:`10.01_B`, sum))
species.RA.df.4 <- mapply(`/`, species.RA.df.2, species.RA.df.3)
Species <- species.RA.df.1 %>% dplyr::select("Species") %>% unique()
species.RA.df.5 <- cbind(Species, species.RA.df.4)
species.RA.df.6 <- species.RA.df.5 %>% mutate_at(vars(Species), ~replace_na(., "z.Not.Assigned")) 
write.csv(species.RA.df.6, "/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/files.for.manipulation/Relative.Abundance.species.intervals.csv")
species.RA.df.6$Species
species.RA.df.7 <- species.RA.df.6 %>% mutate(Species = case_when(
  Species == "Anabaena_PCC7108" ~ "Anabaena",
  Species == "Anabaena" ~ "Anabaena",
  Species == "Aphanizomenon_NIES81" ~ "Aphanizomenon",
  Species == "Caenarcaniphilales_XX" ~ "Caenarcaniphilales",
  Species == "Caldora_VP642b"~ "z.Other",
  Species == "Calothrix_PCC6303"~ "z.Other",
  Species == "Cephalothrix_SAG_75.79"~ "z.Other",    
  Species == "Chroococcidiopsis_PCC_7203"~ "z.Other",
  Species == "Chroococcidiopsis_SAG_2023"~ "Chroococcidiopsis",
  Species == "Cyanobium_PCC6307" ~ "Cyanobium",
  Species == "Cyanothece"~ "z.Other",
  Species == "Dolichospermum_NIES41" ~ "Dolichospermum",
  Species == "Gastranaerophilales_XX"~ "z.Other",
  Species == "Geminocystis_PCC6308"~ "z.Other",
  Species == "Gloeocapsa"~ "z.Other",
  Species == "Haloleptolyngbya"~ "Haloleptolyngbya",
  Species == "Leptolyngbya" ~ "Leptolyngbya",        
  Species == "Microcystis_PCC7914" ~ "Microcystis",
  Species == "Nodularia_PCC9350"~ "z.Other",
  Species == "Nostoc_PCC73102"~ "z.Other",
  Species == "Nostoc_PCC7524"~ "z.Other",
  Species == "Obscuribacterales_XX" ~ "z.Other",
  Species == "Phormidium"~ "z.Other",
  Species == "Phormidium_IAM_M71"~ "z.Other",        
  Species == "Pleurocapsa"~ "z.Other",
  Species == "Pseudanabaena_PCC7429"~ "Pseudanabaena",
  Species == "Synechococcus_CC9902" ~ "z.Other",
  Species == "Synechococcus_PCC7942" ~ "z.Other",
  Species == "Tolypothrix_IAM_M259" ~ "z.Other",
  Species == "Tolypothrix_PCC7601" ~ "Tolypothrix",
  Species == "z.Not.Assigned" ~ "z.Not.Assigned")) %>% group_by(Species) %>% summarise(across(S1_A:`10.01_B`, sum))
species.RA.df.8 <- pivot_longer(species.RA.df.7, cols = c('S1_A', 'S2_A', 'S2_B', 'S3_A', 'S3_B', 'S4_A', 'S4_B', 'S5_A', 'S5_B', 'S6_B', 'S7_A', 'S7_B', 'S8_A', 'S8_B',"06.11_A", "06.11_B", "06.23_A","06.23_B", "07.07_A", "07.07_B", "07.21_A", "07.21_B", "08.12_A", "08.12_B", "08.18_A", "08.18_B","09.01_A", "09.01_B", "09.15_A", "09.15_B", "10.01_A", "10.01_B"), names_to = "Sample") %>% mutate(Sample_Type = case_when(str_detect(Sample, "S") ~ "PM", str_detect(Sample, "0") ~ "Water")) %>% mutate(Site = case_when(str_detect(Sample, "_A") ~ "A", str_detect(Sample, "_B") ~ "B")) %>% mutate_if(is.numeric, ~replace_na(., 0))
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

species.RA.plot.Water.1 <- ggplot(species.RA.df.Water, aes(fill = Species, y=value, x=Sample)) + 
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

species.RA.plot.PM + species.RA.plot.Water.1

## Plotting relative abundance of cyanobacterial genera by CSS, with water values shown as averages over PM sampling periods
species.RA.df.1 <- RA.df.SPs %>% filter(Class == "Cyanobacteria") %>% group_by(Species) %>% summarise_if(is.numeric, ~sum(na.exclude(.)))
species.RA.df.2 <- species.RA.df.1 %>% dplyr::select(S1_A:`W8_B`)
species.RA.df.3 <- species.RA.df.2 %>% summarise(across(S1_A:`W8_B`, sum))
species.RA.df.4 <- mapply(`/`, species.RA.df.2, species.RA.df.3)
Species <- species.RA.df.1 %>% dplyr::select("Species") %>% unique()
species.RA.df.5 <- cbind(Species, species.RA.df.4)
species.RA.df.6 <- species.RA.df.5 %>% mutate_at(vars(Species), ~replace_na(., "z.Not.Assigned")) 
write.csv(species.RA.df.6, "/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/files.for.manipulation/Relative.Abundance.species.csv")
species.RA.df.7 <- species.RA.df.6 %>% mutate(Species = case_when(
  Species == "Anabaena_PCC7108" ~ "Anabaena",
  Species == "Anabaena" ~ "Anabaena",
  Species == "Aphanizomenon_NIES81" ~ "Aphanizomenon",
  Species == "Caenarcaniphilales_XX" ~ "Caenarcaniphilales",
  Species == "Caldora_VP642b"~ "z.Other",
  Species == "Calothrix_PCC6303"~ "z.Other",
  Species == "Cephalothrix_SAG_75.79"~ "z.Other",    
  Species == "Chroococcidiopsis_PCC_7203"~ "z.Other",
  Species == "Chroococcidiopsis_SAG_2023"~ "Chroococcidiopsis",
  Species == "Cyanobium_PCC6307" ~ "Cyanobium",
  Species == "Cyanothece"~ "z.Other",
  Species == "Dolichospermum_NIES41" ~ "Dolichospermum",
  Species == "Gastranaerophilales_XX"~ "z.Other",
  Species == "Geminocystis_PCC6308"~ "z.Other",
  Species == "Gloeocapsa"~ "z.Other",
  Species == "Haloleptolyngbya"~ "Haloleptolyngbya",
  Species == "Leptolyngbya" ~ "Leptolyngbya",    
  Species == "Microcystis_PCC7914" ~ "Microcystis",
  Species == "Nodularia_PCC9350"~ "z.Other",
  Species == "Nostoc_PCC73102"~ "z.Other",
  Species == "Nostoc_PCC7524"~ "z.Other",
  Species == "Obscuribacterales_XX" ~ "z.Other",
  Species == "Phormidium"~ "z.Other",
  Species == "Phormidium_IAM_M71"~ "z.Other",        
  Species == "Pleurocapsa"~ "z.Other",
  Species == "Pseudanabaena_PCC7429"~ "Pseudanabaena",
  Species == "Synechococcus_CC9902" ~ "z.Other",
  Species == "Synechococcus_PCC7942" ~ "z.Other",
  Species == "Tolypothrix_IAM_M259" ~ "z.Other",
  Species == "Tolypothrix_PCC7601" ~ "Tolypothrix",
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

# examining Alpha and Beta Diversity for PM and water samples, for supplementary material
alpha.diversity <- estimate_richness(phyloseq.ob)
Observed <- plot_richness(phyloseq.ob,x = "Sample_Type" , color = "Sample_Type", measures = c("Observed")) + geom_boxplot() + theme_bw() + theme(text = element_text(size=7), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7, angle = 45, hjust = 1), legend.text = element_text(size = 6), legend.key.size = unit(2, "mm"), legend.position = "none") + labs(title = "Alpha Diversity")
Shannon <- plot_richness(phyloseq.ob,x = "Sample_Type" , color = "Sample_Type", measures = c("Shannon")) + geom_boxplot()+ theme_bw() + theme(text = element_text(size=7), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7, angle = 45, hjust = 1), legend.text = element_text(size = 6), legend.key.size = unit(2, "mm"), legend.position = "none") 
InvSimpson <- plot_richness(phyloseq.ob,x = "Sample_Type" , color = "Sample_Type", measures = c("InvSimpson")) + geom_boxplot()+ theme_bw() + theme(text = element_text(size=7), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7, angle = 45, hjust = 1), legend.text = element_text(size = 6), legend.key.size = unit(2, "mm"), legend.position = "none")

Observed + Shannon + InvSimpson

# T-test of Diversity Index results
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
# beta diversity 
wunifrac_dist <- phyloseq::distance(phyloseq.ob, method = "unifrac", weighted = F)
ordination <- ordinate(phyloseq.ob, method="PCoA", distance=wunifrac_dist)
plot_ordination(phyloseq.ob, ordination, color="Sample_Type") + theme(aspect.ratio=1)

#-----------------------------------------------------------------------------------------------------------------------------------------#
# Calculating aerosolization factors (AF) for all classes of bacteria and species of cyanobacteria
RA.df.AF <- pivot_longer(RA.df.6, cols = c('S1_A', 'S2_A', 'S2_B', 'S3_A', 'S3_B', 'S4_A', 'S4_B', 'S5_A', 'S5_B', 'S6_B', 'S7_A', 'S7_B', 'S8_A', 'S8_B', "W1_A", "W2_A" , "W2_B" , "W3_A" , "W3_B" , "W4_A" , "W4_B" , "W5_A" , "W5_B",  "W6_B",  "W7_A"  ,"W7_B"  ,"W8_A",  "W8_B"), names_to = "Sample") %>% mutate(Sample_Type = case_when(str_detect(Sample, "S") ~ "PM", str_detect(Sample, "W") ~ "Water")) %>% mutate(Site = case_when(str_detect(Sample, "_A") ~ "A", str_detect(Sample, "_B") ~ "B")) %>% mutate_if(is.numeric, ~replace_na(., 0)) %>% mutate(Sampling_Period = case_when(str_detect(Sample, "1") ~ "1", str_detect(Sample, "2") ~ "2", str_detect(Sample, "3") ~ "3", str_detect(Sample, "4") ~ "4", str_detect(Sample, "5") ~ "5", str_detect(Sample, "6") ~ "6", str_detect(Sample, "7") ~ "7", str_detect(Sample, "8") ~ "8")) %>% dplyr::select(-Sample) 
class.RAs <- RA.df.AF %>% group_by(Class, Sampling_Period, Site) %>% pivot_wider(names_from = Sample_Type, values_from = value)  
class.AFs <- class.RAs %>% mutate(`r-AF` = PM/Water) %>% mutate(`log.r-AF` = log(`r-AF`)) %>% mutate(`r-AF`= na_if(`r-AF`, "Inf")) %>% mutate(`r-AF`= na_if(`r-AF`, "NaN")) %>% unite(Sampling_Period, Sampling_Period, Site, remove= TRUE) %>% mutate(`log.r-AF` = log(`r-AF`)) %>% mutate_all(~replace(., is.na(.), 0)) %>% mutate_all(~replace(., is.infinite(.), 0)) 

species.RA.df.AF <- pivot_longer(species.RA.df.6, cols = c('S1_A', 'S2_A', 'S2_B', 'S3_A', 'S3_B', 'S4_A', 'S4_B', 'S5_A', 'S5_B', 'S6_B', 'S7_A', 'S7_B', 'S8_A', 'S8_B', "W1_A", "W2_A" , "W2_B" , "W3_A" , "W3_B" , "W4_A" , "W4_B" , "W5_A" , "W5_B",  "W6_B",  "W7_A"  ,"W7_B"  ,"W8_A",  "W8_B"), names_to = "Sample") %>% mutate(Sample_Type = case_when(str_detect(Sample, "S") ~ "PM", str_detect(Sample, "W") ~ "Water")) %>% mutate(Site = case_when(str_detect(Sample, "_A") ~ "A", str_detect(Sample, "_B") ~ "B")) %>% mutate_if(is.numeric, ~replace_na(., 0)) %>% mutate(Sampling_Period = case_when(str_detect(Sample, "1") ~ "1", str_detect(Sample, "2") ~ "2", str_detect(Sample, "3") ~ "3", str_detect(Sample, "4") ~ "4", str_detect(Sample, "5") ~ "5", str_detect(Sample, "6") ~ "6", str_detect(Sample, "7") ~ "7", str_detect(Sample, "8") ~ "8")) %>% dplyr::select(-Sample)
species.RAs <- species.RA.df.AF %>% group_by(Species, Sampling_Period, Site) %>% pivot_wider(names_from = Sample_Type, values_from = value)
species.AFs <- species.RAs %>% mutate(`r-AF` = PM/Water) %>% mutate(`log.r-AF` = log(`r-AF`)) %>% mutate(`r-AF`= na_if(`r-AF`, "Inf")) %>% mutate(`r-AF`= na_if(`r-AF`, "NaN")) %>% unite(Sampling_Period, Sampling_Period, Site, remove= TRUE) %>% mutate_all(~replace(., is.na(.), 0)) %>% mutate_all(~replace(., is.infinite(.), 0)) 

## Q1 FOR STATS CONSULTANT
# the relative aerosolization factor (rAF or r-AF) is the ratio of abundance of each cyanobacteria in air to its abundance in water. This allows us to examine differential species enrichment in spray aerosol. The data has many zeros and is left-skewed. I report the aerosolization factors on average, but for better visualization in the heat-maps, I added 1 and log-transformed. Is this an OK practice? Since it's vizualized this way is it misleading to show the values as non-transformed averages? 
ASV.RA.df.AF <- pivot_longer(overlapping.cyanos.1, cols = c('S1_A', 'S2_A', 'S2_B', 'S3_A', 'S3_B', 'S4_A', 'S4_B', 'S5_A', 'S5_B', 'S6_B', 'S7_A', 'S7_B', 'S8_A', 'S8_B', "W1_A", "W2_A" , "W2_B" , "W3_A" , "W3_B" , "W4_A" , "W4_B" , "W5_A" , "W5_B",  "W6_B",  "W7_A"  ,"W7_B"  ,"W8_A",  "W8_B"), names_to = "Sample")  %>% mutate(Sample_Type = case_when(str_detect(Sample, "S") ~ "PM", str_detect(Sample, "W") ~ "Water")) %>% mutate(Site = case_when(str_detect(Sample, "_A") ~ "A", str_detect(Sample, "_B") ~ "B")) %>% mutate_if(is.numeric, ~replace_na(., 0)) %>% mutate(Sampling_Period = case_when(str_detect(Sample, "1") ~ "1", str_detect(Sample, "2") ~ "2", str_detect(Sample, "3") ~ "3", str_detect(Sample, "4") ~ "4", str_detect(Sample, "5") ~ "5", str_detect(Sample, "6") ~ "6", str_detect(Sample, "7") ~ "7", str_detect(Sample, "8") ~ "8")) %>% dplyr::select(-Sample)
ASV.RAs <- ASV.RA.df.AF %>% group_by(ASV, Sampling_Period, Site) %>% pivot_wider(names_from = Sample_Type, values_from = value) 
ASV.AF.without.log.transformation <- ASV.RAs %>% mutate(`r-AF` = PM/Water) %>% mutate(`r-AF`= na_if(`r-AF`, "Inf")) %>% mutate(`r-AF`= na_if(`r-AF`, "NaN")) %>% unite(Sampling_Period, Sampling_Period, Site, remove= TRUE) %>% group_by(ASV) %>% dplyr::summarise_at(vars(`r-AF`), ~mean(., na.rm = T)) %>% left_join(ASV.w.taxonomy, by = "ASV") %>% dplyr::select(ASV, `r-AF`, Species)
ASV.AFs <- ASV.RAs %>% mutate(`r-AF` = (PM)/(Water)) %>% mutate(`r-AF`= na_if(`r-AF`, "NaN")) %>% mutate(`r-AF`= na_if(`r-AF`, "NA")) %>% mutate(`r-AF`= na_if(`r-AF`, "Inf")) %>% mutate_at(vars(`r-AF`), ~replace_na(., 0)) %>% mutate(`1+r-AF` = (1+PM)/(1+Water)) %>% mutate(`log.r-AF` = log(`r-AF`)) %>% mutate(`log.1+r-AF` = log(`1+r-AF`)) %>% mutate(`r-AF`= na_if(`r-AF`, "Inf")) %>% unite(Sampling_Period, Sampling_Period, Site, remove= TRUE) %>% filter(ASV != "ASV4107" , ASV != "ASV2139" , ASV != "ASV1848" , ASV != "ASV0986") 
AF.summary.stats <- ASV.AFs %>% group_by(ASV) %>% get_summary_stats(`r-AF`, type = "mean_sd") %>% left_join( ASV.w.taxonomy, by = "ASV") %>% select(ASV, variable, n, mean, sd, Species)

#various version of heat map using different data transformations to show enrichment in aerosol
ggplot(ASV.AFs, aes(Sampling_Period, ASV, fill= `r-AF`)) + geom_tile()+ scale_fill_gradientn(colours = c("darkblue", "beige", "darkred"))
ggplot(ASV.AFs, aes(Sampling_Period, ASV, fill= `log.r-AF`)) + geom_tile()+ scale_fill_gradientn(colours = c("darkblue", "beige", "darkred")) + theme_light() + scale_y_discrete(limits = rev(c("ASV1111", "ASV0389", "ASV1999", "ASV0309", "ASV2080","ASV1397", "ASV0546", "ASV0288", "ASV0529", "ASV0413","ASV0108","ASV0196","ASV0889", "ASV2444", "ASV0467")))
ggplot(ASV.AFs, aes(Sampling_Period, ASV, fill= `1+r-AF`)) + geom_tile()+ scale_fill_gradientn(colours = c("darkblue", "beige", "darkred")) + scale_y_discrete(limits = rev(c("ASV1111", "ASV0389", "ASV1999", "ASV0309", "ASV2080","ASV1397", "ASV0546", "ASV0288", "ASV0529", "ASV0413","ASV0108","ASV0196","ASV0889", "ASV2444", "ASV0467")))
ggplot(ASV.AFs, aes(Sampling_Period, ASV, fill= `log.1+r-AF`)) + geom_tile()+ scale_fill_gradientn(colours = c("darkblue", "beige", "darkred")) + theme_light() + scale_y_discrete(limits = rev(c("ASV1111", "ASV0389", "ASV1999", "ASV0309", "ASV2080","ASV1397", "ASV0546", "ASV0288", "ASV0529", "ASV0413","ASV0108","ASV0196","ASV0889", "ASV2444", "ASV0467")))

#r-AF values calculated with START TIMES rather than averages
start.times <- overlapping.cyanos.start.times.ASVs.1 %>% dplyr::rename("W1_A"= "06.11_A", "W1_B" = "06.11_B", "W2_A" = "06.23_A", "W2_B" = "06.23_B", "W3_A" = "07.07_A", "W3_B" = "07.07_B", "W4_A" = "07.21_A", "W4_B" = "07.21_B", "W5_A" = "08.12_A", "W5_B" = "08.12_B", "W6_A" = "08.18_A", "W6_B" = "08.18_B", "W7_A" = "09.01_A", "W7_B" = "09.01_B", "W8_A" = "09.15_A", "W8_B" = "09.15_B")
ASV.RA.df.AF.start <- pivot_longer(start.times, cols = c('S1_A', 'S2_A', 'S2_B', 'S3_A', 'S3_B', 'S4_A', 'S4_B', 'S5_A', 'S5_B', 'S6_B', 'S7_A', 'S7_B', 'S8_A', 'S8_B', "W1_A", "W2_A" , "W2_B" , "W3_A" , "W3_B" , "W4_A" , "W4_B" , "W5_A" , "W5_B",  "W6_B",  "W7_A"  ,"W7_B"  ,"W8_A",  "W8_B"), names_to = "Sample")  %>% mutate(Sample_Type = case_when(str_detect(Sample, "S") ~ "PM", str_detect(Sample, "W") ~ "Water")) %>% mutate(Site = case_when(str_detect(Sample, "_A") ~ "A", str_detect(Sample, "_B") ~ "B")) %>% mutate_if(is.numeric, ~replace_na(., 0)) %>% mutate(Sampling_Period = case_when(str_detect(Sample, "1") ~ "1", str_detect(Sample, "2") ~ "2", str_detect(Sample, "3") ~ "3", str_detect(Sample, "4") ~ "4", str_detect(Sample, "5") ~ "5", str_detect(Sample, "6") ~ "6", str_detect(Sample, "7") ~ "7", str_detect(Sample, "8") ~ "8")) %>% dplyr::select(-Sample)
ASV.RAs.start <- ASV.RA.df.AF.start %>% group_by(ASV, Sampling_Period, Site) %>% pivot_wider(names_from = Sample_Type, values_from = value) 
ASV.AF.without.log.transformation.start <- ASV.RAs.start %>% mutate(`r-AF` = PM/Water) %>% mutate(`r-AF`= na_if(`r-AF`, "Inf")) %>% mutate(`r-AF`= na_if(`r-AF`, "NaN")) %>% unite(Sampling_Period, Sampling_Period, Site, remove= TRUE) %>% group_by(ASV) %>% dplyr::summarise_at(vars(`r-AF`), ~mean(., na.rm = T)) %>% left_join(ASV.w.taxonomy, by = "ASV") %>% dplyr::select(ASV, `r-AF`, Species)
ASV.AFs.start <- ASV.RAs.start %>% mutate(`r-AF` = (PM)/(Water)) %>% mutate(`r-AF`= na_if(`r-AF`, "NaN")) %>% mutate(`r-AF`= na_if(`r-AF`, "NA")) %>% mutate(`r-AF`= na_if(`r-AF`, "Inf")) %>% mutate_at(vars(`r-AF`), ~replace_na(., 0)) %>% mutate(`1+r-AF` = (1+PM)/(1+Water)) %>% mutate(`log.r-AF` = log(`r-AF`)) %>% mutate(`log.1+r-AF` = log(`1+r-AF`)) %>% mutate(`r-AF`= na_if(`r-AF`, "Inf")) %>% unite(Sampling_Period, Sampling_Period, Site, remove= TRUE) %>% filter(ASV != "ASV4107" , ASV != "ASV2139" , ASV != "ASV1848" , ASV != "ASV0986") 
AF.summary.stats.start <- ASV.AFs.start %>% group_by(ASV) %>% get_summary_stats(`r-AF`, type = "mean_sd") %>% left_join( ASV.w.taxonomy, by = "ASV") %>% select(ASV, variable, n, mean, sd, Species)

#r-AF with STOP TIMES instead of averages 
stop.times <- overlapping.cyanos.start.times.ASVs.1 %>% dplyr::rename("W1_A"= "06.23_A", "W1_B" = "06.23_B", "W2_A" = "07.07_A", "W2_B" = "07.07_B", "W3_A" = "07.21_A", "W3_B" = "07.21_B", "W4_A" = "08.12_A", "W4_B" = "08.12_B", "W5_A" = "08.18_A", "W5_B" = "08.18_B", "W6_A" = "09.01_A", "W6_B" = "09.01_B", "W7_A" = "09.15_A", "W7_B" = "09.15_B", "W8_A" = "10.01_A", "W8_B" = "10.01_B")
ASV.RA.df.AF.stop <- pivot_longer(stop.times, cols = c('S1_A', 'S2_A', 'S2_B', 'S3_A', 'S3_B', 'S4_A', 'S4_B', 'S5_A', 'S5_B', 'S6_B', 'S7_A', 'S7_B', 'S8_A', 'S8_B', "W1_A", "W2_A" , "W2_B" , "W3_A" , "W3_B" , "W4_A" , "W4_B" , "W5_A" , "W5_B",  "W6_B",  "W7_A"  ,"W7_B"  ,"W8_A",  "W8_B"), names_to = "Sample")  %>% mutate(Sample_Type = case_when(str_detect(Sample, "S") ~ "PM", str_detect(Sample, "W") ~ "Water")) %>% mutate(Site = case_when(str_detect(Sample, "_A") ~ "A", str_detect(Sample, "_B") ~ "B")) %>% mutate_if(is.numeric, ~replace_na(., 0)) %>% mutate(Sampling_Period = case_when(str_detect(Sample, "1") ~ "1", str_detect(Sample, "2") ~ "2", str_detect(Sample, "3") ~ "3", str_detect(Sample, "4") ~ "4", str_detect(Sample, "5") ~ "5", str_detect(Sample, "6") ~ "6", str_detect(Sample, "7") ~ "7", str_detect(Sample, "8") ~ "8")) %>% dplyr::select(-Sample)
ASV.RAs.stop <- ASV.RA.df.AF.stop %>% group_by(ASV, Sampling_Period, Site) %>% pivot_wider(names_from = Sample_Type, values_from = value) 
ASV.AF.without.log.transformation.stop <- ASV.RAs.stop %>% mutate(`r-AF` = PM/Water) %>% mutate(`r-AF`= na_if(`r-AF`, "Inf")) %>% mutate(`r-AF`= na_if(`r-AF`, "NaN")) %>% unite(Sampling_Period, Sampling_Period, Site, remove= TRUE) %>% group_by(ASV) %>% dplyr::summarise_at(vars(`r-AF`), ~mean(., na.rm = T)) %>% left_join(ASV.w.taxonomy, by = "ASV") %>% dplyr::select(ASV, `r-AF`, Species)
ASV.AFs.stop <- ASV.RAs.stop %>% mutate(`r-AF` = (PM)/(Water)) %>% mutate(`r-AF`= na_if(`r-AF`, "NaN")) %>% mutate(`r-AF`= na_if(`r-AF`, "NA")) %>% mutate(`r-AF`= na_if(`r-AF`, "Inf")) %>% mutate_at(vars(`r-AF`), ~replace_na(., 0)) %>% mutate(`1+r-AF` = (1+PM)/(1+Water)) %>% mutate(`log.r-AF` = log(`r-AF`)) %>% mutate(`log.1+r-AF` = log(`1+r-AF`)) %>% mutate(`r-AF`= na_if(`r-AF`, "Inf")) %>% unite(Sampling_Period, Sampling_Period, Site, remove= TRUE) %>% filter(ASV != "ASV4107" , ASV != "ASV2139" , ASV != "ASV1848" , ASV != "ASV0986") 
AF.summary.stats.stop <- ASV.AFs.stop %>% group_by(ASV) %>% get_summary_stats(`r-AF`, type = "mean_sd") %>% left_join( ASV.w.taxonomy, by = "ASV") %>% select(ASV, variable, n, mean, sd, Species)

## Q2 TO DISCUSS WITH STATS CONSULTANT
# To assess whether the different calculated r-AF values are statistically significant on average between species, I wanted to run some sort of variance testing, but then found that my data were too left skewed for ANOVA. So I alternatively chose to run a Kruskal-Wallis test and paired-Wilcox test. (based off googling)
# 1. I found that I have a couple outliers on single days. Some of which are extreme, but I don't have any explanation as to why I should throw them out, so I decided to keep them in is this OK? 
# 2. I don't understand why my Kruskal-Wallis test resulted in a difference between ASVs, but the pairwise test didn't show any two specific ASVS which had statistically significant different aerosolization factors

#ANOVA assessing difference between average r-AF values 
library(tidyverse);library(ggpubr);library(rstatix)
#quick summary stats and one way ANOVA
AF.summary.stats <- ASV.AFs %>% group_by(ASV) %>% get_summary_stats(`r-AF`, type = "mean_sd") %>% left_join( ASV.w.taxonomy, by = "ASV") %>% select(ASV, variable, n, mean, sd, Species)
write.csv(AF.summary.stats, "/Users/haleyplaas/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/files.for.manipulation/r-AF.summary.stats.csv")
ggboxplot(ASV.AFs, x = "ASV", y = "r-AF")
outliers <- ASV.AFs %>% group_by(ASV) %>% identify_outliers(`r-AF`)
#removing outliers from dataset 
ASV.AFs.1 <- anti_join(ASV.AFs, outliers, by = "r-AF")
model <- lm(`r-AF` ~ ASV, data = ASV.AFs)
ggqqplot(residuals(model))
shapiro_test(residuals(model))
#data is non-normal, pivoting to Kruskal Wallis test 
AF.summary.stats <- ASV.AFs %>% group_by(ASV) %>% summarise(
  count = n(),
  mean = mean(`r-AF`, na.rm = TRUE),
  sd = sd(`r-AF`, na.rm = TRUE),
  median = median(`r-AF`, na.rm = TRUE),
  IQR = IQR(`r-AF`, na.rm = TRUE))
AF.summary.stats
write.csv(AF.summary.stats, "/Users/haleyplaas/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/files.for.manipulation/r-AF.summary.stats.csv")
kruskal.test(`r-AF` ~ ASV, data = ASV.AFs)
pairwise <- pairwise.wilcox.test(ASV.AFs$`r-AF`, ASV.AFs$ASV,
                     p.adjust.method = "BH")
pairwise #none of these yielded were individually 

# ---------------------------- LOADING IN ALL WATER AND METEOROLOGICAL METADATA FOR REGRESSION ANALYSES ----------------------------------
# Loading in the water quality metadata
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
metadata.with.sd <- avg.metadata.1 %>% dplyr::left_join(Sd.metadata.1, by = c("Site", "Date"), keep = FALSE) %>% mutate(Date = as.Date(Date,"%Y-%m-%d"))
str(metadata.with.sd)

# Using start times as the water values rather than averages 
metadata.start.date <- metadata.with.sd %>% mutate(Sampling_Period = case_when(Date == '2020-06-11' ~ "S1", 
                                                                               '2020-06-23' == Date ~ "S2",
                                                                               '2020-07-07' == Date ~ "S3", 
                                                                               '2020-07-21' == Date ~ "S4",
                                                                               "2020-08-12" == Date ~ "S5",
                                                                               "2020-08-18" == Date ~ "S6",
                                                                               "2020-09-01" == Date ~ "S7",
                                                                               "2020-09-15" == Date ~ "S8"))

# Using stop times as the water values rather than averages
metadata.stop.date <- metadata.with.sd %>% mutate(Sampling_Period = case_when(Date == '2020-06-23' ~ "S1", 
                                                                               '2020-07-07' == Date ~ "S2",
                                                                               '2020-07-21' == Date ~ "S3", 
                                                                               '2020-08-12' == Date ~ "S4",
                                                                               "2020-08-18" == Date ~ "S5",
                                                                               "2020-09-01" == Date ~ "S6",
                                                                               "2020-09-15" == Date ~ "S7",
                                                                               "2020-10-01" == Date ~ "S8"))

# Converting water metadata to averages over Sampling Periods 
write.csv(metadata.with.sd, "/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/all.water.metadata.csv") # this data frame was exported to manipulate externally in Excel. To approximate the water conditions during the course of each aerosol sampling, water metadata collected at the start and end of each aerosol sampling period was averaged. 
Sampling_Period_water_metadata.with.sd <- read.csv("/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/all.water.metadata.2.csv", header = TRUE, na.strings = c(""," ", ".", "NA","#DIV/0!"))
Sampling_Period_water_metadata <- Sampling_Period_water_metadata.with.sd %>% select(1:29)

# loading in the atmosphere metadata
# Karsten's cleaned Met station data for wind vectors
wind.data <- read.csv("/Users/haleyplaas/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/DaytimeMsmts_kb.csv", header = TRUE, na.strings = c(""," ", ".", "NA", "NaN", "NAN")) 
wind.data.1 <- separate(wind.data, col = Date_Time, into = c("Date", "Time"), sep = " ", remove = T)
wind.data.2 <- wind.data.1 %>% mutate(NS.vec = WS*cos(WD)) %>% mutate(EW.vec = WS*sin(WD)) %>% mutate(Date = as.Date(Date,"%m/%d/%y"))
wind.data.3 <- wind.data.2 %>% mutate(Sampling_Period = case_when('2020-06-12' <= Date & Date <= '2020-06-23' ~ "S1", 
                                                                  '2020-06-24' <= Date & Date <= '2020-07-07' ~ "S2",
                                                                  '2020-07-08' <= Date & Date <= '2020-07-21' ~ "S3",
                                                                  '2020-07-22' <= Date & Date <= '2020-08-03' ~ "S4",
                                                                  "2020-08-13" <= Date & Date <= '2020-08-18' ~ "S5",
                                                                  "2020-08-19" <= Date & Date <= '2020-09-01' ~ "S6",
                                                                  "2020-09-02" <= Date & Date <= '2020-09-14' ~ "S7",
                                                                  "2020-09-15" <= Date & Date <= '2020-10-01' ~  "S8")) %>% select(-Time)
wind.data.avg.Sampling.Period <- wind.data.3 %>% group_by(Sampling_Period) %>% dplyr::summarise_all((~mean(na.omit(.)))) %>% mutate(Air.Mass = atan(EW.vec/NS.vec)) %>% mutate(Site = "A")
wind.data.avg.Sampling.Period.B <- wind.data.avg.Sampling.Period %>% mutate(Site = "B") 
wind.data.avg.Sampling.Period.1 <- rbind(wind.data.avg.Sampling.Period, wind.data.avg.Sampling.Period.B) %>% dplyr::rename(Wind.Speed = WS) %>% dplyr::rename(Wind.Direction = WD)

# Loading in the other meteorological station data
met.station.data <- read.csv("/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/met.station.data.csv", header = TRUE, na.strings = c(""," ", ".", "NA", "NaN")) 
met.station.data.1 <- met.station.data  %>% mutate(Date = as.Date(Date, "%Y-%m-%d")) 
met.station.data.2 <- met.station.data.1 %>% group_by(Sampling_Period, Site) %>% dplyr::summarise_all((~mean(na.omit(.)))) %>% select(-Date)
met.station.data.2[is.na(met.station.data.2)] <- NA

# combining metadata
# AVERAGED
all.metadata <- Sampling_Period_water_metadata %>% dplyr::left_join(met.station.data.2, by = c("Site", "Sampling_Period"), keep = F) %>% left_join(wind.data.avg.Sampling.Period.1, by = c("Site", "Sampling_Period"), keep = F)

# adding the relative abundances of bacteria and cyanobacteria in PM and water samples to the metadata 
class.RAs$Sampling_Period <- sub("^", "S", class.RAs$Sampling_Period) #don't run twice
class.RAs.2 <- class.RAs %>% filter(Class == "Cyanobacteria") %>% pivot_wider(names_from = Class, values_from = c(PM, Water))
species.RAs$Sampling_Period <- sub("^", "S", species.RAs$Sampling_Period) #don't run twice
species.RAs.2 <- species.RAs %>% pivot_wider(names_from = Species, values_from = c(PM, Water)) 
relative.abundances <- class.RAs.2 %>% left_join(species.RAs.2, by = c("Site", "Sampling_Period"))
non.zeroes <- relative.abundances %>% ungroup(Sampling_Period, Site) %>% select(-Sampling_Period, -Site) %>% select_if(colSums(.) != 0)
relative.abundances.1 <- relative.abundances %>% dplyr::select(Site, Sampling_Period) %>% cbind(non.zeroes) %>% group_by(Site, Sampling_Period)
regression.df.0 <- all.metadata %>% dplyr::left_join(relative.abundances.1, by = c("Site", "Sampling_Period"))

# Loading in individual PM2.5 mass concentrations readings from pDR data to add to metadata 
col.names <- c("Site", "record", "ug/m3", "Temp", "RHumidity", "AtmoPressure", "Flags", "Time", "Date")
june.23.a <- read.delim("/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/pDR/2020-06-23_pDRA.txt", header = FALSE, sep = ",", dec = ".", skip = 25, na.strings = c("NA"))
A <- as.data.frame(rep("A", times = 12529, length.out = NA, each = 1)) 
colnames(A) <- "Site"
june.23.a.1 <- cbind(A, june.23.a)
colnames(june.23.a.1) <- col.names 
june.23.a.1 <- june.23.a.1 %>% dplyr::select(-record)
june.23.b <- read.delim("/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/pDR/2020-06-23_pDRB.txt", header = FALSE, sep = ",", dec = ".", skip = 25, na.strings = c("NA"))
B <- as.data.frame(rep("B", times = 16964, length.out = NA, each = 1)) 
colnames(A) <- "Site"
june.23.b.1 <- cbind(B, june.23.b)
colnames(june.23.b.1) <- col.names 
june.23.b.1 <- june.23.b.1 %>% dplyr::select(-record)
july.07.a.1 <- read.delim("/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/pDR/2020-07-07_pDRA.1.txt", header = FALSE, sep = ",", dec = ".", skip = 25, na.strings = c("NA"))
july.07.a.2 <- read.delim("/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/pDR/2020-07-07_pDRA.2.txt", header = FALSE, sep = ",", dec = ".", skip = 25, na.strings = c("NA"))
july.07.a <- rbind(july.07.a.1, july.07.a.2)
A.1 <- as.data.frame(rep("A", times = 16374, length.out = NA, each = 1)) 
colnames(A.1) <- "Site"
july.07.a.3 <- cbind(A.1, july.07.a)
colnames(july.07.a.3) <- col.names 
july.07.a.3 <- july.07.a.3 %>% dplyr::select(-record)
july.07.b.1 <- read.delim("/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/pDR/2020-07-07_pDRB.1.txt", header = FALSE, sep = ",", dec = ".", skip = 25, na.strings = c("NA"))
july.07.b.2 <- read.delim("/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/pDR/2020-07-07_pDRB.2.txt", header = FALSE, sep = ",", dec = ".", skip = 25, na.strings = c("NA"))
july.07.b <- rbind(july.07.b.1, july.07.b.2)
B.1 <- as.data.frame(rep("A", times = 19648, length.out = NA, each = 1)) 
colnames(B.1) <- "Site"
july.07.b.3 <- cbind(B.1, july.07.b)
colnames(july.07.b.3) <- col.names 
july.07.b.3 <- july.07.b.3 %>% dplyr::select(-record)
july.21.a <- read.delim("/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/pDR/2020-07-21_pDRA.txt", header = FALSE, sep = ",", dec = ".", skip = 25, na.strings = c("NA"))
A.2 <- as.data.frame(rep("A", times = 20098, length.out = NA, each = 1)) 
colnames(A.2) <- "Site"
july.21.a.1 <- cbind(A.2, july.21.a)
colnames(july.21.a.1) <- col.names 
july.21.a.1 <- july.21.a.1 %>% dplyr::select(-record)
july.21.b <- read.delim("/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/pDR/2020-07-21_pDRB.txt", header = FALSE, sep = ",", dec = ".", skip = 25, na.strings = c("NA"))
B.2 <- as.data.frame(rep("B", times = 20041, length.out = NA, each = 1)) 
colnames(B.2) <- "Site"
july.21.b.1 <- cbind(B.2, july.21.b)
colnames(july.21.b.1) <- col.names 
july.21.b.1 <- july.21.b.1 %>% dplyr::select(-record)
aug.03.a <- read.delim("/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/pDR/2020-08-03_pDRA.txt", header = FALSE, sep = ",", dec = ".", skip = 25, na.strings = c("NA"))
A.3 <- as.data.frame(rep("A", times = 3426, length.out = NA, each = 1)) 
colnames(A.3) <- "Site"
aug.03.a.1 <- cbind(A.3, aug.03.a)
colnames(aug.03.a.1) <- col.names 
aug.03.a.1 <- aug.03.a.1 %>% dplyr::select(-record)
aug.03.b <- read.delim("/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/pDR/2020-08-03_pDRB.txt", header = FALSE, sep = ",", dec = ".", skip = 25, na.strings = c("NA"))
B.3 <- as.data.frame(rep("B", times = 18537, length.out = NA, each = 1)) 
colnames(B.3) <- "Site"
aug.03.b.1 <- cbind(B.3, aug.03.b)
colnames(aug.03.b.1) <- col.names 
aug.03.b.1 <- aug.03.b.1 %>% dplyr::select(-record)
aug.18.b <- read.delim("/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/pDR/2020-08-18_pDRB.txt", header = FALSE, sep = ",", dec = ".", skip = 25, na.strings = c("NA"))
B.4 <- as.data.frame(rep("B", times = 8553, length.out = NA, each = 1)) 
colnames(B.4) <- "Site"
aug.18.b.1 <- cbind(B.4, aug.18.b)
colnames(aug.18.b.1) <- col.names 
aug.18.b.1 <- aug.18.b.1 %>% dplyr::select(-record)
sept.01.b <- read.delim("/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/pDR/2020-09-01_pDRB.txt", header = FALSE, sep = ",", dec = ".", skip = 25, na.strings = c("NA"))
B.5 <- as.data.frame(rep("A", times = 34186, length.out = NA, each = 1)) 
colnames(B.5) <- "Site"
sept.01.b.3 <- cbind(B.5, sept.01.b)
colnames(sept.01.b.3) <- col.names 
sept.01.b.3 <- sept.01.b.3 %>% dplyr::select(-record)
sept.15.b <- read.delim("/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/pDR/2020-09-15_pDRB.txt", header = FALSE, sep = ",", dec = ".", skip = 25, na.strings = c("NA"))
B.6 <- as.data.frame(rep("B", times = 20058, length.out = NA, each = 1)) 
colnames(B.6) <- "Site"
sept.15.b.1 <- cbind(B.6, sept.15.b)
colnames(sept.15.b.1) <- col.names 
sept.15.b.1 <- sept.15.b.1 %>% dplyr::select(-record)
oct.01.b <- read.delim("/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/pDR/2020-10-01_pDRB.txt", header = FALSE, sep = ",", dec = ".", skip = 25, na.strings = c("NA"))
B.7 <- as.data.frame(rep("B", times = 23171, length.out = NA, each = 1)) 
colnames(B.7) <- "Site"
oct.01.b.1 <- cbind(B.7, oct.01.b)
colnames(oct.01.b.1) <- col.names 
oct.01.b.1 <- oct.01.b.1 %>% dplyr::select(-record)

# Creating a Time-series of ambient PM2.5 
# cleaning the pDR data 
all.pdr.data <- rbind(june.23.a.1, june.23.b.1, july.07.a.3, july.07.b.3, july.21.a.1, july.21.b.1, aug.03.a.1, aug.03.b.1, aug.18.b.1, sept.01.b.3, sept.15.b.1, oct.01.b.1) 
all.pdr.data <- all.pdr.data %>% dplyr::mutate(Site = as.character(Site), `ug/m3` = as.numeric(`ug/m3`), Temp = as.numeric(Temp), RHumidity = as.numeric(RHumidity), AtmoPressure = as.numeric(AtmoPressure), Time = as.character(Time, "%H%M%S"), Date = as.Date(Date, "%d-%b-%Y")) %>% dplyr::select(-Flags)

# 1. selecting daytime readings only 
all.pdr.data.1 <- all.pdr.data %>% mutate(Time.1 = as.numeric(gsub("[[:punct:]]", "", Time)))
daytime <- all.pdr.data.1 %>% dplyr::filter(Time.1 > 70000 & Time.1 < 190000) 

## Q3. TO DISCUSS WITH STATS CONSULTANT. When is is appropriate to remove outlier readings or entire days that are outliers from a continous dataset? 
# Here, I removed entire days when we had known instances where the readings were extremely high or faulty (found to be extreme outliers), but I only have concrete reasoning for Fourth of July, and 7-11 through 7-20 when readings went beserk and then we found a spider blocking the sensor, and then other days that were found to be "extreme" outliers using the identify_outliers() function.  

# 2. removal of days with known external sources of PM (e.g. fireworks, military activity, spider found in equipment)
cleaned <- daytime %>% dplyr::filter(Date != "2020-07-04", Date != "2020-07-05", Date != "2020-07-06", Date != "2020-06-28", Date != "2020-09-03") 
to.remove <- cleaned %>% dplyr::filter(Site == "B" & Date == "2020-07-11" |
                                         Site == "B" & Date == "2020-07-12"|
                                         Site == "B" & Date == "2020-07-13"|
                                         Site == "B" & Date == "2020-07-14"|
                                         Site == "B" & Date == "2020-07-15"|
                                         Site == "B" & Date == "2020-07-16"|
                                         Site == "B" & Date == "2020-07-17"|
                                         Site == "B" & Date == "2020-07-18"|
                                         Site == "B" & Date == "2020-07-19"|
                                         Site == "B" & Date == "2020-07-20")
cleaned.2 <- anti_join(cleaned, to.remove)

## Q4. TO DISCUSS WITH STATS CONSULTANT
# Here, I removed single readings which were found to be outliers based on calculations of z-scores. Is this OK considering the spread of the data? I plotted the data to observe a normally distribution in order for this, band they were only slightly left skewed, so is there a better way to identify outliers? I note that I am also inconsistent with use of z-scores and the identify_outliers() function because of pasting random code from the internet in attempt to ID outlier data...   

# 3. Outliers calculated and removed via z-scores (z > 3)
all.pdr.data.2 <- daytime %>% mutate(date.time = as.POSIXct(paste(Date, Time), format="%Y-%m-%d %H:%M:%S"))
mass.conc.2 <- all.pdr.data.2 %>% dplyr::select(Site,`ug/m3`, date.time)
only.mass.conc.2 <- mass.conc.2 %>% dplyr::select(`ug/m3`)
z_scores.2 <- mass.conc.2 %>% mutate(z.score = sapply(only.mass.conc.2, function(only.mass.conc.2) (abs(only.mass.conc.2-mean(only.mass.conc.2))/sd(only.mass.conc.2)))) 
no_outliers.2 <- z_scores.2 %>% dplyr::filter(z.score < 3) 
plot.3 <- ggplot(no_outliers.2, aes(x = date.time, y = `ug/m3`, color = Site)) + geom_point(size=0.01) + ylim(c(-1,100))
plot.3 
hist(cleaned.2$`ug/m3`, breaks = seq(from=-100, to=150, by=5)) #normal distribution = fair to use z-scores for outlier identification 
mass.conc <- cleaned.2 %>% dplyr::select(Site,`ug/m3`, Time, Date)
only.mass.conc <- mass.conc %>% dplyr::select(`ug/m3`)
z_scores <- mass.conc %>% mutate(z.score = sapply(only.mass.conc, function(only.mass.conc) (abs(only.mass.conc-mean(only.mass.conc))/sd(only.mass.conc)))) 
no_outliers <- z_scores %>% dplyr::filter(z.score < 3) 

# Collating PM data from both sites due to gaps in data from ambient PM samplers going offline and being inconsistent
# Demonstrating the correlation between readings at both sites to justify integration of data from each sites
correlation <- dplyr::select(no_outliers, Date, Site,`ug/m3`)
correlation.1 <- correlation %>% group_by(Date, Site) %>% summarise_at(vars(c("ug/m3")), ~mean(na.omit(.)))
correlation.2 <- pivot_wider(correlation.1, names_from = Site, values_from = `ug/m3`) %>% na.omit() %>% dplyr::filter(Date < "2020-07-01")
correlation.3 <- lm(correlation.2$A ~ correlation.2$B , data = correlation.2)
summary(correlation.3)
ggplot(correlation.2, aes(x = A, y = B)) + 
  geom_point() + 
  stat_smooth(method = "lm", col = "black")
# Combining readings from both Sites
pdr.together <- no_outliers %>% group_by(Date) %>% summarise_at(vars(c("ug/m3")), ~mean(na.omit(.)))
pdr.together %>% identify_outliers(`ug/m3`) #confirmation of no more outliers
pdr.dates <- as.data.frame(seq(as.Date("2020-06-11"), as.Date("2020-10-01"), by="days"))
name <- "Date"
colnames(pdr.dates) = name
# need to write code here to fix gaps in data on the figure which shows the time series for PM2.5 (for right now I just modified the graph in powerpoint post-plotting and export from R)

# Adding Chlorophyll a and cyanobacterial relative abundance data to visualize alongside PM time series 
Date <- c("2020-06-17", "2020-06-30", "2020-06-30", "2020-07-14","2020-07-14","2020-07-28","2020-07-28","2020-08-15","2020-08-15","2020-08-25","2020-09-08","2020-09-08","2020-09-22","2020-09-22") 
dates <- as.data.frame(Date) %>% mutate(Date = as.Date(Date, "%Y-%m-%d"))
cyano.RA <- regression.df.0 %>% dplyr::select(Site, Sampling_Period, CHLA_avg, Water_Cyanobacteria, Water_Dolichospermum_NIES41) %>% mutate(Water_Cyanobacteria = Water_Cyanobacteria*10^2.5, Water_Dolichospermum_NIES41 = Water_Dolichospermum_NIES41*10^2) %>% mutate(dates) %>% dplyr::select(-Sampling_Period) 
chla.with.PM <- pdr.together %>% left_join(cyano.RA, by = "Date", keep = FALSE) 
chla.with.PM.1 <- chla.with.PM %>% mutate_if(is.numeric, ~replace_na(., NA))
for.avg <- chla.with.PM.1 %>% select(`ug/m3`) %>% na.omit()
average.PM <- as.data.frame(5.81, times = 105)
name <- "avg"
colnames(average.PM) = name
chla.with.PM.2 <- cbind(chla.with.PM.1, average.PM) 
EPA.std <- as.data.frame(10.5, times = 105)
name <- "EPA.std"
colnames(EPA.std) = name
chla.with.PM.3 <- cbind(chla.with.PM.2, EPA.std) 

# plot by sites 
ggplot(chla.with.PM.3, aes(x = Date)) + 
  geom_line(size=0.5, aes(y = `ug/m3`), na.rm=F) + 
  geom_point(size = 3, aes(y=CHLA_avg, color = Site)) +
  geom_line(linetype = 3, aes(y=avg)) + 
  geom_line(linetype = 5, aes(y=EPA.std)) + 
  theme_minimal() + 
  scale_color_grey() + 
  scale_x_date(date_labels="%B") + 
  scale_y_continuous(sec.axis = sec_axis(trans= ~.*1, name = expression(µg~L^-1~Chlorophyll~a))) +
  xlab("Date") + ylab(expression(µg~m^-3~PM[2.5])) +
  theme(plot.title = element_text(face = "bold"), legend.position = "none")

# pooling chlorophyll a data from both sites so it is more representative of entire river between the two sites
cyano.RA.no.sites <- cyano.RA %>% group_by(Date) %>% dplyr::summarise_at(vars(c("CHLA_avg", "Water_Cyanobacteria", "Water_Dolichospermum_NIES41")), ~mean(na.omit(.)))
both.sites.1 <- pdr.together %>%  left_join(cyano.RA.no.sites, by = "Date", keep = FALSE)
both.sites.2 <- both.sites.1 %>% mutate_if(is.numeric, ~replace_na(., NA))
for.avg <- both.sites.1 %>% select(`ug/m3`) %>% na.omit()
average.PM <- as.data.frame(5.81, times = 99)
name <- "avg"
colnames(average.PM) = name
both.sites.3 <- cbind(both.sites.2, average.PM) 
EPA.std <- as.data.frame(10.5, times = 99)
name <- "EPA.std"
colnames(EPA.std) = name
both.sites.4 <- cbind(both.sites.3, EPA.std) 
ggplot(both.sites.4, aes(x = Date)) + 
  geom_line(size=0.5, aes(y = `ug/m3`), na.rm=F) + 
  geom_point(size = 3, aes(y=CHLA_avg)) + 
  geom_point(size = 3, aes(y=Water_Dolichospermum_NIES41), color = "green") + 
  geom_line(linetype = 5, color = "gray", aes(y=avg)) + 
  geom_line(linetype = 5, color = "gray50",aes(y=EPA.std)) + 
  theme_minimal() + 
  scale_color_grey() + 
  scale_x_date(date_labels="%B") + 
  scale_y_continuous(sec.axis = sec_axis(trans= ~.*1, name = expression(µg~L^-1~Chlorophyll~a))) +
  xlab("Date") + ylab(expression(µg~m^-3~PM[2.5])) +
  ylim(c(0,25)) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

# Determining days that exceeded the baseline ambient PM2.5 as a function of bloom conditions 
average.PM <- mean(both.sites.3$`ug/m3`)
median.PM <- median(both.sites.3$`ug/m3`)
pdr.together.1 <- pdr.together %>% mutate(exceeds.baseline = case_when(
  pdr.together$`ug/m3` >= median.PM ~ "1",
  pdr.together$`ug/m3` < median.PM ~ "0")) %>% 
  mutate(exceeds.EPA.standard = case_when(
    pdr.together$`ug/m3` > 12 ~ "1",
    pdr.together$`ug/m3` < 12 ~ "0")) %>% 
  mutate(Bloom = case_when(
    pdr.together$Date >= "2020-06-23" & pdr.together$Date <= "2020-07-21"~ "Bloom",
    pdr.together$Date < "2020-06-23" | pdr.together$Date > "2020-07-21"~ "No Bloom"))
myDate = as.POSIXct(pdr.together.1$Date)
pdr.together.2 <- pdr.together.1 %>% mutate(month = format(myDate,"%b")) #for some reason it's one off. (the first of every month is classified as the month prior--fix this later to save time)
table(pdr.together.2$month)
Jun <- rep("Jun", times = 19)
Jul <- rep("Jul", times = 28)
Aug <- rep("Aug", times = 22)
Sep <- rep("Sep", times = 29)
Oct <- rep("Oct", times = 1)
true.months <- c(Jun,Jul,Aug,Sep,Oct)
pdr.together.3 <- pdr.together.2 %>% mutate(month = true.months)

## TO DISCUSS WITH STATS CONSULTANT
# Q5. Here I used a Wilcoxon test to determine whether bloom vs. non-bloom PM2.5 data were different, but I want to make sure I did this correctly.

# Statistical testing to compare the PM2.5 concentrations on bloom vs. non-bloom days.
pdr.together.test <- pdr.together.3 %>% na.omit()
pdr.together.test %>% shapiro_test(ug/m3)
hist(pdr.together.test$`ug/m3`) #non-normal, use Wilcoxon Signed Ranks Test, the nonparametric alternative to a T test
wilcox.test <- wilcox.test(`ug/m3` ~ Bloom, alternative = "greater", data = pdr.together.test)
wilcox.summary.stats <- pdr.together.test %>% group_by(Bloom) %>% summarise(count = n(), median = median(`ug/m3`), IQR = IQR(`ug/m3`))
wilcox.test
wilcox.summary.stats
ggplot(data = pdr.together.test, aes(x = Bloom, y = `ug/m3`, color=Bloom)) + geom_boxplot() + geom_jitter(aes(color=Bloom, alpha=Bloom)) + theme_bw() + scale_color_manual(values=c("gray50", "gray")) + xlab(" ") + ylab(expression(µg~m^-3~PM[2.5])) + scale_alpha_manual(values=c(1,0.4))

## ------------------------------- Univariate Regressions to examine environmental drivers of aquatic Cyanobacteria in PM -----------------------
# linear regressions
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

#START TIMES ONLY
start.time.metadata <- metadata.start.date %>% dplyr::left_join(met.station.data.2, by = c("Site", "Sampling_Period"), keep = F) %>% left_join(wind.data.avg.Sampling.Period.1, by = c("Site", "Sampling_Period"), keep = F) %>% tidyr::drop_na(Sampling_Period)
start.regression.df.0 <- start.time.metadata %>% dplyr::left_join(relative.abundances.1, by = c("Site", "Sampling_Period"))
#regression.df <- start.regression.df.0 %>% dplyr::select(-Site, -Sampling_Period) %>% mutate_if(is.integer, as.numeric) #only run this when using start times

#STOP TIMES ONLY
stop.time.metadata <- metadata.stop.date %>% dplyr::left_join(met.station.data.2, by = c("Site", "Sampling_Period"), keep = F) %>% left_join(wind.data.avg.Sampling.Period.1, by = c("Site", "Sampling_Period"), keep = F) %>% tidyr::drop_na(Sampling_Period)
stop.regression.df.0 <- stop.time.metadata %>% dplyr::left_join(relative.abundances.1, by = c("Site", "Sampling_Period"))
#regression.df <- stop.regression.df.0 %>% dplyr::select(-Site, -Sampling_Period) %>% mutate_if(is.integer, as.numeric) #only run this when using stop times

## Q6. TO DISCUSS WITH STATS CONSULTANT WHAT ARE THE MOST IMPORTANT VARIABLES TO SHOW FROM A LINEAR REGRESSION OUTPUT AND WHY? I chose adjusted R^2, sigma, F statistic, p value, AIC, and df residuals. Chose univariate regressions due to lack of data points (not enough power for multivariate? right?), plus only interested in association between variables, not using it to predict outcomes. 

# Outcome Variable: Cyanobacteria in PM
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
cyano.1 <- linear.regression.plot(regression.df$PM_Cyanobacteria, regression.df$PM_Pseudanabaena_PCC7429) + labs(x = "Pseudanabaena in PM", y = "Cyanobacteria in PM") 
cyano.2 <- linear.regression.plot(regression.df$PM_Cyanobacteria, regression.df$PM_Anabaena) + labs(x = "Anabaena in PM", y = "Cyanobacteria in PM")
cyano.3 <- linear.regression.plot(regression.df$PM_Cyanobacteria, regression.df$PM_Cyanobium_PCC6307) + labs(x = "Cyanobium in PM", y = "Cyanobacteria in PM") 
cyano.4 <- linear.regression.plot(regression.df$PM_Cyanobacteria, regression.df$PM_Gloeocapsa) + labs(x = "Gloeocapsa in PM", y = "Cyanobacteria in PM") 
cyano.1 + cyano.2 + cyano.3 + cyano.4

# Outcome Variable: Cyanobacteria in water
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

# Outcome Variable: Microcystis in PM
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
micro.1 <- linear.regression.plot(regression.df$PM_Microcystis_PCC7914, regression.df$PM_Dolichospermum_NIES41) + labs(x = "Dolichospermum in PM", y = "Microcystis in PM")
micro.2 <- linear.regression.plot(regression.df$PM_Microcystis_PCC7914, regression.df$PM_Cyanobium_PCC6307) + labs(x = "Cyanobium in PM", y = "Microcystis in PM")
micro.3 <- linear.regression.plot(regression.df$PM_Microcystis_PCC7914, regression.df$Water_Haloleptolyngbya)+ labs(x = "Haloleptolyngbya in Water", y = "Microcystis in PM")
micro.4 <- linear.regression.plot(regression.df$PM_Microcystis_PCC7914, regression.df$Water_Gloeocapsa) + labs(x = "Gloeocapsa in Water", y = "Microcystis in PM")
micro.5 <- linear.regression.plot(regression.df$PM_Microcystis_PCC7914, regression.df$NH4_avg) + labs(x = "NH4 concentration", y = "Microcystis in PM")

micro.1 + micro.2 + micro.3 + micro.4 + micro.5 

# Outcome Variable: Microcystis in water
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

# Outcome Variable: Dolichospermum in PM
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
dolicho.1 <- linear.regression.plot(regression.df$PM_Dolichospermum_NIES41, regression.df$PM_Microcystis_PCC7914) + labs(x = "Microcystis in PM", y = "Dolichospermum in PM")
dolicho.2 <- linear.regression.plot(regression.df$PM_Dolichospermum_NIES41, regression.df$PM_Aphanizomenon_NIES81) + labs(x = "Aphanizomenon in PM", y = "Dolichospermum sp. in PM")
dolicho.3 <- linear.regression.plot(regression.df$PM_Dolichospermum_NIES41, regression.df$C.to.N_avg) + labs(x = "C:N ratio", y = "Dolichospermum in PM")
dolicho.4 <- linear.regression.plot(regression.df$PM_Dolichospermum_NIES41, regression.df$Relative.Humidity) + labs(x = "Relative Humidity", y = "Dolichospermum in PM")

dolicho.1 + dolicho.2 + dolicho.3 + dolicho.4 

#Outcome Variable: Dolichospermum in Water
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
PM.1 <- linear.regression.plot(regression.df$PM.avg, regression.df$Air.Temp) + labs(x = "Air Temperature", y = "PM mass concentration") + ylim(c(0,16))
PM.2 <- linear.regression.plot(regression.df$PM.avg, regression.df$Solar.Rad) + labs(x = "Solar Irradiance", y = "PM mass concentration")+ ylim(c(0,16))
PM.3 <- linear.regression.plot(regression.df$PM.avg, regression.df$Salinity_avg ) + labs(x = "Salinity", y = "PM mass concentration")+ ylim(c(0,16))
PM.4 <- linear.regression.plot(regression.df$PM.avg, regression.df$Relative.Humidity ) + labs(x = "Relative Humidity", y = "PM mass concentration")+ ylim(c(0,16))
PM.5 <- linear.regression.plot(regression.df$PM.avg, regression.df$PM_z.Not.Assigned ) + labs(x = "Unassigned Cyanobacteria in PM", y = "PM mass concentration")+ ylim(c(0,16))
PM.6 <- linear.regression.plot(regression.df$PM.avg, regression.df$Water_Pseudanabaena_PCC7429) + labs(x = "Pseudanabaena in Water", y = "PM mass concentration")+ ylim(c(0,16))
PM.7 <- linear.regression.plot(regression.df$PM.avg, regression.df$C.to.N_avg) + labs(x = "C:N Ratio", y = "Ambient PM concentration")+ ylim(c(0,16))
PM.8 <- linear.regression.plot(regression.df$PM.avg, regression.df$Water_Anabaena) + labs(x = "Anabaena in Water", y = "PM mass concentration")+ ylim(c(0,16))
PM.9 <- linear.regression.plot(regression.df$PM.avg, regression.df$Water_Dolichospermum_NIES41) + labs(x = "Dolichospermum in Water", y = "PM mass concentration")+ ylim(c(0,16))
PM.10 <- linear.regression.plot(regression.df$PM.avg, regression.df$Water_Caenarcaniphilales_XX ) + labs(x = "Caenarcaniphilales in Water", y = "PM mass concentration")+ ylim(c(0,16))
PM.11 <- linear.regression.plot(regression.df$PM.avg, regression.df$Water_Nostoc_PCC7524) + labs(x = "Nostoc in Water", y = "PM mass concentration")+ ylim(c(0,16))
PM.12 <- linear.regression.plot(regression.df$PM.avg, regression.df$PM_Caldora_VP642b) + labs(x = "Caldora in PM", y = "PM mass concentration")+ ylim(c(0,16))
PM.13 <- linear.regression.plot(regression.df$PM.avg, regression.df$Water_Caenarcaniphilales_XX) + labs(x = "Caenarcaniphilales in Water", y = "PM mass concentration")+ ylim(c(0,16))

PM.1 + PM.2 + PM.3 + PM.4 + PM.5 + PM.6 + PM.7 + PM.8 + PM.9 + PM.10 + PM.11 + PM.12 + PM.13 

#Key cyanos in air and water 
Cyano.plot <- linear.regression.plot(regression.df$Water_Cyanobacteria, regression.df$PM_Cyanobacteria) + labs(x = "Water", y = "PM") + ggtitle("Cyanobacteria")
Dolicho.plot <- linear.regression.plot(regression.df$Water_Dolichospermum_NIES41, regression.df$PM_Dolichospermum_NIES41) + labs(x = "Water", y = "PM") + ggtitle("Dolichospermum")
Micro.plot <- linear.regression.plot(regression.df$Water_Microcystis_PCC7914, regression.df$PM_Microcystis_PCC7914) + labs(x = "Water", y = "PM") + ggtitle("Microcystis")
Aphani.plot <- linear.regression.plot(regression.df$Water_Aphanizomenon_NIES81, regression.df$PM_Aphanizomenon_NIES81) + labs(x = "Water", y = "PM") + ggtitle("Aphanizomenon")
Caenar.plot <- linear.regression.plot(regression.df$Water_Caenarcaniphilales_XX, regression.df$PM_Caenarcaniphilales_XX) + labs(x = "Water", y = "PM") + ggtitle("Caenarcaniphilales")

Cyano.plot + Aphani.plot + Dolicho.plot + Micro.plot
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

#showing committee results from using water values as averages between two dates, then only start and only stop dates WRT aerosol sampling
#Exporting Statistically Significant Environmental Drivers from each aerosol linear regression model series of interest
write.csv(significant.predictors.cyanos, "/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/files.for.manipulation/avg.PM.cyanos.csv")
write.csv(significant.predictors.PM_Dolichospermum_NIES41, "/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/files.for.manipulation/avg.PM.dolicho.csv")
write.csv(significant.predictors.PM_Microcystis_PCC7914, "/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/files.for.manipulation/avg.PM.micro.csv")
write.csv(significant.predictors.PM.avg, "/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/files.for.manipulation/avg.PM.csv")

#start only
write.csv(significant.predictors.cyanos, "/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/files.for.manipulation/start.PM.cyanos.csv")
write.csv(significant.predictors.PM_Dolichospermum_NIES41, "/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/files.for.manipulation/start.PM.dolicho.csv")
write.csv(significant.predictors.PM_Microcystis_PCC7914, "/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/files.for.manipulation/start.PM.micro.csv")
write.csv(significant.predictors.PM.avg, "/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/files.for.manipulation/start.PM.csv")

#stop only
write.csv(significant.predictors.cyanos, "/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/files.for.manipulation/stop.PM.cyanos.csv")
write.csv(significant.predictors.PM_Dolichospermum_NIES41, "/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/files.for.manipulation/stop.PM.dolicho.csv")
write.csv(significant.predictors.PM_Microcystis_PCC7914, "/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/files.for.manipulation/stop.PM.micro.csv")
write.csv(significant.predictors.PM.avg, "/Users/haleyplaas/OneDrive - University of North Carolina at Chapel Hill/Coding/R/Chowan Data/CHHEPilotStudy/Chowan_ASVs/files.for.manipulation/stop.PM.csv")
