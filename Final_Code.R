setwd("~/Desktop/R_files")
#### Decontaminate 16S rRNA samples #####
library(phyloseq)
library(ggplot2)
library(decontam)
###creating phyloseq object###
# set working directory 
# read in OTU table 
otu <- read.table(file = "otu_table.txt", header = TRUE)
# read in taxonomy table
tax <- read.table(file = "taxonomy.tsv", sep = '\t', header = TRUE)
# merge files 
merged_file <- merge(otu, tax, by.x = c("OTUID"), by.y=c("OTUID"))
# note: number of rows should equal your shortest file length, drops taxonomy for OTUs that don’t exist in your OTU table

# output merged .txt file
write.table(merged_file, file = "combined_otu_tax", sep = '\t', col.names = TRUE, row.names = FALSE)

# It seems tedious but you need to open the merged .txt file in excel and split into two files: one for taxonomy (containing only the columns OTUID and taxonomic info) and the other for the OTU matrix (containing only OTUID and abundances in each sample). Note: for the taxonomy file, you need to use data —> text-to-columns in Excel and separate on semicolon to get columns for kingdom, phylum, class, etc… once you make these two separate files in excel, save each as a .csv 

# read in otu table
otu_table = read.csv("otu_matrix.csv", sep=",", row.names=1)
otu_table = as.matrix(otu_table)

# read in taxonomy
# seperated by kingdom phylum class order family genus species
taxonomy = read.csv("taxonomy_revised.csv", sep=",", row.names=1)
taxonomy = as.matrix(taxonomy)

# read in metadata 
#SampleID	BarcodeSequence	LinkerPrimerSequence	BarcodePlate	Well	Description	Body_Type	Caste_Type	Species
metadata = read.table("sample-metadata_merged_qpcr.txt", row.names=1, header = TRUE)
# read in tree
phy_tree = read_tree("unrooted_tree.nwk")

# import as phyloseq objects
OTU = otu_table(otu_table, taxa_are_rows = TRUE)
TAX = tax_table(taxonomy)
META = sample_data(metadata)

# (tree was already imported as a phyloseq object)


# merge into one phyloseq object
physeq_con<- phyloseq(OTU, TAX, META, phy_tree)
ps<- physeq_con
head(sample_data(ps))

###Decontam###
df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()


#prevelance method
sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))

ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)

contamcomb<- isContaminant(ps, conc="quant_reading", neg = "is.neg", method = c("combined"),
                           batch = NULL, batch.combine = c("fisher"),
                           threshold = 0.1, normalize = TRUE, detailed = TRUE)
table(contamcomb$contaminant)
head(which(contamcomb$contaminant))
combo <- subset(contamcomb, contaminant=="TRUE")
combo
write.csv(combo, file= "combo.csv")

#remove contaminants
ps.noncontam <- prune_taxa(!contamcomb$contaminant, ps)
ps.noncontam

#remove control samples from dataset
Phylo_decontam = subset_samples(ps.noncontam, Description != "negative1" & Description != "negative2" & Description != "NEGATIVEEL1" & Description != "NEGATIVEEL2")
Phylo_decontam


# rarefying data
library(vegan)
library(ranacapa)
set.seed(100)
Phylo_decontam_nopro <- subset_samples(Phylo_decontam, Species != "Procryptocerus_sp")
rarefied_nopro <- custom_rarefaction(Phylo_decontam_nopro, sample_size =13814, replicates = 100, rngseed = TRUE)
physeq_rarefied_nopro_qiime_rooted <- rarefied_nopro

set.seed(100)
rarefied_withpro_qiime <- custom_rarefaction(Phylo_decontam, sample_size =13335, replicates = 100,  rngseed = TRUE)
physeq_rarefied_withpro_qiime_rooted <- rarefied_withpro_qiime

Phylo_decontam_nopro_noG <- subset_samples(Phylo_decontam_nopro, Body_Type != "G")

#to send back to qiime2
library(biomformat)
otu_biom2 <- as(otu_table(Phylo_decontam), "matrix")
otu_biom2 <- make_biom(data=otu_biom2)
write_biom(otu_biom2, "otu_decontam.biom")

otu_biom <- as(otu_table(Phylo_decontam_nopro), "matrix")
otu_biom <- make_biom(data=otu_biom)
write_biom(otu_biom, "otu_decontam_nopro.biom")


otu_biom <- as(otu_table(Phylo_decontam_nopro_noG), "matrix")
otu_biom <- make_biom(data=otu_biom)
write_biom(otu_biom, "otu_decontam_nopro_nogaster.biom")

otu_forqiime2 <- as(otu_table(physeq_rarefied_withpro_qiime_rooted), "matrix")
otu_biom2 <- make_biom(data=otu_forqiime2)
write_biom(otu_biom2, "physeq_rarefied_withpro_qiime_rooted.biom")

otu_forqiime2 <- as(otu_table(physeq_rarefied_nopro_qiime_rooted), "matrix")
otu_biom2 <- make_biom(data=otu_forqiime2)
write_biom(otu_biom2, "physeq_rarefied_nopro_qiime_rooted.biom") 
#### Relative abundance plot ordered by gut compartment####
#Count # phyla to set color palette
physeq3 = transform_sample_counts(physeq_rarefied_nopro_qiime_rooted, function(x) x/sum(x))
physeq3
glom <- tax_glom(physeq3, taxrank = 'Order')
glom # should list # taxa as # phyla
data_glom<- psmelt(glom) # create dataframe from phyloseq object
data_glom$Order <- as.character(data_glom$Order) #convert to character

#simple way to rename phyla with < 1% abundance
data_glom$Order[data_glom$Abundance < 0.01] <- "< 1% abund."

#Count # phyla to set color palette
Count = length(unique(data_glom$Order))
Count

unique(data_glom$Order)

data_glom$Order <- factor(data_glom$Order, levels = c("Opitutales", "Rhizobiales", "Enterobacteriales", "Acetobacterales", "Rickettsiales", "Xanthomonadales", "Pseudomonadales", "Betaproteobacteriales", "Bacillales", "Micrococcales", "Lactobacillales", "Sphingobacteriales", "Campylobacterales", "Flavobacteriales", "Corynebacteriales", "Cytophagales", "Chloroplast", "Kineosporiales", "Clostridiales", "Pseudonocardiales", "Rhodobacterales", "Bacteroidales", "Bdellovibrionales", "Sphingomonadales", "Propionibacteriales", "Pasteurellales", "Bifidobacteriales", "Actinomycetales", "Erysipelotrichales", "Selenomonadales", "Thermales", "Streptomycetales", "Fusobacteriales", "Caulobacterales", "Chitinophagales", "Blastocatellales", "< 1% abund."))

#plot with condensed phyla into "unknown" category
bodytype_names <- c(
  `G` = "Gaster",
  `I` = "Ileum",
  `M` = "Midgut",
  `R` = "Rectum",
  `C` = "Crop")
data_glom$Body_Type_reorder = factor(data_glom$Body_Type, levels=c("C", "M", "I", "R","G"))
spatial_plot <- ggplot(data= data_glom, aes(x=Sample, y=Abundance, fill=Order)) + facet_grid(~Body_Type_reorder, scales = "free", space="free", labeller = as_labeller(bodytype_names))
spatial_plot + geom_bar(aes(), stat="identity", position="fill") +
  scale_fill_manual(values = c("royalblue4", "deepskyblue", "blue", "cyan2", "darkorchid",
                               "gold1", "forestgreen", "firebrick", "mediumspringgreen", "darkorange1",
                               "saddlebrown", "deeppink", "slategray2", "seagreen", "snow4", "purple4", 
                               "azure", "burlywood4", "chartreuse", "hotpink", "khaki",
                               "darkseagreen1", "yellow2", "orchid2", "paleturquoise1", "deepskyblue4",
                               "gold2", "gray40", "gray74", "peru", "red1", "black", "oldlace", "darkred", "chocolate4", "gold3", "lavenderblush1")) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5)) +
  theme(axis.text.x = element_blank(),
        strip.text.x = element_text(size=20, face= "bold")) 

Rel <- spatial_plot + geom_bar(aes(), stat="identity", position="fill") +
  scale_fill_manual(values = c("royalblue4", "deepskyblue", "blue", "cyan2", "darkorchid",
                               "gold1", "forestgreen", "firebrick", "mediumspringgreen", "darkorange1",
                               "saddlebrown", "deeppink", "slategray2", "seagreen", "snow4", "purple4", 
                               "azure", "burlywood4", "chartreuse", "hotpink", "khaki",
                               "darkseagreen1", "yellow2", "orchid2", "paleturquoise1", "deepskyblue4",
                               "gold2", "gray40", "gray74", "peru", "red1", "black", "oldlace", "darkred", "chocolate4", "gold3", "lavenderblush1")) +
  theme(legend.position="none") + guides(fill=guide_legend(nrow=3)) +  ylab("Relative Abundance") +
  theme(axis.text.x = element_blank(),
        strip.text.x = element_text(size=10, face= "bold")) 

##### Relative Abundance plots ordered by species and seperated by gut compartment ####
library(ggpubr)
library(ggplot2)
#### read abundance plot, with procryptocerus by species by CROP ###
physeq_C <- subset_samples(physeq_rarefied_withpro_qiime_rooted_manuscript, Body_Type == "C")
physeq3 = transform_sample_counts(physeq_C, function(x) x/sum(x))
physeq3
glom <- tax_glom(physeq3, taxrank = 'Order')
glom # should list # taxa as # phyla
data_glom<- psmelt(glom) # create dataframe from phyloseq object
data_glom$Order <- as.character(data_glom$Order) #convert to character

#simple way to rename phyla with < 1% abundance
data_glom$Order[data_glom$Abundance < 0.01] <- "< 1% abund."

#Count # phyla to set color palette
Count = length(unique(data_glom$Order))
Count

unique(data_glom$Order)

data_glom$Order <- factor(data_glom$Order, levels = c("Opitutales", "Rhizobiales", "Enterobacteriales", "Acetobacterales", "Rickettsiales", "Xanthomonadales", "Pseudomonadales", "Betaproteobacteriales", "Bacillales", "Micrococcales", "Lactobacillales", "Sphingobacteriales", "Campylobacterales", "Flavobacteriales", "Corynebacteriales", "Cytophagales", "Chloroplast", "Kineosporiales", "Clostridiales", "Pseudonocardiales", "Rhodobacterales", "Bacteroidales", "Bdellovibrionales", "Sphingomonadales", "Propionibacteriales", "Pasteurellales", "Bifidobacteriales", "Actinomycetales", "Erysipelotrichales", "Selenomonadales", "Thermales", "Streptomycetales", "Fusobacteriales", "Caulobacterales", "Chitinophagales", "Blastocatellales", "< 1% abund."))

#supp.labels <- c("atratus", "auricomus", "cordatus", "grandinosus", "minutus", "multispinosus", "placidus", "setulifer", "simillimus", "texanus", "unimaculatus", "Procryptocerus")
supp.labels <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
names(supp.labels) <- c("C_atratus", "C_auricomus", "C_cordatus", "C_grandinosus", "C_minutus", "C_multispinosus", "C_placidus", "C_setulifer", "C_simillimus", "C_texanus", "C_unimaculatus", "Procryptocerus_sp")
facet_fill <- c("deepskyblue1", "green", "red3", "darkorange1", "darkviolet", "khaki1", "blue", "red1", "plum1", "goldenrod2", "green4")
#plot with condensed phyla into "unknown" category
spatial_plot <- ggplot(data= data_glom, aes(x=Sample, y=Abundance, fill=Order)) + facet_grid(~sample_Species, scales = "free", space= "free", labeller = labeller(sample_Species = supp.labels)) + theme(strip.text.x = element_text(size = 12, face = "bold"))
CTBP <- spatial_plot + geom_bar(aes(), stat="identity", position="fill") +
  scale_fill_manual(values = c("royalblue4", "deepskyblue", "blue", "cyan2", "darkorchid",
                               "gold1", "forestgreen", "firebrick", "mediumspringgreen", "darkorange1",
                               "saddlebrown", "deeppink", "slategray2", "seagreen", "snow4", "purple4", 
                               "azure", "burlywood4", "chartreuse", "hotpink", "khaki",
                               "darkseagreen1", "yellow2", "orchid2", "paleturquoise1", "deepskyblue4",
                               "gold2", "gray40", "gray74", "peru", "red1", "black", "oldlace", "darkred", "chocolate4", "gold3", "lavenderblush1")) +
  theme(legend.position="none") + guides(fill=guide_legend(nrow=5)) + theme(axis.text.x=element_text(size = 8, angle =90)) +  ggtitle("Crop") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size= 20)) 
CTBP
g <- ggplot_gtable(ggplot_build(CTBP))

strips <- which(grepl('strip-', g$layout$name))

pal <- c("deepskyblue1", "green", "darkorange1", "darkviolet", "khaki1", "blue", "red1", "plum1", "goldenrod2", "green4", "grey")


for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- pal[i]
}

CTBP <- g
plot(CTBP)

#### read abundance plot, with procryptocerus by species by MIDGUT ###
physeq_M <- subset_samples(physeq_rarefied_withpro_qiime_rooted_manuscript, Body_Type == "M")
physeq3 = transform_sample_counts(physeq_M, function(x) x/sum(x))
physeq3
glom <- tax_glom(physeq3, taxrank = 'Order')
glom # should list # taxa as # phyla
data_glom<- psmelt(glom) # create dataframe from phyloseq object
data_glom$Order <- as.character(data_glom$Order) #convert to character

#simple way to rename phyla with < 1% abundance
data_glom$Order[data_glom$Abundance < 0.01] <- "< 1% abund."

#Count # phyla to set color palette
Count = length(unique(data_glom$Order))
Count

unique(data_glom$Order)

data_glom$Order <- factor(data_glom$Order, levels = c("Opitutales", "Rhizobiales", "Enterobacteriales", "Acetobacterales", "Rickettsiales", "Xanthomonadales", "Pseudomonadales", "Betaproteobacteriales", "Bacillales", "Micrococcales", "Lactobacillales", "Sphingobacteriales", "Campylobacterales", "Flavobacteriales", "Corynebacteriales", "Cytophagales", "Chloroplast", "Kineosporiales", "Clostridiales", "Pseudonocardiales", "Rhodobacterales", "Bacteroidales", "Bdellovibrionales", "Sphingomonadales", "Propionibacteriales", "Pasteurellales", "Bifidobacteriales", "Actinomycetales", "Erysipelotrichales", "Selenomonadales", "Thermales", "Streptomycetales", "Fusobacteriales", "Caulobacterales", "Chitinophagales", "Blastocatellales", "< 1% abund."))
supp.labels <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
#supp.labels <- c("atratus", "auricomus", "cordatus", "grandinosus", "minutus", "multispinosus", "placidus", "setulifer", "simillimus", "texanus", "unimaculatus", "Procryptocerus")
names(supp.labels) <- c("C_atratus", "C_auricomus", "C_cordatus", "C_grandinosus", "C_minutus", "C_multispinosus", "C_placidus", "C_setulifer", "C_simillimus", "C_texanus", "C_unimaculatus", "Procryptocerus_sp")

#plot with condensed phyla into "unknown" category
spatial_plot <- ggplot(data= data_glom, aes(x=Sample, y=Abundance, fill=Order)) + facet_grid(~sample_Species, scales = "free", space= "free", labeller = labeller(sample_Species = supp.labels)) + theme(strip.text.x = element_text(size = 12, face = "bold"))
MTBP <- spatial_plot + geom_bar(aes(), stat="identity", position="fill") +
  scale_fill_manual(values = c("royalblue4", "deepskyblue", "blue", "cyan2", "darkorchid",
                               "gold1", "forestgreen", "firebrick", "mediumspringgreen", "darkorange1",
                               "black", "peru", "snow4", "deeppink", "purple4", 
                               "azure", "burlywood4", "chartreuse", "hotpink", "khaki",
                               "darkseagreen1", "yellow2", "orchid2", "paleturquoise1", "deepskyblue4",
                               "gold2", "gray40", "gray74", "peru", "red1", "black", "oldlace", "darkred", "chocolate4", "gold3", "lavenderblush1")) +
  theme(legend.position="none") +  theme(axis.text.x=element_text(size = 8, angle =90)) +  ggtitle("Midgut") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size =20))

g <- ggplot_gtable(ggplot_build(MTBP))
strips <- which(grepl('strip-', g$layout$name))
pal <- c("deepskyblue1", "green", "red3", "darkorange1", "darkviolet", "khaki1", "blue", "red1", "plum1", "goldenrod2", "green4", "grey")

for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- pal[i]
}

MTBP <- g
plot(MTBP)

#### read abundance plot, with procryptocerus by species by ILEUM ###
physeq_I <- subset_samples(physeq_rarefied_withpro_qiime_rooted_manuscript, Body_Type == "I")
physeq3 = transform_sample_counts(physeq_I, function(x) x/sum(x))
physeq3
glom <- tax_glom(physeq3, taxrank = 'Order')
glom # should list # taxa as # phyla
data_glom<- psmelt(glom) # create dataframe from phyloseq object
data_glom$Order <- as.character(data_glom$Order) #convert to character

#simple way to rename phyla with < 1% abundance
data_glom$Order[data_glom$Abundance < 0.01] <- "< 1% abund."

#Count # phyla to set color palette
Count = length(unique(data_glom$Order))
Count

unique(data_glom$Order)

data_glom$Order <- factor(data_glom$Order, levels = c("Opitutales", "Rhizobiales", "Enterobacteriales", "Acetobacterales", "Rickettsiales", "Xanthomonadales", "Pseudomonadales", "Betaproteobacteriales", "Bacillales", "Micrococcales", "Lactobacillales", "Sphingobacteriales", "Campylobacterales", "Flavobacteriales", "Corynebacteriales", "Cytophagales", "Chloroplast", "Kineosporiales", "Clostridiales", "Pseudonocardiales", "Rhodobacterales", "Bacteroidales", "Bdellovibrionales", "Sphingomonadales", "Propionibacteriales", "Pasteurellales", "Bifidobacteriales", "Actinomycetales", "Erysipelotrichales", "Selenomonadales", "Thermales", "Streptomycetales", "Fusobacteriales", "Caulobacterales", "Chitinophagales", "Blastocatellales", "< 1% abund."))
supp.labels <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
#supp.labels <- c("atratus", "auricomus", "cordatus", "grandinosus", "minutus", "multispinosus", "placidus", "setulifer", "simillimus", "texanus", "unimaculatus", "Procryptocerus")
names(supp.labels) <- c("C_atratus", "C_auricomus", "C_cordatus", "C_grandinosus", "C_minutus", "C_multispinosus", "C_placidus", "C_setulifer", "C_simillimus", "C_texanus", "C_unimaculatus", "Procryptocerus_sp")

#plot with condensed phyla into "unknown" category
spatial_plot <- ggplot(data= data_glom, aes(x=Sample, y=Abundance, fill=Order)) + facet_grid(~sample_Species, scales = "free", space= "free", labeller = labeller(sample_Species = supp.labels)) + theme(strip.text.x = element_text(size = 12, face = "bold")) 
ITBP <- spatial_plot + geom_bar(aes(), stat="identity", position="fill") +
  scale_fill_manual(values = c("royalblue4", "deepskyblue", "blue", "cyan2", "darkorchid",
                               "gold1", "forestgreen", "firebrick", "mediumspringgreen", "darkorange1",
                               "saddlebrown", "deeppink", "slategray2", "seagreen", "snow4", "purple4", 
                               "azure", "burlywood4", "chartreuse", "hotpink", "khaki",
                               "darkseagreen1", "yellow2", "orchid2", "paleturquoise1", "deepskyblue4",
                               "gold2", "gray40", "gray74", "peru", "red1", "black", "oldlace", "darkred", "chocolate4", "gold3", "lavenderblush1")) +
  theme(legend.position="none") + theme(axis.text.x=element_text(size = 8, angle =90)) +  ggtitle("Ileum") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size =20))
g <- ggplot_gtable(ggplot_build(ITBP))

strips <- which(grepl('strip-', g$layout$name))

pal <- c("deepskyblue1", "green", "red3", "darkorange1", "darkviolet", "khaki1", "blue", "red1", "plum1", "goldenrod2", "green4", "grey")


for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- pal[i]
}

ITBP <- g
plot(ITBP)

#### read abundance plot, with procryptocerus by species by RECTUM ##
physeq_R <- subset_samples(physeq_rarefied_withpro_qiime_rooted_manuscript, Body_Type == "R")
physeq3 = transform_sample_counts(physeq_R, function(x) x/sum(x))
physeq3
glom <- tax_glom(physeq3, taxrank = 'Order')
glom # should list # taxa as # phyla
data_glom<- psmelt(glom) # create dataframe from phyloseq object
data_glom$Order <- as.character(data_glom$Order) #convert to character

#simple way to rename phyla with < 1% abundance
data_glom$Order[data_glom$Abundance < 0.01] <- "< 1% abund."

#Count # phyla to set color palette
Count = length(unique(data_glom$Order))
Count

unique(data_glom$Order)

data_glom$Order <- factor(data_glom$Order, levels = c("Opitutales", "Rhizobiales", "Enterobacteriales", "Acetobacterales", "Rickettsiales", "Xanthomonadales", "Pseudomonadales", "Betaproteobacteriales", "Bacillales", "Micrococcales", "Lactobacillales", "Sphingobacteriales", "Campylobacterales", "Flavobacteriales", "Corynebacteriales", "Cytophagales", "Chloroplast", "Kineosporiales", "Clostridiales", "Pseudonocardiales", "Rhodobacterales", "Bacteroidales", "Bdellovibrionales", "Sphingomonadales", "Propionibacteriales", "Pasteurellales", "Bifidobacteriales", "Actinomycetales", "Erysipelotrichales", "Selenomonadales", "Thermales", "Streptomycetales", "Fusobacteriales", "Caulobacterales", "Chitinophagales", "Blastocatellales", "< 1% abund."))
supp.labels <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
#supp.labels <- c("atratus", "auricomus", "cordatus", "grandinosus", "minutus", "multispinosus", "placidus", "setulifer", "simillimus", "texanus", "unimaculatus", "Procryptocerus")
names(supp.labels) <- c("C_atratus", "C_auricomus", "C_cordatus", "C_grandinosus", "C_minutus", "C_multispinosus", "C_placidus", "C_setulifer", "C_simillimus", "C_texanus", "C_unimaculatus", "Procryptocerus_sp")

#plot with condensed phyla into "unknown" category
spatial_plot <- ggplot(data= data_glom, aes(x=Sample, y=Abundance, fill=Order)) + facet_grid(~sample_Species, scales = "free", space= "free", labeller = labeller(sample_Species = supp.labels)) + theme(strip.text.x = element_text(size = 12, face = "bold")) 
RTBP <- spatial_plot + geom_bar(aes(), stat="identity", position="fill") +
  scale_fill_manual(values = c("royalblue4", "deepskyblue", "blue", "cyan2", "darkorchid",
                               "gold1", "forestgreen", "gold3", "firebrick", "gray74", "purple4", "mediumspringgreen", "darkorange1", "black", "orchid2", "khaki", "peru",
                               "saddlebrown", "deeppink", "slategray2", "seagreen", "snow4", "purple4", 
                               "azure", "burlywood4", "chartreuse", "hotpink", "oldlace",
                               "darkseagreen1", "yellow2", "paleturquoise1", "deepskyblue4",
                               "gold2", "gray40", "peru", "red1", "darkred", "chocolate4", "lavenderblush1")) +
  theme(legend.position="none") + theme(legend.title = element_text(face="bold", size = 16)) + theme(legend.text = element_text(size = 14)) + guides(fill=guide_legend(nrow=5)) + theme(axis.text.x=element_text(size = 8, angle =90)) +  ggtitle("Rectum") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size =20)) 

g <- ggplot_gtable(ggplot_build(RTBP))

strips <- which(grepl('strip-', g$layout$name))

pal <- c("deepskyblue1", "green", "red3", "darkorange1", "darkviolet", "khaki1", "blue", "red1", "plum1", "goldenrod2", "green4", "grey")


for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- pal[i]
}

RTBP <- g
plot(RTBP)
#### read abundance plot, with procryptocerus by species by GASTER ###
physeq_G <- subset_samples(physeq_rarefied_withpro_qiime_rooted_manuscript, Body_Type == "G")
physeq3 = transform_sample_counts(physeq_G, function(x) x/sum(x))
physeq3
glom <- tax_glom(physeq3, taxrank = 'Order')
glom # should list # taxa as # phyla
data_glom<- psmelt(glom) # create dataframe from phyloseq object
data_glom$Order <- as.character(data_glom$Order) #convert to character

#simple way to rename phyla with < 1% abundance
data_glom$Order[data_glom$Abundance < 0.01] <- "< 1% abund."

#Count # phyla to set color palette
Count = length(unique(data_glom$Order))
Count

unique(data_glom$Order)

data_glom$Order <- factor(data_glom$Order, levels = c("Opitutales", "Rhizobiales", "Enterobacteriales", "Acetobacterales", "Rickettsiales", "Xanthomonadales", "Pseudomonadales", "Betaproteobacteriales", "Bacillales", "Micrococcales", "Lactobacillales", "Sphingobacteriales", "Campylobacterales", "Flavobacteriales", "Corynebacteriales", "Cytophagales", "Chloroplast", "Kineosporiales", "Clostridiales", "Pseudonocardiales", "Rhodobacterales", "Bacteroidales", "Bdellovibrionales", "Sphingomonadales", "Propionibacteriales", "Pasteurellales", "Bifidobacteriales", "Actinomycetales", "Erysipelotrichales", "Selenomonadales", "Thermales", "Streptomycetales", "Fusobacteriales", "Caulobacterales", "Chitinophagales", "Blastocatellales", "< 1% abund."))
supp.labels <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
#supp.labels <- c("atratus", "auricomus", "cordatus", "grandinosus", "minutus", "multispinosus", "placidus", "setulifer", "simillimus", "texanus", "unimaculatus", "Procryptocerus")
names(supp.labels) <- c("C_atratus", "C_auricomus", "C_cordatus", "C_grandinosus", "C_minutus", "C_multispinosus", "C_placidus", "C_setulifer", "C_simillimus", "C_texanus", "C_unimaculatus", "Procryptocerus_sp")

#plot with condensed phyla into "unknown" category
spatial_plot <- ggplot(data= data_glom, aes(x=Sample, y=Abundance, fill=Order)) + facet_grid(~sample_Species, scales = "free", space= "free", labeller = labeller(sample_Species = supp.labels)) + theme(strip.text.x = element_text(size = 12, face = "bold"))
GTBP <- spatial_plot + geom_bar(aes(), stat="identity", position="fill") +
  scale_fill_manual(values = c("royalblue4", "deepskyblue", "blue", "cyan2", "darkorchid",
                               "gold1", "forestgreen", "firebrick", "mediumspringgreen", "darkorange1",
                               "saddlebrown", "deeppink", "slategray2", "seagreen", "snow4", "purple4", 
                               "azure", "burlywood4", "chartreuse", "hotpink", "khaki",
                               "darkseagreen1", "yellow2", "orchid2", "paleturquoise1", "deepskyblue4",
                               "gold2", "gray40", "gray74", "peru", "red1", "black", "oldlace", "darkred", "chocolate4", "gold3", "lavenderblush1")) +
  theme(legend.position="none") + theme(axis.text.x=element_text(size = 8, angle =90)) +  ggtitle("Gaster") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20))
g <- ggplot_gtable(ggplot_build(GTBP))

strips <- which(grepl('strip-', g$layout$name))

pal <- c("deepskyblue1", "green", "red3", "darkorange1", "darkviolet", "khaki1", "blue", "red1", "plum1", "goldenrod2", "green4", "grey")


for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- pal[i]
}

GTBP <- g
plot(GTBP)

#### all by body compartment taxa bar plots combined ###
library(ggpubr)
library(gridExtra)
grid.arrange(GTBP, CTBP, MTBP, ITBP, RTBP)
#pdf 13x23


### Relative Abundance of Caste Types ####
#Queen Caste#
physeq_Q <- subset_samples(physeq_rarefied_withpro_qiime_rooted_manuscript, Caste_Type == "Q")
physeq3 = transform_sample_counts(physeq_Q, function(x) x/sum(x))
physeq3
glom <- tax_glom(physeq3, taxrank = 'Order')
glom # should list # taxa as # phyla
data_glom<- psmelt(glom) # create dataframe from phyloseq object
data_glom$Order <- as.character(data_glom$Order) #convert to character

#simple way to rename phyla with < 1% abundance
data_glom$Order[data_glom$Abundance < 0.01] <- "< 1% abund."

#Count # phyla to set color palette
Count = length(unique(data_glom$Order))
Count

unique(data_glom$Order)

data_glom$Order <- factor(data_glom$Order, levels = c("Opitutales", "Rhizobiales", "Enterobacteriales", "Acetobacterales", "Rickettsiales", "Xanthomonadales", "Pseudomonadales", "Betaproteobacteriales", "Bacillales", "Micrococcales", "Lactobacillales", "Sphingobacteriales", "Campylobacterales", "Flavobacteriales", "Corynebacteriales", "Cytophagales", "Chloroplast", "Kineosporiales", "Clostridiales", "Pseudonocardiales", "Rhodobacterales", "Bacteroidales", "Bdellovibrionales", "Sphingomonadales", "Propionibacteriales", "Pasteurellales", "Bifidobacteriales", "Actinomycetales", "Erysipelotrichales", "Selenomonadales", "Thermales", "Streptomycetales", "Fusobacteriales", "Caulobacterales", "Chitinophagales", "Blastocatellales", "< 1% abund."))
supp.labels <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
#supp.labels <- c("atratus", "auricomus", "cordatus", "grandinosus", "minutus", "multispinosus", "placidus", "setulifer", "simillimus", "texanus", "unimaculatus", "Procryptocerus")
names(supp.labels) <- c("C_atratus", "C_auricomus", "C_cordatus", "C_grandinosus", "C_minutus", "C_multispinosus", "C_placidus", "C_setulifer", "C_simillimus", "C_texanus", "C_unimaculatus", "Procryptocerus_sp")

#plot with condensed phyla into "unknown" category
spatial_plot <- ggplot(data= data_glom, aes(x=Sample, y=Abundance, fill=Order)) + facet_grid(~sample_Species, scales = "free", space= "free", labeller = labeller(sample_Species = supp.labels)) + theme(strip.text.x = element_text(size = 12, face = "bold"))
MTBP <- spatial_plot + geom_bar(aes(), stat="identity", position="fill") +
  scale_fill_manual(values = c("royalblue4", "deepskyblue", "cyan2", "darkorchid",
                               "gold1", "forestgreen", "firebrick", "mediumspringgreen",
                               "darkorange1", "khaki", "hotpink", "deeppink", "purple4", 
                               "azure", "burlywood4", "chartreuse", "hotpink", "khaki",
                               "darkseagreen1", "yellow2", "orchid2", "paleturquoise1", "deepskyblue4",
                               "gold2", "gray40", "gray74", "peru", "red1", "black", "oldlace", "darkred", "chocolate4", "gold3", "lavenderblush1"))  + theme(axis.text.x=element_text(size = 8, angle =90)) +  ggtitle("Queen") +   theme(legend.position="none") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size =20))

g <- ggplot_gtable(ggplot_build(MTBP))
strips <- which(grepl('strip-', g$layout$name))
pal <- c("blue", "red1", "plum1", "darkorange1", "darkviolet", "khaki1", "blue", "red1", "plum1", "goldenrod2", "green4", "grey")

for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- pal[i]
}

Queen_Rel <- g
plot(Queen_Rel)

#Worker Caste#
physeq_W <- subset_samples(physeq_rarefied_withpro_qiime_rooted_manuscript, Caste_Type == "W")
physeq3 = transform_sample_counts(physeq_W, function(x) x/sum(x))
physeq3
glom <- tax_glom(physeq3, taxrank = 'Order')
glom # should list # taxa as # phyla
data_glom<- psmelt(glom) # create dataframe from phyloseq object
data_glom$Order <- as.character(data_glom$Order) #convert to character

#simple way to rename phyla with < 1% abundance
data_glom$Order[data_glom$Abundance < 0.01] <- "< 1% abund."

#Count # phyla to set color palette
Count = length(unique(data_glom$Order))
Count

unique(data_glom$Order)

data_glom$Order <- factor(data_glom$Order, levels = c("Opitutales", "Rhizobiales", "Enterobacteriales", "Acetobacterales", "Rickettsiales", "Xanthomonadales", "Pseudomonadales", "Betaproteobacteriales", "Bacillales", "Micrococcales", "Lactobacillales", "Sphingobacteriales", "Campylobacterales", "Flavobacteriales", "Corynebacteriales", "Cytophagales", "Chloroplast", "Kineosporiales", "Clostridiales", "Pseudonocardiales", "Rhodobacterales", "Bacteroidales", "Bdellovibrionales", "Sphingomonadales", "Propionibacteriales", "Pasteurellales", "Bifidobacteriales", "Actinomycetales", "Erysipelotrichales", "Selenomonadales", "Thermales", "Streptomycetales", "Fusobacteriales", "Caulobacterales", "Chitinophagales", "Blastocatellales", "< 1% abund."))
supp.labels <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
#supp.labels <- c("atratus", "auricomus", "cordatus", "grandinosus", "minutus", "multispinosus", "placidus", "setulifer", "simillimus", "texanus", "unimaculatus", "Procryptocerus")
names(supp.labels) <- c("C_atratus", "C_auricomus", "C_cordatus", "C_grandinosus", "C_minutus", "C_multispinosus", "C_placidus", "C_setulifer", "C_simillimus", "C_texanus", "C_unimaculatus", "Procryptocerus_sp")

#plot with condensed phyla into "unknown" category
spatial_plot <- ggplot(data= data_glom, aes(x=Sample, y=Abundance, fill=Order)) + facet_grid(~sample_Species, scales = "free", space= "free", labeller = labeller(sample_Species = supp.labels)) + theme(strip.text.x = element_text(size = 12, face = "bold"))
MTBP <- spatial_plot + geom_bar(aes(), stat="identity", position="fill") +
  scale_fill_manual(values = c("royalblue4", "deepskyblue", "blue", "cyan2", "darkorchid",
                               "gold1", "forestgreen", "firebrick", "mediumspringgreen", "darkorange1",
                               "black", "peru", "snow4", "deeppink", "purple4", 
                               "azure", "burlywood4", "chartreuse", "hotpink", "khaki",
                               "darkseagreen1", "yellow2", "orchid2", "paleturquoise1", "deepskyblue4",
                               "gold2", "gray40", "gray74", "peru", "red1", "black", "oldlace", "darkred", "chocolate4", "gold3", "lavenderblush1")) +
  theme(legend.position="none") +  theme(axis.text.x=element_text(size = 8, angle =90)) +  ggtitle("Worker") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size =20))

g <- ggplot_gtable(ggplot_build(MTBP))
strips <- which(grepl('strip-', g$layout$name))
pal <- c("deepskyblue1", "green", "red3", "darkorange1", "darkviolet", "khaki1", "blue", "red1", "plum1", "goldenrod2", "green4", "grey")

for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- pal[i]
}

Worker_Rel <- g
plot(Worker_Rel)

#Soldier Caste#
physeq_S <- subset_samples(physeq_rarefied_withpro_qiime_rooted_manuscript, Caste_Type == "S")
physeq3 = transform_sample_counts(physeq_S, function(x) x/sum(x))
physeq3
glom <- tax_glom(physeq3, taxrank = 'Order')
glom # should list # taxa as # phyla
data_glom<- psmelt(glom) # create dataframe from phyloseq object
data_glom$Order <- as.character(data_glom$Order) #convert to character

#simple way to rename phyla with < 1% abundance
data_glom$Order[data_glom$Abundance < 0.01] <- "< 1% abund."

#Count # phyla to set color palette
Count = length(unique(data_glom$Order))
Count

unique(data_glom$Order)

data_glom$Order <- factor(data_glom$Order, levels = c("Opitutales", "Rhizobiales", "Enterobacteriales", "Acetobacterales", "Rickettsiales", "Xanthomonadales", "Pseudomonadales", "Betaproteobacteriales", "Bacillales", "Micrococcales", "Lactobacillales", "Sphingobacteriales", "Campylobacterales", "Flavobacteriales", "Corynebacteriales", "Cytophagales", "Chloroplast", "Kineosporiales", "Clostridiales", "Pseudonocardiales", "Rhodobacterales", "Bacteroidales", "Bdellovibrionales", "Sphingomonadales", "Propionibacteriales", "Pasteurellales", "Bifidobacteriales", "Actinomycetales", "Erysipelotrichales", "Selenomonadales", "Thermales", "Streptomycetales", "Fusobacteriales", "Caulobacterales", "Chitinophagales", "Blastocatellales", "< 1% abund."))
supp.labels <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
#supp.labels <- c("atratus", "auricomus", "cordatus", "grandinosus", "minutus", "multispinosus", "placidus", "setulifer", "simillimus", "texanus", "unimaculatus", "Procryptocerus")
names(supp.labels) <- c("C_atratus", "C_auricomus", "C_cordatus", "C_grandinosus", "C_minutus", "C_multispinosus", "C_placidus", "C_setulifer", "C_simillimus", "C_texanus", "C_unimaculatus", "Procryptocerus_sp")

#plot with condensed phyla into "unknown" category
spatial_plot <- ggplot(data= data_glom, aes(x=Sample, y=Abundance, fill=Order)) + facet_grid(~sample_Species, scales = "free", space= "free", labeller = labeller(sample_Species = supp.labels)) + theme(strip.text.x = element_text(size = 12, face = "bold"))
MTBP <- spatial_plot + geom_bar(aes(), stat="identity", position="fill") +
  scale_fill_manual(values = c("royalblue4", "deepskyblue", "blue", "cyan2", "darkorchid",
                               "gold1", "forestgreen", "firebrick", "mediumspringgreen", "darkorange1",
                               "black", "peru", "snow4", "deeppink", "purple4", 
                               "azure", "burlywood4", "chartreuse", "hotpink", "khaki",
                               "darkseagreen1", "yellow2", "orchid2", "paleturquoise1", "deepskyblue4",
                               "gold2", "gray40", "gray74", "peru", "red1", "black", "oldlace", "darkred", "chocolate4", "gold3", "lavenderblush1")) +
  theme(legend.position="none") + theme(axis.text.x=element_text(size = 8, angle =90)) +  ggtitle("Soldier") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size =20))

g <- ggplot_gtable(ggplot_build(MTBP))
strips <- which(grepl('strip-', g$layout$name))
pal <- c("red3", "darkviolet", "red1", "plum1", "darkviolet", "khaki1", "blue", "red1", "plum1", "goldenrod2", "green4", "grey")

for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- pal[i]
}

Soldier_Rel <- g
plot(Soldier_Rel)


#plot all  castes together#
grid.arrange(Queen_Rel, Worker_Rel, Soldier_Rel)
#### Alpha Diversity of Caste Type ####
#alpha diversity metrics by caste type and body type
library(FSA)
library(dplyr)
library(devtools)
library(ggpubr)
library(gridExtra)

Alpha_plots_c<- read.csv("kruskal-wallis-pairwise-Caste_Type1.csv", header= TRUE)
#overall
kruskal.test(shannon ~ Caste_Type, data = Alpha_plots_c)
OS = dunnTest(shannon ~ Caste_Type,data=Alpha_plots_c,method="bh")
OS
kruskal.test(pielou_e ~ Caste_Type, data = Alpha_plots_c)
PO = dunnTest(pielou_e ~ Caste_Type,data=Alpha_plots_c, method="bh")
PO
kruskal.test(mean_ASV_count ~ Caste_Type, data = Alpha_plots_c)
AO = dunnTest(mean_ASV_count ~ Caste_Type,data= Alpha_plots_c, method="bh")
AO
kruskal.test(faith_pd ~ Caste_Type, data = Alpha_plots_c)
FO = dunnTest(faith_pd ~ Caste_Type,data=Alpha_plots_c, method="bh")
FO
#Gaster - no differences in caste
#crop
Crop_Alpha <- filter(Alpha_plots_c, Body_Type == "C")
kruskal.test(shannon ~ Caste_Type, data = Crop_Alpha)
CS = dunnTest(shannon ~ Caste_Type,data=Crop_Alpha,method="bh")
CS
kruskal.test(pielou_e ~ Caste_Type, data = Crop_Alpha)
PC = dunnTest(pielou_e ~ Caste_Type,data=Crop_Alpha, method="bh")
PC
kruskal.test(mean_ASV_count ~ Caste_Type, data = Crop_Alpha)
AC = dunnTest(mean_ASV_count ~ Caste_Type,data= Crop_Alpha, method="bh")
AC
kruskal.test(faith_pd ~ Caste_Type, data = Crop_Alpha)
FC = dunnTest(faith_pd ~ Caste_Type,data=Crop_Alpha, method="bh")
FC
#rectum
Rectum_Alpha <- filter(Alpha_plots_c, Body_Type == "R")
kruskal.test(shannon ~ Caste_Type, data = Rectum_Alpha)
RS = dunnTest(shannon ~ Caste_Type,data=Rectum_Alpha,method="bh")
RS
kruskal.test(pielou_e ~ Caste_Type, data = Rectum_Alpha)
PR = dunnTest(pielou_e ~ Caste_Type,data=Rectum_Alpha, method="bh")
PR
kruskal.test(mean_ASV_count ~ Caste_Type, data = Rectum_Alpha)
AR = dunnTest(mean_ASV_count ~ Caste_Type,data= Rectum_Alpha, method="bh")
AR
kruskal.test(faith_pd ~ Caste_Type, data = Rectum_Alpha)
FR = dunnTest(faith_pd ~ Caste_Type,data=Rectum_Alpha, method="bh")
FR
#midgut
Midgut_Alpha <- filter(Alpha_plots_c, Body_Type == "M")
kruskal.test(shannon ~ Caste_Type, data = Midgut_Alpha)
MS = dunnTest(shannon ~ Caste_Type, data=Midgut_Alpha, method="bh")
MS
kruskal.test(pielou_e ~ Caste_Type, data = Midgut_Alpha)
PM = dunnTest(pielou_e ~ Caste_Type,data=Midgut_Alpha, method="bh")
PM
kruskal.test(mean_ASV_count ~ Caste_Type, data = Midgut_Alpha)
AM = dunnTest(mean_ASV_count ~ Caste_Type,data= Midgut_Alpha, method="bh")
AM
kruskal.test(faith_pd ~ Caste_Type, data = Midgut_Alpha)
FM = dunnTest(faith_pd ~ Caste_Type,data=Midgut_Alpha, method="bh")
FM
#ileum
Ileum_Alpha <- filter(Alpha_plots_c, Body_Type == "I")
kruskal.test(shannon ~ Caste_Type, data = Ileum_Alpha)
SI = dunnTest(shannon ~ Caste_Type,data=Ileum_Alpha, method="bh")
SI
kruskal.test(pielou_e ~ Caste_Type, data = Ileum_Alpha)
PI = dunnTest(pielou_e ~ Caste_Type,data=Ileum_Alpha, method="bh")
PI
kruskal.test(mean_ASV_count ~ Caste_Type, data = Ileum_Alpha)
AI = dunnTest(mean_ASV_count ~ Caste_Type,data=Ileum_Alpha, method="bh")
AI
kruskal.test(faith_pd ~ Caste_Type, data = Ileum_Alpha)
FI = dunnTest(faith_pd ~ Caste_Type,data=Ileum_Alpha, method="bh")
FI

#### Raw  abundance plot####
#abundance plot by body type
physeq_noX <- subset_samples(Phylo_decontam_nopro, Body_Type != "X")
glom <- tax_glom(physeq_noX, taxrank = 'Order')
glom # should list # taxa as # phyla
data_glom<- psmelt(glom) # create dataframe from phyloseq object
data_glom$Order <- as.character(data_glom$Order) #convert to character

#simple way to rename phyla with < 1% abundance
data_glom$Order[data_glom$Abundance < 0.01] <- "< 1% abund."

#Count # phyla to set color palette
Count = length(unique(data_glom$Order))
Count

unique(data_glom$Order)

data_glom$Order <- factor(data_glom$Order, levels = c("Opitutales", "Rhizobiales", "Enterobacteriales", "Acetobacterales", "Rickettsiales", "Xanthomonadales", "Pseudomonadales", "Betaproteobacteriales", "Bacillales", "Micrococcales", "Lactobacillales", "Sphingobacteriales", "Campylobacterales", "Flavobacteriales", "Corynebacteriales", "Cytophagales", "Chloroplast", "Kineosporiales", "Clostridiales", "Pseudonocardiales", "Rhodobacterales", "Bacteroidales", "Bdellovibrionales", "Sphingomonadales", "Propionibacteriales", "Pasteurellales", "Bifidobacteriales", "Actinomycetales", "Erysipelotrichales", "Selenomonadales", "Thermales", "Streptomycetales", "Fusobacteriales", "Caulobacterales", "Chitinophagales", "Blastocatellales", "< 1% abund."))

#plot with condensed phyla into "unknown" category
bodytype_names <- c(
  `G` = "Gaster",
  `I` = "Ileum",
  `M` = "Midgut",
  `R` = "Rectum",
  `C` = "Crop")
data_glom$Body_Type_reorder = factor(data_glom$Body_Type, levels=c("C", "M", "I", "R", "G"))
spatial_plot <- ggplot(data= data_glom, aes(x=Sample, y=Abundance, fill=Order)) + facet_grid(~Body_Type_reorder, scales = "free", space="free", labeller = as_labeller(bodytype_names))
Raw <- spatial_plot + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("royalblue4", "deepskyblue", "blue", "ivory" ,"cyan2", "darkorchid",
                               "gold1", "forestgreen", "firebrick", "mediumspringgreen", "darkorange1",
                               "saddlebrown", "deeppink", "slategray2", "seagreen", "wheat", "snow4", "purple4", "oldlace", "black",
                               "azure", "burlywood4","chocolate4","chartreuse", "khaki",
                               "darkseagreen1", "yellow2", "orchid2", "paleturquoise1", "deepskyblue4",
                               "gold2", "gray40", "gray74", "peru", "red1","hotpink", "darkred", "gold3", "lavenderblush1")) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5)) +
  theme(axis.text.x = element_blank(),
        strip.text.x = element_text(size=20, face= "bold")) 

Raw <- spatial_plot + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("royalblue4", "deepskyblue", "blue", "ivory" ,"cyan2", "darkorchid",
                               "gold1", "forestgreen", "firebrick", "mediumspringgreen", "darkorange1",
                               "saddlebrown", "deeppink", "slategray2", "seagreen", "wheat", "snow4", "purple4", "oldlace", "black",
                               "azure", "burlywood4","chocolate4","chartreuse", "khaki",
                               "darkseagreen1", "yellow2", "orchid2", "paleturquoise1", "deepskyblue4",
                               "gold2", "gray40", "gray74", "peru", "red1","hotpink", "darkred", "gold3", "lavenderblush1")) +
  theme(legend.position="none") + guides(fill=guide_legend(nrow=2)) +  ylab("Raw Read Number") +
  theme(axis.text.x = element_blank(),
        strip.text.x = element_text(size=10, face= "bold")) 
Raw
grid.arrange(Rel, Raw)
#4.5x9 for figure 2
#theme(axis.text.x = element_text(angle=90, vjust= 0.5))
#### qPCR boxplots #####
library(dplyr)
library(devtools)
library(ggpubr)
library(gridExtra)
my_data <- read.csv("qPCR_averages_final.csv")
my_data_gaster <- filter(my_data, Gut.Compartment == "G")
my_data_noX <- filter(my_data, Gut.Compartment != "X")
my_data_noGorX <- filter(my_data, Gut.Compartment != "G" & Gut.Compartment != "X")
my_data_gaster <- filter(my_data, Gut.Compartment == "G")
my_data_crop<- filter(my_data, Gut.Compartment == "C")
my_data_rectum <- filter(my_data, Gut.Compartment == "R")
my_data_midgut <- filter(my_data, Gut.Compartment == "M")
my_data_ileum <- filter(my_data, Gut.Compartment == "I")

#qpcr across all gut compartments
my_data_nopro <- filter(my_data_noX, Species != "Procryptocerus")
my_data_nopro$Gut.Compartment = factor(my_data_nopro$Gut.Compartment, levels=c("C", "M", "I", "R", "G"))
pqpcr = ggboxplot(my_data_nopro, x = "Gut.Compartment", y= "Average", color = "Gut.Compartment") + scale_x_discrete(labels = c("Crop", "Midgut", "Ileum", "Rectum", "Gaster"))
pqpcr = pqpcr + geom_boxplot(fill = c("darkblue", "darkorchid", "darkseagreen", "firebrick", "darkgoldenrod1")) + ggtitle("qPCR values across gut compartments") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")
pqpcr
#pdf is 10x10

pairwise.wilcox.test(my_data_nopro$Average, my_data_nopro$Gut.Compartment, p.adj = "fdr")
pairwise.wilcox.test(my_data_noX$Average, my_data_noX$Species, p.adj = "fdr")
pairwise.wilcox.test(my_data_nopro$Average, my_data_nopro$Caste, p.adj = "fdr")

#boxplots for each gut section all species combined
my_data_nopro <- filter(my_data_noX, Species != "Procryptocerus")
my_data_noX$Gut.Compartment = factor(my_data_noX$Gut.Compartment, levels=c("C", "M", "I", "R", "G"))
pqpcr = ggboxplot(my_data_noX, x = "Gut.Compartment", y= "Average", color = "Gut.Compartment") + scale_x_discrete(labels = c("Crop", "Midgut", "Ileum", "Rectum", "Gaster"))
pqpcr = pqpcr + geom_boxplot(fill = c("darkblue", "darkorchid", "green4", "cadetblue2", "darkgoldenrod1")) + ggtitle("qPCR values across gut compartments") + xlab("Gut Compartment")+ ylab("Mean copies bacterial 16S rRNA gene (rRNA/ul)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")
pqpcr

#boxplots for each gut section and species combined 
my_data_noX$Gut.Compartment = factor(my_data_noX$Gut.Compartment, levels=c("C", "M", "I", "R", "G"))
my_data_noX$Species = factor(my_data_noX$Species, levels=c("atratus", "auricomus", "cordatus", "grandinosus", "minutus", "multispinosus", "placidus", "setulifer", "simillimus", "texanus", "unimaculatus", "Procryptocerus"))
psqpcr = ggboxplot(my_data_noX, x = "Gut.Compartment", y= "Average", fill = "Species") + scale_x_discrete(labels = c("Crop", "Midgut", "Ileum", "Rectum", "Gaster")) 
psqpcr = psqpcr + scale_fill_manual(values = c("deepskyblue1", "green", "red3", "darkorange1", "darkviolet", "khaki1", "blue", "red1", "plum1", "goldenrod2", "green4", "grey")) + ggtitle("qPCR values across gut compartment and species") + xlab("Gut Compartment")+ ylab("Mean copies bacterial 16S rRNA gene (rRNA/ul)") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) + theme(legend.position="bottom")
psqpcr
# 11x19 for pdf

grid.arrange(pqpcr, psqpcr, nrow=1, widths = 1:2)

my_data_gaster$Species = factor(my_data_gaster$Species, levels=c("atratus", "auricomus", "cordatus", "grandinosus", "minutus", "multispinosus", "placidus", "setulifer", "simillimus", "texanus", "unimaculatus", "Procryptocerus"))
pgqpcr = ggboxplot(my_data_gaster, x = "Species", y= "Average", color = "Species")
pgqpcr = pgqpcr + geom_boxplot(fill = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "deeppink", "firebrick", "khaki2", "darkgreen", "brown1", "darkorange1")) + ggtitle("Mean qPCR values for Gaster") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")
pgqpcr
pairwise.wilcox.test(my_data_gaster$Average, my_data_gaster$Species, p.adj = "fdr")

my_data_crop$Species = factor(my_data_crop$Species, levels=c("atratus", "auricomus", "cordatus", "grandinosus", "minutus", "multispinosus", "placidus", "setulifer", "simillimus", "texanus", "unimaculatus","Procryptocerus"))
pcqpcr = ggboxplot(my_data_crop, x = "Species", y= "Average", color = "Species")
pcqpcr = pcqpcr + geom_boxplot(fill = c("deepskyblue1", "green", "red3", "darkorange1", "darkviolet", "khaki1", "blue", "red1", "plum1", "goldenrod2", "green4", "black")) + ggtitle("Mean qPCR values for Crop") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")
pcqpcr
pairwise.wilcox.test(my_data_crop$Average, my_data_crop$Species, p.adj = "fdr")

my_data_rectum$Species = factor(my_data_rectum$Species, levels=c("atratus", "auricomus", "cordatus", "grandinosus", "minutus", "multispinosus", "placidus", "setulifer", "simillimus", "texanus", "unimaculatus","Procryptocerus"))
prqpcr = ggboxplot(my_data_rectum, x = "Species", y= "Average", color = "Species")
prqpcr = prqpcr + geom_boxplot(fill = c("deepskyblue1", "green", "red3", "darkorange1", "darkviolet", "khaki1", "blue", "red1", "plum1", "goldenrod2", "green4", "black")) + ggtitle("Mean qPCR values for Rectum") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")
prqpcr
pairwise.wilcox.test(my_data_rectum$Average, my_data_rectum$Species, p.adj = "fdr")

my_data_ileum$Species = factor(my_data_ileum$Species, levels=c("atratus", "auricomus", "cordatus", "grandinosus", "minutus", "multispinosus", "placidus", "setulifer", "simillimus", "texanus", "unimaculatus","Procryptocerus"))
piqpcr = ggboxplot(my_data_ileum, x = "Species", y= "Average", color = "Species")
piqpcr = piqpcr + geom_boxplot(fill = c("deepskyblue1", "green", "red3", "darkorange1", "darkviolet", "khaki1", "blue", "red1", "plum1", "goldenrod2", "green4", "black")) + ggtitle("Mean qPCR values for Ileum") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")
piqpcr
pairwise.wilcox.test(my_data_ileum$Average, my_data_ileum$Species, p.adj = "fdr")

my_data_midgut$Species = factor(my_data_midgut$Species, levels=c("atratus", "auricomus", "cordatus", "grandinosus", "minutus", "multispinosus", "placidus", "setulifer", "simillimus", "texanus", "unimaculatus", "Procryptocerus"))
p = ggboxplot(my_data_midgut, x = "Species", y= "Average", color = "Species")
p = p + geom_boxplot(fill = c("deepskyblue1", "green", "red3", "darkorange1", "darkviolet", "khaki1", "blue", "red1", "plum1", "goldenrod2", "green4", "black")) + ggtitle("Mean qPCR values for Midgut") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")
p
pairwise.wilcox.test(my_data_midgut$Average, my_data_midgut$Species, p.adj = "fdr")


my_data_noX$Species = factor(my_data_noX$Species, levels=c("atratus", "auricomus", "cordatus", "grandinosus", "minutus", "multispinosus", "placidus", "setulifer", "simillimus", "texanus", "unimaculatus", "Procryptocerus"))
pall = ggboxplot(my_data_noX, x = "Species", y= "Average", color = "Species")
pall = pall + geom_boxplot(fill = c("deepskyblue1", "green", "red3", "darkorange1", "darkviolet", "khaki1", "blue", "red1", "plum1", "goldenrod2", "green4", "black")) + ggtitle("Mean qPCR values for all samples") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")
pall

grid.arrange(pcqpcr,p, piqpcr, prqpcr)
#18x30 for pdf

###STATS for LMM for qPCR data####
mixed.lmer1 <- lmer(Average ~ Gut.Compartment +  (1|Species/Colony), data = my_data_nopro, REML=T)
mixed.lmer2 <- lmer(Average ~ Gut.Compartment + Caste  + (1|Species/Colony), data =my_data_nopro, REML=T)
mixed.lmer3 <- lmer(Average ~ Gut.Compartment + Caste + Gut.Compartment:Caste + (1|Species/Colony), data = my_data_nopro, REML=T)

anova(mixed.lmer1)
anova(mixed.lmer3, mixed.lmer2, mixed.lmer1)

em <- emmeans(mixed.lmer1, "Gut.Compartment")
contrast(em, adjust = "FDR", method = "pairwise")

#### Alpha Diversity Metrics ####
#alpha diversity metrics for gut compartment 
library(gridExtra)
library(ggsignif)
library(ggpubr)
Alpha_plots <- read.csv("kruskal-wallis-pairwise-Body_Type.csv", header= TRUE)
Alpha_plots$Body_Type = factor(Alpha_plots$Body_Type, levels=c("C", "M", "I", "R", "G")) 
ps = ggboxplot(Alpha_plots, x = "Body_Type", y= "mean_shannon", color = "Body_Type") + scale_x_discrete(labels = c("Crop", "Midgut", "Ileum", "Rectum", "Gaster"))
ps= ps + ggtitle("Mean Shannon Diversity") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none") + xlab("Gut Compartment")
ps = ps + geom_boxplot(fill = c("darkblue", "darkorchid", "green4", "cadetblue2", "darkgoldenrod1"))
ps 

pd = ggboxplot(Alpha_plots, x = "Body_Type", y= "mean_faith_pd", color = "Body_Type") + scale_x_discrete(labels = c("Crop", "Midgut", "Ileum", "Rectum", "Gaster"))
pd = pd  + scale_colour_manual(values = c("darkblue", "darkorchid", "darkseagreen", "firebrick","darkgoldenrod1")) + ggtitle("Mean Faith's PD") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none") 
pd = pd + geom_boxplot(fill = c("darkblue", "darkorchid", "green4", "cadetblue2", "darkgoldenrod1")) + xlab("Gut Compartment")
pd

pe = ggboxplot(Alpha_plots, x = "Body_Type", y= "mean_evenness", color = "Body_Type") + scale_x_discrete(labels = c("Crop", "Midgut", "Ileum", "Rectum", "Gaster"))
pe = pe + scale_colour_manual(values = c("darkblue", "darkorchid", "darkseagreen", "firebrick","darkgoldenrod1")) + ggtitle("Mean Pielou's Evenness") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none") 
pe = pe + geom_boxplot(fill = c("darkblue", "darkorchid", "green4", "cadetblue2", "darkgoldenrod1")) + xlab("Gut Compartment")
pe

pa = ggboxplot(Alpha_plots, x = "Body_Type", y= "mean_ASV_count", color = "Body_Type") + scale_x_discrete(labels = c("Crop", "Midgut", "Ileum", "Rectum", "Gaster"))
pa = pa + scale_colour_manual(values = c("darkblue", "darkorchid", "darkseagreen", "firebrick","darkgoldenrod1")) + ggtitle("Mean ASV Count") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none") 
pa = pa + geom_boxplot(fill = c("darkblue", "darkorchid", "green4", "cadetblue2", "darkgoldenrod1")) + xlab("Gut Compartment")
pa

grid.arrange(ps, pd, pe, pa)
# pdf is 10x14

#alpha diversity metrics by caste type 
Alpha_plots_c<- read.csv("kruskal-wallis-pairwise-Caste_Type1.csv", header= TRUE)
Alpha_plots_c$Caste_Type = factor(Alpha_plots_c$Caste_Type, levels=c("Q", "S", "W")) 
psc = ggboxplot(Alpha_plots_c, x = "Caste_Type", y= "shannon", color = "Caste_Type") + scale_x_discrete(labels = c("Queen", "Soldier", "Worker"))
psc= psc + ggtitle("Mean Shannon Diversity") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none") + xlab("Caste")
psc = psc + geom_boxplot(fill = c("darkblue", "darkorchid", "darkseagreen")) + xlab("Caste")
psc

pdc = ggboxplot(Alpha_plots_c, x = "Caste_Type", y= "faith_pd", color = "Caste_Type") + scale_x_discrete(labels = c("Queen", "Soldier", "Worker"))
pdc = pdc  + scale_colour_manual(values = c("darkblue", "darkorchid", "darkseagreen", "darkolivegreen1","darkgoldenrod1")) + ggtitle("Mean Faith's PD") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none") 
pdc = pdc + geom_boxplot(fill = c("darkblue", "darkorchid", "darkseagreen")) + xlab("Caste")
pdc

pec = ggboxplot(Alpha_plots, x = "Caste_Type", y= "pielou_e", color = "Caste_Type") + scale_x_discrete(labels = c("Queen", "Soldier", "Worker"))
pec = pec + scale_colour_manual(values = c("darkblue", "darkorchid", "darkseagreen", "darkolivegreen1","darkgoldenrod1")) + ggtitle("Mean Pielou's Evenness") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none") 
pec = pec + geom_boxplot(fill = c("darkblue", "darkorchid", "darkseagreen")) + xlab("Caste")
pec

pac = ggboxplot(Alpha_plots_c, x = "Caste_Type", y= "mean_ASV_count", color = "Caste_Type") + scale_x_discrete(labels = c("Queen", "Soldier", "Worker"))
pac = pac + scale_colour_manual(values = c("darkblue", "darkorchid", "darkseagreen", "darkolivegreen1","darkgoldenrod1")) + ggtitle("Mean ASV Count") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none") 
pac = pac + geom_boxplot(fill = c("darkblue", "darkorchid", "darkseagreen")) + xlab("Caste")
pac

grid.arrange(psc, pdc, pec, pac)


#### Beta Diversity Metrics ####

#bray-curtis all samples no pro #rooted
plot.brayr = phyloseq::ordinate(physeq_rarefied_nopro_qiime_rooted, "PCoA", "bray")
pbr = phyloseq::plot_ordination(physeq_rarefied_nopro_qiime_rooted, plot.brayr, color = "Species") 
#p = p + scale_colour_brewer(type="qual", palette="Set2")
pbr = pbr + scale_colour_manual(values = c("deepskyblue1", "green", "red3", "darkorange1", "darkviolet", "khaki1", "blue", "red1", "plum1", "goldenrod2", "green4"))
pbr + geom_point() + stat_ellipse()
#pdf is 10x10

#bray curtis by gut compartment without procryptocerus
#bray-curtis just gaster #rooted
physeq_rarefied_nopro_qiime_rooted_justgaster1 <- phyloseq::subset_samples(physeq_rarefied_nopro_qiime_rooted, Body_Type == "G")
plot.bray = phyloseq::ordinate(physeq_rarefied_nopro_qiime_rooted_justgaster1, "PCoA", "bray")
pbg = phyloseq::plot_ordination(physeq_rarefied_nopro_qiime_rooted_justgaster1, plot.bray, color = "Species") 
#p = p + scale_colour_brewer(type="qual", palette="Set2")
pbg = pbg + scale_colour_manual(values = c("deepskyblue1", "green", "red3", "darkorange1", "darkviolet", "khaki1", "blue", "red1", "plum1", "goldenrod2", "green4", "black"))
pbg = pbg + geom_point() + stat_ellipse() + ggtitle("Gaster") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size =22), axis.title.x = element_text(size=12, face = "bold"), axis.title.y = element_text(size=12, face = "bold")) + theme(legend.position="none") + guides(fill=guide_legend(nrow=5)) + geom_point(size =3)
pbg
#bray-curtis just midgut #rooted
physeq_rarefied_nopro_qiime_rooted_justmidgut <- phyloseq::subset_samples(physeq_rarefied_nopro_qiime_rooted, Body_Type == "M")
plot.bray = phyloseq::ordinate(physeq_rarefied_nopro_qiime_rooted_justmidgut, "PCoA", "bray")
pbgm = phyloseq::plot_ordination(physeq_rarefied_nopro_qiime_rooted_justmidgut, plot.bray, color = "Species") 
#p = p + scale_colour_brewer(type="qual", palette="Set2")
pbgm = pbgm + scale_colour_manual(values = c("deepskyblue1", "green", "red3", "darkorange1", "darkviolet", "khaki1", "blue", "red1", "plum1", "goldenrod2", "green4", "black"))
pbgm = pbgm + geom_point() + stat_ellipse() + ggtitle("Midgut") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size =22), axis.title.x = element_text(size=12, face = "bold"), axis.title.y = element_text(size=12, face = "bold")) + theme(legend.position="none") + guides(fill=guide_legend(nrow=5)) + geom_point(size =3)
pbgm
#bray-curtis just crop #rooted
physeq_rarefied_nopro_qiime_rooted_justcrop <- subset_samples(physeq_rarefied_nopro_qiime_rooted, Body_Type == "C")
plot.bray = phyloseq::ordinate(physeq_rarefied_nopro_qiime_rooted_justcrop, "PCoA", "bray")
pbgc = phyloseq::plot_ordination(physeq_rarefied_nopro_qiime_rooted_justcrop, plot.bray, color = "Species") 
#p = p + scale_colour_brewer(type="qual", palette="Set2")
pbgc = pbgc + scale_colour_manual(values = c("deepskyblue1", "green", "red3", "darkorange1", "darkviolet", "khaki1", "blue", "red1", "plum1", "goldenrod2", "green4", "black"))
pbgc = pbgc + geom_point() + stat_ellipse() + ggtitle("Crop") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size =22), axis.title.x = element_text(size=12, face = "bold"), axis.title.y = element_text(size=12, face = "bold")) + theme(legend.position="none") + guides(fill=guide_legend(nrow=5)) + geom_point(size =3)
pbgc


#bray-curtis just ileum #rooted #with phylogeny coloring
physeq_rarefied_nopro_qiime_rooted_justileum <- subset_samples(physeq_rarefied_nopro_qiime_rooted, Body_Type == "I")
plot.bray = phyloseq::ordinate(physeq_rarefied_nopro_qiime_rooted_justileum, "PCoA", "bray")
pbgi = phyloseq::plot_ordination(physeq_rarefied_nopro_qiime_rooted_justileum, plot.bray, color = "Species") 
#p = p + scale_colour_brewer(type="qual", palette="Set2")
pbgi = pbgi + scale_colour_manual(values = c("deepskyblue1", "green", "red3", "darkorange1", "darkviolet", "khaki2", "blue", "red1", "plum1", "goldenrod2", "green4")) + ggtitle("Ileum") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size =15))
pbgi = pbgi + geom_point() + stat_ellipse() + ggtitle("Ileum") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size =22), axis.title.x = element_text(size=12, face = "bold"), axis.title.y = element_text(size=12, face = "bold")) + theme(legend.position="center")  + geom_point(size =3)
pbgi 

#bray-curtis just rectum #rooted
physeq_rarefied_nopro_qiime_rooted_justrectum <- subset_samples(physeq_rarefied_nopro_qiime_rooted, Body_Type == "R")
plot.bray =phyloseq:: ordinate(physeq_rarefied_nopro_qiime_rooted_justrectum, "PCoA", "bray")
pbgr = phyloseq::plot_ordination(physeq_rarefied_nopro_qiime_rooted_justrectum, plot.bray, color = "Species") 
#p = p + scale_colour_brewer(type="qual", palette="Set2")
pbgr = pbgr + scale_colour_manual(values = c("deepskyblue1", "green", "red3", "darkorange1", "darkviolet", "khaki1", "blue", "red1", "plum1", "goldenrod2", "green4"))
pbgr = pbgr + geom_point() + stat_ellipse() + ggtitle("Rectum") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size =22), axis.title.x = element_text(size=12, face = "bold"), axis.title.y = element_text(size=12, face = "bold")) + theme(legend.position="none") + guides(fill=guide_legend(nrow=5)) + geom_point(size =3) 
pbgr 

grid.arrange(pbg, pbgc, pbgm, pbgi, pbgr)
#pdf 20x27

#LEGEND
physeq_rarefied_withpro_qiime_rooted_justrectum <- subset_samples(physeq_rarefied_withpro_qiime_rooted, Body_Type == "R")
plot.bray = phyloseq::ordinate(physeq_rarefied_withpro_qiime_rooted_justrectum, "PCoA", "bray")
pbgrl = phyloseq::plot_ordination(physeq_rarefied_withpro_qiime_rooted_justrectum, plot.bray, color = "Species") 
#p = p + scale_colour_brewer(type="qual", palette="Set2")
pbgrl = pbgrl + scale_colour_manual(values = c("deepskyblue1", "green", "red3", "darkorange1", "darkviolet", "khaki1", "blue", "red1", "plum1", "goldenrod2", "green4", "black"))
legendbray = pbgrl + geom_point() + stat_ellipse() + ggtitle("Rectum") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size =22), axis.text.x = element_text(size=6)) + theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))
legendbray

#bray-curtis by gut section with procryptocerus 
#bray-curtis just gaster #rooted
physeq_rarefied_withpro_qiime_rooted_justgaster <- phyloseq::subset_samples(physeq_rarefied_withpro_qiime_rooted, Body_Type == "G")
plot.bray = phyloseq::ordinate(physeq_rarefied_withpro_qiime_rooted_justgaster, "PCoA", "jaccard")
pbg = phyloseq::plot_ordination(physeq_rarefied_withpro_qiime_rooted_justgaster, plot.bray, color = "Species") 
#p = p + scale_colour_brewer(type="qual", palette="Set2")
pbg = pbg + scale_colour_manual(values = c("deepskyblue1", "green", "red3", "darkorange1", "darkviolet", "khaki1", "blue", "red1", "plum1", "goldenrod2", "green4", "grey59"))
pbg = pbg + geom_point() + stat_ellipse() + ggtitle("Gaster") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size =22)) + theme(legend.position="none") + guides(fill=guide_legend(nrow=5)) + geom_point(size =3)
pbg
#bray-curtis just midgut #rooted
physeq_rarefied_withpro_qiime_rooted_justmidgut <- subset_samples(physeq_rarefied_withpro_qiime_rooted, Body_Type == "M")
plot.bray = phyloseq::ordinate(physeq_rarefied_withpro_qiime_rooted_justmidgut, "PCoA", "bray")
pbgm = phyloseq::plot_ordination(physeq_rarefied_withpro_qiime_rooted_justmidgut, plot.bray, color = "Species") 
#p = p + scale_colour_brewer(type="qual", palette="Set2")
pbgm = pbgm + scale_colour_manual(values = c("deepskyblue1", "green", "red3", "darkorange1", "darkviolet", "khaki1", "blue", "red1", "plum1", "goldenrod2", "green4", "grey70"))
pbgm = pbgm + geom_point() + stat_ellipse() + ggtitle("Midgut") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size =22)) + theme(legend.position="none") + guides(fill=guide_legend(nrow=5)) + geom_point(size =3)
pbgm
#bray-curtis just crop #rooted
physeq_rarefied_withpro_qiime_rooted_justcrop <- subset_samples(physeq_rarefied_withpro_qiime_rooted, Body_Type == "C")
plot.bray = phyloseq::ordinate(physeq_rarefied_withpro_qiime_rooted_justcrop, "PCoA", "bray")
pbgc = phyloseq::plot_ordination(physeq_rarefied_withpro_qiime_rooted_justcrop, plot.bray, color = "Species") 
#p = p + scale_colour_brewer(type="qual", palette="Set2")
pbgc = pbgc + scale_colour_manual(values = c("deepskyblue1", "green", "red3", "darkorange1", "darkviolet", "khaki1", "blue", "red1", "plum1", "goldenrod2", "green4", "grey59"))
pbgc = pbgc + geom_point() + stat_ellipse() + ggtitle("Crop") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size =20)) + theme(legend.position="none") + guides(fill=guide_legend(nrow=5)) + geom_point(size =3)
pbgc
#bray-curtis just ileum #rooted
phyloseq::physeq_rarefied_withpro_qiime_rooted_justileum <- subset_samples(physeq_rarefied_withpro_qiime_rooted, Body_Type == "I")
plot.bray = phyloseq::ordinate(physeq_rarefied_withpro_qiime_rooted_justileum, "PCoA", "bray")
pbgi = phyloseq::plot_ordination(physeq_rarefied_withpro_qiime_rooted_justileum, plot.bray, color = "Species") 
#p = p + scale_colour_brewer(type="qual", palette="Set2")
pbgi = pbgi + scale_colour_manual(values = c("deepskyblue1", "green", "red3", "darkorange1", "darkviolet", "khaki1", "blue", "red1", "plum1", "goldenrod2", "green4", "grey59")) + ggtitle("Ileum") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size =15))
pbgi = pbgi + geom_point() + stat_ellipse() + ggtitle("Ileum") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size =22)) + theme(legend.position="center")  + geom_point(size =3)
pbgi 


#bray-curtis just ileum #rooted #with phylogeny coloring
physeq_rarefied_nopro_qiime_rooted_justileum <- subset_samples(physeq_rarefied_withpro_qiime_rooted, Body_Type == "I")
plot.bray = phyloseq::ordinate(physeq_rarefied_nopro_qiime_rooted_justileum, "PCoA", "bray")
pbgi = phyloseq::plot_ordination(physeq_rarefied_nopro_qiime_rooted_justileum, plot.bray, color = "Species") 
#p = p + scale_colour_brewer(type="qual", palette="Set2")
pbgi = pbgi + scale_colour_manual(values = c("deepskyblue1", "green", "red3", "darkorange1", "darkviolet", "khaki2", "blue", "red1", "plum1", "goldenrod2", "green4", "grey59")) + ggtitle("Ileum") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size =15))
pbgi = pbgi + geom_point() + stat_ellipse() + ggtitle("Ileum") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size =22)) + theme(legend.position="center")  + geom_point(size =3)
pbgi 

#bray-curtis just rectum #rooted
physeq_rarefied_withpro_qiime_rooted_justrectum <- subset_samples(physeq_rarefied_withpro_qiime_rooted, Body_Type == "R")
plot.bray =phyloseq:: ordinate(physeq_rarefied_withpro_qiime_rooted_justrectum, "PCoA", "bray")
pbgr = phyloseq::plot_ordination(physeq_rarefied_withpro_qiime_rooted_justrectum, plot.bray, color = "Species") 
#p = p + scale_colour_brewer(type="qual", palette="Set2")
pbgr = pbgr + scale_colour_manual(values = c("deepskyblue1", "green", "red3", "darkorange1", "darkviolet", "khaki1", "blue", "red1", "plum1", "goldenrod2", "green4", "grey59"))
pbgr = pbgr + geom_point() + stat_ellipse() + ggtitle("Rectum") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size =22)) + theme(legend.position="none") + guides(fill=guide_legend(nrow=5)) + geom_point(size =3) 
pbgr 

grid.arrange(pbg, pbgc, pbgm, pbgi, pbgr)
#pdf 20x27

#LEGEND
physeq_rarefied_withpro_qiime_rooted_justrectum <- subset_samples(physeq_rarefied_withpro_qiime_rooted, Body_Type == "R")
plot.bray = phyloseq::ordinate(physeq_rarefied_withpro_qiime_rooted_justrectum, "PCoA", "bray")
pbgrl = phyloseq::plot_ordination(physeq_rarefied_withpro_qiime_rooted_justrectum, plot.bray, color = "Species") 
#p = p + scale_colour_brewer(type="qual", palette="Set2")
pbgrl = pbgrl + scale_colour_manual(values = c("deepskyblue1", "green", "red3", "darkorange1", "darkviolet", "khaki1", "blue", "red1", "plum1", "goldenrod2", "green4", "grey59"))
legendbray = pbgrl + geom_point() + stat_ellipse() + ggtitle("Rectum") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size =15))  + theme(legend.text=element_text(size=10, face = "bold")) + guides(fill=guide_legend(nrow=5, byrow= TRUE) + theme(legend.position="bottom"))
legendbray


#wunifrac all samples (including gaster) #rooted
plot.wunifracnp = ordinate(physeq_rarefied_nopro_qiime_rooted, "PCoA", "wunifrac")
pwnp = plot_ordination(physeq_rarefied_nopro_qiime_rooted, plot.wunifracnp, color = "Body_Type") 
#p = p + scale_colour_brewer(type="qual", palette="Set2")
pwnp = pwnp + scale_colour_manual(values = c("darkblue", "darkgoldenrod1", "green4", "darkorchid", "cadetblue2", "lightskyblue", "deeppink", "firebrick", "khaki2", "darkgreen", "brown1", "darkorange1", "cyan1", "darksalmon", "royalblue4"))
pwnp = pwnp +  geom_point() + stat_ellipse()
pwnp
#wunifrac without gaster samples added #rooted #looks awesome
physeq_rarefied_nopro_qiime_rooted_nogaster <- subset_samples(physeq_rarefied_nopro_qiime_rooted, Body_Type != "G")
plot.wunifracnpng = ordinate(physeq_rarefied_nopro_qiime_rooted_nogaster, "PCoA", "wunifrac")
pwngnp = plot_ordination(physeq_rarefied_nopro_qiime_rooted_nogaster, plot.wunifracnpng, color = "Body_Type") 
pwngnp = pwngnp + scale_colour_manual(values = c("darkblue", "green4", "darkorchid", "cadetblue2", "lightskyblue", "deeppink", "darkolivegreen1", "khaki2", "darkgreen", "brown1", "darkorange1", "cyan1", "darksalmon", "royalblue4"))
pwngnp = pwngnp +  geom_point() + stat_ellipse() 

grid.arrange(pwngnp, pwnp, nrow = 1)

###### SIMPER boxplots ####
library(ggpubr)

SIMPER_data_rare <- read.csv("physeq_rarefied_nopro_qiime_rooted_SIMPER.csv", header= TRUE)
SIMPER_data_rare_withpro <- read.csv("physeq_rarefied_withpro_qiime_rooted_SIMPER_transpo.csv", header= TRUE)

# SIMPER plot for "b395480d9912c42a0fd8f73dbcec81dc"
SIMPER_data_rare$Body_Type = factor(SIMPER_data_rare$Body_Type, levels=c("C", "M", "I", "R", "G")) 
p = ggboxplot(SIMPER_data_rare, x = "Body_Type", y= "b395480d9912c42a0fd8f73dbcec81dc", color = "Body_Type", ylim = c(0,13000)) + scale_x_discrete(labels = c("Crop", "Midgut", "Ileum", "Rectum", "Gaster"))
p = p + geom_boxplot(fill = c("darkblue", "darkorchid", "darkseagreen", "darkolivegreen1", "darkgoldenrod1")) + ggtitle("Opitutales ASV2 by gut compartment") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none") + ylab("Read Number")
p

SIMPER_data_rare$Body_Type = factor(SIMPER_data_rare$Body_Type, levels=c("C", "M", "I", "R", "G")) 
SIMPER_data_rare_texanus <- subset(SIMPER_data_rare, Species == "C_texanus")
ggg <- ggboxplot(SIMPER_data_rare_texanus, x = "Body_Type", y= "b395480d9912c42a0fd8f73dbcec81dc", ylim = c(0,13000)) +
  scale_x_discrete(labels = c("Crop", "Midgut", "Ileum", "Rectum", "Gaster")) + 
  geom_boxplot(fill = c("darkblue", "darkorchid", "green4", "cadetblue2", "darkgoldenrod1")) + 
  ggtitle("Opitutales ASV2 by gut compartment for C. texanus") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none") +
  ylab("Read Number for C. texanus samples") + xlab("Gut Compartment")
ggg
pairwise.wilcox.test(SIMPER_data_rare_texanus$b395480d9912c42a0fd8f73dbcec81dc, SIMPER_data_rare_texanus$Body_Type, p.adj = "fdr")

SIMPER_data_rare_withpro$Species <- factor(SIMPER_data_rare_withpro$Species, levels=c("C_atratus", "C_auricomus", "C_cordatus", "C_grandinosus", "C_minutus", "C_multispinosus", "C_placidus",  "C_setulifer", "C_simillimus",  "C_texanus",  "C_unimaculatus",  "Procryptocerus_sp"))
ggB <- ggboxplot(SIMPER_data_rare_withpro, x = "Species", y= "b395480d9912c42a0fd8f73dbcec81dc", color = "Species", ylim = c(0,15000)) + ggtitle("Opitutales ASV2 by species") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none") + ylab("Read Number") 
ggB <- ggB + scale_x_discrete(labels=c("C_atratus" = "atratus", "C_auricomus" = "auricomus", "C_cordatus" = "cordatus", "C_grandinosus" = "grandinosus", "C_minutus" = "minutus", "C_multispinosus" = "multispinosus", "C_placidus" = "placidus", "C_setulifer" = "setulifer", "C_simillimus" = "simillimus", "C_texanus" = "texanus", "C_unimaculatus" = "unimaculatus", "Procryptocerus_sp" = "Procryptocerus")) +
  geom_boxplot(fill = c("goldenrod2", "green", "red3", "darkorange1", "darkviolet", "khaki1", "blue", "red1", "plum1", "goldenrod2", "green4", "black"))
ggB
pairwise.wilcox.test(SIMPER_data_rare_withpro$b395480d9912c42a0fd8f73dbcec81dc, SIMPER_data_rare_withpro$Species,
                     p.adj = "fdr")

grid.arrange(ggg, ggB)
#pdf is 10x15

# SIMPER plot for "X19fec32f5ac2990cd5512922f099546a"
SIMPER_data_rare$Body_Type = factor(SIMPER_data_rare$Body_Type, levels=c("C", "M", "I", "R", "G")) 
p = ggboxplot(SIMPER_data_rare, x = "Body_Type", y= "X19fec32f5ac2990cd5512922f099546a", color = "Body_Type", ylim = c(0,13000)) + scale_x_discrete(labels = c("Crop", "Midgut", "Ileum", "Rectum", "Gaster"))
p = p + geom_boxplot(fill = c("darkblue", "darkorchid", "green4", "cadetblue2", "darkgoldenrod1")) +
  ggtitle("Opitutales ASV1 by gut compartment") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none") +
  ylab("Read Number") + 
  xlab("Gut Compartment")
p
pairwise.wilcox.test(SIMPER_data_rare$X19fec32f5ac2990cd5512922f099546a, SIMPER_data_rare$Body_Type, p.adj = "fdr")

SIMPER_data_rare_withpro$Species <- factor(SIMPER_data_rare_withpro$Species, levels=c("C_atratus", "C_auricomus", "C_cordatus", "C_grandinosus", "C_minutus", "C_multispinosus", "C_placidus",  "C_setulifer", "C_simillimus",  "C_texanus",  "C_unimaculatus",  "Procryptocerus_sp"))
ggx <- ggboxplot(SIMPER_data_rare_withpro, x = "Species", y= "X19fec32f5ac2990cd5512922f099546a", color = "Species", ylim = c(0,15000)) + ggtitle("Opitutales ASV1 by species") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none") + ylab("Read Number") 
ggx <- ggx + scale_x_discrete(labels=c("C_atratus" = "atratus", "C_auricomus" = "auricomus", "C_cordatus" = "cordatus", "C_grandinosus" = "grandinosus", "C_minutus" = "minutus", "C_multispinosus" = "multispinosus", "C_placidus" = "placidus", "C_setulifer" = "setulifer", "C_simillimus" = "simillimus", "C_texanus" = "texanus", "C_unimaculatus" = "unimaculatus", "Procryptocerus_sp" = "Procryptocerus")) +
  geom_boxplot(fill = c("deepskyblue1", "green", "red3", "darkorange1", "darkviolet", "khaki1", "blue", "red1", "plum1", "goldenrod2", "green4", "gray")) 
ggx
pairwise.wilcox.test(SIMPER_data_rare_withpro$X19fec32f5ac2990cd5512922f099546a, SIMPER_data_rare_withpro$Species, p.adj = "fdr")

SIMPER_data_rare_M <- subset(SIMPER_data_rare_withpro, Body_Type == "M")
SIMPER_data_rare_M$Species <- factor(SIMPER_data_rare_M$Species, levels=c("C_atratus", "C_auricomus", "C_cordatus", "C_grandinosus", "C_minutus", "C_multispinosus", "C_placidus",  "C_setulifer", "C_simillimus",  "C_texanus",  "C_unimaculatus",  "Procryptocerus_sp"))
ggM <- ggboxplot(SIMPER_data_rare_M, x = "Species", y= "X19fec32f5ac2990cd5512922f099546a", ylim = c(0,15000)) + ggtitle("Opitutales ASV1 by species in Midgut") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none") + ylab("Read Number in Midgut") 
ggM <- ggM + scale_x_discrete(labels=c("C_atratus" = "atratus", "C_auricomus" = "auricomus", "C_cordatus" = "cordatus", "C_grandinosus" = "grandinosus", "C_minutus" = "minutus", "C_multispinosus" = "multispinosus", "C_placidus" = "placidus", "C_setulifer" = "setulifer", "C_simillimus" = "simillimus", "C_texanus" = "texanus", "C_unimaculatus" = "unimaculatus", "Procryptocerus_sp" = "Procryptocerus")) + 
  geom_boxplot(fill = c("deepskyblue1", "green", "red3", "darkorange1", "darkviolet", "khaki1", "blue", "red1", "plum1", "goldenrod2", "green4", "gray")) 
ggM
pairwise.wilcox.test(SIMPER_data_rare_M$X19fec32f5ac2990cd5512922f099546a, SIMPER_data_rare_M$Species, p.adj = "fdr")

grid.arrange(p, ggx, ggM)
#pdf is 10x15

# SIMPER plot for "d7f90b5fb735570e661ecfacbfec037a"
SIMPER_data_rare_withpro$Species <- factor(SIMPER_data_rare_withpro$Species, levels=c("C_atratus", "C_auricomus", "C_cordatus", "C_grandinosus", "C_minutus", "C_multispinosus", "C_placidus",  "C_setulifer", "C_simillimus",  "C_texanus",  "C_unimaculatus",  "Procryptocerus_sp"))
d7 <- ggboxplot(SIMPER_data_rare_withpro, x = "Species", y= "d7f90b5fb735570e661ecfacbfec037a", color = "Species", ylim = c(0,7000))  + ggtitle("Xanthomonadales ASV2 by species") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none") + ylab("Read Number") 
d7 <- d7 + scale_x_discrete(labels=c("C_atratus" = "atratus", "C_auricomus" = "auricomus", "C_cordatus" = "cordatus", "C_grandinosus" = "grandinosus", "C_minutus" = "minutus", "C_multispinosus" = "multispinosus", "C_placidus" = "placidus", "C_setulifer" = "setulifer", "C_simillimus" = "simillimus", "C_texanus" = "texanus", "C_unimaculatus" = "unimaculatus", "Procryptocerus_sp" = "Procryptocerus")) + 
  geom_boxplot(fill = c("deepskyblue1", "green", "red3", "darkorange1", "darkviolet", "khaki1", "blue", "red1", "plum1", "goldenrod2", "green4", "black")) 
d7
#pdf is 10x15
pairwise.wilcox.test(SIMPER_data_rare_withpro$d7f90b5fb735570e661ecfacbfec037a, SIMPER_data_rare_withpro$Species, p.adj = "fdr")



#### Cephalotes Phylogeny ####
library("phytools")
library("phylotools")
library("ape")
#original phylogeny input from Price et al. 2016 
Cephalotes_Tree <- read.nexus("cephalotes.tree")
Cephalotes_Tree
#choose the tips you want to keep
species <- c("SP44_atratus_atr", "MLA01_auricomus_ham", "13_cordatus_dep","GB15_grandinosus_grand","JSC06_minutus_lam", "29_multispinosus_multi","SP81_placidus_atr", "JSC05_sp_Pro", "SY04_texanus_tex", "JTL06_setulifer_coff","TRS06_simillimus_lam", "MI01_unimaculatus_ham")
#
new_tips <- c("C_atratus", "C_auricomus", "C_cordatus", "C_grandinosus", "C_minutus", "C_multispinosus", "C_placidus", "Procryptocerus_sp", "C_texanus", "C_setulifer", "C_simillimus", "C_unimaculatus")

new_mantel_tips <- c("atratus", "auricomus", "cordatus", "grandinosus", "minutus", "multispinosus", "placidus", "Procryptocerus", "texanus", "setulifer", "simillimus", "unimaculatus")

#this command will drop all tips that are not the ones selected
pruned_ceph.tree <- drop.tip(Cephalotes_Tree, Cephalotes_Tree$tip.label[-match(species, Cephalotes_Tree$tip.label)])
dat <- data.frame(species, new_tips)
ntree <- sub.taxa.label(pruned_ceph.tree, dat)
ntree
#writes this phylogeny to a tree file
write.tree(ntree)
#writes this tree file to a nexus file
write.tree(ntree, file ="cephalotes1.tre")
ntree$edge.length<-NULL
write.tree(ntree, file = "pruned_ceph_tree_good1.tre")
#####Mantel Tests #####
library(vegan)
library(ecodist)
library(ape)

####Cephalotes Host Phylogenies###

#cephalotes trees#
ceph <- read.tree(file = "cephalotes.tre")
PatristicDistMatrix <-cophenetic(ceph)
PatristicDistMatrix1 <- as.dist(PatristicDistMatrix)
cephtree <-PatristicDistMatrix1 
cephtree 

#cephalotes crop: no cordatus
cephC <- drop.tip(ceph, "C_cordatus" )
cephC <-cophenetic(cephC)
cephC <- as.dist(cephC)
cephC

tips <- c("C_placidus", "C_auricomus", "C_cordatus")
cephRhC <- drop.tip(ceph, tips)
cephRhC <-cophenetic(cephRhC)
cephRhC <- as.dist(cephRhC)
cephRhC

#cephalotes enterobacteriales crop
tips <- c("C_cordatus", "C_auricomus", "C_grandinosus", "C_unimaculatus")
cephE <- drop.tip(ceph, tips)
cephE <-cophenetic(cephE)
cephE <- as.dist(cephE)
cephE

tips <- c("C_cordatus", "C_multispinosus", "C_unimaculatus")
cephEC <- drop.tip(ceph, tips)
cephEC <-cophenetic(cephEC)
cephEC <- as.dist(cephEC)
cephEC

tips <- c("C_cordatus", "C_placidus", "C_grandinosus", "C_unimaculatus")
cephEI <- drop.tip(ceph, tips)
cephEI <-cophenetic(cephEI)
cephEI <- as.dist(cephEI)
cephEI

tips <- c("C_cordatus", "C_placidus", "C_unimaculatus")
cephER <- drop.tip(ceph, tips)
cephER <-cophenetic(cephER)
cephER <- as.dist(cephER)
cephER

tips <- c("C_atratus", "C_multispinosus", "C_setulifer", "C_texanus", "C_grandinosus", "C_unimaculatus")
cephR <- drop.tip(ceph, tips)
cephR <-cophenetic(cephR)
cephR <- as.dist(cephR)
cephR

tips <- c("C_atratus", "C_grandinosus")
cephR <- drop.tip(ceph, tips)
cephR <-cophenetic(cephR)
cephR <- as.dist(cephR)
cephR

tips <- c("C_multispinosus", "C_unimaculatus")
cephE <- drop.tip(ceph, tips)
cephE <-cophenetic(cephE)
cephE <- as.dist(cephE)
cephE

tips <- c("C_grandinosus", "C_unimaculatus", "C_simillimus", "C_setulifer")
cephRG <- drop.tip(ceph, tips)
cephRG <-cophenetic(cephRG)
cephRG <- as.dist(cephRG)
cephRG

tips <- c("C_texanus", "C_atratus", "C_placidus", "C_unimaculatus", "C_grandinosus")
cephRM <- drop.tip(ceph, tips)
cephRM <-cophenetic(cephRM)
cephRM <- as.dist(cephRM)
cephRM

tips <- c("C_auricomus", "C_cordatus", "C_placidus", "C_unimaculatus", "C_grandinosus", "C_minutus", "C_setulifer", "C_multispinosus")
cephEM <- drop.tip(ceph, tips)
cephEM <-cophenetic(cephEM)
cephEM <- as.dist(cephEM)
cephEM


tips <- c("Procryptocerus_sp", "C_cordatus", "C_placidus", "C_grandinosus", "C_setulifer", "C_multispinosus")
cephEG <- drop.tip(ceph, tips)
cephEG <-cophenetic(cephEG)
cephEG <- as.dist(cephEG)
cephEG

tips <- c("C_atratus", "C_multispinosus", "C_grandinosus")
cephRR <- drop.tip(ceph, tips)
cephRR <-cophenetic(cephRR)
cephRR <- as.dist(cephRR)
cephRR

tips <- c("C_atratus", "C_cordatus", "C_grandinosus")
cephRC <- drop.tip(ceph, tips)
cephRC <-cophenetic(cephRC)
cephRC <- as.dist(cephRC)
cephRC

tips <- c("C_atratus", "C_multispinosus", "C_grandinosus")
cephRR <- drop.tip(ceph, tips)
cephRR <-cophenetic(cephRR)
cephRR <- as.dist(cephRR)
cephRR

tips <- c("C_texanus", "C_atratus", "C_multispinosus", "C_unimaculatus", "C_simillimus", "C_setulifer")
cephRI <- drop.tip(ceph, tips)
cephRI <-cophenetic(cephRI)
cephRI <- as.dist(cephRI)
cephRI

tips <- c("C_grandinosus", "C_unimaculatus", "C_simillimus", "C_setulifer")
cephRG <- drop.tip(ceph, tips)
cephRG <-cophenetic(cephRG)
cephRG <- as.dist(cephRG)
cephRG

###MANTEL TESTS for different gut compartments from 16S rRNA amplicon sequencing###

#overall
merged<- merge_samples(physeq_rarefied_withpro_qiime_rooted, "Species")
otu_mean <- as(otu_table(merged), "matrix")
otu_mean_dist_T <- distance(otu_mean, method = "bray")
vegan::mantel(otu_mean_dist_T, cephtree, method = "pearson", permutations = 1000000)

#midgut
merged_Midgut <- subset_samples(physeq_rarefied_withpro_qiime_rooted, Body_Type == "M")
merged_M_mean<- merge_samples(merged_Midgut, "Species")
M_otu_mean <- as(otu_table(merged_M_mean), "matrix")
M_otu_mean_dist_T <- distance(M_otu_mean, method = "bray")
vegan::mantel(M_otu_mean_dist_T, cephtree, method = "pearson", permutations = 1000000)

#ileum
merged_Ileum <- subset_samples(physeq_rarefied_withpro_qiime_rooted, Body_Type == "I")
merged_withproI <- merge_samples(merged_Ileum, "Species")
tryI_OTU <- as(otu_table(merged_withproI), "matrix")
tryI_OTU_T <- t(tryI_OTU)
tryI.dist_OTU <- distance(tryI_OTU, method = "bray")
vegan::mantel(tryI.dist_OTU, cephtree, method = "pearson", permutations = 1000000)

#rectum
merged_Rectum <- subset_samples(physeq_rarefied_withpro_qiime_rooted, Body_Type == "R")
merged_withproR <- merge_samples(merged_Rectum, "Species")
R_otu_mean <- as(otu_table(merged_withproR), "matrix")
tryR.dist_OTU <- distance(R_otu_mean, method = "bray")
vegan::mantel(tryR.dist_OTU, cephtree, method = "pearson", permutations = 1000000)

#gaster
merged_Gaster <- subset_samples(physeq_rarefied_withpro_qiime_rooted, Body_Type == "G")
merged_withproG <- merge_samples(merged_Gaster, "Species")
tryG_OTU <- as(otu_table(merged_withproG), "matrix")
tryG.dist_OTU <- distance(tryG_OTU, method = "bray")
vegan::mantel(tryG.dist_OTU, cephtree, method = "pearson", permutations = 1000000)

#crop
merged_Crop <- subset_samples(physeq_rarefied_withpro_qiime_rooted, Body_Type == "C")
merged_withproC <- merge_samples(merged_Crop, "Species")
tryC_OTU <- as(otu_table(merged_withproC), "matrix")
tryC.dist_OTU <- distance(tryC_OTU, method = "bray")
vegan::mantel(tryC.dist_OTU, cephC, method = "pearson", permutations = 1000000)

###### Betaproteobacteriales ###
#betaproteobacteriales
merged_beta <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Betaproteobacteriales")
merged_beta_mean<- merge_samples(merged_beta, "Species")
beta_otu_mean <- as(otu_table(merged_beta_mean), "matrix")
beta_otu_mean_dist_T <- distance(beta_otu_mean, method = "bray")
vegan::mantel(beta_otu_mean_dist_T, cephtree, method = "pearson", permutations = 1000000)

#betaproteobacteriales and crop
merged_beta <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Betaproteobacteriales")
merged_beta_C <- subset_samples(merged_beta, Body_Type == "C")
merged_beta_C_mean<- merge_samples(merged_beta_C, "Species")
beta_otu_C_mean <- as(otu_table(merged_beta_C_mean), "matrix")
rowSums(beta_otu_C_mean)
beta_otu_C_mean_dist <- distance(beta_otu_C_mean, method = "bray")
vegan::mantel(beta_otu_C_mean_dist, cephC, method = "pearson", permutations =1000000)


#betaproteobacteriales and midgut
merged_beta <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Betaproteobacteriales")
merged_beta_M <- subset_samples(merged_beta, Body_Type == "M")
merged_beta_M_mean<- merge_samples(merged_beta_M, "Species")
beta_otu_M_mean <- as(otu_table(merged_beta_M_mean), "matrix")
beta_otu_M_mean_dist <- distance(beta_otu_M_mean, method = "bray")
vegan::mantel(beta_otu_M_mean_dist, cephtree, method = "pearson")

#betaproteobacteriales and ileum
merged_beta <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Betaproteobacteriales")
merged_beta_I <- subset_samples(merged_beta, Body_Type == "I")
merged_beta_I_mean<- merge_samples(merged_beta_I, "Species")
beta_otu_I_mean <- as(otu_table(merged_beta_I_mean), "matrix")
beta_otu_I_mean_dist <- distance(beta_otu_I_mean, method = "bray")
vegan::mantel(beta_otu_I_mean_dist, cephtree, method = "pearson", permutations = 100000)

#betaproteobacteriales and rectum
merged_beta <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Betaproteobacteriales")
merged_beta_R <- subset_samples(merged_beta, Body_Type == "R")
merged_beta_R_mean<- merge_samples(merged_beta_R, "Species")
beta_otu_R_mean <- as(otu_table(merged_beta_R_mean), "matrix")
beta_otu_R_mean_dist <- distance(beta_otu_R_mean, method = "bray")
vegan::mantel(beta_otu_R_mean_dist, cephtree, method = "pearson", permutations = 1000000)

#betaproteobacteriales and gaster
merged_beta <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Betaproteobacteriales")
merged_beta_G <- subset_samples(merged_beta, Body_Type == "G")
merged_beta_G_mean<- merge_samples(merged_beta_G, "Species")
beta_otu_G_mean <- as(otu_table(merged_beta_G_mean), "matrix")
beta_otu_G_mean_dist <- distance(beta_otu_G_mean, method = "bray")
vegan::mantel(beta_otu_G_mean_dist, cephtree, method = "pearson", permutations = 1000000)


###### Rhizobiales ###
#rhizobiales
merged_Rhizo <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Rhizobiales")
merged_Rhizo_mean<- merge_samples(merged_Rhizo, "Species")
Rhizo_otu_mean <- as(otu_table(merged_Rhizo_mean), "matrix")
Rhizo_otu_mean_dist <- distance(Rhizo_otu_mean, method = "bray")
vegan::mantel(Rhizo_otu_mean_dist, cephtree, method = "pearson", permutations = 1000000)

#rhizobiales and crop
merged_Rhizo <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Rhizobiales")
merged_Rhizo_Crop <- subset_samples(merged_Rhizo, Body_Type == "C")
merged_Rhizo_meanC <- merge_samples(merged_Rhizo_Crop, "Species")
Rhizo_otu_meanC <- as(otu_table(merged_Rhizo_meanC), "matrix")
rowSums(Rhizo_otu_meanC)
Rhizo_otu_mean_TC <- Rhizo_otu_meanC[-c(2, 6), ]
Rhizo_otu_mean_dist_TC <- distance(Rhizo_otu_mean_TC, method = "bray")
vegan::mantel(Rhizo_otu_mean_dist_TC, cephRhC, method = "pearson", permutations = 1000000)

#rhizobiales and midgut
merged_Rhizo <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Rhizobiales")
merged_Rhizo_M <- subset_samples(merged_Rhizo, Body_Type == "M")
merged_Rhizo_meanM <- merge_samples(merged_Rhizo_M, "Species")
Rhizo_otu_meanM <- as(otu_table(merged_Rhizo_meanM), "matrix")
Rhizo_otu_mean_distM<- distance(Rhizo_otu_meanM, method = "bray")
vegan::mantel(Rhizo_otu_mean_distM, cephtree, method = "pearson", permutations = 1000000)

#rhizobiales and ileum
merged_Rhizo <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Rhizobiales")
merged_Rhizo_I <- subset_samples(merged_Rhizo, Body_Type == "I")
merged_Rhizo_meanI <- merge_samples(merged_Rhizo_I, "Species")
Rhizo_otu_meanI <- as(otu_table(merged_Rhizo_meanI), "matrix")
Rhizo_otu_mean_dist_TI <- distance(Rhizo_otu_meanI, method = "bray")
vegan::mantel(Rhizo_otu_mean_dist_TI, cephtree, method = "pearson", permutations = 1000000)

#rhizobiales and rectum
merged_Rhizo <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Rhizobiales")
merged_Rhizo_R <- subset_samples(merged_Rhizo, Body_Type == "R")
merged_Rhizo_meanR <- merge_samples(merged_Rhizo_R, "Species")
Rhizo_otu_meanR <- as(otu_table(merged_Rhizo_meanR), "matrix")
Rhizo_otu_mean_dist_TR <- distance(Rhizo_otu_meanR, method = "bray")
vegan::mantel(Rhizo_otu_mean_dist_TR, cephtree, method = "pearson", permutations = 1000000)

#rhizobiales and gaster
merged_Rhizo <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Rhizobiales")
merged_Rhizo_G <- subset_samples(merged_Rhizo, Body_Type == "G")
merged_Rhizo_meanG <- merge_samples(merged_Rhizo_G, "Species")
Rhizo_otu_meanG <- as(otu_table(merged_Rhizo_meanG), "matrix")
Rhizo_otu_mean_dist_TG <- distance(Rhizo_otu_meanG, method = "bray")
vegan::mantel(Rhizo_otu_mean_dist_TG, cephtree, method = "pearson", permutations = 1000000)

#####Opitutales ###

#opitutales
merged_Opit <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Opitutales")
merged_Opit_mean<- merge_samples(merged_Opit, "Species")
Opit_otu_mean <- as(otu_table(merged_Opit_mean), "matrix")
Opit_otu_mean_dist <- distance(Opit_otu_mean, method = "bray")
vegan::mantel(Opit_otu_mean_dist, cephtree, method = "pearson", permutations = 1000000)

#opitutales and crop
merged_Opit <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Opitutales")
merged_Opit_C <- subset_samples(merged_Opit, Body_Type == "C")
merged_Opit_meanC<- merge_samples(merged_Opit_C, "Species")
Opit_otu_meanC <- as(otu_table(merged_Opit_meanC), "matrix")
Opit_otu_mean_distC <- distance(Opit_otu_meanC, method = "bray")
vegan::mantel(Opit_otu_mean_distC, cephC, method = "pearson", permutations =1000000)

#opitutales and midgut
merged_Opit <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Opitutales")
merged_Opit_M <- subset_samples(merged_Opit, Body_Type == "M")
merged_Opit_meanM<- merge_samples(merged_Opit_M, "Species")
Opit_otu_mean_M <- as(otu_table(merged_Opit_meanM), "matrix")
Opit_otu_mean_dist_M <- distance(Opit_otu_mean_M, method = "bray")
vegan::mantel(Opit_otu_mean_dist_M, cephtree, method = "pearson", permutations =1000000)

#opitutales and ileum
merged_Opit <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Opitutales")
merged_Opit_I <- subset_samples(merged_Opit, Body_Type == "I")
merged_Opit_meanI<- merge_samples(merged_Opit_I, "Species")
Opit_otu_meanI <- as(otu_table(merged_Opit_meanI), "matrix")
Opit_otu_mean_dist_I <- distance(Opit_otu_meanI, method = "bray")
vegan::mantel(Opit_otu_mean_dist_I, cephtree, method = "pearson", permutations = 1000000)

#opitutales and rectum
merged_Opit <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Opitutales")
merged_Opit_R <- subset_samples(merged_Opit, Body_Type == "R")
merged_Opit_meanR<- merge_samples(merged_Opit_R, "Species")
Opit_otu_meanR <- as(otu_table(merged_Opit_meanR), "matrix")
Opit_otu_mean_dist_R <- distance(Opit_otu_meanR, method = "bray")
vegan::mantel(Opit_otu_mean_dist_R, cephtree, method = "pearson", permutations = 1000000)

#opitutales and gaster
merged_Opit <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Opitutales")
merged_Opit_G <- subset_samples(merged_Opit, Body_Type == "G")
merged_Opit_meanG <- merge_samples(merged_Opit_G, "Species")
Opit_otu_meanG <- as(otu_table(merged_Opit_meanG), "matrix")
Opit_otu_mean_dist_G <- distance(Opit_otu_meanG, method = "bray")
vegan::mantel(Opit_otu_mean_dist_G, cephtree, method = "pearson", permutations = 1000000)

##### Xanthomonadales ###

#xanthomonadales
merged_Xantho <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Xanthomonadales")
merged_Xantho_mean<- merge_samples(merged_Xantho, "Species")
Xantho_otu_mean <- as(otu_table(merged_Xantho_mean), "matrix")
Xantho_otu_mean_dist <- distance(Xantho_otu_mean, method = "bray")
vegan::mantel(Xantho_otu_mean_dist, cephtree, method = "pearson", permutations = 1000000)

#xanthomonadales and crop
merged_Xantho <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Xanthomonadales")
merged_Xantho_C <- subset_samples(merged_Xantho, Body_Type == "C")
merged_Xantho_meanC <- merge_samples(merged_Xantho_C, "Species")
Xantho_otu_meanC <- as(otu_table(merged_Xantho_meanC), "matrix")
Xantho_otu_mean_dist_C <- distance(Xantho_otu_meanC, method = "bray")
vegan::mantel(Xantho_otu_mean_dist_C, cephC, method = "pearson", permutations = 100000)

#xanthomonadales and midgut
merged_Xantho <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Xanthomonadales")
merged_Xantho_M <- subset_samples(merged_Xantho, Body_Type == "M")
merged_Xantho_meanM <- merge_samples(merged_Xantho_M, "Species")
Xantho_otu_meanM <- as(otu_table(merged_Xantho_meanM), "matrix")
Xantho_otu_mean_dist_TM <- distance(Xantho_otu_meanM, method = "bray")
vegan::mantel(Xantho_otu_mean_dist_TM, cephtree, method = "pearson", permutations = 100000)

#xanthomonadales and ileum
merged_Xantho <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Xanthomonadales")
merged_Xantho_I <- subset_samples(merged_Xantho, Body_Type == "I")
merged_Xantho_meanI <- merge_samples(merged_Xantho_I, "Species")
Xantho_otu_meanI<- as(otu_table(merged_Xantho_meanI), "matrix")
Xantho_otu_mean_dist_TI <- distance(Xantho_otu_meanI, method = "bray")
vegan::mantel(Xantho_otu_mean_dist_TI, cephtree, method = "pearson", permutations = 1000000)

#xanthomonadales and rectum
merged_Xantho <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Xanthomonadales")
merged_Xantho_R <- subset_samples(merged_Xantho, Body_Type == "R")
merged_Xantho_meanR <- merge_samples(merged_Xantho_R, "Species")
Xantho_otu_meanR<- as(otu_table(merged_Xantho_meanR), "matrix")
Xantho_otu_mean_dist_TR <- distance(Xantho_otu_meanR, method = "bray")
vegan::mantel(Xantho_otu_mean_dist_TR, cephtree, method = "pearson", permutations = 1000000)

#xanthomonadales and gaster
merged_Xantho <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Xanthomonadales")
merged_Xantho_G <- subset_samples(merged_Xantho, Body_Type == "G")
merged_Xantho_meanG <- merge_samples(merged_Xantho_G, "Species")
Xantho_otu_meanG <- as(otu_table(merged_Xantho_meanG), "matrix")
Xantho_otu_mean_dist_TG <- distance(Xantho_otu_meanG, method = "bray")
vegan::mantel(Xantho_otu_mean_dist_TG, cephtree, method = "pearson", permutations = 1000000)

##### Enterobacteriales ###

#enterobacteriales
merged_Entero <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Enterobacteriales")
merged_Entero_mean<- merge_samples(merged_Entero, "Species")
Entero_otu_mean <- as(otu_table(merged_Entero_mean), "matrix")
rowSums(Entero_otu_mean)
Entero_otu_mean_row <- Entero_otu_mean[-c(3), ]
Entero_otu_mean_dist <- distance(Entero_otu_mean_row, method = "bray")
vegan::mantel(Entero_otu_mean_dist, cephC, method = "pearson", permutations =1000000)

#enterobacteriales and crop
merged_Entero <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Enterobacteriales")
merged_Entero_C <- subset_samples(merged_Entero, Body_Type == "C")
merged_Entero_meanC <- merge_samples(merged_Entero_C, "Species")
Entero_otu_meanC <- as(otu_table(merged_Entero_meanC), "matrix")
rowSums(Entero_otu_meanC)
Entero_otu_mean_rowC <- Entero_otu_meanC[-c(5,10), ]
Entero_otu_mean_dist_TC <- distance(Entero_otu_mean_rowC, method = "bray")
vegan::mantel(Entero_otu_mean_dist_TC, cephEC, method = "pearson", permutations = 1000000)


#enterobacteriales and midgut
merged_Entero <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Enterobacteriales")
merged_Entero_M <- subset_samples(merged_Entero, Body_Type == "M")
merged_Entero_meanM <- merge_samples(merged_Entero_M, "Species")
Entero_otu_meanM <- as(otu_table(merged_Entero_meanM), "matrix")
rowSums(Entero_otu_meanM)
Entero_otu_mean_rowM <- Entero_otu_meanM[-c(2,3,4,5,6,7,8,11), ]
Entero_otu_mean_dist_TM <- distance(Entero_otu_mean_rowM, method = "bray")
vegan::mantel(Entero_otu_mean_dist_TM, cephEM, method = "pearson", permutations = 1000000)


#enterobacteriales and ileum
merged_Entero <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Enterobacteriales")
merged_Entero_I <- subset_samples(merged_Entero, Body_Type == "I")
merged_Entero_meanI <- merge_samples(merged_Entero_I, "Species")
Entero_otu_meanI <- as(otu_table(merged_Entero_meanI), "matrix")
rowSums(Entero_otu_meanI)
Entero_otu_mean_rowI <- Entero_otu_meanI[-c(3,4,7,11), ]
Entero_otu_mean_dist_TI <- distance(Entero_otu_mean_rowI, method = "bray")
vegan::mantel(Entero_otu_mean_dist_TI, cephEI, method = "pearson", permutations = 1000000)

#enterobacteriales and rectum
merged_Entero <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Enterobacteriales")
merged_Entero_R <- subset_samples(merged_Entero, Body_Type == "R")
merged_Entero_meanR <- merge_samples(merged_Entero_R, "Species")
Entero_otu_meanR <- as(otu_table(merged_Entero_meanR), "matrix")
rowSums(Entero_otu_meanR)
Entero_otu_mean_rowR <- Entero_otu_meanR[-c(3,7,11), ]
Entero_otu_mean_dist_TR <- distance(Entero_otu_mean_rowR, method = "bray")
vegan::mantel(Entero_otu_mean_dist_TR, cephER, method = "pearson", permutations = 1000000)

#enterobacteriales and gaster
merged_Entero <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Enterobacteriales")
merged_Entero_G <- subset_samples(merged_Entero, Body_Type == "G")
merged_Entero_meanG <- merge_samples(merged_Entero_G, "Species")
Entero_otu_meanG <- as(otu_table(merged_Entero_meanG), "matrix")
rowSums(Entero_otu_meanG)
Entero_otu_mean_rowG <- Entero_otu_meanG[-c(3,4,6,7,8,12), ]
Entero_otu_mean_dist_TG <- distance(Entero_otu_mean_rowG, method = "bray")
vegan::mantel(Entero_otu_mean_dist_TG, cephEG, method = "pearson", permutations = 1000000)

######Ricketssiales ###
#rickettsiales
merged_ric <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Rickettsiales")
merged_ric_mean<- merge_samples(merged_ric, "Species")
ric_otu_mean <- as(otu_table(merged_ric_mean), "matrix")
rowSums(ric_otu_mean)
ric_otu_mean_dist<- distance(ric_otu_mean, method = "bray")
vegan::mantel(ric_otu_mean_dist, cephtree, method = "pearson", permutations =1000000)

#rickettsiales and crop
merged_ric <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Rickettsiales")
merged_ric_C <- subset_samples(merged_ric, Body_Type == "C")
merged_ric_meanC <- merge_samples(merged_ric_C, "Species")
ric_otu_meanC <- as(otu_table(merged_ric_meanC), "matrix")
rowSums(ric_otu_meanC)
ric_otu_meanC_row <- ric_otu_meanC[-c(1,3), ]
ric_otu_mean_dist_TC <- distance(ric_otu_meanC_row, method = "bray")
vegan::mantel(ric_otu_mean_dist_TC, cephRC, method = "pearson", permutations =1000000)

#rickettsiales and midgut
merged_ric <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Rickettsiales")
merged_ric_M <- subset_samples(merged_ric, Body_Type == "M")
merged_ric_meanM<- merge_samples(merged_ric_M, "Species")
ric_otu_meanM <- as(otu_table(merged_ric_meanM), "matrix")
rowSums(ric_otu_meanM)
ric_otu_meanM_row <- ric_otu_meanM[-c(1,4,6,7,10), ]
ric_otu_mean_dist_TM <- distance(ric_otu_meanM_row, method = "bray")
vegan::mantel(ric_otu_mean_dist_TM, cephRM, method = "pearson", permutations =1000000)

#rickettsiales and ileum
merged_ric <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Rickettsiales")
merged_ric_I<- subset_samples(merged_ric, Body_Type == "I")
merged_ric_meanI<- merge_samples(merged_ric_I, "Species")
ric_otu_meanI <- as(otu_table(merged_ric_meanI), "matrix")
rowSums(ric_otu_meanI)
ric_otu_meanI_row <- ric_otu_meanI[-c(1,6,8,9,10,11), ]
ric_otu_mean_dist_TI <- distance(ric_otu_meanI_row, method = "bray")
vegan::mantel(ric_otu_mean_dist_TI, cephRI, method = "pearson", permutations =1000000)

#rickettsiales and rectum
merged_ric <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Rickettsiales")
merged_ric_R<- subset_samples(merged_ric, Body_Type == "R")
merged_ric_meanR<- merge_samples(merged_ric_R, "Species")
ric_otu_meanR <- as(otu_table(merged_ric_meanR), "matrix")
rowSums(ric_otu_meanR)
ric_otu_meanR_row <- ric_otu_meanR[-c(1,4,6), ]
ric_otu_mean_dist_TR <- distance(ric_otu_meanR_row, method = "bray")
vegan::mantel(ric_otu_mean_dist_TR, cephRR, method = "pearson", permutations =1000000)

#rickettsiales and gaster
merged_ric <- subset_taxa(physeq_rarefied_withpro_qiime_rooted, Order == "Rickettsiales")
merged_ric_G <- subset_samples(merged_ric, Body_Type == "G")
merged_ric_meanG<- merge_samples(merged_ric_G, "Species")
ric_otu_meanG <- as(otu_table(merged_ric_meanG), "matrix")
rowSums(ric_otu_meanG)
ric_otu_meanG_row <- ric_otu_meanG[-c(4,8,9,11), ]
ric_otu_mean_dist_TG <- distance(ric_otu_meanG_row, method = "bray")
vegan::mantel(ric_otu_mean_dist_TG, cephRG, method = "pearson", permutations =1000000)

#### Beta Diversity Statistics ####

#### permanova for GASTER###
physeq_G <- subset_samples(physeq_rarefied_withpro_qiime_rooted, Body_Type == "G")
ps.prop_allvend_ctrl <- transform_sample_counts(physeq_G, function(otu) {otu/sum(otu)})
set.seed(100)
ps_allvend_ctrl_pcoa.bray <- phyloseq::distance(ps.prop_allvend_ctrl, method = "bray")
sampledf <- data.frame(sample_data(physeq_G))
adonis(ps_allvend_ctrl_pcoa.bray ~ Species, data = sampledf)
#####permanova for crop###
physeq_C <- subset_samples(physeq_rarefied_withpro_qiime_rooted, Body_Type == "C")
ps.prop_allvend_ctrl <- transform_sample_counts(physeq_C, function(otu) {otu/sum(otu)})
set.seed(100)
ps_allvend_ctrl_pcoa.bray <- phyloseq::distance(ps.prop_allvend_ctrl, method = "bray")
sampledf <- data.frame(sample_data(physeq_C))
adonis(ps_allvend_ctrl_pcoa.bray ~ Species, data = sampledf)
####permanova for midgut ###
physeq_M <- subset_samples(physeq_rarefied_withpro_qiime_rooted, Body_Type == "M")
ps.prop_allvend_ctrl <- transform_sample_counts(physeq_M, function(otu) {otu/sum(otu)})
set.seed(100)
ps_allvend_ctrl_pcoa.bray <- phyloseq::distance(ps.prop_allvend_ctrl, method = "bray")
sampledf <- data.frame(sample_data(physeq_M))
adonis(ps_allvend_ctrl_pcoa.bray ~ Species, data = sampledf)
#####permanova for ileum###
physeq_I <- subset_samples(physeq_rarefied_withpro_qiime_rooted, Body_Type == "I")
ps.prop_allvend_ctrl <- transform_sample_counts(physeq_I, function(otu) {otu/sum(otu)})
set.seed(100)
ps_allvend_ctrl_pcoa.bray <- phyloseq::distance(ps.prop_allvend_ctrl, method = "bray")
sampledf <- data.frame(sample_data(physeq_I))
adonis(ps_allvend_ctrl_pcoa.bray ~ Species, data = sampledf)
sampledf <- data.frame(sample_data(physeq_G))
groups_GP <- sampledf[["Species"]]
bray_dispGP <-betadisper(ps_allvend_ctrl_pcoa.bray, groups_GP, type=c("median"))
anova(bray_dispGP)
permutest(bray_dispGP, pairwise= "TRUE")
####permanova for rectum###
physeq_R <- subset_samples(physeq_rarefied_withpro_qiime_rooted, Body_Type == "R")
ps.prop_allvend_ctrl <- transform_sample_counts(physeq_R, function(otu) {otu/sum(otu)})
set.seed(100)
ps_allvend_ctrl_pcoa.bray <- phyloseq::distance(ps.prop_allvend_ctrl, method = "bray")
sampledf <- data.frame(sample_data(physeq_R))
adonis2(ps_allvend_ctrl_pcoa.bray ~ Species, data = sampledf)
