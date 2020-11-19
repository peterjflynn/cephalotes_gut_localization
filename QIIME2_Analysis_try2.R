#QIIME2 Analysis 

##### 16S rRNA amplicon sequencing run ####
# Make paired end sequence artifact
qiime tools import \
--type EMPPairedEndSequences \
--input-path paired-end-sequences \
--output-path paired-end-sequences.qza

# Demultiplex the sequence reads
qiime demux emp-paired \
--m-barcodes-file sample-metadata_try3.tsv \
--m-barcodes-column BarcodeSequence \
--i-seqs emp-paired-end-sequences.qza \
--o-per-sample-sequences demux.qza \ 

qiime demux summarize \
--i-data demux.qza \
--o-visualization demux.qzv

# View demultiplexed reads plot

qiime tools view demux.qzv


#chose trim left at 12 because first time that base pair quality score went to 38
qiime dada2 denoise-paired \ 
--i-demultiplexed-seqs demux.qza \ 
--o-table table_run1.qza \ 
--o-representative-sequences rep-seqs-run1.qza \ 
--p-trim-left-f 12 \ 
--p-trim-left-r 12 \ 
--p-trunc-len-f 150 \ 
--p-trunc-len-r 150 


#### Run 2 (with extra Cephalotes texanus samples) ###
# Make paired end sequence artifact

qiime tools import \
--type EMPPairedEndSequences \
--input-path paired-end-sequences-run2 \
--output-path paired-end-sequences_run2.qza

qiime demux emp-paired-sequences-run2 \ 
--m-barcodes-file argonne-fixed-metadata-try3.tsv \ 
--m-barcodes-column BarcodeSequence \ 
--i-seqs emp-paired-end-sequences_run2.qza \ 
--o-per-sample-sequences demux_try2.qza  

qiime demux summarize \
--i-data demux_try2.qza \
--o-visualization demux_try2.qzv

qiime dada2 denoise-paired \ 
--i-demultiplexed-seqs demux_try2.qza \ 
--p-trim-left-f 12 \ 
--p-trim-left-r 12 \ 
--p-trunc-len-f 150 \ 
--p-trunc-len-r 150 \ 
--o-table table-run2.qza \ 
--o-representative-sequences rep-seqs-run2.qza \ 
--o-denoising-stats denoising-stats-run2.qza 

Filter  
qiime feature-table filter-samples \ 
--i-table table-run2.qza \ 
--m-metadata-file samples-to-keep.tsv \ 
--o-filtered-table id-filtered-run2.qza 

qiime feature-table filter-seqs \ 
--i-data rep-seqs-run2.qza \ 
--i-table id-filtered-run2.qza  \ 
--p-exclude-ids \ 
--o-filtered-data rep-seqs-filtered-run2.qza 

#visualize feature table from run2 
qiime feature-table summarize \
--i-table table_run1.qza \
--o-visualization table_run1.qzv \
--m-sample-metadata-file sample-metadata_try3.tsv 

#visualize table from run1
qiime feature-table summarize \
--i-table id-filtered-run2.qza \
--o-visualization id-filtered-run2.qzv \
--m-sample-metadata-file argonne-fixed-metadata-try3.tsv

qiime feature-table tabulate-seqs \
--i-data rep-seqs-filtered-run2.qza \
--o-visualization rep-seqs-filtered-run2.qzv

qiime feature-table merge \ 
--i-tables id-filtered-run2.qza \ 
--i-tables table_run1.qza \ 
--o-merged-table merged-table.qza 

qiime feature-table merge-seqs \ 
--i-data rep-seqs-run1.qza \ 
--i-data rep-seqs-filtered-run2.qza \ 
--o-merged-data merged-rep-seqs.qza 

FOR FILTERED AND MERGED DATA 
qiime tools import \ 
--type 'FeatureData[Sequence]' \ 
--input-path silva_132_99_16S.fna \ 
--output-path 99_silva_otus.qza 

qiime tools import \ 
--type 'FeatureData[Taxonomy]' \ 
--input-format HeaderlessTSVTaxonomyFormat \ 
--input-path taxonomy_7_levels.txt \ 
--output-path ref-taxonomy.qza 

#Extract Reference Reads
qiime feature-classifier extract-reads \ 
--i-sequences 99_silva_otus.qza \ 
--p-f-primer GTGCCAGCMGCCGCGGTAA \ 
--p-r-primer GGACTACHVGGGTWTCTAAT \ 
--p-trunc-len 120 \ 
--o-reads ref-seqs.qza 

#Train the Classifier
qiime feature-classifier fit-classifier-naive-bayes \ 
--i-reference-reads ref-seqs.qza \ 
--i-reference-taxonomy ref-taxonomy.qza \ 
--o-classifier classifier.qza 

#Test the classifier
qiime feature-classifier classify-sklearn \ 
--i-classifier classifier.qza \ 
--i-reads merged-rep-seqs.qza \ 
--o-classification taxonomy.qza 

#Extract Reference Reads
qiime feature-classifier extract-reads \ 
--i-sequences 99_silva_otus.qza \ 
--p-f-primer GTGCCAGCMGCCGCGGTAA \ 
--p-r-primer GGACTACHVGGGTWTCTAAT \ 
--p-trunc-len 120 \ 
--o-reads ref-seqs.qza 

qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ref-seqs.qza \
--i-reference-taxonomy ref-taxonomy.qza \
--o-classifier classifier.qza

qiime feature-classifier classify-sklearn \ 
--i-classifier classifier.qza \ 
--i-reads merged-rep-seqs.qza \ 
--o-classification taxonomy.qza 



qiime feature-table summarize \
--i-table merged-table.qza \
--o-visualization merged-table.qzv \
--m-sample-metadata-file sample-metadata_merged.tsv

qiime feature-table tabulate-seqs \
--i-data merged-rep-seqs.qza \
--o-visualization merged-rep-seqs.qzv

qiime phylogeny align-to-tree-mafft-fasttree \ 
--i-sequences merged-rep-seqs.qza \ 
--o-alignment aligned-rep-seqs.qza \ 
--o-masked-alignment masked-aligned-rep-seqs.qza \ 
--o-tree unrooted-tree.qza \ 
--o-rooted-tree rooted-tree.qza


### convert QIIME to phyloseq in R ###

# How to export a feature (OTU) table and convert from biom to .tsv 
# Step 1, export OTU table 
qiime tools export \
merged-table.qza \
--output-dir phyloseq

# OTU table exports as feature-table.biom so convert to .tsv
biom convert -i phyloseq/feature-table.biom -o phyloseq/otu_table.txt --to-tsv

# now you have otu_table.txt
# open it up in text edit and change #OTUID to OTUID

# Step 2, export taxonomy table 
qiime tools export \
taxonomy.qza\
--output-dir phyloseq

# now you have taxonomy.tsv
# open it up in text edit and change Feature ID to OTUID

# Step 3, export tree 
qiime tools export \
unrooted-tree.qza \
--output-dir phyloseq

# Step 4 QIIME2 doesn’t filter out taxonomy, so you have to merge the two files in R and output a merged file
(the rest of the code is in the R file)

#without procryptocerus
qiime tools import \ 
--input-path otu_decontam_nopro.biom \ 
--type 'FeatureTable[Frequency]' \ 
--input-format BIOMV100Format \ 
--output-path feature-table_phylo_nopro_decontam.qza 

#filter rep-seqs file for decontaminated OTUs 
qiime feature-table filter-seqs \ 
--i-data merged-rep-seqs.qza \ 
--i-table feature-table_phylo_nopro_decontam.qza   \ 
--p-no-exclude-ids \ 
--o-filtered-data rep-seqs-decontam-nopro.qza 

qiime phylogeny align-to-tree-mafft-fasttree \ 
--i-sequences rep-seqs-decontam-nopro.qza \ 
--o-alignment aligned-rep-seqs-decontam-nopro.qza \ 
--o-masked-alignment masked-aligned-rep-seqs-decontam-nopro.qza \ 
--o-tree unrooted-tree-decontam-nopro.qza \ 
--o-rooted-tree rooted-tree-decontam-nopro.qza 

#with procryptocerus
qiime tools import \ 
--input-path otu_decontam.biom \ 
--type 'FeatureTable[Frequency]' \ 
--input-format BIOMV100Format \ 
--output-path feature-table_phylo_decontam.qza 

qiime feature-table filter-seqs \ 
--i-data merged-rep-seqs.qza \ 
--i-table feature-table_phylo_decontam.qza   \ 
--p-no-exclude-ids \ 
--o-filtered-data rep-seqs-decontam-withpro.qza 

qiime phylogeny align-to-tree-mafft-fasttree \ 
--i-sequences rep-seqs-decontam-withpro.qza \ 
--o-alignment aligned-rep-seqs-decontam-withpro.qza \ 
--o-masked-alignment masked-aligned-rep-seqs-decontam-withpro.qza \ 
--o-tree unrooted-tree-decontam-withpro.qza \ 
--o-rooted-tree rooted-tree-decontam-withpro.qza 

####

qiime diversity core-metrics-phylogenetic \
--i-phylogeny rooted-tree-decontam-nopro.qza \ 
--i-table feature-table_phylo_nopro_decontam.qza \ 
--p-sampling-depth 13814 \ 
--m-metadata-file sample-metadata_merged_qpcr1_nopro_again.tsv \ 
--output-dir core-metrics-results-rarefied-nopro-rooted-try2

qiime diversity core-metrics-phylogenetic \
--i-phylogeny rooted-tree-decontam-withpro.qza \ 
--i-table feature-table_phylo_decontam.qza \ 
--p-sampling-depth 13335 \ 
--m-metadata-file sample-metadata_rarefied_withpro_qiime_caste.tsv \ 
--output-dir core-metrics-results-rarefied-withpro-rooted-try2 


#alpha diversity statistics

qiime diversity alpha-group-significance \ 
--i-alpha-diversity shannon_vector-nopro.qza \ 
--m-metadata-file sample-metadata_merged_qpcr1_nopro_again.tsv \ 
--o-visualization shannon-group-significance-nopro.qzv 

qiime diversity alpha-group-significance \ 
--i-alpha-diversity evenness_vector-nopro.qza \ 
--m-metadata-file sample-metadata_merged_qpcr1_nopro_again.tsv \ 
--o-visualization evenness-group-significance-nopro.qzv 

qiime diversity alpha-group-significance \ 
--i-alpha-diversity faith_pd_vector-nopro.qza \ 
--m-metadata-file sample-metadata_merged_qpcr1_nopro_again.tsv \ 
--o-visualization faith-pd-group-significance-nopro.qzv 

qiime diversity alpha-group-significance \ 
--i-alpha-diversity observed_otus_vector-nopro.qza \ 
--m-metadata-file sample-metadata_merged_qpcr1_nopro_again.tsv \ 
--o-visualization observed-otus-group-significance-nopro.qzv 

#beta diversity statistics 
qiime diversity adonis \
--i-distance-matrix weighted_unifrac_distance_matrix.qza \
--m-metadata-file sample-metadata_merged_qpcr1_nopro_again.tsv \
--p-formula "Body_Type*Species*Caste_Type*Colony" \
--p-permutations 999 \
--o-visualization ADONIS-weighted_unifrac_colony.qzv


qiime diversity adonis \
--i-distance-matrix unweighted_unifrac_distance_matrix.qza \
--m-metadata-file sample-metadata_merged_qpcr1_nopro_again.tsv \
--p-formula "Body_Type*Species*Caste_Type*Colony" \
--p-permutations 999 \
--o-visualization ADONIS-unweighted_unifrac_colony.qzv


qiime diversity adonis \
--i-distance-matrix bray_curtis_distance_matrix.qza \
--m-metadata-file sample-metadata_merged_qpcr1_nopro_again.tsv \
--p-formula "Body_Type*Species*Caste_Type*Colony" \
--p-permutations 999 \
--o-visualization ADONIS-bray_colony.qzv

qiime diversity adonis \
--i-distance-matrix bray_curtis_distance_matrix.qza \
--m-metadata-file sample-metadata_merged_qpcr1_nopro_again.tsv \
--p-formula "Body_Type*Species*Caste_Type" \
--p-permutations 999 \
--o-visualization ADONIS-bray_nocolony.qzv

qiime diversity adonis \
--i-distance-matrix bray_curtis_distance_matrix.qza \
--m-metadata-file sample-metadata_merged_qpcr1_nopro_again.tsv \
--p-formula “Species" \
--p-permutations 999 \
--o-visualization ADONIS-bray_species.qzv

qiime diversity adonis \
--i-distance-matrix bray_curtis_distance_matrix.qza \
--m-metadata-file sample-metadata_merged_qpcr1_nopro_again.tsv \
--p-formula “Species" \
--p-permutations 999 \
--o-visualization ADONIS-bray_species.qzv

#permdisp

qiime diversity beta-group-significance \
--i-distance-matrix bray_curtis_distance_matrix.qza \
--m-metadata-file sample-metadata_merged_qpcr1_nopro_again.tsv \
--m-metadata-column Species \
--p-method 'permdisp' \
--p-pairwise \
--o-visualization bray__species_permdisp.qzv

qiime diversity beta-group-significance \
--i-distance-matrix core-metrics-results-rarefied-nopro-rooted-try2/weighted_unifrac_distance_matrix.qza \
--m-metadata-file sample-metadata_merged_qpcr1_nopro_again.tsv \
--m-metadata-column Body_Type \
--p-method 'permdisp' \
--p-pairwise \
--o-visualization core-metrics-results-rarefied-nopro-rooted-try2/wunifrac_body_permdisp.qzv


#beta diversity statistics: pairwise

qiime diversity beta-group-significance \ 
--i-distance-matrix unweighted_unifrac_distance_matrix.qza \ 
--m-metadata-file sample-metadata_merged_qpcr1_nopro_again.tsv \ 
--m-metadata-column Body_Type \
--o-visualization unweighted-unifrac-compartment-significance.qzv \ 
--p-pairwise
qiime diversity beta-group-significance \ 
--i-distance-matrix weighted_unifrac_distance_matrix.qza \ 
--m-metadata-file sample-metadata_merged_qpcr1_nopro_again.tsv \ 
--m-metadata-column Body_Type \
--o-visualization weighted-unifrac-compartment-significance.qzv \ 
--p-pairwise

qiime diversity beta-group-significance \ 
--i-distance-matrix weighted_unifrac_distance_matrix.qza \ 
--m-metadata-file sample-metadata_merged_qpcr1_nopro_again.tsv \ 
--m-metadata-column Species \
--o-visualization weighted-unifrac-species-significance.qzv \ 
--p-pairwise

qiime diversity beta-group-significance \ 
--i-distance-matrix unweighted_unifrac_distance_matrix.qza \ 
--m-metadata-file sample-metadata_merged_qpcr1_nopro_again.tsv \ 
--m-metadata-column Species \
--o-visualization unweighted-unifrac-species-significance.qzv \ 
--p-pairwise

qiime diversity beta-group-significance \ 
--i-distance-matrix weighted_unifrac_distance_matrix.qza \ 
--m-metadata-file sample-metadata_merged_qpcr1_nopro_again.tsv \ 
--m-metadata-column Caste_Type \
--o-visualization weighted-unifrac-caste-significance.qzv \ 
--p-pairwise

qiime diversity beta-group-significance \ 
--i-distance-matrix unweighted_unifrac_distance_matrix.qza \ 
--m-metadata-file sample-metadata_merged_qpcr1_nopro_again.tsv \ 
--m-metadata-column Caste_Type \
--o-visualization unweighted-unifrac-caste-significance.qzv \ 
--p-pairwise


qiime diversity beta-group-significance \ 
--i-distance-matrix bray_curtis_distance_matrix.qza \
--m-metadata-file sample-metadata_merged_qpcr1_nopro_again.tsv \ 
--m-metadata-column Body_Type \
--o-visualization bray-compartment-significance.qzv \ 
--p-pairwise

qiime diversity beta-group-significance \ 
--i-distance-matrix bray_curtis_distance_matrix.qza \
--m-metadata-file sample-metadata_merged_qpcr1_nopro_again.tsv \ 
--m-metadata-column Species \
--o-visualization bray-species-significance.qzv \ 
--p-pairwise

qiime diversity beta-group-significance \ 
--i-distance-matrix bray_curtis_distance_matrix.qza \
--m-metadata-file sample-metadata_merged_qpcr1_nopro_again.tsv \ 
--m-metadata-column Caste_Type \
--o-visualization bray-caste-significance.qzv \ 
--p-pairwise

qiime diversity beta-group-significance \ 
--i-distance-matrix bray_curtis_distance_matrix.qza \
--m-metadata-file sample-metadata_merged_qpcr1_nopro_again.tsv \ 
--m-metadata-column Colony \
--o-visualization bray-colony-significance.qzv \ 
--p-pairwise

qiime diversity beta-group-significance \ 
--i-distance-matrix weighted_unifrac_distance_matrix.qza \
--m-metadata-file sample-metadata_rarefied_nopro_qiime_colony.tsv \ 
--m-metadata-column Colony \
--o-visualization weighted_unifrac-colony-significance.qzv \ 
--p-pairwise

qiime diversity beta-group-significance \ 
--i-distance-matrix unweighted_unifrac_distance_matrix.qza \
--m-metadata-file sample-metadata_rarefied_nopro_qiime_colony.tsv \ 
--m-metadata-column Colony \
--o-visualization unweighted_unifrac-colony-significance.qzv \ 
--p-pairwise


#no gaster

qiime tools import \ 
--type 'FeatureTable[Frequency]'  \ 
--input-path otu_decontam_nopro_nogaster.biom \ 
--input-format BIOMV100Format \
--output-path feature-table_phylo_nopro_nogaster_decontam.qza


qiime feature-table filter-seqs \
--i-data merged-rep-seqs.qza  \
--i-table feature-table_phylo_nopro_nogaster_decontam.qza \
--p-no-exclude-ids  \
--o-filtered-data rep-seqs-decontam-nopro-nogaster.qza


qiime phylogeny align-to-tree-mafft-fasttree \ 
--i-sequences rep-seqs-decontam-nopro_nogaster.qza \ 
--o-alignment aligned-rep-seqs-decontam-nopro_noG.qza \ 
--o-masked-alignment masked-aligned-rep-seqs-decontam-nopro_noG.qza \ 
--o-tree unrooted-tree-decontam-nopro_noG.qza \ 
--o-rooted-tree  rooted-tree-decontam-nopro_noG.qza  

qiime diversity core-metrics-phylogenetic \ 
--i-phylogeny rooted-tree-decontam-nopro_noG.qza \ 
--i-table feature-table_phylo_nopro_nogaster_decontam.qza \ 
--p-sampling-depth 13814 \ 
--m-metadata-file sample-metadata_merged_qpcr1_nopro_again.tsv \ 
--output-dir core-metrics-results-rarefied-nopro-rooted-nogaster


qiime diversity adonis \
--i-distance-matrix core-metrics-results-rarefied-nopro-rooted-nogaster/weighted_unifrac_distance_matrix_nogaster.qza  \
--m-metadata-file sample-metadata_merged_qpcr1_nopro_again.tsv \
--p-formula "Body_Type*Species*Caste_Type*Colony" \
--p-permutations 999 \
--o-visualization core-metrics-results-rarefied-nopro-rooted-nogaster/ADONIS-weighted_unifrac_nogaster.qzv 



#### Bacterial Phylogenies for JANE Analysis #####

#importing the feature table from phyloseq
qiime tools import \ 
--input-path physeq_rarefied_withpro_qiime_rooted.biom \ 
--type 'FeatureTable[Frequency]' \ 
--input-format BIOMV100Format \ 
--output-path feature-table_physeq_rarefied_withpro_qiime_rooted.qza 

#For grouping by species 
qiime feature-table group \ 
--i-table feature-table_physeq_rarefied_withpro_qiime_rooted.qza\ 
--p-axis sample \ 
--m-metadata-file sample-metadata_rarefied_withpro_qiime.tsv\ 
--m-metadata-column Species \ 
--p-mode sum \ 
--output-dir species-rarefied-withpro-rooted-merged-table.qza

qiime feature-table summarize \ 
--i-table species-rarefied-withpro-rooted-merged-table.qza \ 
--o-visualization species-rarefied-withpro-rooted-merged-table.qzv 

#filter rep-seqs file for decontaminated OTUs 
qiime feature-table filter-seqs \ 
--i-data merged-rep-seqs.qza \ 
--i-table species-rarefied-withpro-rooted-merged-table.qza \ 
--p-no-exclude-ids \ 
--o-filtered-data rep-seqs-species-rarefied-withpro-rooted-merged-table.qza 


#betaproteobacteriales #
Betaproteo Filtering Sequences and new subset tree 
qiime taxa filter-table \
--i-table species-rarefied-withpro-rooted-merged-table.qza \ 
--i-taxonomy taxonomy.qza \
--p-include  'D_0__Bacteria;D_1__Proteobacteria;D_2__Gammaproteobacteria;D_3__Betaproteobacteriales' \
--o-filtered-table betaproteo-table.qza

qiime feature-table summarize \ 
--i-table betaproteo-table.qza \ 
--o-visualization betaproteo-table.qzv 

qiime feature-table filter-seqs \
--i-data rep-seqs-species-rarefied-withpro-rooted-merged-table.qza \
--m-metadata-file betaproteo-ASVs-to-keep.tsv \ 
--p-no-exclude-ids \
--o-filtered-data betaproteo-merged-rep-seqs.qza 

qiime feature-table tabulate-seqs \ 
--i-data betaproteo-merged-rep-seqs.qza \ 
--o-visualization betaproteo-merged-rep-seqs.qzv 

qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences betaproteo-merged-rep-seqs.qza\ 
--o-alignment aligned-betaproteo-rep-seqs.qza \
--o-masked-alignment masked-betaproteo-aligned-rep-seqs.qza \
--o-tree unrooted-betaproteo-tree.qza \
--o-rooted-tree rooted-betaproteo-tree.qza

qiime tools export \
--input-path unrooted-betaproteo-tree.qza\
--output-path betaproteo_species 

qiime tools export \
--input-path rooted-betaproteo-tree.qza\
--output-path betaproteo_species 

qiime tools extract \ 
--input-path betaproteo-table.qza \ 
--output-path extracted-feature-table-betaproteo 

biom convert -i feature-table.biom -o betaproteo_otu_table.txt --to-tsv 


#RHIZOBIALES#

#Rhizobiales Filtering Sequences and new subset tree #
qiime taxa filter-table \
--i-table species-rarefied-withpro-rooted-merged-table.qza \ 
--i-taxonomy taxonomy.qza \
--p-include 'D_0__Bacteria;D_1__Proteobacteria;D_2__Alphaproteobacteria;D_3__Rhizobiales' \
--o-filtered-table rhizobiales-table1.qza

qiime feature-table summarize \ 
--i-table rhizobiales-table1.qza \
--o-visualization rhizobiales-table1.qzv

qiime feature-table filter-seqs \
--i-data rep-seqs-species-rarefied-withpro-rooted-merged-table.qza \ 
--m-metadata-file rhizobiales-ASVs-to-keep.tsv \ 
--p-no-exclude-ids \
--o-filtered-data rhizobiales-merged-rep-seqs1.qza

qiime feature-table tabulate-seqs \ 
--i-data rhizobiales-merged-rep-seqs1.qza \ 
--o-visualization rhizobiales-merged-rep-seqs1.qzv 

qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences rhizobiales-merged-rep-seqs1.qza \ 
--o-alignment aligned-rhizobiales-rep-seqs1.qza \
--o-masked-alignment masked-rhizobiales-aligned-rep-seqs1.qza \
--o-tree unrooted-rhizobiales-tree1.qza \
--o-rooted-tree rooted-rhizobiales-tree1.qza

qiime tools export \
--input-path unrooted-rhizobiales-tree1.qza\
--output-path exported-feature-table

qiime tools export \
--input-path rooted-rhizobiales-tree1.qza\ 
--output-path rooted_rhizo

qiime tools export \
--input-path unrooted-rhizobiales-tree1.qza\
--output-path Rhizobiales_species

mkdir extracted-feature-table 
qiime tools extract \ 
--input-path Rhizobiales_species/rhizobiales-table1.qza \ 
--output-path Rhizobiales_species/extracted-feature-table-rhizo 

qiime tools export \ 
--input-path Rhizobiales_species/rhizobiales-table1.qza \ 
--output-path exported-feature-table 

biom convert -i feature-table.biom -o otu_table.txt --to-tsv

#OPITUTALES#
#Opitutales Filtering Sequences and new subset tree #
qiime taxa filter-table \
--i-table species-rarefied-withpro-rooted-merged-table.qza \ 
--i-taxonomy taxonomy.qza \
--p-include  'D_0__Bacteria;D_1__Verrucomicrobia;D_2__Verrucomicrobiae;D_3__Opitutales' \
--o-filtered-table opitutales-table-rare.qza

qiime feature-table summarize \ 
--i-table opitutales-table-rare.qza \ 
--o-visualization opitutales-table-rare.qzv 

qiime feature-table filter-seqs \
--i-data rep-seqs-species-rarefied-withpro-rooted-merged-table.qza \ 
--m-metadata-file opitutales-ASVs-to-keep-rare.tsv \ 
--p-no-exclude-ids \
--o-filtered-data opitutales-rare-merged-rep-seqs.qza

qiime feature-table tabulate-seqs \ 
--i-data opitutales-rare-merged-rep-seqs.qza \ 
--o-visualization opitutales-rare-merged-rep-seqs.qzv 

qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences opitutales-rare-merged-rep-seqs.qza\ 
--o-alignment aligned-opitutales-rare-rep-seqs.qza \
--o-masked-alignment masked-opitutales-rare-aligned-rep-seqs.qza \
--o-tree unrooted-opitutales-rare-tree.qza \
--o-rooted-tree rooted-opitutales-rare-tree.qza

qiime tools export \
--input-path unrooted-opitutales-rare-tree.qza\
--output-path Opitutales1

qiime tools export \
--input-path rooted-opitutales-rare-tree.qza\
--output-path Opitutales1

qiime tools extract \ 
--input-path opitutales-table-rare.qza \ 
--output-path extracted-feature-table-opit 

biom convert -i feature-table.biom -o opit_otu_table.txt --to-tsv

#RICKETTSIALES#
Rickettsiales Filtering Sequences and new subset tree 
qiime taxa filter-table \
--i-table species-rarefied-withpro-rooted-merged-table.qza \ 
--i-taxonomy taxonomy.qza \
--p-include  'D_0__Bacteria;D_1__Proteobacteria;D_2__Alphaproteobacteria;D_3__Rickettsiales' \
--o-filtered-table rickettsiales-table.qza

qiime feature-table summarize \ 
--i-table rickettsiales-table.qza \ 
--o-visualization rickettsiales-table.qzv 

qiime feature-table filter-seqs \
--i-data rep-seqs-species-rarefied-withpro-rooted-merged-table.qza \
--m-metadata-file rickettsiales-ASVs-to-keep.tsv \ 
--p-no-exclude-ids \
--o-filtered-data rickettsiales-merged-rep-seqs.qza 

qiime feature-table tabulate-seqs \ 
--i-data rickettsiales-merged-rep-seqs.qza \ 
--o-visualization rickettsiales-merged-rep-seqs.qzv 

qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences rickettsiales-merged-rep-seqs.qza\ 
--o-alignment aligned-rickettsiales-rep-seqs.qza \
--o-masked-alignment masked-rickettsiales-aligned-rep-seqs.qza \
--o-tree unrooted-rickettsiales-tree.qza \
--o-rooted-tree rooted-rickettsiales-tree.qza

qiime tools export \
--input-path unrooted-rickettsiales-tree.qza\
--output-path Rickettsiales_species

qiime tools export \
--input-path rooted-rickettsiales-tree.qza\
--output-path Rickettsiales_species

qiime tools extract \ 
--input-path rickettsiales-table.qza \ 
--output-path extracted-feature-table-rickettsiales 

biom convert -i feature-table.biom -o rickettsiales_otu_table.txt --to-tsv

#Enterobacteriales #
Enterobacteriales Filtering Sequences and new subset tree 
qiime taxa filter-table \
--i-table species-rarefied-withpro-rooted-merged-table.qza \ 
--i-taxonomy taxonomy.qza \
--p-include  'D_0__Bacteria;D_1__Proteobacteria;D_2__Gammaproteobacteria;D_3__Enterobacteriales' \
--o-filtered-table enterobacteriales-table.qza

qiime feature-table summarize \ 
--i-table enterobacteriales-table.qza \ 
--o-visualization enterobacteriales-table.qzv 

qiime feature-table filter-seqs \
--i-data rep-seqs-species-rarefied-withpro-rooted-merged-table.qza \ 
--m-metadata-file enterobacteriales-ASVs-to-keep.tsv \ 
--p-no-exclude-ids \
--o-filtered-data enterobacteriales-merged-rep-seqs.qza 

qiime feature-table tabulate-seqs \ 
--i-data enterobacteriales-merged-rep-seqs.qza \ 
--o-visualization enterobacteriales-merged-rep-seqs.qzv 

qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences enterobacteriales-merged-rep-seqs.qza\ 
--o-alignment aligned-enterobacteriales-rep-seqs.qza \
--o-masked-alignment masked-enterobacteriales-aligned-rep-seqs.qza \
--o-tree unrooted-enterobacteriales-tree.qza \
--o-rooted-tree rooted-enterobacteriales-tree.qza

qiime tools export \
--input-path unrooted-enterobacteriales-tree.qza\
--output-path Enterobacteriales_species

qiime tools export \
--input-path rooted-enterobacteriales-tree.qza\
--output-path Enterobacteriales_species

qiime tools extract \ 
--input-path Enterobacteriales-table.qza \ 
--output-path extracted-feature-table-entero 

biom convert -i feature-table.biom -o entero_otu_table.txt --to-tsv

#Xanthomonadales #
#Xanthomonadales Filtering Sequences and new subset tree #
qiime taxa filter-table \
--i-table species-rarefied-withpro-rooted-merged-table.qza \ 
--i-taxonomy taxonomy.qza \
--p-include  'D_0__Bacteria;D_1__Proteobacteria;D_2__Gammaproteobacteria;D_3__Xanthomonadales' \
--o-filtered-table Xanthomonadales-table.qza 

qiime feature-table summarize \ 
--i-table Xanthomonadales-table.qza \ 
--o-visualization Xanthomonadales-table.qzv 

qiime feature-table filter-seqs \
--i-data rep-seqs-species-rarefied-withpro-rooted-merged-table.qza \ 
--m-metadata-file Xanthomonadales-ASVs-to-keep.tsv \ 
--p-no-exclude-ids \
--o-filtered-data Xanthomonadales-merged-rep-seqs.qza 

qiime feature-table tabulate-seqs \ 
--i-data Xanthomonadales-merged-rep-seqs.qza \ 
--o-visualization Xanthomonadales-merged-rep-seqs.qzv 

qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences Xanthomonadales-merged-rep-seqs.qza\ 
--o-alignment aligned-Xanthomonadales-rep-seqs.qza \
--o-masked-alignment masked-Xanthomonadales-aligned-rep-seqs.qza \
--o-tree unrooted-Xanthomonadales-tree.qza \
--o-rooted-tree rooted-Xanthomonadales-tree.qza

qiime tools export \
--input-path unrooted-Xanthomonadales-tree.qza\
--output-path Xanthomonadales_species 

qiime tools export \
--input-path rooted-Xanthomonadales-tree.qza\
--output-path Xanthomonadales_species 

qiime tools extract \ 
--input-path Xanthomonadales-table.qza \ 
--output-path extracted-feature-table-xantha 

biom convert -i feature-table.biom -o otu_table.txt --to-tsv


