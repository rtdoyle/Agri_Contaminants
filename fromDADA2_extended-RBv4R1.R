#
#title: "from DADA2"
#author: "Laura Rossi"
#date: "Sept26, 2018"
#
setwd("C:/Users/Laura/Documents/labwork/2023/RB/")
path ="C:/Users/Laura/Documents/labwork/2023/RB//"
getwd()
# load libraries for ordination and clustering
library("cluster")
library(data.table)
library("phyloseq")
packageDescription("phyloseq")$Version
library("ggplot2")
library("plyr")
library("grid")
library("ape")
library("phangorn")
library("phytools")
library("vegan")
library("ShortRead")
library("seqinr")
#install.packages("seqRFLP")
library("seqRFLP")
#install.packages("pipeR")
#install.packages('devtools')
devtools::load_all("C:/Users/Laura/Documents/labwork/R_scripts/aftersl1p-master/aftersl1p-master")
theme_set(theme_bw()) ## Default is gray, In case you want the white
#
#import files to make the phyloseq object 
otufile = "seqtab_nochimtransposed_R1v4_RBtest_l200.csv"
otu_df = read.csv(otufile, row.names = 1)   # if there are non-otu table lines at the top, use the skip = argument.
#rownames(otu_df)=otu_df[,1]  #assign contents of first column to rownames
#otu_df=otu_df[,-1]  #removed first column
seqs = rownames(otu_df)
rownames(otu_df) = NULL
#
taxfile = "taxa_R1v4_RBtest_silva138wsp.csv"
tax_df = read.csv(taxfile, row.names = 1)

#rownames(tax_df)=tax_df$X.TAXONOMY #assign contents of first column to rownames
#tax_df=tax_df[,-1] #remove first column
all(seqs == rownames(tax_df))
rownames(tax_df) = NULL
#
metfile = "mapfilev4R1.csv"
met_df = read.csv(metfile) # read tab-delimited instead of comma
#the following commands help to clean up the map file, they are not necessary but may
#help avoid errors in creating the phyloseq object
#for the following commands to work the column label of your sample names, cellA1
#needs to be #NAME - or - you need to change X.NAME (after the $) to whatever you have 
#labelled the first column
met_df$NAME = gsub('-', '.', as.character(met_df$NAME))#
met_df$NAME = gsub(' ', '.', met_df$NAME)
met_df$NAME = as.factor(met_df$NAME)
rownames(met_df) = met_df$NAME   # make sure the rownames in the metadata data frame are the same as the row (or column, whichever is samples) in the otu data frame
#create the phyloseq object
dat = phyloseq(otu_table(otu_df, taxa_are_rows = TRUE), # or FALSE if false
               tax_table(as.matrix(tax_df)),
               sample_data(met_df))

#preprocessing
#write ASV sequences to file
write.csv(seqs, file="ASV_sequencesv4R1.csv")
#try writing fasta file with new package **this works!
seqs.fasta = dataframe2fas(seqs, file="ASVseqs.fasta")
#simple sums and write to file
samplesums = sort(sample_sums(dat))
write.csv(samplesums, file="readcounts_per_sample_sorted.csv")
plot(samplesums)
head(samplesums)
#filter table of non-bacterial taxa
dat_lessHOST = subset_taxa(dat, Kingdom=="Bacteria", Family!="Mitochondria")
#write.csv(otu_table(dat_lessHOST),file='seqtab_nochim_lessHOST_AKLRCC.csv')
#write.csv(tax_table(dat_lessHOST),file='taxa_AKLRCC_lessHOST.csv')
#remove less than 0.001% ?
rel_abun_all = transform_sample_counts(dat_lessHOST, function(x) x/sum(x))
rel_abun_all_prune = prune_taxa(taxa_sums(rel_abun_all) > 0.001,
                                rel_abun_all)
#remove less than 5 reads
dat_lessHOST_n5 = prune_taxa(taxa_sums(dat_lessHOST)>5, dat_lessHOST)
#rarefy
min_lib <- min(sample_sums(dat))
dat_r =rarefy_even_depth(dat_lessHOST,sample.size=50000,verbose=FALSE,replace = FALSE)
#dat_r1000=rarefy_even_depth(dat_lessHOST, sample.size=1000, verbose=FALSE, replace=FALSE)
#write.csv(otu_table(dat_r),file="rarefy_otu_table.csv")
#glom to Genus level (change "Genus" if you want the other levels)
dat_lessHOST_genus <- tax_glom(dat_lessHOST, taxrank="Genus")
write.csv(otu_table(dat_lessHOST_genus),file="otu_table_noHOST_genus.csv")
#set colour palettes
colours <- c("#330000", "#594c16", "#6cd9b5", "#0000ff", "#b32d86", "#d91d00", "#f2e6b6", "#006652", "#000066", "#d9a3bf", "#cc7466", "#eeff00", "#394d4b", "#333366", "#7f0033", "#59392d", "#6b7300", "#6cd2d9", "#8979f2", "#e5003d", "#ffd0bf", "#caf279", "#1d6273", "#732699", "#f27999", "#f2aa79", "#5a664d", "#3db6f2", "#524359", "#7f0011", "#a65800", "#004d00", "#001b33", "#cc00ff", "#592d33", "#a68a53", "#00ff22", "#296ca6", "#2b0d33", "#e5b800", "#2db33e", "#a3aad9", "#f780ff", "#330000", "#594c16", "#6cd9b5", "#0000ff", "#b32d86", "#d91d00", "#f2e6b6")
colours2 <- c("#4363d8", "#f58231", "#911eb4",
              "#e6194B", "#3cb44b", "#ffe119",
              "#42d4f4", "#f032e6", "#bfef45",
              "#fabebe", "#469990", "#e6beff",
              "#9A6324", "#fffac8", "#800000",
              "#aaffc3", "#808000", "#ffd8b1",
              "#000075", "#a9a9a9", "#ffffff",
              "#000000", "#330000", "#594c16",
              "#6cd9b5", "#0000ff", "#b32d86",
              "#d91d00", "#f2e6b6", "#006652", "#000066","#4363d8", "#f58231", "#911eb4",
              "#e6194B", "#3cb44b", "#ffe119",
              "#42d4f4", "#f032e6", "#bfef45",
              "#fabebe", "#469990", "#e6beff",
              "#9A6324", "#fffac8", "#800000",
              "#aaffc3", "#808000", "#ffd8b1",
              "#000075", "#a9a9a9", "#ffffff",
              "#000000", "#330000", "#594c16",
              "#6cd9b5", "#0000ff", "#b32d86",
              "#d91d00", "#f2e6b6", "#006652", "#000066")
#
#alpha diversity
plot_richness(dat_lessHOST, x="SampleID", measures=c("Shannon", "Simpson"))
plot_richness(dat_r, x="SampleID", measures=c("Shannon", "Simpson"))
#alpha = estimate_richness(dat_lessHOST, measures=c("Shannon", "Simpson"))
#write.csv(alpha, file="alpha_Shannon_Simpson_lessHOST.csv")
alpha1 = plot_richness(dat_r, x="SampleID", measures=c("Shannon"), color = "Region")
alpha1 + scale_color_manual(values = colours2) + geom_point(size=3) + scale_y_continuous(limits = c(1, 5))
alpha2 = plot_richness(dat_r, x="SampleID", measures=c("Simpson"), color = "Region")
alpha2 + scale_color_manual(values = colours2) + geom_point(size=3) + scale_y_continuous(limits = c(0, 1))
#
#PCoA-relative abundance OTUS
relabunbray.ord <- ordinate(rel_abun_all_prune, method = "PCoA", distance = "bray")
relabunbray.plot <- plot_ordination(rel_abun_all_prune, relabunbray.ord,
                                    color = "SampleID",
                                    shape = "PCR",
                                    axes = c(1,2),
                                    title = "Bray-Curtis rel abund")
relabunbray.plot+geom_point(size=5)+scale_color_manual(values=colours2)
rarefy.ord = ordinate(dat_r, method="PCoA", distance = "bray")
prare.plot <- plot_ordination(dat_r, rarefy.ord,
                           color = "SampleID",
                           shape = "PCR",
                           axes = c(1,2),
                           title = "Bray-Curtis rarefy")
prare.plot+geom_point(size=5)+scale_color_manual(values=colours2)
#
#barcharts
top20 <- names(sort(taxa_sums(dat_lessHOST), decreasing=TRUE))[1:20]
dat.top20 <- transform_sample_counts(dat, function(OTU) OTU/sum(OTU))
dat.top20 <- prune_taxa(top20, dat.top20)
top50 <- names(sort(taxa_sums(dat), decreasing=TRUE))[1:50]
dat.top50 <- transform_sample_counts(dat, function(OTU) OTU/sum(OTU))
dat.top50 <- prune_taxa(top50, dat.top50)
top100 <- names(sort(taxa_sums(dat), decreasing=TRUE))[1:100]
dat.top100 <- transform_sample_counts(dat, function(OTU) OTU/sum(OTU))
dat.top100 <- prune_taxa(top100, dat.top100)
plot_bar(dat.top50, fill='Phylum')#too busy
#try new colours
tax_plot=plot_bar(dat.top50, fill='Genus')
tax_plot + geom_bar(aes(), stat="identity", position="stack") + scale_fill_manual(values = colours) + theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))
#try some things
tax_plot3 = plot_bar(dat.top50, fill="Genus", title="top 50 Genus")
tax_plot3 + geom_bar(aes(), stat="identity", position="fill")+
  facet_wrap("SampleID", scales="free_x")+
  scale_fill_manual(values = colours2)+
  theme(legend.position="right", legend.key.size = unit(0.3, "cm"))+
  guides(fill=guide_legend(nrow=35))
tax_plot7 = plot_bar(dat.top100, fill="Genus", title="top 100 Genus")
tax_plot7 + geom_bar(aes(), stat="identity", position="fill")+
  facet_wrap("SampleID", scales="free_x")+
  scale_fill_manual(values = colours2)+
  theme(legend.position="right", legend.key.size = unit(0.5, "cm"))+
  guides(fill=guide_legend(nrow=25))
tax_plot8 = plot_bar(dat.top100, fill="Family", title="top 100 Family")
tax_plot8 + geom_bar(aes(), stat="identity", position="fill")+
  facet_wrap("SampleID", scales="free_x")+
  scale_fill_manual(values = colours2)+
  theme(legend.position="right", legend.key.size = unit(0.5, "cm"))+
  guides(fill=guide_legend(nrow=25))
tax_plot5 = plot_bar(dat.top50, x="Media", fill="Family", title="top 50 Family")
tax_plot5 + geom_bar(aes(), stat="identity", position="fill")+
  facet_grid(rows = vars(SampleNum), cols = vars(Condition), scales="free_x")+
  scale_fill_manual(values = colours2) +
  theme(legend.position="bottom", legend.key.size = unit(0.5, "cm"))+
  guides(fill=guide_legend(nrow=10))
tax_plot8 = plot_bar(dat.top100, x="Media", fill="Genus", title="top 100 Genus")
tax_plot8 + geom_bar(aes(), stat="identity", position="fill")+ facet_grid(rows = vars(SampleNum), cols = vars(Condition), scales="free_x") + scale_fill_manual(values = colours2) + theme(legend.position="right", legend.key.size = unit(0.2, "cm"))+ guides(fill=guide_legend(nrow=35))
tax_plot6 = plot_bar(dat.top100, x="Media", fill="Genus", title="top 100 Genus")
tax_plot6 + geom_bar(aes(), stat="identity", position="fill")+ facet_wrap("Condition", scales="free_x")+ scale_fill_manual(values = colours) + theme(legend.position="right", legend.key.size = unit(0.5, "cm"))+ guides(fill=guide_legend(nrow=25))

#Jake's aftersl1p
phy_df = make_phy_df(rel_abun_all_prune, rank="Family", cutoff=0.01)
genus_plot = plot_tax_bar(phy_df, rank ='Family', sample="Media", )
genus_plot + facet_wrap(~Condition, scales='free_x')+ theme(legend.position="bottom", legend.key.size = unit(0.2, "cm"))+ guides(fill=guide_legend(nrow=8))
genus_plot + facet_grid(rows = vars(SampleNum), cols = vars(Condition), scales="free_x")+ theme(legend.position="bottom", legend.key.size = unit(0.2, "cm"))+ guides(fill=guide_legend(nrow=6))
#for PCoA's 
make_ord_df()
plt_ord()
samreaddepth = plot_read_depth(dat)
samreaddepth