
#Rithik Castelino (Walters Lab, NCI)

#Loading in packages if already installed. If not already installed-rerun code again to load.
if(!require(tidyverse)) { install.packages("tidyverse"); library(tidyverse)}
if(!require(ggrepel)) { install.packages("ggrepel"); library(ggrepel)}
if(!require(scales)) { install.packages("scales"); library(scales)}
if(!require(ggh4x)) { install.packages("ggh4x"); library(ggh4x)}

#Installing BiocManager to get TissueEnrich from Bioconductor. Un-comment lower-line to install.
if (!require(BiocManager)) { install.packages("BiocManager"); library(BiocManager)}
#BiocManager::install(version = "3.19")

#Loading in TissueEnrich if already installed. Un-comment lower-line to install.
if (!require(TissueEnrich)) { install.packages("TissueEnrich"); library(TissueEnrich)}
#BiocManager::install("TissueEnrich", force=TRUE)

#Loading in data and selecting appropriate gene-edited cell lines.
TMTMS_rawdata <- read.csv("./Input/Abundance_vs_pvalues.csv", skip = 1)
TMTMS_data <- TMTMS_rawdata %>% 
  dplyr::select(., 1, 4, 5) %>% 
  dplyr::rename(., Log2Abundance = Log2_.Ab) %>% 
  dplyr::rename(., Neglog10pvalue = Neg_log10..pvalue..1)

#Extracting gene names.
GeneNames <- unique(tail(TMTMS_data$Gene, -1))

#Processing set of Genes through the 'teEnrichment' function of TissueEnrich package.
gs <- GeneSet(geneIds=GeneNames, organism="Homo Sapiens", geneIdType=SymbolIdentifier())
output<-teEnrichment(inputGenes=gs, #using list of gene names extracted previously
	rnaSeqDataset = 1, #Using value of '1' to specify Human Protein Atlas'
	tissueSpecificGeneType = 1) #Using '1' to group all tissue specific genes under an aggregate category

#Initializing empty list
TissueSpecificGenes <- list()

#Brute force population of list with gene names of tissue specific proteins 
for (i in 1:length(output[[3]])) {
  enrichmentOutput<-data.frame(assay(output[[3]][[i]]))$Gene
  TissueSpecificGenes <- c(TissueSpecificGenes, enrichmentOutput)
}

#Extracting uncategorized gene Ids for reference
Uncategorized <- geneIds(output[[4]])

#Adding on a new column to the 'TMTMS_data' defining if a protein is Tissue Specific, not Tissue Specific, or uncategorized.
TMTMS_data_labelled <- TMTMS_data %>% 
  mutate(., TissueSpecificity = ifelse(Gene %in% TissueSpecificGenes, "Yes", "No")) %>% 
  mutate(., TissueSpecificity = ifelse(Gene %in% Uncategorized, "NA", TissueSpecificity))

#Generating graph
graph <- ggplot(data=TMTMS_data_labelled, aes(x=Log2Abundance, y=Neglog10pvalue, colour = TissueSpecificity, label=Gene)) +
  theme_classic() +
  geom_point(data = TMTMS_data_labelled[which(abs(TMTMS_data_labelled$Log2Abundance) >= 1 & TMTMS_data_labelled$Neglog10pvalue >= 1.301),], size = 2.5) +
  geom_point(data = TMTMS_data_labelled[which(abs(TMTMS_data_labelled$Log2Abundance) < 1 & TMTMS_data_labelled$Neglog10pvalue >= 1.301),], size = 1) +
  geom_point(data = TMTMS_data_labelled[which(TMTMS_data_labelled$Neglog10pvalue < 1.301),], size = 1, colour = "grey") +
  geom_text_repel(data = TMTMS_data_labelled[which((TMTMS_data_labelled$Log2Abundance) > 1 & TMTMS_data_labelled$Neglog10pvalue > 1.301),], size = 3.5, max.overlaps = 3, min.segment.length = 1e-6, force_pull = 0.3, aes(family = "Arial Narrow"), segment.size = unit(0.08, "cm"), box.padding = 0.2, segment.curvature = -0.3, segment.ncp = 10, segment.angle = 20, nudge_y = 0.2, nudge_x = 0.2) +
  geom_text_repel(data = TMTMS_data_labelled[which((TMTMS_data_labelled$Log2Abundance) < -1 & TMTMS_data_labelled$Neglog10pvalue > 1.301),], size = 3.5, max.overlaps = 3, min.segment.length = 1e-6, force_pull = 0.3, aes(family = "Arial Narrow"), segment.size = unit(0.08, "cm"), box.padding = 0.2, segment.curvature = -0.3, segment.ncp = 10, segment.angle = 20, nudge_y = 0.2, nudge_x = -0.2) +
  geom_vline(xintercept = c(-1, 1), linetype="dotted", colour = "black") +
  geom_hline(yintercept = c(1.301), linetype="dotted", colour = "black") +
  scale_x_continuous(limits = c(-6, 6), breaks=seq(-6, 6, 1), guide="axis_minor", expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 7), breaks=seq(0, 7, 1), guide="axis_minor", expand = c(0,0)) +
  labs(x= "Log2(∆UIM/WT)", y="-Log10(p-value)") +
  theme(legend.position = "None") + 
  scale_colour_discrete(breaks=c("NA", "No", "Yes"), type = c("black", "cornflowerblue", "maroon"))
graph

#Saving graph as a '.tiff' file.
ggsave("./Output/tissuespecificity_∆UIM.tiff", graph, width = 3456/250, height = 2234/(250/1.547), dpi = 300, bg = "transparent")

#Calculating relative percentage of Tissue Specific proteins within certain thresholds. See Table S1 in Lu & Osei-Amponsa et al., 2025.
Num_Specific_less1 <- nrow(TMTMS_data_labelled[which(abs(TMTMS_data_labelled$Log2Abundance) <= 1 & TMTMS_data_labelled$TissueSpecificity == "Yes" & TMTMS_data_labelled$Neglog10pvalue > 1.301),])
Num_Total_less1 <- nrow(TMTMS_data_labelled[which(abs(TMTMS_data_labelled$Log2Abundance) <= 1 & TMTMS_data_labelled$TissueSpecificity != "NA" & TMTMS_data_labelled$Neglog10pvalue > 1.301),])
Percent_Specific_Total_less1 = Num_Specific_less1/Num_Total_less1

Num_Specific_more1 <- nrow(TMTMS_data_labelled[which(abs(TMTMS_data_labelled$Log2Abundance) >= 1 & TMTMS_data_labelled$TissueSpecificity == "Yes" & TMTMS_data_labelled$Neglog10pvalue > 1.301),])
Num_Total_more1 <- nrow(TMTMS_data_labelled[which(abs(TMTMS_data_labelled$Log2Abundance) >= 1 & TMTMS_data_labelled$TissueSpecificity != "NA" & TMTMS_data_labelled$Neglog10pvalue > 1.301),])
Percent_Specific_Total_more1 = Num_Specific_more1/Num_Total_more1

Num_Specific_more2 <- nrow(TMTMS_data_labelled[which(abs(TMTMS_data_labelled$Log2Abundance) >= 2 & TMTMS_data_labelled$TissueSpecificity == "Yes" & TMTMS_data_labelled$Neglog10pvalue > 1.301),])
Num_Total_more2 <- nrow(TMTMS_data_labelled[which(abs(TMTMS_data_labelled$Log2Abundance) >= 2 & TMTMS_data_labelled$TissueSpecificity != "NA" & TMTMS_data_labelled$Neglog10pvalue > 1.301),])
Percent_Specific_Total_more2 = Num_Specific_more2/Num_Total_more2

#_____________________________________________
