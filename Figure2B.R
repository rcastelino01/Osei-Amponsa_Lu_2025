
#Rithik Castelino (Walters Lab, NCI)

#Loading in packages if already installed. If not already installed-rerun code again to load.
if(!require(tidyverse)) { install.packages("tidyverse"); library(tidyverse)}
if(!require(readr)) { install.packages("readr"); library(readr)}

#Loading in Consense RNAseq from The Human Protein Atlas v24.0 (https://www.proteinatlas.org/about/download)
consensus_rnaseq <- readr::read_tsv("./Input/rna_tissue_consensus.tsv")
#Loading in TMT-MS data
TMTMS_rawdata <- read.csv("./Input/Abundance_vs_pvalues.csv", skip = 1)

#Loading in TMT-MS data and selecting TMT-MS data for ∆UIM.
TMTMS_data <- TMTMS_rawdata %>%
  dplyr::select(., 1, 4, 5) %>% 
  dplyr::rename(., Log2Abundance = Log2_.Ab) %>% 
  dplyr::rename(., Neglog10pvalue = Neg_log10..pvalue..1) %>% 
  filter(abs(.$Log2Abundance) >= 2 & .$Neglog10pvalue >= 1.3) #Filtering for proteins which are significant and have 4x log fold change.
  
#Filtering for proteins present in both data sets
consensus_rnaseq_filtered <- consensus_rnaseq %>% 
  filter(.$`Gene name` %in% TMTMS_data$Gene)

#Extracting the max nTPM value across all tissues for each protein
consensus_rnaseq_summary <- consensus_rnaseq_filtered %>% 
  group_by(`Gene name`) %>% 
  summarise(max = max(nTPM))

#Normalizing each nTPM value as a function of tissue relative to the max nTPM across-tissues for each protein
consensus_rnaseq_normtomax <- consensus_rnaseq_filtered %>% 
  inner_join(., consensus_rnaseq_summary) %>% 
  mutate(norm_nTPM = nTPM/max)

#Generating graph
graph <- ggplot(data = consensus_rnaseq_normtomax, aes(x=reorder(`Gene name`, norm_nTPM), y=fct_rev(Tissue), fill=norm_nTPM)) +
  theme_classic() +
  geom_tile(color = "black", linewidth = 0.5) +
  scale_fill_gradient2(low = '#fdfefe', mid = "#fdca69", high="#ff0000", midpoint=0.5) +
  theme(axis.text = element_text(size=5, colour = "black"), 
        axis.title = element_text(size=5, colour = "black"), 
        axis.title.x = element_text(margin = margin(t = 2)),
        axis.title.y = element_text(margin = margin(r = 2)),
        axis.line = element_line(size = 0.25, colour = "black"),
        axis.ticks.length = unit(0.087, "cm"),
        axis.ticks = element_line(size = 0.087),
        text = element_text(family="Arial Narrow"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.background = element_blank(), panel.background = element_rect(fill = 'transparent'), 
        legend.background = element_rect(fill='transparent'))

#Displaying Graph
graph

#Saving graph as a '.tiff' file. Duplicate and modify as needed for other file types (e.g. '.svg')
ggsave("./Output/Figure2B.tiff", graph)

#----------------------------------------
