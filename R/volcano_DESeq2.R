if(!any(rownames(installed.packages()) == "pacman")){
      install.packages("pacman")
}

pacman::p_load(tidyverse,
               RColorBrewer,
               ggrepel,
               here)

DE <- read.csv(here("output","DESeq2_results.csv")) 
rownames(DE) <- DE$X


# Theme 
theme_set(theme_classic(base_size = 20) +
                theme(
                      axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = "black"),
                      axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = "black"),
                      plot.title = element_text(hjust = 0.5)
                ))


## Set logFC and pvalue cut-off
DE$diffexpressed <- "NO"
DE$diffexpressed[DE$log2FoldChange > 1.5 & DE$padj < 0.05] <- "UP"
DE$diffexpressed[DE$log2FoldChange < -1.5 & DE$padj < 0.05] <- "DOWN"

## GGrepel top 30 lipid species
top30deg <- head(DE[order(DE$padj), "X"], 30) 


DE$delabel <- ifelse(DE$X %in% top30deg, DE$X, NA)
## Plot
volcano_plot <- ggplot(data = DE, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed, label = delabel)) +
      geom_vline(xintercept = c(-1.5, 1.5), col = "gray", linetype = "dashed") +
      geom_hline(yintercept = c(-log10(0.05)), col = "gray", linetype = "dashed") +
      geom_point(size = 2) +
      scale_color_manual(values = c("#00AFBB", "grey","#bb0c00"),
                         labels = c("Downregulated", "Not Significant", "Upregulated")) +
      coord_cartesian(ylim = c(0, 10), xlim = c(-10, 10)) +
      scale_x_continuous(breaks = seq(-10, 10, 1)) +
      labs(color = "Treatment1 vs Control", x = "logFC", y = expression("-log"[10]*"(Adjusted p-value)")) +
      ggtitle(NULL) +
      geom_text_repel(max.overlaps = Inf)
volcano_plot

DE %>%
      count(diffexpressed)

ggsave(filename = here("figures","volcano_plot_DESeq2.jpg"),
       dpi = 300, width = 12, height = 10,
       plot = volcano_plot)

