if(!any(rownames(installed.packages()) == "pacman")){
      install.packages("pacman")
}

if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

pacman::p_load(here,
               org.Mm.eg.db,
               clusterProfiler,
               enrichplot,
               tidyverse,
               readxl,
               forcats,
               DOSE)


tt <- as.data.frame(read_csv(here("output","DESeq2_results.csv")) |>
                          rename(gene = ...1))
tt$fcsign <- sign(tt$log2FoldChange)
tt$logP<- -log10(tt$pvalue)
tt$metric <- tt$logP/tt$fcsign
y<-tt[,c("gene", "metric")]
filtered <- na.omit(y)

######  How to prepare your own geneList #####
d <- as.data.frame(read_csv(here("output","DESeq2_results.csv")) |>
                         rename(gene = ...1) |>
                         select(gene, log2FoldChange))
filtered$gene <- mapIds(org.Mm.eg.db, filtered$gene, 'ENTREZID', 'SYMBOL')
names(filtered$gene) <- NULL
filtered <- na.omit(filtered)
## assume 1st column is ID
## 2nd column is FC

## feature 1: numeric vector
geneList = filtered[,2]

## feature 2: named vector
names(geneList) <- filtered[,1]

## feature 3: decreasing order (Higher FC in Treatment to lower FC)
geneList = geneList[order(geneList, decreasing = TRUE)]

################################################

## GSEA GO : BP
BPgo <- gseGO(geneList     = geneList,
              OrgDb        = org.Mm.eg.db,
              ont          = "BP",
              pvalueCutoff = 0.05,
              verbose      = FALSE)
head(BPgo)
dotplot(BPgo, showCategory=20)
goplot(BPgo, showCategory = 10)



## NES
BPgo2 <- arrange(BPgo, desc(abs(NES))) %>%
      group_by(sign(NES)) %>%
      slice(1:5)

BPgo_plot <- ggplot(BPgo2, showCategory = 10,
                    aes(NES,
                        fct_reorder(Description, NES))) +
      geom_segment(aes(xend=0, yend = Description)) +
      geom_point(aes(color=qvalue, size = Count)) +
      scale_color_gradientn(colors=c("#f7ca64", "#46bac2",
                                     "#7e62a3"),
                            trans = "log10",
                            guide=guide_colorbar(reverse=TRUE,
                                                 order=1)) +
      scale_size_continuous(range=c(2, 10)) +
      theme_dose(12) +
      xlab("Normalized Enrichment Score") +
      ylab(NULL) +
      ggtitle("GO: Biological Process")
BPgo_plot

ggsave(plot = BPgo_plot, filename = here("figures","GO_biological_process_10min.jpg"),
       dpi = 300,
       width = 8,
       height = 8,
       scale =1.3
)


## GO BP Network
## convert gene ID to Symbol
BPgo_net <- setReadable(BPgo, 'org.Mm.eg.db', 'ENTREZID')
p1 <- cnetplot(BPgo_net, categorySize="pvalue", color.params = list(foldChange=geneList))
ggsave(plot = p1, here("figures","GSEA_biological_process_network.jpg"),
       dpi = 300,
       width = 8,
       height = 8,
       scale = 2)


### Heat Plot
p1_heat <- heatplot(BPgo_net, foldChange=geneList, showCategory=5)

### Tree Plot
BPgo_net2 <- pairwise_termsim(BPgo_net)
treeplot(BPgo_net2, hclust_method = "average")
#### GO Molecular Function
## GSEA GO : MF
MFgo <- gseGO(geneList     = geneList,
              OrgDb        = org.Mm.eg.db,
              ont          = "MF",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
head(MFgo)
dotplot(MFgo, showCategory=20)
goplot(MFgo, showCategory = 10)



## NES
MFgo2 <- arrange(MFgo, desc(abs(NES))) %>%
      group_by(sign(NES)) %>%
      slice(1:5)

MFgo_plot <- ggplot(MFgo2, showCategory = 10,
                    aes(NES,
                        fct_reorder(Description, NES))) +
      geom_segment(aes(xend=0, yend = Description)) +
      geom_point(aes(color=qvalue, size = Count)) +
      scale_color_gradientn(colors=c("#f7ca64", "#46bac2",
                                     "#7e62a3"),
                            trans = "log10",
                            guide=guide_colorbar(reverse=TRUE,
                                                 order=1)) +
      scale_size_continuous(range=c(2, 10)) +
      theme_dose(12) +
      xlab("Normalized Enrichment Score") +
      ylab(NULL) +
      ggtitle("GO: Molecular Function")
MFgo_plot

ggsave(plot = MFgo_plot, filename = here("figures","GO_molecular_function.jpg"),
       dpi = 300,
       width = 8,
       height = 8,
       scale =1.3
)


## convert gene ID to Symbol
MFgo_net <- setReadable(MFgo, 'org.Mm.eg.db', 'ENTREZID')
p2 <- cnetplot(MFgo_net, categorySize="pvalue", color.params = list(foldChange=geneList))
ggsave(plot = p2 , here("figures","GSEA_molecular_function_network.jpg"),
       dpi = 300,
       width = 8,
       height = 8,
       scale = 2.5)

## KEGG pathway gene set enrichment analysis
kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'mmu',
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)
dotplot(kk2, showCategory=20)
goplot(kk2, showCategory = 10)


## NES
kk2_nes <- arrange(kk2, desc(abs(NES))) %>%
      group_by(sign(NES)) %>%
      slice(1:5)

kegg_plot <- ggplot(kk2_nes, showCategory = 10,
                    aes(NES,
                        fct_reorder(Description, NES))) +
      geom_segment(aes(xend=0, yend = Description)) +
      geom_point(aes(color=qvalue, size = Count)) +
      scale_color_gradientn(colors=c("#f7ca64", "#46bac2",
                                     "#7e62a3"),
                            trans = "log10",
                            guide=guide_colorbar(reverse=TRUE,
                                                 order=1)) +
      scale_size_continuous(range=c(2, 10)) +
      theme_dose(12) +
      xlab("Normalized Enrichment Score") +
      ylab(NULL) +
      ggtitle("KEGG")
kegg_plot

ggsave(plot = kegg_plot, filename = here("figures","kegg.jpg"),
       dpi = 300,
       width = 8,
       height = 8,
       scale =1.3
)


## GO BP Network
## convert gene ID to Symbol
kk2_net <- setReadable(kk2, 'org.Mm.eg.db', 'ENTREZID')
p1 <- cnetplot(kk2_net, categorySize="pvalue", color.params = list(foldChange=geneList))
ggsave(plot = p1, here("figures","kegg_network.jpg"),
       dpi = 300,
       width = 8,
       height = 8,
       scale = 2)
