if(!any(rownames(installed.packages()) == "pacman")){
      install.packages("pacman")
}

if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

pacman::p_load(here,
               DESeq2,
               apeglm,
               Mus.musculus,
               tidyverse,
               gplots)



## Reading in count-data
x <- read_delim(here("counts", "counts_hiseq.tabular"),  col_names = TRUE) |>
      rename(EntrezID = `#KEY` ) |>
      filter(!EntrezID %in% "Geneid") |>
      column_to_rownames(var ="EntrezID")


colnames(x) <- gsub(".fq", "", colnames(x))

### Extracting Gene and mice group
geneid <- rownames(x)

coldata <- read_csv(here("factor","factordata_hidden.csv"))


x_matrix <- x |>
      as.matrix() 



counts_numeric <- apply(x_matrix, 2, as.numeric)
rownames(counts_numeric) <- rownames(x)
dds <- DESeqDataSetFromMatrix(countData = counts_numeric,
                              colData = coldata,
                              design= ~ group)

#Add genes
featureData <- data.frame(gene=geneid)
mcols(dds) <- DataFrame(mcols(dds), featureData)

## Filter 
#keep <- rowSums(counts(dds) >= 10) >= 3
#dds <- dds[keep,]

dds$group <- factor(dds$group, levels = c("control","treatment1"))
dds$group <- relevel(dds$group, ref = "control")

## DE Analysis
dds <- DESeq(dds)
res <- results(dds)
res

## Specific Contrast
res <- results(dds, name="group_treatment1_vs_control")
res <- results(dds, contrast=c("group","treatment1","control"))

resultsNames(dds)

# Reorder by p-value
resOrdered <- res[order(res$pvalue),]

summary(res)
sum(res$padj < 0.1, na.rm=TRUE)

res05 <- results(dds, alpha=0.05)
summary(res05)

sum(res05$padj < 0.05, na.rm=TRUE)
df <- as.data.frame(resOrdered)

write.csv(as.data.frame(resOrdered), 
          file=here("output","DESeq2_results.csv"))


# Data transformations and visualization
#In order to test for differential expression, we operate on raw counts and use discrete distributions as described in the previous section on differential expression. However for other downstream analyses – e.g. for visualization or clustering – it might be useful to work with transformed versions of the count data.
vsd <- vst(dds, blind=FALSE) # Similar to putting data on log2 scale while dealing with sampling variability of low counts


## PCA

pcaData <- plotPCA(vsd, intgroup=c("group"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=group)) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
      coord_fixed()

ggsave(plot = last_plot(), here("figures","PCA_deseq2.jpg"),
       dpi = 300,
       width = 8,
       height = 8)
