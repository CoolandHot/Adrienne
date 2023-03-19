# BiocManager::install(c("DESeq2", "ReportingTools", "S4Vectors", "org.Mm.eg.db", "AnnotationDbi"))
# http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#standard-workflow

library(S4Vectors)
# otherwise it will occur `could not find function "mcols"`

setwd("/home/huan/projects/Adrienne")

counts_data <- read.table("Mapped_reads_counts.STAR.txt",
  row.names = 1, sep = ",",
  header = TRUE
)
rownames(counts_data) <- make.names(counts_data$GeneSymbol, unique = TRUE)

sample_info <- read.table("meta-data.csv",
  row.names = 1, sep = ",",
  header = TRUE
) |>
  {
    # remove F group temporary because it is still processing in HTCondor
    \(x) subset(x, row.names(x) != "p22303.s019_GL261.luc2.BHB5.4day.a_S19")
  }() |>
  dplyr::mutate(condition = as.factor(condition))

##################################################
##### Start differential expression analysis
##################################################
counts_data_sub <- subset(counts_data,
  select = row.names(sample_info)
)
# making sure both matches
# identical(names(counts_data_sub), rownames(sample_info))

#########################
##### create DESeq class
#########################
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_data_sub,
  colData = sample_info,
  design = ~ 0 + condition
) |>
  {
    # QC
    \(x) x[rowSums(DESeq2::counts(x)) >= 10, ]
  }()
dds$condition <- relevel(dds$condition, ref = "A")
dds_diff <- DESeq2::DESeq(dds)

# coef(dds_diff) |> head()

vsd <- DESeq2::vst(dds_diff, blind = FALSE)
p1 <- BiocGenerics::plotPCA(vsd, intgroup = c("condition")) +
  ggplot2::labs(title = "PCA for variance")

ggpubr::ggexport(
  filename = "PCA_of_samples_for_overall_variances.pdf",
  plotlist = list(p1), nrow = 1, ncol = 1,
  width = 15, height = 9
)


############# export result
# https://support.bioconductor.org/p/117940/#126943

### Dose response of long term treatment

res <- DESeq2::results(dds_diff,
  contrast = c("condition", "A", "B"),
  alpha = 0.05
) |>
  {
    \(x) x[order(x$padj), ]
  }()

write.csv(as.data.frame(res),
  file = paste0("condition_", "A", "_vs_", "B", ".csv")
)

######################
##### HTML report
######################
# des2Report <- ReportingTools::HTMLReport(
#   shortName = "RNAseq_analysis_with_DESeq2",
#   title = paste0("RNA-seq analysis of differential expression: ", "A", "_vs_", "F"),
#   reportDirectory = paste0("./condition_", "A", "_vs_", "F")
# )
# ReportingTools::publish(dds_diff, des2Report,
#   pvalueCutoff = 0.05,
#   annotation.db = "org.Mm.eg.db",
#   factor = sample_info$condition,
#   reportDir = paste0("./condition_", "A", "_vs_", "F")
# )
# ReportingTools::finish(des2Report)

##################################################
##### plot
##################################################
## Blue dots are genes with p.adjust < alpha(5%)
## triangle on the edges of top & bottom are the genes with logFC > axis.limit
## genes on top right and bottom right are considered to be interested (more counts, higher fold changes)
# DESeq2::plotMA(res)

## plot counts of a certain gene between groups
# library("ggplot2")
# idx <- which.min(res$padj)
# DESeq2::plotCounts(dds_diff,
#   gene = idx,
#   intgroup = "condition",
#   returnData = TRUE
# ) |>
#   dplyr::filter(condition %in% c("A", "F")) |>
#   ggplot(aes(x = condition, y = count)) +
#   geom_point(position = position_jitter(w = 0.1, h = 0)) +
#   labs(
#     x = "sample group",
#     y = "normalised counts",
#     title = paste(
#       "the most significant gene:",
#       counts_data[row.names(res)[idx], "GeneSymbol"]
#     )
#   )
