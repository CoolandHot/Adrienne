# BiocManager::install(c("DESeq2", "ReportingTools", "S4Vectors", "org.Mm.eg.db", "AnnotationDbi", "apeglm"))
# devtools::install_github("stephens999/ashr")

# http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#standard-workflow

library(S4Vectors)
# otherwise it will occur `could not find function "mcols"`
library(BiocParallel)
register(MulticoreParam(2))

setwd("~/projects/Arizona_brain_tumor_research")

##################### steps ###################
# 1. read gene count table and meta-data for each sample
#   make sure the row names for meta-data are identical to the column names in gene count table
# 2. create a DESeq2 class and run through `DESeq2::DESeq`
#    `DESeq2::design` is for GLM formula:
# https://docs.google.com/presentation/d/1B9zW1_F-kBqQEu4xqxIJrudYP5DecytYMRR6bY4H6aM/edit#slide=id.gaddde84de4_0_211
#            `~condition`: β0+β1*x1+β2*x2+...
#            `~1+condition`: β0+β1*x1+β2*x2+...
#            `~0+condition`: β1*x1+β2*x2+...
# 3. PCA to see variance of each group (how good are the samples sequenced)
##################### steps ###################

counts_data <- read.table("Mapped_reads_counts.STAR.txt",
    row.names = 1, sep = ",",
    header = TRUE
)
rownames(counts_data) <- make.names(counts_data$GeneSymbol, unique = TRUE)

sample_info <- read.table("meta-data.csv",
    row.names = 1, sep = ",",
    header = TRUE
) |>
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
dds$condition <- relevel(dds$condition, ref = ctrl_sampl)
dds_diff <- DESeq2::DESeq(dds, parallel = TRUE)
saveRDS(dds_diff, file = "./output/diff_expr_A-J_groups_individual.RDS")

# DESeq2::resultsNames(dds_diff)
# coef(dds_diff) |> head()

vsd <- DESeq2::vst(dds_diff, blind = FALSE)
p1 <- BiocGenerics::plotPCA(vsd, intgroup = c("condition")) +
    ggplot2::labs(title = "PCA for variance")
ggpubr::ggexport(
    filename = "./output/PCA_of_samples_for_overall_variances.pdf",
    plotlist = list(p1), nrow = 1, ncol = 1,
    width = 15, height = 9
)

DESeq2::design(dds) <- ~condition
dds_diff <- DESeq2::DESeq(dds, parallel = TRUE)
saveRDS(dds_diff, file = "./output/diff_expr_A-J_groups_intercept.RDS")


####################################
##### HTML report for all genes
####################################
# des2Report <- ReportingTools::HTMLReport(
#   shortName = "RNAseq_analysis_with_DESeq2",
#   title = paste0("RNA-seq analysis of differential expression: ", ctrl_sampl, "_vs_", treat_sampl),
#   reportDirectory = paste0("./output/condition_", ctrl_sampl, "_vs_", treat_sampl)
# )
# ReportingTools::publish(dds_diff, des2Report,
#   pvalueCutoff = 0.05,
#   annotation.db = "org.Mm.eg.db",
#   factor = sample_info$condition,
#   reportDir = paste0("./output/condition_", ctrl_sampl, "_vs_", treat_sampl)
# )
# ReportingTools::finish(des2Report)
