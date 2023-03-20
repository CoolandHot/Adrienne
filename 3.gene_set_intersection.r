# install.packages("ggVennDiagram")
setwd("~/projects/Arizona_brain_tumor_research")

# dds_diff <- readRDS(file = "./output/diff_expr_A-J_groups_individual.RDS")
dds_diff <- readRDS(file = "./output/diff_expr_A-J_groups_intercept.RDS")
# DESeq2::resultsNames(dds_diff)

# res <- DESeq2::results(dds_diff,
#     contrast = c("condition", "A", "B"),
#     alpha = 0.05,
#     lfcThreshold = 1,
#     altHypothesis = "greaterAbs"
# )
# DESeq2::plotCounts(dds_diff,
#     gene = which(row.names(res) == "B3galt2"),
#     # gene = which.min(res$padj),
#     intgroup = "condition"
# )
# res[which.min(res$padj), ]

#######################################
### Dose difference, long term
### list(c("A", "B"), c("A", "C"), c("A", "D"), c("B", "C"), c("B", "D"))
#######################################
dose.long.ab <- read.csv(file = paste0(
    "./output/doseDiff_longTerm_",
    "A", "_vs_", "B", ".csv"
)) |>
    dplyr::filter(padj < 0.05) |>
    dplyr::mutate(log2FoldChange = -log2FoldChange, stat = -stat)
dose.long.ab.up <- dplyr::filter(dose.long.ab, log2FoldChange > 0)
dose.long.ab.down <- dplyr::filter(dose.long.ab, log2FoldChange < 0)

dose.long.ac <- read.csv(file = paste0(
    "./output/doseDiff_longTerm_",
    "A", "_vs_", "C", ".csv"
)) |>
    dplyr::filter(padj < 0.05) |>
    dplyr::mutate(log2FoldChange = -log2FoldChange, stat = -stat)
dose.long.ac.up <- dplyr::filter(dose.long.ac, log2FoldChange > 0)
dose.long.ac.down <- dplyr::filter(dose.long.ac, log2FoldChange < 0)

dose.long.ad <- read.csv(file = paste0(
    "./output/doseDiff_longTerm_",
    "A", "_vs_", "D", ".csv"
)) |>
    dplyr::filter(padj < 0.05) |>
    dplyr::mutate(log2FoldChange = -log2FoldChange, stat = -stat)
dose.long.ad.up <- dplyr::filter(dose.long.ad, log2FoldChange > 0)
dose.long.ad.down <- dplyr::filter(dose.long.ad, log2FoldChange < 0)

dose.long.bc <- read.csv(file = paste0(
    "./output/doseDiff_longTerm_",
    "B", "_vs_", "C", ".csv"
)) |>
    dplyr::filter(padj < 0.05) |>
    dplyr::mutate(log2FoldChange = -log2FoldChange, stat = -stat)
dose.long.bc.up <- dplyr::filter(dose.long.bc, log2FoldChange > 0)
dose.long.bc.down <- dplyr::filter(dose.long.bc, log2FoldChange < 0)

dose.long.bd <- read.csv(file = paste0(
    "./output/doseDiff_longTerm_",
    "B", "_vs_", "D", ".csv"
)) |>
    dplyr::filter(padj < 0.05) |>
    dplyr::mutate(log2FoldChange = -log2FoldChange, stat = -stat)
dose.long.bd.up <- dplyr::filter(dose.long.bd, log2FoldChange > 0)
dose.long.bd.down <- dplyr::filter(dose.long.bd, log2FoldChange < 0)

########
ggVennDiagram::ggVennDiagram(list(
    "D vs B" = dose.long.bd.up$X,
    "C vs B" = dose.long.bc.up$X,
    "C vs A" = dose.long.ac.up$X,
    "D vs A" = dose.long.ad.up$X
))
Reduce(intersect, list(
    dose.long.bd.up$X, dose.long.bc.up$X, dose.long.ac.up$X, dose.long.ad.up$X
))
dose.long.a2bcd <- Reduce(intersect, list(
    dose.long.ab.up$X, dose.long.ac.up$X, dose.long.ad.up$X
))

intersect(dose.long.bd.up$X, dose.long.bc.up$X)


#######################################
### Dose difference, short term
#######################################
list(c("A", "E"), c("A", "F"), c("A", "G"), c("E", "F"), c("E", "G"))
read.csv(
    file = paste0(
        "./output/doseDiff_shortTerm_",
        ctrl_sampl, "_vs_", treat_sampl, ".csv"
    )
)
#######################################
### same dose, long vs short term
#######################################
list(c("B", "E"), c("C", "F"), c("D", "G"))
read.csv(
    file = paste0(
        "./output/termDiff_doseSame_",
        ctrl_sampl, "_vs_", treat_sampl, ".csv"
    )
)


#######################################
### long term, dose removal
#######################################
list(c("B", "H"), c("C", "I"), c("D", "J"))
read.csv(
    file = paste0(
        "./output/doseRemoval_longTerm_",
        ctrl_sampl, "_vs_", treat_sampl, ".csv"
    )
)


#######################################
### long term, dose removal, gene expression return to baseline/short term?
#######################################
list(
    c("A", "H"), c("A", "I"), c("A", "J"),
    c("E", "H"), c("F", "I"), c("G", "J")
)
read.csv(
    file = paste0(
        "./output/doseRemoval_expressReturn_",
        ctrl_sampl, "_vs_", treat_sampl, ".csv"
    )
)
