# https://support.bioconductor.org/p/117940/#126943

# library(S4Vectors)
# # otherwise it will occur `could not find function "mcols"`
# library(BiocParallel)
# register(MulticoreParam(2))
setwd("~/projects/Arizona_brain_tumor_research")

# dds_diff <- readRDS(file = "./output/diff_expr_A-J_groups_individual.RDS")
dds_diff <- readRDS(file = "./output/diff_expr_A-J_groups_intercept.RDS")
# DESeq2::resultsNames(dds_diff)

# ctrl_sampl <- "A"
# treat_sampl <- "B"

#######################################
### Dose difference, long term
#######################################
lapply(
  list(c("A", "B"), c("A", "C"), c("A", "D"), c("B", "C"), c("B", "D")),
  function(x) {
    ctrl_sampl <- x[1]
    treat_sampl <- x[2]
    res <- DESeq2::results(dds_diff,
      contrast = c("condition", ctrl_sampl, treat_sampl),
      alpha = 0.05,
      lfcThreshold = 1,
      altHypothesis = "greaterAbs"
    ) |>
      {
        \(x) x[order(x$padj), ]
      }()

    write.csv(as.data.frame(res),
      file = paste0(
        "./output/doseDiff_longTerm_",
        ctrl_sampl, "_vs_", treat_sampl, ".csv"
      )
    )
  }
)

#######################################
### Dose difference, short term
#######################################
lapply(
  list(c("A", "E"), c("A", "F"), c("A", "G"), c("E", "F"), c("E", "G")),
  function(x) {
    ctrl_sampl <- x[1]
    treat_sampl <- x[2]
    res <- DESeq2::results(dds_diff,
      contrast = c("condition", ctrl_sampl, treat_sampl),
      alpha = 0.05,
      lfcThreshold = 1,
      altHypothesis = "greaterAbs"
    ) |>
      {
        \(x) x[order(x$padj), ]
      }()

    write.csv(as.data.frame(res),
      file = paste0(
        "./output/doseDiff_shortTerm_",
        ctrl_sampl, "_vs_", treat_sampl, ".csv"
      )
    )
  }
)

#######################################
### same dose, long vs short term
#######################################
lapply(
  list(c("B", "E"), c("C", "F"), c("D", "G")),
  function(x) {
    ctrl_sampl <- x[1]
    treat_sampl <- x[2]
    res <- DESeq2::results(dds_diff,
      contrast = c("condition", ctrl_sampl, treat_sampl),
      alpha = 0.05,
      lfcThreshold = 1,
      altHypothesis = "greaterAbs"
    ) |>
      {
        \(x) x[order(x$padj), ]
      }()

    write.csv(as.data.frame(res),
      file = paste0(
        "./output/termDiff_doseSame_",
        ctrl_sampl, "_vs_", treat_sampl, ".csv"
      )
    )
  }
)


#######################################
### long term, dose removal
#######################################
lapply(
  list(c("B", "H"), c("C", "I"), c("D", "J")),
  function(x) {
    ctrl_sampl <- x[1]
    treat_sampl <- x[2]
    res <- DESeq2::results(dds_diff,
      contrast = c("condition", ctrl_sampl, treat_sampl),
      alpha = 0.05,
      lfcThreshold = 1,
      altHypothesis = "greaterAbs"
    ) |>
      {
        \(x) x[order(x$padj), ]
      }()

    write.csv(as.data.frame(res),
      file = paste0(
        "./output/doseRemoval_longTerm_",
        ctrl_sampl, "_vs_", treat_sampl, ".csv"
      )
    )
  }
)


#######################################
### long term, dose removal, gene expression return to baseline/short term?
#######################################
lapply(
  list(
    c("A", "H"), c("A", "I"), c("A", "J"),
    c("E", "H"), c("F", "I"), c("G", "J")
  ),
  function(x) {
    ctrl_sampl <- x[1]
    treat_sampl <- x[2]
    res <- DESeq2::results(dds_diff,
      contrast = c("condition", ctrl_sampl, treat_sampl),
      alpha = 0.05,
      lfcThreshold = 1,
      altHypothesis = "greaterAbs"
    ) |>
      {
        \(x) x[order(x$padj), ]
      }()

    write.csv(as.data.frame(res),
      file = paste0(
        "./output/doseRemoval_expressReturn_",
        ctrl_sampl, "_vs_", treat_sampl, ".csv"
      )
    )
  }
)


# DESeq2::summary(res)

# ## Blue dots are genes with p.adjust < alpha(5%)
# ## triangle on the edges of top & bottom are the genes with logFC > axis.limit
# ## genes on top right and bottom right are considered to be interested (more counts, higher fold changes)
# DESeq2::plotMA(res, main = "Log fold change vs normalised counts")
# abline(h = c(-1, 1), col = "dodgerblue", lwd = 2)


# resLFC <- DESeq2::lfcShrink(dds_diff,
#   coef = "condition_B_vs_A",
#   type = "apeglm",
#   parallel = TRUE
# )
# resAsh <- DESeq2::lfcShrink(dds_diff,
#   contrast = c("condition", "A", "B"),
#   type = "ashr",
#   parallel = TRUE
# )
# resNorm <- DESeq2::lfcShrink(dds_diff,
#   contrast = c("condition", "A", "B"),
#   type = "normal",
#   parallel = TRUE
# )

# DESeq2::plotMA(resLFC, main = "apeglm shrinkage")
# abline(h = c(-1, 1), col = "dodgerblue", lwd = 2)

##################################################
############## plot counts of a certain gene between groups
##################################################
# library("ggplot2")
# idx <- which.min(res$padj)
# DESeq2::plotCounts(dds_diff,
#   gene = idx,
#   intgroup = "condition",
#   returnData = TRUE
# ) |>
#   dplyr::filter(condition %in% c(ctrl_sampl, treat_sampl)) |>
#   ggplot(aes(x = condition, y = count)) +
#   geom_point(position = position_jitter(w = 0.1, h = 0)) +
#   labs(
#     x = "sample group",
#     y = "normalised counts",
#     title = paste(
#       "the most significant differentially expressed gene:",
#       # in case this gene has a duplicated name which was changed in make.names()
#       sub("\\.\\d+", "", row.names(res)[idx])
#     )
#   )
