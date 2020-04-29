#' Run differential expression
#'
#' @param counts A matrix of counts
#' @param groups conditions
#' @param method DE method to use
#' @param row_filter filter for rows
#' @param col_filter filter for columns
#' @param ... plotting parameters
#' @return \code{combined_out}
#' @examples
#'
#'
#' #load("data/counts_data.Rdata")
#' counts <- matrix(  2 ^ rpois(1e5, 5) , ncol=10, nrow=1000)
#' rownames(counts) <- paste0('gene', 1:1000 )
#' colnames(counts) <- paste0('sample', 1:10 )
#'
#' groups <- c( rep(1,5), rep(2,5) )
#' method <- 'wilcox'
#'
#' degs <- calc_DE(counts, groups, method)
#'
#' @import edgeR DESeq2 stats utils
#' @export
#'

calc_DE <- function(counts, groups, method, row_filter = FALSE,
                    col_filter = FALSE, ...) {

    cpm <- edgeR::cpm(counts)

    if (row_filter) {
        cpm <- cpm[row_filter, ]
        if (counts) {
            counts <- counts[row_filter, ]
        }
    }
    if (col_filter) {
        cpm <- cpm[, col_filter]
        groups <- groups[col_filter]
        if (counts) {
            counts <- counts[, col_filter]
        }
    }

    if (method == "edgeR") {
        groups_adj <- groups
        groups_adj[groups == 0] <- 2
        groups_adj[groups == 2] <- 0
        y <- edgeR::DGEList(counts = counts, group = groups_adj )
        y <- edgeR::estimateGLMCommonDisp(y)
        design <- model.matrix(~groups_adj )
        fit <- edgeR::glmFit(y, design)
        output <- edgeR::glmLRT(fit, coef = 2)
        padj <- p.adjust(output$table$PValue )
        degs <- data.frame(  mean_cpm = 2^output$table$logCPM ,
                             log2_fc = output$table$logFC,
                             pvals = output$table$PValue,
                             padj = padj )
    }
    if (method == "DESeq2") {
        conditions <- groups
        samples <- colnames(cpm)
        col_data <- as.data.frame(cbind(samples, conditions))
        colnames(col_data) <- c("samples", "conditions")
        dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                              colData = col_data,
                                              design = ~conditions)
        dds <- DESeq2::DESeq(dds)
        output <- DESeq2::results(dds, contrast = c("conditions", "1", "2"))
        degs <- data.frame(  mean_cpm = output$baseMean,
                             log2_fc = output$log2FoldChange,
                             pvals = output$pvalue,
                             padj =  output$padj )
    }

    if (method == "wilcox" | method == "") {
        if (sum(groups == 1) < 2) {
            m_cpm1 <- (cpm[, groups == 1])
        } else {
            m_cpm1 <- rowMeans(cpm[, groups == 1])
        }
        if (sum(groups == 2) < 2) {
            m_cpm2 <- (cpm[, groups == 2])
        } else {
            m_cpm2 <- rowMeans(cpm[, groups == 2])
        }
        m_cpm <- rowMeans(cpm)
        logfc <- log2(m_cpm1 / m_cpm2)
        n_genes <- dim(cpm)[1]

        cpm_ps_g <- sapply(1:n_genes, function(k)
            wilcox.test(cpm[k, groups == 1], cpm[k, groups == 2],
                        alt = "g")$p.val)

        cpm_padj_g <- p.adjust(cpm_ps_g, method = "BH")

        cpm_ps_l <- sapply(1:n_genes, function(k)
            suppressWarnings( wilcox.test(cpm[k, groups == 1], cpm[k, groups == 2],
                        alt = "l")$p.val))

        cpm_padj_l <- p.adjust(cpm_ps_l, method = "BH")

        cpm_ps <- sapply(1:n_genes, function(k)
            suppressWarnings(wilcox.test(cpm[k, groups == 1], cpm[k, groups == 2])$p.val) )

        cpm_padj <- p.adjust(cpm_ps, method = "BH")

        degs <- data.frame( mean_cpm = m_cpm,
                            log2_fc = logfc,
                            pvals = cpm_ps,
                            padj = cpm_padj)

        output <- data.frame(cbind(m_cpm, logfc, cpm_ps, cpm_padj,
                      cpm_ps_g, cpm_padj_g,
                      cpm_ps_l, cpm_padj_l,
                      m_cpm1, m_cpm2) )

        rownames(degs) <- rownames(cpm)
        rownames(output) <- rownames(cpm)

    }
    combined_out <- list(degs, output)
    names(combined_out) = c("degs", "output")
    return( combined_out )
}
