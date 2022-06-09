#' Get sub-network
#'
#' @param deg A matrix (differential expression analysis)
#' @param network_type name of pre-built network to use (generic, blood, brain)
#' @param method top (fold change only) or reg (fold change and FDR cutoff)
#' @param n_top number of genes to consider
#' @param q q-value
#' @param fct abs fold change threshold
#' @param flag_occr boolean to select occurence network or coexpression network
#' @param dir dummy directory for the networks
#' @return \code{combined_output}
#' @examples
#'
#' i <- 1
#' studies <- paste0('study', 1:10 )
#' nettype <- 'random'
#' dir <- 'out_test'
#' label <- studies[i]
#' filename <- 'output'
#' genes <- rhdf5::h5read('/directflow/SCCGGroupShare/projects/sarba2/data/agg_coexp/generic.genes.h5',
#'  'genes' )
#'
#'
#' counts <- matrix(  2 ^ rpois(1e5, 5) , ncol=10, nrow=1000)
#' rownames(counts) <- sample(genes, 1000 )
#' colnames(counts) <- paste0('sample', 1:10 )
#'
#' groups <- c( rep(1,5), rep(2,5) )
#' method <- 'wilcox'
#'
#' deg <- calc_DE(counts, groups, method)
#' dir <- '/directflow/SCCGGroupShare/projects/sarba2/data/agg_coexp/'
#' dir <- 'S:/data/agg_coexp/'
#'
#' network_type='generic'
#' subset_network_hdf5(deg$degs, network_type, dir=dir)
#'
#' @import rhdf5 stats utils
#' @export
#'

subset_network_hdf5 <- function(deg, network_type = "generic",
                                method = "top", n_top = 100, min_cpm = 1,
                                q = 0.05, fct = 2, flag_occr = TRUE,
                                dir = "") {


    if (network_type == "generic" | network_type == "generic230"
        | network_type == "generic75neg"
        | network_type == "blood" | network_type == "brain") {
        if (flag_occr == T) {
            network_type <- paste0(network_type, ".occr")
        }

        genes_hdf5 <- paste0(dir, network_type, ".genes.h5")
        median_hdf5 <- paste0(dir, network_type, ".med.h5")

        genes <- rhdf5::h5read(genes_hdf5, "genes")
        colnames(genes) <- c("entrezID", "name", "ensemblID")

        median_net <- rhdf5::h5read(median_hdf5, "median")

        net_hdf5 <- paste0(dir, network_type, ".net.h5")

        # Check genes (probably should make this a sep function)
        grep_res <- grep ("^ENSG", head ( rownames(deg)  )  )

        if( length(grep_res ) >0 ) {
            # set genes to ensemblID
            genes <- genes[,3]

        } else {
            # check if entrez IDs
             grep_res <- grep ("^\\d+", head ( rownames(deg)  )  )

             # if true, set to entrez ID
             if( length(grep_res ) >0 ) {
                 genes <- genes[,1]
             # if not, set to gene symbols/names
             } else {
                genes <- genes[,2]
             }
        }
    }

    deg <- deg[is.finite(deg$log2_fc) ,]



    m <- match(rownames(deg), genes)
    f_m <- !is.na(m)
    f_am <- m[f_m]

    deg_sub <- deg[f_m,]

    m_cpm <- log2(deg_sub[, 1])
    fc    <- deg_sub[, 2]
    padj  <- p.adjust(deg_sub[,3], method="fdr" )

    sub_net <- list()
    deg_sig <- list()
    fc_sig <- list()
    node_degrees <- list()

    if (method == "top") {
        fc_temp = fc
        fc_temp[ m_cpm <= log2(min_cpm)  ] = 0
        deg_sig[["up"]] <- tail(order(fc_temp), n = n_top)
        deg_sig[["down"]] <- head(order(fc_temp), n = n_top)
    }
    if (method == "reg") {
        deg_sig[["down"]] <- which(fc <= -fct & padj <= q & m_cpm > log2(min_cpm) )
        deg_sig[["up"]] <- which(fc >= fct & padj <= q & m_cpm > log2(min_cpm) )
    }

    fc_sig[["down"]] <- cbind(m_cpm, fc, padj)[deg_sig[["down"]], ]
    fc_sig[["up"]] <- cbind(m_cpm, fc, padj)[deg_sig[["up"]], ]

    genes_index_down <- f_am[deg_sig[["down"]]]
    genes_index_up <- f_am[deg_sig[["up"]]]

    temp_net_down <- rhdf5::h5read(net_hdf5, "net",
                           index = list(genes_index_down, NULL))
    temp_net_up <- rhdf5::h5read(net_hdf5, "net",
                           index = list(genes_index_up, NULL))

    rownames(temp_net_down) <-  rownames(deg_sub)[deg_sig[["down"]]]
    colnames(temp_net_down) <-  genes
    rownames(temp_net_up)   <-  rownames(deg_sub)[deg_sig[["up"]]]
    colnames(temp_net_up)   <-  genes


    node_degrees[["down"]] <- cbind( rowSums(temp_net_down), colSums(temp_net_down)[genes_index_down] )
    node_degrees[["up"]] <- cbind( rowSums(temp_net_up), colSums(temp_net_up)[genes_index_up] )
    node_degrees[["n_genes_total"]] <- dim(temp_net_down)[2]
    node_degrees[["n_genes_down"]]  <- dim(temp_net_down)[1]
    node_degrees[["n_genes_up"]]    <- dim(temp_net_up)[1]


    sub_net[["down"]] <- temp_net_down[,genes_index_down]
    sub_net[["up"]] <- temp_net_up[,genes_index_up]


    combined_output <- list(deg_sig, fc_sig, sub_net, node_degrees, median_net)
    names(combined_output) <- c("deg_sig", "fc_sig", "sub_net", "node_degrees", "median")
    return(combined_output)
}





#' Get sub-network
#'
#' @param gene_list A list of genes (entrez ID)
#' @param network_type name of pre-built network to use (generic, blood, brain)
#' @param flag_occr boolean to select occurence network or coexpression network
#' @param dir dummy directory for the networks
#' @return \code{combined_output}
#' @examples
#'
#' dir <- '/directflow/SCCGGroupShare/projects/sarba2/data/agg_coexp/'
#' dir <- 'S:/data/agg_coexp/'
#'
#' network_type='generic'
#' subset_network_hdf5_gene_list(gene_list, network_type, dir=dir)
#'
#' @import rhdf5 stats utils
#' @export
#'

subset_network_hdf5_gene_list <- function(gene_list, network_type = "generic",
                                          flag_occr = TRUE, dir = "") {


  if (network_type == "generic" | network_type == "generic230" | network_type == "generic75neg"  | network_type == "blood" | network_type == "brain") {
    if (flag_occr == TRUE) {
      network_type <- paste0(network_type, ".occr")
    }

    genes_hdf5 <- paste0(dir, network_type, ".genes.h5")
    median_hdf5 <- paste0(dir, network_type, ".med.h5")

    genes <- rhdf5::h5read(genes_hdf5, "genes")
    colnames(genes) <- c("entrezID", "name", "ensemblID")
    median_net <- rhdf5::h5read(median_hdf5, "median")

    net_hdf5 <- paste0(dir, network_type, ".net.h5")



    # Check genes (probably should make this a separate function)
    grep_res <- grep ("^ENSG", head (gene_list)  )

    if( length(grep_res ) >0 ) {
      # set genes to ensemblID
      genes <- genes[,3]

    } else {
      # check if entrez IDs
      grep_res <- grep ("^\\d+", head (gene_list) )

      # if true, set to entrez ID
      if( length(grep_res ) >0 ) {
        genes <- genes[,1]
        # if not, set to gene symbols/names
      } else {
        genes <- genes[,2]
      }
    }

  }

  m <- match(gene_list, genes)
  f_m <- !is.na(m)
  f_am <- m[f_m]

  gene_list_match <- gene_list[f_m]

  sub_net <- list()
  node_degrees <- list()
  genes_index <- f_am

  temp_net <- rhdf5::h5read(net_hdf5, "net",
                            index = list(genes_index, NULL))

  rownames(temp_net) <-  gene_list_match
  colnames(temp_net) <-  genes


  node_degrees[["genes"]] <- cbind( rowSums(temp_net ), colSums(temp_net )[genes_index ] )
  node_degrees[["n_genes_total"]] <- dim(temp_net)[2]
  node_degrees[["n_genes"]]  <- dim(temp_net)[1]

  sub_net[["genes"]] <- temp_net[,genes_index]

  combined_output <- list(gene_list_match, sub_net, node_degrees, median_net)
  names(combined_output) <- c("gene_list_match",  "sub_net", "node_degrees", "median")
  return(combined_output)
}

