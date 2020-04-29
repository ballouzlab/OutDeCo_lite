#' Run the recurrence and filtering pipeline
#'
#' @param subgenesets DE genes
#' @param flag up or down regulated genes
#' @param net network to use
#' @param flag_out boolean to save intermediate output files
#' @return \code{clusterids}
#' @examples
#'
#' cpm <- matrix(  abs(rnorm(100)*10) , ncol=1, nrow=100)
#' cpm <- sapply(1:10, function(i) cpm[,1])
#' cpm <- jitter(cpm, amount = 4)
#' network <- matrix(rank(cor(t(cpm))), nrow=100) - 1
#' network <- network/max(network)
#' diag(network) <- 1
#' rownames(network) <-paste0('gene', 1:100 )
#' colnames(network) <-paste0('gene', 1:100 )
#' nettype <- 'random'
#' dir <- 'out_test'
#' cluster_coexp(network)
#'
#' @import EGAD dynamicTreeCut viridis stats utils
#' @importFrom gplots heatmap.2
#' @export
#'
#'
run_filtering <- function( subgenesets, flag, network_type, flag_out=FALSE,
                           filt_min=6, flag_med=TRUE, flag_dist = FALSE){

  res.prec <- list()
  sub_net <- list()
  clust_net <- list()

  n_studies <- dim(subgenesets)[2]

  gene_names <- rownames(subgenesets)
  gene_names_entrez <- EGAD::attr.human$entrezID[match( rownames(subgenesets), EGAD::attr.human$name) ]

  f_e <- !is.na(gene_names_entrez)
  subgenesets <- subgenesets[f_e,]
  gene_names_entrez <- gene_names_entrez[f_e]

  recur <-  rowSums(subgenesets)
  fdrs  <- calc_fdrs_recur( subgenesets )
  gene_list_sub <- gene_names_entrez[recur>0]
  gene_list_sub <- gene_list_sub[!is.na(gene_list_sub)]

  # get coexpression data once, could cause memory problems if gene list is large, needs a warning
  temp_sub_net <- subset_network_hdf5_gene_list(gene_list_sub , network_type, dir=GLOBAL_DIR)


  subgenesets_filt <- subgenesets * 0

  for(i in 1:n_studies ){

    m <- match ( temp_sub_net$gene_list_match , gene_names_entrez[which(subgenesets[,i] == 1 )] )
    f_m <- !is.na(m)
    f_am <- m[f_m]

    sub_net[[i]]  <-  temp_sub_net$sub_net$genes[f_m,f_m]
    n <- dim(sub_net[[i]])[1]
    if( length(n) == 0 ) { print(n); next; } # no genes overlapping with coexpression network


    clust_net[[i]] <- cluster_coexp( sub_net[[i]] , medK = as.numeric(temp_sub_net$median),
                                     flag_med=flag_med, flag_dist = flag_dist )
    clust_size <- plyr::count(clust_net[[i]]$clusters$labels )
    clust_keep <-  clust_size[clust_size[,2] < filt_min ,1]
    genes_keep <- !is.na(match( clust_net[[i]]$clusters$labels, clust_keep))

    tempfilt = rep(0,n)
    tempfilt[genes_keep] <- 1:sum(genes_keep)

    names(tempfilt) <- rownames(sub_net[[i]])
    res.prec[[i]] <- tempfilt

    m <- match(  clust_net[[i]]$clusters$genes[genes_keep] , gene_names_entrez)
    subgenesets_filt[m,i] <- 1

  }

  pre_post <- cbind(colSums(subgenesets), colSums(subgenesets_filt))
  recur_filt <- rowSums(subgenesets_filt)
  fdrs_filt <- calc_fdrs_recur( subgenesets_filt )

  #if(flag_out==T){
   # save(res.prec, file=paste(disease,"res.prec", flag, "Rdata", sep="."))
  #  save(fdrs.filt, subgenesets.filt, file=paste("fdr_calcs", flag, "Rdata", sep="."))
  #  save(fdrs.filt, subgenesets.filt, file=paste("fdr_calcs.filt", flag, "Rdata", sep="."))
  #  save(recur, recur.filt, subgenesets.filt, subgenesets, pre.post, file=paste(disease,"recur.filt.genes", flag, "Rdata", sep="."))
  #}
  combined_output <- list(recur, recur_filt,
                          fdrs, fdrs_filt,
                          subgenesets, subgenesets_filt,
                          pre_post,
                          sub_net, clust_net)

  names(combined_output) <- c("Recurrence", "Recurrence_filtered",
                              "FDRs", "FDRs_filtered",
                              "genesets", "genesets_filtered",
                              "pre_post_totals", "sub_networks", "clustering_results")
  return( combined_output )
}

