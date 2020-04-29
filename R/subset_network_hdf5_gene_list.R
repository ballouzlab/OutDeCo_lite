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
