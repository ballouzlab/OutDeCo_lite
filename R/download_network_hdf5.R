#' Download networks from dropbox
#'
#' @param network_type name of pre-built network to use (generic, blood, brain)
#' @param flag_occr boolean to select occurence network or coexpression network
#' @param dir dummy directory for the networks
#' @examples
#' network_type='generic'
#'
#' @import utils
#' @export
#'

download_network_hdf5 <- function(network_type = "generic", flag_occr = TRUE, dir = "") {

    net_files_keys <- c("b0v6405hz5zlmv8",  "qkoenzheon8nafj",  "299y0pnwewv9ee6", "2np3e78gjnvoe10",  "wsqrji519uyh03k", "tayd6axapwt29ck")
    gene_files_keys <- c("8fo67lvq6jemjs4", "bs232ltz50yez7o",  "mi25kj1dtxubzw7",  "s4865kljzg5p8pv", "waqkeem6agg05ve", "fyuq0xkhq4s0ars")
    med_files_keys <- c("ucc8uj6p6gc14hu", "gr9ghxp17pe1gaf", "1qwcvvdjr92o22o", "xkioxsl7989ems6", "uykxeie4qgz6dns", "2xoutlukp4cv29x")

    i <- 0
    if (network_type == "generic" ) { i <-  1 }
    if (network_type == "blood" ) { i <- 3 }
    if (network_type == "brain" ) { i <- 5 }
    if (flag_occr == FALSE) { i <- i + 1 }
    if (flag_occr == TRUE ) { network_type <- paste0(network_type, ".occr") }

    genes_hdf5  <- paste0(network_type, ".genes.h5")
    median_hdf5 <- paste0(network_type, ".med.h5")
    net_hdf5    <- paste0(network_type, ".net.h5")
    url <- "https://www.dropbox.com/s/"


    if( i >  0 ) {
        genes_hdf5_dl  <- paste0(url, gene_files_keys[i], "/", genes_hdf5, "?raw=1")
        median_hdf5_dl <- paste0(url, med_files_keys[i], "/", median_hdf5, "?raw=1")
        net_hdf5_dl    <- paste0(url, net_files_keys[i], "/", net_hdf5, "?raw=1")


        if(!file.exists(genes_hdf5)){
            tryCatch( download.file(genes_hdf5_dl,   destfile=genes_hdf5) )
        }

        if(!file.exists(median_hdf5)){
            tryCatch( download.file(median_hdf5_dl,   destfile=median_hdf5) )
        }

        if(!file.exists(net_hdf5)){
            tryCatch( download.file(net_hdf5_dl,   destfile=net_hdf5) )
        }

    }

}


