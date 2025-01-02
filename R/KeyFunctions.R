#' Calculates Spearman dissimilarity and t-SNE from a given dataset
#'
#' This function computes the Spearman dissimilarity matrix from the input dataset,
#' processes it for missing values, and performs t-SNE for dimensionality reduction.
#'
#' @param allptmtable A data frame containing numeric data for post-translational modifications.
#' @return A matrix containing t-SNE coordinates (3D).
#' @export
#'
#' @examples
#' SpearmanDissimilarity(allptmtable)
SpearmanDissimilarity <- function(allptmtable) {
    # Add if statement here to make sure functions are formatted correctly #
    # Ensure allptmtable is a data frame with numeric values #
    allptmtable <- as.data.frame(lapply(allptmtable, as.numeric))
    print("Converting Data Types...")

    # Calculate Spearman correlation #
    allptmtable.cor <- cor(t(allptmtable), use = "pairwise.complete.obs", method = "spearman")
    print("Calculating Spearman Correlation...")

    # Replace diagonal with NA #
    diag(allptmtable.cor) <- NA

    # Calculate dissimilarity #
    dissimilarity.allptmtable <- 1 - abs(allptmtable.cor)
    print("Calculating Spearman Dissimilarity...")

    # Handle any remaining NA values by setting them to the maximum dissimilarity #
    max_dissimilarity <- max(dissimilarity.allptmtable, na.rm = TRUE)
    dissimilarity.allptmtable[is.na(dissimilarity.allptmtable)] <- max_dissimilarity
    print("Filtering missing values...")

    # Make sure the dissimilarity matrix is numeric and suitable for t-SNE #
    dissimilarity.allptmtable <- as.matrix(dissimilarity.allptmtable)

    # Run t-SNE #
    tsne_results <- Rtsne(dissimilarity.allptmtable, dims = 3, perplexity = 15, theta = 0.25, max_iter = 5000, check_duplicates = FALSE, pca = FALSE)
    print("Mapping Data Points...")
    # Return t-SNE results #
    return(tsne_results$Y)
}

#' Calculates Euclidean Distance between
#'
#' @param
#'Post translation modification data frame
#' @return
#' TSNE Euclidean Distance
#' @export
#'
#' @examples
#' (EuclideanDistance(allptmtable.df)
EuclideanDistance <- function(allptmtable.df) {
    # Add if statement here to make sure functions are formatted correctly #
    # Convert the dataframe to a distance matrix using Euclidean distance #
    allptmtable.df.dist = as.matrix(dist(allptmtable.df, method = "euclidean"))
    print("Converting Data Types...")

    # Compute the maximum distance in the matrix, excluding NA values #
    max_dist = max(allptmtable.df.dist, na.rm = TRUE)
    print("Finding maximum distance...")

    # Replace NA values in the distance matrix with 100 times the maximum distance #
    allptmtable.df.dist[is.na(allptmtable.df.dist)] <- 100 * max_dist
    print("Filtering missing values...")

    # Normalize the distance matrix by scaling it to a range from 0 to 100 #
    allptmtable.df.dist.1 <- 100 * allptmtable.df.dist / max_dist
    print("Normalizing distances...")

    # Apply t-SNE to the distance matrix to reduce dimensions to 3 #
    # Parameters: dims = 3 (3D output), perplexity = 15, theta = 0.25 (speed/accuracy trade-off) #
    # max_iter = 5000 (number of iterations), check_duplicates = FALSE (treat rows as unique) #
    # pca = FALSE (no initial PCA) #
    eu.allptms.tsne.list <- Rtsne(as.matrix(allptmtable.df.dist.1), dims = 3, perplexity = 15, theta = 0.25, max_iter = 5000, check_duplicates = FALSE, pca = FALSE)

    # Extract the t-SNE results from the output list #
    eu.allptms.tsne <- eu.allptms.tsne.list$Y
    print("Mapping Data Points...")

    # Return the t-SNE results #
    return(eu.allptms.tsne)
}

#' Finds the difference between Euclidean Distance and Spearman Dissimilarity
#'
#' @param
#' PTM data set & data frame
#' @return
#' TSNE coordinates
#' @export
#'
#' @examples
#' CombinedPar(allptmtable.df, allptmtable)
CombinedPar <- function(allptmtable.df, allptmtable) {
    # Creates a cluster
    cl <- makeCluster(2)  # Uses two cores, may increase later #
    # Using makecluster & not parLapply so that this works with Windows machines as well as Unix based ones #
    registerDoParallel(cl)

    # Export necessary functions and data to each cluster node #
    clusterExport(cl, list("SpearmanDissimilarity", "EuclideanDistance", "allptmtable", "allptmtable.df"))
    clusterEvalQ(cl, {
        library(Rtsne)
        library(parallel)
        library(foreach)
    })

    # Run SpearmanDissimilarity and EuclideanDistance in parallel #
    results <- foreach(i = 1:2, .combine = 'list', .packages = c("Rtsne")) %dopar% {
        if (i == 1) {
            return(SpearmanDissimilarity(allptmtable))
        } else {
            return(EuclideanDistance(allptmtable.df))
        }
    }

    # Extract results #
    spearman_result <- results[[1]]
    euclidean_result <- results[[2]]

    # Continue with the rest of the function #
    combined_distance <- (spearman_result + euclidean_result) / 2

    # Perform t-SNE on the combined distances #
    tsne_result <- Rtsne(as.matrix(combined_distance), dims = 3, perplexity = 15, theta = 0.25, check_duplicates = FALSE, pca = FALSE)
    tsne_coordinates <- tsne_result$Y

    # Stop the cluster #
    stopCluster(cl)
    return(tsne_coordinates)
}

#' Makes a list of cluster groupings
#'
#' @param
#' TSNE Data, Cluster Distance Number, TBD
#' @return
#' TSNE list (filtered)
#' @export
#'
#' @examples
MakeClusterList <- function(tsnedata, toolong, tbl.sc)	{ # Run for all three not just one
    toolong = 3.5
    tsne.span2 <- spantree(dist(tsnedata), toolong=toolong)
    tsnedata.disc2 <-  distconnected(dist(tsnedata), toolong = toolong, trace = TRUE)  # test
    cat ("threshold dissimilarity", toolong, "\n", max(tsnedata.disc2), " groups","\n")
    ordiplot(tsnedata)
    #lines(tsne.span2, tsnedata)
    ordihull(tsnedata, tsnedata.disc2, col="red", lwd=2)
    # Find groups
    tsnedata.span2.df <- data.frame(rownames(tbl.sc))
    names(tsnedata.span2.df) <- "Gene.Name"
    tsnedata.span2.df$group <- tsnedata.disc2
    tsnedata.span2.list <- dlply(tsnedata.span2.df, .(group))  # GROUP LIST  !
    return(tsnedata.span2.list)
}

#' Finds correlations between clusters
#'
#' @param eu_allptms_tsne
#' @param sp_allptms_tsne
#' @param sed_allptms_tsne
#' @param allptmtable_df
#' @param output_dir
#'
#' @return list(eu_allptms_list = eu_allptms_list, sp_allptms_list = sp_allptms_list, sed_allptms_list = sed_allptms_list
#' @export
#'
#' @examples
FindCommonCluster <- function(eu_allptms_tsne, sp_allptms_tsne, sed_allptms_tsne, allptmtable_df, output_dir = "plots") {
    if (!exists("MakeClusterList")) {
        stop("The function 'MakeClusterList' is not defined.")
    }

    # Create output directory if it doesn't exist #
    if (!dir.exists(output_dir)) {
        dir.create(output_dir)
    }

    # Create cluster lists (To be changed) #
    eu_allptms_list <- MakeClusterList(eu_allptms_tsne, 3.8, allptmtable_df)
    sp_allptms_list <- MakeClusterList(sp_allptms_tsne, 3.8, allptmtable_df)  # sp.groups
    sed_allptms_list <- MakeClusterList(sed_allptms_tsne, 3.0, allptmtable_df)  # sed.groups

    # Calculate cluster sizes #
    spsizes_allptms <- sapply(sp_allptms_list, function(x) dim(x)[1])
    sedsizes_allptms <- sapply(sed_allptms_list, function(x) dim(x)[1])
    esizes_allptms <- sapply(eu_allptms_list, function(x) dim(x)[1])

    # Plot and save histograms #
    plot_names <- c("Euclidean_tSNE_Cluster_Sizes.png",
                    "Spearman_tSNE_Cluster_Sizes.png",
                    "Combined_tSNE_Cluster_Sizes.png")

    plot_data <- list(esizes_allptms, spsizes_allptms, sedsizes_allptms)
    plot_colors <- c("yellow", "purple", "brown")
    plot_titles <- c("Euclidean t-SNE Cluster Sizes",
                     "Spearman t-SNE Cluster Sizes",
                     "Combined t-SNE Cluster Sizes")

    for (i in 1:3) {
        png(file.path(output_dir, plot_names[i]), width = 800, height = 600)
        hist(plot_data[[i]], breaks = 100, col = plot_colors[i],
             main = plot_titles[i], xlab = "Cluster Size", ylab = "Frequency")
        dev.off()
        print(paste("Saved plot:", plot_names[i]))
    }

    # Return the cluster lists for further use if needed
    return(list(eu_allptms_list = eu_allptms_list,
                sp_allptms_list = sp_allptms_list,
                sed_allptms_list = sed_allptms_list))
}

# Helper function to find intersections of clusters #
list.common <- function(list1, list2, keeplength = 3) {
  parse <- lapply(list1, function(y) sapply(list2, function(x) intersect(x, y)))
  dims <- lapply(parse, function(x) sapply(x, length))
  keep <- which(sapply(dims, sum) > keeplength)
  pare <- parse[keep]
  prune <- lapply(pare, function(y) return(y[which(sapply(y, function(x) which(length(x) > keeplength)) > 0)]))
  newlist <- unlist(prune, recursive = FALSE)
  return(newlist)
}

# Made into one function instead of two, so that the user has less to manage #
#' Title
#'
#' @param eu.allptms.list
#' @param sp.allptms.list
#' @param sed.allptms.list
#' @param allptmtable.df
#' @param keeplength
#' @param output_dir
#'
#' @return
#' @export
#'
#' @examples
GenerateAndConstructAllptmsNetwork <- function(eu.allptms.list, sp.allptms.list, sed.allptms.list,
                                               allptmtable.df, keeplength = 2, output_dir = "plots") {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  # Mark's Functions #
  "%w/o%" <- function(x, y) x[!x %in% y] #--  x without y
  without <- function(x, y) x[!x %in% y] #--  x without y
  nmissing <- function(x) sum(is.na(x))
  filled <- function (x) {length(x) - nmissing(x)}
  fractNA <- function(df) {
    result <- nmissing(df)/(dim(df)[1]*dim(df)[2])
    return(result)
  }
  mean.na <- function(x) mean(x, na.rm=TRUE)
  max.na <- function(x) max(x, na.rm=TRUE)
  min.na <- function(x) min(x, na.rm=TRUE)
  sd.na <- function(x) sd(x, na.rm=TRUE)
  outersect <- function(x,y){sort(c(setdiff(x,y), setdiff(y,x)))}

  # Convert lists to data frames #
  eu.allptms.df <- ldply(eu.allptms.list)[, 2:3]
  sp.allptms.df <- ldply(sp.allptms.list)[, 2:3]
  sed.allptms.df <- ldply(sed.allptms.list)[, 2:3]

  # Make group names unique #
  eu.allptms.df$group <- paste(eu.allptms.df$group, "e", sep = "")
  sp.allptms.df$group <- paste(sp.allptms.df$group, "s", sep = "")
  sed.allptms.df$group <- paste(sed.allptms.df$group, "sed", sep = "")

  # Group everything together #
  allptmsgroups.df <- rbind(eu.allptms.df, sed.allptms.df, sp.allptms.df)

  # Functions to extract gene names and PTMs #
  extract.genes.from.clist <- function(clusterlist.element) {
    element <- clusterlist.element[1]
    genes <- unique(sapply(as.character(element$Gene.Name), function(x) unlist(strsplit(x, " ", fixed = TRUE))[1]))
    return(genes)
  }

  extract.peps.from.clist <- function(clusterlist.element) {
    element <- clusterlist.element[1]
    return(as.character(element$Gene.Name))
  }

  eu.allptms.genes <- lapply(eu.allptms.list, extract.genes.from.clist)
  sp.allptms.genes <- lapply(sp.allptms.list, extract.genes.from.clist)
  sed.allptms.genes <- lapply(sed.allptms.list, extract.genes.from.clist)

  eu.allptms.peps <- lapply(eu.allptms.list, extract.peps.from.clist)
  sp.allptms.peps <- lapply(sp.allptms.list, extract.peps.from.clist)
  sed.allptms.peps <- lapply(sed.allptms.list, extract.peps.from.clist)


  eu.sp.allptms <- list.common(eu.allptms.peps, sp.allptms.peps, keeplength)
  eu.sp.allptms.sizes <- sapply(eu.sp.allptms, length)
  eu.sp.sed.allptms <- list.common(eu.sp.allptms, sed.allptms.peps, keeplength)
  eu.sp.sed.allptms.sizes <- sapply(eu.sp.sed.allptms, length)

  # Function to generate data frames for heatmaps and evaluations #
  clust.data.from.vec <- function(vec, tbl) {
    if (class(vec) == "list") {
      vec <- unlist(vec)
    }
    at <- tbl[vec, ]
    acol <- names(at[, which(numcolwise(filled)(at) != 0)])
    if (length(acol) == 1) {
      ats <- data.frame(cbind(rownames(at), as.numeric(at[, acol])))
      names(ats) <- c("Gene.Name", acol)
    } else if (length(acol) >= 2) {
      ats <- cbind(rownames(at), at[, acol])
      names(ats)[1] <- "Gene.Name"
    }
    clust.data <- ats
    return(clust.data)
  }

  # Generate data lists for evaluations #
  eu.sp.sed.allptms.data <- list()
  for (i in 1:length(eu.sp.sed.allptms)) {
    if (length(intersect(eu.sp.sed.allptms[[i]], rownames(allptmtable.df))) == 0) next
    at <- allptmtable.df[unlist(eu.sp.sed.allptms[[i]]), ]
    if (dim(at)[1] < 2 | dim(at)[2] < 2) next
    eu.sp.sed.allptms.data[[i]] <- clust.data.from.vec(eu.sp.sed.allptms[[i]], tbl = allptmtable.df)

    # Save the plot
    plot_file <- file.path(output_dir, paste0("plot_", i, ".png"))
    png(plot_file, width = 800, height = 600)
    plot(eu.sp.sed.allptms.data[[i]])
    dev.off()

    print(paste("Saved plot", i, "to", plot_file))
  }

  # Trim datasets #
  alltrimmedsamples <- apply(allptmtable.df, 1, filled)
  allptms.t <- allptmtable.df[which(alltrimmedsamples > 2), ]
  allptmtable.df <- allptms.t

  # Repair bad clusters #
  bad.clusterlist <- list()
  badptms <- unique(outersect(rownames(allptmtable.df), rownames(allptmtable.df)))

  return(list(allptmtable.df = allptmtable.df, eu.sp.sed.allptms.data = eu.sp.sed.allptms.data))
}

#' Creates an adjacency matrix from the given data set
#'
#' @param list.element
#'
#' @return list.el.mat
#' @export
#'
#' @examples
MakeAdjMatrix <- function(list.element) {
  list.el.mat <- matrix(1, nrow = length(list.element), ncol = length(list.element))
  rownames(list.el.mat) <- list.element
  colnames(list.el.mat) <- list.element
  return(list.el.mat)
}

#' Function to bind matrices and return the necessary outputs
#'
#' @param cluster_list
#' @param correlation_matrix
#'
#' @return
#' @export
#'
#' @examples
BindMatrices <- function(cluster_list, correlation_matrix) {
  # Generate the combined adjacency matrix
  adj_matrix <- rbind.fill.matrix(llply(cluster_list, MakeAdjMatrix))
  rownames(adj_matrix) <- colnames(adj_matrix)

  # Order the adjacency matrix by row and column names
  adj_matrix_ordered <- adj_matrix[order(rownames(adj_matrix)), order(colnames(adj_matrix))]

  # Align the correlation matrix with the ordered adjacency matrix
  matched_rows <- intersect(rownames(adj_matrix_ordered), rownames(correlation_matrix))
  matched_cols <- intersect(colnames(adj_matrix_ordered), colnames(correlation_matrix))
  cccn_matrix <- correlation_matrix[matched_rows, matched_cols]

  # Replace NA values in the correlation matrix
  na_indices <- which(is.na(adj_matrix_ordered), arr.ind = TRUE)
  cccn_matrix <- replace(cccn_matrix, na_indices, NA)

  # Remove self-loops by setting diagonal to NA
  if (any(!is.na(diag(cccn_matrix)))) {
    diag(cccn_matrix) <- NA
  }

  # Return the adjacency and CCCN matrices as a list
  return(list(adj_matrix = adj_matrix_ordered, cccn_matrix = cccn_matrix))
}

#' Function to create a correlation network
#'
#' @param bind_result
#'
#' @return
#' @export
#'
#' @examples
CorrelationNetwork <- function(bind_result) {
  adj_matrix <- bind_result$adj_matrix
  cccn_matrix <- bind_result$cccn_matrix

  # Make igraph object, replacing NA with 0
  cccn_matrix0 <- cccn_matrix
  cccn_matrix0[is.na(cccn_matrix0)] <- 0
  graph <- graph_from_adjacency_matrix(as.matrix(cccn_matrix0), mode = "lower", diag = FALSE, weighted = "Weight")

  # Return the graph object
  return(graph)
}

zero.to.NA <- function(df) {
  cf <- df
  zer0 <- which(cf==0, arr.ind = TRUE)
  cfNA <- as.matrix(cf)
  cfNA[zer0] <- NA
  cfNA <- data.frame(cfNA)
  return(cfNA)
}

#' Title
#'
#' @param eu.sp.sed.allptms
#' @param sed.allptms.peps
#' @param AlldataPTMs_cor
#'
#' @return
#' @export
#'
#' @examples
process_ptms_data <- function(eu.sp.sed.allptms, sed.allptms.peps, AlldataPTMs_cor) {
  # Set variables
  eu_sp_sed_allptms <- list.common(eu.sp.sed.allptms, sed.allptms.peps, keeplength = 2)

  # Create adjacency matrices
  allptms_adj <- rbind.fill.matrix(llply(eu_sp_sed_allptms, MakeAdjMatrix))
  rownames(allptms_adj) <- colnames(allptms_adj)

  # Order and align matrices
  allptms_adj_o <- allptms_adj[order(rownames(allptms_adj)), order(colnames(allptms_adj))]

  allptms_cccn_1 <- AlldataPTMs_cor[rownames(AlldataPTMs_cor) %in% rownames(allptms_adj_o), colnames(AlldataPTMs_cor) %in% colnames(allptms_adj_o)]

  # Check matrices
  if(length(setdiff(rownames(allptms_adj), rownames(allptms_cccn_1))) != 0) stop("Mismatch in rownames")
  if(length(intersect(rownames(allptms_adj), rownames(AlldataPTMs_cor))) != nrow(allptms_adj)) stop("Mismatch in intersect rownames")

  # Add correlation as edge values in adjacency matrix
  allptms_cccn <- AlldataPTMs_cor[intersect(rownames(allptms_adj_o), rownames(AlldataPTMs_cor)), intersect(colnames(allptms_adj_o), colnames(AlldataPTMs_cor))]

  # Replace NA values
  allptms_NA <- which(is.na(allptms_adj_o), arr.ind = TRUE)
  allptms_cccn <- replace(allptms_cccn, allptms_NA, NA)
  if (any(!is.na(diag(allptms_cccn)))) diag(allptms_cccn) <- NA

  # Make igraph objects
  allptms_cccn0 <- allptms_cccn
  allptms_cccn0[is.na(allptms_cccn0)] <- 0
  allptms_cccn_g <- graph_from_adjacency_matrix(as.matrix(allptms_cccn0), mode = "lower", diag = FALSE, weighted = "Weight")

  # Gene CCCN construction
  allptms_gene_cccn <- data.frame(allptms_cccn, row.names = rownames(allptms_cccn), check.rows = TRUE, check.names = FALSE, fix.empty.names = FALSE)
  allptms_gene_cccn$Gene_Name <- sapply(rownames(allptms_gene_cccn), function(x) unlist(strsplit(x, " ", fixed = TRUE))[1])

  allptms_gene_cccn[lower.tri(allptms_gene_cccn)] <- NA

  allptms_gene_cccn2 <- ddply(allptms_gene_cccn, .(Gene_Name), numcolwise(function(x) sum(x, na.rm = TRUE)), .progress = "tk")

  rownames(allptms_gene_cccn2) <- allptms_gene_cccn2$Gene_Name
  allptms_gene_cccn2 <- allptms_gene_cccn2[, 2:ncol(allptms_gene_cccn2)]
  allptms_gene_cccn2 <- data.frame(t(allptms_gene_cccn2))
  allptms_gene_cccn2$Gene <- sapply(rownames(allptms_gene_cccn2), function(x) unlist(strsplit(x, " ", fixed = TRUE))[1])

  allptms_gene_cccn3 <- ddply(allptms_gene_cccn2, .(Gene), numcolwise(function(x) sum(x, na.rm = TRUE)), .progress = "tk")

  names(allptms_gene_cccn3)[2:ncol(allptms_gene_cccn3)] <- allptms_gene_cccn3$Gene
  rownames(allptms_gene_cccn3) <- allptms_gene_cccn3$Gene

  allptms_gene_cccn0 <- allptms_gene_cccn3[, 2:ncol(allptms_gene_cccn3)]
  allptms_gene_cccn_na <- zero_to_NA_func(allptms_gene_cccn0)

  allptms_gene_cccn_g <- graph.adjacency(as.matrix(allptms_gene_cccn0), mode = "lower", diag = FALSE, weighted = "Weight")

  allptms_gene_cccn_edges <- data.frame(as_edgelist(allptms_gene_cccn_g))
  names(allptms_gene_cccn_edges) <- c("Gene.1", "Gene.2")
  allptms_gene_cccn_edges$Weight <- edge_attr(allptms_gene_cccn_g)[[1]]
  allptms_gene_cccn_edges$interaction <- "correlation"
  allptms_gene_cccn_edges$interaction[allptms_gene_cccn_edges$Weight <= -0.5] <- "negative correlation"
  allptms_gene_cccn_edges$interaction[allptms_gene_cccn_edges$Weight >= 0.5] <- "positive correlation"

  return(allptms_gene_cccn_edges)
}


get.gene.names.from.peps <- function(pepvec, pepsep="; ") {
  genevec=NULL
  for(i in 1:length(pepvec)) {
    x <- unlist(strsplit(as.character(pepvec[i]), pepsep, fixed=TRUE))
    genes <- unique(sapply(as.character(x),  function (x) unlist(strsplit(x, " ",  fixed=TRUE))[1]))
    genevec <- c(genevec, genes)
  }
  return(genevec)
}

find_ppi_edges <- function(input_dataset, gmfilename, nodenames) {
  # Load PPI edges from other databases
  load("PPIEdges.RData")

  # Initialize the STRING database object
  string_db <- STRINGdb$new(version="12.0", species=9606, score_threshold=0, link_data="detailed", input_directory="")

  # Retrieve the proteins from the STRING database
  string_proteins <- string_db$get_proteins()
  print(dim(string_proteins))

  # Read the dataset that you want to combine with the STRING database
  filter_db <- read.table(input_dataset, header = TRUE, sep = "\t")
  print(colnames(filter_db))

  if (!"experimental" %in% colnames(filter_db)) {
    stop("Column 'experimental' not found in input dataset.")
  }

  # Map the genes to STRING IDs
  mapped_genes <- string_db$map(filter_db, "experimental", removeUnmappedRows = TRUE)
  print(head(mapped_genes))

  # Retrieve the interactions for the mapped genes
  interactions <- string_db$get_interactions(mapped_genes$STRING_id)

  # Convert protein IDs to gene names
  interactions$Gene.1 <- sapply(interactions$from, function(x) string_proteins[match(x, string_proteins$protein_external_id), "preferred_name"])
  interactions$Gene.2 <- sapply(interactions$to, function(x) string_proteins[match(x, string_proteins$protein_external_id), "preferred_name"])

  # Filter interactions based on evidence types
  str.e <- interactions[interactions$experiments > 0, ]
  str.et <- interactions[interactions$experiments_transferred > 0, ]
  str.d <- interactions[interactions$database > 0, ]
  str.dt <- interactions[interactions$database_transferred > 0, ]

  # Combine filtered interactions
  combined_interactions <- unique(rbind(str.e, str.et, str.d, str.dt))

  # Assign edge types
  combined_interactions$edgeType <- "STRINGdb"
  combined_interactions[combined_interactions$database > 0, "edgeType"] <- "database"
  combined_interactions[combined_interactions$database_transferred > 0, "edgeType"] <- "database"
  combined_interactions[combined_interactions$experiments > 0, "edgeType"] <- "experiments"
  combined_interactions[combined_interactions$experiments_transferred > 0, "edgeType"] <- "experiments"

  # Calculate weights
  combined_interactions$Weight <- rowSums(combined_interactions[, c("experiments", "experiments_transferred", "database", "database_transferred")])
  combined_interactions$Weight <- combined_interactions$Weight / 1000

  # Create the final edges dataframe from STRINGdb
  combined_edges <- combined_interactions[, c("Gene.1", "Gene.2", "Weight", "edgeType")]

  # Get GeneMANIA edges
  gm_edges <- get.GM.edgefile(gmfilename, nodenames)

  # Combine STRINGdb and GeneMANIA edges
  final_edges <- rbind(combined_edges, gm_edges)

  return(final_edges)
}

# Function to extract gene names from peptide names
pepgene <- function(peps) {
  unique(sapply(peps, function(x) unlist(strsplit(x, " ", fixed=TRUE))[1]))
}

# Function to extract gene names from peptide edge file
extract.gene.names <- function(peptide.edgefile) {
  peps <- c(peptide.edgefile[,1], peptide.edgefile[,2])
  genes <- unique(sapply(peps, function(x) unlist(strsplit(x, " ", fixed=TRUE))[1]))
  return(genes)
}

# Function to create peptide edges
genepep.edges.3 <- function(nodelist, pepkey=ld.key) {
  nodelist <- unique(nodelist)
  gpedges <- pepkey[pepkey$Gene.Name %in% nodelist, 1:2]
  names(gpedges)[1:2] <- c("Gene.1", "Gene.2")
  gpedges$edgeType <- "peptide"
  gpedges$Weight <- 1
  gpedges$Alt.Weight <- 100
  gpedges$Directed <- FALSE
  return(unique(gpedges))
}

# Function to process correlation edges
process_correlation_edges <- function(cor_matrix, mode="lower") {
  g <- graph_from_adjacency_matrix(as.matrix(cor_matrix), mode=mode, diag=FALSE, weighted="Weight")
  edges <- data.frame(as_edgelist(g))
  edges$Weight <- edge_attr(g)[[1]]
  edges$edgeType <- "correlation"
  edges$edgeType[edges$Weight <= -0.5] <- "negative correlation"
  edges$edgeType[edges$Weight >= 0.5] <- "positive correlation"
  edges <- edges[!is.na(edges$Weight),]
  names(edges)[1:2] <- c("Peptide.1", "Peptide.2")
  edges$Gene.1 <- sapply(edges$Peptide.1, pepgene)
  edges$Gene.2 <- sapply(edges$Peptide.2, pepgene)
  return(edges)
}

# Function to filter dual modifications
filter_dual_modifications <- function(edges, mod1, mod2) {
  dual_mod <- edges[intersect(grep(mod1, edges$Peptide.1), grep(mod2, edges$Peptide.2)), ]
  return(dual_mod)
}

# Function to analyze negative correlations
analyze_negative_correlations <- function(edges) {
  neg <- edges[edges$Weight < 0, ]
  vneg <- neg[abs(neg$Weight) >= 0.5, ]
  vvneg <- neg[abs(neg$Weight) > 0.543, ]

  neg_genes <- unique(neg$Gene.1)
  vneg_genes <- unique(vneg$Gene.1)
  vvneg_genes <- unique(vvneg$Gene.1)

  return(list(neg=neg, vneg=vneg, vvneg=vvneg,
              neg_genes=neg_genes, vneg_genes=vneg_genes, vvneg_genes=vvneg_genes))
}

