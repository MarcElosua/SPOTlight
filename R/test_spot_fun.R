#' This functions takes in a seurat object and creates different mixtures resembling spots at different proportions
#'
#' @param se_obj seurat object; The seurat object NEEDS to have the variable features calculated
#' @param clust_vr Name of the variable containing the cell clustering
#' @param n number of spots to generate
#' @param verbose Name of the variable containing the cell clustering
#' @return This function returnsa dataframe where each column is the proportion of each spot and a sparse matrix with the synthetic spots with the genenames as rownames and colnames being the spot names
#' @export
#' @examples
#'

test_spot_fun <- function(se_obj, clust_vr, n=1000, verbose=TRUE){

  # Check variables
  if(is(se_obj)!="Seurat") {stop("ERROR: se_obj must be a Seurat object!")}
  if(!is.character(clust_vr)){stop("ERROR: clust_vr must be a character string!")}
  if(!is.numeric(n)){stop("ERROR: n must be an integer!")}
  if(!is.logical(verbose)){stop("ERROR: verbose must be a logical object!")}

  cat("Loading packages...", sep="\n")
  suppressMessages(require(DropletUtils)) # For the downsampling
  suppressMessages(require(dtplyr)) # To use dplyr commands with DT speed
  suppressMessages(require(dplyr)) # To use dplyr commands with DT speed
  suppressMessages(require(tidyr)) # To use dplyr commands with DT speed


  se_obj$seurat_clusters <- droplevels(factor(se_obj@meta.data[,clust_vr]))

  print('Generating synthetic test spots...')
  start_gen <- Sys.time()
  # create progress bar
  pb <- txtProgressBar(min = 0, max = n, style = 3)

  # Defining a vector with the highly variable genes
  # hvg <- VariableFeatures(se_obj)

  # Save count matrix
  count_mtrx <- as.matrix(se_obj@assays$RNA@counts)
  # Remove genes not expressed in at least 0.5% of the cells
  min_cells <- 0.005*ncol(se_obj)
  count_mtrx <- count_mtrx[rowSums(count_mtrx != 0) > min_cells,]


  # Iterator to define the name of the spots generated
  # id <- 1

  ds_spots <- lapply(1:n, function(i){

    # Select between 2 and 10 cells randomly from the count matrix
    cell_pool <- sample(colnames(count_mtrx), sample(x = 2:10, size = 1))

    # Determine the weight each cell will have on the synthetic spot
    # weigh <-runif(length(cell_pool))
    # weigh <- weigh/sum(weigh)

    # We're not going to sum the reads as is bc spots are **enriched**
    # so we'll add up the counts and downsample to the ~depth of a typical spot.

    # Create a name for the spot with the necessary info to deconvolute it
    pos <- which(colnames(count_mtrx) %in% cell_pool)
    tmp_ds <- se_obj@meta.data[pos,] %>% mutate(weight = 1)
    # tmp_ds[,'weight'] <- weigh
    name_simp <- paste('spot_',i,sep='')
    # name_long <- paste(with(tmp_ds, paste(rownames(tmp_ds), 'cluster', seurat_clusters, 'weight', weight, sep = '_')), collapse = '-')

    spot_ds <- tmp_ds %>%
      dplyr::select(seurat_clusters, weight) %>%
      dplyr::mutate(seurat_clusters = paste('clust_',seurat_clusters,sep = '')) %>%
      dplyr::group_by(seurat_clusters) %>%
      dplyr::summarise(sum_weights = sum(weight)) %>%
      dplyr::ungroup() %>%
      tidyr::pivot_wider(names_from = seurat_clusters, values_from = sum_weights) %>%
      dplyr::mutate(name = name_simp)

    # Generate synthetic spot

    ## Here we multiply each vector by its weight
    # weighted <- lapply(1:length(cell_pool), function(ii) expr <- as.integer(round(weigh[[ii]]*count_mtrx[hvg,cell_pool[[ii]]],0)))

    ## Next step we add all the vectors by position
    # syn_spot <- Reduce(`+`,weighted)
    # ret_ds <- data.frame(gene=hvg, tmp=syn_spot)

    ## Here we add up the counts of each cell
    syn_spot <- rowSums(as.matrix(count_mtrx[,cell_pool])); sum(syn_spot)
    names_genes <- names(syn_spot)
    ## Downsample
    ### 25k is a bit above average 20k UMIs observed in spatial transcriptomics data then downsample to 20k
    if(sum(syn_spot) > 25000){
      syn_spot_sparse <- downsampleMatrix(Matrix::Matrix(syn_spot,sparse = T), prop = 20000/sum(syn_spot))
    } else {
      syn_spot_sparse <- Matrix::Matrix(syn_spot,sparse = T)
    }

    rownames(syn_spot_sparse) <- names_genes
    colnames(syn_spot_sparse) <- name_simp
    # ret_ds <- data.frame(syn_spot) %>% mutate(gene=names_hvg)

    # colnames(ret_ds) <- c('gene', i)
    # id <<- id+1

    # update progress bar
    setTxtProgressBar(pb, i)
    # Return the transpose so that documents are on the rows and genes are on the columns
    return(list(syn_spot_sparse,spot_ds))
  })

  # Generate sparse matrix of spots
  # ds_syn_spots <- map(ds_spots, 1) %>%
  # Reduce(function(...) merge(..., by='gene', all.x=TRUE), .) %>%
  # column_to_rownames(var="gene") %>%
  # as.matrix() %>%
  # Matrix(sparse = TRUE) %>%
  # t()

  ds_syn_spots <- map(ds_spots, 1) %>%
    base::Reduce(function(m1,m2) cbind(unlist(m1),unlist(m2)), .)

  # Generate dataframe of spot characteristic
  ds_spots_metadata <- map(ds_spots, 2) %>%
    dplyr::bind_rows() %>%
    dtplyr::lazy_dt(.) %>%  # Convert to lazy data so it can be avaluated as a data.table
    # mutate_all(~replace_na(.,0)) %>%
    data.frame()

  ds_spots_metadata[is.na(ds_spots_metadata)] <- 0

  # change column order so that its progressive
  lev_mod <- gsub("[\\+|\\ ]", ".", levels(se_obj$seurat_clusters))
  all_cn <- c(paste('clust_',lev_mod,sep = ''),'name')
  if( sum(all_cn %in% colnames(ds_spots_metadata)) == (nlevels(se_obj$seurat_clusters)+1) ){
    ds_spots_metadata <- ds_spots_metadata[,all_cn]
  } else {

    # stringr::str_replace(levels(se_obj$seurat_clusters), pattern = '[\\+]', replacement = '.')
    # lev_mod <- stringr::str_replace(levels(se_obj$seurat_clusters), pattern = '[\\+]', replacement = '.')
    # lev_mod <- gsub(" ", ".", lev_mod, fixed = TRUE)
    # lev_mod <- stringr::str_replace(lev_mod, pattern = '\\s', replacement = '.')
    # lev_mod <- paste('clust',lev_mod,sep = '_')

    missing_cols <- all_cn[ which( ! all_cn %in% colnames(ds_spots_metadata) ) ]
    ds_spots_metadata[missing_cols] <- 0
    ds_spots_metadata <- ds_spots_metadata[,all_cn]
  }

  # Close progress bar
  close(pb)

  print(sprintf('Generation of %s test spots took %s mins', n, difftime(Sys.time(), start_gen, units = 'mins')))
  print('output consists of a list with two dataframes, this first one has the weighted count matrix and the second has the metadata for each spot')
  return(list(topic_profiles=ds_syn_spots,cell_composition=ds_spots_metadata))
}
