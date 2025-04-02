
find_control_gene_sets <- function(gene_exp_mat,gene_sets,num_bins=10,gene_universe=NULL) {
    mean_gene_exp_vec <- rowMeans(gene_exp_mat)
    genes_by_bin <- list()
    if (is.null(gene_universe)) {
        gene_universe <- rownames(gene_exp_mat)#rownames(gene_exp_mat)
    }
    gene_exp_df <- tibble::enframe(mean_gene_exp_vec,name="gene",value="gene_exp") %>% 
    mutate(bin=ntile(gene_exp,num_bins)) %>%
    dplyr::filter(gene %in% gene_universe)
    for (exp_bin in 1:num_bins) {
        genes_by_bin[[exp_bin]] <- gene_exp_df %>% dplyr::filter(bin == exp_bin) %>% pull(gene)
        #genes_by_bin[[exp_bin]] <- genes_by_bin[[exp_bin]]
    }

    control_gene_sets <- list()
    for (gene_set_name in names(gene_sets)) {
        gene_set <- gene_sets[[gene_set_name]]
        bin_dist_df <- gene_exp_df %>% dplyr::filter(gene %in% gene_set) %>% group_by(bin) %>% dplyr::count()
        control_genes <- c()
        for (idx in 1:nrow(bin_dist_df)) {
            bin_num <- bin_dist_df[idx,] %>% pull(bin)
            num_genes <- bin_dist_df[idx,] %>% pull(n)
            control_genes <- c(control_genes,sample(genes_by_bin[[bin_num]],size=num_genes))
        }
        control_gene_sets[[paste(gene_set_name,"Control",sep="_")]] <- control_genes
    }
    
    return(control_gene_sets) 
}

compute_bulk_raw_gene_set_scores <- function(gene_exp_mat,gene_sets) {
    gene_set_names <- names(gene_sets)
    gene_set_score_mat <- matrix(nrow=length(gene_set_names),ncol=ncol(gene_exp_mat), dimnames=list(gene_set_names,colnames(gene_exp_mat)))

    for (gene_set_name in gene_set_names)  {
        genes <- gene_sets[[gene_set_name]]
        genes <- genes[genes %in% rownames(gene_exp_mat)]
        gene_set_score_mat[gene_set_name,] <- colMeans(gene_exp_mat[genes,])
    }

    return(gene_set_score_mat)
}

compute_bulk_normalized_gene_set_scores <-
function(gene_exp_mat,gene_sets,num_controls=100,num_bins=10,q_thresh=0.95,gene_universe=NULL,
shuffle_type="gene sets") {
    gene_set_names <- names(gene_sets)
    control_gene_set_score_mat <- matrix(nrow=length(gene_set_names),
    ncol=ncol(gene_exp_mat), dimnames=list(gene_set_names,colnames(gene_exp_mat)))
    foreground_score_mat <- compute_bulk_raw_gene_set_scores(gene_exp_mat,gene_sets)


    background_score_df <- data.frame()
    if (shuffle_type == "gene sets") {
        for (control in 1:num_controls) {
            print(control)
            flush.console()
            control_gene_sets <- find_control_gene_sets(gene_exp_mat,gene_sets,num_bins=num_bins,gene_universe=gene_universe)
            background_score_mat <- compute_bulk_raw_gene_set_scores(gene_exp_mat,control_gene_sets)
            control_score_df <- reshape2::melt( background_score_mat ) %>% mutate(control=control) %>%
            dplyr::rename(gene_set=Var1,sample_name=Var2,score=value)
            background_score_df <- rbind(background_score_df,control_score_df)
        }
        background_score_df <- group_by(background_score_df,sample_name,gene_set) %>% dplyr::summarize(score_thresh=quantile(score,q_thresh))
        background_score_mat <- tidyr::pivot_wider(background_score_df,names_from="sample_name",values_from="score_thresh") %>%
        tibble::column_to_rownames("gene_set") %>% as.matrix
    } else if (shuffle_type == "gene exp") {
        all_genes <- intersect( unique(unlist(gene_sets)), rownames(gene_exp_mat) )
        reduced_exp_mat <- gene_exp_mat[all_genes,]
        for (control in 1:num_controls) {
            print(control)
            flush.console()
            shuffled_mat <- apply(reduced_exp_mat,1,sample)
            rownames(shuffled_mat) <- colnames(reduced_exp_mat)
            shuffled_mat <- t(shuffled_mat)
            background_score_mat <- compute_bulk_raw_gene_set_scores(shuffled_mat,gene_sets)
            control_score_df <- reshape2::melt( background_score_mat ) %>% mutate(control=control) %>%
            dplyr::rename(gene_set=Var1,sample_name=Var2,score=value)
            background_score_df <- rbind(background_score_df,control_score_df)
        }

    }

    return(list("fg"=foreground_score_mat,"bg"=background_score_mat))
}
