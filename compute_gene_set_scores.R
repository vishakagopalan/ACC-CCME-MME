library(data.table)
library(dplyr)
library(AUCell)

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

compute_AUCell_scores <-
function(seurat_obj=NULL,gene_sets=NULL,compute_thresholds=T,threshold_type="Global_k1",rankings_obj=NULL,assay_to_use="RNA",
nCores=3,aucMaxRank=0.05) {
    gene_set_names <- names(gene_sets)
    to_return <- list()
    if (!is.null(rankings_obj)) {
        cells_rankings <- rankings_obj
    } else {
        if ("counts" %in% slotNames(seurat_obj[[assay_to_use]])) {
            mat <- seurat_obj[[assay_to_use]]@counts
        } else {
            mat <- seurat_obj[[assay_to_use]]@layers$counts 
            rownames(mat) <- rownames(seurat_obj)
            colnames(mat) <- Cells(seurat_obj)
        }

        cells_rankings <- AUCell_buildRankings(mat, nCores=nCores, plotStats=F,verbose=F)
        to_return$rankings <- cells_rankings
    }
    cells_AUC <- AUCell_calcAUC(gene_sets,
    cells_rankings,verbose=F,nCores=nCores,aucMaxRank=ceiling(aucMaxRank*nrow(cells_rankings)))
    aucell_scores_mat <- t(getAUC(cells_AUC))
    to_return[["auc_mat"]] =aucell_scores_mat

    if (compute_thresholds == T) {
        cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=F )
        auc_thr <- sapply(cells_assignment, function(x){return(x$aucThr$thresholds[threshold_type,"threshold"])})
        to_return[["thresholds"]] <- auc_thr
    }

    return(to_return)
}

compute_shuffled_gene_set_AUCell_scores <-
function(seurat_obj_,gene_sets,sample_info_column="orig.ident",do_sample_wise=T, num_controls=100,
num_bins=10,assay_to_use="RNA",
q_thresh=1.0,nCores=3,gene_universe=NULL, aucMaxRank=0.05 ) {
    control_sd_df <- data.frame()
    rankings_obj_list <- list()
    control_gene_set_list <- list()

    samples <- unique(seurat_obj_@meta.data[[sample_info_column]])
    for (control in 1:num_controls) {
        print(control)
        flush.console()

        if (do_sample_wise) {
            sd_df <- data.frame()
            for (sample_name in samples) {
                cells <- rownames(seurat_obj_@meta.data %>% dplyr::filter(!!sym(sample_info_column) == sample_name))
                if (!sample_name %in% names(control_gene_set_list)) {
                    control_gene_set_list[[sample_name]] <- find_control_gene_sets(seurat_obj_[["RNA"]]@data[,cells],gene_sets)
                }

                if (control == 1) {
                    auc_output <- compute_AUCell_scores(seurat_obj_[,cells],
                                                        compute_thresholds = F,assay_to_use=assay_to_use,
                                                        control_gene_set_list[[sample_name]],rankings_obj=NULL,
                                                        nCores=nCores,aucMaxRank=aucMaxRank )
                    rankings_obj_list[[sample_name]] <- auc_output$rankings
                } else {
                    auc_output <-
                    compute_AUCell_scores(rankings_obj=rankings_obj_list[[sample_name]],assay_to_use=assay_to_use,
                                                        compute_thresholds = F,gene_sets=control_gene_set_list[[sample_name]],
                                                        nCores=nCores,aucMaxRank=aucMaxRank)
                }

                temp_df <- apply(as.matrix(auc_output$auc_mat[cells,]),2,sd) %>%
                tibble::enframe(.,name="gene_set",value="stdev") %>%
                mutate(shuffle=control,sample=sample_name,mean=apply(as.matrix(auc_output$auc_mat[cells,]),2,
                function(x){return(quantile(x,q_thresh))}))
                sd_df <- rbind(sd_df,temp_df)
            }
        } else {
                control_gene_set_list <-
                find_control_gene_sets(seurat_obj_[["RNA"]]@data,gene_sets,gene_universe=gene_universe)
                if (control == 1) {
                    auc_output <- compute_AUCell_scores(seurat_obj_,
                                                        compute_thresholds = F,
                                                        control_gene_set_list,rankings_obj=NULL,
                                                        nCores=nCores,aucMaxRank=aucMaxRank )
                    rankings_obj <- auc_output$rankings
                } else {
                    auc_output <- compute_AUCell_scores(rankings_obj=rankings_obj,
                                                        compute_thresholds = F,gene_sets=control_gene_set_list,
                                                        nCores=nCores,aucMaxRank=aucMaxRank)
                }
                sd_df <- apply(as.matrix(auc_output$auc_mat),2,sd) %>%
                tibble::enframe(.,name="gene_set",value="stdev") %>% 
                mutate(mean=apply(as.matrix(auc_output$auc_mat),2,function(x){return(quantile(x,q_thresh))}),shuffle=control)
        }
        control_sd_df <- rbind(control_sd_df,sd_df)
    }

    return(control_sd_df)
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

