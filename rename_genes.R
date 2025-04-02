library(data.table)
library(dplyr)
library(stringr)

create_gene_name_update_map <- function() {
    hgnc_alias_dt <- fread("HGNC_Aliases.tsv") %>% dplyr::select(-`Previous name`) %>% unique
    prev_to_curr_gene_map <- list()
    curr_symbol_prev_symbol_df <- dplyr::select(hgnc_alias_dt,`Approved symbol`,prev_sym=`Previous symbol`) %>%
                                               mutate(symbol_type="Previous symbol")
    curr_symbol_alias_symbol_df <- dplyr::select(hgnc_alias_dt,`Approved symbol`,prev_sym=`Alias symbol`) %>%
    mutate(symbol_type="Alias symbol")
    curr_alias_prev_map_df <- rbind(curr_symbol_prev_symbol_df,curr_symbol_alias_symbol_df) %>% unique

    prev_name_curr_name_map_vec <- tibble::deframe( curr_alias_prev_map_df %>% dplyr::filter(prev_sym != "") %>%
                   dplyr::select(prev_sym,`Approved symbol`)) 
    
    return(prev_name_curr_name_map_vec)
}

update_gene_names <- function(gene_names_to_update,reference_gene_names,gene_name_map_vec) {
    already_present_genes <- c()
    new_gene_names <- c()
    renamed_genes <- c()
    for (gene in gene_names_to_update) {
        if (!gene %in% reference_gene_names) {
            new_gene_name <- gene_name_map_vec[gene]
            if (new_gene_name %in% reference_gene_names) {
                new_gene_names <- c(new_gene_names,new_gene_name)
                renamed_genes <- c(renamed_genes,new_gene_name)
            } else {
                new_gene_names <- c(new_gene_names,NA)
            }
        } else {
            new_gene_names <- c(new_gene_names,gene)
            already_present_genes <- c(already_present_genes,gene)
        }
    }
    names(new_gene_names) <- gene_names_to_update
    
    return(new_gene_names)
}

rename_seurat_object_genes <- function(obj,reference_gene_names) {
    old_exp_mat <- obj[["RNA"]]@counts
    gene_name_map_vec <- create_gene_name_update_map()
    updated_obj_gene_names <- update_gene_names( gene_names_to_update=rownames(old_exp_mat), 
                                                  reference_gene_names=reference_gene_names, gene_name_map_vec ) %>%
                                                  na.omit

    updated_gene_name_df <- data.frame(new_name=updated_obj_gene_names,
              old_name=names(updated_obj_gene_names))

    num_mismatched_names <- dplyr::filter(updated_gene_name_df,new_name!=old_name) %>% nrow
    if (num_mismatched_names > 0) {
        duplicated_gene_names <- dplyr::count(updated_gene_name_df,new_name) %>% dplyr::filter(n>1) %>% pull(new_name)
        gene_aliases <- list()
        for (gene in duplicated_gene_names) {
            gene_aliases[[gene]] <- dplyr::filter(updated_gene_name_df,new_name == gene) %>% pull(old_name)
        }

        unique_gene_names <- updated_obj_gene_names[!updated_obj_gene_names %in% duplicated_gene_names]
        new_exp_mat <- old_exp_mat[names(unique_gene_names),]
        rownames(new_exp_mat) <- unique_gene_names
        
        mat <- NULL
        for (gene in duplicated_gene_names) {
            summed_mat <- as.matrix(Matrix::colSums(old_exp_mat[gene_aliases[[gene]],])) %>% t
            rownames(summed_mat) <- gene

            if (is.null(mat)) {
                mat <- summed_mat
            } else {
                mat <- rbind(mat,summed_mat)
            }
        }
        mat <- as(mat,"dgCMatrix")
        new_exp_mat <- rbind(new_exp_mat,mat)
    } else {
        new_exp_mat <- old_exp_mat
    }
    meta_data_df <- obj@meta.data
    new_obj <- CreateSeuratObject( new_exp_mat, meta.data=meta_data_df)

    old_lib_sizes <- colSums(obj[["RNA"]]@counts)
    new_lib_sizes <- colSums(new_obj[["RNA"]]@counts)

    lib_size_df <- data.frame(old_lib_size=old_lib_sizes,new_lib_size=new_lib_sizes,diff=new_lib_sizes-old_lib_sizes)
    
    return(list("obj"=new_obj,"conv_table"=updated_gene_name_df,"lib_size_change_table"=lib_size_df))
}


#
# nLiver_corrected_obj <- rename_seurat_object_genes(nLiver_obj,acc_genes)
# nLung_corrected_obj <- rename_seurat_object_genes(nLung_obj,acc_genes)
