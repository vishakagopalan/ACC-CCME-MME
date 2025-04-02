compute_met_adj_healthy_de_genes <- function(seurat_obj,compartment_name=NULL,meta_data_col="cluster",prefix="",
                                            latent_vars_to_use=c("patient_code","site", "nCount_RNA"),reuse_old_results=T,min_cells=1) {
    
    if (is.null(compartment_name)) {
        print("Specify name for compartment. This is only meant for saving DE results to file.")
        flush.console()
        return(NULL)
    }
    file_name <- paste(prefix,compartment_name,"Markers_Met_Adj_Normal_By",meta_data_col,sep="_") %>% paste(.,"rds",sep=".")
    if (file.exists(file_name) && reuse_old_results) {
        print(paste("Reloading saved DE results stored in",file_name))
        flush.console()
        markers_df_list <- readRDS(file_name)
    } else {
        markers_df_list <- list()
    }

    pairs <-  list(c("Adj_Tissue","Healthy"),c("ACC_Met","Adj_Tissue"),c("ACC_Met","Healthy"))
    DefaultAssay(seurat_obj) <- "RNA"
    seurat_obj <- NormalizeData(seurat_obj)
    for (category in unique(seurat_obj$cluster)) {
        print(paste("Cell type :",category))
        flush.console()
        obj <- subset( seurat_obj, subset = cluster == category)
        min_counts <- dplyr::count(obj@meta.data,cluster,sample_type) %>% arrange(n) %>% head(1) %>% pull(n)
        if (min_counts < min_cells) {
            print(min_counts)
            next
        }

        if (length(Cells(obj)) < 100) 
            next
        
        DefaultAssay(obj) <- "RNA"
        obj <- SetIdent(obj,value="sample_type")
        if (!category %in% names(markers_df_list)) {
            markers_df_list[[category]] <- data.frame()
        }  else {
            print("Already run.Skipping.")
            flush.console()
            next
        }
        for (marker_pair in pairs) {
            pair_key <- paste(marker_pair,collapse="_vs_")
            print(pair_key)
            flush.console()

            min_cells_in_smallest_site <- obj@meta.data %>% dplyr::filter(cluster == category) %>% dplyr::count(site) %>% pull(n) %>% min
            num_sites <- obj@meta.data %>% dplyr::filter(cluster == category) %>% dplyr::count(site) %>% nrow
            min_cells_in_sample_type <- obj@meta.data %>% dplyr::filter(cluster == category & sample_type %in% marker_pair) %>% 
            dplyr::count(sample_type) %>% pull(n) %>% min
            count_df <- obj@meta.data %>% dplyr::filter(cluster == category & sample_type %in% marker_pair) %>% dplyr::count()

            if (nrow(count_df) < 2)
                next

            min_cells_in_sample_type <- min(count_df$n)
            if (min_cells_in_sample_type < min_cells)
                next

            if (min_cells_in_smallest_site == 1 || num_sites == 1) {
                latent_vars <- c("patient_code")
            } else {
                latent_vars <- c("patient_code","site")
            }
            
            temp_df <- FindMarkers( obj, ident.1=marker_pair[1], ident.2=marker_pair[2],
                                  test.use="MAST", latent.var=latent_vars,
                                  assay="RNA") %>% tibble::rownames_to_column("gene") %>%
            mutate(comparison=pair_key,latent_vars=paste(latent_vars,collapse=","))
            markers_df_list[[category]] <- rbind(markers_df_list[[category]],temp_df)
            saveRDS(markers_df_list,file_name)
        }
    }

    return(markers_df_list)
}

run_gsea <- function(markers_df_list,pathway_list,nperm=2000,minSize=20) {
    fgsea_combined_list <- list()
    fgsea_combined_dt <- data.table()

    for (cell_type in names(markers_df_list)) {
        df <- markers_df_list[[cell_type]]
        comparisons <- unique(df$comparison)
        for (comp in comparisons) {
            if (length(names(fgsea_combined_list)) == 0 ) {
                fgsea_combined_list[[comp]] <- data.table()
            }
            logfc_df <- dplyr::filter(df,comparison == comp) %>% arrange(-avg_log2FC)
            vec <- dplyr::select(logfc_df,gene,avg_log2FC) %>% tibble::deframe(.)
            fgsea_dt <- fgseaSimple(pathway_list,vec,nperm=nperm,minSize=minSize) %>%
            #mutate(cell_type=cell_type) %>% dplyr::select(-c(pval,ES,nMoreExtreme,leadingEdge))
            mutate(cell_type=cell_type) %>% dplyr::select(-c(pval,ES,nMoreExtreme))
            fgsea_combined_list[[comp]] <- rbind(fgsea_combined_list[[comp]],fgsea_dt)
            fgsea_combined_dt <- rbind(fgsea_combined_dt,fgsea_dt %>% mutate(comparison=comp))
        }
    }

    return(fgsea_combined_dt)
}

make_gsea_heatmap <- function(fgsea_output,min_comparisons_enriched=3,p_val_thresh=0.1,group_heatmap="clustering",num_row_clusters=3,num_column_clusters=3) {
    fgsea_output %>% 
    mutate(comparison=case_when(comparison == "Adj_Tissue_vs_Healthy" ~ "A_vs_H",
                               comparison == "ACC_Met_vs_Healthy" ~ "M_vs_H",
                               comparison == "ACC_Met_vs_Adj_Tissue" ~ "M_vs_A")) %>%
#      dplyr::filter(comparison == "A_vs_H") %>%
#    dplyr::filter(comparison != "M_vs_H") %>%
#      dplyr::filter(comparison == "M_vs_A") %>%
     mutate(column_name=paste(cell_type,comparison,sep="|")) %>%
    dplyr::filter(padj < p_val_thresh) %>% dplyr::select(pathway,NES,column_name) %>% 
    tidyr::pivot_wider(names_from="column_name",values_from="NES",values_fill=0) %>%
    #mutate(pathway=gsub("^KEGG_|^REACTOME_|^WP_|^GSE|^HALLMARK_","",pathway)) %>%
    tibble::column_to_rownames("pathway") %>% 
    as.matrix -> fgsea_mat

    temp <- str_split( colnames(fgsea_mat), "\\|" )
    anno_df <- data.frame( cell_type=sapply(temp,function(x){x[1]}),
                          comparison=sapply(temp,function(x){x[2]}))
    rownames(anno_df) <- colnames(fgsea_mat)

    colors <- circlize::colorRamp2(c(min(fgsea_mat),0,max(fgsea_mat)),c("blue","white","red"))
    #colors <- circlize::colorRamp2(c(min(fgsea_mat),max(fgsea_mat)),c("green","red"))


    comparisons <- unique(anno_df$comparison)
    num_comparisons <- length(comparisons)
    comparison_colors <- brewer.pal(num_comparisons,"YlGnBu") %>% setNames(.,comparisons)
    comparison_colors <- comparison_colors[1:num_comparisons]
    print(comparison_colors)
    anno_colors  <- list(comparison=comparison_colors)#,site=site_colors)
    num_cell_types <- length(unique(anno_df$cell_type))
    if (num_cell_types <= 12) {
        cell_type_colors <- brewer.pal(num_cell_types,"Paired") %>% setNames(.,unique(anno_df$cell_type))
        anno_colors[["cell_type"]] <- cell_type_colors
    }

    # site_colors <- brewer.pal(2,"RdPu")[1:2] %>% setNames(.,unique(anno_df$site))

    #colnames(fgsea_mat) <- gsub("\\|",":",colnames(fgsea_mat))

    if (group_heatmap == "clustering") {
        top_anno <- HeatmapAnnotation(df=anno_df,gp=gpar(lwd=1),
                                      col=anno_colors)
        ht <- Heatmap( fgsea_mat[rowSums(abs(fgsea_mat) > 0) >= min_comparisons_enriched,], col=colors, row_names_gp=gpar(fontsize=8), name="NES",#legend_label_gp=gpar(fontsize=12),
        #                heatmap_legend_param=list(legend_label_gp=gpar(fontsize=12)),
                #column_km = num_column_clusters,
                row_km=num_row_clusters,
                rect_gp=gpar(lwd=1),
                column_names_gp=gpar(fontsize=12),
                top_annotation=top_anno,
                row_gap=unit(5,"mm"),
                column_gap=unit(5,"mm"),
                column_split=anno_df$comparison
        #        cluster_columns =F
        #         cell_fun = function(j, i,"" x, y, width, height, fill) {
        #         grid.text(sprintf("%.2f", fgsea_mat[i, j]), x, y, gp = gpar(fontsize = 12)) }
               )
    } else {

        a_vs_h_columns <- grep("A_vs_H",colnames(fgsea_mat),value=T)
        m_vs_h_columns <- grep("M_vs_H",colnames(fgsea_mat),value=T)
        m_vs_a_columns <- grep("M_vs_A",colnames(fgsea_mat),value=T)

        vec <- c(rep(1,times=length(a_vs_h_columns)),
                rep(2,times=length(m_vs_h_columns)),
                rep(3,times=length(m_vs_a_columns)))
        names(vec) <- c(a_vs_h_columns,m_vs_h_columns,m_vs_a_columns)
        anno_row_names <- gsub(":","|",names(vec))

        top_anno <- HeatmapAnnotation(df=anno_df[anno_row_names,],gp=gpar(lwd=1),
                                      col=anno_colors)
        fgsea_mat <- fgsea_mat[,c(a_vs_h_columns,m_vs_a_columns,m_vs_h_columns)]
        fgsea_mat <- fgsea_mat[rowSums(abs(fgsea_mat) > 0) >= min_comparisons_enriched,]
        ht <- Heatmap( fgsea_mat, col=colors, row_names_gp=gpar(fontsize=8), name="NES",#legend_label_gp=gpar(fontsize=12),
        #                heatmap_legend_param=list(legend_label_gp=gpar(fontsize=12)),
                      heatmap_legend_param = list(direction = "horizontal"),
        #        column_split = vec,
                cluster_columns=F,
                row_km=num_row_clusters,
                rect_gp=gpar(lwd=1),
                column_names_gp=gpar(fontsize=12),
                top_annotation=top_anno,
                row_gap=unit(5,"mm"),
                row_names_max_width = unit(12, "cm"),
                column_gap=unit(5,"mm"),
                column_split=anno_df$comparison
        #         cell_fun = function(j, i,"" x, y, width, height, fill) {
        #         grid.text(sprintf("%.2f", fgsea_mat[i, j]), x, y, gp = gpar(fontsize = 12)) }
               )
    }
    
    draw(ht, merge_legend = TRUE, heatmap_legend_side = "top", 
    annotation_legend_side = "top")

    #draw(ht)
    return(grid.grabExpr(draw(ht)))
}

make_deg_heatmap <- function(markers_df_list,freq_thresh=0.4,p_val_thresh=0.1,
num_row_clusters=3,num_column_clusters=3) {
    final_markers_df <- data.frame()
    for (category in names(markers_df_list)) {
        final_markers_df <- rbind(final_markers_df,
                                  markers_df_list[[category]] %>% mutate(cell_type=category) %>% filter(comparison != "ACC_Met_vs_Adj_Norm") %>%
        mutate(column_id=paste(cell_type,comparison,sep="-")) )
    }
    
    num_de_cols <- unique(final_markers_df$column_id) %>% length
    freq_df <- final_markers_df %>% dplyr::filter(p_val_adj < p_val_thresh) %>% dplyr::count(gene) %>% mutate(freq=n/num_de_cols)
    genes_to_keep <- dplyr::filter(freq_df,freq >= freq_thresh) %>% pull(gene)

    deg_mat <- final_markers_df %>% dplyr::filter(gene %in% genes_to_keep & p_val_adj < p_val_thresh) %>%
    mutate(column_id=gsub("Adj_Tissue","A",column_id) %>% gsub("Healthy","H",.) %>%
           gsub("ACC_Met","M",.)) %>%
    dplyr::select(gene,avg_log2FC,column_id) %>% 
    tidyr::pivot_wider(names_from="column_id",values_from="avg_log2FC",values_fill=0) %>%
    tibble::column_to_rownames("gene") %>% as.matrix -> deg_mat

    cell_types <- str_split(colnames(deg_mat),"-") %>% sapply(.,function(x){x[1]})
    comparisons <- str_split(colnames(deg_mat),"-") %>% sapply(.,function(x){tail(x,1)}) 
    

    anno_df <- data.frame(cell_type=cell_types,comparison=comparisons,row.names=colnames(deg_mat)) %>% filter(comparison != "")

    colors <- circlize::colorRamp2(c(min(deg_mat),0,max(deg_mat)),c("blue","white","red"))
    num_cell_types <- length(unique(anno_df$cell_type))
    comparisons <- unique(anno_df$comparison)
    num_comparisons <- length(comparisons)
    comparison_colors <- brewer.pal(num_comparisons,"YlGnBu") %>% setNames(.,comparisons)
    comparison_colors <- comparison_colors[1:num_comparisons]
    anno_colors  <- list(comparison=comparison_colors)
    if (num_cell_types <= 12) {
        cell_type_colors <- brewer.pal(num_cell_types,"Paired") %>% setNames(.,unique(anno_df$cell_type))
        anno_colors[["cell_type"]] <- cell_type_colors
    }


    top_anno <- HeatmapAnnotation(df=anno_df,gp=gpar(lwd=1),
                                  col=anno_colors)
    ht <- Heatmap( deg_mat, col=colors, row_names_gp=gpar(fontsize=12), name="Log2FC",
                  #legend_label_gp=gpar(fontsize=12),
    #                heatmap_legend_param=list(legend_label_gp=gpar(fontsize=12)),
            column_km = num_column_clusters,
            row_km=num_row_clusters,
            rect_gp=gpar(lwd=1),
            column_names_gp=gpar(fontsize=12),
            row_gap=unit(5,"mm"),
            column_gap=unit(5,"mm"),
            top_annotation=top_anno,
    #         cell_fun = function(j, i,"" x, y, width, height, fill) {
    #         grid.text(sprintf("%.2f", fgsea_mat[i, j]), x, y, gp = gpar(fontsize = 12)) }
           )

#     ggsave("Endothelial_Type_DEG_Heatmap.svg",grid.grabExpr(draw(ht)),height=20,width=12,limitsize=F)

#     options(repr.plot.height=20,repr.plot.width=12)
    draw(ht)
    return(list("heatmap"=grid.grabExpr(draw(ht),"heatmap_genes"=rownames(deg_mat))))
#     options(repr.plot.height=7,repr.plot.width=7)
}

