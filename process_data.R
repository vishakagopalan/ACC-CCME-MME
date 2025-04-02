load_generic_cancer_data <- function(cancer_type) {
    obj_path <- file.path(base_path,paste(cancer_type,"Seurat_Object.rds",sep="_"))
    if (!file.exists(obj_path)) {
        mat <- Read10X(file.path(base_path,"export",paste(cancer_type,"counts",sep="_")))
        meta_data_df <- fread(file.path(base_path,paste(cancer_type,"metadata.csv.gz",sep="_"))) %>% tibble::column_to_rownames("Cell")
        obj <- CreateSeuratObject(mat,project=cancer_type,meta.data=meta_data_df)
        saveRDS(obj,obj_path)
    }  else {
        obj <- readRDS(obj_path) 
    }

    return(obj)
}

load_colon_cancer_data <- function() {
    return(load_generic_cancer_data("CRC"))
}

TCGAtranslateID = function(file_ids, legacy = FALSE) {
    info = files(legacy = legacy) %>%
        filter( ~ file_id %in% file_ids) %>%
        GenomicDataCommons::select('cases.samples.submitter_id') %>%
        results_all()
    # The mess of code below is to extract TCGA barcodes
    # id_list will contain a list (one item for each file_id)
    # of TCGA barcodes of the form 'TCGA-XX-YYYY-ZZZ'
    id_list = lapply(info$cases,function(a) {
        a[[1]][[1]][[1]]})
    # so we can later expand to a data.frame of the right size
    barcodes_per_file = sapply(id_list,length)
    # And build the data.frame
    return(data.frame(file_id = rep(ids(info),barcodes_per_file),
                      submitter_id = unlist(id_list)))
}

create_TCGA_matrices <- function(project_name,purity_quantiles=c("low_purity"=0.2,"high_purity"=0.8)) {
    tcga_output_file_path <- sprintf("../public_data/%s_Project_Files.rds",project_name)
    purity_df <- fread("../public_data/pancan_tumor_purity.tsv")   
    if (file.exists(tcga_output_file_path)) {
        tcga_var_list <- readRDS(tcga_output_file_path)
        return(tcga_var_list)
    }
   tcga_var_list <- list()
   project_case_ids <- cases() %>% filter(project.project_id == project_name) %>% ids()
   project_clinical_info <- gdc_clinical(project_case_ids)

   project_filter <- files() %>%
        filter( cases.project.project_id == project_name)

    project_adj_norm_manifest_df <- project_filter %>% filter(cases.samples.sample_type %in% c("solid tissue normal") &
              analysis.workflow_type == "STAR - Counts") %>% 
    manifest %>% as.data.frame 

    project_mal_manifest_df <- project_filter %>% filter(cases.samples.sample_type %in% c("primary tumor") &
              analysis.workflow_type == "STAR - Counts") %>% 
    manifest %>% as.data.frame 

    project_adj_norm_translate_df <- TCGAtranslateID(project_adj_norm_manifest_df$id) %>%
    mutate(sample_type="Adj_Norm",purity=NA)
    project_mal_translate_df <- TCGAtranslateID(project_mal_manifest_df$id) %>%
    mutate(sample_type="Tumor",array=gsub("[A-Z]$","",submitter_id)) %>%
    merge(.,dplyr::select(purity_df,array,purity),by="array") %>%
    dplyr::select(-array)
     
    project_translate_df <- rbind(project_adj_norm_translate_df,project_mal_translate_df)

    tcga_var_list[["sample_type"]] <- project_translate_df

    read_count_files <- list("Log2_Counts"="../data/TcgaTargetGtex_gene_expected_count.gz",
                             "Log2_FPKM"="../data/TcgaTargetGtex_rsem_gene_fpkm.gz")

    v23_gtf <- rtracklayer::readGFF("../data/gencode.v23.annotation.gtf.gz") %>%
    dplyr::filter(type=="gene") %>% dplyr::filter(!grepl("pseudogene",gene_type)) %>%
    dplyr::select(sample=gene_id,gene_name)
    tcga_samples <- gsub("[A-Z]$","",project_translate_df$submitter_id)

    normal_list <- list("TCGA-UCEC"="Uterus","TCGA-ESCA"="Esophagus","TCGA-BLCA"="Bladder",
                       "TCGA-KICH"="Kidney","TCGA-KIRC"="Kidney","TCGA-KIRP"="Kidney","TCGA-LIHC"="Liver",
                       "TCGA-STAD"="Stomach","TCGA-COAD"="Colon","TCGA-LUAD"="Lung","TCGA-LUSC"="Lung",
                       "TCGA-PRAD"="Prostate","TCGA-THCA"="Thyroid","TCGA-BRCA"="Breast",
                       "TCGA-GBM"="Brain","TCGA-SKCM"="Skin")

    gtex_sample_info_df <- fread("../data/GTEx_v7_Annotations_SampleAttributesDS.txt") %>%
    mutate(SUBJID=as.character(str_match(SAMPID,".*?-.*?-")) %>% gsub("-$","",.),
           tissue=SMTS)
    
    normal_samples <- gtex_sample_info_df %>% dplyr::filter(grepl(normal_list[[project_name]],tissue,ignore.case = T)) %>% pull(SAMPID)
    if (length(normal_samples) == 0) {
        print(sprintf("Zero normal samples found for %s. Exiting.",project_name))
        return(NULL)
    }
    all_samples <- c(tcga_samples,normal_samples)

    for (file_type in names(read_count_files)) {
        print(file_type)
        flush.console()
        file_path <- read_count_files[[file_type]]
        dt <- fread(file_path, nrows=0)
        cols_to_read <- c("sample")
        sample_names <- colnames(dt)
        for (col in all_samples) {
            cols_to_read <- c(cols_to_read,sample_names[grepl(col,sample_names)])
        }

        gene_exp_dt <- fread(file_path,select=unique(cols_to_read)) %>% merge(v23_gtf,.,by="sample") %>% dplyr::select(-sample)
        sample_names <- setdiff(colnames(gene_exp_dt),"gene_name")
        mat <- gene_exp_dt %>% group_by(gene_name) %>% summarize_at(sample_names,sum) %>% 
        tibble::column_to_rownames("gene_name") %>% as.matrix
        tcga_var_list[[file_type]] <- mat
    }

    sample_names <- setdiff(colnames(tcga_var_list[["Log2_Counts"]]),"gene_name")
    
    purity_quantiles <- c("all_samples"=NA,purity_quantiles)
    tcga_var_list[["DE_Adj_Norm_vs_Tumor"]] <- list()
    tcga_var_list[["DE_Healthy_vs_Tumor"]] <- list()
    for (purity_quantile in names(purity_quantiles)) {
        tcga_adj_norm_anno_df <- dplyr::select(project_translate_df,submitter_id,sample_type) %>%
        dplyr::filter(sample_type == "Adj_Norm") %>%
        mutate(submitter_id=gsub("[A-Z]$","",submitter_id),purity=NA) %>%
        dplyr::filter(submitter_id %in% sample_names) %>% unique
            
        tcga_mal_anno_df <- dplyr::select(project_translate_df,submitter_id,sample_type,purity) %>%
        dplyr::filter(sample_type == "Tumor") %>%
        mutate(submitter_id=gsub("[A-Z]$","",submitter_id)) %>%
        dplyr::filter(submitter_id %in% sample_names) %>% unique

        if (purity_quantile == "low_purity") {
            tcga_mal_anno_df <- dplyr::filter(tcga_mal_anno_df,purity < quantile(purity,purity_quantiles["low_purity"],na.rm=T))
        } else if (purity_quantile == "high_purity") {
            tcga_mal_anno_df <- dplyr::filter(tcga_mal_anno_df,purity > quantile(purity,purity_quantiles["high_purity"],na.rm=T))
        }
    
        tcga_anno_df <- rbind(tcga_adj_norm_anno_df,tcga_mal_anno_df)
        rownames(tcga_anno_df) <- NULL
        gtex_anno_df <- data.frame(submitter_id=intersect(normal_samples,cols_to_read),sample_type="Healthy",purity=NA)
        anno_df <- rbind(tcga_anno_df,gtex_anno_df)

        read_count_mat <- floor(2**tcga_var_list[["Log2_Counts"]] - 1)
        mode(read_count_mat) <- "integer"
        read_count_mat <- na.omit(read_count_mat)
        anno_df <- dplyr::filter(anno_df,submitter_id %in% colnames(read_count_mat),) %>% tibble::column_to_rownames("submitter_id") %>%
        mutate(sample_type=factor(sample_type,levels=c("Tumor","Adj_Norm","Healthy")))
        read_count_mat <- read_count_mat[,rownames(anno_df)]
        print(purity_quantile)
        print(nrow(anno_df))
        flush.console()
        obj <- DESeqDataSetFromMatrix(read_count_mat,colData=anno_df,design=~sample_type)

        dds <- DESeq(obj)
        adj_norm_diff_exp_df <- DESeq2::results(dds,contrast=c("sample_type","Adj_Norm","Tumor")) %>% as.data.frame %>% na.omit %>% arrange(log2FoldChange) %>%
        tibble::rownames_to_column("gene") %>% mutate(comparison="Adj_Norm_vs_Tumor")
        healthy_diff_exp_df <- DESeq2::results(dds,contrast=c("sample_type","Healthy","Tumor")) %>% as.data.frame %>% na.omit %>% arrange(log2FoldChange) %>%
        tibble::rownames_to_column("gene") %>% mutate(comparison="Healthy_vs_Tumor")
         
        tcga_var_list[["DE_Adj_Norm_vs_Tumor"]][[purity_quantile]] <- adj_norm_diff_exp_df
        tcga_var_list[["DE_Healthy_vs_Tumor"]][[purity_quantile]] <- healthy_diff_exp_df
    }

    project_case_ids <- cases() %>% filter(project.project_id == project_name) %>% ids()
    project_clinical_info <- gdc_clinical(project_case_ids)
    tcga_var_list[["purity_quantiles"]] <- purity_quantiles[2:length(purity_quantiles)]
    tcga_var_list[["Clinical"]] <- merge( project_clinical_info$diagnoses %>% 
                                mutate(submitter_id=gsub("_diagnosis","",submitter_id)),
      project_clinical_info$demographic %>% 
                                dplyr::select(case_id,days_to_death,gender,age_at_index),
      by="case_id") 

    saveRDS(tcga_var_list,tcga_output_file_path)
    return(tcga_var_list)
}
