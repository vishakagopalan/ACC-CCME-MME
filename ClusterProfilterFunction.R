# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .R
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.2
#   kernelspec:
#     display_name: R/4.2
#     language: R
#     name: ir42
# ---

# +
gene.df_function <- function(markers_filtered){
    
gene <- markers_filtered[, 8] 
gene
#gene <- gene$gene

gene.df <- bitr(gene, fromType = "ALIAS",
        toType = c("ENSEMBL","ENTREZID","SYMBOL"),
        OrgDb = org.Hs.eg.db)
dim(gene.df)


gene.df
}

# +
gene.df_function_I <- function(markers_filtered){
gene <- markers_filtered[, 1]
#gene

gene.df <- bitr(gene, fromType = "ALIAS",
        toType = c("ENSEMBL","ENTREZID","SYMBOL"),
        OrgDb = org.Hs.eg.db)
dim(gene.df)
gene.df
}


# +

geneID_converter <- function(markers_filtered){
gene <- markers_filtered[, 8] 
gene
#gene <- gene$gene

gene.df <- bitr(gene, fromType = "ALIAS",
        toType = c("ENSEMBL","ENTREZID","SYMBOL"),
        OrgDb = org.Hs.eg.db)
dim(gene.df)


gene.df


colnames(gene.df) <- c('gene','ENSEMBL','ENTREZID', "SYMBOL")
head(gene.df)
markers_filtered <- merge( markers_filtered,gene.df,by="gene")

markers_filtered <- markers_filtered %>% 
distinct(gene, .keep_all = TRUE) 

head(markers_filtered %>% arrange(desc(avg_log2FC)))


#d <- read.csv(your_csv_file)
## assume that 1st column is ID
## 2nd column is fold change

## feature 1: numeric vector
GeneList <- markers_filtered[,4]
#GeneList

# ## feature 2: named vector

 names(GeneList) <- as.character(markers_filtered[,10])

# # ## feature 3: decreasing order
 GeneList <- sort(GeneList, decreasing = TRUE)



 GeneList

    

}



# +
geneID_converter_I <- function(markers_filtered){
gene <- markers_filtered[, 1]
#gene

gene.df <- bitr(gene, fromType = "ALIAS",
        toType = c("ENSEMBL","ENTREZID","SYMBOL"),
        OrgDb = org.Hs.eg.db)
dim(gene.df)
#gene.df

colnames(gene.df) <- c('gene','ENSEMBL','ENTREZID', "SYMBOL")
#head(gene.df)
markers_filtered <- merge( markers_filtered,gene.df,by="gene")

markers_filtered <- markers_filtered %>% 
distinct(gene, .keep_all = TRUE) 
#head(markers_filtered %>% arrange(desc(avg_log2FC)), 20)

## feature 1: numeric vector
GeneList <- markers_filtered[,3]

## feature 2: named vector
names(GeneList) <- as.character(markers_filtered [,10])

## feature 3: decreasing order
GeneList <- sort(GeneList, decreasing = TRUE)

GeneList
}
# -



# +
## function for gene name conversion
gene_list <- function(df){
    
gene <- rownames(df) ## row names of each of the list elements
    
gene.df <- bitr(gene, fromType = "ALIAS",
        toType = c("ENSEMBL","ENTREZID","SYMBOL"),
        OrgDb = org.Hs.eg.db)
    
#dim(gene.df)
#head(gene.df)

colnames(gene.df) <- c('gene','ENSEMBL','ENTREZID', "SYMBOL")
head(gene.df)
    
gene.df <- gene.df %>% 
distinct(gene, .keep_all = TRUE) 
gene.df$group <- colnames(df)
gene.df <- gene.df %>% dplyr::select(-c("ENSEMBL", "SYMBOL"))
head(gene.df)
dim(gene.df)
gene.df

} 

# +
back_gene_func <- function(object){
    
    genes.to.keep <- Matrix::rowSums(object@assays$RNA@counts > 0) >= floor(0.1 * ncol(object@assays$RNA@counts))
counts.sub <- object@assays$RNA@counts[genes.to.keep,]

length(rownames(counts.sub))
background_genes <- as.data.frame(rownames(counts.sub))
head(background_genes)

rownames(background_genes) <- background_genes$`rownames(counts.sub)`
head(background_genes)

mydf_background <- gene_list(background_genes) ### use gene_list function to convert gene symbol to ENtrez
head(mydf_background)
mydf_background
    
}


