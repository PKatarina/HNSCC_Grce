library(dplyr)
# library(biomaRt)
library(stringr)
# library(multiMiR)
library(RColorBrewer)
library(pheatmap)


idMap <- function(ensembl.id, genes_file){
  # @description: Mapping from ENSEMBL id to hgnc with GRanges object provided
  # by Anja.
  # @ensembl.id - vector of ensembl 
  # @file - path to the GRanges file
  
  genes_granges <- readRDS(file=genes_file) %>%
    mcols()

  return(genes_granges[ensembl.id,])
}


integrateOmicsData <- function(mrna_res_obj, mir_res_obj, dnam_res_obj, mir_targets_cor, genes_file){
  # @description: This data takes res objects and creates a joined data frame
  # that contains significant mRNA and its respected mir targets and dnam 
  # data (which are not filtered for significance)
  
  # Transforming mrna res object to a dataframe
  mrna_res_obj <- as.data.frame(mrna_res_obj) %>%
    dplyr::filter(padj < 0.05)
  
  # Transforming mir data to dataframe and renaming columns
  mir_res_obj <- 
    tibble::rownames_to_column(as.data.frame(mir_res_obj), var = "mir") %>%
    dplyr::select(mir,
                  mir_mean = "baseMean",
                  mir_FC = "log2FoldChange",
                  mir_pval = "pvalue",
                  mir_padj = "padj")
  
  # Transformin dnam data to df and renaming columns
  dnam_res_obj = unique(dnam_res_obj) %>%
    as.data.frame() %>%
    dplyr::mutate(dnam_condensed_ranges = str_c(seqnames,":",start,"-",end, sep = "")) %>%
    dplyr::select(gene = "ensembl",
                  dnam_quot = "log.quot",
                  dnam_pval = "pvalue",
                  dnam_padj = "padj",
                  dnam_num_sites = "num.sites",
                  dnam_condensed_ranges)

  
  # Integrated (joined) table stored as res
  res <- 
    data.frame(gene=rownames(mrna_res_obj),
               gene_mean = mrna_res_obj$baseMean,
               gene_FC=mrna_res_obj$log2FoldChange,
               gene_pval=mrna_res_obj$pvalue,
               gene_pajd=mrna_res_obj$padj) %>%
    # Adding hgnc symbols
    dplyr::mutate(symbol= idMap(gene,genes_file = genes_file)) %>%
    # Reaaranging columns
    dplyr::select(gene, symbol, everything()) %>%
    # Adding mir correlation data
    dplyr::left_join(dplyr::select(mir_targets_cor, gene, mir, mir_r = r, mir_r_pval = pvalue), by = "gene") %>%
    # Adding mir differential expression data
    dplyr::left_join(mir_res_obj, by="mir") %>%
    # Adding dnam quotient of intesity of DNA metylation in cgi promoters
    dplyr::left_join(dnam_res_obj, by = "gene")
  
  #print(res)
  
  return(res)
}

plot_pheatmap <- function(dat, colVec, colBreaks, annotRows, gaps, prefix){
  f <- paste0(prefix, ".png") # prefix defines plot filename
  
  # various colour schemes for annotation columns
  clPurp <- brewer.pal(4, "Purples")
  clPurp[1] <- "white"
  clGr <- brewer.pal(4, "Greens")
  # clGr[1] <- "white"
  clPink <- brewer.pal(4, "RdPu")
  clPink[1] <- "white"
  clPanel <- brewer.pal(8, "Set2")
  clClus <- brewer.pal(12, "Paired")
  
  # assigning schemes to columns in annotRows 
  annotColors=list(
    Cluster=clClus[1:length(unique(annotRows$Cluster))],
    Cpg_num = brewer.pal(8,"Greys"))
  #Abundance_F=clPink,
  #Abundance_T= clClus,
  #bundance_F=clPink,
  #Abundance= c(Rare=clGr[1], Low=clGr[2], Intermediate=clGr[3], Abundant=clGr[4]),
  
  #    Disease=c(noDisease="white", CD=clPurp[1], `CD,UC`=clPurp[2], UC=clPurp[3]),
  #  names(annotColors$Cluster) <- levels(annotRows$Cluster)
  # names(annotColors$abundance) <- levels(annotRows$abundance)
  #  names(annotColors$Panel) <- levels(annotRows$Panel)
  names(annotColors$Cluster) <- unique(annotRows$Cluster)
  print(gaps)
  
  pheatmap(
    mat=dat,
    color=colVec,
    #breaks=colBreaks,
    #legend_breaks=colBreaks,
    # #annotation_row=tb$kmeans,
    annotation_row=annotRows,
    #annotation_col=annotCols,
    annotation_colors=annotColors,
    show_rownames=F,
    fontsize_row=5,
    gaps_row=gaps,
    cluster_rows=F, 
    cluster_cols=F,
    border_color='white',    
    filename=f)
}

kmeansHeatmap <- function(tb, ann_tb = FALSE, clusters, nstart, GSEA=FALSE){
  # @ description
  
  
  timestamp <- Sys.time()
  timestamp <- format(timestamp, "%Y%m%d_%H%M_")
  heatmap.col_vec2 <- brewer.pal(n=9, name="RdBu")
  
  # Coloring try 2
  heatmap.col_vec2 <-colorRampPalette(brewer.pal(n=9, name="RdBu"))(100)
  
  
  for(c in clusters){
    print(paste0("Processing kmeans with ", c, " clusters."))
    
    # nstart ensures stability of results, otherwise you get a different kmeans     # every time
    k <- kmeans(tb, centers = c, nstart = nstart)
    
    # reorder rows by clustering
    tb <- tb %>%
      dplyr::mutate(kmeans= as.numeric(k$cluster)) %>%
      dplyr::arrange(kmeans)
    
    # specify white space between clusters
    x <- table(as.numeric(tb$kmeans))
    gps <- numeric(length(x))
    for(i in 1:length(x)){
      gps[i:length(x)] <- gps[i:length(x)]+x[i]
    }
    
    
    rg <- max(abs(select(tb, -kmeans)))
    
    # Annotation row beside the heatmap
    annotRows <- data.frame(Cluster = as.numeric(tb$kmeans)) %>%
      dplyr::arrange(Cluster) %>%
      dplyr::mutate(Cluster = as.factor(as.character(Cluster)))
    #print(str(annotRows))
    
    annotRows[is.na(annotRows)]<- 0
    rownames(annotRows) <- rownames(tb)

    
    clClus <- brewer.pal(12, "Paired")
    # assigning schemes to columns in annotRows 
    clColors <- clClus[1:c]
    names(clColors) <- 1:c
    annotColors=list(Cluster=clColors)
        #Cpg_num = brewer.pal(8,"Greys"))



    pheatmap(
      mat=select(tb, -kmeans),
      color=heatmap.col_vec2,
      breaks = seq(-rg, rg, length.out = 100),
      annotation_row=annotRows,
      annotation_colors=annotColors,
      show_rownames=F,
      fontsize_row=5,
      gaps_row=gps,
      cluster_rows=F, 
      cluster_cols=F,
      border_color='white',    
      filename=paste0("clustering/",timestamp,"_", c,
                      "clusters","_",nstart,"nstart",".png"))

  
  }
  
}
