library(dplyr)
# library(biomaRt)
library(stringr)
# library(multiMiR)
library(RColorBrewer)
library(pheatmap)
library(gprofiler2)

GSEA_SOURCES <- c("GO:BP","KEGG","REAC", "MIRNA", "TF")


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

  names(dnam_res_obj) <-  1:length(dnam_res_obj)
  # dnam_res_obj = unique(dnam_res_obj) %>%
  dnam_res_obj <- dnam_res_obj %>%
    as.data.frame() %>%
    dplyr::mutate(dnam_condensed_ranges = stringr::str_c(seqnames,":",start,"-",end, sep = "")) %>%
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

# Do not use
plot_pheatmap_ <- function(dat, colVec, colBreaks, annotRows, gaps, prefix){
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

calcGaps <- function(vec){
  # @description: returns a vector of gaps separating clusters for plotting 
  # cluster heatmap
  
  # specify white space between clusters
  # using factor to keep the order of clusters from vec
  x <- table(factor(vec, levels=unique(vec)))
  gps <- numeric(length(x))
  for(i in 1:length(x)){
    gps[i:length(x)] <- gps[i:length(x)]+x[i]
  }

  return(gps)
}

plotKmeansPheatmap <- function(tb, additional_annotation = NULL, custom_cluster_names = NULL, file=NA, breaks = NULL ){
  
  # Number of clusters
  c <- length(unique(tb$kmeans))
  
  # Colors
  heat_colors <- colorRampPalette(brewer.pal(n=9, name="RdBu"))(7)
  clusters_annot_colors <- brewer.pal(12, "Paired")[1:c]
  # Named vectors for clusters
  if(is.null(custom_cluster_names)){
    names(clusters_annot_colors) <- 1:c
  }

  
  # For continuous colors (passed in breaks in pheatmap)
  rg <- max(abs(dplyr::select(tb, -kmeans)))
  
  # Specify white space between clusters
  gps <- calcGaps(tb$kmeans)

  # Annotation row beside the heatmap
  annotRows <- data.frame(row.names = rownames(tb),
                          Cluster = as.numeric(tb$kmeans)) %>%
    #dplyr::arrange(Cluster) %>%
    dplyr::mutate(Cluster = as.factor(Cluster))
  
  # Saving colors scheme to list (has to have the SAME name as in annotRows)
  annotColors=list(Cluster=clusters_annot_colors)
  #rownames(tb) <- 1:nrow(tb)

  # Additional annotations
  if(!is.null(additional_annotation)){pass}


  # Main plotting function
  pheatmap(
    mat=dplyr::select(tb, -kmeans),
    color=heat_colors,
    #breaks = seq(-rg, rg, length.out = 100),
    annotation_row=annotRows,
    annotation_colors=annotColors,
    breaks = breaks,
    show_rownames=F,
    fontsize_row=5,
    gaps_row=gps,
    cluster_rows=F, 
    cluster_cols=F,
    border_color='white',    
    filename=file)

}

gseaOnKmeans <- function(tb, custom_sources=NULL){
  # @description: gsea on kmeans table and saving it in xlsx format
  
  if(is.null(custom_sources)){
    GSEA_SOURCES <- c("GO:BP","KEGG","REAC", "MIRNA", "TF")
  }
  else{
    GSEA_SOURCES <- custom_sources
  }

  # Preparing list 
  query <- list()
  for (c in unique(tb$kmeans)) {
    query[[c]] <- rownames(dplyr::filter(tb, kmeans == c))  
  }
  
  print(paste0("GSEA on kmeans table with ", c, " clusters."))

  gost_k <-gost(query = query,
               organism= "hsapiens",
               exclude_iea=TRUE,
               domain_scope="annotated",
               ordered_query=F,
               sources=GSEA_SOURCES)
  
  return(gost_k)
  
}

kmeansHeatmap <- function(tb, ann_tb=FALSE, clusters, nstart, GSEA=FALSE, scaled, directory, plot_cost=FALSE){
  # @ description
  
  timestamp <- Sys.time()
  timestamp <- format(timestamp, "%Y%m%d_%H%M_")
  
  
  # COlors for clusters
  #kmeans_color <- colorRampPalette(brewer.pal(n=9, name="RdBu"))(7)
  
  # Scaling as part of kmeans func
  if(scaled == "dnam2") {
    tb$dnam_quot <- 2*(tb$dnam_quot)
    # kmeans color breaks for scaled
    pheatmap_breaks <- c(-6,-4,-2,-0.000001,0.000001,2,4,6)
  }
  if(scaled == "scaled"){
    tb$dnam_quot <- scale(tb$dnam_quot)
    tb$gene_FC <- scale(tb$gene_FC)
    tb$mir_FC <- scale(tb$mir_FC)
    # kmeans color breaks for scaled
    pheatmap_breaks <- c(-2.25,-1.5,-0.75,-0.000001,0.000001,0.75,1.5,2.25)

  }
  else{
    pheatmap_breaks <- c(-6,-4,-2,-0.000001,0.000001,2,4,6)
  }
  # setting NAs to 0
  tb[is.na(tb)] <- 0
  

  cost_vec <- vector(mode = "numeric")

  
  for(c in clusters){
    print(paste0("Processing kmeans with ", c, " clusters."))
    
    # nstart ensures stability of results, otherwise you get a different kmeans
    # every time
    k <- kmeans(tb, centers = c, nstart = nstart, iter.max=10)
    
    
    # reorder rows by clustering
    tb <- tb %>%
      as.data.frame() %>%
      dplyr::mutate(kmeans= as.numeric(k$cluster)) %>%
      # Saving the rownames to column
      dplyr::mutate(gene = rownames(tb))
    

    
    
    # ordering by mean values
    tb_order <- tb %>%
      group_by(kmeans) %>%
      dplyr::summarise(avg_kmeans = mean(gene_FC)) %>%
      dplyr::arrange(desc(avg_kmeans)) %>%
      dplyr::mutate(k = 1:c)
    
    tb <- left_join(tb_order, tb, by="kmeans") %>%
      dplyr::select(gene, gene_FC, mir_FC, dnam_quot,  -avg_kmeans, kmeans = k) %>%
      as.data.frame()
    
    # Explicitly adding back rownames and dropping the gene column
    rownames(tb) <- tb$gene
    tb <- dplyr::select(tb, dplyr::everything(), -gene)

    # Timestamped filename
    filename = 
      paste0(directory ,timestamp,"_", c, "clusters","_",nstart,"nstart","_", scaled)
    
    
    # Passing kmeans results to the plotting function
    plotKmeansPheatmap(tb, file = paste0(filename,".png"), breaks = pheatmap_breaks)
    
    # Saving tables for later use
    saveRDS(tb, file = paste0(filename,".rds"))
    
    if(GSEA){
      write.csv(dplyr::select(gseaOnKmeans(tb)$result, -parents), file = paste0(filename,".csv"))
    }
    if(plot_cost){
      cost_vec <- append(cost_vec, k$tot.withinss)
    }
  
  }
  
  if(plot_cost){
    cost_data <- data.frame(cl = clusters, cost=cost_vec)
    print(cost_data)
    
    write.table(cost_data, file = paste0(filename, "_cost_data.csv"), row.names = FALSE)
    
    ggplot(cost_data, aes(x=cl, y=cost))+
      geom_line()+
      geom_point(color="blue")+
      theme_bw()+
      theme(panel.background = element_rect(fill = "white"))
    
    ggsave(filename = paste0(filename, "_cost_plot.png"), )
    
  }
  
  
}


barPlotData <- function(tb, cluster, unscaled_data = NULL){
  sub_tb <- dplyr::filter(tb, kmeans == cluster)
  
  if(!is.null(unscaled_data)){
    unscaled_data[is.na(unscaled_data)] <- 0
    sub_tb <- unscaled_data[rownames(sub_tb),]
  }
    
  output <- data.frame(
    groups = c("gene_FC", "mir_FC", "dnam_quot"),
    mean = c(mean(sub_tb$gene_FC), mean(sub_tb$mir_FC), mean(sub_tb$dnam_quot)),
    sd = c(sd(sub_tb$gene_FC), sd(sub_tb$mir_FC), sd(sub_tb$dnam_quot))
  )
  
  output$groups <- factor(output$groups, levels = output$groups)
  return(output)
}
