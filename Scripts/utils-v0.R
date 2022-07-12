library(dplyr)
library(biomaRt)
library(stringr)
library(multiMiR)


idMap <- function(ensembl.id, genes_file){
  # @description: Mapping from ENSEMBL id to hgnc with
  # GRanges object provided by Anja.
  # @ensembl.id - vector of ensembl 
  # @file - path to the GRanges file
  
  genes_granges <- readRDS(file=genes_file) %>%
    mcols()

  return(genes_granges[ensembl.id,])
}

#idMap("ENSG00000082898", genes_file = "01-Exploratory_analysis_and_differential_expression/data/genes.hg19.ensembl.rds")

#this function creates the IDTABLES FOR idMapper
createIdTable <- function(assembly = "hg19"){
  att <- c("ensembl_gene_id","hgnc_symbol","mirbase_id",
           "uniprot_gn_symbol", "description")
  
  if(assembly == "hg19"){
    ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                      host="grch37.ensembl.org",
                      path="/biomart/martservice", 
                      dataset="hsapiens_gene_ensembl")
  }
   
  else if(assembly == "hg38"){
    ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                       dataset = "hsapiens_gene_ensembl")
  }
  else{
    stop("Argument assembly is invalid, use hg19 or hg38")
  }
  filename <- paste0("IDTABLE_",assembly)
  IDTABLE <- getBM(attributes = att, mart = ensembl)
  saveRDS(IDTABLE, file = filename)
  print(paste0("IDTABLE saved as filename: ", filename))
  
  
}


# IDTABLE_hg38 <- readRDS("IDTABLE")
# IDTABLE_hg19 <- readRDS("IDTABLE_hg19")

idMapper <- function(input.list,
                     input.id = c("ensembl_gene_id"), 
                     output.id = "hgnc_symbol",
                     assembly = "hg38"){
  
  # if(!file.exists(file.path(getwd(), paste0("IDTABLE_", assembly))))
  
  if(assembly == "hg38"){
    IDTABLE <- IDTABLE_hg38
  }
  else if(assembly == "hg19"){
    IDTABLE <- IDTABLE_hg19
  }
  else{
    stop("Argument assembly is invalid, use hg19 or hg38")
    }
  
  res <- NULL
  for(id in input.list){
    #print(id)
    #print(IDTABLE[IDTABLE[,input.id] == id, output.id])
    if(id %in% IDTABLE[,input.id]){
      if(!is.na(IDTABLE[IDTABLE[,input.id] == id, output.id]) && 
         IDTABLE[IDTABLE[,input.id] == id, output.id] != ""){
          res <- append(res,IDTABLE[IDTABLE[,input.id] ==id,output.id][1])
        }
      else{
        res <- append(res,id)
        }
      }
  

    else{
      res <- append(res,id)
    }
  
  }  
  #print(res)
  
  return(res)
  
}

idMapper0.2 <- function(input.list,
                     input.id = "ensembl_gene_id", 
                     output.id = "hgnc_symbol",
                     assembly = "hg19"){
  
  if(assembly == "hg38"){
    IDTABLE <- unique(dplyr::select(IDTABLE_hg38, input.id,output.id))
  }
  else if(assembly == "hg19"){
    IDTABLE <- unique(dplyr::select(IDTABLE_hg19, input.id,output.id))
  }
  else{
    stop("Argument assembly is invalid, use hg19 or hg38")
  }
  IDTABLE <- IDTABLE[!duplicated(dplyr::select(IDTABLE, input.id)),]
  
  res <- data.frame(input.list,NA)
  res <- left_join(res,IDTABLE, by=c("input.list"=input.id))
  #res[res$output.id =="",output.id] <- filter(res, output.id == "")$input.list
  
  #print(unique(dplyr::select(res,input.list,output.id)))
  
    
  return(res[,output.id])
}

mapTester <- function(){

dds <- readRDS("RDS/dds.mRNA")
map.test <- rownames(dds)[sample(1:length(rownames(dds)), 10, replace = T)]
#print(map.test)

map.test.res <- idMapper0.2(map.test, input.id="ensembl_gene_id", assembly = "hg19")
#print(map.test.res)

}

#mapTester()

#map.test <- res.export$row
#map.test <- idMapper(input.id = "ensembl_gene_id",input.list = map.test)
#map.test2 <- idMapper(input.id = "hgnc_symbol",output.id = "ensembl_gene_id",input.list = map.test)

result_export <- function(dds.object,contrast.vector,
                          p.value = 0.05, change.rownames = FALSE){
  res.export <- results(dds.object,tidy = T, contrast = contrast.vector)
  
  #filtering out p- value less than 0.05
  res.export <- filter (res.export, padj != is.na(padj))
  res.export <- res.export[res.export$padj < p.value,]
  
  if(change.rownames){
    res.export$row <- idMapper(input.id = "ensembl_gene_id",
                               input.list = res.export$row)
    print()
  }
  print(res.export$row)
  
  write_xlsx(res.export, path = str_c("Results/3_result_renamed_", contrast.vector[1],"_",contrast.vector[2],"_vs_",contrast.vector[3],".xlsx"))
  
  return(length(rownames(res.export)))
}



filter.mRNA <- function(dds.object, counts = 10){
  #' mRNA filter function
  #' Give a mRNA dds object and the min row count sum.
  
  
  keep <- rowSums(counts(dds.object))> counts
  dds.object <- dds.object[keep,]
  return(dds.object) 
}


filter.microRNA <- function(dds.object, filter.isomirs = T, row.sum = 10){
  #' microRNA filter function
  #' Give a microRNA dds object and boolean isomir filter option. Filters out rows which count sum is less than 10 and isomirs.
  
  if(filter.isomirs){
    numbers <- c("1","2","3","4","5","6","7","8","9","0")
    keep <- rowSums(counts(dds.object))> row.sum & 
      is.element(str_sub(rownames(dds.object),1,1), numbers)==FALSE
  }
  else{
    keep <- rowSums(counts(dds.object))> row.counts 
  }
  dds.object <- dds.object[keep,]
  return(dds.object)  
  
}

filterPvalue <- function(dds.res, padj = T, p.treshold =0.05){
  if(padj){
    dds.res <- dds.res[!is.na(dds.res$padj) & dds.res$padj < 0.05,]
  }
  else{
    dds.res <- dds.res[!is.na(dds.res$pval) & dds.res$pval < 0.05,]
    
  }
  
  return (dds.res)
}

countsPlot <- function(dds.object, regulated = c("up","down"),
                       gene.name.type=c("miRbase","ensembl"), n=9){
  res <- results(dds.object, contrast = c("Hist","cancer","normal"))
  #filtering only significant res
  res <- filterPvalue(res)
  # ordering according to significance
  res <- res[order(res$padj),]
  #print(res)
  if(regulated == "up"){
    res <-res[res$log2FoldChange > 0,]
  }
  else{
    res <-res[res$log2FoldChange < 0,]
  }
  
  gene.symbol <- idMapper(rownames(res[1:n,]), assembly = "hg19")
  
  plot_list <- list()
  
  for(i in (1:n)){
    
    gene <- rownames(res[i,])
    plot.counts <- plotCounts(dds.object, gene = gene, returnData = T,
                              intgroup = "Hist")

    ens.title <- paste0(gene.symbol[i],sep= ", p: ",signif(res[gene,"padj"],3))
    
    p <- ggplot(plot.counts, aes(x=Hist, y=count))+
      {if(max(plot.counts$count > 1000))scale_y_continuous(trans = "log10", labels = scientific)}+
      geom_jitter(width=0.1)+
      {if(gene.name.type=="miRbase")ggtitle(paste0(str_split(gene,"R", simplify = T)[2], sep = ", p: ",signif(res[gene,"padj"],3)))}+
      {if(gene.name.type=="ensembl")ggtitle(ens.title)}+
      theme_classic(base_size = 10)+
      {if(gene.name.type=="ensembl" & str_length(ens.title)>25)theme(plot.title = element_text(size=8))}
    
    
    plot_list[[i]]<- p
  }
  gridExtra::grid.arrange(grobs = plot_list)
  
  
}

filterMRNA <- function(dds.object, row.counts = 10){
  #' mRNA filter function
  #' Give a mRNA dds object and the min row count sum.
  
  
  keep <- rowSums(counts(dds.object))> row.counts
  dds.object <- dds.object[keep,]
  return(dds.object)  
}


filterMicroRNA <- function(dds.object, filter.isomirs = F, row.sum = 10){
  #' microRNA filter function
  #' Give a microRNA dds object and boolean isomir filter option. Filters out rows which count sum is less than 10 and isomirs.
  
  if(filter.isomirs){
    numbers <- c("1","2","3","4","5","6","7","8","9","0")
    keep <- rowSums(counts(dds.object))> row.sum & 
      is.element(str_sub(rownames(dds.object),1,1), numbers)==FALSE
  }
  else{
    keep <- rowSums(counts(dds.object))> row.sum 
  }
  dds.object <- dds.object[keep,]
  return(dds.object)  
  
}

miRNA_target_table <- function(microRNA.res, mRNA.res, database ="mirtarbase", keys = "validated", functional.mti.weak = T){
  
  mRNA.res$row <- str_split( mRNA.res$row,  "\\.", simplify = T)[,1]
  
  #getting target genes from the mirtarbase db
  multi.db <- get_multimir(mirna = microRNA.res$row, 
                           target = mRNA.res$row,
                           table = database, summary = TRUE, 
                           use.tibble = TRUE, predicted.cutoff = 10, 
                           predicted.cutoff.type = "p")
  print(multi.db)
  
  
  result <- multiMiR::select(multi.db, 
                             keytype = "type",
                             keys = keys,
                             columns = columns(multi.db))
  if(length(result)== 0){
    warning("No microRNA target genes in the mRNA gene list")
    break()
  }
  #View(result)
  if (!functional.mti.weak){
    result <- filter(result, support_type == "Functional MTI")
    
  }
  
  #selecting the mirbase microRNA ids and respective target genes 
  result <- result[,c("mature_mirna_id", "target_ensembl","target_symbol")]
  
  
  #merging result table with the miRNA result table to identify fold change for microRNA
  result <- merge(result, microRNA.res, by.x = "mature_mirna_id", by.y = "row")
  print(result)
  
  #merging result table with mRNA res table to identify fold change for mRNA
  result <- merge(result, mRNA.res, by.x ="target_ensembl", by.y ="row")
  print(result)
  
  #getting rid of duplicates
  result <- result[ !duplicated(result) ,]
  
  return(result)
}

