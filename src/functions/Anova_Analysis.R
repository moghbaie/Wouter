# Mehrnoosh Oghbaie
## self
## 09/21/2020
## Preparing data from MaxQuant or ...

## Data preparation consists of four stages:
##  1. Remove contaminants and reverse proteins
##  2. Log transformation
## t-test and fold change analysis ...


set.seed(123)

###########################################################################
### Proteomic template
###########################################################################
TemplateProtein <- R6Class("TemplateProtein",
                           list(
                             input.dir = NA,
                             input = NA,
                             df = NA,
                             df_log = NA,
                             df_log_revised = NA,
                             df_imputed = list(),
                             df_significant = list(),
                             df_significant_woI = list(),
                             condition = NA,
                             dp = NA,
                             dp_log = NA,
                             count_records = NA,
                             regulation_table= list(),
                             enrichment_table =list()
                             
                           )
)


TemplateProtein$set("public","importInput", function(input.dir){
  print("Importing Maxquant output")
  self$input.dir <- input.dir
  self$input <- read.delim(file.path(self$input.dir, 'proteinGroups.txt'), header = TRUE)
  self$condition <- gsub("LFQ.intensity.","",unique(colnames(self$input)[grepl("LFQ.intensity",colnames(self$input))]))
  self$dp <- read.delim(file.path(self$input.dir,"Phospho (STY)Sites.txt"))
  
  invisible(self)
}
)

TemplateProtein$set("public","removeContaminant", function(){
  print("Removing contaminant proteins and reversed sequences")
  cols =colnames(self$input)
  dl = self[["input"]] 
  
  dl$Peptide.counts..unique. <-  apply(dl,1,function(x) strsplit(as.character(x["Peptide.counts..unique."]),";")[[1]][1])
  
  
  dl$Protein.ID <- apply(dl,1,function(x) strsplit(as.character(x["Protein.IDs"]),";")[[1]][1])
  dl$Gene.name <- apply(dl,1,function(x) strsplit(as.character(x["Gene.names"]),";")[[1]][1])
  dl$Gene.name[is.na(dl$Gene.name)] <- dl$Protein.ID[is.na(dl$Gene.name)]
  
  self$df <- dl[!(dl$Potential.contaminant=="+"|dl$Reverse=="+"|dl[["Q.value"]]>0.05) ,
                c("Protein.ID","Gene.name","Peptide.counts..unique.","Fasta.headers",cols[grepl("LFQ.intensity",cols)])] 
  colp <- colnames(self$dp)
  dlp <-  self$dp
  dlp <- dlp %>% dplyr::mutate(Protein = as.character(Leading.proteins),
                               Gene.names = as.character(Gene.names),
                               Potential.contaminant = as.character(Potential.contaminant))
  np <- dim(dlp)[1]
  self$dp <- dlp[!(dlp$Potential.contaminant=="+") ,]
  print(paste((np-dim(self$dp)[1]),"records out of",np,"records were contaminants or reversed sequences and were removed from monoclonal phospho table."))
  
  
  invisible(self)
})



TemplateProtein$set("public","transformData", function(){
  print("Log transformation from intensities")
  
  self[["df_log"]] <- self$df
  self[["df_log"]][,colnames(self[["df_log"]])[grepl("LFQ.intensity.",colnames(self[["df_log"]]))]] <- data.frame(apply(self$df[,colnames(self$df)[grepl("LFQ.intensity.",colnames(self$df))]],2, function(x) ifelse(x>0,log2(as.numeric(x)),0)))
  
  dlp <- self$dp[,c("Gene.names", "Protein" ,"PEP","Fasta.headers",
                    colnames(self$dp)[grepl(paste0("Intensity."), colnames(self$dp))])]
  
  colp <- colnames(dlp) 
  dlp[,grepl("Intensity.",colp)&!grepl("___", colp)] <- apply(dlp[,grepl("Intensity.",colp)&!grepl("___", colp)],2,function(x) ifelse(x!=0, log2(x),0))
  
  self[["dp_log"]] <- dlp
  
  invisible(self)
})

TemplateProtein$set("public","choosingConditions", function(Comparison){
  print("Removing redundant replicates")
  dx <- self$df_log 
  self$count_records <- data.frame(apply(dx[,grepl("LFQ.intensity",colnames(dx))],2, function(x) length(x[x!=0])))
  
  
  
  for(i in 1:dim(Comparison)[1]){
    for(j in 1:4){
      col <- paste0(c(Comparison[i,c(1,3,2)],j), collapse="_")
      cols <- colnames(dx)[grepl(col, colnames(dx))]
      if(length(cols)>1){
        tableCount <- data.frame(apply(dx[,cols],2, function(x) length(x[x!=0])))
        
        dx <- dx[,!colnames(dx) %in% rownames(tableCount)[tableCount[,1]!=max(tableCount[,1])]]
      }
    }
  }
  
  colnames(dx) <- gsub("_revB|_revC|revD","",colnames(dx))
  self$df_log_revised <- dx
  invisible(self)
})

TemplateProtein$set("public","visualize", function(Comparison){
  print("Performing imputation and t-test")
  df <- self$df_log_revised 
  cells <- unique(Comparison[,1])
  extract <- unique(Comparison[,3])
  grid <- expand.grid(cells,extract)
  
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  doParallel::registerDoParallel(cl)
  
  for(i in 1:dim(grid)[1]){
    grid[i,]
    comp <- Comparison[Comparison[,1]==grid[i,1] &Comparison[,3]==grid[i,2],]
    case <- paste0(comp[1,c(1,3,2)], collapse="_")
    control <- paste0(comp[2,c(1,3,2)], collapse="_")
    case_reps <- colnames(df)[grepl(case, colnames(df))]
    control_reps <- colnames(df)[grepl(control, colnames(df))]
    ## Removing records with zero intensities in all replicates
    dfnz <- df[apply(df[,c(case_reps, control_reps)],1,sum)!=0,]
    print(paste('imputing',paste(grid[i,])))
    
    
    
    ddd <- foreach(seed=1:100, .combine='rbind', .packages = c("doParallel"), .export =c("impute", "calculate_stats_nonzeros","impute_all_zeros","min_col","count_zeros","na_zeros_impute","impute_partial_zeros")) %dopar% {
      set.seed(seed)
      dfi <- dfnz[,!grepl("LFQ.intensity", colnames( dfnz))]
      cols2 <- c(colnames(dfi), case_reps, control_reps)
      dfi <- cbind(dfi, matrix(NA, ncol=length(c(case_reps, control_reps)), nrow= dim(dfi)[1]))
      colnames(dfi) <- cols2
      dfi[,case_reps] <- impute(dfnz[case_reps], amm = "2", pmm = "6")
      dfi[,control_reps] <- impute(dfnz[control_reps], amm = "2", pmm = "6")
      return(dfi)
    }
    
    dfi <- ddd %>% group_by(Protein.ID, Gene.name, Peptide.counts..unique., Fasta.headers) %>%
      summarise(across(c(control_reps, case_reps), ~ mean(.x)))
    
    
    self$df_imputed[[paste0(as.character(unlist(grid[i,])), collapse="_")]] <- dfi
    
    dff <- dfnz[,!grepl("LFQ.intensity", colnames( dfnz))]
    dff[["avg.Case.log.intensity"]] <- apply(dfnz[,case_reps],1, function(x) ifelse(!is.na(mean(x[x!=0])), mean(x[x!=0]),0)) 
    dff[["avg.Control.log.intensity"]] <-apply(dfnz[,control_reps],1, function(x) ifelse(!is.na(mean(x[x!=0])), mean(x[x!=0]),0))
    dff[["avg.Case.num.reps"]] <- apply(dfnz[,case_reps],1, function(x) length(x[x!=0])) 
    dff[["avg.Control.num.reps"]] <-apply(dfnz[,control_reps],1, function(x) length(x[x!=0]))
    dff[["log2foldChange"]] <- as.numeric(dff[["avg.Case.log.intensity"]]) - as.numeric(dff[["avg.Control.log.intensity"]])
    dff["p.value"] <- NA
    
    tt <- foreach(k=1:dim(dff)[1], .combine=rbind)%do%{
      t.test(dfnz[k,case_reps],dfnz[k,control_reps])$p.value
    }
    dff[["p.value"]] <- data.frame(tt)$tt
    dff[["p.adjust.value"]] <- p.adjust(dff[["p.value"]] ,method="BH")
    dff[["Significant"]] <- ifelse(((dff[['log2foldChange']] >1 &dff[["avg.Case.num.reps"]]>1)|(dff[['log2foldChange']] < -1 &dff[["avg.Control.num.reps"]]>1))&dff[["p.adjust.value"]]<0.05,'Yes','No')
    
    self$df_significant_woI[[paste0(as.character(unlist(grid[i,])), collapse="_")]] <- dff
    
    dfx <- dfi[,!grepl("LFQ.intensity", colnames( dfi))]
    dfx[["avg.Case.log.intensity"]] <- apply(dfi[,case_reps],1, function(x) ifelse(!is.na(mean(x[x!=0])), mean(x[x!=0]),0)) 
    dfx[["avg.Control.log.intensity"]] <-apply(dfi[,control_reps],1, function(x) ifelse(!is.na(mean(x[x!=0])), mean(x[x!=0]),0))
    dfx[["avg.Case.num.reps"]] <- apply(dfnz[,case_reps],1, function(x) length(x[x!=0])) 
    dfx[["avg.Control.num.reps"]] <-apply(dfnz[,control_reps],1, function(x) length(x[x!=0]))
    dfx[["log2foldChange"]] <- as.numeric(dfx[["avg.Case.log.intensity"]]) - as.numeric(dfx[["avg.Control.log.intensity"]])
    dfx[["log2foldChange.adjusted"]] <- as.numeric(dfx[["avg.Case.log.intensity"]]) - as.numeric(dfx[["avg.Control.log.intensity"]])*(median(as.numeric(dfx[["avg.Case.log.intensity"]]))/median(as.numeric(dfx[["avg.Control.log.intensity"]])))
    
    dfx["p.value"] <- NA
    
    tt2 <- foreach(k=1:dim(dff)[1], .combine=rbind)%do%{
      t.test(dfi[k,case_reps],dfi[k,control_reps])$p.value
    }
    dfx[["p.value"]] <- data.frame(tt2)$tt2
    dfx[["p.adjust.value"]] <- p.adjust(dfx[["p.value"]] ,method="BH")
    dfx[["Significant"]] <- ifelse(((dfx[['log2foldChange']] >1 &dfx[["avg.Case.num.reps"]]>1)|(dfx[['log2foldChange']] < -1 &dfx[["avg.Control.num.reps"]]>1))&dfx[["p.adjust.value"]]<0.05,'Yes','No')
    dfx[["Significant.adjusted"]] <- ifelse(((dfx[['log2foldChange.adjusted']] >1 &dfx[["avg.Case.num.reps"]]>1)|(dfx[['log2foldChange.adjusted']] < -1 &dfx[["avg.Control.num.reps"]]>1))&dfx[["p.adjust.value"]]<0.05,'Yes','No')
    
    
    bfe <- foreach(k=1:dim(dfi)[1], .combine=rbind)%do%{
      dh <- data.frame(rbind(cbind(unlist(dfi[k,case_reps]),rep("case",length(dfi[k,case_reps]))),
                             cbind(unlist(dfi[k,control_reps]),rep("control",length(dfi[k,control_reps])))))
      colnames(dh) <- c("intensity", "condition")
      dh$intensity <- as.numeric(dh$intensity)
      return(abs(1/ttestBF(formula = intensity ~ condition, data = dh)@bayesFactor$bf))
    }
    
    dfx[["bf.error"]] <- data.frame(bfe)$bfe
    dfx[["SignificantB"]] <- ifelse((dfx$bf.error < (1/3)) & (dfx$log2foldChange > 1 & dfx$avg.Case.num.reps > 1), "Yes","No")
    
    dfx[["SignificantB.adjusted"]] <- ifelse((dfx$bf.error < (1/3)) & (dfx$log2foldChange.adjusted > 1 & dfx$avg.Case.num.reps > 1), "Yes","No")
    
    self$df_significant[[paste0(as.character(unlist(grid[i,])), collapse="_")]] <- dfx
    
  }
  
  parallel::stopCluster(cl)
  
  invisible(self)
})





TemplateProtein$set("public","anovaAnalysis", function(Comparison){
  dt <- self$count_records
  colnames(dt) <- 'count'
  dt$replicate <- gsub("_revB|_revC|_revD","",rownames(dt))
  dt1 <- dt%>%group_by(replicate)%>% summarize(count = max(count))
  
  dt1$replicate <- gsub("LFQ.intensity.","",dt1$replicate)
  
  
  p1 <- ggplot(dt1, (aes(x = replicate, y = count))) +
    geom_bar(aes(fill = count), stat = "identity") +
    scale_color_gradient(low="blue", high="green")+
    xlab("sample") +
    ylab("Number of proteins with non zero intensity") +
    ggtitle("Number of proteins with intensity > 0") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    scale_fill_viridis_b()
  
  ggsave(file= file.path("../img/proteomics/intensity.barplot.pdf"), p1, width=9, height=9, dpi=100)
  ggsave(file= file.path("../img/proteomics/intensity.barplot.png"), p1, width=9, height=9, dpi=100)
  
  
  for(cell in names(self$df_imputed)){
    cols <- colnames(self$df_imputed[[cell]])[grepl("LFQ.",colnames(self$df_imputed[[cell]]))]
    dz <- self$df_log_revised[apply(self$df_log_revised[,cols], 1,sum )>0, cols]
    dw <- self$df_imputed[[cell]][,grepl("LFQ", colnames(self$df_imputed[[cell]]))]
    
    dzm <- melt(dz)
    dwm <- melt(dw)
    dzm[["category"]] <- "PreImputation"
    dwm[["category"]] <- "PostImputation"
    
    dm <- rbind(dzm,dwm)
    colnames(dm) <- c("condition","value", "category")
    dm$condition <- gsub("LFQ.intensity.","", dm$condition)
    
    pd <- ggplot(dm, aes(x = value, y = condition, color = category, point_color = category, fill = category)) +
      geom_density_ridges(
        jittered_points = TRUE, scale = .95, rel_min_height = .01,
        point_shape = "|", point_size = 3, size = 0.25,
        position = position_points_jitter(height = 0)
      ) +
      scale_y_discrete(expand = c(0, 0)) +
      scale_x_continuous(expand = c(0, 0), name = "log LFQ intensity") +
      scale_fill_manual(values = c("#D55E0050", "#0072B250"), labels = c("Post Imputation", "Pre Imputation")) +
      scale_color_manual(values = c("#D55E00", "#0072B2"), guide = "none") +
      scale_discrete_manual("point_color", values = c("#D55E00", "#0072B2"), guide = "none") +
      coord_cartesian(clip = "off") +
      guides(fill = guide_legend(
        override.aes = list(
          fill = c("#D55E00A0", "#0072B2A0"),
          color = NA, point_color = NA)
      )
      ) +
      ggtitle(paste("Density plot before and after imputation", cell)) +
      theme_ridges(center = TRUE)+
      theme(axis.text.y = element_text(size = 10))
    
    ggsave(file= file.path("../img/proteomics/",paste0(cell,".density.imputation.pdf")), pd, width=8, height=9, dpi=100)
    ggsave(file= file.path("../img/proteomics/",paste0(cell,".density.imputation.png")), pd, width=8, height=9, dpi=100)
    
    
  }
  
  
  
  for(cell in names(self$df_significant)){
    dx <- self$df_significant[[cell]]
    shade = data.frame(x1=c(-1, 1), 
                       x2=c(-Inf, Inf),
                       y1=c(-log10(0.05), -log10(0.05)), 
                       y2=c(Inf, Inf))
    fold_cutoff = 1
    pvalue_cutoff = 0.05
    
    
    order <- c(Comparison[Comparison[,1]==strsplit(cell,"_")[[1]][1]&Comparison[,3]==strsplit(cell,"_")[[1]][2],2])
    
    p <- ggplot(dx) +
      theme_bw()+
      geom_rect(data=shade, 
                mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill='grey95')+
      geom_point(data  = subset(dx, (Significant=="No")), aes(x =log2foldChange, y = -log10(p.adjust.value),colour=Significant), alpha = 0.5,size=2)+
      geom_point(data  = subset(dx, (Significant=="Yes")), aes(x =log2foldChange, y = -log10(p.adjust.value),colour=Significant), alpha = 0.8,size=3)+
      geom_vline(xintercept = fold_cutoff, col = "blue")+
      geom_vline(xintercept = -fold_cutoff, col = "blue")+
      geom_hline(yintercept = -log10(pvalue_cutoff), col = "green")+
      ggtitle(paste(cell, "(",order[[2]],"~",order[[1]],")"))+
      scale_colour_manual(values=c("gray","red"))+
      geom_text_repel(data  = subset(dx, (Significant=="Yes")),
                      aes(x=log2foldChange, y=-log10(p.adjust.value),label=Gene.name),
                      segment.alpha =0.35,
                      size = 2.5 )
    
    ggsave(file=paste0("../img/proteomics/",paste0(cell,"_", paste0(order, collapse=".")),"_volcanoplot.pdf"), p, width=11, height=9, dpi=200)
    ggsave(file=paste0("../img/proteomics/",paste0(cell,"_", paste0(order, collapse=".")),"_volcanoplot.png"), p, width=11, height=9, dpi=200)
  }
  
  
  for(cell in names(self$df_significant_woI)){
    dx <- self$df_significant_woI[[cell]]
    shade = data.frame(x1=c(-1, 1), 
                       x2=c(-Inf, Inf),
                       y1=c(-log10(0.05), -log10(0.05)), 
                       y2=c(Inf, Inf))
    fold_cutoff = 1
    pvalue_cutoff = 0.05
    
    
    order <- c(Comparison[Comparison[,1]==strsplit(cell,"_")[[1]][1]&Comparison[,3]==strsplit(cell,"_")[[1]][2],2])
    
    p <- ggplot(dx) +
      theme_bw()+
      geom_rect(data=shade, 
                mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill='grey95')+
      geom_point(data  = subset(dx, (Significant=="No")), aes(x =log2foldChange, y = -log10(p.adjust.value),colour=Significant), alpha = 0.5,size=2)+
      geom_point(data  = subset(dx, (Significant=="Yes")), aes(x =log2foldChange, y = -log10(p.adjust.value),colour=Significant), alpha = 0.8,size=3)+
      geom_vline(xintercept = fold_cutoff, col = "blue")+
      geom_vline(xintercept = -fold_cutoff, col = "blue")+
      geom_hline(yintercept = -log10(pvalue_cutoff), col = "green")+
      ggtitle(paste(cell, "(",order[[2]],"~",order[[1]],")"))+
      scale_colour_manual(values=c("gray","red"))+
      #geom_text_repel(data  = subset(dx, (Significant=="Yes")),
      #                aes(x=log2foldChange, y=-log10(p.adjust.value),label=Gene.name),
      #                segment.alpha =0.35,
      #                size = 2.5 )
    
    ggsave(file=paste0("../img/proteomics/",paste0(cell,"_", paste0(order, collapse=".")),"_volcanoplotWOI.pdf"), p, width=11, height=9, dpi=200)
    ggsave(file=paste0("../img/proteomics/",paste0(cell,"_", paste0(order, collapse=".")),"_volcanoplotWOI.png"), p, width=11, height=9, dpi=200)
  }
  
  
  
  png(file=paste0("../img/proteomics/pvalue.distribution.png"), width=800, height=900)
  
  par(mfrow=c(4,2))
  for(cell in names(self$df_significant)){
    
    order <- c(Comparison[Comparison[,1]==strsplit(cell,"_")[[1]][1]&Comparison[,3]==strsplit(cell,"_")[[1]][2],2])
    
    hist(self$df_significant_woI[[cell]]$p.value, breaks="FD", main=paste("Pre imputation", cell,paste0(order, collapse=".vs.")), xlab="")
    hist(self$df_significant_woI[[cell]]$p.adjust.value, breaks="FD", add=T, col='red')
    
    hist(self$df_significant[[cell]]$p.value, breaks="FD", main=paste("Post imputation",cell,paste0(order, collapse=".vs.")), xlab="")
    hist(self$df_significant[[cell]]$p.adjust.value, breaks="FD", add=T, col='red')
    legend("top", inset=.02, title="",
           c("Pre Imputation","Post Imputation"), fill=c("grey","red"), horiz=TRUE, cex=1.5)
  }
  
  dev.off()
  
  invisible(self)
})

TemplateProtein$set("public","compareRnaProtein", function(RNAInput){
  self$enrichment_table[["CC"]] <- list()
  self$enrichment_table[["BP"]] <- list()
  for(name in names(RNAInput$rna_toptable)){
    print(name)
    cell <- unique(Comparison[Comparison[,2] %in% strsplit(name, ".vs.")[[1]],1])
    
    dd <- RNAInput$rna_toptable[[name]]
    dd <- dd[!grepl("NA.",dd$genes),]
    dd$regulate <- ifelse((dd$logFC>1 &dd$FDR<0.05),1,ifelse((dd$logFC< -1 &dd$FDR<0.05),-1,0))
    dd <- dd[dd$regulate!=0,]
    
    dp1 <- self$df_significant[[paste0(c(cell,"insol"), collapse="_")]]
    dp1$regulate <- ifelse((dp1$Significant=="Yes"&dp1$log2foldChange>= 1),1,ifelse((dp1$Significant=="Yes"&dp1$log2foldChange < -1),-1,0))
    
    dp2 <- self$df_significant[[paste0(c(cell,"wcl"), collapse="_")]]
    dp2$regulate <- ifelse((dp2$Significant=="Yes"&dp2$log2foldChange>= 1),1,ifelse((dp2$Significant=="Yes"&dp2$log2foldChange < -1),-1,0))
    
    dt <- data.frame(Gene = dd$genes,
                     RNA= dd$regulate,
                     ProtinInsol=dp1$regulate[match(dd$genes,dp1$Gene.name)],
                     ProtinWcl= dp2$regulate[match(dd$genes,dp2$Gene.name)])
    
    self$regulation_table[[name]] <- dt
    
    dt1 <- dt[(!is.na(dt$ProtinInsol))|(!is.na(dt$ProtinWcl)),]
    dt1 <- dt1[apply(dt1[,c("ProtinInsol","ProtinWcl")],1, function(x) (!is.na(x["ProtinInsol"]))|(!is.na(x["ProtinWcl"]))),]
    dt1.m <- melt(dt1)
    dt1.m$value <- ifelse(dt1.m$value==1,"Upregulated", ifelse(dt1.m$value==-1,"Downregulated","Unchanged"))
    dt1.m$value <- as.factor(dt1.m$value)
    
    qd <-  ggplot(data=dt1.m, aes( variable,Gene, fill= value)) + 
      geom_tile() +
      scale_fill_manual(values=c("seagreen3","white","red"), limits=c("Downregulated","Unchanged",'Upregulated'))+
      ggtitle(paste(cell, name))+
      xlab("Experiment")
    
    ggsave(file=paste0("../img/proteomics/",cell,".",name,".diff_tilePlot.pdf"), qd, width=7, height=7+round((dim(dt1.m)[1]-50)/30), dpi=200)
    ggsave(file=paste0("../img/proteomics/",cell,".",name,".diff_tilePlot.png"), qd, width=7, height=7+round((dim(dt1.m)[1]-50)/30), dpi=200)
    
    
    gene <- unique(c(dd$genes,dp1$Gene.name, dp2$Gene.name))
    countTable <- data.frame(gene = gene)
    countTable$RNA <- dd$regulate[match(countTable$gene,dd$genes)]
    countTable$ProteinInsol <- dp1$regulate[match(countTable$gene,dp1$Gene.name)]
    countTable$ProteinWcl <- dp2$regulate[match(countTable$gene,dp2$Gene.name)]
    
    countTable[is.na(countTable)] <- 0
    countTable[countTable== -1] <- 1
    
    png(paste0("../img/proteomics/","Significant_Venndiagram_",cell, name,".png"),width = 500, height = 500)
    venn(countTable[,-1],zcolor = "style", sncs =0.9, ilcs =0.85)
    dev.off()
    
    
    for(pro in c("CC","BP")){
      geneID <- mapIds(org.Hs.eg.db, countTable$gene[countTable$RNA==1], 'ENTREZID', 'SYMBOL')
      gene.df <- bitr(geneID , fromType = "ENTREZID",
                      toType = c("ENSEMBL", "SYMBOL"),
                      OrgDb = org.Hs.eg.db)
      
      ego <- enrichGO(gene          = gene.df$ENTREZID,
                      # universe      = gene.df$SYMBOL,
                      OrgDb         = org.Hs.eg.db,
                      ont           = pro,
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05)
      
      egos <-  clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
      
      #####
      geneID2 <- mapIds(org.Hs.eg.db, countTable$gene[countTable$ProteinInsol==1], 'ENTREZID', 'SYMBOL')
      gene.df2 <- bitr(geneID2 , fromType = "ENTREZID",
                       toType = c("ENSEMBL", "SYMBOL"),
                       OrgDb = org.Hs.eg.db)
      
      ego2 <- enrichGO(gene          = gene.df2$ENTREZID,
                       # universe      = gene.df$SYMBOL,
                       OrgDb         = org.Hs.eg.db,
                       ont           = pro,
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05)
      
      egos2 <-  clusterProfiler::simplify(ego2, cutoff=0.7, by="p.adjust", select_fun=min)
      
      ####
      geneID3 <- mapIds(org.Hs.eg.db, countTable$gene[countTable$ProteinWcl==1], 'ENTREZID', 'SYMBOL')
      gene.df3 <- bitr(geneID3 , fromType = "ENTREZID",
                       toType = c("ENSEMBL", "SYMBOL"),
                       OrgDb = org.Hs.eg.db)
      
      ego3 <- enrichGO(gene          = gene.df3$ENTREZID,
                       # universe      = gene.df$SYMBOL,
                       OrgDb         = org.Hs.eg.db,
                       ont           = pro,
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05)
      
      egos3 <-  clusterProfiler::simplify(ego3, cutoff=0.7, by="p.adjust", select_fun=min)
      
      egos@result$Cluster <- "RNA"
      egos2@result$Cluster <- "ProteinInsol"
      if(dim(egos3@result)[1]>0){
        egos3@result$Cluster <- "ProteinWcl"
      }
      
      egoList <- rbind(egos@result, egos2@result, egos3@result)%>% filter(qvalue<=0.05&p.adjust<=0.01)
      egoList$GeneRatio <- round(unlist(lapply(egoList$GeneRatio, function(x) as.numeric(strsplit(x,"/")[[1]][1])/as.numeric(strsplit(x,"/")[[1]][2]))),2)
      
      
      egoList$Symbol <- NA
      for(i in 1:dim(egoList)[1]){
        
        geneID <- strsplit(egoList$geneID[i],"/")[[1]]
        symbols <- mapIds(org.Hs.eg.db, geneID, 'SYMBOL', 'ENTREZID')
        egoList$Symbol[i] <- paste(unname(symbols), collapse="|")
      }
      
      self$enrichment_table[[pro]][[name]] <- egoList
      
      qe <- ggplot(egoList, aes(x=Description, y=Cluster)) + 
        geom_point(stat='identity', aes(color=p.adjust, size=GeneRatio)) +
        labs(title=paste("Comparison between RNA-Seq andproteomic data"), 
             subtitle=ifelse(pro=="CC", "Cellular component","Biological Process")) + 
        scale_color_gradient(low="blue", high="red")+
        coord_flip()+
        theme_bw()+
        theme(axis.text.x = element_text(angle = 90))
      
      ggsave(file=paste0("../img/proteomics/",cell,".",name,".Enrichment.Comparison.",pro,".pdf"), qe, width=5+(dim(egoList)[1]/30), height=7+(dim(egoList)[1]/10), dpi=200)
      ggsave(file=paste0("../img/proteomics/",cell,".",name,".Enrichment.Comparison.",pro,".png"), qe, width=5+(dim(egoList)[1]/30), height=7+(dim(egoList)[1]/10), dpi=200)
      
    }
    
  }
invisible(self)
})


TemplateProtein$set("public","writeFiles", function(){
  
  
  if(file.exists("../table/significant_table.xlsx")){
    file.remove("../table/significant_table.xlsx")}
  names <- names(self$df_significant)
  
  for(name in names){
    xlsx::write.xlsx(as.data.frame(self$df_significant[[name]]), file="../table/significant_table.xlsx",sheetName=name, row.names=FALSE, append=TRUE)
  }
  
  if(file.exists("../table/significant_tableWOI.xlsx")){
    file.remove("../table/significant_tableWOI.xlsx")}
  
  for(name in names){
    xlsx::write.xlsx(as.data.frame(self$df_significant_woI[[name]]), file="../table/significant_tableWOI.xlsx",sheetName=name, row.names=FALSE, append=TRUE)
  }
  
  if(file.exists("../table/enrichment_tableCC.xlsx")){
    file.remove("../table/enrichment_tableCC.xlsx")}
  names <- names(self$enrichment_table$CC)
  for(name in names){
    xlsx::write.xlsx(as.data.frame(self$enrichment_table$CC[[name]]), file="../table/enrichment_tableCC.xlsx",sheetName=name, row.names=FALSE, append=TRUE)
  }
  
  if(file.exists("../table/enrichment_tableBP.xlsx")){
    file.remove("../table/enrichment_tableBP.xlsx")}
  names <- names(self$enrichment_table$BP)
  for(name in names){
    xlsx::write.xlsx(as.data.frame(self$enrichment_table$BP[[name]]), file="../table/enrichment_tableBP.xlsx",sheetName=name, row.names=FALSE, append=TRUE)
  }
  
  if(file.exists("../table/regulation_table.xlsx")){
    file.remove("../table/regulation_table.xlsx")}
  names <- names(self$regulation_table)
  for(name in names){
    xlsx::write.xlsx(as.data.frame(self$regulation_table[[name]]), file="../table/regulation_table.xlsx",sheetName=name, row.names=FALSE, append=TRUE)
  }
  
  invisible(self)
})


TemplateProtein$set("public","drawScatterplot", function(RNAInput){
  sig <- c(self$df_significant$HEK293T_insol$Gene.name[self$df_significant$HEK293T_insol$Significant=="Yes"],
           self$df_significant$HEK293T_wcl$Gene.name[self$df_significant$HEK293T_wcl$Significant=="Yes"])
  
  
  
  dx <- merge(self$df_significant$HEK293T_insol[,c("Gene.name","log2foldChange")], 
              self$df_significant$HEK293T_wcl[,c("Gene.name","log2foldChange")], by="Gene.name", all=T)
  
  colnames(dx) <- c("Gene.name","insol","wcl")
  
  pS1 <- ggplot(dx, aes(x = insol, y = wcl)) +
    geom_point(data =subset(dx, !Gene.name %in% sig) ,colour="gray", size=2)+
    geom_point(data =subset(dx, Gene.name %in% sig) ,colour="red", size=3)+
    geom_text_repel(data  = subset(dx, Gene.name %in% sig),
                  aes(x = insol, y = wcl,label = Gene.name),
                  segment.alpha =0.35,
                  size = 2.5 ) +
    ggtitle("HEK293T")
  
  ggsave(file=paste0("../img/proteomics/Hek293T_scatterplot.pdf"), pS1, width=11, height=9, dpi=200)
  ggsave(file=paste0("../img/proteomics/Hek293T_scatterplot.png"), pS1, width=11, height=9, dpi=200)
  
  
  
  #self$df_significant$U2OS_insol
  #self$df_significant$U2OS_wcl
  
  
  sig <- c(self$df_significant$U2OS_insol$Gene.name[self$df_significant$U2OS_insol$Significant=="Yes"],
           self$df_significant$U2OS_wcl$Gene.name[self$df_significant$U2OS_wcl$Significant=="Yes"])
  
  
  
  dx <- merge(self$df_significant$U2OS_insol[,c("Gene.name","log2foldChange")], 
              self$df_significant$U2OS_wcl[,c("Gene.name","log2foldChange")], by="Gene.name", all=T)
  
  colnames(dx) <- c("Gene.name","insol","wcl")
  
  pS2 <- ggplot(dx, aes(x = insol, y = wcl)) +
    geom_point(data =subset(dx, !Gene.name %in% sig) ,colour="gray", size=2)+
    geom_point(data =subset(dx, Gene.name %in% sig) ,colour="red", size=3)+
    geom_text_repel(data  = subset(dx, Gene.name %in% sig),
                    aes(x = insol, y = wcl,label = Gene.name),
                    segment.alpha =0.35,
                    size = 2.5 ) +
    ggtitle("U2OS")
  
  ggsave(file=paste0("../img/proteomics/U2OS_scatterplot.pdf"), pS2, width=11, height=9, dpi=200)
  ggsave(file=paste0("../img/proteomics/U2OS_scatterplot.png"), pS2, width=11, height=9, dpi=200)
  
  
  
  ###########################################################
  #self$df_significant$U2OS_insol
  #self$df_significant$U2OS_wcl
  #RNAInput$rna_toptable$ATM.KO.vs.WT
  
  sig <- c(self$df_significant$U2OS_insol$Gene.name[self$df_significant$U2OS_insol$Significant=="Yes"],
           self$df_significant$U2OS_wcl$Gene.name[self$df_significant$U2OS_wcl$Significant=="Yes"],
           RNAInput$rna_toptable$ATM.KO.vs.WT$genes[abs(RNAInput$rna_toptable$ATM.KO.vs.WT$logFC)>1&RNAInput$rna_toptable$ATM.KO.vs.WT$FDR<0.05])
  
  sig_wcl <- c(self$df_significant$U2OS_wcl$Gene.name[self$df_significant$U2OS_wcl$Significant=="Yes"],
               RNAInput$rna_toptable$ATM.KO.vs.WT$genes[abs(RNAInput$rna_toptable$ATM.KO.vs.WT$logFC)>1&RNAInput$rna_toptable$ATM.KO.vs.WT$FDR<0.05])
  
  sig_wcl_common <- intersect(self$df_significant$U2OS_wcl$Gene.name[self$df_significant$U2OS_wcl$Significant=="Yes"],
                              RNAInput$rna_toptable$ATM.KO.vs.WT$genes[abs(RNAInput$rna_toptable$ATM.KO.vs.WT$logFC)>1&RNAInput$rna_toptable$ATM.KO.vs.WT$FDR<0.05])
  
  
  sig_insol <- c(self$df_significant$U2OS_insol$Gene.name[self$df_significant$U2OS_insol$Significant=="Yes"],
                 RNAInput$rna_toptable$ATM.KO.vs.WT$genes[abs(RNAInput$rna_toptable$ATM.KO.vs.WT$logFC)>1&RNAInput$rna_toptable$ATM.KO.vs.WT$FDR<0.05])
  
  sig_insol_common <- intersect(self$df_significant$U2OS_insol$Gene.name[self$df_significant$U2OS_insol$Significant=="Yes"],
                                RNAInput$rna_toptable$ATM.KO.vs.WT$genes[abs(RNAInput$rna_toptable$ATM.KO.vs.WT$logFC)>1&RNAInput$rna_toptable$ATM.KO.vs.WT$FDR<0.05])
  
  
  dx <- merge(self$df_significant$U2OS_insol[,c("Gene.name","log2foldChange")], 
              self$df_significant$U2OS_wcl[,c("Gene.name","log2foldChange")], by="Gene.name", all=T)
  
  dxx <- merge(dx, RNAInput$rna_toptable$ATM.KO.vs.WT[,c("genes","logFC")], by.x="Gene.name", by.y = "genes", all=T)
  
  colnames(dxx) <- c("Gene.name","insol","wcl", "RNA")
  
  shade = data.frame(x1=c(-1, 1,-1,1), 
                     x2=c(-Inf, Inf, -Inf, Inf),
                     y1=c(-1,1,1,-1), 
                     y2=c(-Inf, Inf, Inf, -Inf))
  
  pS3 <- ggplot(data=dxx[dxx$Gene.name %in% sig_wcl, ]) +
    geom_rect(data=shade, 
              mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill='grey95')+
    geom_point(data =subset(dxx, !Gene.name %in% sig_wcl), aes(x = RNA, y = wcl) ,colour="gray", size=2, alpha= 0.5)+
    geom_point(data =subset(dxx, Gene.name %in% sig_wcl) , aes(x = RNA, y = wcl),colour="pink", size=3, alpha=1)+
    geom_point(data =subset(dxx, Gene.name %in% sig_wcl_common), aes(x = RNA, y = wcl) ,colour="red", size=3)+
    geom_text_repel(data  = subset(dxx, Gene.name %in% sig_wcl_common),
                    aes(x = RNA, y = wcl,label = Gene.name),
                    segment.alpha =0.35,
                    size = 4 ) +
    geom_vline(xintercept = 1, col = "blue")+
    geom_vline(xintercept = -1, col = "blue")+
    geom_hline(yintercept = -1, col = "green")+
    geom_hline(yintercept = 1, col = "green")+
    ggtitle("U2OS Wcl")+
    theme_classic()
  
  
  pS4 <- ggplot(dxx[dxx$Gene.name %in% sig_insol, ]) +
    geom_rect(data=shade, 
              mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill='grey95')+
    geom_point(data =subset(dxx, !Gene.name %in% sig_insol), aes(x = RNA, y = insol) ,colour="gray", size=2, alpha= 0.5)+
    geom_point(data =subset(dxx, Gene.name %in% sig_insol) , aes(x = RNA, y = insol),colour="pink", size=3, alpha=1)+
    geom_point(data =subset(dxx, Gene.name %in% sig_insol_common), aes(x = RNA, y = insol) ,colour="red", size=3)+
    geom_text_repel(data  = subset(dxx, Gene.name %in% sig_insol_common),
                    aes(x = RNA, y = insol,label = Gene.name),
                    segment.alpha =0.35,
                    size = 4 ) +
    geom_vline(xintercept = 1, col = "blue")+
    geom_vline(xintercept = -1, col = "blue")+
    geom_hline(yintercept = -1, col = "green")+
    geom_hline(yintercept = 1, col = "green")+
    scale_y_reverse()+
    ggtitle("U2OS Insol")+
    theme_classic()
  
  ps <- grid.arrange(pS3, pS4 , 
                     ncol = 1, nrow = 2)
  ggsave(file=paste0("../img/proteomics/U2OS_RNA_scatterplot.pdf"), ps, width=11, height=11, dpi=200)
  ggsave(file=paste0("../img/proteomics/U2OS_RNA_scatterplot.png"), ps, width=11, height=11, dpi=200)
  
  ################################
  
  #self$df_significant$HEK293T_insol
  #self$df_significant$HEK293T_wcl
  #RNAInput$rna_toptable$CPT.vs.DMSO
  
  sig <- c(self$df_significant$HEK293T_insol$Gene.name[self$df_significant$HEK293T_insol$Significant=="Yes"],
           self$df_significant$HEK293T_wcl$Gene.name[self$df_significant$HEK293T_wcl$Significant=="Yes"],
           RNAInput$rna_toptable$CPT.vs.DMSO$genes[abs(RNAInput$rna_toptable$CPT.vs.DMSO$logFC)>1&RNAInput$rna_toptable$CPT.vs.DMSO$FDR<0.05])
  
  sig_wcl <- c(self$df_significant$HEK293T_wcl$Gene.name[self$df_significant$HEK293T_wcl$Significant=="Yes"],
               RNAInput$rna_toptable$CPT.vs.DMSO$genes[abs(RNAInput$rna_toptable$CPT.vs.DMSO$logFC)>1&RNAInput$rna_toptable$CPT.vs.DMSO$FDR<0.05])
  
  sig_wcl_common <- intersect(self$df_significant$HEK293T_wcl$Gene.name[self$df_significant$HEK293T_wcl$Significant=="Yes"],
                              RNAInput$rna_toptable$CPT.vs.DMSO$genes[abs(RNAInput$rna_toptable$CPT.vs.DMSO$logFC)>1&RNAInput$rna_toptable$CPT.vs.DMSO$FDR<0.05])
  
  
  sig_insol <- c(self$df_significant$HEK293T_insol$Gene.name[self$df_significant$HEK293T_insol$Significant=="Yes"],
                 RNAInput$rna_toptable$CPT.vs.DMSO$genes[abs(RNAInput$rna_toptable$CPT.vs.DMSO$logFC)>1&RNAInput$rna_toptable$CPT.vs.DMSO$FDR<0.05])
  
  sig_insol_common <- intersect(self$df_significant$HEK293T_insol$Gene.name[self$df_significant$HEK293T_insol$Significant=="Yes"],
                                RNAInput$rna_toptable$CPT.vs.DMSO$genes[abs(RNAInput$rna_toptable$CPT.vs.DMSO$logFC)>1&RNAInput$rna_toptable$CPT.vs.DMSO$FDR<0.05])
  
  
  dx <- merge(self$df_significant$HEK293T_insol[,c("Gene.name","log2foldChange")], 
              self$df_significant$HEK293T_wcl[,c("Gene.name","log2foldChange")], by="Gene.name", all=T)
  
  dxx <- merge(dx, RNAInput$rna_toptable$CPT.vs.DMSO[,c("genes","logFC")], by.x="Gene.name", by.y = "genes", all=T)
  
  colnames(dxx) <- c("Gene.name","insol","wcl", "RNA")
  
  shade = data.frame(x1=c(-1, 1,-1,1), 
                     x2=c(-Inf, Inf, -Inf, Inf),
                     y1=c(-1,1,1,-1), 
                     y2=c(-Inf, Inf, Inf, -Inf))
  
  pS3 <- ggplot(data=dxx[dxx$Gene.name %in% sig_wcl, ]) +
    geom_rect(data=shade, 
              mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill='grey95')+
    geom_point(data =subset(dxx, !Gene.name %in% sig_wcl), aes(x = RNA, y = wcl) ,colour="gray", size=2, alpha= 0.5)+
    geom_point(data =subset(dxx, Gene.name %in% sig_wcl) , aes(x = RNA, y = wcl),colour="pink", size=3, alpha=1)+
    geom_point(data =subset(dxx, Gene.name %in% sig_wcl_common), aes(x = RNA, y = wcl) ,colour="red", size=3)+
    geom_text_repel(data  = subset(dxx, Gene.name %in% sig_wcl_common),
                    aes(x = RNA, y = wcl,label = Gene.name),
                    segment.alpha =0.35,
                    size = 4) +
    geom_vline(xintercept = 1, col = "blue")+
    geom_vline(xintercept = -1, col = "blue")+
    geom_hline(yintercept = -1, col = "green")+
    geom_hline(yintercept = 1, col = "green")+
    ggtitle("HEK293T Wcl")+
    theme_classic()
  
  
  pS4 <- ggplot(dxx[dxx$Gene.name %in% sig_insol, ]) +
    geom_rect(data=shade, 
              mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill='grey95')+
    geom_point(data =subset(dxx, !Gene.name %in% sig_insol), aes(x = RNA, y = insol) ,colour="gray", size=2, alpha= 0.5)+
    geom_point(data =subset(dxx, Gene.name %in% sig_insol) , aes(x = RNA, y = insol),colour="pink", size=3, alpha=1)+
    geom_point(data =subset(dxx, Gene.name %in% sig_insol_common), aes(x = RNA, y = insol) ,colour="red", size=3)+
    geom_text_repel(data  = subset(dxx, Gene.name %in% sig_insol_common),
                    aes(x = RNA, y = insol,label = Gene.name),
                    segment.alpha =0.35,
                    size = 4 ) +
    geom_vline(xintercept = 1, col = "blue")+
    geom_vline(xintercept = -1, col = "blue")+
    geom_hline(yintercept = -1, col = "green")+
    geom_hline(yintercept = 1, col = "green")+
    scale_y_reverse()+
    ggtitle("HEK293T Insol")+
    theme_classic()
  
  ps <- grid.arrange(pS3, pS4 , 
                     ncol = 1, nrow = 2)
  ggsave(file=paste0("../img/proteomics/HEK293T_RNA_scatterplot.pdf"), ps, width=11, height=11, dpi=200)
  ggsave(file=paste0("../img/proteomics/HEK293T_RNA_scatterplot.png"), ps, width=11, height=11, dpi=200)
  
  invisible(self)
  
})