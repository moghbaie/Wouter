# Mehrnoosh Oghbaie
## self
## 09/24/2020
## Differential Expression analysis ...

################################################################################
### RNA -Seq template
################################################################################

TemplateRna <- R6Class("TemplateRna",
                       list(
                         input.dir = NA,
                         count_records = NA,
                         rna_count_reads = NA,
                         rna_target = NA,
                         geneExp = NA,
                         edgeList = list(),
                         edgeListNormalized = list(),
                         diffExpress = list(),
                         rna_toptable = list(),
                         rna_toptags = list()
                         
                       )
)


TemplateRna$set("public","importInput", function(){
  
  print("Importing RNA_Seq reads")
  countmtx <- read.csv("../input/countmtx.csv", row.names=1)
  
  # import the target file
  target <- as.data.frame(read_excel("../input/targetfile.xlsx")) # add more information to the targetfile in excel
  rownames(target) <- target$sampleID
  target$group <- apply(target,1, function(x) paste0(c(x["Cell"],x["condition"], x["BiologicalReplicate"]),collapse="_"))
  colnames(countmtx) <- target$group[match(colnames(countmtx),target$sampleID)]
  
  geneExp <- read.delim("../input/Homo_sapiens.GRCh38.91_shorter.txt", header=FALSE)
  colnames(geneExp) <- c("Ensembl_ID", "Gene_sym")
  
  self$rna_count_reads <- countmtx
  self$rna_target <- target
  self$geneExp <- geneExp
  
  invisible(self)
}
)

TemplateRna$set("public","makeEdgeList", function(){
  print("Making the Edge list and normalize it.")
  target <- self$rna_target
  geneExp <- self$geneExp
  countmtx <- self$rna_count_reads
    
  for(cell in unique(target$Cell)){
    split <- countmtx[,colnames(countmtx) %in% target$group[target$Cell==cell]]
    geneMat <- geneExp[match(rownames(split), geneExp$Ensembl_ID),]
    rownames(split) <- make.names(geneMat$Gene_sym, unique=TRUE)
    dgList <- obtainFilteredDGE(split, target[target$Cell==cell,], rownames(split))
    self$edgeList[[cell]] <- dgList
    logCPMcounts <- cpm(dgList,log=TRUE)
    colnames(logCPMcounts) <- target$group[target$Cell==cell]
    self$edgeListNormalized[[cell]] <- logCPMcounts
    
    df <- as.data.frame(dgList$samples)
    
    p1 <- ggplot(df, aes(x = as.factor(group), y = lib.size)) + # basic of ggplot, define x and y
      geom_bar(stat="identity", fill = "darksalmon") + # add a barplot
      ggtitle("Barplot of library sizes") + # add title
      theme_classic() + # select color theme
      ylab("Library size") + # y axis label
      xlab("Sample") + # x axis label
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) # rotate x labels
    
    
    # barplot of number of expressed genes per sample
    nGene <- data.frame(matrix(NA, nrow = length(colnames(split)), ncol = 1)) # maak leeg dataframe
    rownames(nGene) <- target$group[target$Cell==cell] # rownames word sample names
    colnames(nGene) <- "genes"
    for (i in 1:ncol(split)){ # per sample
      nGene[i,1] <- length(which(split[,i] > 0)) # tellen hoeveel genen een expressie hebben van > 0
    }
    
    p2 <- ggplot(nGene, (aes(x = rownames(nGene), y = nGene$genes))) +
      geom_bar(stat = "identity", fill = "green4") +
      xlab("sample") +
      ylab("Number of genes expressed") +
      ggtitle("Number of genes with count > 0") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
    
    p <-  grid.arrange(arrangeGrob(p1,
                                    p2,
                                    ncol=2, widths=c(1,1)))
    
    ggsave(file= file.path("../img/rna-seq/",paste0(cell,".library.size.pdf")), p, width=15, height=8, dpi=100)
    ggsave(file= file.path("../img/rna-seq/",paste0(cell,".library.size.png")), p, width=15, height=8, dpi=100)
 
    
    df <- melt(logCPMcounts)
    colnames(df) <- c("Gene.name","condition","logCPM")
    df$condition <- as.character(df$condition)
    q <- ggplot(df, aes(x=condition, y=logCPM)) + 
      geom_boxplot()+
      geom_hline(yintercept=median(logCPMcounts),col="blue")+
      theme_classic()+
      theme(axis.text.x = element_text(angle = 90))+
      ggtitle("Boxplots of logCPMs")
    
    ggsave(file= file.path("../img/rna-seq/",paste0(cell,".boxplot.pdf")), q, width=5, height=9, dpi=100)
    ggsave(file= file.path("../img/rna-seq/",paste0(cell,".boxplot.png")), q, width=5, height=9, dpi=100)
  
    
    project.pca <- prcomp(t(logCPMcounts)) # logCPM matrix moet getransposed worden, vandaar de t
    summary(project.pca)
    
    # plot
    dg <- project.pca$x
    dg <- dg[,1:2]
    dg <- as.data.frame(dg)
    dg <- cbind(dg, dgList$samples)

    # PCA plot colored by library size shows that the separation of the samples is not caused by the library size
    # So that is good because that indicates that the samples segregate based on biological variation instead of technical variation
    # so i hope these 3 groups represent your experimental groups :) 
    qq <- ggplot(dg, aes(PC1, PC2, colour = lib.size)) + # you can change lib.size by any column in your target file
      geom_point(size = 3) + 
      scale_shape_manual(values =  c(21, 23)) +
      geom_text_repel(aes(label = group), nudge_x = 5) +
      xlab(paste("PC1 (", (round((project.pca$sdev[1]^2)/(sum(project.pca$sdev^2))*100, digits = 2)), "%)", sep = "")) +
      ylab(paste("PC2 (", (round((project.pca$sdev[2]^2)/(sum(project.pca$sdev^2))*100, digits = 2)), "%)", sep = "")) +
      theme_classic() +
      ggtitle("PCA logCPM")
    
    ggsave(file= file.path("../img/rna-seq/",paste0(cell,".pcaplot.pdf")), qq, width=9, height=9, dpi=100)
    ggsave(file= file.path("../img/rna-seq/",paste0(cell,".pcaplot.png")), qq, width=9, height=9, dpi=100)
    
    #split1 <- split[,!colnames(split) %in% c("Hek293t_DMSO_1","self$rna_target")]
    #target1 <- target[target$group!=c("Hek293t_DMSO_1","self$rna_target"),]
    #dgList1 <- obtainFilteredDGE(split1, target1[target1$Cell==cell,], rownames(split1))
    #self$edgeList[[cell]] <- dgList1
    #logCPMcounts1 <- cpm(dgList1,log=TRUE)
    #colnames(logCPMcounts1) <- target1$group[target1$Cell==cell]
    #self$edgeListNormalized[[cell]] <- logCPMcounts1
    #self$rna_target <- target1
    
     }
  
  invisible(self)
})


TemplateRna$set("public","diffExpression", function(){
  print("Calculate differential expression")
  target <- self$rna_target
  
  for(cell in unique(target$Cell)){
    dgList <- self$edgeList[[cell]] 
    logCPMcounts <- self$edgeListNormalized[[cell]]
    
    # Creating a design 
    design <- model.matrix(~0 + dgList$samples$condition, data=dgList$samples) # basically making a binary matrix that contains the information which sample is in which experimental group
    
    # estimate dispersions 
    dge.full <- estimateDisp(dgList, design)
    plotBCV(dge.full) 
    fit <- glmFit(dge.full, design) # fitting of the glm
    
    name <- paste0(sort(unique(fit$samples$condition)),collapse=".vs.")
    
    # differential gene expression function, contains all settings you need
    # change p- and logFC thresholds if you like
    DiffExpr <- function(contrast){
      name <- paste0(sort(unique(fit$samples$condition)),collapse=".vs.")
      # do the actual DE analysis
      lrt <- glmLRT(fit, contrast = contrast)
      tags <- topTags(lrt, n=dim(dge.full[[1]])[1], adjust.method="BH", sort.by="logFC") #this is the full table for all genes
      #assign(x = paste("tags", name, sep = ""), value = tags, .GlobalEnv)
      self$rna_toptags[[name]] <- tags
      # extract significant genes (logFC > 1 and p < 0.05)
      toptable <- tags[[1]][which(abs(tags[[1]]$logFC) > 0.5 & tags[[1]]$FDR < 0.05),] # this is the table with only the significant genes
      plotSmear(lrt, de.tags = toptable$genes)
      abline(h=c(-1, 1), col=2)
      title(paste("DEGS for ", name, sep = ""))
      
      #assign(x = paste("topTable", name, sep = ""), value = toptable, .GlobalEnv)
      self$rna_toptable[[name]] <- toptable
      # print info
      return(summary(decideTests(lrt, p = 0.05, lfc = 1, adjust.method = "BH")))
    }
    
    # setting of the contrast is based on the colnames of your design. In this example there are 3 groups.
    # but if you have more than you should also add more values in de contrast vector. 
    colnames(design)
    contrast <- c(1, -1)
    
    
    
    name <- paste0(sort(unique(fit$samples$condition)),collapse=".vs.")
    # now you have defined everything you need so you can use the DiffExpr function
    DiffExpr(contrast) # alles wat rood is in the plot is een DEG
    toptable <- self$rna_toptable[[name]]
    # plot DE genes in a heatmap
    n_allgenes <- logCPMcounts[rownames( toptable), ] # here you can also opt to extract the values from the countspermillion matrix, but logCPM is usually nicer
    n_allgenes <- n_allgenes[order(-abs(toptable$logFC)),]
    n_allgenes <- n_allgenes[!grepl("NA.",rownames(n_allgenes)),]
    
    colnames(n_allgenes) <- target$group[target$Cell==cell]
    colfunc <- colorRampPalette(c("steelblue", "slategray2", "white", "tan1", "firebrick3")) # these are nice colors
    
    pdf( file= file.path("../img/rna-seq/",paste0(cell,".heatmap.pdf")),  width=10, height=10)
    heatmap.2(n_allgenes, 
              col=colfunc(50), # define number of sections in color scale
              keysize = 1,
              scale="row", # does a row z-score
              trace="none", # is really ugly
              cexRow = 0.7, # font size of the row labels
              cexCol = 0.9, # font size of the col labels
              #Colv = "NA", # do unsupervised, meaning you don't do column clustering
              dendrogram = 'both', # or choose 'row' or 'col or 'none
              key.title = NA,
              main = "Differentially expressed genes")
    
    dev.off()
    
    # volcano plot
    # all dots that have a color are significantly differentially expressed between the groups
    # grey dots are genes that are not significant DE (p > 0,05 or logFC > 1)
    qp <- ggplot( toptable, aes(x = logFC, y = -log(FDR))) +
      geom_point(colour = "grey80", size = 0.5) + 
      geom_hline(yintercept = -log(0.05), colour = "grey80") + # these are the sifnigicane thresholds (p < 0.05)
      geom_vline(xintercept = c(-1,1), colour = "grey80") + # logFC thresholds (logFC > 1 or < -1)
      geom_point(data = toptable[[name]], aes(x = logFC, y = -log(FDR)), shape = 21, 
                 colour = "black", 
                 fill = ifelse(toptable$logFC > 1, "violetred3", "mediumaquamarine")) +
      geom_text_repel(data = toptable, label = ifelse(-log(toptable$FDR) > 50, toptable$genes, "")) + # set treshold from when you want to show gene names, or put a # for this line if you dont want to show gene names
      theme_classic() +
      xlim(-15,15) +
      ggtitle(paste0(sort(unique(fit$samples$condition)),collapse=".vs.")) 
    
    ggsave(file= file.path("../img/rna-seq/",paste0(cell,".volcanoplot.pdf")), qp, width=12, height=12, dpi=100)
    ggsave(file= file.path("../img/rna-seq/",paste0(cell,".volcanoplot.png")), qp, width=12, height=12, dpi=100)
  }
  
  invisible(self)
})