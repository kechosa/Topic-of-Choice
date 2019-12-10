Downloaded the 14 tumor counts table from GEOAccession, Accession number GSE136661

     tar xvf GSE136661_RAW.tar
Comment

    gzip -d GSM4054970_17-01-031R_htseq.counts.txt.gz
    gzip -d GSM4054971_17-01-035R_htseq.counts.txt.gz
    gzip -d GSM4054972_18-01-003R_htseq.counts.txt.gz
    gzip -d GSM4054973_18-01-010R_htseq.counts.txt.gz
    gzip -d GSM4054975_18-01-033R_htseq.counts.txt.gz
    gzip -d GSM4054976_18-01-034R_htseq.counts.txt.gz
    gzip -d GSM4054980_18-01-047R_htseq.counts.txt.gz
    gzip -d GSM4054983_DCC04-15-01-029R-T_htseq.counts.txt.gz
    gzip -d GSM4054984_DCC04-15-01-030R-T_htseq.counts.txt.gz
    gzip -d GSM4054986_DCC14-15-01-028R-T_htseq.counts.txt.gz
    gzip -d GSM4054987_17-01-014-R1_htseq.counts.txt.gz
    gzip -d GSM4054988_18-01-004R_htseq.counts.txt.gz
    gzip -d GSM4054989_18-01-001R_htseq.counts.txt.gz
    gzip -d GSM4054990_18-01-002R_htseq.counts.txt.gz
Comment

     library(tidyverse)
     library(readr)
     library(DESeq2)
     library(ggplot2)
Comment

    colnames(WHO_I_a)[1] <- 'Gene'
    colnames(WHO_I_a)[2] <- 'Ia Expression'
    colnames(WHO_I_b)[1] <- 'Gene'
    colnames(WHO_I_b)[2] <- 'Ib Expression'
    colnames(WHO_I_c)[1] <- 'Gene'
    colnames(WHO_I_c)[2] <- 'Ic Expression'
    colnames(WHO_I_d)[1] <- 'Gene'
    colnames(WHO_I_d)[2] <- 'Id Expression'
    colnames(WHO_I_e)[1] <- 'Gene'
    colnames(WHO_I_e)[2] <- 'Ie Expression'
    colnames(WHO_I_f)[1] <- 'Gene'
    colnames(WHO_I_f)[2] <- 'If Expression'
    colnames(WHO_I_g)[1] <- 'Gene'
    colnames(WHO_I_g)[2] <- 'Ig Expression'
    colnames(WHO_III_a)[1] <- 'Gene'
    colnames(WHO_III_a)[2] <- 'IIIa Expression'
    colnames(WHO_III_b)[1] <- 'Gene'
    colnames(WHO_III_b)[2] <- 'IIIb Expression'
    colnames(WHO_III_c)[1] <- 'Gene'
    colnames(WHO_III_c)[2] <- 'IIIc Expression'
    colnames(WHO_III_d)[1] <- 'Gene'
    colnames(WHO_III_d)[2] <- 'IIId Expression'
    colnames(WHO_III_e)[1] <- 'Gene'
    colnames(WHO_III_e)[2] <- 'IIIe Expression'
    colnames(WHO_III_f)[1] <- 'Gene'
    colnames(WHO_III_f)[2] <- 'IIIf Expression'
    colnames(WHO_III_g)[1] <- 'Gene'
    colnames(WHO_III_g)[2] <- 'IIIg Expression'
Comment

    total <- merge(WHO_I_a,WHO_I_b)
    total2 <-merge(WHO_I_c,WHO_I_d)
    total3 <-merge(total, total2)
    total4 <-merge(WHO_I_e,WHO_I_f)
    total5 <-merge(total3,total4)
    WHO_I <-merge(total5, WHO_I_g)

    total6 <-merge(WHO_III_a,WHO_III_b)
    total7 <-merge(WHO_III_c,WHO_III_d)
    total8 <-merge(WHO_III_e,WHO_III_f)
    total9 <-merge(total6,total7)
    total10 <-merge(total9,total8)
    WHO_III <-merge(total10,WHO_III_g)
Comment

    WHO_III = WHO_III[-(1:2),]
    WHO_I = WHO_I[-(1:2),]
Comment

Imported excel table with id columns for all 14 tumors, named Ia-Ig Expression and IIIa-IIIg, to match the counts tables. dex were columns named WHO I and WHO III, depending on the type of tumor

    tumorcounts <- merge(WHO_I,WHO_III)
    tumorcounts <- as.data.frame(tumorcounts)
    metadata <- as.data.frame(metadata)
    all(names(tumorcounts)[-1]==metadata$id)
Comment

     dds.data <- DESeqDataSetFromMatrix(countData=tumorcounts, 
                                   colData=metadata, 
                                   design=~dex, 
                                   tidy=TRUE)
Comment

    dds <- DESeq(dds.data)
    res <- results(dds, contrast=c("dex","WHO I","WHO III"), tidy=TRUE)
    res2 <- arrange(res, padj)
    res3 <- filter(res2, padj<=0.05 & log2FoldChange>=2)
    res4 <- filter(res2, padj<=0.05 & log2FoldChange<=-2)
    res5 <- arrange(res3, log2FoldChange)
    res6 <- arrange(res4, log2FoldChange)
    res7 <- arrange(res3, baseMean)
    res8 <- filter(res2, padj<=0.05, log2FoldChange>=2 & baseMean>=100)
    res9 <- arrange(res8, log2FoldChange)
    res10 <- filter(res2, padj<=0.05, log2FoldChange<=-2 & baseMean>=100)
    res11 <- arrange(res10, log2FoldChange)
    res12 <- filter(res2, padj<=0.05 &log2FoldChange>=0)
    res13 <- filter(res2, padj<=0.05 & log2FoldChange<=0)
Comment

#Volcano Plot

    ggplot() + 
      geom_point(data=res8, aes(x=log2FoldChange, y=-log10(padj)), color='salmon') + 
      geom_point(data=res10, aes(x=log2FoldChange, y=-log10(padj)), color='darkturquoise') +
      ggtitle("WHO Type I vs WHO Type III Differential Expression") +
      xlab("Log2 Fold Change") +
      ylab("-Log10 Adjusted P-value")+
      theme(axis.text=element_text(size=14),
       axis.title=element_text(size=14,face="bold"))
Comment

#WHO Type I Boxplot

      plotCounts(dds, gene= "ENSG00000179344", intgroup="dex", returnData=TRUE) %>%
       ggplot(aes(dex,count))+
        geom_boxplot(aes(fill=dex))+
        scale_y_log10() + 
         ggtitle("HLA-DQB1")+
        xlab("WHO Typing")+
        ylab("Log10 Expression")+
        theme(axis.text=element_text(size=14),
              axis.title=element_text(size=14,face="bold"))+
        theme(legend.position = "none")

      plotCounts(dds, gene= "ENSG00000196735", intgroup="dex", returnData=TRUE) %>%
      ggplot(aes(dex,count))+
       geom_boxplot(aes(fill=dex))+
        scale_y_log10() + 
        ggtitle("HLA-DQA1")+
        xlab("WHO Typing")+
        ylab("Log10 Expression")+
        theme(axis.text=element_text(size=14),
         axis.title=element_text(size=14,face="bold"))+
        theme(legend.position = "none")

      plotCounts(dds, gene= "ENSG00000204287", intgroup="dex", returnData=TRUE) %>%
        ggplot(aes(dex,count))+
        geom_boxplot(aes(fill=dex))+
        scale_y_log10() + 
         ggtitle("HLA-DRA")+
         xlab("WHO Typing")+
         ylab("Log10 Expression")+
          theme(axis.text=element_text(size=14),
               axis.title=element_text(size=14,face="bold"))+
          theme(legend.position = "none") 
Comment

#WHO III

      plotCounts(dds, gene= "ENSG00000167654", intgroup="dex", returnData=TRUE) %>%
       ggplot(aes(dex,count))+
        geom_boxplot(aes(fill=dex))+
        scale_y_log10() + 
       ggtitle("ATCAY")+
       xlab("WHO Typing")+
       ylab("Log10 Expression")+
       theme(axis.text=element_text(size=14),
             axis.title=element_text(size=14,face="bold"))+
        theme(legend.position = "none")

      plotCounts(dds, gene= "ENSG00000180818", intgroup="dex", returnData=TRUE) %>%
      ggplot(aes(dex,count))+
       geom_boxplot(aes(fill=dex))+
       scale_y_log10() + 
       ggtitle("HOXC10")+
       xlab("WHO Typing")+
       ylab("Log10 Expression")+
       theme(axis.text=element_text(size=14),
              axis.title=element_text(size=14,face="bold"))+
       theme(legend.position = "none")

    plotCounts(dds, gene= "ENSG00000250978", intgroup="dex", returnData=TRUE) %>%
     ggplot(aes(dex,count))+
     geom_boxplot(aes(fill=dex))+
     scale_y_log10() + 
      ggtitle("Novel Transcript")+
      xlab("WHO Typing")+
      ylab("Log10 Expression")+
      theme(axis.text=element_text(size=14),
           axis.title=element_text(size=14,face="bold"))+
     theme(legend.position = "none")
Comment

    citation("ggplot2")
    citation("DESeq2")
Comment
