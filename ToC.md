Downloaded the 14 tumor counts table from GEOAccession, Accession number GSE136661

### On PuTTY

     tar xvf GSE136661_RAW.tar
tar xvf untars and extracts the selected tumor count files as gz files.

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
gzip -d was used to decompress all the tumor counts files extracted from the tarred file. Decompressing these files allows them to be opened and used in other programs.

### In RStudio

     library(tidyverse)
     library(readr)
     library(DESeq2)
     library(ggplot2)
The library command allows us to access and use these packages in RStudio.

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
After importing all the  counts tables into RStudio by tranferring them to the Techtmann server using WinSCP, colnames was used to rename the columns to prepare the tables for DESeq2. [1] is used to signify the first column and [2] is used to signify the second.

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
In order to have all the tumor counts for each tumor type in one table the merge function was used. This commands allows for the merging of two tables at a time, thus it was used multiple times until all the counts for each tumor grade were in one file, named WHO_I and WHO_III.

    WHO_III = WHO_III[-(1:2),]
    WHO_I = WHO_I[-(1:2),]
This command was to remove two unnecessary rows from both counts files.

Imported excel table with id columns for all 14 tumors, named Ia-Ig Expression and IIIa-IIIg, to match the counts tables. dex were columns named WHO I or WHO III, depending on the type of tumor

    tumorcounts <- merge(WHO_I,WHO_III)
    tumorcounts <- as.data.frame(tumorcounts)
    metadata <- as.data.frame(metadata)
    all(names(tumorcounts)[-1]==metadata$id)
These commands were run to prepare the data for DESeq2. The merge command was again used to merge our two counts files into one named tumor counts. as.data.frame is used to convert both our tumorcounts table and metadata table into data frames which are necessary to run DESeq2. The last command is used to make sure that all our column names in tumorcounts (except the first, which is the gene identification) are the same as our id columns in our metadata.

     dds.data <- DESeqDataSetFromMatrix(countData=tumorcounts, 
                                   colData=metadata, 
                                   design=~dex, 
                                   tidy=TRUE)
This command is used to create a DESeq2 object that is needed to run DESeq2. We specify that our countData is our tumorcounts file, our colData is our metadata file, design=~dex was used to specify what we're comparing, which in this case is WHO grade I tumors to WHO grade III tumors. The last part of the command, tidy=TRUE, tells DESeq2 that the first column should be the rownames (our genes), as a column named row.

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
dds <- DESeq(dds.data) runs DESeq2 on the file created in the last command. The res command creates different versions of the results, which are used for plotting in later commands. In the first res, contrast is used to tell our results which comparison to make. In this case, we want our results to compare expression our dex, WHO I and WHO III. arrange does just as it says and arranges the results based on what is specified. In res9 and res11, it was used to find the three most differentially expressed genes for both tumor grades. filter is used to remove results deemed not significant, in this case anything with a log2FoldChange greater than or equal to 2 or less than or equal to -2, an adjusted p-value of 0.05, and later a basemean greater than or equal to 100.

### Volcano Plot

    ggplot() + 
      geom_point(data=res8, aes(x=log2FoldChange, y=-log10(padj)), color='salmon') + 
      geom_point(data=res10, aes(x=log2FoldChange, y=-log10(padj)), color='darkturquoise') +
      ggtitle("WHO Type I vs WHO Type III Differential Expression") +
      xlab("Log2 Fold Change") +
      ylab("-Log10 Adjusted P-value")+
      theme(axis.text=element_text(size=14),
       axis.title=element_text(size=14,face="bold"))
To display all the significant differentially expressed genes from the two tumor types ggplot2 was used to make a volcano plot. For this command, we use geom_point to specify that we want to use points and which results table should be plotted. aes is used to state which part of the results table will be used for our x and y axes and color allows us to pick the color of the points, here we have selected salmon for WHO I tumors and dark turquoise for WHO III tumors. ggtitle is used to give the plot a name and xlab and ylab are used to create labels for the x and y axes. The theme part of this command is used to increase the font size on the plot and bold the title.

### WHO Type I Boxplots

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

### WHO Type III Boxplots

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
Using plotCounts the expression of a single gene was able to be plotted in a boxplot format. This command was used to create plots of our three most significantly expressed genes in the WHO I and WHO III tumors. dds signifies using our dds data, gene= is used to state the gene we wish to plot, ingroup= is used to determine what comparison we want to make in this case it is between our dex groups, WHO I and WHO III, and returnData=TRUE allows us to return the raw data for this gene. We combined the basic plotCounts command with ggplot2 using %>%, allowing us to make a better plot. aes is used to specify that we want to compare our two dex groups and have them represent our x-axis, while our y-axis has our counts for the gene. The addition of scale_y_log10() is used so that the difference in expression can be more easily seen. Like in the volcano plot, ggtitle, xlab, ylab, and theme are used to name the graph, label the axes, and change the font size. We also added theme(legend.position = "none") to state that we do not want a legend. 

    citation("ggplot2")
    citation("DESeq2")
The citation command was used to obtain the information needed to cite both ggplot2 and DESeq2.
