
library(clusterProfiler)
library(org.Hs.eg.db)

############Gene Ontology Enrichemnt Analysis###########
#to identify whether certain functional categories or biological processes are significantly associated with a set of genes of interest. 
#This is particularly useful when you have a list of genes resulting from experiments such as differential gene expression analysis

####1- Go classification###

# Load the geneList data from the DOSE package
data(geneList, package = "DOSE")

# Select genes with absolute values greater than 2, This type of filtering is often used 
#in differential gene expression analysis to focus on genes that have undergone substantial changes in expression levels,
#regardless of whether the change is upregulation (positive) or downregulation (negative)
#By choosing a threshold of absolute value greater than 2, 
#you are focusing on genes that have undergone at least a twofold change 
#in expression (either up or down) between the two conditions. This threshold helps identify genes that 
#have potentially meaningful and biologically significant changes in expression levels.
gene <- names(geneList)[abs(geneList) > 2]


# Perform gene ID conversion using bitr
gene.df <- bitr(gene, fromType = "ENTREZID", toType = c("ENSEMBL", "SYMBOL"), OrgDb = org.Hs.eg.db)

# Display the first few rows of the resulting gene.df data frame
head(gene.df)



#readable: When set to TRUE, it indicates that the GO term names should be converted into more readable forms.
#ont: Specifies the ontology to perform enrichment analysis on. 
#In this case, it's set to "CC" for the Cellular Component ontology.
ggo <- groupGO(gene     = gene,
               OrgDb    = org.Hs.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)


#the output provide information about the enriched GO terms and their statistical significance in relation to your selected gene list.
head(ggo)




######2- GO over-representation test #######
# to determine whether specific categories of genes or gene products are overrepresented in a set of genes of interest 
#compared to what would be expected by chance to determine whether specific categories of genes or gene products are overrepresented 
#in a set of genes of interest compared to what would be expected by chance

ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)


gene_symbols <- bitr(gene.df$ENSEMBL, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)








#1-keytype: Specifies the type of gene identifiers you are using. In this case, 
#it's set to 'ENSEMBL' since you are using Ensembl IDs
#2-pAdjustMethod: Specifies the method for adjusting p-values for multiple testing. 
#"BH" refers to the Benjamini-Hochberg method.
ego2 <- enrichGO(gene         = gene_symbols$SYMBOL,
                 OrgDb         = org.Hs.eg.db,
                 #keytype       = 'ENSEMBL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)


#This function maps the Ensembl IDs associated with enriched GO terms to their corresponding gene symbols, 
#making the results more interpretable.
ego2 <- setReadable(ego2, OrgDb = org.Hs.eg.db)





