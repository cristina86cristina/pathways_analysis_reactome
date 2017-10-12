##########################################################################################################################
####### file_name : reactome_overlap_heatmap.r
####### description : this script will produce a heatmap of enrichment z-scores with a structured ordered of Reactome pathways 
####### based on Jaccardi index
####### author : Cristina Venturini
###############################################################################################################################

###packages 
library(dplyr)
library(XGR)
library(pheatmap)
library(RColorBrewer)
library(viridis)
### files and parameters
## 2 versions: from gene list or from matrix of z-scores
ontology <- "MsigdbC2REACTOME" ##choose ontology, here Reactome
#for this we just need the miru_clusters_annotated - 2 col: 1 column cluster/comparison 1 column genes 
### i.e.Cluster Gene_name
###     group1 gene1
###     group1 gene2
###     group2 gene3 
mydir <- "~/Downloads/sepsis_variance2/" ##change
myfile <- "miru_annotated_check.csv"     ##change


myannotatedfile <- read.csv(paste(mydir,myfile,sep=""), as.is = TRUE) ##important to have as.is, no factors! 

##this is for the pathway enrichment analysis for Reactome by Group (comparison, cluster, whatever) 
##takes a moment - go a grab some coffee! 
test_pathways <- myannotatedfile %>% 
  group_by(Cluster) %>% 
  do(as.data.frame(tryCatch(xEnrichViewer(xEnricherGenes(data=.$Gene_name,ontology=ontology,p.adjust.method="BH")),error=function(e)NA)))


## we want a nice summary table with the name of the pathways and z-scores for each group
## i.e. name Group1 Group2
##      Immune system 3 NA 
##      Cell cycle NA 2
tidier <- test_pathways  %>% select(name,Cluster,zscore) %>%
  spread(Cluster, zscore)  %>% group_by(name) %>% arrange(name)
#tidier %>% head(8)  #Uncomment this if you want to check

tidier[sapply(tidier, function(x) all(is.na(x)))] <- NULL # in my case I had clusters with no terms, so I need to get rid of them to make my heatmap tidier and prettier

write.csv(tidier,paste(mydir,"tidier_pathway_forheatmap.csv",sep="") #we can write this as .csv so we dont need to rerun the enrichment analysis everytime

##getting ready for the overlap matrix
#tidier <- read.csv("tidier_pathway_forheatmap.csv",as.is = TRUE)

##load Reactome data 
GS <- xRDataLoader(RData.customised=paste('org.Hs.eg', ontology, sep='')) 

new_list <- inner_join(tidier,GS$set_info) ##get names REACTOME_IMMUNE_SYSTEM


######this part is to build the Jaccardi index among Reactome pathways 
#it's here for referenceI have saved the results in a .csv file and can just use that
#this includes ALL Reactome pathways
##source("https://bioconductor.org/biocLite.R")
#biocLite("GeneOverlap")
#library(GeneOverlap)


#gom.obj <- newGOM(GS$gs) ##this takes a while - time for tea 
#jaccard_overlap_m <-  getMatrix(gom.obj, "Jaccard")
#jacc_df <- as.data.frame(jaccard_overlap_m)
#write.csv(jacc_df,"Jaccardi_index_all_reactome.csv")

#I have created Jaccardi_index_all_reactome.csv - much quicker
jacc_df <- read.csv(paste(mydir,"Jaccardi_index_all_reactome.csv",sep=""), row.names = 1) #read file
jaccard_sel1 <- jacc_df[rownames(jacc_df) %in% new_list$setID, ]  #first step sel
jaccard_final <- jaccard_sel1[,colnames(jaccard_sel1) %in% new_list$setID ] #second step selection - we filter the pathways we need 

jacc_m <- as.matrix(jaccard_final)
res <- pheatmap(jacc_m,cluster_rows =T,cluster_cols = F,fontsize = 6) #get structure based on Jacc index

#let's create the matrix of z scores for heatmap
pathways_nice_order <- new_list
#let's make a bit of order - some columns are not numeric and we dont need them
pathways_nice_order$cluster<-NULL; pathways_nice_order$setID<-NULL;pathways_nice_order$namespace<-NULL;pathways_nice_order$distance<-NULL; pathways_nice_order$`NA`<-NULL

#pathways_nice_order[sapply(pathways_nice_order, function(x) all(is.na(x)))] <- NULL
pathways_nice_order[is.na(pathways_nice_order)] <- 0.000000000001 #change NA into 0s
row.names(pathways_nice_order)<-pathways_nice_order$name #make col "name" into rownames
pathways_nice_order$name<-NULL

#voila' - the heatmap! 

pathways_nice_order_col <- colnames(pathways_nice_order)


pdf(paste(mydir,"heatmap_pretty_order.pdf",sep=""),width = 10,height=9)
heatord <- pheatmap(pathways_nice_order,
                    cluster_rows =res$tree_row$order,  ##this is the most important part, do not change this
                    cluster_cols = T,
                    fontsize = 8)
dev.off()





          