library(readr)

#Import recon3D reaction network in sif format
reaction_network_recon3_no_cofact <- as.data.frame(read_csv("Dropbox/Meta_PKN/recon3D_netowrk/reaction_network_recon3_no_cofact.csv"))

#Import STITCH allosteric interactions in SIF format
STITCH_900_sif <- as.data.frame(read_csv("Dropbox/Meta_PKN/STITCH_network/STITCH_900_sif.csv"))

#Import omnipath
url <- "http://omnipathdb.org/interactions?types=post_translational,transcriptional&datasets=omnipath,pathwayextra,dorothea&fields=sources,references,curation_effort,dorothea_level,type&genesymbols=yes"

download_omnipath <- function(){
  
  read.table(url, sep = '\t', header = TRUE)
  
}

omni_network <- download_omnipath()
omni_network <- omni_network[omni_network$consensus_stimulation != 0 | omni_network$consensus_inhibition != 0,]
omni_network <- omni_network[,c(3,4,9,10)]
omni_network$Interaction <- omni_network$consensus_stimulation - omni_network$consensus_inhibition
cons_0 <- omni_network[omni_network$Interaction == 0,]
omni_network <- omni_network[omni_network$Interaction != 0,]
cons_0$consensus_inhibition <- 0
cons_0$Interaction <- 1
omni_network <- as.data.frame(rbind(omni_network,cons_0))
cons_0$consensus_inhibition <- 1
cons_0$consensus_stimulation <- 0
cons_0$Interaction <- -1
omni_network <- as.data.frame(rbind(omni_network,cons_0))
omni_network <- omni_network[,c(1,5,2)]

#Get a vector of omnipath proteins
omni_prots <- unique(c(omni_network$source, omni_network$target))

#Map omnipath proteins to NCBI gene ids (entrez)
library(biomaRt)

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

G_list <- getBM(filters = 'hgnc_symbol', 
                attributes = c("ensembl_peptide_id",'hgnc_symbol','entrezgene_id', "description"),
                values = omni_prots, mart = ensembl)

G_list <- G_list[,c(2,3)]
G_list <- unique(G_list)
G_list_vec <- G_list$entrezgene
names(G_list_vec) <- G_list$hgnc_symbol

#Convert Gene symbols of omnipath to NCBI Ids
for(i in c(1,3))
{
  for(j in 1:length(omni_network[,1]))
  {
    if(omni_network[j,i] %in% names(G_list_vec))
    {
      omni_network[j,i] <- G_list_vec[omni_network[j,i]]
    }
  }
}

recon_nodes <- unique(c(reaction_network_recon3_no_cofact$V1, reaction_network_recon3_no_cofact$V2))
recon_nodes <- recon_nodes[!grepl("Metab",recon_nodes)]
recon_nodes <- recon_nodes[!grepl("_reverse",recon_nodes)]
recon_nodes <- recon_nodes[!grepl("__[A-Za-z]",recon_nodes)]
recon_nodes <- recon_nodes[!grepl("[A-Za-z]$",recon_nodes)]

recon_nodes <- gsub("Gene[0-9]+__","",recon_nodes)

omni_network <- omni_network[!(omni_network$source_genesymbol %in% recon_nodes),]
### We need to ad comprtments to the metabolites in STITCH so we can link them to recon3D

#make a vector of metaoblites in STITCH
metabs <- unique(c(reaction_network_recon3_no_cofact$V1, reaction_network_recon3_no_cofact$V2))
metabs <- metabs[grepl("Metab__",metabs)]

#We create a dataframe mapping metabolites names with compartments
metabs <- as.data.frame(cbind(metabs, gsub("[[][a-z][]]","",metabs)))

#We create a new STITCH SIF dataframe where the metabolites avec the compartment suffixes
row_list <- list()
for(i in 1:length(STITCH_900_sif[,1]))
{
  sub_metabs <- metabs[metabs$V2 == STITCH_900_sif[i,1],]
  if(length(sub_metabs[,1]) > 0)
  {
    sub_STITCH <- as.data.frame(matrix(NA,length(sub_metabs[,1]),3))
    sub_STITCH[,1] <- sub_metabs[,1]
    sub_STITCH[,2] <- STITCH_900_sif[i,2]
    sub_STITCH[,3] <- STITCH_900_sif[i,3]
    # print(sub_STITCH)
    row_list[[i]] <- sub_STITCH 
  }
  else
  {
    sub_STITCH <- as.data.frame(matrix(NA,0,3))
    row_list[[i]] <- sub_STITCH 
  }
}

STITCH_900_sif_compartiments <- as.data.frame(do.call(rbind,row_list))

### Genes of the reaction network are uniquelly identified with respect to the reaction they belong too. There are also gene complexes. Thus we need to map them to their original gene name to conect them to the other networks.

# We create a dataframe conecting unique gene_reactio nidentifiers to genes.
genes <- unique(c(reaction_network_recon3_no_cofact$V1, reaction_network_recon3_no_cofact$V2))
genes <- genes[grepl("Gene",genes)]
genes <- as.data.frame(cbind(genes,gsub("Gene[0-9]+__","",genes)))
genes$V2 <- gsub("_reverse","",genes$V2)
genes <- unique(genes)

# We expand this dataframe by splitting the gene complexes. We assume that each complexe only one of it's gene to be activated for the complexe to be activated. Cohenrently, only one of the member need to be inhibited for the complexe to be inhibited.
row_list <- list()
for(i in 1:length(genes[,1]))
{
  if(grepl("[0-9]_[0-9]",genes[i,2]))
  {
    complexe_members <- unlist(strsplit(genes[i,2], split = "_"))
    complexe_df <- as.data.frame(matrix(NA,length(complexe_members),2))
    complexe_df[,1] <- genes[i,1]
    complexe_df[,2] <- complexe_members
    names(complexe_df) <- names(genes[i,])
    row_list[[i]] <- complexe_df
  }
  else
  {
    row_list[[i]] <- genes[i,]
  }
}

meta_nodes <- as.data.frame(do.call(rbind,row_list))

#We filter out that are not in omnipath
meta_nodes <- meta_nodes[(meta_nodes$V2 %in% omni_network$source | meta_nodes$V2 %in% omni_network$target) | meta_nodes$V2 %in% STITCH_900_sif_compartiments$V3,]
meta_nodes$sign <- 1
meta_nodes <- meta_nodes[,c(2,3,1)]

#### We just need to filter out nodes of STITCH that are not in omnipath and then connect everything together

# Filter out nodes of STITCH that are not connected to the other networks
STITCH_900_sif_compartiments_filtered <- STITCH_900_sif_compartiments[
  STITCH_900_sif_compartiments$V3 %in% meta_nodes[1,] |
    STITCH_900_sif_compartiments$V3 %in% omni_network$source |
    STITCH_900_sif_compartiments$V3 %in% omni_network$target,
  ]

# Make names consistent and crbind everything
names(STITCH_900_sif_compartiments_filtered) <- c("source","interaction","target")

names(meta_nodes) <- c("source","interaction","target")

reaction_network_recon3_no_cofact$interaction <- 1
reaction_network_recon3_no_cofact <- reaction_network_recon3_no_cofact[,c(1,3,2)]
names(reaction_network_recon3_no_cofact) <- c("source","interaction","target")

names(omni_network) <- c("source","interaction","target")

meta_network <- do.call(rbind,list(STITCH_900_sif_compartiments_filtered,meta_nodes, reaction_network_recon3_no_cofact, omni_network))
meta_network <- meta_network[complete.cases(meta_network),]
# Export network
meta_network <- meta_network[meta_network$source != meta_network$target,]
meta_network <- meta_network[c(length(meta_network[,1]),1:(length(meta_network[,1])-1)),]
write_csv(meta_network, "Dropbox/Meta_PKN/result/meta_network_fullomni_metfiltered.csv")

meta_network_carnival_ready <- meta_network
meta_network_carnival_ready$source <- gsub("[[]","___",meta_network_carnival_ready$source)
meta_network_carnival_ready$target <- gsub("[[]","___",meta_network_carnival_ready$target)

meta_network_carnival_ready$source <- gsub("[]]","____",meta_network_carnival_ready$source)
meta_network_carnival_ready$target <- gsub("[]]","____",meta_network_carnival_ready$target)

meta_network_carnival_ready <- unique(as.data.frame(
  rbind(meta_network_carnival_ready[12158,],meta_network_carnival_ready)
))

write_csv(meta_network_carnival_ready, "Dropbox/Meta_PKN/result/meta_network_carnival_ready_fullomni_metfiltered.csv")
