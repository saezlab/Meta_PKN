library(readr)

reaction_network_recon3_no_cofact <- as.data.frame(read_csv("Dropbox/Meta_PKN/recon3D_netowrk/reaction_network_recon3_no_cofact.csv"))

STITCH_900_sif <- as.data.frame(read_csv("Dropbox/Meta_PKN/STITCH_network/STITCH_900_sif.csv"))

omni_network <- as.data.frame(read_delim("Dropbox/Meta_PKN/omnipath_netowrk/omni_PPI_clean.sif", 
                      "\t", escape_double = FALSE, trim_ws = TRUE))

omni_prots <- unique(c(omni_network$source, omni_network$target))

library(biomaRt)

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

G_list <- getBM(filters = 'hgnc_symbol', 
                attributes = c("ensembl_peptide_id",'hgnc_symbol','entrezgene', "description"),
                values = omni_prots, mart = ensembl)

G_list <- G_list[,c(2,3)]
G_list <- unique(G_list)
G_list_vec <- G_list$entrezgene
names(G_list_vec) <- G_list$hgnc_symbol

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

#######################

metabs <- unique(c(reaction_network_recon3_no_cofact$V1, reaction_network_recon3_no_cofact$V2))
metabs <- metabs[grepl("Metab__",metabs)]

metabs <- as.data.frame(cbind(metabs, gsub("[[][a-z][]]","",metabs)))

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

#################

genes <- unique(c(reaction_network_recon3_no_cofact$V1, reaction_network_recon3_no_cofact$V2))
genes <- genes[grepl("Gene",genes)]
genes <- as.data.frame(cbind(genes,gsub("Gene[0-9]+__","",genes)))
genes$V2 <- gsub("_reverse","",genes$V2)
genes <- unique(genes)

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
meta_nodes <- meta_nodes[meta_nodes$V2 %in% omni_network$source | meta_nodes$V2 %in% omni_network$target,]
meta_nodes$sign <- 1
meta_nodes <- meta_nodes[,c(2,3,1)]

###############
STITCH_900_sif_compartiments_filtered <- STITCH_900_sif_compartiments[
  STITCH_900_sif_compartiments$V3 %in% meta_nodes[1,] |
  STITCH_900_sif_compartiments$V3 %in% omni_network$source |
    STITCH_900_sif_compartiments$V3 %in% omni_network$target,
]

names(STITCH_900_sif_compartiments_filtered) <- c("source","interaction","target")

names(meta_nodes) <- c("source","interaction","target")

reaction_network_recon3_no_cofact$interaction <- 1
reaction_network_recon3_no_cofact <- reaction_network_recon3_no_cofact[,c(1,3,2)]
names(reaction_network_recon3_no_cofact) <- c("source","interaction","target")

names(omni_network) <- c("source","interaction","target")

meta_network <- do.call(rbind,list(STITCH_900_sif_compartiments_filtered,meta_nodes, reaction_network_recon3_no_cofact, omni_network))

write_csv(meta_network, "Dropbox/Meta_PKN/result/meta_network.csv")
