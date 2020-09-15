library(R.matlab)
library(stringr)
library(readr)

recon3D <- readMat("Dropbox/Meta_PKN/recon3D_netowrk/Recon3DModel_301.mat")

#get the stochio matrix
S <- as.matrix(recon3D$Recon3DModel[[1]])

#get the list of gene rules for reactions
reaction_list <- recon3D$Recon3DModel[[21]]

#get reversible reactions
lbs <- recon3D$Recon3DModel[[6]]
reversible <- ifelse(lbs < 0, TRUE, FALSE)

#get the reaction ids
reaction_ids <- unlist(recon3D$Recon3DModel[[5]])

#create a dataframe to map reaction indexes with genes
reaction_to_genes <- list()
for(i in 1:length(reaction_list))
{
  if(length(reaction_list[[i]][[1]] > 0))
  {
    genes <- unique(gsub(" and ","_",gsub("[()]","",gsub("[.][0-9]","",strsplit(reaction_list[[i]][[1]], split = " or ")[[1]]))))
    df <- as.data.frame(matrix(NA,length(genes), 2))
    df[,1] <- i
    df[,2] <- genes
    reaction_to_genes[[i]] <- df
  } else
  {
    df <- as.data.frame(matrix(NA,1, 2))
    df[,1] <- i
    df[,2] <- reaction_ids[i]
    reaction_to_genes[[i]] <- df
  }
}
reaction_to_genes <- as.data.frame(do.call(rbind,reaction_to_genes))
# reaction_to_genes$V2 <- paste("Gene_",reaction_to_genes$V2, sep = "")
reaction_to_genes_original <- reaction_to_genes

#get metabolites
metabolites <- unlist(recon3D$Recon3DModel[[2]])
metabolites_names <- unlist(recon3D$Recon3DModel[[15]])
# KEGG <- recon3D$Recon3DModel[[18]]
PubChem <- recon3D$Recon3DModel[[19]]
# for(i in 1:length(metabolites))
# {
#   if(length(KEGG[[i]][[1]]) > 0)
#   {
#     comp <- gsub(".*[[]","",metabolites[i])
#     metabolites[i] <- paste(KEGG[[i]][[1]],comp, sep = "[")
#   }
# }
metab_to_pubchem <- as.data.frame(matrix(NA,length(metabolites), 2))
names(metab_to_pubchem) <- c("name","pubchem")

for(i in 1:length(metabolites))
{
  metab_to_pubchem[i,1] <- metabolites_names[i]
  if(length(PubChem[[i]][[1]]) > 0)
  {
    metab_to_pubchem[i,2] <- PubChem[[i]][[1]]
    comp <- gsub(".*[[]","",metabolites[i])
    metabolites[i] <- paste(PubChem[[i]][[1]],comp, sep = "[")
  }
}
metab_to_pubchem <- metab_to_pubchem[complete.cases(metab_to_pubchem),]
metab_to_pubchem <- unique(metab_to_pubchem)
write_csv(metab_to_pubchem, "~/Dropbox/Meta_PKN/support/metab_to_pubchem.csv")

#We will modify the name of genes as we go so we save the original names
reaction_to_genes <- reaction_to_genes_original

#create the 2 coulmn format network between metabolites and enzymes
reactions_df <- list()

#we do reactions 1 by 1 (colunms of stochiometric matrix)
for(i in 1:length(S[1,]))
{
  print(i)
  #get the reactions stochiometry
  reaction <- S[,i]
  
  #modify gene name so reactions that are catalised by same enzyme stay separated
  reaction_to_genes[reaction_to_genes$V1 == i,2] <- paste(paste("Gene",i, sep = ""),reaction_to_genes[reaction_to_genes$V1 == i,2], sep = "__")
  
  #get the enzymes associated with reaction
  genes <- reaction_to_genes[reaction_to_genes$V1 == i,2]
  

  #get reactant
  reactants <- paste("Metab__",metabolites[reaction == -1],sep = "")
  
  #get the products
  products <- paste("Metab__",metabolites[reaction == 1], sep = "")
  
  #check how many rows(interactions) will be necessary in the 2 column format of this reaction
  number_of_interations <- length(reactants) + length(products)
  
  #now for each enzymes, we create a two column datframe recapitulating the interactions between the metabolites and this enzyme
  reaction_df <- list()
  j <- 1
  for(gene in genes)
  {
    gene_df <- as.data.frame(matrix(NA,number_of_interations,2))
    gene_df$V1 <- c(reactants,rep(gene,number_of_interations-length(reactants))) #reactants followed by the enzyme (the enzyme is repeated asmany time as they are products)
    gene_df$V2 <- c(rep(gene,number_of_interations-length(products)),products) #enzyme(repeated asmany time as they are reactants) followed by products
    
    if(reversible[i]) #if reaction is reversible, we do the same but inverse reactant and products
    {
      gene_df_reverse <- as.data.frame(matrix(NA,number_of_interations,2))
      gene_df_reverse$V1 <- c(rep(paste(gene,"_reverse",sep = ""),number_of_interations-length(products)),products)
      gene_df_reverse$V2 <- c(reactants,rep(paste(gene,"_reverse",sep = ""),number_of_interations-length(reactants)))
      gene_df <- as.data.frame(rbind(gene_df,gene_df_reverse))
    }
    reaction_df[[j]] <- gene_df 
    j <- j+1
  }
  #the individual enzyme dataframes of this reaction are combined into one reaction dataframe
  reaction_df <- as.data.frame(do.call(rbind,reaction_df))
  
  #the reaction dataframe is added to the list of all reaction reaction dataframes
  reactions_df[[i]] <- reaction_df
}
#the individual reaction dataframesare combined into one
reactions_df <- as.data.frame(do.call(rbind,reactions_df))

write_csv(reactions_df, "Dropbox/Meta_PKN/recon3D_netowrk/reaction_network_recon3.csv")

###############################
###############################
###############################

gene_metab_network <- reactions_df
gene_metab_network$V1 <- gsub("[[][a-z][]]","",gene_metab_network$V1)
gene_metab_network$V2 <- gsub("[[][a-z][]]","",gene_metab_network$V2)

nodes <- rep(1,length(unique(c(gene_metab_network$V1, gene_metab_network$V2))))
names(nodes) <- unique(c(gene_metab_network$V1, gene_metab_network$V2))

metabs <- nodes[grepl("Metab__",names(nodes))]
# names(metabs) <- gsub("[[][a-z][]]","",names(metabs))

i <- 1
for(metab in names(metabs))
{
  print(i)
  metabs[metab] <- sum(metab == gene_metab_network$V1)
  metabs[metab] <- metabs[metab] + sum(metab == gene_metab_network$V2)
  i <- i+1
}

metabs_sorted <- sort(metabs, decreasing = T)

# write(as.character(gsub("Metab__","",names(metabs_sorted))),"~/Dropbox/Meta_PKN/recon3D_netowrk/pubchem_compound_list.txt")

coenzymes <- as.character(as.data.frame(read_csv("Dropbox/Meta_PKN/recon3D_netowrk/coenzymes.txt", 
                      col_names = FALSE))[,1]) #see 01b_list_coenzymes.R

metabs_sorted <- metabs_sorted[!(gsub("Metab__","",gsub("[[][a-z][]]","",names(metabs_sorted))) %in% coenzymes)]
to_write <- as.data.frame(metabs_sorted)
to_write$pubchem <- row.names(to_write)
write_csv(to_write[,c(2,1)], "~/Dropbox/Meta_PKN/support/metab_cardinality.csv")
metabs_sorted <- metabs_sorted[metabs_sorted < 350] #smallest number of connections before we find important metabolties

gene_metab_network_no_cofac <- reactions_df[gsub("[[][a-z][]]","",reactions_df$V1) %in% names(metabs_sorted) | gsub("[[][a-z][]]","",reactions_df$V2) %in% names(metabs_sorted),]

write_csv(gene_metab_network_no_cofac, "Dropbox/Meta_PKN/recon3D_netowrk/reaction_network_recon3_no_cofact.csv")

metabs <- unique(c(gene_metab_network_no_cofac$V1, gene_metab_network_no_cofac$V2))
metabs <- metabs[grepl("Metab",metabs)]
metabs <- gsub("Metab__","",metabs)
metabs <- gsub("[[].*","",metabs)
metabs <- unique(metabs)

metab_to_pubchem <- metab_to_pubchem[metab_to_pubchem$pubchem %in% metabs,]
metab_to_pubchem <- unique(metab_to_pubchem)

write_csv(metab_to_pubchem, "~/Dropbox/Meta_PKN/support/metab_to_pubchem_cofactfiltered.csv")
