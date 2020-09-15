library(readr)

### Import the high confidence interactions of STITCH
STITCH_900 <- as.data.frame(read_delim("Dropbox/Meta_PKN/STITCH_network/STITCH_900.tsv", 
                         "\t", escape_double = FALSE, trim_ws = TRUE))

####THIS SECTION TO REMOVE TEXTMINING
X9606_protein_chemical_links_detailed_v5_0 <- as.data.frame(read_delim("Dropbox/Meta_PKN/STITCH_network/9606.protein_chemical.links.detailed.v5.0.tsv", 
                                                         "\t", escape_double = FALSE, trim_ws = TRUE))

not_text_mining <- X9606_protein_chemical_links_detailed_v5_0[X9606_protein_chemical_links_detailed_v5_0$combined_score >= 900,]
rm(X9606_protein_chemical_links_detailed_v5_0)

not_text_mining$ID <- paste(not_text_mining$chemical, not_text_mining$protein , sep = "_")
not_text_mining$ID_reverse <- paste(not_text_mining$protein, not_text_mining$chemical, sep = "_")

not_text_mining <- not_text_mining[(not_text_mining$experimental + not_text_mining$prediction + not_text_mining$prediction) >= 900,]

### We only care about allosteric interactions
STITCH_900 <- STITCH_900[STITCH_900$a_is_acting,]
STITCH_900$ID <- paste(STITCH_900$item_id_a, STITCH_900$item_id_b, sep = "_")
STITCH_900 <- STITCH_900[STITCH_900$ID %in% not_text_mining$ID | STITCH_900$ID %in% not_text_mining$ID_reverse,]
STITCH_900 <- STITCH_900[,-7]
####END OF REMOVE TEXTMINING

### We need to map the Ensembl Ids to NCBI gene ids

## We make a vector of uniue ensembl IDs
prots <- unique(c(STITCH_900$item_id_a, STITCH_900$item_id_b))
prots <- prots[grepl("9606[.]ENSP", prots)]

## We remove the taxon ID so the IDs can be mapped
prots <- as.data.frame(cbind(prots, gsub("9606[.]","",prots)))

##We use biomart package to map the ensembl ids to NCBI ones
library(biomaRt)

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

G_list <- getBM(filters = "ensembl_peptide_id", 
                attributes = c("ensembl_peptide_id",'hgnc_symbol','entrezgene_id', "description"),
                values = prots$V2, mart = ensembl)
names(G_list)[1] <- "V2"

## We put all the ID versions into one dataframe
prots <- merge(prots, G_list, by = "V2")
prots <- prots[prots$entrezgene != "",]

## We create a named vector to make the id conversion more efficient
prots_vec <- prots$entrezgene
names(prots_vec) <- prots$prots

## Ids are converted in the STITCH interaciton dataframe
for(i in 1:2)
{
  for(j in 1:length(STITCH_900[,1]))
  {
    if(STITCH_900[j,i] %in% names(prots_vec))
    {
      STITCH_900[j,i] <- prots_vec[STITCH_900[j,i]]
    }
    else
    {
      STITCH_900[j,i] <- gsub("CID[a-z]0*","Metab__",STITCH_900[j,i])
    }
  }
}

### We converted the interactions of STITCH into a SIF format
STITCH_900$sign <- ifelse(STITCH_900$action == "inhibition", -1, 1)
STITCH_900 <- STITCH_900[grepl("Metab__",STITCH_900$item_id_a),]

STITCH_900 <- STITCH_900[,c(1,7,2)]
names(STITCH_900) <- c("source","sign","target")

### Save the SIF network as csv
write_csv(STITCH_900,"Dropbox/Meta_PKN/STITCH_network/STITCH_900_sif.csv")
