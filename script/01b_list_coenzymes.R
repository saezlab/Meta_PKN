library(RJSONIO)
library(httr)
library(readr)
library(snowfall)
library(parallel)
library(rlist)

get_MeSH_class <- function(cid)
{
  prolog <- 'https://pubchem.ncbi.nlm.nih.gov/rest/pug'
  input <- '/compound/cid'
  output <- '/classification/JSON'
  
  qurl <- paste0(prolog, input, output)
  
  cont <- content(POST(qurl,body = list("cid" = paste(cid, collapse = ','))), type = 'text', encoding = 'UTF-8')
  
  cont <- fromJSON(cont)
  kegg <- FALSE
  if(names(cont) != "Fault")
  {
    cont <- cont[[1]][[1]]
    
    for(i in 1:length(cont))
    {
      if(cont[[i]]$SourceName == "MeSH")
      {
        kegg_class <- cont[[i]]
        kegg <- TRUE
      }
    }
  }
  if(kegg)
  {
    return(list(cid,kegg_class))
  }
  else
  {
    return(cid)
  }
}

metab_list <- as.data.frame(
  read_csv("Dropbox/Meta_PKN/recon3D_netowrk/pubchem_compound_list.txt", 
                       col_names = FALSE))

metab_vec <- as.character(metab_list[,1])
metab_vec <- metab_vec[!grepl("[a-zA-Z]",metab_vec)]

classification_list <- mclapply(metab_vec,get_MeSH_class, mc.cores = 3)

# classification_list <- list()
# for(i in 1:length(metab_list[,1]))
# {
#   print(i)
#   classification_list[[i]] <- get_kegg_class(metab_list[i,])
# }

classification_list <- classification_list[unlist(lapply(classification_list,is.list))]
classification_list <- lapply(classification_list,list.flatten)

cofactors <- unlist(lapply(classification_list,function(x) {sum(grepl("Coenzymes",x)) > 0}))
cofactors <- classification_list[cofactors]
cofactors <- lapply(cofactors, function(x){x[[1]]})
cofactors <- unlist(cofactors)

write(as.character(cofactors),"~/Dropbox/Meta_PKN/recon3D_netowrk/coenzymes.txt")
