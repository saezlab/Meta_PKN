library(readr)

meta_network_carnival_ready <- as.data.frame(
  read_csv("Dropbox/Meta_PKN/result/meta_network_carnival_ready_fullomni.csv"))

meta_network_carnival_ready$source <- paste("X",meta_network_carnival_ready$source, sep = "")
meta_network_carnival_ready$target <- paste("X",meta_network_carnival_ready$target, sep = "")

meta_network_carnival_ready$source <- gsub("[-+{},;() ]","______",meta_network_carnival_ready$source)
meta_network_carnival_ready$target <- gsub("[-+{},;() ]","______",meta_network_carnival_ready$target)

metab_enzyme_reactions <- unique(c(
  meta_network_carnival_ready$source,
  meta_network_carnival_ready$target))

metab_enzyme_reactions <- unique(metab_enzyme_reactions[grepl("XGene.*", metab_enzyme_reactions)])

exchanges <- rep(F,length(metab_enzyme_reactions))

e <- 1
k <- 1
df_exchanges <- list()
for(metab_enzyme_reaction in metab_enzyme_reactions) #in this loop, we separate exchange reaction into distinct path to ensure different metabolites cannot be transformed into one another through exchanges
{
  # print(e)
  df <- meta_network_carnival_ready[
    meta_network_carnival_ready$source == metab_enzyme_reaction |
      meta_network_carnival_ready$target == metab_enzyme_reaction,
    ]
  df_inverse <- df[,c(3,2,1)]
  if(metab_enzyme_reaction == "XGene5599__3767_6833")
  {
    print(e)
  }
  for(i in 1:length(df[,1]))
  {
    for(j in i:length(df_inverse[,1]))
    {
      if(sum(gsub("___[a-z]____","",df[i,]) == gsub("___[a-z]____","",df_inverse[j,])) == 3)
      {
        df[i,grep("XGene.*",df[i,])] <- paste(df[i,grep("XGene.*",df[i,])], paste("EXCHANGE",i,sep = ""), sep = "")
        df[j,grep("XGene.*",df[j,])] <- paste(df[j,grep("XGene.*",df[j,])], paste("EXCHANGE",i,sep = ""), sep = "")
        
        if(sum(grepl("X[0-9]+$",df[,1])) != 0)
        {
          new_row <- df[grepl("X[0-9]+$",df[,1]) & !grepl("EXCHANGE",df[,3]),]
          if(length(new_row[,1]) == 1){
            new_row[3] <- paste(new_row[3], paste("EXCHANGE",i,sep = ""),sep = "")
            df <- as.data.frame(rbind(df,new_row))
          } else
          {
            for(m in 1:length(new_row[,1]))
            {
              new_row[m,3] <- paste(new_row[m,3], paste("EXCHANGE",i,sep = ""),sep = "")
              df <- as.data.frame(rbind(df,new_row))
            }
          }
          
        } else
        {
          if(sum(grepl("X[0-9]+$",df[,3])) != 0)
          {
            new_row <- df[grepl("X[0-9]+$",df[,3]) & !grepl("EXCHANGE",df[,1]),]
            if(length(new_row[,1]) == 1){
              new_row[3] <- paste(new_row[1], paste("EXCHANGE",i,sep = ""),sep = "")
              df <- as.data.frame(rbind(df,new_row))
            } else
            {
              for(m in 1:length(new_row[,1]))
              {
                new_row[m,1] <- paste(new_row[m,1], paste("EXCHANGE",i,sep = ""),sep = "")
                df <- as.data.frame(rbind(df,new_row))
              }
            }
          }
        }
        # print(df)
        exchanges[e] <- TRUE
      }
    }
  }
  if(exchanges[e] & sum(grepl("X[0-9]+$",df[,1]) | grepl("X[0-9]+$",df[,3])) != 0)
  {
    df <- df[unique(c(grep("EXCHANGE[0-9]+$",df[,1]), grep("EXCHANGE[0-9]+$",df[,3]))),]
  }
  if(exchanges[e])
  {
    df_exchanges[[k]] <- df
    k <- k+1
  }
  e <- e +1
}

df_exchanges <- as.data.frame(do.call(rbind,df_exchanges))

exchange_reactions <- metab_enzyme_reactions[exchanges]

meta_network_carnival_ready <- meta_network_carnival_ready[!(
  meta_network_carnival_ready$source %in% exchange_reactions |
    meta_network_carnival_ready$target %in% exchange_reactions),
  ]

meta_network_carnival_ready <- as.data.frame(rbind(meta_network_carnival_ready, df_exchanges))
row.names(meta_network_carnival_ready) <- c(1:length(meta_network_carnival_ready[,1]))


meta_network_carnival_ready <- unique(meta_network_carnival_ready)
write_csv(meta_network_carnival_ready, "~/Dropbox/Meta_PKN/result/meta_network_carnival_ready_exch_solved_fullomni.csv")
write_tsv(meta_network_carnival_ready, "~/Dropbox/Meta_PKN/result/meta_network_carnival_ready_exch_solved_fullomni.tsv")


