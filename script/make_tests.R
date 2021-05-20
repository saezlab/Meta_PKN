library(cosmos)
library(readr)

RNA_ttop_tumorvshealthy <- as.data.frame(
  read_csv("Dropbox/COSMOS_MSB/data/RNA_ttop_tumorvshealthy.csv"))

toy_RNA <- RNA_ttop_tumorvshealthy$t
names(toy_RNA) <- paste("X",RNA_ttop_tumorvshealthy$ID, sep ="")

toy_network <- toy_network
toy_signaling_input_carnival_vec <- toy_signaling_input_carnival_vec
toy_metab_input_carnival_vec <- toy_metab_input_carnival_vec

library(igraph)

gnet <- graph_from_data_frame(toy_network[,c(1,3,2)])

shortest_paths(gnet, from = "XMetab__190___c____",
               to = "X6733")

p1 <- shortest_paths(gnet, from = "X3725", 
               to = "XMetab__6426851___c____")

p2 <- shortest_paths(gnet, from = "X2305", 
               to = "XMetab__124886___c____")

p3 <- shortest_paths(gnet, from = "XMetab__60961___c____", 
               to = "X3725")

p4 <- shortest_paths(gnet, from = "XMetab__190___c____",
               to = "X6733")

#for short
nodes <- unique(c(names(unlist(p1$vpath)),
                  names(unlist(p2$vpath)),
                  names(unlist(p3$vpath)),
                  names(unlist(p4$vpath))))

#for super_short
# nodes <- unique(c(names(unlist(p1$vpath))))

toy_network_lpSolve <- toy_network[toy_network$source %in% nodes & toy_network$target %in% nodes,]

toy_signaling_input_carnival_vec_lpSolve <- toy_signaling_input_carnival_vec[names(toy_signaling_input_carnival_vec) %in% nodes]
toy_metab_input_carnival_vec_lpSove <- toy_metab_input_carnival_vec[names(toy_metab_input_carnival_vec) %in% nodes]

toy_network_lpSolve$interaction <- abs(toy_network_lpSolve$interaction)

toy_signaling_input_carnival_vec_lpSolve <- abs(toy_signaling_input_carnival_vec_lpSolve)
toy_metab_input_carnival_vec_lpSove <- abs(toy_metab_input_carnival_vec_lpSove)

new_toys <- list("network" = toy_network_lpSolve,
                 "metab" = toy_metab_input_carnival_vec_lpSove,
                 "sig" = toy_signaling_input_carnival_vec_lpSolve,
                 "RNA" = toy_RNA)

save(new_toys, file = "Desktop/new_toys.RData")
