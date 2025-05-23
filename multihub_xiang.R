############################################################################
#How does the presence of another major, connected hub prevent the cascade# 
######### from becoming total, compared to a single BA graph? ##############
############################################################################

library(igraph)

cascade.Xiang.multi = function(g, seed, beta = 1, T_cap = 1.2, alpha = 1, max.steps = 10000) {
  N = vcount(g)
  deg = igraph::degree(g)
  load0 = beta * (deg^alpha)
  cap = T_cap * load0
  
  load = load0
  failed = logical(N)
  failed[seed] = TRUE
  
  for(s in seed){
    
    nbrs = as.integer(neighbors(g, s))
    active_nbrs = nbrs[!failed[nbrs]] # neighbors that are not themselves initial targets
    
    if (length(active_nbrs) > 0) {
      # use a small epsilon for degrees if they are 0 to avoid 0 weights if alpha > 0
      # or sum(wts) = 0 if all active_nbrs have deg=0
      deg_active_nbrs = deg[active_nbrs]
      wts = (deg_active_nbrs + ifelse(deg_active_nbrs==0, 1e-9, 0))^alpha
      
      if(sum(wts) > 0){
        load[active_nbrs] = load[active_nbrs] + load0[s] * (wts / sum(wts))
      } else if (length(active_nbrs) > 0) { # if sum(wts) is 0 but there are neighbors (e.g. all deg=0)
        load[active_nbrs] = load[active_nbrs] + load0[s] / length(active_nbrs) # Distribute evenly
      }
    }
    load[s] = 0
  }
  
  series = c(sum(failed) / N) 
  step = 0
  
  while (step < max.steps) {
    step = step + 1
    
    # identify nodes that will fail in this step
    to_fail = which(!failed & (load > cap))
    if (!length(to_fail)) break
    
    failed_nodes = c() # track nodes that actually fail in this iteration
    
    for (v in to_fail) {
      if(failed[v]) next # should not happen due to `which` condition, but safe check
      
      failed[v] = TRUE
      failed_nodes = c(failed_nodes, v)
      
      load_to_redistribute = load[v] # this is the load that causes failure
      load[v] = 0 
      
      nbrs2 = as.integer(neighbors(g, v))
      # redistribute only to neighbors that are not already marked as failed
      active_nbrs2 = nbrs2[!failed[nbrs2]]
      
      if (length(active_nbrs2) > 0) {
        deg_active_nbrs2 = deg[active_nbrs2]
        w2 = (deg_active_nbrs2 + ifelse(deg_active_nbrs2==0, 1e-9, 0))^alpha
        
        if(sum(w2) > 0){
          load[active_nbrs2] = load[active_nbrs2] + load_to_redistribute * (w2 / sum(w2))
        } else if (length(active_nbrs2) > 0) {
          load[active_nbrs2] = load[active_nbrs2] + load_to_redistribute / length(active_nbrs2)
        }
      }
    }
    series = c(series, sum(failed) / N)
  }
  return(series)
  
}

set.seed(11123)
n1 = 50 
n2 = 50 
m_param = 2 
num_additional_inter_edges = 0 # number of EXTRA edges to connect hub neighborhoods

# 1. Create two BA graphs
g1 = sample_pa(n=n1, power=1, m=m_param, directed=FALSE)
V(g1)$name = paste0("g1v", 1:vcount(g1)) 

g2 = sample_pa(n=n2, power=1, m=m_param, directed=FALSE)
V(g2)$name = paste0("g2v", 1:vcount(g2))

# 2. Find the hub (highest degree node) in each graph
hub1 = which.max(degree(g1))
hub2 = which.max(degree(g2))

# Get neighbors of these hubs in their original graphs
nbrs_hub1 = neighbors(g1, hub1)
nbrs_hub2 = neighbors(g2, hub2)

# 3. Combine the graphs
combined_g = disjoint_union(g1, g2)

# 4. Identify hubs and their neighbors in the combined graph
hub1_combined = hub1
hub2_combined = n1 + hub2

# Map original neighbor indices to combined graph indices
neighbors_hub1_combined = as.integer(nbrs_hub1) # these are already correct as they are from g1
neighbors_hub2_combined = as.integer(nbrs_hub2) + n1 # offset for g2 nodes

# 5. Add edges to connect the hub regions
edges_to_add = c()

# 5a. Direct hub-to-hub connection
edges_to_add = c(edges_to_add, hub1_combined, hub2_combined)

# 5b. Additional sparse connections between hub neighborhoods
if (length(neighbors_hub1_combined) > 0 && length(neighbors_hub2_combined) > 0 && num_additional_inter_edges > 0) {
  
  # Ensure we don't try to sample more unique pairs than available
  max_edges = length(neighbors_hub1_combined) * length(neighbors_hub2_combined)
  edges_to_attempt = min(num_additional_inter_edges, max_edges)
  
  added_count = 0
  attempt_limit = edges_to_attempt * 5 # try a few times to get unique edges
  current_attempts = 0
  
  # store existing edges to check for duplicates
  existing_edges_mat = as_edgelist(combined_g, names=FALSE)
  if(length(edges_to_add) > 0) { # include the hub-hub edge if already defined
    existing_edges_mat = rbind(existing_edges_mat, matrix(edges_to_add, ncol=2, byrow=TRUE))
  }
  
  while(added_count < edges_to_attempt && current_attempts < attempt_limit){
    current_attempts = current_attempts + 1
    # randomly pick one neighbor from each hub's neighborhood
    n1_sample = sample(neighbors_hub1_combined, 1)
    n2_sample = sample(neighbors_hub2_combined, 1)
    
    new_edge = sort(c(n1_sample, n2_sample)) # sort to handle undirected nature for duplicate check
    
    # check if this edge (or its reverse) already exists or is planned
    is_duplicate = FALSE
    if (nrow(existing_edges_mat) > 0) {
      for(k in 1:nrow(existing_edges_mat)){
        if(all(sort(existing_edges_mat[k,]) == new_edge)){
          is_duplicate = TRUE
          break
        }
      }
    }
    # check against already selected additional edges for this loop
    if (!is_duplicate && length(edges_to_add) > 0) {
      temp_added_edges = matrix(edges_to_add, ncol=2, byrow=TRUE)
      for(k in 1:nrow(temp_added_edges)){
        if(all(sort(temp_added_edges[k,]) == new_edge)){
          is_duplicate = TRUE
          break
        }
      }
    }
    
    if(!is_duplicate){
      edges_to_add = c(edges_to_add, n1_sample, n2_sample)
      added_count = added_count + 1
      cat("Adding inter-neighborhood edge between", n1_sample, "and", n2_sample, "\n")
    }
  }
  
}

# create the final connected graph
combined_g_connected_multi = add_edges(combined_g, edges_to_add)
V(combined_g_connected_multi)$id = 1:vcount(combined_g_connected_multi)

# 6. Create a single BA graph of comparable size for comparison (as before)
n_total = n1 + n2
g_single_ba = sample_pa(n=n_total, power=1, m=m_param, directed=FALSE)
V(g_single_ba)$id = 1:vcount(g_single_ba)

# 7. Run Cascades and Compare
# fail one of the original hubs in the multi-connected graph
seed_multi_connected = hub1
cascade_multi_connected_xiang = cascade.Xiang.multi(
  g = combined_g_connected_multi,
  seed = seed_multi_connected,
  T_cap = 1.2, alpha = 1.5)
# Cascade in the single BA graph, failing its main hub
hub_single_ba = which.max(degree(g_single_ba))
cascade_single_ba_xiang = cascade.Xiang.multi(
  g = g_single_ba,
  seed = hub_single_ba,
  T_cap = 1.2, alpha = 1.5)


max_len = max(length(cascade_multi_connected_xiang), 
              length(cascade_single_ba_xiang), n_total)


if(length(cascade_multi_connected_xiang) < max_len) {
  padded_multi_connected_xiang = c(cascade_multi_connected_xiang, 
                                   rep(tail(cascade_multi_connected_xiang,1), 
                                       max_len - length(cascade_multi_connected_xiang)))
} else {
  padded_multi_connected_xiang = cascade_multi_connected_xiang[1:max_len]
}

padded_single_ba_xiang = cascade_single_ba_xiang
if(length(cascade_single_ba_xiang) < max_len) {
  padded_single_ba_xiang = c(cascade_single_ba_xiang, rep(tail(cascade_single_ba_xiang,1), 
                                                          max_len - length(padded_single_ba_xiang)))
} else {
  padded_single_ba_xiang = cascade_single_ba_xiang[1:max_len]
}


plot(padded_multi_connected_xiang, type="o", col="blue", pch=1, ylim=c(0,1),
     xlab="Step", ylab="Proportion Failed", main="Multi Hub Connection (Xiang)")
lines(padded_single_ba_xiang, type="o", col="red", pch=2)
legend("bottomright", 
       legend=c("Multi-Connected BA",
                "Single BA"),
       col=c("blue", "red"), pch=c(1,2), lty=1, bty="n")
