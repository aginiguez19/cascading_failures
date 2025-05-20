# Cascading Failure concept by Watts (2002)
# Intellectual property of Jonathan J. Park, PhD @UC Davis

## The logic of this failure system is that neighbors are tolerant to failures of their neighbors
## If enough of your neighbors fail, you will too
### Failure thresholds for each $i^{th}$ individual (\phi_i) are assumed to be random
### They are drawn from a uniform distribution from e.g., 0.10 - 0.20
### Thus, if 10- to 20-percent of neighbors fail, you will too

library(ggplot2)
library(tidyr)
library(dplyr)
library(igraph)

cascade.Watts = function(g, thrsh.lim = c(0.1, 0.20), max.steps = 10000){
  Cascades = list()
  metrics = igraph::degree(g)
  # Identify the high, average, and low degree vertices
  max_node = which.max(metrics)
  if(min(metrics, na.rm = TRUE) == 0){
    min_node = which.min(ifelse(metrics != 0 & !is.na(metrics), metrics, 999))
  }else(min_node = which.min(ifelse(!is.na(metrics), metrics, 999)))
  avg_node = sample(which(abs(metrics - floor(mean(metrics, na.rm = TRUE))) == min(abs(metrics - floor(mean(metrics, na.rm = TRUE))), na.rm = TRUE) & !is.na(metrics)), 1)
  seeds = list(high = max_node, 
               avg = avg_node,
               low = min_node)
  # For the high, middle, and low vertices initiate a cascade
  for(name in names(seeds)) {
    target = seeds[[name]]
    if(thrsh.lim[1] > thrsh.lim[2]){
      message("First limit should be the lower limit! Flipping for you")
      thrsh.lim = c(thrsh.lim[2], thrsh.lim[1])
    }
    # Uniformly apply thresholds to all vertices between limits
    V(g)$threshold = runif(vcount(g), min = thrsh.lim[1], max = thrsh.lim[2])
    V(g)$state = 0
    V(g)$state[target] = 1
    flip.count = c(sum(V(g)$state))
    step = 0
    # Initiate Cascade
    while(step < max.steps) {
      step = step + 1
      new.state = V(g)$state
      # For each vertex, check if it failed. If so, skip
      # If not, then identify its neighbors
      # Check how many of the current vertex's neighbors have failed relative to total
      # Compare to threshold; if bigger then fail
      # Repeat
      for (v in V(g)) {
        if (V(g)$state[v] == 1) next
        neighbors.v = neighbors(g, v)
        if (length(neighbors.v) == 0) next
          active_fraction = sum(V(g)$state[neighbors.v]) / length(neighbors.v)
        if (active_fraction >= V(g)$threshold[v]) {
          new.state[v] = 1
        }
      }
      if(all(new.state == V(g)$state)) break
      V(g)$state = new.state
      flip.count = c(flip.count, sum(V(g)$state))
    }
    Cascades[[name]] = flip.count / vcount(g)
  }
  names(Cascades) = c("High", "Average", "Low")
  return(Cascades)
}

# Initial Parameters for Simulation
n_nodes = 30

# Simulation Function
sim.funk = function(rep_id) {
  set.seed(3887891 + rep_id) 
  g_ba = sample_pa(n = n_nodes, power = 1, m = 1, directed = FALSE)
  g_er = sample_gnm(n = n_nodes, m = ecount(g_ba), directed = FALSE)
  g_sw = sample_smallworld(dim = 1, size = n_nodes, nei = floor(ecount(g_ba)/n_nodes), p = 0.05)
  ER.set = cascade.Watts(g_er, c(0.1, 0.20), 10000)
  BA.set = cascade.Watts(g_ba, c(0.1, 0.20), 10000)
  SW.set = cascade.Watts(g_sw, c(0.1, 0.20), 10000)
  list(
    ER_high = ER.set$High,
    ER_avg  = ER.set$Average,
    ER_low  = ER.set$Low,
    BA_high = BA.set$High,
    BA_avg  = BA.set$Average,
    BA_low  = BA.set$Low,
    SW_high = SW.set$High,
    SW_avg  = SW.set$Average,
    SW_low  = SW.set$Low
  )
}

results = sim.funk(1245)
models = c("ER","BA","SW")
max_len = 30
op = par(mfrow=c(3,1), mar=c(4,4,2,1))

for(m in models) {
  idx = grep(paste0("^", m, "_"), names(results))
  raw = results[idx]
  mat = sapply(raw, function(v) {
    L = length(v)
    if (L < max_len) {
      c(v, rep(v[L], max_len - L))
    } else {
      v[1:max_len]
    }
  })
  
  matplot(1:max_len, mat,
          type = "o", pch = 1:3,
          col = 1:3, xlab = "Step",
          ylab = "Prop. Failed",
          main = paste(m, "model"),
          xlim = c(1, max_len),
          ylim = c(0, 1))
  legend("bottomright", legend = names(results)[idx],
         col = 1:3, pch = 1:3, bty = "n")
}
par(op) 
 

### "Complex" Xiang Cascade 

cascade.Xiang = function(g, beta = 1, T = 1.2, alpha = 1, max.steps = 10000) {
  N = vcount(g)
  deg = igraph::degree(g)
  # Load is based on beta, alpha, and degree
  load0 = beta * (deg^alpha)
  cap = T * load0
  # Choose high, average, and lowest degree vertices
  high = which.max(deg)
  low = which.min(deg)
  avg = which.min(abs(deg - mean(deg)))[1]
  seeds = list(High=high, Average=avg, Low=low)
  results = lapply(seeds, function(seed) {
    load = load0
    failed = logical(N)
    # Kill the target
    # Find neighbors
    # Distribute failure to other neighbors
    # High degree neighbors take more load
    failed[seed] = TRUE
    nbrs = as.integer(neighbors(g,seed))
    if (length(nbrs)) {
      wts = deg[nbrs]^alpha
      load[nbrs] = load[nbrs] + load[seed]*(wts/sum(wts))
    }
    load[seed] = 0
    
    series = sum(failed)/N
    t = 1
    while (t < max.steps) {
      to_fail = which(!failed & (load > cap))
      if (!length(to_fail)) break
      for (v in to_fail) {
        failed[v] = TRUE
        nbrs2 = setdiff(as.integer(neighbors(g,v)), which(failed))
        if (length(nbrs2)) {
          w2 = deg[nbrs2]^alpha
          load[nbrs2] = load[nbrs2] + load[v]*(w2/sum(w2))
        }
        load[v] = 0
      }
      t = t+1
      series[t] = sum(failed)/N
    }
    series
  })
  
  names(results) = names(seeds)
  results
}


sim.funk = function(rep_id, n_nodes) {
  set.seed(15846 + rep_id) 
  g_ba = sample_pa(n = n_nodes, power = 1, m = 2, directed = FALSE)
  g_er = sample_gnm(n = n_nodes, m = ecount(g_ba), directed = FALSE)
  g_sw = sample_smallworld(dim = 1, size = n_nodes, nei = floor(ecount(g_ba)/n_nodes), p = 0.05)
  ER.set = cascade.Xiang(g_er, beta = 1, T = 1.2, alpha = 1, max.steps = 1000)
  BA.set = cascade.Xiang(g_ba, beta = 1, T = 1.2, alpha = 1, max.steps = 1000)
  SW.set = cascade.Xiang(g_sw, beta = 1, T = 1.2, alpha = 1, max.steps = 1000)
  list(
    ER_high = ER.set$High,
    ER_avg  = ER.set$Average,
    ER_low  = ER.set$Low,
    BA_high = BA.set$High,
    BA_avg  = BA.set$Average,
    BA_low  = BA.set$Low,
    SW_high = SW.set$High,
    SW_avg  = SW.set$Average,
    SW_low  = SW.set$Low
  )
}
n_nodes = 30
results = sim.funk(34, n_nodes = n_nodes)
models = c("ER","BA","SW")
max_len = 30
op = par(mfrow=c(3,1), mar=c(4,4,2,1))

for(m in models) {
  idx = grep(paste0("^", m, "_"), names(results))
  raw = results[idx]
  mat = sapply(raw, function(v) {
    L = length(v)
    if (L < max_len) {
      c(v, rep(v[L], max_len - L))
    } else {
      v[1:max_len]
    }
  })
  
  matplot(1:max_len, mat,
          type = "o", pch = 1:3,
          col = 1:3, xlab = "Step",
          ylab = "Prop. Failed",
          main = paste(m, "model"),
          xlim = c(1, max_len),
          ylim = c(0, 1))
  legend("bottomright", legend = names(results)[idx],
         col = 1:3, pch = 1:3, bty = "n")
}
par(op)


library(visNetwork)
g = barabasi.game(n = 20, power = 1)
V(g)$name = as.character(1:vcount(g))
nodes = data.frame(id = V(g)$name, label = V(g)$name)
edges = as_data_frame(g, what = "edges")
visNetwork(nodes, edges) %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visPhysics(stabilization = TRUE) %>%
  visInteraction(dragNodes = TRUE, dragView = TRUE)
