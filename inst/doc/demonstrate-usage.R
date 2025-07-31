## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align ='center', 
  fig.width = 7, 
  cache = FALSE
)

## ----setup, echo = FALSE------------------------------------------------------
library(Pv3Rs)

## -----------------------------------------------------------------------------
y <- list("Enrolment" = list(m1 = c('b','c','d'), 
                              m2 = c('a','b'),
                              m3 = c('b','c','d')),
          "Recurrence 1" = list(m1 = c('b','d'), 
                                m2 = c('a'), 
                                m3 = c('a','b')),
          "Recurrence 2" = list(m1 = c('d'), 
                                m2 = c('a'), 
                                m3 = c('a')))

fs <- list(m1 = c(a = 0.27, b = 0.35, c = 0.18, d = 0.20),
           m2 = c(a = 0.78, b = 0.14, c = 0.07, d = 0.01),
           m3 = c(a = 0.21, b = 0.45, c = 0.26, d = 0.08))

## ----fig.align='center', fig.height=6-----------------------------------------
plot_data(ys = list("Participant data" = y), fs = fs)

## ----cache=TRUE---------------------------------------------------------------
post <- compute_posterior(y, fs)

## -----------------------------------------------------------------------------
sort(post$joint, decreasing = T)

## -----------------------------------------------------------------------------
post$marg

## -----------------------------------------------------------------------------
post$joint["LC"] + post$joint["LL"] + post$joint["LI"]

## -----------------------------------------------------------------------------
oldpar <- par(no.readonly = TRUE)
par(mar = c(0,0,0,0))
plot_simplex(p.coords = rbind(post$marg, Prior = rep(1/3, 3)), pch = 20)
par(oldpar)

## -----------------------------------------------------------------------------
post <- compute_posterior(y, fs, return.RG = TRUE, return.logp = TRUE)

## -----------------------------------------------------------------------------
sort(post$joint, decreasing = T)

## -----------------------------------------------------------------------------


# Extract all log likelihoods
llikes <- sapply(post$RGs, function(RG) RG$logp)

# Get maximum log likelihood
mllikes <- max(llikes)

# Extract the relationship graphs (RGs) with the maximum log likelihood
RGs <- post$RGs[which(abs(llikes - mllikes) < .Machine$double.eps^0.5)]

# Plot RGs with maximum log likelihoods
oldpar <- par(no.readonly = TRUE) # Store user's options
par(mar = rep(0.1,4), mfrow = c(1,2))
for(i in 1:length(RGs)) {
  plot_RG(RG_to_igraph(RGs[[i]], determine_MOIs(y)), vertex.size = 20)
  box()
}

# Add a legend
legend("bottomright", pch = 21, 
       pt.bg =  RColorBrewer::brewer.pal(n = 8, "Set2") [1:length(y)], 
       bty = "n", legend = names(y), title = "Episode")

par(oldpar) # Restore user's options

## -----------------------------------------------------------------------------
# In the following code, we place two graphs in the same equivalence class if
# they share the same likelihood. This is not ideal (two graphs that are not
# isomorphic up to permutation could share the same likelihood), but it works
# here: the plot shows only isomorphic graphs within the equivalence class.

sorted_llikes <- sort(llikes, decreasing = T) # Sort log likelihoods
adj_equal <- abs(diff(sorted_llikes, lag = 1)) < .Machine$double.eps^0.5 # Find matches
decr_idxs <- which(adj_equal == FALSE) # Change points: 2, 8, 14, 20, 32, ...
class_sizes <- c(decr_idxs[1], diff(decr_idxs)) # Number of graphs per class

# log likelihood of representative from each 'equivalence class' (EC)
llikes_unique <- sorted_llikes[decr_idxs]

# EC likelihood
class_ps <- exp(llikes_unique)*class_sizes
max_class_p <- which(class_ps == max(class_ps)) # ML EC index 
max_idx <- decr_idxs[max_class_p] # Index of last graph in ML EC
max_size <- class_sizes[max_class_p] # Number of graphs in ML EC

# Plot all graphs within the ML EC 
oldpar <- par(no.readonly = TRUE) # Store user's options
par(mar = rep(0.1,4), mfrow = c(3,4))
RG_order <- order(llikes, decreasing = T) # order RGs by logl
for(i in (max_idx-max_size+1):max_idx) { # EC consists of the RGs with logl rank 21-32
  RG <- post$RGs[[RG_order[i]]]
  RG_igraph <- RG_to_igraph(RG, determine_MOIs(y))
  plot_RG(RG_igraph, vertex.size = 25, vertex.label = NA)
  box()
}
par(oldpar) # Restore user's options

