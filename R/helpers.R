library('ggm')
library('pcalg')
library('ggpubr')
library('qgraph')
library('igraph')
library('ggplot2')
library('mvtnorm')
library('reshape2')
library('doParallel')
# library('JuliaCall')


#' Estimates a (Gaussian) Markov Random Field using glasso
#'
#' @param dat data matrix
#' @returns a matrix specifying the MRF
get_mrf <- function(dat) {
  return(suppressWarnings(EBICglasso(cov(dat), n = nrow(dat))))
}


#' Estimates a (Gaussian) Markov Random Field using CI tests
#'
#' @param dat data matrix
#' @param alpha test level
#' @returns a matrix specifying the MRF
get_mrf2 <- function(dat, alpha = .01) {
  n <- nrow(dat)
  p <- ncol(dat)
  G <- diag(p)
  S <- var(dat)
  
  for (i in seq(p)) {
    for (j in seq(i, p)) {
      if (i != j) {
        r <- pcor(c(i, j), S)
        test <- pcor.test(r, p - 2, n)
        G[i, j] <- G[j, i] <- ifelse(test$pvalue < alpha, r, 0)
      }
    }
  }
  
  G
}


#' Topologically sorts a DAG generated from pcalg::randDAG
#'
#' @param g a pcalg DAG object
#' @returns a topologically sorted DAG
topsort <- function(g) {
  amat <- wgtMatrix(g)
  toporder <- topo_sort(igraph.from.graphNEL(g), mode = 'in')
  amat_ord <- amat[toporder, toporder]
  igraph.to.graphNEL(
    graph_from_adjacency_matrix(amat_ord, mode = 'directed', weighted = TRUE)
  )
}


#' Computes the KL-based Causal Effect for each all in a DAG
#'
#' @param dag a pcalg DAG object
#' @param error the variance of the error
#' @returns a numeric vector indicating the KL-based Causal Effect for each node
get_causal_effect <- function(dag, error = 1) {
  A <- t(wgtMatrix(dag))
  nnodes <- colnames(A)
  p <- length(nnodes)
  res <- numeric(p)
  
  I <- diag(1, p)
  E <- diag(error, p)
  
  for (i in seq(p)) {
    node <- nnodes[i]
    D <- qr.solve(I - A)
    S <- D %*% E %*% t(D)
    
    As <- A
    As[node == nnodes, ] <- 0
    Dcut <- qr.solve(I - As)
    
    Ecut <- E + As %*% diag(diag(S), p) %*% t(As)
    Scut <- Dcut %*% Ecut %*% t(Dcut)
    
    # change to Julia on the server
    # julia_assign('Scut', Scut)
    # Scut.inv <- julia_eval('inv(Scut)')
    
    Scut.inv <- qr.solve(Scut)
    KL <- .5 * (sum(diag(Scut.inv %*% S)) - (determinant(S)$modulus - determinant(Scut)$modulus) - p)
    res[as.numeric(node)] <- KL
  }
  
  res
}


#' Computes the ACE-based Causal Effect for all nodes in a DAG
#' (We do not use this in the paper, but compute ACE only locally)
#'
#' @param g a pcalg DAG object
#' @returns a numeric vector indicating the ACE-based Causal Effect for each node
get_total_effect <- function(g) {
  ig <- igraph.from.graphNEL(g)
  wmat <- t(wgtMatrix(g))
  p <- ncol(wmat)
  res <- numeric(p)
  
  for (i in seq(p)) {
    total_effect <- 0
    for (j in seq(p)) {
      effect_on_j <- 0
      paths <- all_simple_paths(ig, from = i, to = j, mode = 'out')
      
      n <- length(paths)
      if (n != 0) {
        for (k in seq(n)) {
          path <- paths[[k]]
          
          w <- 1
          l <- 1
          while (l < length(path)) {
            w <- w*wmat[path[l], path[l + 1]]
            l <- l + 1
          }
          effect_on_j <- effect_on_j + w
        }
      }
      total_effect <- total_effect + abs(effect_on_j)
    }
    res[i] <- total_effect
  }
  res
}


#' Computes eigenvector centrality
#' 
#' @param G a precision matrix
#' @returns a numeric vector of eigenvector centralities for each node
eigencentrality <- function(G) {
  ig <- graph_from_adjacency_matrix(G, weighted = TRUE, mode = 'undirected')
  eigen_centrality(ig)$vector
}


#' Computes the local causal effect based on the true DAG
#' 
#' @param g is a true DAG
#' @returns a numeric vector of local causal effects
get_local_causal_effect <- function(g) {
  wMat <- wgtMatrix(g)
  apply(wMat, 2, function(x) sum(abs(x)))
}


#' Converts node connectivity to average node degree
#'
#' @param nodes the number of nodes
#' @param conn the network connectivity/density
#' @returns a number indicating the average degree of each node
conn2degree <- function(nodes, conn) {
  conn*(nodes - 1)
}


#' Samples from a DAG
#'
#' This is the same as the pcalg::rmvDAG function, except that I fixed a (terrible) bug
#' I also wrote dag@nodes instead of nodes(dag), as nodes(.) is not exported from pcalg
rmvDAG_modified <- function (n, dag, errDist = c("normal", "cauchy", "t4", "mix", 
                              "mixt3", "mixN100"), mix = 0.1, errMat = NULL, back.compatible = FALSE, 
          use.node.names = !back.compatible) 
{
  stopifnot(is(dag, "graph"), (p <- length(dag@nodes)) >= 
              2)
  weightMatrix <- if (back.compatible) 
    wgtMatrix.0(dag)
  else wgtMatrix(dag)
  nonZeros <- which(weightMatrix != 0, arr.ind = TRUE)
  if (nrow(nonZeros) > 0) {
    if (any(nonZeros[, 1] - nonZeros[, 2] < 0) || any(diag(weightMatrix) != 
                                                      0)) 
      stop("Input DAG must be topologically ordered!")
  }
  errDist <- match.arg(errDist)
  if (grepl("^mix", errDist)) 
    eMat <- function(outs) {
      X <- c(rnorm(n * p - length(outs)), outs)
      matrix(sample(X), nrow = n)
    }
  if (is.null(errMat)) {
    errMat <- switch(errDist,
                     normal = matrix(rnorm(n * p), nrow = n),
                     cauchy = matrix(rcauchy(n * p), nrow = n), 
                     t4 = matrix(rt(n * p, df = 4), nrow = n),
                     mix = eMat(rcauchy(round(mix * n * p))),
                     mixt3 = eMat(rt(round(mix * n * p), df = 3)),
                     mixN100 = eMat(rnorm(round(mix * n * p), sd = 10)))
  }
  else {
    stopifnot(!is.null(dim.eM <- dim(errMat)), dim.eM == c(n, p), is.numeric(errMat))
  }
  if (use.node.names) 
    colnames(errMat) <- dag@nodes
  if (sum(abs(weightMatrix)) > 0) { # bug fixed! (see https://github.com/cran/pcalg/pull/2)
    X <- errMat
    for (j in 2:p) {
      ij <- 1:(j - 1)
      X[, j] <- X[, j] + X[, ij, drop = FALSE] %*% weightMatrix[j, ij]
    }
    X
  }
  else errMat
}


#' Generates data from a DAG
#'
#' @param n the number of data points to be generated
#' @param dag a pcalg DAG object
#' @returns a p x n data matrix where p is the number of nodes; error variance is 1
rDAG <- function(n, dag) {
  tsdag <- topsort(dag)
  dat <- rmvDAG_modified(n, tsdag, errDist = 'normal', mix = 0)
  dat[, order(as.numeric(colnames(dat)))]
}


#' Runs one iteration of the simulation study
#'
#' @param nobs the number of observations to generate
#' @param nodes the number of nodes
#' @param conn the network connectivity
#' @param method the type of network to generate
#' @param times the number of simulation iterations
#' @param ... more arguments to pcalg::randDAG
#' @returns a p x n data matrix where p is the number of nodes; error variance is 1
simulation_study <- function(
  nobs = 500, nodes = 5, conn = .5,
  method = 'er', times = 100, ...
  ) {
  
  degree <- conn2degree(nodes, conn)
  
  measures <- c(
    'ce', 'lace', 'strength', 'closeness', 'betweenness', 'degree', 'eigen'
  )
  
  m <- length(measures)
  res <- matrix(NA, ncol = nodes + 1, nrow = times * m)
  colnames(res) <- c(as.character(seq(nodes)), 'measure')
  get_abs_sum <- function(x) apply(x, 1, function(row) sum(abs(row)))
  
  k <- 1
  for (i in seq(times)) {
    
    ERROR <- FALSE
    dag <- tryCatch({
      randDAG(n = nodes, d = degree, method = method, wFUN = list(rnorm, mean = 0, sd = .5), ...)
    }, error = function(e) {
      ERROR <<- TRUE
    })
    
    if (!ERROR) {
      dat <- rDAG(n = nobs, dag)
      
      ace <- get_local_causal_effect(dag)
      ce <- get_causal_effect(dag, error = 1)
      
      mrf <- get_mrf(dat, alpha = alpha)
      central <- centrality(mrf)
      strength <- get_abs_sum(mrf)
      eigencent <- eigencentrality(mrf)
      
      r <- rbind(
        ce,
        ace,
        strength,
        central$Closeness,
        central$Betweenness,
        central$InDegree,
        eigencent
      )
      
      r <- cbind(r, measures)
    
    } else {
      na <- rep(NA, m)
      r <- matrix(NA, nrow = m, ncol = nodes)
      r <- cbind(r, measures)
    }
    
    res[seq(k, k + m - 1), ] <- r
    k <- k + m
  }
  
  d <- data.frame(res)
  d[, seq(nodes)] <- apply(res[, seq(nodes)], 2, as.numeric)
  d$nodes <- nodes
  d$nobs <- nobs
  d$avg_degree <- degree
  d$type <- method
  d$conn <- conn
  
  d
}


#' Computes a particular measure for each simulation run
#'
#' @param res a matrix of lists with the simulation results
#' @param fun a function specifying the measure
#' @returns a list including the output of applying the measure function
#'          with information about the simulation run
get_measure <- function(res, fun, measure = 'ce') {
  dat <- foreach(i = seq(nrow(res)), .combine = rbind) %dopar% {
  # for (i in seq(nrow(res))) {
    simdat <- res[i, ][[1]]
    p <- simdat$nodes[1]
    nobs <- simdat$nobs[1]
    type <- simdat$type[1]
    avg_degree <- simdat$avg_degree[1]
    conn <- simdat$conn[1]
    
    cordat <- fun(simdat, p, measure)
    cordat$nodes <- p
    cordat$nobs <- nobs
    cordat$type <- type
    cordat$conn <- conn
    cordat
  }
  
  mapping <- list(
    'regular' = 'Regular',
    'watts' = 'Watts-Strogatz',
    'er' = 'Erdos-Renyi',
    'power' = 'Power-law',
    'bipartite' = 'Bipartite',
    'barabasi' = 'Barabasi',
    'geometric' = 'Geometric',
    'interEr' = 'Islands-Graph'
  )
  
  mapping2 <- list(
    'ce' = 'CE',
    'lace' = 'ACE',
    'ace' = 'ACE',
    'strength' = 'Strength',
    'closeness' = 'Closeness',
    'betweenness' = 'Betweenness',
    'degree' = 'Degree',
    'eigen' = 'Eigenvector'
  )
  
 d <- dat %>% 
   dplyr::mutate(
     type = unname(unlist(mapping[type])),
     type = factor(type),
     measure = as.character(measure),
     measure = unname(unlist(mapping2[measure])),
     measure = factor(measure),
     nodes = factor(nodes),
     conn = factor(conn)
   ) %>% as_tibble
 
 d$nodes <- factor(d$nodes,
                   labels = sapply(levels(d$nodes),
                                   function(i) paste0('Nodes: ', as.numeric(as.character(i)))))
 d
}


#' Computes Spearman's rho for a simulation run
#'
#' @param simdat a data.frame with the simulation results
#' @param p the number of nodes
#' @returns a list with mean and standard deviation of Spearman's rho
compute_corr <- function(simdat, p, measure = 'ce') {
  measures <- as.character(unique(simdat[, p + 1]))
  mapping <- list()
  
  for (i in seq(length(measures))) {
    mapping[[measures[i]]] <- i
  }
  
  dat <- simdat[, seq(p + 1)]
  n <- length(measures)
  times <- table(simdat$measure)[1]
  res <- array(dim = c(n, n, times))
  
  for (i in seq(n)) {
    for (j in seq(i, n)) {
      m1 <- dat[dat[, p + 1] == measures[i], -(p + 1)]
      m2 <- dat[dat[, p + 1] == measures[j], -(p + 1)]
      
      # compute correlations across simulation runs
      for (k in seq(times)) {
        if (i == j) {
          res[i, j, k] <- 1
        } else {
          res[i, j, k] <- res[j, i, k] <- cor(as.numeric(m1[k, ]), as.numeric(m2[k, ]), method = 'spearman')
        }
      }
    }
  }
  
  if (measure == 'ce') {
    m <- measures[-1]
  } else if (measure == 'betweenness') {
    corrs <- c(res[5, 1, ], res[5, 2, ], res[5, 3, ], res[5, 4, ], res[5, 6, ], res[5, 7, ])
    m <- measures[-5]
  } else if (measure == 'lace') {
    corrs <- c(res[2, 1, ], res[2, 3, ], res[2, 4, ], res[2, 5, ], res[2, 6, ], res[2, 7, ])
    m <- measures[-2]
  }
  
  data.frame(value = corrs, measure = rep(m, each = times))
}


#' Computes how successful the classification of the top k nodes is
#'
#' @param simdat a data.frame with the simulation results
#' @param p the number of nodes
#' @param topk the number of top nodes to consider
#' @returns a list with mean and standard deviation of Spearman's rho
compute_class <- function(simdat, p, topk = 1, measure = 'ce') {
  measures <- as.character(unique(simdat[, p + 1]))
  
  dat <- simdat[, seq(p + 1)]
  n <- length(measures)
  times <- table(simdat$measure)[1]
  res <- array(dim = c(n, n, times))
  
  get_pos <- function(r) sapply(seq(topk), function(i) which(r == i))
  mavek <- function(r1, r2) sum(which.min(r1) %in% get_pos(r2))
  
  for (i in seq(n)) {
    for (j in seq(i, n)) {
      m1 <- dat[dat[, p + 1] == measures[i], -(p + 1)]
      m2 <- dat[dat[, p + 1] == measures[j], -(p + 1)]
      
      for (k in seq(times)) {
        if (i == j) {
          res[i, j, k] <- 1
        } else {
          r1 <- order(m1[k, ], decreasing = TRUE)
          r2 <- order(m2[k, ], decreasing = TRUE)
          
          res[i, j, k] <- res[j, i, k] <- ifelse(
            any(is.na(m1[k, ])) || any(is.na(m2[k, ])), NA,
            mavek(r1, r2)
          )
        }
      }
    }
  }
  
  if (measure == 'ce') {
    corrs <- c(res[1, 2, ], res[1, 3, ], res[1, 4, ], res[1, 5, ], res[1, 6, ], res[1, 7, ])
    m <- measures[-1]
  } else if (measure == 'betweenness') {
    corrs <- c(res[5, 1, ], res[5, 2, ], res[5, 3, ], res[5, 4, ], res[5, 6, ], res[5, 7, ])
    m <- measures[-5]
  }
  
  data.frame(value = corrs, measure = rep(m, each = times))
}



#' Wrapper function to run the simulation
#'
#' @param nodes vector of the number of nodes
#' @param gtypes vector of the types of networks to simulate from
#' @param nobs vector of the number of observations to generate
#' @param connectivity vector of the connectivity/density levels
#' @param times number of times to run the simulation for each configuration
#' @returns a matrix of lists holding the simulation results
run_simulation_study <- function(nodes, gtypes, nobs = 1000, connectivity = .5, times = 500) {
  
  comb <- expand.grid(nodes = nodes, type = gtypes, nobs = nobs, conn = connectivity)
  ncomb <- nrow(comb)
  comb$type <- as.character(comb$type)
  
  res <- foreach(i = seq(ncomb), .combine = rbind) %dopar% {
  # for(i in seq(ncomb)) {
    nodes <- comb[i, 1]
    type <- comb[i, 2]
    nobs <- comb[i, 3]
    conn <- comb[i, 4]
    
    par1 <- NULL
    par2 <- NULL
    
    if (type == 'watts' || type == 'bipartite') { par1 <- .5 }
    if (type == 'barabasi') { par1 <- 1 }
    if (type == 'geometric') { par1 <- 2 }
    if (type == 'interEr') { par1 <- 2; par2 <- .25 }
    
    simdat <- simulation_study(
      nobs = nobs, nodes = nodes, conn = conn,
      method = type, times = times, par1 = par1, par2 = par2
    )
    
    list(
      'simdat' = simdat
    )
  }
  
  res
}


#' Z-standardises a numeric input vector
#' 
#' @param x a numeric vector
#' @returns a z-standardised numeric vector
std <- function(x) {
  if (sd(x) == 0) {
    return(rep(NA, length(x)))
  }
  (x - mean(x)) / sd(x)
}


#' Plots a DAG and MRF next to each other
#' 
#' @param dag a pcalg::randDAG object
#' @param mrf a matrix specifying the (Gaussian) MRF
#' @param network_measure the network_measure
#' @param ... additional arguments to qgraph
#' @returns NULL
compplot <- function(dag, mrf, network_measure = 'Betweenness', causal_measure = 'KL', ...) {
  nm <- ifelse(network_measure == 'Degree', 'InDegree', network_measure)
  
  get_abs_sum <- function(x) apply(x, 1, function(row) sum(abs(row)))
  
  if (causal_measure == 'KL') {
    ce <- std(get_causal_effect(dag))
  } else {
    ce <- std(get_local_causal_effect(dag))
  }
  
  if (nm == 'Strength') {
    cent <- std(get_abs_sum(mrf))
  } else if (nm == 'Eigenvector') {
    cent <- std(eigencentrality(mrf))
  } else {
    cent <- std(centrality(mrf)[[network_measure]])
  }
  
  par(mfrow = c(1, 2))
  qgraph(dag, vsize = (ce - min(ce) + 1)*2, layout = 'circle', asize = 12, aspect = TRUE, ...)
  qgraph(mrf, vsize = (cent - min(cent) + 1)*2, layout = 'circle', aspect = TRUE, ...)
}


#' Plots the distribution of node scores
#' 
#' @param dag a pcalg::randDAG object
#' @param mrf a matrix specifying the (Gaussian) MRF
#' @param network_measure the network_measure
#' @returns NULL
plot_single <- function(dag, mrf, network_measure = 'Betweenness', causal_measure = 'KL') {
  p <- length(dag@nodes)
  nm <- ifelse(network_measure == 'Degree', 'InDegree', network_measure)
  
  if (causal_measure == 'KL') {
    ce <- std(get_causal_effect(dag))
  } else {
    ce <- std(get_local_causal_effect(dag))
  }
  
  get_abs_sum <- function(x) apply(x, 1, function(row) sum(abs(row)))
  
  if (nm == 'Strength') {
    cent <- std(get_abs_sum(mrf))
  } else if (nm == 'Eigenvector') {
    cent <- std(eigencentrality(mrf))
  } else {
    cent <- std(centrality(mrf)[[network_measure]])
  }
  
  r <- round(cor(ce, cent, method = 'spearman'), 3)
  title <- paste0('Distribution of Node Scores (Rank Correlation: ', r, ')')
  network_measure <- ifelse(network_measure == 'InDegree', 'Degree', network_measure)
  
  dat <- data.frame(
    Measure = rep(c('Causal Effect', network_measure), each = p),
    Value = c(ce, cent),
    Node = rep(seq(p), 2)
  )
  
  dat$Node <- with(dat, reorder(Node, Value, function(x) -x[1]))
  
  graph <- ggplot(dat, aes(x = Node, y = Value, colour = Measure, group = Measure)) +
    geom_point() +
    geom_line() +
    scale_colour_discrete(name='') +
    scale_y_continuous(
      breaks = scales::pretty_breaks(n = 11)
      ) +
    xlab('\nNode') +
    ylab('z-Value\n') +
    ggtitle(title) +
    theme_pubclean() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.major.y = element_line(size = .1, color = 'black') ,
          plot.title = element_text(hjust = .5),
          text = element_text(size = 16))
  
  graph
}


#' Plots the (mean) simulation results
#'
#' @param dat is the prepared data, ready to be plotted
#' @param title is the title of the plot
#' @param n is the number of simulation repetitions per configuration
#' @param class is a boolean indicating whether to plot correlation or classification
#' @returns NULL
plot_mean_results <- function(dat, title, n = 500, class = FALSE) {
  nodes <- sapply(strsplit(as.character(dat$nodes), ':'), function(x) as.numeric(x[2]))
  dat$hline <- 1/nodes
  
  ytitle <- ifelse(class, 'Probability Correct\n', 'Rank Correlation\n')
  
  p <- ggplot(dat, aes(x = conn, y = mean, colour = measure, order = measure))
    
  if (class) {
    p <- p + geom_hline(aes(yintercept = hline))
  }
  
  p <- p +
    geom_point(position = position_dodge(.9), size = 1.5) +
    scale_x_discrete(labels = paste0('.', seq(1, 9, 1))) +
    xlab('\nConnectivity') +
    ylab(ytitle) + 
    facet_wrap(~ type + nodes, ncol = length(unique(dat$nodes))) +
    ggtitle(title) +
    theme_pubclean() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.major.y = element_line( size=.1, color="black") ,
          plot.title = element_text(hjust = .5),
          text = element_text(size = 16)) +
    scale_colour_manual(name = '', values = RColorBrewer::brewer.pal(n = 6, 'Set1'))
    
  if (class) {
    p <- p + 
      scale_y_continuous(breaks = scales::pretty_breaks(n = 11), limits = c(0, .35)) +
      geom_errorbar(aes(ymin = pmax(mean - 1.96 * sqrt(mean * (1 - mean) / n), 0),
                        ymax = pmin(mean + 1.96 * sqrt(mean * (1 - mean) / n), 1)),
                    width = .2, size = .5, position = position_dodge(.9))
    
  } else {
    p <- p + 
      scale_y_continuous(breaks = scales::pretty_breaks(n = 11), limits = c(-1, 1)) +
      geom_errorbar(aes(ymin = pmax(mean - 1.96 * sd/sqrt(n), -.5),
                        ymax = pmin(mean + 1.96 * sd/sqrt(n), 1)),
                    width = .2, size = .5, position = position_dodge(.9))
  }
  
  p
}


#' Plots the simulation results
#'
#' @param dat is the prepared data, ready to be plotted
#' @param title is the title of the plot
#' @param n is the number of simulation repetitions per configuration
#' @param class is a boolean indicating whether to plot correlation or classification
#' @returns NULL
plot_results <- function(dat, title, n = 500, class = FALSE) {
  nodes <- sapply(strsplit(as.character(dat$nodes), ':'), function(x) as.numeric(x[2]))
  dat$hline <- 1/nodes
  
  sd.box <- function(d) {
    data.frame(
      middle = mean(d),
      ymin = quantile(d, .1), ymax = quantile(d, .9),
      upper = quantile(d, .75), lower = quantile(d, .25)
    )
  }
  
  p <- ggplot(dat, aes(x = conn, y = value, fill = measure, order = measure)) +
    stat_summary(fun.data = sd.box, geom = 'boxplot', position = position_dodge(.9), size = .4)
    
  if (class) {
    p <- p + 
      scale_y_continuous(breaks = scales::pretty_breaks(n = 11), limits = c(0, .4)) +
      geom_hline(aes(yintercept = hline))
    
  } else {
    p <- p + 
      scale_y_continuous(breaks = scales::pretty_breaks(n = 11), limits = c(-1, 1))
  }
  
  ytitle <- ifelse(class, 'Probability Correct\n', 'Rank Correlation\n')
  p <- p + 
    scale_x_discrete(labels = paste0('.', seq(1, 9, 1))) +
    xlab('\nConnectivity') +
    ylab(ytitle) + 
    facet_wrap(~ type + nodes, ncol = length(unique(dat$nodes))) +
    # scale_colour_discrete(name = '') +
    ggtitle(title) +
    theme_pubclean() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.major.y = element_line(size = .1, color="black") ,
          plot.title = element_text(hjust = .5),
          text = element_text(size = 16)) +
    scale_fill_manual(name = '', values = RColorBrewer::brewer.pal(n = 6, 'Set1'))
  p
}




# Various functions to check if everything works right
######################################################
check <- function(n, d, ...) {
  g <- randDAG(n = n, d = d, wFUN = list(rnorm, mean = 0, sd = .5), ...)
  dat <- rDAG(2000, g)
  G <- get_mrf(dat)
  cent <- centrality(G)
  
  std <- function(x) {
    if (sd(x) == 0) {
      return(rep(0, length(x)))
    }
    (x - mean(x)) / sd(x)
  }
  
  degree <- cent$InDegree
  # g <- topsort(g)
  ce <- std(get_causal_effect(g))
  ace <- std(get_total_effect(g))
  
  
  S <- cbind(ce, ace, degree, betweenness = cent$Betweenness, closeness = cent$Closeness)
  list(
    G,
    (cor(S, method = 'spearman'))
  )
}

check2 <- function(p, conn, ...) {
  g <- randDAG(p, conn2degree(p, conn), wFUN = list(rnorm, mean = .3, sd = .5), ...)
  dat <- rDAG(1000, g)
  dat
  # mrf <- get_mrf(dat)
  # sum(mrf) == 0
}

check3 <- function(p, conn, ...) {
  g <- randDAG(p, conn2degree(p, conn), wFUN = list(rnorm, mean = 0.1, sd = .5), ...)
  d <- rDAG(1000, g)
  C <- cor(d)
  diag(C) <- 0
  as.numeric(C)
}