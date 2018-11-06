library('dplyr')
source('helpers.R')
registerDoParallel(10)


# Read in simulation results and prepare them for plotting
##########################################################

res <- readRDS('Simulation-Results/simres-mean-0.0-5-80.RDS')
res <- readRDS('Simulation-Results/simres-mean-0.3-5-50.RDS')

dcorr <- get_measure(res, compute_corr, 'lace')
dclass <- get_measure(res,
                       function(simdat, p, measure) compute_class(simdat, p, topk = 1, measure),
                       measure = 'ce')

dcorr <- readRDS('Simulation-Results/corr-mean-0.0-5-80.RDS')
dclass <- readRDS('Simulation-Results/class-mean-0.0-5-80.RDS')

dcorr <- readRDS('Simulation-Results/corr-mean-0.3-5-50.RDS')
dclass <- readRDS('Simulation-Results/class-mean-0.3-5-50.RDS')

remove <- paste0('Nodes: ', c(5, 40, 60, 70))
remove <- paste0('Nodes: ', c(5))

dcorr <- dcorr %>% 
  filter(
    !(nodes %in% remove)
    # measure != 'ACE'
  ) %>% 
  mutate(
    nodes = factor(nodes),
    measure = factor(measure, levels(measure)[c(seq(2, 6), 1)]) # put ACE back
  )

dclass <- dclass %>% 
  filter(
    !(nodes %in% remove),
    measure != 'ACE'
  ) %>% 
  mutate(
    nodes = factor(nodes),
    measure = factor(measure, levels = levels(dcorr$measure))
  )

dclass_mean <- dclass %>% 
  group_by(nodes, conn, type, measure) %>% 
  summarize(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE)
  )

# recovers Figures 5 & 6 (or 7 & 8)
plot_mean_results(dclass_mean, 'Probability of Successful Classification', class = TRUE)
plot_results(dcorr, 'Correlation between Causal Effect and Centrality Measures', class = FALSE)


# recovers Figure 3 in the paper
set.seed(7)
p <- 10
g <- randDAG(p, conn2degree(p, .3), method = 'er', wFUN = list(rnorm, mean = .1, sd = .5))
dat <- rDAG(1000, g)
mrf <- get_mrf(dat)

p <- plot_single(g, mrf, 'Betweenness')
compplot(g, mrf, 'Betweenness')
