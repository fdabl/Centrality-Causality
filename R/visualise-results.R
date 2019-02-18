library('dplyr')
source('helpers.R')
registerDoParallel(10)


# Recovers Figure 3A and 3B from the paper
##########################################
set.seed(7)
p <- 10
g <- randDAG(p, conn2degree(p, .3), method = 'er', wFUN = list(rnorm, mean = .1, sd = .5))
dat <- rDAG(1000, g)
mrf <- get_mrf(dat)

p <- plot_single(g, mrf, 'Betweenness')
compplot(g, mrf, 'Betweenness')


##########################################################
# Read in simulation results and prepare them for plotting
# These are too large to store them on Github
# We therefore only store the already prepared data
##########################################################

# res <- readRDS('../Simulation-Results/simres-mean-0.0-5-80.RDS')
# res <- readRDS('../Simulation-Results/simres-mean-0.3-5-50.RDS')
# 
# dcorr_lace <- get_measure(res, compute_corr, 'lace')
# dclass_lace <- get_measure(res,
#                        function(simdat, p, measure) compute_class(simdat, p, topk = 1, measure),
#                        measure = 'lace')
# dranking_lace <- get_measure(res, compute_relative_ranking, measure = 'lace')

prepare_dat <- function(dat, type = 'correlation', correlate_with = 'CE', max_size = 50) {
  if (max_size == 50) {
    to_remove <- paste0('Nodes: ', c(5))
  } else {
    to_remove <- paste0('Nodes: ', c(5, 40, 60, 70))
  }
  
  if (correlate_with == 'CE') {
    reshuffle <- c(seq(2, 6), 1) # put ACE last
  } else {
    reshuffle <- c(1, 3, 4, 5, 6, 2) # put CE last
  }
  
  dat_ind <- dat %>% 
    filter(
      !(nodes %in% to_remove)
    ) %>% 
    mutate(
      nodes = factor(nodes),
      measure = factor(measure, levels(measure)[reshuffle]) # put (A)CE in last position
    )
  
  if (type != 'correlation') {
    dat_mean <- dat_ind %>% 
      group_by(nodes, conn, type, measure) %>% 
      summarize(
        mean = mean(value, na.rm = TRUE),
        sd = sd(value, na.rm = TRUE)
      )
    return(dat_mean)
  }
  
  dat_ind
}


## Data for correlations with CE Measure

# true mean 0.3
dcorr_CE <- prepare_dat(
  type = 'correlation', correlate_with = 'CE', max_size = 50,
  dat = readRDS('../Simulation-Results/CE/CE-corr-mean-0.3-5-50.RDS')
)
dclass_CE <- prepare_dat(
  type = 'classification', correlate_with = 'CE', max_size = 50,
  dat = readRDS('../Simulation-Results/CE/CE-class-mean-0.3-5-50.RDS')
)
dranking_CE <- prepare_dat(
  type = 'ranking', correlate_with = 'CE', max_size = 50,
  dat = readRDS('../Simulation-Results/CE/CE-ranking-mean-0.3-5-50.RDS')
)


# true mean 0.0
dcorr_CE0 <- prepare_dat(
  type = 'correlation', correlate_with = 'CE', max_size = 80,
  dat = readRDS('../Simulation-Results/CE/CE-corr-mean-0.0-5-80.RDS')
)
dclass_CE0 <- prepare_dat(
  type = 'classification', correlate_with = 'CE', max_size = 80,
  dat = readRDS('../Simulation-Results/CE/CE-class-mean-0.0-5-80.RDS')
)
dranking_CE0 <- prepare_dat(
  type = 'ranking', correlate_with = 'CE', max_size = 80,
  dat = readRDS('../Simulation-Results/CE/CE-ranking-mean-0.0-5-80.RDS')
)


## Data for correlations with ACE Measure [not reported in the paper]
# true mean 0.3
dcorr_ACE <- prepare_dat(
  type = 'correlation', correlate_with = 'ACE', max_size = 50,
  dat = readRDS('../Simulation-Results/ACE/ACE-corr-mean-0.3-5-50.RDS')
)
dclass_ACE <- prepare_dat(
  type = 'classification', correlate_with = 'ACE', max_size = 50,
  dat = readRDS('../Simulation-Results/ACE/ACE-class-mean-0.3-5-50.RDS')
)
dranking_ACE <- prepare_dat(
  type = 'ranking', correlate_with = 'ACE', max_size = 50,
  dat = readRDS('../Simulation-Results/ACE/ACE-ranking-mean-0.3-5-50.RDS')
)


# true mean 0.0
dcorr_ACE0 <- prepare_dat(
  type = 'correlation', correlate_with = 'ACE', max_size = 80,
  dat = readRDS('../Simulation-Results/ACE/ACE-corr-mean-0.0-5-80.RDS')
)
dclass_ACE0 <- prepare_dat(
  type = 'classification', correlate_with = 'ACE', max_size = 80,
  dat = readRDS('../Simulation-Results/ACE/ACE-class-mean-0.0-5-80.RDS')
)
dranking_ACE0 <- prepare_dat(
  type = 'ranking', correlate_with = 'ACE', max_size = 80,
  dat = readRDS('../Simulation-Results/ACE/ACE-ranking-mean-0.0-5-80.RDS')
)


## Correlation with CE Measure [Main Paper]
# Recovers Figures 5 and 6 (main text) and a relative ranking plot (not in paper)
##############################################################################
pdf('../Figures/Fig-CE-Corr-3.pdf', width = 10, height = 10)
plot_results(dcorr_CE, 'Correlation between Causal Effect and Centrality Measures', class = FALSE)
dev.off()

pdf('../Figures/Fig-CE-Class-3.pdf', width = 10, height = 10)
plot_mean_results(dclass_CE, 'Probability of Successful Classification', type = 'classification')
dev.off()

pdf('../Figures/Fig-CE-Ranking-3.pdf', width = 10, height = 10)
plot_mean_results(dranking_CE, 'Relative Ranking of Causally Most Effective Node', type = 'ranking')
dev.off()


## Correlation with CE Measure [Appendix]
# Recovers Figures 7 and 8 (appendix) and a relative ranking plot (not in paper)
##############################################################################
pdf('../Figures/Fig-CE-Corr-0.pdf', width = 10, height = 10)
plot_results(dcorr_CE0, 'Correlation between Causal Effect and Centrality Measures', class = FALSE)
dev.off()

pdf('../Figures/Fig-CE-Class-0.pdf', width = 10, height = 10)
plot_mean_results(dclass_CE0, 'Probability of Successful Classification', type = 'classification')
dev.off()

pdf('../Figures/Fig-CE-Ranking-0.pdf', width = 10, height = 10)
plot_mean_results(dranking_CE0, 'Relative Ranking of Causally Most Effective Node', type = 'ranking', max_size = 80)
dev.off()



## Correlation with ACE Measure with true mean = 0.3 [Not reported in paper]
##############################################################################
pdf('../Figures/Fig-ACE-Corr-3.pdf', width = 10, height = 10)
plot_results(dcorr_ACE, 'Correlation between Causal Effect and Centrality Measures', class = FALSE)
dev.off()

pdf('../Figures/Fig-ACE-Class-3.pdf', width = 10, height = 10)
plot_mean_results(dclass_ACE, 'Probability of Successful Classification', type = 'classification')
dev.off()

pdf('../Figures/Fig-ACE-Ranking-3.pdf', width = 10, height = 10)
plot_mean_results(dranking_ACE, 'Relative Ranking of Causally Most Effective Node', type = 'ranking')
dev.off()


## Correlation with ACE Measure with true mean = 0.0 [Not reported in paper]
##############################################################################
pdf('../Figures/Fig-ACE-Corr-0.pdf', width = 10, height = 10)
plot_results(dcorr_ACE0, 'Correlation between Causal Effect and Centrality Measures', class = FALSE)
dev.off()

pdf('../Figures/Fig-ACE-Class-0.pdf', width = 10, height = 10)
plot_mean_results(dclass_ACE0, 'Probability of Successful Classification', type = 'classification')
dev.off()

pdf('../Figures/Fig-ACE-Ranking-0.pdf', width = 10, height = 10)
plot_mean_results(dranking_ACE0, 'Relative Ranking of Causally Most Effective Node', type = 'ranking', max_size = 80)
dev.off()