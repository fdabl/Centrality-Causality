source('helpers.R')
registerDoParallel(10)


# nodes <- c(5, 10, 20, 40, 50, 60, 70, 80)
# connectivity <- seq(.1, .9, .1)
nobs <- c(2000)
nodes <- c(5, 10)
connectivity <- c(.5)
gtypes <- c('watts', 'er', 'power', 'geometric')


run_sim <- function(gtypes, fname) {
  tryCatch({
    start <- Sys.time()
    res <- run_simulation_study(nodes, gtypes, nobs = nobs, connectivity = connectivity, times = 5)
    end <- Sys.time()
    print(end - start)
    saveRDS(res, fname)
  }, error = function(e) {
    print(e)
  })
}

run_sim(gtypes, 'res-test.RDS')