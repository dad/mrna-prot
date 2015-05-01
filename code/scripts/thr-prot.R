#! /usr/bin/env Rscript
#SBATCH -n 1
#SBATCH -o out-%a.txt
#SBATCH -e err-%a.txt
#SBATCH -p airoldi
#SBATCH --mem-per-cpu=4096
#SBATCH --mail-user=csardi.gabor@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -a 1-33

if (Sys.getenv("SLURM_JOB_ID") != "") {
  library(SCM)
  no_drop <- seq(1800, 5000, by = 100)

  my_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  my_no_drop <- no_drop[my_id]

  state <- get_start_state()
  no_to_drop <- ifelse(state$lj == "abund", my_no_drop, 0)
  start <- threshold_state(state, no_to_drop)

  no_steps <- 20
  for (i in 1:1000) {
    start <- BayesCFA(start, noSteps = no_steps)[[no_steps]]
    out_file <- sprintf("thr-prot/run-%i-%i/sample-%i.Rdata",
                        my_id, my_no_drop, i)
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    save(start, file = out_file)
  }
}

load_from_file <- function(file, what) {
  cat(".")
  env <- new.env()
  load(file, envir = env)
  get(ls(env)[1], envir = env)[[what]]
}

post_mean <- function(run, what = "G", samples = 801:1000,
                      in_dir = "thr-prot") {
  library(magrittr)
  sample_files <- paste0("-", run, "-") %>%
    grep(dir(in_dir), value = TRUE) %>%
    grep(pattern = run, value = TRUE) %>%
    file.path(in_dir, ., "sample-%i.Rdata") %>%
    sprintf(samples)

  res <- sapply(sample_files, load_from_file, what = what) %>%
    rowMeans()

  cat("\n")

  res
}

post_mean_runs <- function(runs, what = "G", samples = 801:1000,
                           in_dir = "thr-prot") {
  sapply(runs, post_mean, what = what, samples =samples, in_dir = in_dir)
}

nice_g_plot <- function() {

  library(ggplot2)
  library(reshape2)

  if (file.exists("thr-prot/G.Rdata")) {
    load("thr-prot/G.Rdata")
  } else {
    G <- post_mean_runs(1:33)
    save(G, file = "thr-prot/G.Rdata")
  }

  abund <- c("futcher", "gygi", "newman", "degodoy", "lee", "lu",
             "nagaraj", "peng", "thakur", "washburn", "ghaem")

  ms <- setdiff(abund, c("futcher", "ghaem", "gygi", "newman"))

  to_plot <- G[ms, ] %>%
    as.data.frame() %>%
    set_colnames(5854 - seq(1800, 5000, by = 100)) %>%
    cbind(experiment = ms) %>%
    melt(is.vars = "experiment") %>%
    set_colnames(c("experiment", "max genes allowed", "G"))

  to_plot$`max genes allowed` <- to_plot$`max genes allowed` %>%
    as.character() %>%
    as.numeric()

  load("thr-prot/run-1-1800/sample-1.Rdata")
  X <- start$X
  X[start$I == 0] <- NA
  to_plot_pts <- colSums(is.na(X)) %>%
    tapply(start$kj, min) %>%
    extract(ms) %>%
    subtract(5854, .) %>%
    data.frame() %>%
    set_names("genes_observed") %>%
    cbind(experiment = ms)

  ## Round to closest point measured
  to_plot_pts$`max genes allowed` <-
    sapply(to_plot_pts$genes_observed, function(x) {
      to_plot[,2][which.min(abs(x - to_plot[,2]))]
    })

  to_plot <- merge(to_plot_pts, to_plot, all = TRUE)

  (
    ggplot(to_plot, aes(x = `max genes allowed`, y = G, group = experiment,
                        colour = experiment)) +
      scale_x_reverse() +
      geom_line() +
      geom_point(aes(x = genes_observed, y = G), size = 3) +
      scale_size(2)
  ) %>%
    print()

  ggsave("G-vs-miss.pdf", width = 14)
}
