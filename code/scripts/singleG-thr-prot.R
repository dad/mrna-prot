#! /usr/bin/env Rscript
#SBATCH -n 1
#SBATCH -o out-%a.txt
#SBATCH -e err-%a.txt
#SBATCH -p airoldi
#SBATCH --mem-per-cpu=4096
#SBATCH --mail-user=afranks@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -a 1-34

library(methods)
library(magrittr)
library(falsy)

if (Sys.getenv("SLURM_JOB_ID") != "") {
  library(SCM)
  no_drop <- c(seq(1800, 5000, by = 100), 0)

  my_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  my_no_drop <- no_drop[my_id]

  state <- get_start_state()
  no_to_drop <- ifelse(state$lj == "abund", my_no_drop, 0)
  start <- threshold_state(state, no_to_drop)

  no_steps <- 20
  for (i in 1:1000) {
    start <- BayesCFA(start, noSteps = no_steps, singleG = TRUE)[[no_steps]]
    out_file <- sprintf("singleG-thr-prot/run-%i-%i/sample-%i.Rdata",
                        my_id, my_no_drop, i)
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    save(start, file = out_file)
  }
}

total_no_samples <- dir("singleG-thr-prot/run-1-1800") %>%
  gsub(pattern = "[^0-9]", replacement ="") %>%
  as.numeric() %>%
  max()

no_samples <- 100
sample_nos <- seq(
  total_no_samples - no_samples + 1,
  total_no_samples
)

sample_files <- sprintf("singleG-thr-prot/run-1-1800/sample-%i.Rdata",
                        sample_nos)

samples <- lapply(sample_files, function(file) {
  cat(".")
  e <- new.env()
  load(file, envir = e)
  get(ls(e)[1], envir = e, inherits = FALSE)
})

corrs <- list(sapply(samples, "[[", "psi")[2,])
save(corrs, file = "~/g2-nonpooled-tech-sample-corrs.RData")

L1 <- sapply(samples, function(x) x$L[,1])
L2 <- sapply(samples, function(x) x$L[,2])

save(L1, L2, file = "~/g2-nonpooled-tech-sample-L.Rdata")

df <- data.frame(
  stringsAsFactors = FALSE,
  orf = rownames(samples[[1]]$X),
  sd.mrna = apply(L2, 1, sd),
  sd.prot = apply(L1, 1, sd),
  mean.mrna = rowMeans(L2),
  mean.prot = rowMeans(L1)
)

write.table(df, file = "~/g2-nonpooled-tech-sample-summary.txt",
            quote = FALSE, row.names = FALSE)

get_err <- function(samp) {
  err <- with(samp,
              (diag(Tau)[tj] + diag(Xi)[kj] + diag(Theta)) / G[kj]^2
              )
  names(err) <- samp$Rnames
  err
}

err <- sapply(samples, get_err) %>% rowMeans()
samp <- samples[[no_samples]]
save(err, samp, file = "~/g2-nonpooled-tech-sample.Rdata")

## ------------------------------------------------------------------------

trace_data <- function(dir, what = "G", pattern = "sample-.*.Rdata") {

  sort_num <- function(x) {
    idx <- gsub("[^0-9]*", "", x) %>%
      as.numeric() %>%
      order()
    x[idx]
  }

  load_from_file <- function(file, what, dot = TRUE) {
    dot %&&% cat(".")
    tmp_env <- new.env()
    load(file, envir = tmp_env)
    tmp_env %>%
      ls() %>%
      extract(1) %>%
      get(envir = tmp_env, inherits = FALSE) %>%
      extract2(what)
  }

  dir(dir, full.names = TRUE) %>%
    grep(pattern = pattern, value = TRUE) %>%
    sort_num() %>%
    sapply(load_from_file, what = what)
}

trace_plot <- function(dir, what = "G", pattern = "sample-.*.Rdata",
                       data = NULL) {

  data <- data %||% trace_data(dir, what = what, pattern = pattern)
  data %>%
    t() %>%
    matplot()

}

load_from_file <- function(file, what) {
  cat(".")
  env <- new.env()
  load(file, envir = env)
  get(ls(env)[1], envir = env)[[what]]
}

post_plot <- function(dir, what = "G", pattern = "sample-.*.Rdata",
                      samples = -(1:100), data = NULL, layout = NULL) {

  data <- data %||% trace_data(dir, what = what, pattern = pattern)

  if (what == "G") {
    data <- data[c(1, 20), ]

  } else if (what == "psi") {
    data <- data[2, , drop = FALSE]

  }

  if (samples[1] < 0) samples <- rev(ncol(data) + samples + 1)

  if (ncol(data) < length(samples)) warning("Not enough samples")

  data <- data[, samples]

  layout <- layout %||% matrix(seq_len(nrow(data)), ncol = min(4, nrow(data)))

  layout(layout)

  for (i in seq_len(nrow(data))) {
    hist(data[i,])
  }

}

post_mean <- function(dir, what = "G", samples = -(1:100)) {

  sort_num <- function(x) {
    idx <- gsub("[^0-9]*", "", x) %>%
      as.numeric() %>%
      order()
    x[idx]
  }

  sample_files <- dir(dir, full.names = TRUE) %>%
    grep(pattern = "sample-", value = TRUE) %>%
    sort_num()

  if (samples[1] < 0) samples <- rev(length(sample_files) + samples + 1)

  sample_files <- sample_files[samples]

  res <- sapply(sample_files, load_from_file, what = what) %>%
    rowMeans()

  cat("\n")

  res
}

post_mean_runs <- function(dirs, what = "G", samples = -(1:100)) {
  sapply(dirs, post_mean, what = what, samples = samples)
}

no_drop <- seq(1800, 5000, by = 100)
dirs <- sprintf("singleG-thr-prot/run-%i-%i", seq_along(no_drop),
                no_drop)
