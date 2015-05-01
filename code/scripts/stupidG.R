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
  my_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
} else {
  warning("Job is running locally, _not_ on Odyssey")
  my_id <- 1
}

library(SCM)
no_drop <- seq(1800, 5000, by = 100)
my_no_drop <- no_drop[my_id]

load("/n/airoldifs2/lab/afranks/covariance/singleG-thr-prot/run-1-1800/sample-780.Rdata")
state <- start

no_to_drop <- ifelse(state$lj == "abund", my_no_drop, 0)
start <- threshold_state(state, no_to_drop)

no_steps <- 20
for (i in 1:1000) {
  start <- BayesCFA(start, noSteps = no_steps, singleG = TRUE,
                    stupidG = TRUE)[[no_steps]]
  out_file <- sprintf("stupidG-thr-prot/run-%i-%i/sample-%i.Rdata",
                      my_id, my_no_drop, i)
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  save(start, file = out_file)
}

## -----------------------------------------------------------------------

get_err <- function(samp) {
  err <- with(samp,
              (diag(Tau)[tj] + diag(Xi)[kj] + diag(Theta)) / G[kj]^2
              )
  names(err) <- samp$Rnames
  err
}

load("stupidG-thr-prot/run-1-1800/sample-167.Rdata")
stupid_g_err <- get_err(start)

save(stupid_g_err, file = "~/stupid-g-error.Rdata")

## -----------------------------------------------------------------------

library(ggplot2)

load("data/stupid-g-error.Rdata")
load("data/g2-nonpooled-tech-sample.Rdata")

## Average within experiments
err2 <- tapply(err, samp$kj, mean)
stupid_g_err2 <- tapply(stupid_g_err, samp$kj, mean)

lj <- ifelse(samp$lj == "abund", "protein", "mRNA")
  
df <- data.frame(
  stringsAsFactors = FALSE,
  check.names = FALSE,
  "1 / signal to noise ratio" = c(err2, stupid_g_err2),
  which = rep(c("different scale", "common scale"), each = length(err2)),
  type = tapply(lj, samp$kj, unique) %>% rep(2),
  experiment = names(err2) %>% rep(2)
)

pdf("common-scale.pdf", width = 7.34, height = 4.17)
ggplot(df, aes(x = which, y = `1 / signal to noise ratio`, group = experiment,
               colour = experiment)) +
  geom_line() +
  geom_point() +
  facet_grid(~ type) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  guides(colour = guide_legend(ncol=2))
dev.off()

