
library(ggplot2)
library(magrittr)
library(gridExtra)
library(reshape2)
library(scales)

load("data/g2-nonpooled-tech-sample.Rdata")

state <- samp
col <- "expr.comm.causton.peroxidea"

invlogit <- function(x, eta0, eta1) {
  1 / (1 + exp( -(eta0 + eta1 * x)))
}

nmar_figure3 <- function(state, superset, subsets, nms) {
  cols <- c(superset,subsets)
  Isuper <- state$I[, c(superset, subsets)] *
    state$I[, rep(superset, length(cols))]

  X <- exp(state$X)
  X <- X[,rep(superset,length(cols))]
  colnames(X) <- colnames(Isuper) <- nms

  df <- data.frame(stringsAsFactors=FALSE,
                   melt(X,as.factor=FALSE))
  colnames(df) <- c("Gene","Experiment","Protein level (mol./cell)")
  df$I <- as.vector(Isuper)

  df.obs <- df[df$I==1,]
  ggplot(df.obs,
    aes(x=`Protein level (mol./cell)`, y = ..count.., fill = `Experiment`)) +
    xlim(c(3,13.5))+
    geom_density(alpha = .3, size = .1) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
#    geom_vline(
#      data=aggregate(df.obs[3], df.obs[2], mean),
#      mapping=aes(xintercept=`Protein level (mol./cell)`,
#        color=Experiment), linetype="dashed") +
    theme_bw() +
    theme(
      plot.margin = unit(c(.2, .5, .2, .5), "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      aspect.ratio=1,
      legend.position="bottom",
      legend.title = element_blank()) +
    guides(fill = guide_legend(ncol = 3))
}

hist_miss_1 <- function(state, col) {

  kj <- state$kj %>%
    set_names(state$Rnames)
  eta0 <- state$eta0[kj[col]]
  eta1 <- state$eta1[kj[col]]

  df <- data.frame(
    stringsAsFactors = FALSE,
    m = state$X[,col] %>% exp(),
    o = ifelse(state$I[,col], "observed", "missing") %>%
           factor(levels = c("observed", "missing", "total"))
  ) %>% set_names(c("mRNA level", "Causton"))

  df2 <- df[, 1, drop = FALSE] %>%
    cbind(Causton = factor("total"))
  
  p2 <- ggplot(df, aes(x = `mRNA level`, y = ..count.., fill = Causton)) +
    geom_density(alpha = .3, size = .1) +
    geom_density(data = df2, fill = "#FFFFFF", alpha = 0) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    theme_bw() +
    theme(plot.margin = unit(c(.2, .5, .2, .5), "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          aspect.ratio = 1,
          legend.position = "bottom")

  dfcens <- data.frame(
    check.names = FALSE,
    `mRNA level` = exp(seq(min(log(df[,1])), max(log(df[,1])), length.out = 100))
  )

  dfcens$`P(obs)` <- sapply(log(dfcens[,1]), invlogit, eta0 = eta0, eta1 = eta1)
  dfcens$Causton <- factor("")
  
  p1 <- ggplot(dfcens, aes(x = `mRNA level`, y = `P(obs)`)) +
     geom_line() +
     scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                   labels = trans_format("log10", math_format(10^.x))) +
     theme_bw() +
     theme(axis.title.x = element_blank(),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           plot.margin = unit(c(.2,.5,.2,.5), "cm"),
           aspect.ratio = .32)


  cols <- c("abund.wblot.ghaem.1",
            "abund.ms.peng.1",
            "abund.gfp.newman.1")
  p3 <- nmar_figure3(state, cols[1], cols[2:3],
                     c("All Ghaemmaghami", "Detected in Peng",
                       "Detected in Newman"))
  
  pdf("miss-hist.pdf", width = 9*1.2, height = 6*1.2)
  grid.arrange(arrangeGrob(virtualGrob, p1, p3, p2, heights=c(1/4, 3/4), ncol=2))
  dev.off()
}

hist_miss_1(state, col)







df_row_names <- function(x, name="name") {
  x <- data.frame(rownames(x), x, stringsAsFactors = FALSE)
  rownames(x) <- NULL
  colnames(x)[1] <- name
  x
}

nmar_figure <- function(state) {

  X <- state$X
  X[state$I == 0] <- NA
  
  X <- scale(X) %>%
    data.frame(stringsAsFactors = FALSE, orf = rownames(X), .) %>%
    melt(id.vars = "orf") %>%
    set_names(c("orf", "experiment", "value"))
  
  X$type <- grepl("^abund", X$experiment) %>% ifelse("protein", "mRNA")

  miss <- tapply(is.na(X$value), X$orf, FUN=sum) %>%
    cut(breaks=c(-1,5,10,15,20,30,40,60)) %>%
    data.frame(orf = X$orf, missing = ., stringsAsFactors = FALSE)

  X$missing <- miss$missing[match(miss$orf, X$orf)]

  names(X) <- c("orf", "experiment", "scaled raw value", "type", "missing")
  
  ggplot(X, aes(missing, `scaled raw value`)) +
    geom_boxplot() +
    facet_wrap(~ type, ncol = 1) +
    theme_bw()

  ggsave("nmar-miss.pdf", height = 5.34, width = 5.34)
}

nmar_figure2 <- function(state,cols,nms=cols) {
    
    p1 <- with(state,{
        X <- X[,cols]
        colnames(X) <- nms
        xCenter <- X-rep(nu[cols],each=N)
        df <- data.frame(stringsAsFactors=FALSE,
                         melt(xCenter,as.factor=FALSE))
        colnames(df) <- c("Gene","Experiment","Protein Abundance")
        df$I <- as.vector(I[,cols])

        df.obs <- df[df$I==1,]
        p1 <- ggplot(df.obs, aes(x=`Protein Abundance`,
                            y = ..count..,
                            fill = `Experiment`), xlim=c(-3,15)) +
              geom_density(alpha = .3, size = .1) +
              geom_vline(data=aggregate(df.obs[3],df.obs[2], mean),
                         mapping=aes(xintercept=`Protein Abundance`, color=Experiment),linetype="dashed") +
              theme_bw() +
                  theme(plot.margin = unit(c(.2, .5, .2, .5), "cm"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank())
        p1
    })
    ggsave("nmar-old.pdf")
    print(p1)
    p1
}
