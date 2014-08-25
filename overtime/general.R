#!/usr/bin/Rscript
library(ggplot2)
library(scales)
library(gridExtra)

load_data <- function(filename) {
  data <- read.csv(filename, header=TRUE)
  data <- subset(data, problem != "Deceptive3")
  data$solver <- factor(data$solver, levels = c("lambdalambda", "hc", "ltga", "p3", "hboa", "phboa"))
  return(data)
}

median_100 <- function(data) {
  worst <- max(data, na.rm=TRUE)
  data <- c(data, rep(worst+1, 100 - length(data)))
  result <- median(data)
  if(result > worst) {
    result <- NA
  }
  return(result)
}


blind_colors <- c("#CC79A7", "#E69F00","#56B4E9", "#000000", "#D55E00", "#009E73", "#F0E442", "#0072B2")
printer_colors <- c("#a6cee3", "#1f78b4", "#fb9a99", "#000000", "#b2df8a", "#33a02c", "#e31a1c")


log_x <- scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x)))
log_y <- scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x)))

clean <- theme_bw() + theme(
    plot.background = element_blank()
   ,panel.grid.major = element_blank()
   ,panel.grid.minor = element_blank()
   ,panel.border = element_blank()
  ) + theme(axis.line = element_line(color = 'black'))

seven_plot <- function(p1, p2, p3, p4, p5, p6, p7, xlabel, ylabel) {
  tmp <- ggplot_gtable(ggplot_build(p1))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  group_legend <- tmp$grobs[[leg]]

  nl <- theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank()) 
  result <- arrangeGrob(arrangeGrob(p1 + nl),
                       arrangeGrob(p2 + nl),
                       arrangeGrob(p3 + nl),
                       arrangeGrob(p4 + nl),
                       arrangeGrob(p5 + nl),
                       arrangeGrob(p6 + nl),
                       arrangeGrob(p7 + nl),
                       group_legend,
                       sub=paste(xlabel, "\n", sep=""), left=paste("\n", ylabel, sep=""), ncol=3)
  return(result)
}
