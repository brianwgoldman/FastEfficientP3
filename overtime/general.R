#!/usr/bin/Rscript
library(ggplot2)
library(scales)
library(gridExtra)

load_data <- function(filename) {
  data <- read.csv(filename, header=TRUE)
  # Ignores preliminary experiments on Deceptive 3
  data <- subset(data, problem != "Deceptive3")
  # Ensure solvers and problems appear in the order I like
  data$solver <- factor(data$solver, levels = c("hc", "lambdalambda", "hboa", "phboa", "ltga", "p3"))
  data$problem <- factor(data$problem, levels = c("DeceptiveTrap", "DeceptiveStepTrap", "HIFF", "Rastrigin",
                                                  "NearestNeighborNK", "IsingSpinGlass", "MAXSAT"))
  return(data)
}

# Treat missing runs as higher than the highest run.
# Returns NA if the median run is then missing.
median_100 <- function(data) {
  worst <- max(data, na.rm=TRUE)
  data <- c(data, rep(worst+1, 100 - length(data)))
  result <- median(data)
  if(result > worst) {
    result <- NA
  }
  return(result)
}

# Treat missing runs as lower than the lowest run.
# Returns NA if the median run is then missing.
median_100_low <- function(data) {
  worst <- min(data, na.rm=TRUE)
  data <- c(data, rep(worst-1, 100 - length(data)))
  result <- median(data)
  if(result < worst) {
    result <- NA
  }
  return(result)
}


problem_colors <- c("#E69F00","#CC79A7", "#D55E00", "#009E73", "#56B4E9", "#000000", "#0072B2", "#F0E442")
solver_colors <- c("#1f78b4", "#a6cee3", "#b2df8a", "#33a02c", "#e31a1c", "#000000", "#fb9a99")

# Convert raw solver names into prettier labels
solver_labels <- c('Hill Climber', '1+(Lambda,Lambda)', 'hBOA', 'Parameter-less hBOA', 'LTGA', 'P3')

# Lookup table which coverts raw problem names into prettier labels
problem_labels <- c("Deceptive Trap", "Deceptive Step Trap", "HIFF", "Rastrigin", "Nearest Neighbor NK", "Ising Spin Glass", "MAX-SAT")
names(problem_labels) <- c("DeceptiveTrap", "DeceptiveStepTrap", "HIFF", "Rastrigin", "NearestNeighborNK", "IsingSpinGlass", "MAXSAT")

# Largest problem size in which LTGA was successfully tuned
largest <- c(805, 805, 2048, 800, 700, 784, 60)
names(largest) <- c("DeceptiveTrap", "DeceptiveStepTrap", "HIFF", "Rastrigin",
                    "NearestNeighborNK", "IsingSpinGlass", "MAXSAT")

# Largest problem size in which hBOA was successfully tuned
hboa_largest <- c(203, 203, 512, 350, 200, 484, 40)
names(hboa_largest) <- c("DeceptiveTrap", "DeceptiveStepTrap", "HIFF", "Rastrigin",
                    "NearestNeighborNK", "IsingSpinGlass", "MAXSAT")

# Symbols to use when plotting problems
problem_shapes <- c(24, 25, 21, 22, 23, 3, 4)

# Make pretty logarithmic scales
log_x <- scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x)))
log_y <- scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x)))

# Remove excess chart junk
clean <- theme_bw() + theme(
    plot.background = element_blank()
   ,panel.grid.major = element_blank()
   ,panel.grid.minor = element_blank()
   ,panel.border = element_blank()
   ,strip.background=element_blank()
  ) + theme(axis.line = element_line(color = 'black'))

# Garbage method for plotting 7 things on one plot. Should be replaced with facet_wrap
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
