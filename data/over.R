#!/usr/bin/Rscript
source("general.R")

data <- load_data("over.csv")

make_plot <- function(all_data, problem_name, length_choice) {
  plt <- ggplot(data = subset(all_data, problem==problem_name & length==length_choice),
                aes(x = evals, y=fitness, color=solver, shape=solver))
  plt <- plt + geom_line(stat="summary", fun.y=median_100_low, show_guide=FALSE)
  plt <- plt + geom_point(stat="summary", fun.y=median_100_low, size=3)
  plt <- plt + log_x
  plt <- plt + scale_color_manual(name="Optimizer", values=solver_colors, labels=solver_labels)
  plt <- plt + scale_shape_manual(name="Optimizer", values=c(16, 17, 15, 7, 3, 8), labels=solver_labels)
  plt <- plt + clean + ggtitle(paste(problem_labels[problem_name], length_choice))
  return(plt)
}

result <- seven_plot(make_plot(data, "DeceptiveTrap", 203),
                     make_plot(data, "DeceptiveStepTrap",203),
                     make_plot(data, "HIFF", 512),
                     make_plot(data, "Rastrigin", 350),
                     make_plot(data, "NearestNeighborNK", 200),
                     make_plot(data, "IsingSpinGlass", 484),
                     make_plot(data, "MAXSAT", 40),
                     "Evaluations", "Median Best Fitness")
ggsave("fitness-over-time.eps", result, width=10, height=12)

low_pop <- read.csv("overlow.csv", header=TRUE)

make_plot <- function(all_data, problem_name, length_choice) {
  plt <- ggplot(data = subset(all_data, problem==problem_name & length==length_choice),
                aes(x = evals, y=fitness, color=solver, shape=solver))
  plt <- plt + geom_line(stat="summary", fun.y=median_100_low, show_guide=FALSE)
  plt <- plt + geom_point(stat="summary", fun.y=median_100_low, size=3)
  plt <- plt + log_x
  plt <- plt + scale_color_manual(name="Optimizer", values=solver_colors[c(5, 7, 6)], labels=c("LTGA_Original", "LTGA_Tenth", "P3"))
  plt <- plt + scale_shape_manual(name="Optimizer", values=c(3, 4, 8), labels=c("LTGA_Original", "LTGA_Tenth", "P3"))
  plt <- plt + clean + theme(legend.position="bottom") + xlab("Evaluations") + ylab("Median Best Fitness")
  return(plt)
}
ggsave("small-pop-dst.eps", make_plot(low_pop, "DeceptiveStepTrap",203), width=6, height=6)
ggsave("small-pop-nk.eps", make_plot(low_pop, "NearestNeighborNK", 200), width=6, height=6)

medians <- aggregate(fitness~problem+solver+length+evals, subset(data, length==hboa_largest[problem]), FUN=median_100)
# How often P3 had the best of all
nrow(subset(merge(aggregate(fitness~problem+evals, medians, FUN=max), medians), solver=="p3"))
# HC better than P3
181 - nrow(subset(merge(aggregate(fitness~problem+evals, subset(medians, solver=="p3" | solver=="hc"), FUN=max), medians), solver=="p3"))
# Lambda Lambda better than P3
181 - nrow(subset(merge(aggregate(fitness~problem+evals, subset(medians, solver=="p3" | solver=="lambdalambda"), FUN=max), medians), solver=="p3"))
# hBOA better than P3
181 - nrow(subset(merge(aggregate(fitness~problem+evals, subset(medians, solver=="p3" | solver=="hboa"), FUN=max), medians), solver=="p3"))
# Parameter-less hBOA better than P3
181 - nrow(subset(merge(aggregate(fitness~problem+evals, subset(medians, solver=="p3" | solver=="phboa"), FUN=max), medians), solver=="p3"))
# LTGA better than P3
181 - nrow(subset(merge(aggregate(fitness~problem+evals, subset(medians, solver=="p3" | solver=="ltga"), FUN=max), medians), solver=="p3"))

# Statistical test
binom.test(181-50, 181, .5, alternative="greater")

# hBOA worst
nrow(subset(merge(aggregate(fitness~problem+evals, medians, FUN=min), medians), solver=="hboa"))
# Parameter-less hBOA better than hBOA
nrow(subset(merge(aggregate(fitness~problem+evals, subset(medians, solver=="hboa" | solver=="phboa"), FUN=max), medians), solver=="phboa"))
