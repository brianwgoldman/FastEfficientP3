#!/usr/bin/Rscript
source("general.R")

data <- load_data("time-old.csv")

make_plot <- function(all_data, problem_name, key) {
  plt <- ggplot(data = subset(all_data, problem==problem_name),
                aes_string(x = "length", y=key, color="solver", shape="solver"))
  plt <- plt + geom_line(stat="summary", fun.y=median_100, show_guide=FALSE)
  plt <- plt + geom_point(stat="summary", fun.y=median_100, size=3)
  plt <- plt + log_x + log_y
  plt <- plt + scale_color_manual(name="Optimizer", values=solver_colors, labels=solver_labels)
  plt <- plt + scale_shape_manual(name="Optimizer", values=c(16, 17, 15, 7, 3, 8), labels=solver_labels)
  plt <- plt + clean + ggtitle(problem_labels[problem_name])
  return(plt)
}

result <- seven_plot(make_plot(data, "DeceptiveTrap", "evals"),
                     make_plot(data, "DeceptiveStepTrap","evals"),
                     make_plot(data, "HIFF", "evals"),
                     make_plot(data, "Rastrigin", "evals"),
                     make_plot(data, "NearestNeighborNK", "evals"),
                     make_plot(data, "IsingSpinGlass", "evals"),
                     make_plot(data, "MAXSAT", "evals"),
                     "Problem Size", "Median Evaluations to Success")
ggsave("evals-to-success.eps", result, width=10, height=12)

data_correct_seconds <- load_data("time.csv")
result <- seven_plot(make_plot(data_correct_seconds, "DeceptiveTrap", "seconds"),
                     make_plot(data_correct_seconds, "DeceptiveStepTrap","seconds"),
                     make_plot(data_correct_seconds, "HIFF", "seconds"),
                     make_plot(data_correct_seconds, "Rastrigin", "seconds"),
                     make_plot(data_correct_seconds, "NearestNeighborNK", "seconds"),
                     make_plot(data_correct_seconds, "IsingSpinGlass", "seconds"),
                     make_plot(data_correct_seconds, "MAXSAT", "seconds"),
                     "Problem Size", "Median Seconds to Success")
ggsave("seconds-to-success.eps", result, width=10, height=12)

make_plot <- function(all_data, solver_name) {
  plt <- ggplot(data = subset(all_data, solver==solver_name),
                aes(x=length, y=pop_size, color=problem, shape=problem, fill=problem))
  plt <- plt + geom_line(stat="summary", fun.y=median_100, show_guide=FALSE)
  plt <- plt + geom_point(stat="summary", fun.y=median_100, size=3)
  plt <- plt + scale_shape_manual(name="Problem", values=problem_shapes, labels=problem_labels)
  plt <- plt + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                             labels = trans_format("log10", math_format(10^.x)), limits=c(min(data$length), max(data$length)))
  pop_sizes <- subset(data, solver=="ltga" | solver=="hboa")$pop_size
  plt <- plt + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x)), limits=c(min(pop_sizes), max(pop_sizes)))
  plt <- plt + scale_color_manual(name="Problem", values=problem_colors, labels=problem_labels)
  plt <- plt + scale_fill_manual(name="Problem", values=problem_colors, labels=problem_labels)
  plt <- plt + clean
  plt <- plt + theme(legend.position="bottom") + guides(shape=guide_legend(nrow=2, byrow=TRUE)) + xlab("Problem Size")
  return(plt)
}

ggsave("pop-ltga.eps", make_plot(data, "ltga") + ylab("Population Size"), width=6, height=6)
ggsave("pop-hboa.eps", make_plot(data, "hboa") + ylab("Population Size"), width=6, height=6)

# Comparison of algorithms
medians <- aggregate(evals~problem+solver+length, data, FUN=median_100)
# P3 best
nrow(subset(merge(aggregate(evals~problem+length, medians, FUN=min), medians), solver=="p3"))
subset(merge(aggregate(evals~problem+length, subset(medians, solver=="p3" | solver=="hc"), FUN=min), medians), solver=="hc")
subset(merge(aggregate(evals~problem+length, subset(medians, solver=="p3" | solver=="lambdalambda"), FUN=min), medians), solver=="lambdalambda")
subset(merge(aggregate(evals~problem+length, subset(medians, solver=="p3" | solver=="hboa"), FUN=min), medians), solver=="hboa")
subset(merge(aggregate(evals~problem+length, subset(medians, solver=="p3" | solver=="phboa"), FUN=min), medians), solver=="phboa")
subset(merge(aggregate(evals~problem+length, subset(medians, solver=="p3" | solver=="ltga"), FUN=min), medians), solver=="ltga")

############# Looking at largest variance ###################
bigtwo <- subset(data_correct_seconds, length==largest[problem] & (solver=="p3" | solver=="ltga"))
bigtwo$problem <- factor(problem_labels[bigtwo$problem], levels = problem_labels)
plt <- ggplot(bigtwo, aes(solver, evals)) + geom_boxplot() + facet_wrap(~problem, scales="free", ncol=3) + clean
plt <- plt + scale_x_discrete(labels=c("LTGA", "P3")) + xlab("") + ylab("Evaluations to Success")
ggsave("evals-to-success-boxplot.eps", plt + theme(strip.text = element_text(size=14)), width=10, height=12)

############## Looking at aggregate statistics
binom.test(130-10, 130, .5, alternative="greater")
