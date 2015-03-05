#!/usr/bin/Rscript
source("general.R")

data <- load_data("meta.csv")

make_plot <- function(all_data, key) {
  plt <- ggplot(data = subset(all_data, solver=="p3"),
                aes_string(x="length", y=key, color="problem", shape="problem", fill="problem"))
  plt <- plt + geom_line(stat="summary", fun.y=median_100, show_guide=FALSE)
  plt <- plt + geom_point(stat="summary", fun.y=median_100, size=3)
  plt <- plt + scale_shape_manual(name="Problem", values=problem_shapes, labels=problem_labels)
  plt <- plt + coord_cartesian(xlim = c(0, 850))
  plt <- plt + scale_color_manual(name="Problem", values=problem_colors, labels=problem_labels)
  plt <- plt + scale_fill_manual(name="Problem", values=problem_colors, labels=problem_labels)
  plt <- plt + clean
  plt <- plt + theme(legend.position="bottom") + guides(shape=guide_legend(nrow=2, byrow=TRUE)) + xlab("Problem Size")
  return(plt)
}

ggsave("rebuilds.eps", make_plot(data, "rebuilds*length/(cross+local)") + ylab("Model Rebuild Cost per Evaluation"), width=6, height=6)
ggsave("donations.eps", make_plot(data, "donation_attempts/(cross+local)") + ylab("Donations per Evaluation"), width=6, height=6)
ggsave("cross.eps", make_plot(data, "cross/(cross+local)") + ylab("Percentage of Evaluations Spent on Crossover"), width=6, height=6)
ggsave("cross-success.eps", make_plot(data, "successes/(ties+successes+failures)") + ylab("Percentage of Crosovers Resulting in Fitness Improvement") + log_y, width=6, height=6)

