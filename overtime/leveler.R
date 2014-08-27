#!/usr/bin/Rscript
source("general.R")

data <- load_data("leveler.csv")

largest <- c(805, 805, 2048, 800, 700, 784, 60)
names(largest) <- c("DeceptiveTrap", "DeceptiveStepTrap", "HIFF", "Rastrigin",
                    "NearestNeighborNK", "IsingSpinGlass", "MAXSAT")

make_plot <- function(all_data, problem_name, key) {
  plt <- ggplot(data = subset(all_data, problem==problem_name & solver=="p3"),
                aes_string(x="level", y=key, group="factor(length)", color="factor(length)"))
  plt <- plt + geom_line(stat="summary", fun.y=median_100_low, show_guide=FALSE)
  #plt <- plt + geom_point(stat="summary", fun.y=median_100_low, size=3)
  # plt <- plt + scale_color_manual(values=blind_colors)
  # plt <- plt + scale_color_manual(values=printer_colors)
  plt <- plt + clean + ggtitle(problem_name)
  return(plt)
}
# make_plot(data, "DeceptiveStepTrap", "size")

make_plot <- function(all_data, length_choices, key) {
  plt <- ggplot(data = subset(all_data, solver=="p3" & length==length_choices[problem]),
                aes_string(x="level", y=key, color="problem", shape="problem", fill="problem"))
  plt <- plt + geom_line(stat="summary", fun.y=median_100_low, show_guide=FALSE)
  plt <- plt + geom_point(stat="summary", fun.y=median_100_low, size=3)
  plt <- plt + scale_shape_manual(values=c(24, 25, 21, 22, 23, 3, 4))
  # plt <- plt + coord_cartesian(xlim = c(0, 850))
  plt <- plt + scale_color_manual(values=blind_colors) + scale_fill_manual(values=blind_colors)
  # plt <- plt + scale_color_manual(values=printer_colors) + scale_fill_manual(values=printer_colors)
  plt <- plt + clean
  plt <- plt + theme(legend.position="bottom") + guides(shape=guide_legend(nrow=2, byrow=TRUE))
  return(plt)
}

ggsave("level-size.eps", make_plot(data, largest, "size"), width=6, height=6)
ggsave("level-success.eps", make_plot(data, largest, "successes/(successes+ties+failures)") + log_y, width=6, height=6)






