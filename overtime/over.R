#!/usr/bin/Rscript
source("general.R")

data <- load_data("over.csv")

make_plot <- function(all_data, problem_name, length_choice) {
  plt <- ggplot(data = subset(all_data, problem==problem_name & length==length_choice),
                aes(x = evals, y=fitness, color=solver, shape=solver))
  plt <- plt + geom_line(stat="summary", fun.y=median_100, show_guide=FALSE)
  plt <- plt + geom_point(stat="summary", fun.y=median_100, size=3)
  plt <- plt + log_x
  # plt <- plt + scale_color_manual(values=blind_colors)
  plt <- plt + scale_color_manual(values=printer_colors)
  plt <- plt + clean + ggtitle(problem_name)
  return(plt)
}

result <- seven_plot(make_plot(data, "DeceptiveTrap", 154),
                     make_plot(data, "DeceptiveStepTrap",154),
                     make_plot(data, "HIFF", 128),
                     make_plot(data, "Rastrigin", 150),
                     make_plot(data, "NearestNeighborNK", 150),
                     make_plot(data, "IsingSpinGlass", 144),
                     make_plot(data, "MAXSAT", 40),
                     "Evaluations", "Median Best Fitness")
ggsave("fitness-over-time.eps", result, width=10, height=12)

