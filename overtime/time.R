#!/usr/bin/Rscript
source("general.R")

data <- load_data("time.csv")

make_plot <- function(all_data, problem_name, key) {
  plt <- ggplot(data = subset(all_data, problem==problem_name),
                aes_string(x = "length", y=key, color="solver", shape="solver"))
  plt <- plt + geom_line(stat="summary", fun.y=median_100, show_guide=FALSE)
  plt <- plt + geom_point(stat="summary", fun.y=median_100, size=3)
  plt <- plt + log_x + log_y
  # plt <- plt + scale_color_manual(values=blind_colors)
  plt <- plt + scale_color_manual(values=printer_colors)
  plt <- plt + clean + ggtitle(problem_name)
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

major_data <- subset(data, solver=="p3"| solver=="ltga")

make_plot <- function(all_data, problem_name, key) {
  plt <- ggplot(data = subset(all_data, problem==problem_name),
                aes_string(x = "length", y=key, color="solver", linetype="solver"))
  plt <- plt + geom_ribbon(stat="summary", fun.data="median_hilow", alpha=0, size=1.3)
  plt <- plt + log_x + log_y
  # plt <- plt + scale_color_manual(values=blind_colors[c(3, 4, 5)])
  plt <- plt + scale_color_manual(values=printer_colors[c(3, 4, 5)])
  #plt <- plt + scale_linetype_manual(values=c(2,1,4))
  plt <- plt + clean + ggtitle(problem_name)
  return(plt)
}
result <- seven_plot(make_plot(major_data, "DeceptiveTrap", "evals"),
                     make_plot(major_data, "DeceptiveStepTrap","evals"),
                     make_plot(major_data, "HIFF", "evals"),
                     make_plot(major_data, "Rastrigin", "evals"),
                     make_plot(major_data, "NearestNeighborNK", "evals"),
                     make_plot(major_data, "IsingSpinGlass", "evals"),
                     make_plot(major_data, "MAXSAT", "evals"),
                     "Problem Size", "Median Evaluations to Success")
ggsave("evals-to-success-range.eps", result, width=10, height=12)

