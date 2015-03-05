#!/usr/bin/Rscript
source("general.R")

data <- load_data("leveler.csv")

########## P3 Level based information #############
make_plot <- function(all_data, length_choices, key) {
  plt <- ggplot(data = subset(all_data, solver=="p3" & length==length_choices[problem]),
                aes_string(x="level", y=key, color="problem", shape="problem", fill="problem"))
  plt <- plt + geom_line(stat="summary", fun.y=median_100_low, show_guide=FALSE)
  plt <- plt + geom_point(stat="summary", fun.y=median_100_low, size=3)
  plt <- plt + scale_shape_manual(name="Problem", values=problem_shapes, labels=problem_labels)
  plt <- plt + scale_color_manual(name="Problem", values=problem_colors, labels=problem_labels)
  plt <- plt + scale_fill_manual(name="Problem", values=problem_colors, labels=problem_labels)
  plt <- plt + clean
  plt <- plt + theme(legend.position="bottom") + guides(shape=guide_legend(nrow=2, byrow=TRUE)) + xlab("Pyramid Level")
  return(plt)
}

ggsave("level-size.eps", make_plot(data, largest, "size") + ylab("Number of Solutions"), width=6, height=6)
ggsave("level-success.eps", make_plot(data, largest, "successes/(successes+ties+failures)") + ylab("Percentage of Crosovers Resulting in Fitness Improvement") + log_y, width=6, height=6)

############## P3 aggregate sizes ####################
p3pops <- aggregate(size~problem+length+solver+seed, subset(data, solver=="p3"), FUN=sum)
p3largest <- subset(p3pops, length==largest[problem])
colnames(p3largest)[5] <- "pop_size"
ltgadata <- aggregate(pop_size~problem+length+solver+seed, subset(data, length==largest[problem] & solver=="ltga"), FUN=mean)

plt <- ggplot(p3largest, aes(problem, pop_size)) + geom_boxplot() + log_y + clean + geom_point(data=ltgadata, shape='+', color="#e31a1c", size=10)
plt <- plt + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0))
plt <- plt + xlab("") + ylab("Concurrent Solutions") + scale_x_discrete(labels=problem_labels)
ggsave("solutions-stored.eps", plt, width=6, height=6)

plt <- ggplot(p3pops, aes(length, size, color=problem, shape=problem, fill=problem))
plt <- plt + geom_line(stat="summary", fun.y=median_100, show_guide=FALSE) + geom_point(stat="summary", fun.y=median_100, size=3)
plt <- plt + scale_shape_manual(name="Problem", values=problem_shapes, labels=problem_labels)
plt <- plt + scale_color_manual(name="Problem", values=problem_colors, labels=problem_labels)
plt <- plt + scale_fill_manual(name="Problem", values=problem_colors, labels=problem_labels)
plt <- plt + clean + log_x + log_y + ylab("Concurrent Solutions") + xlab("Problem Size")
plt <- plt + theme(legend.position="bottom") + guides(shape=guide_legend(nrow=2, byrow=TRUE)) 

ggsave("pop-p3.eps", plt, width=6, height=6)

