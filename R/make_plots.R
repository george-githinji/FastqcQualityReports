require(cowplot)
require(reshape2)
source("R/get_data.R")

plot.per.base.sequence.content <- function(data_frame){
  plot <- ggplot(melt(data_frame,id.vars="base"), aes(x = base, y = value, fill =variable )) +
    geom_bar(stat = "identity",width = .7) +
    scale_fill_manual(name="",values=c("red","blue","purple","grey30")) +
    theme(axis.text.x  = element_text(size=10),
          legend.position="") +
    ylab("proportion") +
    xlab("")

  return(plot)
}

plot.per.base.N.count <- function(data_frame){
  plot <- ggplot(data_frame,aes(base,N_count)) +
    geom_bar(stat = "identity",width=.7,fill="red") +
    xlab("position in read(bp)") +
    ylab("% ambigous bases")

  return(plot)
}

plot.seq.duplication.levels <- function(data_frame){
  plot <- ggplot(melt(data_frame,id.vars=c("duplication_level")),aes(duplication_level,value,colour=as.factor(variable),ymax=max(value)*1.05)) +
    #geom_bar(position="dodge",stat="identity",width=.7) +
    geom_line(position=position_dodge(.1),aes(group=variable)) +
    #scale_y_discrete(breaks=seq(0, 100, 10)) +
    scale_colour_manual(name="",labels=c("Distinct","Total"),values = c("red","blue")) +
    ylim(0,100) +
    theme(axis.text.x  = element_text(size=10)) +
    xlab("sequence deduplication level") +
    ylab("")

  return(plot)
}

plot.per.sequence.quality.scores <- function(data_frame){
  plot <- ggplot(data_frame,aes(quality,count)) +
    geom_bar(stat="identity",width=.5)

  return(plot)
}

plot.per.base.sequence.quality <- function(data_frame){
  plot <- ggplot(data_frame,aes(base,mean)) +
    geom_smooth(colour="red",method="loess",se=F) +
    geom_smooth(aes(base,median),colour="blue",method="loess",se=F) +
    ylim(1,40) +
    xlim(1,150) +
    ylab("phred scores")

  return(plot)
}


json_data1 <- read.fastqc("~/RSV_analysis/test/fastqc/ERR303259_1_fastqc/ERR303259_1_fastqc.json")
json_data2 <- read.fastqc("~/RSV_analysis/test/fastqc/ERR303259_2_fastqc/ERR303259_2_fastqc.json")

json_data1 <- read.fastqc("~/RSV_analysis/test/fastqc/ERR438932_1_fastqc/ERR438932_1_fastqc.json")
json_data2 <- read.fastqc("~/RSV_analysis/test/fastqc/ERR438932_2_fastqc/ERR438932_2_fastqc.json")


  basic.info1 <- basic_stats(json_data1)
  basic.info2 <- basic_stats(json_data2)
  file1_name <- tools::file_path_sans_ext(basic.info1$file_name)
  file2_name <- tools::file_path_sans_ext(basic.info2$file_name)


  plot_grid(
    plot.per.base.sequence.quality(per_base_sequence_quality(json_data1)) + ggtitle(file1_name),
    plot.per.base.sequence.quality(per_base_sequence_quality(json_data2)) + ggtitle(file2_name)
  )

  plot_grid(
    plot.per.sequence.quality.scores(per_sequence_quality_scores(json_data1)) + ggtitle(file1_name),
    plot.per.sequence.quality.scores(per_sequence_quality_scores(json_data2)) + ggtitle(file2_name)
  )

  plot_grid(
    plot.per.base.sequence.content(per_base_sequence_content(json_data1)) + ggtitle(file1_name),
    plot.per.base.sequence.content(per_base_sequence_content(json_data2)) + ggtitle(file2_name),
    nrow = 2
  )

  plot_grid(
    plot.per.base.N.count(per_base_N_count(json_data1)) + ggtitle(file1_name),
    plot.per.base.N.count(per_base_N_count(json_data2)) + ggtitle(file2_name),
    nrow=2
  )

  plot_grid(
    plot.seq.duplication.levels(sequence_duplication_levels(json_data1)) + ggtitle(file1_name),
    plot.seq.duplication.levels(sequence_duplication_levels(json_data2)) + ggtitle(file2_name),
    nrow = 2
  )



