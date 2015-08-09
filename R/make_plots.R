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
    geom_bar(stat = "identity",width=.7) +
    xlab("position in read(bp)") +
    ylab("% ambigous bases")

  return(plot)
}

plot.seq.duplication.levels <- function(data_frame){
  plot <- ggplot(melt(data_frame,id.vars=c("duplication_level")),aes(duplication_level,value,fill=as.factor(variable))) +
    geom_bar(position="dodge",stat="identity") +
    scale_fill_manual(name="",labels=c("% Deduplicated","% Total"),values = c("red","blue")) +
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


#examples
json_data <- read.fastqc("~/RSV_analysis/test/fastqc/ERR303259_1_fastqc.json")
basic.info <- basic_stats(json_data)

plot.per.base.sequence.quality(per_base_sequence_quality(json_data))

plot.per.sequence.quality.scores(per_sequence_quality_scores(json_data))

plot.per.base.sequence.quality(per_base_sequence_quality(json_data))

plot.per.base.sequence.content(per_base_sequence_content(json_data))

plot.per.base.N.count(per_base_N_count(json_data))

plot.seq.duplication.levels(sequence_duplication_levels(json_data))
