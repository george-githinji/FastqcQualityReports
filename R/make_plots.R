require(cowplot)
require(reshape2)

source("/Users/george/Code/R/FastqcQualityReports/R/get_data.R")

plot.per.base.sequence.content <- function(data_frame){
  plot <- ggplot(melt(data_frame,id.vars="base"), aes(x = base, y = value, colour =variable )) +
    geom_line(aes(group=variable)) +
    #geom_bar(stat = "identity",width = .7) +
    scale_colour_manual(name="bases",values=c("red","blue","green","grey30")) +
    theme(axis.text.x  = element_text(size=10),
          legend.position="") +
    ylim(0,100) +
    ylab("proportion") +
    xlab("")

  return(plot)
}

plot.per.base.N.count <- function(data_frame){
  plot <- ggplot(data_frame,aes(base,N_count)) +
    geom_bar(stat = "identity",width=.5,fill="red") +
    ylim(0,1) +
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
    #geom_bar(stat="identity",width=.4) +
    geom_point(colour="red") +
    geom_line() +
  xlim(0,40)

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

plot.seq.length.distribution <- function(data_frame){
  plot <- ggplot(data_frame,aes(factor(length),count)) +
    geom_bar(width=.5,stat="identity",fill="grey70") +
    #geom_freqpoly(stat = "identity") +
   # scale_x_discrete(breaks=seq(0,150,50)) +
    #xlim(0,150) +
    ylab("count") +
    xlab("length") 
  
  return(plot)
}
