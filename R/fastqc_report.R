library(rjson)
library(reshape2)
library(dplyr)
library(ggplot2)

#' read a json formated fastqc data file and return a list object
#' @param filename A path to the json formatted fastqc file
#' @return A list
#' @examples
#' read.fastqc("path_to_file",format="json")
read.fastqc <- function(filename,format="json"){
  return(fromJSON(file=filename))
}

#' get basic statistics
#' This function returns the basic statics in the fastqc results module.
#' @param jsonlist A list object. The list can be derived from read.fastqc function
#' @return A dataframe
#' @examples
#' basic_stats(listobject)
basic_stats <- function(jsonlist){
  dta <- data.frame(matrix(unlist(jsonlist$`Basic Statistics`$contents), ncol = 7, byrow =T))
  dta <- dta %>%
    transmute(
      file_name  = X4,
      encoding = X2,
      file_type = X3,
      sequence_length = X5,
      poor_quality_seqs = X6,
      total_sequences = X7,
      gc_percent = X1
    )
  dta[] <- lapply(dta, as.character)
  return(dta)
}

#get the per base sequence quality
#' This function returns the per base sequence quality scores
#' @param
per_base_sequence_quality <- function(json_file){
  dta <- data.frame(matrix(unlist(json_file$`per_base_sequence_quality`$contents), ncol = 7, byrow =T))
  dta <- dta %>%
    transmute(
         base = X3,
         p10  = X1,
         p90  = X2,
         lower_quartile = X4,
         mean = X5,
         median = X6,
         upper_quartile = X7
         )
  dta[] <- lapply(dta, as.character)
  dta[] <- lapply(dta, as.numeric)
  return(dta)
}

#get per base sequence content
per_base_sequence_content <- function(json_file){
  dta <- data.frame(matrix(unlist(json_file$`Per base sequence content`$contents),ncol = 5,byrow =T))
  dta <- dta%>%
    transmute(base=X2,
              A=X1,
              C=X3,
              G=X4,
              T=X5)

  dta[] <- lapply(dta, as.character)
  dta[] <- lapply(dta, as.numeric)
  return(dta)
}

per_base_N_count <- function(json_file){
  dta <- data.frame(matrix(unlist(json_file$`Per base N content`$contents), ncol = 2, byrow =T))
  dta <- dta %>%
    transmute(
      base = X1,
      N_count  = X2
    )
  dta[] <- lapply(dta, as.character)
  dta[] <- lapply(dta, as.numeric)
  return(dta)
}


#sequence duplication levels
sequence_duplication_levels <- function(json_file){
  dta <- data.frame(matrix(unlist(json_file$`Sequence Duplication Levels`$contents), ncol = 3, byrow =T))
  dta <- dta %>%
    transmute(
      duplication_level = X1,
      percent_deduplicated  = X2,
      percent_total = X3
    )
  dta[] <- lapply(dta, as.character)
  dta   <- transmute(dta,
                     duplication_level = duplication_level,
                     percent_deduplicated = as.numeric(percent_deduplicated),
                     percent_total = as.numeric(percent_total))

  dta <- dta$duplication_level <- factor(as.character( dta$duplication_level),
                                                     levels=c("1","2","3","4","5","6","7","8","9",">10",">50",">100",">500",">1k",">5k",">10k+"))
  return(dta)
}

#sequence quality scores
per_sequence_quality_scores <- function(json_file){
  dta <- data.frame(matrix(unlist(json_file$`Per sequence quality scores`$contents), ncol = 2, byrow =T))
  dta <- dta %>%
    transmute(
      quality = X2,
      count  = X1
    )
  dta[] <- lapply(dta, as.character)
  dta[] <- lapply(dta, as.numeric)
  return(dta)
}

# length distribution
sequence_length_distribution <- function(json_file){
  dta <- data.frame(matrix(unlist(json_file$`Sequence Length Distribution`$contents), ncol = 2, byrow =T))
  dta <- dta %>%
    transmute(
      count  = X1,
      length = X2
    )
  dta[] <- lapply(dta, as.character)
  dta[] <- lapply(dta, as.numeric)
  return(dta)
}


per.base.sequence.content <- per_base_sequence_content(json_data)

per_base_sequence_content.m <- melt(per.base.sequence.content,id.vars = "base")
per.base.sequence.content.plot <- ggplot(per_base_sequence_content.m, aes(x = base, y = value, fill =variable )) +
  geom_bar(stat = "identity",width = .7) +
  scale_fill_manual(name="",values=c("red","blue","purple","grey30")) +
  theme(axis.text.x  = element_text(size=10),
        legend.position="") +
  ylab("proportion") +
  xlab("")


#per base N count
per.base.N.count <- per_base_N_count(json_data)

per.base.N.count.plot <- ggplot(per.base.N.count,aes(base,N_count)) +
  geom_bar(stat = "identity",width=.7) +
  xlab("position in read(bp)") +
  ylab("% ambigous bases")

seq.duplication.levels <- sequence_duplication_levels(json_data)

seq.duplication.levels.long <- melt(seq.duplication.levels,id.vars=c("duplication_level"))

seq.duplication.levels.plot <- ggplot(seq.duplication.levels.long,aes(duplication_level,value,fill=as.factor(variable))) +
  geom_bar(position="dodge",stat="identity") +
  scale_fill_manual(name="",labels=c("% Deduplicated","% Total"),values = c("red","blue")) +
  ylim(0,100) +
  theme(axis.text.x  = element_text(size=10)) +
  xlab("sequence deduplication level") +
  ylab("")


json_data <- read.fastqc("~/RSV_analysis/test/fastqc/ERR303259_1_fastqc.json")
basic.info <- basic_stats(json_data)

per.sequence.quality.scores <- per_sequence_quality_scores(json_data)

per.base.sequence.quality <- per_base_sequence_quality(json_data)

per.sequence.quality.scores.plot <- ggplot(per.sequence.quality.scores,aes(quality,count)) +
  geom_bar(stat="identity",width=.5)

per.base.sequence.quality.plot <- ggplot(per.base.sequence.quality,aes(base,mean)) +
  geom_smooth(colour="red",method="loess",se=F) +
  geom_smooth(aes(base,median),colour="blue",method="loess",se=F) +
  ylim(1,40) +
  xlim(1,150) +
  ylab("phred scores")

