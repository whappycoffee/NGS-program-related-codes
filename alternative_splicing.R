library(ggplot2)
iPS = read.csv("iPS/ToFU/AS_summary.csv", header = TRUE)
S1 = read.csv("S1/ToFU/AS_summary.csv", header = TRUE)
S2 = read.csv("S2/ToFU/AS_summary.csv", header = TRUE)
S3 = read.csv("S3/ToFU/AS_summary.csv", header = TRUE)
ISM.D0 = read.csv("ISM.D0/ToFU/AS_summary.csv", header = TRUE)
ISM.D4 = read.csv("ISM.D4/ToFU/AS_summary.csv", header = TRUE)
ADM.D0 = read.csv("ADM.D0/ToFU/AS_summary.csv", header = TRUE)
ADM.D4 = read.csv("ADM.D4/ToFU/AS_summary.csv", header = TRUE)
View(iPS)

iPS = c(sum(iPS$Alt5), sum(iPS$Alt3), sum(iPS$ES), sum(iPS$IR), sum(iPS$Total))
S1 = c(sum(S1$Alt5), sum(S1$Alt3), sum(S1$ES), sum(S1$IR), sum(S1$Total))
S2 = c(sum(S2$Alt5), sum(S2$Alt3), sum(S2$ES), sum(S2$IR), sum(S2$Total))
S3 = c(sum(S3$Alt5), sum(S3$Alt3), sum(S3$ES), sum(S3$IR), sum(S3$Total))
ISM.D0 = c(sum(ISM.D0$Alt5), sum(ISM.D0$Alt3), sum(ISM.D0$ES), sum(ISM.D0$IR), sum(ISM.D0$Total))
ISM.D4 = c(sum(ISM.D4$Alt5), sum(ISM.D4$Alt3), sum(ISM.D4$ES), sum(ISM.D4$IR), sum(ISM.D4$Total))
ADM.D0 = c(sum(ADM.D0$Alt5), sum(ADM.D0$Alt3), sum(ADM.D0$ES), sum(ADM.D0$IR), sum(ADM.D0$Total))
ADM.D4 = c(sum(ADM.D4$Alt5), sum(ADM.D4$Alt3), sum(ADM.D4$ES), sum(ADM.D4$IR), sum(ADM.D4$Total))

AS_events = c("Alt5", "Alt3", "ES", "IR", "Total")
AS = as.data.frame(iPS, AS_events)
AS$S1 = S1
AS$S2 = S2
AS$S3 = S3
AS$ISM.D0 = ISM.D0
AS$ISM.D4 = ISM.D4
AS$ADM.D0 = ADM.D0
AS$ADM.D4 = ADM.D4
AS$events = AS_events

numbers = c(iPS, S1, S2, S3, ISM.D0, ISM.D4, ADM.D0, ADM.D4)
stage = c(rep("iPS", time = 5), rep("S1", time = 5), rep("S2", time = 5), rep("S3", time = 5), rep("ISM.D0", time = 5), rep("ISM.D4", time = 5), rep("ADM.D0", time = 5), rep("ADM.D4", time = 5))
stage
events = c(rep(c("Alt5", "Alt3", "ES", "IR", "Total"), times = 8))
AS = as.data.frame(numbers, stage)
AS$events = events
View(AS)
AS$stage = stage
row.names(AS) = seq(1,40)

ggplot(AS, aes(x = events,y = numbers, fill = stage)) + geom_bar(stat = "identity",position = "dodge") + 
  xlab("Alternative Splicing events") + ylab("Numbers") 
  
### iPS cell line
iPSline = subset(AS, stage == "iPS"|stage == "S1"| stage == "S2" | stage=="S3")

ggplot(iPSline, aes(x = events,y = numbers, fill = stage)) + geom_bar(stat = "identity",position = "dodge") + 
  xlab("Alternative Splicing events") + ylab("Numbers") 

### ISM cell line
ISMline = subset(AS, stage == "ISM.D0"| stage == "ISM.D4")
ggplot(ISMline , aes(x = events,y = numbers, fill = stage)) + geom_bar(stat = "identity",position = "dodge") + 
  xlab("Alternative Splicing events") + ylab("Numbers")
View(ISMline)
### ADM cell line
ADMline = subset(AS, stage == "ADM.D0"| stage == "ADM.D4")
ggplot(ADMline, aes(x = events,y = numbers, fill = stage)) + geom_bar(stat = "identity",position = "dodge") + 
  xlab("Alternative Splicing events") + ylab("Numbers") 

  