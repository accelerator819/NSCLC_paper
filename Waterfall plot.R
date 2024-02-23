library(ggplot2)
library(ggbreak)
library(dplyr)
library(cowplot)

setwd("D:/aPKU/Research/2023/20231016/")
data = read.csv("D:/aPKU/Research/2023/20231016/10例肺癌临床数据.csv",encoding ="UTF-8")

data$'-log2FC(药敏/0.7)' <- -log(data$药敏/0.7,2)
colnames(data)[13] <- 'scaled_药敏'


data$敏感与否 <- 'Sensitive'
for (i in 1:11) 
{
  if(data[i,11]>0.7)
  {
    data$敏感与否[i] <- 'Resistant'
  }
}

data$文章编号 <- substr(data$编号,nchar(data$编号)-3,nchar(data$编号))



p1 <- ggplot(data, aes(x=reorder(文章编号,-药敏))) +
  geom_bar(aes(y=药敏,fill=最佳疗效),stat='identity') +
  #scale_y_break(c(4,130), scale = 0.2) +
  #geom_hline(yintercept = 0.75) +
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))+
  scale_x_discrete(name = "Patient") +
  scale_y_continuous(name = "Cell viability")+
  labs(fill = "Treatment Effect")+
  scale_fill_manual(values=c("PD"='red','SD'='#d9d9d9','PR'='#2297E6'))

p1


p2 <- ggplot(data, aes(x=reorder(文章编号,-IFN变化率))) +
  geom_bar(aes(y=IFN变化率,fill=最佳疗效),stat='identity') +
  # scale_y_break(c(4,130), scale = 0.2) +
  #geom_hline(yintercept = 1) +
  #geom_hline(yintercept = 0) +
  #geom_hline(yintercept = -0.3) +
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))+
  scale_x_discrete(name = "Patient") +
  scale_y_continuous(name = "IFN-γ Fold Change")+
  labs(fill = "Treatment Effect")+
  #scale_fill_discrete(labels = c("Resistant","Sensitive"))+
  scale_fill_manual(values=c("PD"='red','SD'='#d9d9d9','PR'='#2297E6'))+
  theme(text = element_text(size = 15))
  
p2