#绘制10例肺癌样本生存曲线
library(survminer)
library(survival)
library(showtext)
library(ggplot2)
library(ggsignif)

font_add("myfont","C:/Windows/Fonts/simhei.ttf")
showtext_auto()

setwd("D:/aPKU/Research/20231016/")
data = read.csv("D:/aPKU/Research/20231016/10例肺癌临床数据.csv",encoding ="UTF-8")

plotdata <- data.frame(as.numeric(c(3.5,9,7,17.7,8,4.7,4,3,2.7,0.6,6)))
colnames(plotdata) <- "PFS"

plotdata$status <- rep(1,11)
plotdata$ID <- paste0('Patient-',seq(1,11))
plotdata$class <- data$最佳疗效
plotdata$IFN <- data$IFN变化率

plotdata$twoclass <- c("PR&SD","PR&SD","PR&SD","PR&SD","PR&SD","PR&SD","PD","PD","PD","PD","PD")

#plotdata <- plotdata[-12,]

fit <- survfit(Surv(PFS,status) ~ twoclass,  
               data = plotdata)


fit

custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5)
    )
}


ggsurvplot(fit, 
           data = plotdata,
           legend.labs=c("PD", "PR&SD"),
           legend.title="Best effect",
           pval = TRUE,pval.size=4,
           pval.coord=c(15,0.65),
           #ggtheme=custom_theme(),
           #title = "PFS",
           legend=c(0.8,0.8))+
  labs(x="Months elapsed",
       #y = "Percentage of survival",
       y="Progression-free Survival")

  
  


plotdata$药敏 <- data$抑制剂.白介素刺激后杀伤
plotdata$药敏 <- data$药敏

ggplot(data = plotdata,aes(x=药敏,y=PFS,colour=class,size=IFN))+
  geom_point()+
  #ggtitle("药敏与PFS") +
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="药敏",y = "PFS",title = "药敏与PFS")+
  theme(plot.title = element_text( size=14, face="bold.italic"),
        axis.title.x = element_text( size=14, face="bold"),
        axis.title.y = element_text( size=14, face="bold"))+
  scale_size(name="IFN-γ变化率")+
  scale_color_discrete(name="最佳疗效")



ggplot(data = plotdata,aes(x=IFN,y=PFS,colour=class))+
  geom_point()

ggplot(data = plotdata,aes(x=药敏,y=IFN,fill=class))+
  #geom_hline(yintercept = 1.2,col="grey",lwd=1,linetype = "dashed") +
  #geom_vline(xintercept = 0.75,col="grey",lwd=1,linetype = "dashed") + 
  geom_point(size = 13,shape=21)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="Cell viability",y = "IFN-γ Fold Change")+
  theme(plot.title = element_text( size=30, face="bold.italic"),
        axis.title.x = element_text( size=30, face="bold"),
        axis.title.y = element_text( size=30, face="bold"),
        axis.text.x=element_text(size=20,face = "bold"),
        axis.text.y=element_text(size=20,face = "bold"))+
  guides(fill =guide_legend(title.theme = element_text(size = 30,face = "bold"),
                            label.theme =element_text(size = 25,face = "bold")))+
  scale_color_discrete(name="Treatment Effect")

highlight <- plotdata[12,]

ggplot(data = plotdata,aes(x=药敏,y=IFN,fill=class))+
  #geom_hline(yintercept = 1.2,col="grey",lwd=1,linetype = "dashed") +
  #geom_vline(xintercept = 0.75,col="grey",lwd=1,linetype = "dashed") + 
  geom_point(size = 3,shape=21)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="Cell viability",y = "IFN-γ Fold Change")+
  theme(plot.title = element_text( size=30, face="bold.italic"),
        axis.title.x = element_text( size=20, face="bold"),
        axis.title.y = element_text( size=20, face="bold"))
  #scale_color_discrete(name="Treatment Effect")+
  #coord_equal()+
  geom_point(data = highlight,aes(x=药敏,y=IFN),fill=NA,size=5,colour='red',stroke=2,pch=21)


a=0.5

plotdata$score <- a*(1-plotdata$药敏)+(1-a)*plotdata$IFN

ggplot(plotdata, aes(x=reorder(ID,-score))) +
  geom_bar(aes(y=score,fill=class),stat='identity') +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))+
  scale_x_discrete(name = "Patient") +
  scale_y_continuous(name = "Score")+
  labs(fill = "Treatment Effect",title="Score and Effect")+
  scale_fill_manual(values=c("PD"='red','SD'='#d9d9d9','PR'='#2297E6'))

