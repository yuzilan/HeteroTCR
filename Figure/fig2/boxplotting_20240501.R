getwd()

df = read.csv("df_exp1.csv")
head(df)

library(ggplot2)
library(ggsignif)

df$Model=factor(df$Model,level=c("NetTCR-1.0","NetTCR-2.0","ERGO_LSTM","ERGO_AE","DLpTCR","HeteroTCR"))

p = ggplot(data=df,aes(x=Model,y=AUC,fill=factor(Score)))+
    labs(x="",y="AUC",fill="VDJdb confidence score")+
    stat_boxplot(geom="errorbar",width=0.3,position=position_dodge(0.7),linewidth=0.3)+ 
    geom_boxplot(position=position_dodge(0.7),width=0.6,alpha=0.5,outlier.size=1,linewidth=0.3)+
    geom_point(position = position_dodge(0.7), size = 1, shape = 21, stroke=0.3) +
    theme_bw()+
    theme(panel.grid=element_blank(),
          line=element_line(linewidth=0.3),
          axis.text=element_text(family="Helvetica", size=6),
          axis.title=element_text(family="Helvetica", size=7),
          legend.text=element_text(family="Helvetica", size=6),
          legend.title=element_text(family="Helvetica", size=7)) +
    scale_fill_manual(values=c('cornflowerblue','darkorange','forestgreen'))+  # 箱线图的颜色
    scale_y_continuous(limits=c(0.52,0.82),breaks=seq(0.52,0.82,0.04))+
    geom_signif(comparisons = list(c("HeteroTCR","DLpTCR"),
                                   c("HeteroTCR","ERGO_AE"),
                                   c("HeteroTCR","ERGO_LSTM"),
                                   c("HeteroTCR","NetTCR-2.0"),
                                   c("HeteroTCR","NetTCR-1.0")),
                test="t.test",
                test.args=list(paired=TRUE,alternative="greater",var.equal=FALSE),
                map_signif_level=TRUE,
                step_increase=0.1, size=0.3, textsize = 2.5)
                
ggsave("Figure2_20240501.pdf",p,device = "pdf",width=18,height=7.4,units="cm")

# install.packages("svglite")
# library(svglite)
# ggsave("Figure2_20240501.svg",p,device = "svg",width=18,height=7.4,units="cm")