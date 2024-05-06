getwd()
library(ggplot2)
library(ggsignif)

# vdj0
df = read.csv("df_exp2_vdj0.csv")
head(df)
df$Model=factor(df$Model,level=c("NetTCR-1.0","NetTCR-2.0","ERGO_LSTM","ERGO_AE","DLpTCR","HeteroTCR"))
df$based=factor(df$based,level=c("Pair-based","TCR-based","Antigen-based","Strict-based"))

p = ggplot(data=df,aes(x=Model,y=AUC,fill=factor(based)))+
    labs(x="",y="AUC",fill="Four types of data splitting methods")+
    stat_boxplot(geom="errorbar",width=0.3,position=position_dodge(0.7),linewidth=0.3)+
    geom_boxplot(position=position_dodge(0.7),width=0.6,alpha=0.5,outlier.size=1,linewidth=0.3)+
    geom_point(position = position_dodge(0.7), size=1, shape=21, stroke=0.3) +
    theme_bw()+
    theme(panel.grid=element_blank(),
          line=element_line(linewidth=0.3),
          axis.text=element_text(family="Helvetica", size=6),
          axis.title=element_text(family="Helvetica", size=7),
          legend.text=element_text(family="Helvetica", size=6),
          legend.title=element_text(family="Helvetica", size=7)) +
    scale_fill_manual(values=c('cornflowerblue','darkorange','forestgreen','firebrick'))+
    scale_y_continuous(limits=c(0.48,0.84),breaks=seq(0.48,0.84,0.04))+
    geom_signif(comparisons = list(c("HeteroTCR","DLpTCR"),
                                   c("HeteroTCR","ERGO_AE"),
                                   c("HeteroTCR","ERGO_LSTM"),
                                   c("HeteroTCR","NetTCR-2.0"),
                                   c("HeteroTCR","NetTCR-1.0")),
                test="t.test",
                test.args=list(paired=TRUE,alternative="greater",var.equal=FALSE),
                map_signif_level=TRUE,
                step_increase=0.1, size=0.3, textsize = 2.5)

ggsave("Figure3_20240501.pdf",p,device = "pdf",width=18,height=7.4,units="cm")


#vdj1
df = read.csv("df_exp2_vdj1.csv")
df$Model=factor(df$Model,level=c("NetTCR-1.0","NetTCR-2.0","ERGO_LSTM","ERGO_AE","DLpTCR","HeteroTCR"))
df$based=factor(df$based,level=c("Pair-based","TCR-based","Antigen-based","Strict-based"))

p = ggplot(data=df,aes(x=Model,y=AUC,fill=factor(based)))+
    labs(x="",y="AUC",fill="Four types of data splitting methods")+
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
    scale_fill_manual(values=c('cornflowerblue','darkorange','forestgreen','firebrick'))+
    scale_y_continuous(limits=c(0.36,0.88),breaks=seq(0.36,0.88,0.08))+
    geom_signif(comparisons = list(c("HeteroTCR","DLpTCR"),
                                   c("HeteroTCR","ERGO_AE"),
                                   c("HeteroTCR","ERGO_LSTM"),
                                   c("HeteroTCR","NetTCR-2.0"),
                                   c("HeteroTCR","NetTCR-1.0")),
                test="t.test",
                test.args=list(paired=TRUE,alternative="greater",var.equal=FALSE),
                map_signif_level=TRUE,
                step_increase=0.1, size=0.3, textsize = 2.5)

ggsave("supplementaryFig1_20240501.pdf",p,device = "pdf",width=18,height=7.4,units="cm")


#vdj2
df = read.csv("df_exp2_vdj2.csv")
df$Model=factor(df$Model,level=c("NetTCR-1.0","NetTCR-2.0","ERGO_LSTM","ERGO_AE","DLpTCR","HeteroTCR"))
df$based=factor(df$based,level=c("Pair-based","TCR-based","Antigen-based","Strict-based"))

p = ggplot(data=df,aes(x=Model,y=AUC,fill=factor(based)))+
    labs(x="",y="AUC",fill="Four types of data splitting methods")+
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
    scale_fill_manual(values=c('cornflowerblue','darkorange','forestgreen','firebrick'))+
    scale_y_continuous(limits=c(0.1,1.3),breaks=seq(0.1,1.0,0.1))+
    geom_signif(comparisons = list(c("HeteroTCR","DLpTCR"),
                                   c("HeteroTCR","ERGO_AE"),
                                   c("HeteroTCR","ERGO_LSTM"),
                                   c("HeteroTCR","NetTCR-2.0"),
                                   c("HeteroTCR","NetTCR-1.0")),
                test="t.test",
                test.args=list(paired=TRUE,alternative="greater",var.equal=FALSE),
                map_signif_level=TRUE,
                step_increase=0.1, size=0.3, textsize = 2.5)

ggsave("supplementaryFig2_20240501.pdf",p,device = "pdf",width=18,height=7.4,units="cm")


#vdj3
 df = read.csv("df_exp2_vdj3.csv")
 df$Model=factor(df$Model,level=c("NetTCR-1.0","NetTCR-2.0","ERGO_LSTM","ERGO_AE","DLpTCR","HeteroTCR"))
df$based=factor(df$based,level=c("Pair-based","TCR-based","Antigen-based","Strict-based"))

p = ggplot(data=df,aes(x=Model,y=AUC,fill=factor(based)))+
    labs(x="",y="AUC",fill="Four types of data splitting methods")+
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
    scale_fill_manual(values=c('cornflowerblue','darkorange','forestgreen','firebrick'))+
    scale_y_continuous(limits=c(0.2,1.4),breaks=seq(0.2,1.0,0.1))+
    geom_signif(comparisons = list(c("HeteroTCR","DLpTCR"),
                                   c("HeteroTCR","ERGO_AE"),
                                   c("HeteroTCR","ERGO_LSTM"),
                                   c("HeteroTCR","NetTCR-2.0"),
                                   c("HeteroTCR","NetTCR-1.0")),
                test="t.test",
                test.args=list(paired=TRUE,alternative="greater",var.equal=FALSE),
                map_signif_level=TRUE,
                step_increase=0.1, size=0.3, textsize = 2.5)

ggsave("supplementaryFig3_20240501.pdf",p,device = "pdf",width=18,height=7.4,units="cm")

