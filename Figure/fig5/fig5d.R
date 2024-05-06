getwd()

df = read.csv("df_exp4_vdj1.csv")

head(df)

library(ggplot2)
library(ggsignif)
library(ggpubr)

df$Model=factor(df$Model,level=c("Without Heterogeneous GNN Module","With Heterogeneous GNN Module"))
df$Peptide=factor(df$Peptide,level=c("GTSGSPIVNR","FLKEKGGL","GTSGSPIINR","KLVALGINAV","GLCTLVAML","FPRPWLHGL","GILGFVFTL","NLVPMVATV","KAFSPEVIPMF","TPQDLNTML","ELAGIGILTV","KRWIILGLNK"))


p = ggplot(data=df,aes(x=Peptide,y=AUC,fill=Model))+
#     stat_summary(geom="bar",fun=mean,aes(fill=Model),width=0.6,position=position_dodge(0.7),colour="black")+
#     stat_summary(geom="errorbar",fun.min=min,fun.max=max,width=0.3,position=position_dodge(0.7))+
    labs(x="",y="AUC",fill="Models")+
    stat_boxplot(geom="errorbar",width=0.3,position=position_dodge(0.7),linewidth=0.3)+ 
    geom_boxplot(position=position_dodge(0.7),width=0.6,alpha=0.5,outlier.size=1,linewidth=0.3)+
    geom_point(position = position_dodge(0.7), size = 1, shape = 21, stroke=0.3) +
    theme_bw()+
    theme(panel.grid=element_blank(),
        axis.text=element_text(family="Helvetica", size=6),
        axis.title=element_text(family="Helvetica", size=7),
        legend.text=element_text(family="Helvetica", size=6),
        legend.title=element_text(family="Helvetica", size=7),
        axis.text.x=element_text(angle=45,hjust = 0.5,vjust = 0.5))+
    scale_fill_manual(values=c('firebrick','cornflowerblue'))+
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.1))
    
ggsave("Figure5d_20240502.pdf",p,device = "pdf",width=18,height=7.4,units="cm")

