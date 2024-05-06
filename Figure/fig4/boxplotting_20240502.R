getwd()

library(ggplot2)
library(ggsignif)
library(ggpubr)
# library(scales)

# squash_axis <- function(from, to, factor) { 
#   trans <- function(x) {    
#     # get indices for the relevant regions
#     isq <- !is.na(x) & x > from & x < to
#     ito <- !is.na(x) & x >= to
#     
#     # apply transformation
#     x[isq] <- from + (x[isq] - from)/factor
#     x[ito] <- from + (to - from)/factor + (x[ito] - to)
#     
#     return(x)
#   }
#   
#   inv <- function(x) {
#     # get indices for the relevant regions
#     isq <- !is.na(x) & x > from & x < from + (to - from)/factor
#     ito <- !is.na(x) & x >= from + (to - from)/factor
#     
#     # apply inverse transformation
#     x[isq] <- from + (x[isq] - from) * factor
#     x[ito] <- to + (x[ito] - (from + (to - from)/factor))
#     
#     return(x)
#   }
#   
#   # return the transformation
#   return(trans_new("squash_axis", trans, inv))
# }


df = read.csv("df_exp4_vdj0.csv")
head(df)

df$Model=factor(df$Model,level=c("Baseline","HeteroTCR"))
df$based=factor(df$based,level=c("Pair-based","TCR-based","Antigen-based","Strict-based"))

p = ggplot(data=df,aes(x=based,y=AUC,fill=factor(Model)))+
    labs(x="",y="AUC",fill="Models")+
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
    scale_fill_manual(values=c('firebrick','cornflowerblue'))+
    scale_y_continuous(limits=c(0.46,0.78),breaks=seq(0.46,0.78,0.04))+
#     coord_trans(y = squash_axis(0.52, 0.64, 5))+
    geom_signif(annotations=c("***","***","***","***"),
                         y_position=c(0.76,0.76,0.76,0.76),
                         xmin=c(0.83,1.83,2.83,3.83),
                         xmax=c(1.17,2.17,3.17,4.17),
                         tip_length=c(0.03,0.03,0.03,0.03),
                         size=0.3, textsize = 2.5)

ggsave("Figure4_20240502.pdf",p,device = "pdf",width=18,height=7.4,units="cm")
