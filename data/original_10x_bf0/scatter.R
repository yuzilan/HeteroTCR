library(ggpubr)
library(ggplot2)
library(ggsci)
library(ggsignif)
library(gghalves)
library(tidyverse)

# Donor1 new clone fraction
pred_donor1 = read.table('donor1/pred.tsv', header=TRUE, sep='\t')
test_donor1 = read.table('donor1/test.tsv', header=TRUE, sep='\t')
Hetero_donor1 = left_join(test_donor1, pred_donor1, by=c("cdr3"="cdr3", "peptide"="peptide"))
p_Hetero_UMI_donor1 = ggplot(Hetero_donor1, aes(x = fraction, y = probability)) +
            labs(x='Binding Fraction', y = 'Predicted Probability', title = 'Donor 1')+ 
            geom_point(size = 0.5,color='#3E4A7B') +
            stat_cor(data = Hetero_donor1, method = "spearman", label.x.npc = 0.55)+
            theme_bw()+
            theme(panel.grid=element_blank())+
            theme(plot.title = element_text(size = 16, family = "Arial", hjust = 0.5,face = "bold"),
                  plot.subtitle = element_text(size = 14, family = "Arial", hjust = 0.5,face = "bold"),
                  axis.title.x = element_text(size = 14,family = "Arial", face = "bold"),
                  axis.title.y = element_text(size = 14,family = "Arial",face = "bold"),
                  axis.text.x = element_text(size = 13, family = "Arial",face = "bold"),
                  axis.text.y = element_text(size = 13, family = "Arial",face = "bold"),
                  strip.text.x = element_text(size = 13, family = "Arial",face = "bold"),
                  legend.position = "None",
                  )
ggsave("donor1.png",p_Hetero_UMI_donor1,width=5,height=4.5,dpi=300)

# Donor2
pred_donor2 = read.table('donor2/pred.tsv', header=TRUE, sep='\t')
test_donor2 = read.table('donor2/test.tsv', header=TRUE, sep='\t')
Hetero_donor2 = left_join(test_donor2, pred_donor2, by=c("cdr3"="cdr3", "peptide"="peptide"))
p_Hetero_UMI_donor2 = ggplot(Hetero_donor2, aes(x = fraction, y = probability)) +
            labs(x='Binding Fraction', y = 'Predicted Probability', title = 'Donor 2')+ 
            geom_point(size = 0.5,color='#3E4A7B') + 
            stat_cor(data = Hetero_donor2, method = "spearman", label.x.npc = 0.55)+
            theme_bw()+
            theme(panel.grid=element_blank())+
            theme(plot.title = element_text(size = 16, family = "Arial", hjust = 0.5,face = "bold"),
                  plot.subtitle = element_text(size = 14, family = "Arial", hjust = 0.5,face = "bold"),
                  axis.title.x = element_text(size = 14,family = "Arial", face = "bold"),
                  axis.title.y = element_text(size = 14,family = "Arial",face = "bold"),
                  axis.text.x = element_text(size = 13, family = "Arial",face = "bold"),
                  axis.text.y = element_text(size = 13, family = "Arial",face = "bold"),
                  strip.text.x = element_text(size = 13, family = "Arial",face = "bold"),
                  legend.position = "None",
                  )
ggsave("donor2.png",p_Hetero_UMI_donor2,width=5,height=4.5,dpi=300)

# Donor3
pred_donor3 = read.table('donor3/pred.tsv', header=TRUE, sep='\t')
test_donor3 = read.table('donor3/test.tsv', header=TRUE, sep='\t')
Hetero_donor3 = left_join(test_donor3, pred_donor3, by=c("cdr3"="cdr3", "peptide"="peptide"))
Hetero_donor3$log_UMI = log10(Hetero_donor3$UMI)
p_Hetero_UMI_donor3 = ggplot(Hetero_donor3, aes(x = fraction, y = probability)) +
            labs(x='Binding Fraction', y = 'Predicted Probability', title = 'Donor 3')+ 
            geom_point(size = 0.5,color='#3E4A7B') +
            stat_cor(data = Hetero_donor3, method = "spearman", label.x.npc = 0.55)+
            theme_bw()+
            theme(panel.grid=element_blank())+
            theme(plot.title = element_text(size = 16, family = "Arial", hjust = 0.5,face = "bold"),
                  plot.subtitle = element_text(size = 14, family = "Arial", hjust = 0.5,face = "bold"),
                  axis.title.x = element_text(size = 14,family = "Arial", face = "bold"),
                  axis.title.y = element_text(size = 14,family = "Arial",face = "bold"),
                  axis.text.x = element_text(size = 13, family = "Arial",face = "bold"),
                  axis.text.y = element_text(size = 13, family = "Arial",face = "bold"),
                  strip.text.x = element_text(size = 13, family = "Arial",face = "bold"),
                  legend.position = "None",
                  )
ggsave("donor3.png",p_Hetero_UMI_donor3,width=5,height=4.5,dpi=300)

# Donor4
pred_donor4 = read.table('donor4/pred.tsv', header=TRUE, sep='\t')
test_donor4 = read.table('donor4/test.tsv', header=TRUE, sep='\t')
Hetero_donor4 = left_join(test_donor4, pred_donor4, by=c("cdr3"="cdr3", "peptide"="peptide"))
p_Hetero_UMI_donor4 = ggplot(Hetero_donor4, aes(x = fraction, y = probability)) +
            labs(x='Binding Fraction', y = 'Predicted Probability', title = 'Donor 4')+ 
            geom_point(size = 0.5,color='#3E4A7B') + 
            stat_cor(data = Hetero_donor4, method = "spearman", label.x.npc = 0.55)+
            theme_bw()+
            theme(panel.grid=element_blank())+
            theme(plot.title = element_text(size = 16, family = "Arial", hjust = 0.5,face = "bold"),
                  plot.subtitle = element_text(size = 14, family = "Arial", hjust = 0.5,face = "bold"),
                  axis.title.x = element_text(size = 14,family = "Arial", face = "bold"),
                  axis.title.y = element_text(size = 14,family = "Arial",face = "bold"),
                  axis.text.x = element_text(size = 13, family = "Arial",face = "bold"),
                  axis.text.y = element_text(size = 13, family = "Arial",face = "bold"),
                  strip.text.x = element_text(size = 13, family = "Arial",face = "bold"),
                  legend.position = "None",
                  )
ggsave("donor4.png",p_Hetero_UMI_donor4,width=5,height=4.5,dpi=300)