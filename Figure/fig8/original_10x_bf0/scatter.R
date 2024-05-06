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
            geom_point(size = 0.3,color='#3E4A7B') +
            stat_cor(data = Hetero_donor1, method = "spearman", label.x.npc = 0.7, size=2)+
            theme_bw()+
            theme(panel.grid=element_blank(),
                plot.title = element_text(size = 7, family = "Helvetica", hjust = 0.5),
                axis.title = element_text(size = 7,family = "Helvetica"),
                axis.text = element_text(size = 6, family = "Helvetica"))+
            scale_y_continuous(limits=c(0,1))
ggsave("fig8a.pdf",p_Hetero_UMI_donor1,device = "pdf",width=10,height=9,units="cm")

# Donor2
pred_donor2 = read.table('donor2/pred.tsv', header=TRUE, sep='\t')
test_donor2 = read.table('donor2/test.tsv', header=TRUE, sep='\t')
Hetero_donor2 = left_join(test_donor2, pred_donor2, by=c("cdr3"="cdr3", "peptide"="peptide"))
p_Hetero_UMI_donor2 = ggplot(Hetero_donor2, aes(x = fraction, y = probability)) +
            labs(x='Binding Fraction', y = 'Predicted Probability', title = 'Donor 2')+ 
            geom_point(size = 0.3,color='#3E4A7B') + 
            stat_cor(data = Hetero_donor2, method = "spearman", label.x.npc = 0.7, size=2)+
            theme_bw()+
            theme(panel.grid=element_blank(),
                plot.title = element_text(size = 7, family = "Helvetica", hjust = 0.5),
                axis.title = element_text(size = 7,family = "Helvetica"),
                axis.text = element_text(size = 6, family = "Helvetica"))+
            scale_y_continuous(limits=c(0,1))
ggsave("fig8b.pdf",p_Hetero_UMI_donor2,device = "pdf",width=10,height=9,units="cm")

# Donor3
pred_donor3 = read.table('donor3/pred.tsv', header=TRUE, sep='\t')
test_donor3 = read.table('donor3/test.tsv', header=TRUE, sep='\t')
Hetero_donor3 = left_join(test_donor3, pred_donor3, by=c("cdr3"="cdr3", "peptide"="peptide"))
p_Hetero_UMI_donor3 = ggplot(Hetero_donor3, aes(x = fraction, y = probability)) +
            labs(x='Binding Fraction', y = 'Predicted Probability', title = 'Donor 3')+ 
            geom_point(size = 0.3,color='#3E4A7B') +
            stat_cor(data = Hetero_donor3, method = "spearman", label.x.npc = 0.7, size=2)+
            theme_bw()+
            theme(panel.grid=element_blank(),
                plot.title = element_text(size = 7, family = "Helvetica", hjust = 0.5),
                axis.title = element_text(size = 7,family = "Helvetica"),
                axis.text = element_text(size = 6, family = "Helvetica"))+
            scale_y_continuous(limits=c(0,1))
ggsave("fig8c.pdf",p_Hetero_UMI_donor3,device = "pdf",width=10,height=9,units="cm")

# Donor4
pred_donor4 = read.table('donor4/pred.tsv', header=TRUE, sep='\t')
test_donor4 = read.table('donor4/test.tsv', header=TRUE, sep='\t')
Hetero_donor4 = left_join(test_donor4, pred_donor4, by=c("cdr3"="cdr3", "peptide"="peptide"))
p_Hetero_UMI_donor4 = ggplot(Hetero_donor4, aes(x = fraction, y = probability)) +
            labs(x='Binding Fraction', y = 'Predicted Probability', title = 'Donor 4')+ 
            geom_point(size = 0.3,color='#3E4A7B') + 
            stat_cor(data = Hetero_donor4, method = "spearman", label.x.npc = 0.7, size=2)+
            theme_bw()+
            theme(panel.grid=element_blank(),
                plot.title = element_text(size = 7, family = "Helvetica", hjust = 0.5),
                axis.title = element_text(size = 7,family = "Helvetica"),
                axis.text = element_text(size = 6, family = "Helvetica"))+
            scale_y_continuous(limits=c(0,1))
ggsave("fig8d.pdf",p_Hetero_UMI_donor4,device = "pdf",width=10,height=9,units="cm")