

source("Simulation_experiments/simulation_experiment_England.R")
# palette = colorRampPalette(c('#1f77b4', '#d62728'))(100)
idx1 = (case.df$SAMPLE_DT_numeric > 40) & ( case.df$SAMPLE_DT_numeric < 60) & (case.df$warning.score < 0.95)
#colors = palette[as.integer(case.df[idx,]$warning.score * 100)+1]
idx2 = (case.df$SAMPLE_DT_numeric > 40) & ( case.df$SAMPLE_DT_numeric < 60) & (case.df$warning.score > 0.95)

# png(paste0('Fig/Fig_7AB', label.str,'.png'), width = 5 * 2.8 / 3 * 2, height = 4, units = 'in', res=600, pointsize = 15)
par(mfrow=c(1, 2), mar= c(4, 4, 0.5, 0.5) )

palette=inferno(length(40:59),begin=0.15, end=0.85, direction = -1)
color = ifelse(colSums(epi.matrix) > 0,palette,'white')


plot(40:59, colSums(epi.matrix), col=color, xlab='t', pch=20, ylab='no. epidemic cases')
mtext('A', side=3, padj=2, at=42, cex=1.5)

plot(case.df$SAMPLE_DT_numeric, case.df$warning.score, xlab='t', ylab=expression(italic(w)),
     col=ifelse(case.df$true_positive, "#d6272855", "#1f77b455"), pch=ifelse(case.df$true_positive, 20, s1))
legend('bottomright', c('epidemic','endemic'), col=c("#d6272855", "#1f77b455"), pch=c(20,1))


mtext('B', side=3, padj=2, at=95, cex=1.5)

# plot(c(41,59), range(case.df[case.df$true_positive,11:30]), col='white', xlab = expression(tau),
#      ylab=expression(italic(w)), axes = FALSE)
# for(i in 1:NROW(case.df[case.df$true_positive, 11:30])){
#   t1 = case.df[case.df$true_positive,][i,'SAMPLE_DT_numeric'] +1
#   id = t1 - 40 + 11
#   if (t1 <= 59){
#     lines(t1:59, case.df[case.df$true_positive,][i, id:30], col='#d62728')
#     points(t1:59, case.df[case.df$true_positive,][i, id:30], pch=20, cex=0.4, col='#d62728')
#   }
# }
# axis(2, at=c(0,0.25,0.5,0.75,1))
# axis(1, at = seq(40,59))
#mtext('C', side=3, padj=2, at=42, cex=1.5)
# dev.off()
