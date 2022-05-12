
```{r}
# png(paste0('Fig/Fig_7DE', label.str,'.png'), width = 5 * 2.8 / 3, height = 4 * 2, units = 'in', res=600, pointsize = 15)
par(mfrow=c(2, 1), mar= c(4, 4, 0.5, 0.5) )
library(viridisLite)
palette=inferno(length(40:59),begin=0.15, end=0.85, direction=-1)

plot(case.df$warning.score, case.df$warning.score2,
     xlab = expression(paste(italic(w),' (rep. 1)')), ylab = expression(paste(italic(w)," (rep. 2)")), col=colors, pch=ifelse(case.df$true_positive, 20,1))
legend('bottomright', c('epidemic','endemic'), col=c("#d6272855", "#1f77b455"), pch=c(20,1))

mtext('D', side=3, padj=2, at=(42-40)/(58.5-40), cex=1.5)

source("Simulation_experiments/simulation_experiment_England.R")
#png(paste0('Fig/timelines_', label.str, '.png'), width = 4 * 1.2, height = 3, units = 'in', res=400, pointsize = 13)
#par(mar=c(4,4,1,1)+0.1)

plot(c(40,59), range(case.df[case.df$true_positive,11:30]), col='white', xlab = expression(tau),
     ylab=expression(italic(w)), axes = FALSE)
#abline(h=0.95, col='#7f7f7f')

for(i in 1:NROW(case.df[case.df$true_positive,  ])){
  t1 = case.df[case.df$true_positive,][i,'SAMPLE_DT_numeric']
  id = t1 - 40 + 18
  # 40 is the week of the start of the outbreak
  # 18 is the number of columns before the column `warning.score40`
  color = palette[t1-39]
  if (t1 <= 59){
    #  print(c(dim(case.df[case.df$true_positive,]),  i, id:31))
    #    print( unlist(case.df[case.df$true_positive,][i, id:31]))
    lines(t1:59, case.df[case.df$true_positive,][i, id:37], col=color)
    points(t1:59, case.df[case.df$true_positive,][i, id:37], pch=20, cex=0.4, col=color)
  }
}
axis(2, at=c(0,0.25,0.5,0.75,1))
axis(1, at = seq(40,59))
mtext('E', side=3, padj=3, at=42, cex=1.5)

# dev.off()
```