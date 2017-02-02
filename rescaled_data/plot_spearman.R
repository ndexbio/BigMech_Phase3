library('ggplot2')

prior_data=read.table('prior_spearman_normalized_results_rescaled_unix.txt',header=TRUE)
high_conf_data=read.table('high_confidence_spearman_results.txt',header=TRUE)

drugs<-c('AK','Tm','ZS','901','ST','HN','PLX','RO','RY','P6','SR','NT')

for (d in drugs){
M=mean(prior_data$Spearman[grepl(d,prior_data$Drug)])
print(paste(d," average Spearman ",M))
}

to_plot=data.frame(corr=c(high_conf_data$Spearman,prior_data$Spearman),
model=c(rep('prior plus literature',
length(high_conf_data$Spearman)),rep('prior',
length(prior_data$Spearman))))

m<-ggplot(to_plot,aes(corr,color=model)) + geom_density()

ggsave(m,file='Spearman_corr.png')