rm(list=ls())
library(ggplot2)
library(RColorBrewer)
load("~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_betas_with_freq.RData")

#The beta dist. for the hermaphrodites

i=1
plot(density(beta[i,],bw=.01),col=brewer.pal(n = 6, name = 'Dark2')[i],bty="n",las=1,xlab=expression(paste("Net directional selection gradient (", beta[net], ")")),main="",xlim=c(-.1,0.05),ylab="Density",ylim=c(0,30))
#mtext(side=2,"Density",padj=-2.4)
mtext(side=3,"Transition rates")
for(i in 2:6) lines(density(beta[i,],bw=.01),col=brewer.pal(n = 6, name = 'Dark2')[i])
abline(v=rowMedians(beta),col=brewer.pal(n = 6, name = 'Dark2'))
legend(-.1,28,c("SF","SB","FS","FB","BS","BF"),lwd=2,col=brewer.pal(n = 6, name = 'Dark2'))

#The beta dist. for the males

load("~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_betas_with_freq_Males.RData")

i=1
plot(density(beta[i,],bw=.01),col=brewer.pal(n = 6, name = 'Dark2')[i],bty="n",las=1,xlab=expression(paste("Net directional selection gradient (", beta[net], ")")),main="",xlim=c(-.1,0.05),ylab="Density",ylim=c(0,30))
#mtext(side=2,"Density",padj=-2.4)
mtext(side=3,"Transition rates")
for(i in 2:6) lines(density(beta[i,],bw=.01),col=brewer.pal(n = 6, name = 'Dark2')[i])
abline(v=rowMedians(beta),col=brewer.pal(n = 6, name = 'Dark2'))
legend(-.1,28,c("SF","SB","FS","FB","BS","BF"),lwd=2,col=brewer.pal(n = 6, name = 'Dark2'))


