library(tidyverse)
library(R.matlab)

fpath   = './'
PngPath = './png/'
fname = 'RFPc_rawData.rds'
ds = readRDS(file = paste(fpath, fname, sep = ''))


vpnames = unique(ds$vp)

# get fraction 'spiky' responses  
prz1all = ds %>% group_by(vp, nodd, rfpbase) %>%
      summarize(f_spikey = mean(isRound == 0)) %>%
      filter(rfpbase == '0.5') %>% 
      ungroup()

prz2all = ds %>% group_by(vp, nodd, rfpbase) %>%
       summarize(f_spikey = mean(isRound == 0)) %>%
       filter(rfpbase == '1.5') %>% 
       ungroup()

for (k in 1:length(vpnames)){
PngFileName = paste(PngPath, 'RFPc-', vpnames[k], '-linear.png', sep = '')
png(filename = PngFileName,
    width  = 360,
    height = 360)

prz1  = prz1all %>% filter(vp == vpnames[k])  
spike = prz1$f_spikey
erg = lm(spike ~ prz1$nodd)
plot(spike ~ prz1$nodd, col = 'red', xlim=c(0,12), ylim=c(0,1), ylab = 'p spitz', xlab = 'n odd-symmetrics', lwd=2)
abline(erg, col= 'red', lwd=2)

prz2  = prz2all %>% filter(vp == vpnames[k])
spike = prz2$f_spikey
erg = lm(spike ~ prz2$nodd)
points(spike ~ prz1$nodd, col = 'blue', xlim=c(0,12), ylim=c(0,1), ylab = 'p spitz', xlab = 'n odd-symmetrics', lwd=2)
abline(erg, col= 'blue', lwd=2)

text(2.5,0.1, vpnames[k])
legend(1,0.4,legend = c('0.5', '1.5'), lty = c('solid', 'solid'),col = c('red', 'blue'), bty='n')
dev.off()
}

