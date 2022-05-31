library(tidyverse)
library(R.matlab)


flist = sort(dir('./raw/', pattern = 'S*_*.mat', full.names=T))

tmp = data.frame()
for (k in 1:length(flist)){
df    = readMat(flist[k]) 
tmp2 = data.frame(df$RFP[[5]], df$RFP[[1]], df$RFP[[2]], df$RFP[[3]], df$RFP[[4]])
tmp = rbind(tmp, tmp2)
}

ds = as_tibble(tmp) %>%
     mutate(vp       = as_factor(df.RFP..1..),
            mapping  = fct_recode(as_factor(df.RFP..2..), rund_links = '1', rund_rechts = '2'),
            date     = as.POSIXct(df.RFP..3.., format='%Y%m%dT%H%M%S')) %>%
     rename(trial       = X1,
            block       = X2,
            origTrial   = X3,
            nodd        = X4,
            rfpbase     = X5,
            rfpfreq     = X6,
            audiofreq   = X7,
            phi         = X8,
            rt          = X9,
            keycode     = X10,
            isRound     = X11,
            ISI         = X12,
            size        = X13) %>%
    select(vp, block, trial, rt, isRound, rfpbase,
           rfpfreq, nodd, audiofreq , size,
           origTrial, mapping, date)

ds = ds %>%
     filter(!vp %in% c('S02', 'S10')) # exclude bad data 

# S19 had one additional block, remove block 13
ds[ds$vp=='S19',]$block = rep(1:13, each = 48) 
ds  = ds %>% filter(block < 13) 


fname = 'RFPc_rawData.rds'
saveRDS(ds, file = fname)
