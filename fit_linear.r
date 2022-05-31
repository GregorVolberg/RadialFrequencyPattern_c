library(tidyverse)

source('./RFPc_functions.r')

fname <- 'RFPc_rawData.rds'
ds    <- readRDS(fname) 

# sanity check
n <- n_distinct(ds$vp)
trialsPerPerson <- ds %>%
  group_by(vp) %>%
  summarize(nperP = length(vp)) %>%
  ungroup()
cat('\n******\n', n, 'participants in data set\n')


k1 = ds %>% group_by(vp, nodd, rfpbase) %>%
        summarize(N = n(), pSpitz = sum(isRound == 0)/n()) %>%
        filter(rfpbase == 0.5) %>%
        ungroup() %>%
        group_by(vp) %>%
        group_map(~ lm(pSpitz~nodd, data=.)) %>%
        map_df('coefficients') %>%
        rename(rfp05b0 = '(Intercept)',
               rfp05b1 = 'nodd')
k2 = ds %>% group_by(vp, nodd, rfpbase) %>%
        summarize(N = n(), pSpitz = sum(isRound == 0)/n()) %>%
        filter(rfpbase == 1.5) %>%
        ungroup() %>%
        group_by(vp) %>%
        group_map(~ lm(pSpitz~nodd, data=.)) %>%
        map_df('coefficients') %>%
        rename(rfp15b0 = '(Intercept)',
               rfp15b1 = 'nodd')

erg = cbind(vp = unique(ds$vp), k1,k2)         

# export raw data
write_excel_csv2(erg, 'fittedParameters_RFPc.csv')

dsFull = ds %>% group_by(vp, nodd, rfpbase, audiofreq) %>%
  summarize(N = n(),
            p = mean(isRound == 0),
            RT = mean(rt)) %>%
  ungroup()

write_excel_csv2(dsFull, 'fullData.csv')


