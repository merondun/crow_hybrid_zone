# ~~ CGI ~~

## CG 5x

```
library(methylKit)
library(lme4)
library(Matrix)
library(emmeans)
library(multcomp)
library(matrixStats)
library(ggplot2)
library(gridExtra)
library(HDInterval)
library(dplyr)
library(stringr)

setwd("E:/Research/scratch/crow_hybrid_paper/CGI/CG/dplyr/")

#import all the data and average by CGI
tab <- read.table('D_Ko_C29_BL_ADL_F.CGI.bed',header=FALSE);  tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);C29_ADL_F <- tab2[,c(7,1,2,3,8,9,10)]; names(C29_ADL_F) <- c('site','chr','start','end','numCs1','numTs1','Sites1');C29_ADL_F$C29_ADL_F <- (C29_ADL_F$numCs1/(C29_ADL_F$numCs1+C29_ADL_F$numTs1))
tab <- read.table('D_Ko_C29_BL_CHK_F.CGI.bed',header=FALSE);  tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);C29_CHK_F <- tab2[,c(7,1,2,3,8,9,10)]; names(C29_CHK_F) <- c('site','chr','start','end','numCs2','numTs2','Sites2');C29_CHK_F$C29_CHK_F <- (C29_CHK_F$numCs2/(C29_CHK_F$numCs2+C29_CHK_F$numTs2))
tab <- read.table('D_Ko_C29_BL_YRL_F.CGI.bed',header=FALSE);  tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);C29_YRL_F <- tab2[,c(7,1,2,3,8,9,10)]; names(C29_YRL_F) <- c('site','chr','start','end','numCs3','numTs3','Sites3');C29_YRL_F$C29_YRL_F <- (C29_YRL_F$numCs3/(C29_YRL_F$numCs3+C29_YRL_F$numTs3))
tab <- read.table('D_Ko_C31_BL_ADL_M.CGI.bed',header=FALSE);  tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);C31_ADL_M <- tab2[,c(7,1,2,3,8,9,10)]; names(C31_ADL_M) <- c('site','chr','start','end','numCs4','numTs4','Sites4');C31_ADL_M$C31_ADL_M <- (C31_ADL_M$numCs4/(C31_ADL_M$numCs4+C31_ADL_M$numTs4))
tab <- read.table('D_Ko_C31_BL_YRL_M.CGI.bed',header=FALSE);  tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);C31_YRL_M <- tab2[,c(7,1,2,3,8,9,10)]; names(C31_YRL_M) <- c('site','chr','start','end','numCs5','numTs5','Sites5');C31_YRL_M$C31_YRL_M <- (C31_YRL_M$numCs5/(C31_YRL_M$numCs5+C31_YRL_M$numTs5))
tab <- read.table('D_Ko_C36_BL_ADL_M.CGI.bed',header=FALSE);  tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);C36_ADL_M <- tab2[,c(7,1,2,3,8,9,10)]; names(C36_ADL_M) <- c('site','chr','start','end','numCs6','numTs6','Sites6');C36_ADL_M$C36_ADL_M <- (C36_ADL_M$numCs6/(C36_ADL_M$numCs6+C36_ADL_M$numTs6))
tab <- read.table('D_Ko_C36_BL_CHK_M.CGI.bed',header=FALSE);  tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);C36_CHK_M <- tab2[,c(7,1,2,3,8,9,10)]; names(C36_CHK_M) <- c('site','chr','start','end','numCs7','numTs7','Sites7');C36_CHK_M$C36_CHK_M <- (C36_CHK_M$numCs7/(C36_CHK_M$numCs7+C36_CHK_M$numTs7))
tab <- read.table('D_Ko_C36_BL_YRL_M.CGI.bed',header=FALSE);  tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);C36_YRL_M <- tab2[,c(7,1,2,3,8,9,10)]; names(C36_YRL_M) <- c('site','chr','start','end','numCs8','numTs8','Sites8');C36_YRL_M$C36_YRL_M <- (C36_YRL_M$numCs8/(C36_YRL_M$numCs8+C36_YRL_M$numTs8))
tab <- read.table('D_Ko_C45_BL_ADL_F.CGI.bed',header=FALSE);  tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);C45_ADL_F <- tab2[,c(7,1,2,3,8,9,10)]; names(C45_ADL_F) <- c('site','chr','start','end','numCs9','numTs9','Sites9');C45_ADL_F$C45_ADL_F <- (C45_ADL_F$numCs9/(C45_ADL_F$numCs9+C45_ADL_F$numTs9))
tab <- read.table('D_Ko_C45_BL_CHK_F.CGI.bed',header=FALSE);  tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);C45_CHK_F <- tab2[,c(7,1,2,3,8,9,10)]; names(C45_CHK_F) <- c('site','chr','start','end','numCs10','numTs10','Sites10');C45_CHK_F$C45_CHK_F <- (C45_CHK_F$numCs10/(C45_CHK_F$numCs10+C45_CHK_F$numTs10))
tab <- read.table('D_Ko_C45_BL_YRL_F.CGI.bed',header=FALSE);  tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);C45_YRL_F <- tab2[,c(7,1,2,3,8,9,10)]; names(C45_YRL_F) <- c('site','chr','start','end','numCs11','numTs11','Sites11');C45_YRL_F$C45_YRL_F <- (C45_YRL_F$numCs11/(C45_YRL_F$numCs11+C45_YRL_F$numTs11))
tab <- read.table('S_Up_H59_BL_ADL_M.CGI.bed',header=FALSE);  tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);H59_ADL_M <- tab2[,c(7,1,2,3,8,9,10)]; names(H59_ADL_M) <- c('site','chr','start','end','numCs12','numTs12','Sites12');H59_ADL_M$H59_ADL_M <- (H59_ADL_M$numCs12/(H59_ADL_M$numCs12+H59_ADL_M$numTs12))
tab <- read.table('S_Up_H59_BL_CHK_M.CGI.bed',header=FALSE);  tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);H59_CHK_M <- tab2[,c(7,1,2,3,8,9,10)]; names(H59_CHK_M) <- c('site','chr','start','end','numCs13','numTs13','Sites13');H59_CHK_M$H59_CHK_M <- (H59_CHK_M$numCs13/(H59_CHK_M$numCs13+H59_CHK_M$numTs13))
tab <- read.table('S_Up_H59_BL_YRL_M.CGI.bed',header=FALSE);  tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);H59_YRL_M <- tab2[,c(7,1,2,3,8,9,10)]; names(H59_YRL_M) <- c('site','chr','start','end','numCs14','numTs14','Sites14');H59_YRL_M$H59_YRL_M <- (H59_YRL_M$numCs14/(H59_YRL_M$numCs14+H59_YRL_M$numTs14))
tab <- read.table('S_Up_H60_BL_ADL_F.CGI.bed',header=FALSE);  tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);H60_ADL_F <- tab2[,c(7,1,2,3,8,9,10)]; names(H60_ADL_F) <- c('site','chr','start','end','numCs15','numTs15','Sites15');H60_ADL_F$H60_ADL_F <- (H60_ADL_F$numCs15/(H60_ADL_F$numCs15+H60_ADL_F$numTs15))
tab <- read.table('S_Up_H60_BL_CHK_F.CGI.bed',header=FALSE);  tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);H60_CHK_F <- tab2[,c(7,1,2,3,8,9,10)]; names(H60_CHK_F) <- c('site','chr','start','end','numCs16','numTs16','Sites16');H60_CHK_F$H60_CHK_F <- (H60_CHK_F$numCs16/(H60_CHK_F$numCs16+H60_CHK_F$numTs16))
tab <- read.table('S_Up_H60_BL_YRL_F.CGI.bed',header=FALSE);  tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);H60_YRL_F <- tab2[,c(7,1,2,3,8,9,10)]; names(H60_YRL_F) <- c('site','chr','start','end','numCs17','numTs17','Sites17');H60_YRL_F$H60_YRL_F <- (H60_YRL_F$numCs17/(H60_YRL_F$numCs17+H60_YRL_F$numTs17))
tab <- read.table('S_Up_H65_BL_ADL_F.CGI.bed',header=FALSE);  tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);H65_ADL_F <- tab2[,c(7,1,2,3,8,9,10)]; names(H65_ADL_F) <- c('site','chr','start','end','numCs18','numTs18','Sites18');H65_ADL_F$H65_ADL_F <- (H65_ADL_F$numCs18/(H65_ADL_F$numCs18+H65_ADL_F$numTs18))
tab <- read.table('S_Up_H65_BL_CHK_F.CGI.bed',header=FALSE);  tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);H65_CHK_F <- tab2[,c(7,1,2,3,8,9,10)]; names(H65_CHK_F) <- c('site','chr','start','end','numCs19','numTs19','Sites19');H65_CHK_F$H65_CHK_F <- (H65_CHK_F$numCs19/(H65_CHK_F$numCs19+H65_CHK_F$numTs19))
tab <- read.table('S_Up_H65_BL_YRL_F.CGI.bed',header=FALSE);  tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);H65_YRL_F <- tab2[,c(7,1,2,3,8,9,10)]; names(H65_YRL_F) <- c('site','chr','start','end','numCs20','numTs20','Sites20');H65_YRL_F$H65_YRL_F <- (H65_YRL_F$numCs20/(H65_YRL_F$numCs20+H65_YRL_F$numTs20))
tab <- read.table('S_Up_H75_BL_ADL_M.CGI.bed',header=FALSE);  tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);H75_ADL_M <- tab2[,c(7,1,2,3,8,9,10)]; names(H75_ADL_M) <- c('site','chr','start','end','numCs21','numTs21','Sites21');H75_ADL_M$H75_ADL_M <- (H75_ADL_M$numCs21/(H75_ADL_M$numCs21+H75_ADL_M$numTs21))
tab <- read.table('S_Up_H75_BL_CHK_M.CGI.bed',header=FALSE);  tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);H75_CHK_M <- tab2[,c(7,1,2,3,8,9,10)]; names(H75_CHK_M) <- c('site','chr','start','end','numCs22','numTs22','Sites22');H75_CHK_M$H75_CHK_M <- (H75_CHK_M$numCs22/(H75_CHK_M$numCs22+H75_CHK_M$numTs22))
tab <- read.table('S_Up_H75_BL_YRL_M.CGI.bed',header=FALSE);  tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);H75_YRL_M <- tab2[,c(7,1,2,3,8,9,10)]; names(H75_YRL_M) <- c('site','chr','start','end','numCs23','numTs23','Sites23');H75_YRL_M$H75_YRL_M <- (H75_YRL_M$numCs23/(H75_YRL_M$numCs23+H75_YRL_M$numTs23))

head(H75_YRL_M)

master <- Reduce(function(x,y) merge(x = x, y = y, by=c('site','chr','start','end')),
                 list(C29_ADL_F,C29_CHK_F,C29_YRL_F,C31_ADL_M,C31_YRL_M,C36_ADL_M,C36_CHK_M,C36_YRL_M,C45_ADL_F,C45_CHK_F,C45_YRL_F,H59_ADL_M,H59_CHK_M,H59_YRL_M,H60_ADL_F,H60_CHK_F,H60_YRL_F,H65_ADL_F,H65_CHK_F,H65_YRL_F,H75_ADL_M,H75_CHK_M,H75_YRL_M))
head(master)

#calculate coverage and CpGs
cov <- master %>% mutate(minsite1 = ifelse(Sites1 < 3, NA, Sites1),
                         minsite2 = ifelse(Sites2 < 3, NA, Sites2),
                         minsite3 = ifelse(Sites3 < 3, NA, Sites3),
                         minsite4 = ifelse(Sites4 < 3, NA, Sites4),
                         minsite5 = ifelse(Sites5 < 3, NA, Sites5),
                         minsite6 = ifelse(Sites6 < 3, NA, Sites6),
                         minsite7 = ifelse(Sites7 < 3, NA, Sites7),
                         minsite8 = ifelse(Sites8 < 3, NA, Sites8),
                         minsite9 = ifelse(Sites9 < 3, NA, Sites9),
                         minsite10 = ifelse(Sites10 < 3, NA, Sites10),
                         minsite11 = ifelse(Sites11 < 3, NA, Sites11),
                         minsite12 = ifelse(Sites12 < 3, NA, Sites12),
                         minsite13 = ifelse(Sites13 < 3, NA, Sites13),
                         minsite14 = ifelse(Sites14 < 3, NA, Sites14),
                         minsite15 = ifelse(Sites15 < 3, NA, Sites15),
                         minsite16 = ifelse(Sites16 < 3, NA, Sites16),
                         minsite17 = ifelse(Sites17 < 3, NA, Sites17),
                         minsite18 = ifelse(Sites18 < 3, NA, Sites18),
                         minsite19 = ifelse(Sites19 < 3, NA, Sites19),
                         minsite20 = ifelse(Sites20 < 3, NA, Sites20),
                         minsite21 = ifelse(Sites21 < 3, NA, Sites21),
                         minsite22 = ifelse(Sites22 < 3, NA, Sites22),
                         minsite23 = ifelse(Sites23 < 3, NA, Sites23))

#count NAs
colSums(is.na(cov))

#filter according to missingness
sitedat <- cov[rowSums(is.na(cov[grepl('^minsite', names(cov))])) <= 4, ]
nrow(cov)
nrow(sitedat)

master <- sitedat

#define groups
master$MALE <- rowMeans(master[grepl('_M$', names(master))],na.rm=TRUE) #0=F  
master$FEMALE <- rowMeans(master[grepl('_F$', names(master))],na.rm=TRUE) #1=M
master$ADULT <- rowMeans(master[grepl('_ADL|_YRL', names(master))],na.rm=TRUE) #0=A
master$CHICK <- rowMeans(master[grepl('_CHK_', names(master))],na.rm=TRUE) #1=C
master$HOODED <- rowMeans(master[grepl('^H', names(master))],na.rm=TRUE) #1=C
master$CARRION <- rowMeans(master[grepl('^C', names(master))],na.rm=TRUE) #0=H
compare <- master[,c(1,(ncol(master)-4-1):ncol(master))]
head(compare)
#now we must do pairwise comparisons between our 3 treatments for % 5mC
compare[8] <- compare[2]-compare[3]
compare[9] <- compare[4]-compare[5]
compare[10] <- compare[6]-compare[7]
head(compare)
win_diff <- compare[,c(8,9,10)]

#merge the total percent diff with the original frame
master_meth <- cbind(master,win_diff)
colnames(master_meth)[colnames(master_meth) == 'MALE.1'] <- 'sex_M.F'
colnames(master_meth)[colnames(master_meth) == 'ADULT.1'] <- 'stage_ADL.CHK'
colnames(master_meth)[colnames(master_meth) == 'HOODED.1'] <- 'species_H.C'
nrow(master_meth)
names(master_meth)


#now subset for GLMM # C=1 H=0 / F=0 M=1 / ADL+YRL=0 CHK=1
IND1 <- as.data.frame(c(master_meth[1],"C29",master_meth[5], master_meth[6],"1","0","0"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Species","Sex","Stage"))
IND2 <- as.data.frame(c(master_meth[1],"C29",master_meth[9], master_meth[10],"1","0","1"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Species","Sex","Stage"))
IND3 <- as.data.frame(c(master_meth[1],"C29",master_meth[13], master_meth[14],"1","0","0"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Species","Sex","Stage"))
IND4 <- as.data.frame(c(master_meth[1],"C31",master_meth[17], master_meth[18],"1","1","0"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Species","Sex","Stage"))
IND5 <- as.data.frame(c(master_meth[1],"C31",master_meth[21], master_meth[22],"1","1","0"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Species","Sex","Stage"))
IND6 <- as.data.frame(c(master_meth[1],"C36",master_meth[25], master_meth[26],"1","1","0"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Species","Sex","Stage"))
IND7 <- as.data.frame(c(master_meth[1],"C36",master_meth[29], master_meth[30],"1","1","1"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Species","Sex","Stage"))
IND8 <- as.data.frame(c(master_meth[1],"C36",master_meth[33], master_meth[34],"1","1","0"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Species","Sex","Stage"))
IND9 <- as.data.frame(c(master_meth[1],"C45",master_meth[37], master_meth[38],"1","0","0"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Species","Sex","Stage"))
IND10 <- as.data.frame(c(master_meth[1],"C45",master_meth[41], master_meth[42],"1","0","1"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Species","Sex","Stage"))
IND11 <- as.data.frame(c(master_meth[1],"C45",master_meth[45], master_meth[46],"1","0","0"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Species","Sex","Stage"))
IND12 <- as.data.frame(c(master_meth[1],"H59",master_meth[49], master_meth[50],"0","1","0"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Species","Sex","Stage"))
IND13 <- as.data.frame(c(master_meth[1],"H59",master_meth[53], master_meth[54],"0","1","1"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Species","Sex","Stage"))
IND14 <- as.data.frame(c(master_meth[1],"H59",master_meth[57], master_meth[58],"0","1","0"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Species","Sex","Stage"))
IND15 <- as.data.frame(c(master_meth[1],"H60",master_meth[61], master_meth[62],"0","0","0"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Species","Sex","Stage"))
IND16 <- as.data.frame(c(master_meth[1],"H60",master_meth[65], master_meth[66],"0","0","1"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Species","Sex","Stage"))
IND17 <- as.data.frame(c(master_meth[1],"H60",master_meth[69], master_meth[70],"0","0","0"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Species","Sex","Stage"))
IND18 <- as.data.frame(c(master_meth[1],"H65",master_meth[73], master_meth[74],"0","0","0"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Species","Sex","Stage"))
IND19 <- as.data.frame(c(master_meth[1],"H65",master_meth[77], master_meth[78],"0","0","1"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Species","Sex","Stage"))
IND20 <- as.data.frame(c(master_meth[1],"H65",master_meth[81], master_meth[82],"0","0","0"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Species","Sex","Stage"))
IND21 <- as.data.frame(c(master_meth[1],"H75",master_meth[85], master_meth[86],"0","1","0"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Species","Sex","Stage"))
IND22 <- as.data.frame(c(master_meth[1],"H75",master_meth[89], master_meth[90],"0","1","1"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Species","Sex","Stage"))
IND23 <- as.data.frame(c(master_meth[1],"H75",master_meth[93], master_meth[94],"0","1","0"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Species","Sex","Stage"))

head(IND4)
head(IND23)

#Bind them all together and order it by window
datmeth <- rbind(IND1,IND2,IND3,IND4,IND5,IND6,IND7,IND8,IND9,IND10,IND11,IND12,IND13,IND14,IND15,IND16,IND17,IND18,IND19,IND20,IND21,IND22,IND23)
datmeth <- datmeth[order(datmeth$CGI),]

head(datmeth)

#sanity check 
nrow(datmeth)/23

#save some files we might need later 
write.csv(datmeth,file="CG_GLMM_Input_5x.csv",row.names=FALSE)
write.csv(master_meth,file="CG_FULL_5x.csv",row.names=FALSE)

#now we can start here 
datmeth <- read.csv("CG_GLMM_Input_5x.csv",header=TRUE)
head(datmeth)
nrow(datmeth)/23
#set these as factor variables 
datmeth$Species <- as.factor(datmeth$Species)
datmeth$Stage <- as.factor(datmeth$Stage)
datmeth$Sex <- as.factor(datmeth$Sex)
datmeth$ID <- as.factor(datmeth$ID)

##### MODEL #####
DMA_GLMM_pValues <- NULL
DMA_GLMM_Estimates <- NULL
DMA_GLMM_SEs <- NULL
DMA_GLMM_RandomEffect_Var <- NULL
DMA_GLMM_DispersionStat <- NULL
DMA_GLMM_Singular <- NULL

for (i in seq(1, nrow(datmeth), by = 23)) {
  tryCatch({
    cat('\n Running model for CpG site: ',((i-1)/23))
    # model each site
    CpGsite <- datmeth[i:(i + 22), ]
    # run the model and get the contrasts between sex, stage, species
    glmm_CpG <- glmer(cbind(MethCounts, UnMethCounts) ~
                        Sex + Stage + Species + (1 | ID) ,data = CpGsite, family = "binomial",
                      control = glmerControl(optimizer = "bobyqa", boundary.tol = 0.01,
                                             optCtrl = list(maxfun = 2e+08)))
    DMA_GLMM_Singularf <- isSingular(glmm_CpG)
    sum_mod <- summary(glmm_CpG)
    
    ## 1 random effects ( (1 | ID) )
    GLMM_RandomEffect_Var <- as.data.frame(VarCorr(glmm_CpG))$sdcor
    
    # calculate dispersion statistic (following AF Zuur, JM Hilbe
    # & EN Ieno. A beginner’s guide to GLM and GLMM with R - a
    # frequentist and Bayesian perspective for
    # ecologists(Highland Statistics Ltd., 2013).)
    E1 <- residuals(glmm_CpG)
    # number of parameters: fixed effects + 1 random effect
    p1 <- length(fixef(glmm_CpG)) + 1
    GLMM_DispersionStat <- sum(E1^2)/(nrow(CpGsite) - p1)
    # format estimates to add estimates for each loop to
    # estimates from previous loops
    
    # TO GET THE COLUMN HEADERS WITH THE CONTRASTS, check the excel file. First I look at the glmm_CpG_contrast object, which lists the models in order of 1-8, for this comparison. From here, make the 28 pairwise comparison charts, and substitute the 0 and 1 with your variables
    GLMM_Estimates_f <- data.frame(site = as.character(CpGsite[1,1]),Sex1 = sum_mod$coefficients[2,1], Stage1 = sum_mod$coefficients[3,1], Species1 = sum_mod$coefficients[4,1])
    GLMM_SEs_f <- data.frame(site = as.character(CpGsite[1,1]),Sex1 = sum_mod$coefficients[2,2], Stage1 = sum_mod$coefficients[3,2], Species1 = sum_mod$coefficients[4,2])
    GLMM_pValues_f <- data.frame(site = as.character(CpGsite[1,1]),Sex1 = sum_mod$coefficients[2,4], Stage1 = sum_mod$coefficients[3,4], Species1 = sum_mod$coefficients[4,4])
    GLMM_RandomEffect_Var_f <- data.frame(site = as.character(CpGsite[1,
                                                                      1]), ID = GLMM_RandomEffect_Var[1])
    GLMM_DispersionStat_f <- data.frame(site = as.character(CpGsite[1,
                                                                    1]), DispStat = GLMM_DispersionStat)
    # add estimates for each loop to estimates from previous
    # loops
    DMA_GLMM_pValues <- rbind(DMA_GLMM_pValues, GLMM_pValues_f)
    DMA_GLMM_Estimates <- rbind(DMA_GLMM_Estimates, GLMM_Estimates_f)
    DMA_GLMM_SEs <- rbind(DMA_GLMM_SEs, GLMM_SEs_f)
    DMA_GLMM_RandomEffect_Var <- rbind(DMA_GLMM_RandomEffect_Var,GLMM_RandomEffect_Var_f)
    DMA_GLMM_DispersionStat <- rbind(DMA_GLMM_DispersionStat,GLMM_DispersionStat_f)
    DMA_GLMM_Singular <- rbind(DMA_GLMM_Singular,DMA_GLMM_Singularf)
    
    #Save errors to this output
  }, error=function(e){cat(unique(CpGsite[,1]), "\n")
    
  }
  )
}

# 3. Save loop results
save(DMA_GLMM_pValues, DMA_GLMM_Estimates, DMA_GLMM_SEs, DMA_GLMM_RandomEffect_Var,
     DMA_GLMM_DispersionStat,DMA_GLMM_Singular, file = "CG_DMA_CrowResults_5x.RData")

#should be the samem otherwise we have models that didn't converge 
nrow(DMA_GLMM_pValues)
nrow(datmeth)/23
```

### Volcanos

```#
library(methylKit)
library(lme4)
library(Matrix)
library(emmeans)
library(multcomp)
library(tidyr)
library(ggplot2)
library(psych)
library(gridExtra)
library(HDInterval)
library(dplyr)
library(data.table)
library(stringr)
library(LICORS)
library(matrixStats)
library(viridis)
library(RColorBrewer)

setwd("E:/Research/scratch/crow_hybrid_paper/CGI/CG/dplyr")

load("CG_DMA_CrowResults_5x.RData")

#check out some of the data
head(DMA_GLMM_Estimates)
dim(DMA_GLMM_Estimates)
head(DMA_GLMM_pValues)
head(DMA_GLMM_SEs)
head(DMA_GLMM_Singular)


#let's check out the distributio of 5mC contrasts
multi.hist(DMA_GLMM_Estimates[,2:4],nrow=1,ncol=NULL,breaks=21)

#plot dispersion boxplot, we will exclude any sites outside the 99% interval
disp_limits <- hdi(DMA_GLMM_DispersionStat$DispStat,credMass = 0.99)
disp_limits

# lower     upper 
# 0.0251642 1.3674338 

jpeg("plots/CG_dispersion_plot.jpg",units="in",res=300, height=6, width=5)
theme_set(theme_classic(base_size = 16))
dplot <- ggplot(DMA_GLMM_DispersionStat,aes(y=DispStat,x=0))+geom_violin()+theme_classic(base_size = 16)+
  geom_jitter(width=0.1)+
  geom_text(aes(y = (disp_limits[2]+0.05),x=0.2, label = "Upper Limit"))+
  geom_hline(yintercept = disp_limits[2],linetype="dashed",col="red") +
  geom_text(aes(y = (disp_limits[1]+0.05),x=0.2, label = "Lower Limit"))+
  geom_hline(yintercept = disp_limits[1],linetype="dashed",col="red") +
  xlab("")+ylab("Dispersion Statistic")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
grid.arrange(dplot)
dev.off()

#subset bad sites
bad_sites <- subset(DMA_GLMM_DispersionStat,DispStat < disp_limits[1] | DispStat > disp_limits[2])
nrow(bad_sites)

# calculate 10% differences, now we care about direction so we re-do this, since before it was absolute value 
meth <- read.csv("CG_FULL_5x.csv",header=TRUE)
head(meth)

names(DMA_GLMM_pValues) <- c('site','pSex','pStage','pSpecies')
names(DMA_GLMM_Estimates) <- c('site','eSex','eStage','eSpecies')
names(DMA_GLMM_SEs) <- c('site','sSex','sStage','sSpecies')
DMA_RF_S <- cbind(DMA_GLMM_RandomEffect_Var,DMA_GLMM_Singular)
names(DMA_RF_S) <- c('site','RandomEffect','IsSingular')
master_meth1 <- merge(meth, DMA_GLMM_pValues, by='site')
master_meth2 <- merge(master_meth1, DMA_GLMM_Estimates, by='site')
master_meth3 <- merge(master_meth2, DMA_GLMM_SEs, by='site')
master_meth4 <- merge(master_meth3, DMA_GLMM_DispersionStat, by='site')
master_meth <- merge(master_meth4, DMA_RF_S, by='site')
head(master_meth)
nrow(master_meth)

#### Volcano with bonferonni
pthresh <- (0.05/nrow(master_meth))
master_meth <- master_meth %>% mutate(OD = ifelse(DispStat < disp_limits[1] | DispStat > disp_limits[2], "OVERDISPERSED","OKAY"),
                                      SING = ifelse(IsSingular == FALSE, "OKAY","SINGULAR"),
                                      DMR = ifelse(pSex < pthresh & abs(sex_M.F) > 0.05 & pSpecies > pthresh & pStage > pthresh, "SEX", 
                                                   ifelse(pSpecies < pthresh & abs(species_H.C) > 0.05 &  pSex > pthresh & pStage > pthresh,"SPECIES",
                                                          ifelse(pStage < pthresh & abs(stage_ADL.CHK) > 0.05  & pSex > pthresh & pSpecies > pthresh ,"STAGE",
                                                                 ifelse(pSex < pthresh & pSpecies < pthresh & abs(sex_M.F) > 0.05 & abs(species_H.C) > 0.05 &  pStage > pthresh ,"SEX_SPECIES", 
                                                                        ifelse(pSex < pthresh & pStage < pthresh & abs(sex_M.F) > 0.05 & abs(stage_ADL.CHK) > 0.05 & pSex > pthresh ,"SEX_STAGE",
                                                                               ifelse(pSpecies < pthresh & pStage < pthresh & abs(species_H.C) > 0.05 & abs(stage_ADL.CHK) > 0.05 &  pSex > pthresh ,"SPECIES_STAGE",
                                                                                      ifelse(pSpecies < pthresh & pStage < pthresh &  pSex < pthresh ,"ALL","UNKNOWN"))))))))
master_meth %>% dplyr::count(OD)
master_meth %>% dplyr::count(SING)
master_meth <- subset(master_meth, OD == 'OKAY')
master_meth %>% dplyr::count(DMR)

### Or with FDR
pval <- 0.05
master_meth$fSex <- p.adjust(master_meth$pSex,method='BH',n=nrow(master_meth))
master_meth$fStage <- p.adjust(master_meth$pStage,method='BH',n=nrow(master_meth))
master_meth$fSpecies <- p.adjust(master_meth$pSpecies,method='BH',n=nrow(master_meth))
master_meth <- master_meth %>% mutate(OD = ifelse(DispStat < disp_limits[1] | DispStat > disp_limits[2], "OVERDISPERSED","OKAY"),
                                      SING = ifelse(IsSingular == FALSE, "OKAY","SINGULAR"),
                                      DMR = ifelse(fSex < pval & abs(sex_M.F) > 0.05 & fSpecies > pval & fStage > pval, "SEX", 
                                                   ifelse(fSpecies < pval & abs(species_H.C) > 0.05 &  fSex > pval & fStage > pval,"SPECIES",
                                                          ifelse(fStage < pval & abs(stage_ADL.CHK) > 0.05  & fSex > pval & fSpecies > pval ,"STAGE",
                                                                 ifelse(fSex < pval & fSpecies < pval & abs(sex_M.F) > 0.05 & abs(species_H.C) > 0.05 &  fStage > pval ,"SEX_SPECIES", 
                                                                        ifelse(fSex < pval & fStage < pval & abs(sex_M.F) > 0.05 & abs(stage_ADL.CHK) > 0.05 & fSex > pval ,"SEX_STAGE",
                                                                               ifelse(fSpecies < pval & fStage < pval & abs(species_H.C) > 0.05 & abs(stage_ADL.CHK) > 0.05 &  fSex > pval ,"SPECIES_STAGE",
                                                                                      ifelse(fSpecies < pval & fStage < pval &  fSex < pval ,"ALL","UNKNOWN"))))))))


master_meth %>% dplyr::count(OD)
master_meth %>% dplyr::count(SING)
#master_meth <- subset(master_meth, OD == 'OKAY')
master_meth %>% dplyr::count(DMR)

#save for manhattan
write.csv(master_meth,file='CG_Manhattan_5x_bonferonni.csv')

#hunt for colors
#display.brewer.all(colorblindFriendly=TRUE) 
cols <- brewer.pal(8,'Dark2')
#cols <- inferno(5)
shapes <- c(15,16,17,18,16,16,1)

#plot
p1 <- ggplot(master_meth, aes(x=sex_M.F,y=threshold(-log10(pSex),max=quantile(-log10(pSex),0.99)))) +
  geom_point(aes(colour = DMR,pch= DMR),size=4,show.legend = FALSE) +
  scale_colour_manual(values = c("UNKNOWN"= "black","SEX"=cols[5],"SPECIES"=cols[3],"SEX_SPECIES"=cols[4], "STAGE"=cols[2],"SEX_STAGE"=cols[6],"SPECIES_STAGE"=cols[7],"ALL"=cols[8]))+
  theme_classic(base_size = 16)+xlab("5mC % Difference")+
  scale_shape_manual(values=shapes)+
  geom_hline(yintercept=-log10(pthresh),col="maroon",lwd=1.5,lty='dashed')+
  geom_vline(xintercept=c(-0.05,0.05),lwd=1.5,col="black")+
  ggtitle("Sex")+
  ylab('-log10(p)')
p1

p2 <- ggplot(master_meth, aes(x=stage_ADL.CHK,y=threshold(-log10(pStage),max=quantile(-log10(pStage),0.99)))) +
  geom_point(aes(colour = DMR,pch= DMR), size=4,show.legend = FALSE) +
  scale_colour_manual(values = c("UNKNOWN"= "black","SEX"=cols[5],"SPECIES"=cols[3],"SEX_SPECIES"=cols[4], "STAGE"=cols[2],"SEX_STAGE"=cols[6],"SPECIES_STAGE"=cols[7],"ALL"=cols[8]))+
  theme_classic(base_size = 16)+xlab("5mC % Difference")+
  scale_shape_manual(values=shapes)+
  geom_hline(yintercept=-log10(pthresh),col="maroon",lwd=1.5,lty='dashed')+
  geom_vline(xintercept=c(-0.05,0.05),lwd=1.5,col="black")+
  ggtitle("Stage")+
  ylab('-log10(p)')
p2

p3 <- ggplot(master_meth, aes(x=species_H.C,y=threshold(-log10(pSpecies),max=quantile(-log10(pSpecies),0.99)))) +
  geom_point(aes(colour = DMR,pch= DMR), size=4,show.legend = FALSE) +
  scale_colour_manual(values = c("UNKNOWN"= "black","SEX"=cols[5],"SPECIES"=cols[3],"SEX_SPECIES"=cols[4], "STAGE"=cols[2],"SEX_STAGE"=cols[6],"SPECIES_STAGE"=cols[7],"ALL"=cols[8]))+
  theme_classic(base_size = 16)+xlab("5mC % Difference")+
  geom_hline(yintercept=-log10(pthresh),col="maroon",lwd=1.5,lty='dashed')+
  geom_vline(xintercept=c(-0.05,0.05),lwd=1.5,col="black")+
  scale_shape_manual(values=shapes)+
  ggtitle("Species")+
  ylab('-log10(p)')
p3

#save it 
jpeg("plots/CG_DMA_Volcano_5x_bonferonni.jpg",units="in",res=300, height=5, width=12)
theme_set(theme_classic(base_size = 16))
grid.arrange(p1,p3,p2,nrow=1)
dev.off()


#save file
write.csv(master_meth,file='CG_DMA_DMRS_5x_bonferonni.csv',row.names=FALSE)

```

## HZ 5x

```
library(methylKit)
library(lme4)
library(Matrix)
library(emmeans)
library(multcomp)
library(matrixStats)
library(ggplot2)
library(gridExtra)
library(HDInterval)
library(dplyr)
library(stringr)

setwd("E:/Research/scratch/crow_hybrid_paper/CGI/HZ/dplyr/")

#import all the data and average by CGI

tab <- read.table('D_Ba_H02_BL_CHK_M.CGI.bed',header=FALSE); tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);D_Ba_H02 <- tab2[,c(7,1,2,3,8,9,10)]; names(D_Ba_H02) <- c('site','chr','start','end','numCs1','numTs1','Sites1');D_Ba_H02$D_Ba_H02 <- (D_Ba_H02$numCs1/(D_Ba_H02$numCs1+D_Ba_H02$numTs1))
tab <- read.table('D_Ba_H09_BL_CHK_M.CGI.bed',header=FALSE); tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);D_Ba_H09 <- tab2[,c(7,1,2,3,8,9,10)]; names(D_Ba_H09) <- c('site','chr','start','end','numCs2','numTs2','Sites2');D_Ba_H09$D_Ba_H09 <- (D_Ba_H09$numCs2/(D_Ba_H09$numCs2+D_Ba_H09$numTs2))
tab <- read.table('D_Ba_H19_BL_CHK_M.CGI.bed',header=FALSE); tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);D_Ba_H19 <- tab2[,c(7,1,2,3,8,9,10)]; names(D_Ba_H19) <- c('site','chr','start','end','numCs3','numTs3','Sites3');D_Ba_H19$D_Ba_H19 <- (D_Ba_H19$numCs3/(D_Ba_H19$numCs3+D_Ba_H19$numTs3))
tab <- read.table('D_Ba_H21_BL_CHK_M.CGI.bed',header=FALSE); tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);D_Ba_H21 <- tab2[,c(7,1,2,3,8,9,10)]; names(D_Ba_H21) <- c('site','chr','start','end','numCs4','numTs4','Sites4');D_Ba_H21$D_Ba_H21 <- (D_Ba_H21$numCs4/(D_Ba_H21$numCs4+D_Ba_H21$numTs4))
tab <- read.table('D_Hi_C03_BL_CHK_M.CGI.bed',header=FALSE); tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);D_Hi_C03 <- tab2[,c(7,1,2,3,8,9,10)]; names(D_Hi_C03) <- c('site','chr','start','end','numCs5','numTs5','Sites5');D_Hi_C03$D_Hi_C03 <- (D_Hi_C03$numCs5/(D_Hi_C03$numCs5+D_Hi_C03$numTs5))
tab <- read.table('D_Ko_C22_BL_CHK_M.CGI.bed',header=FALSE); tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);D_Ko_C22 <- tab2[,c(7,1,2,3,8,9,10)]; names(D_Ko_C22) <- c('site','chr','start','end','numCs6','numTs6','Sites6');D_Ko_C22$D_Ko_C22 <- (D_Ko_C22$numCs6/(D_Ko_C22$numCs6+D_Ko_C22$numTs6))
tab <- read.table('D_Ko_C40_BL_CHK_M.CGI.bed',header=FALSE); tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);D_Ko_C40 <- tab2[,c(7,1,2,3,8,9,10)]; names(D_Ko_C40) <- c('site','chr','start','end','numCs7','numTs7','Sites7');D_Ko_C40$D_Ko_C40 <- (D_Ko_C40$numCs7/(D_Ko_C40$numCs7+D_Ko_C40$numTs7))
tab <- read.table('D_Ko_C42_BL_CHK_M.CGI.bed',header=FALSE); tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);D_Ko_C42 <- tab2[,c(7,1,2,3,8,9,10)]; names(D_Ko_C42) <- c('site','chr','start','end','numCs8','numTs8','Sites8');D_Ko_C42$D_Ko_C42 <- (D_Ko_C42$numCs8/(D_Ko_C42$numCs8+D_Ko_C42$numTs8))
tab <- read.table('D_Ku_H02_BL_CHK_M.CGI.bed',header=FALSE); tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);D_Ku_H02 <- tab2[,c(7,1,2,3,8,9,10)]; names(D_Ku_H02) <- c('site','chr','start','end','numCs9','numTs9','Sites9');D_Ku_H02$D_Ku_H02 <- (D_Ku_H02$numCs9/(D_Ku_H02$numCs9+D_Ku_H02$numTs9))
tab <- read.table('D_Lo_C19_BL_CHK_M.CGI.bed',header=FALSE); tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);D_Lo_C19 <- tab2[,c(7,1,2,3,8,9,10)]; names(D_Lo_C19) <- c('site','chr','start','end','numCs10','numTs10','Sites10');D_Lo_C19$D_Lo_C19 <- (D_Lo_C19$numCs10/(D_Lo_C19$numCs10+D_Lo_C19$numTs10))
tab <- read.table('D_Lo_C20_BL_CHK_M.CGI.bed',header=FALSE); tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);D_Lo_C20 <- tab2[,c(7,1,2,3,8,9,10)]; names(D_Lo_C20) <- c('site','chr','start','end','numCs11','numTs11','Sites11');D_Lo_C20$D_Lo_C20 <- (D_Lo_C20$numCs11/(D_Lo_C20$numCs11+D_Lo_C20$numTs11))
tab <- read.table('D_Ne_Y03_BL_CHK_M.CGI.bed',header=FALSE); tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);D_Ne_Y03 <- tab2[,c(7,1,2,3,8,9,10)]; names(D_Ne_Y03) <- c('site','chr','start','end','numCs12','numTs12','Sites12');D_Ne_Y03$D_Ne_Y03 <- (D_Ne_Y03$numCs12/(D_Ne_Y03$numCs12+D_Ne_Y03$numTs12))
tab <- read.table('D_Ne_Y14_BL_CHK_M.CGI.bed',header=FALSE); tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);D_Ne_Y14 <- tab2[,c(7,1,2,3,8,9,10)]; names(D_Ne_Y14) <- c('site','chr','start','end','numCs13','numTs13','Sites13');D_Ne_Y14$D_Ne_Y14 <- (D_Ne_Y14$numCs13/(D_Ne_Y14$numCs13+D_Ne_Y14$numTs13))
tab <- read.table('D_Ne_Y31_BL_CHK_M.CGI.bed',header=FALSE); tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);D_Ne_Y31 <- tab2[,c(7,1,2,3,8,9,10)]; names(D_Ne_Y31) <- c('site','chr','start','end','numCs14','numTs14','Sites14');D_Ne_Y31$D_Ne_Y31 <- (D_Ne_Y31$numCs14/(D_Ne_Y31$numCs14+D_Ne_Y31$numTs14))
tab <- read.table('D_Ne_Y36_BL_CHK_M.CGI.bed',header=FALSE); tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);D_Ne_Y36 <- tab2[,c(7,1,2,3,8,9,10)]; names(D_Ne_Y36) <- c('site','chr','start','end','numCs15','numTs15','Sites15');D_Ne_Y36$D_Ne_Y36 <- (D_Ne_Y36$numCs15/(D_Ne_Y36$numCs15+D_Ne_Y36$numTs15))
tab <- read.table('D_Ne_Y40_BL_CHK_M.CGI.bed',header=FALSE); tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);D_Ne_Y40 <- tab2[,c(7,1,2,3,8,9,10)]; names(D_Ne_Y40) <- c('site','chr','start','end','numCs16','numTs16','Sites16');D_Ne_Y40$D_Ne_Y40 <- (D_Ne_Y40$numCs16/(D_Ne_Y40$numCs16+D_Ne_Y40$numTs16))
tab <- read.table('D_Ne_Y42_BL_CHK_M.CGI.bed',header=FALSE); tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);D_Ne_Y42 <- tab2[,c(7,1,2,3,8,9,10)]; names(D_Ne_Y42) <- c('site','chr','start','end','numCs17','numTs17','Sites17');D_Ne_Y42$D_Ne_Y42 <- (D_Ne_Y42$numCs17/(D_Ne_Y42$numCs17+D_Ne_Y42$numTs17))
tab <- read.table('D_Rb_Y07_BL_CHK_M.CGI.bed',header=FALSE); tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);D_Rb_Y07 <- tab2[,c(7,1,2,3,8,9,10)]; names(D_Rb_Y07) <- c('site','chr','start','end','numCs18','numTs18','Sites18');D_Rb_Y07$D_Rb_Y07 <- (D_Rb_Y07$numCs18/(D_Rb_Y07$numCs18+D_Rb_Y07$numTs18))
tab <- read.table('D_Rb_Y08_BL_CHK_M.CGI.bed',header=FALSE); tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);D_Rb_Y08 <- tab2[,c(7,1,2,3,8,9,10)]; names(D_Rb_Y08) <- c('site','chr','start','end','numCs19','numTs19','Sites19');D_Rb_Y08$D_Rb_Y08 <- (D_Rb_Y08$numCs19/(D_Rb_Y08$numCs19+D_Rb_Y08$numTs19))
tab <- read.table('D_Rb_Y13_BL_CHK_M.CGI.bed',header=FALSE); tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);D_Rb_Y13 <- tab2[,c(7,1,2,3,8,9,10)]; names(D_Rb_Y13) <- c('site','chr','start','end','numCs20','numTs20','Sites20');D_Rb_Y13$D_Rb_Y13 <- (D_Rb_Y13$numCs20/(D_Rb_Y13$numCs20+D_Rb_Y13$numTs20))
tab <- read.table('D_Rb_Y15_BL_CHK_M.CGI.bed',header=FALSE); tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);D_Rb_Y15 <- tab2[,c(7,1,2,3,8,9,10)]; names(D_Rb_Y15) <- c('site','chr','start','end','numCs21','numTs21','Sites21');D_Rb_Y15$D_Rb_Y15 <- (D_Rb_Y15$numCs21/(D_Rb_Y15$numCs21+D_Rb_Y15$numTs21))
tab <- read.table('S_Up_H77_BL_CHK_M.CGI.bed',header=FALSE); tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);S_Up_H77 <- tab2[,c(7,1,2,3,8,9,10)]; names(S_Up_H77) <- c('site','chr','start','end','numCs22','numTs22','Sites22');S_Up_H77$S_Up_H77 <- (S_Up_H77$numCs22/(S_Up_H77$numCs22+S_Up_H77$numTs22))
tab <- read.table('S_Up_H80_BL_CHK_M.CGI.bed',header=FALSE); tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);S_Up_H80 <- tab2[,c(7,1,2,3,8,9,10)]; names(S_Up_H80) <- c('site','chr','start','end','numCs23','numTs23','Sites23');S_Up_H80$S_Up_H80 <- (S_Up_H80$numCs23/(S_Up_H80$numCs23+S_Up_H80$numTs23))
tab <- read.table('S_Up_H81_BL_CHK_M.CGI.bed',header=FALSE); tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);S_Up_H81 <- tab2[,c(7,1,2,3,8,9,10)]; names(S_Up_H81) <- c('site','chr','start','end','numCs24','numTs24','Sites24');S_Up_H81$S_Up_H81 <- (S_Up_H81$numCs24/(S_Up_H81$numCs24+S_Up_H81$numTs24))

head(S_Up_H81)

master <- Reduce(function(x,y) merge(x = x, y = y, by=c('site','chr','start','end')),
                 list(D_Ba_H02,D_Ba_H09,D_Ba_H19,D_Ba_H21,D_Hi_C03,D_Ko_C22,D_Ko_C40,D_Ko_C42,D_Ku_H02,D_Lo_C19,D_Lo_C20,D_Ne_Y03,D_Ne_Y14,D_Ne_Y31,D_Ne_Y36,D_Ne_Y40,D_Ne_Y42,D_Rb_Y07,D_Rb_Y08,D_Rb_Y13,D_Rb_Y15,S_Up_H77,S_Up_H80,S_Up_H81))
head(master)

#calculate coverage and CpGs
cov <- master %>% mutate(minsite1 = ifelse(Sites1 < 3, NA, Sites1),
                         minsite2 = ifelse(Sites2 < 3, NA, Sites2),
                         minsite3 = ifelse(Sites3 < 3, NA, Sites3),
                         minsite4 = ifelse(Sites4 < 3, NA, Sites4),
                         minsite5 = ifelse(Sites5 < 3, NA, Sites5),
                         minsite6 = ifelse(Sites6 < 3, NA, Sites6),
                         minsite7 = ifelse(Sites7 < 3, NA, Sites7),
                         minsite8 = ifelse(Sites8 < 3, NA, Sites8),
                         minsite9 = ifelse(Sites9 < 3, NA, Sites9),
                         minsite10 = ifelse(Sites10 < 3, NA, Sites10),
                         minsite11 = ifelse(Sites11 < 3, NA, Sites11),
                         minsite12 = ifelse(Sites12 < 3, NA, Sites12),
                         minsite13 = ifelse(Sites13 < 3, NA, Sites13),
                         minsite14 = ifelse(Sites14 < 3, NA, Sites14),
                         minsite15 = ifelse(Sites15 < 3, NA, Sites15),
                         minsite16 = ifelse(Sites16 < 3, NA, Sites16),
                         minsite17 = ifelse(Sites17 < 3, NA, Sites17),
                         minsite18 = ifelse(Sites18 < 3, NA, Sites18),
                         minsite19 = ifelse(Sites19 < 3, NA, Sites19),
                         minsite20 = ifelse(Sites20 < 3, NA, Sites20),
                         minsite21 = ifelse(Sites21 < 3, NA, Sites21),
                         minsite22 = ifelse(Sites22 < 3, NA, Sites22),
                         minsite23 = ifelse(Sites23 < 3, NA, Sites23),
                         minsite24 = ifelse(Sites24 < 3, NA, Sites24))

#count NAs
colSums(is.na(cov))

#filter according to missingness
sitedat <- cov[rowSums(is.na(cov[grepl('^minsite', names(cov))])) <= 4, ]
nrow(cov)
nrow(sitedat)

master <- sitedat
nrow(master)

#count NAs
colSums(is.na(master))

#define groups
master$HOODED <- rowMeans(master[grepl('^S_Up', names(master))],na.rm=TRUE)
master$CARRION <- rowMeans(master[grepl('^D_Ko', names(master))],na.rm=TRUE)

compare <- master[,c(1,(ncol(master)-0-1):ncol(master))]
head(compare)
#now we must do 15 pairwise comparisons between our 3 treatments for % 5mC, keep track of this in the template file
compare[4] <- compare[2]-compare[3]

head(compare)
win_diff <- compare[,c(4)]

#merge the total percent diff with the original frame
master_meth <- cbind(master,win_diff)
colnames(master_meth)[colnames(master_meth) == 'win_diff'] <- 'species_H.C'
nrow(master_meth)
names(master_meth)

#create the GLMM input #C=0 Hybrid=1 H=2 / Years / etc
IND1 <- as.data.frame(c(master_meth[1],"D_Ba_H02",master_meth[5], master_meth[6],"0.574193","0","588","1","Ba_0"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Hybrid_Index","Year","Distance","Classification","Locality"))
IND2 <- as.data.frame(c(master_meth[1],"D_Ba_H09",master_meth[9], master_meth[10],"0.5466904","0","588","1","Ba_0"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Hybrid_Index","Year","Distance","Classification","Locality"))
IND3 <- as.data.frame(c(master_meth[1],"D_Ba_H19",master_meth[13], master_meth[14],"0.4055915","0","588","1","Ba_0"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Hybrid_Index","Year","Distance","Classification","Locality"))
IND4 <- as.data.frame(c(master_meth[1],"D_Ba_H21",master_meth[17], master_meth[18],"0.5258975","0","588","1","Ba_0"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Hybrid_Index","Year","Distance","Classification","Locality"))
IND5 <- as.data.frame(c(master_meth[1],"D_Hi_C03",master_meth[21], master_meth[22],"0.4279232","1","466","1","Hi_0"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Hybrid_Index","Year","Distance","Classification","Locality"))
IND6 <- as.data.frame(c(master_meth[1],"D_Ko_C22",master_meth[25], master_meth[26],"0","1","0","0","Ko_1"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Hybrid_Index","Year","Distance","Classification","Locality"))
IND7 <- as.data.frame(c(master_meth[1],"D_Ko_C40",master_meth[29], master_meth[30],"0","1","0","0","Ko_1"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Hybrid_Index","Year","Distance","Classification","Locality"))
IND8 <- as.data.frame(c(master_meth[1],"D_Ko_C42",master_meth[33], master_meth[34],"0","1","0","0","Ko_1"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Hybrid_Index","Year","Distance","Classification","Locality"))
IND9 <- as.data.frame(c(master_meth[1],"D_Ku_H02",master_meth[37], master_meth[38],"0.7063256","0","589","1","Ku_2"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Hybrid_Index","Year","Distance","Classification","Locality"))
IND10 <- as.data.frame(c(master_meth[1],"D_Lo_C19",master_meth[41], master_meth[42],"0.2260644","1","440","1","Lo_2"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Hybrid_Index","Year","Distance","Classification","Locality"))
IND11 <- as.data.frame(c(master_meth[1],"D_Lo_C20",master_meth[45], master_meth[46],"0.3307552","1","440","1","Lo_2"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Hybrid_Index","Year","Distance","Classification","Locality"))
IND12 <- as.data.frame(c(master_meth[1],"D_Ne_Y03",master_meth[49], master_meth[50],"0.4295311","2","546","1","Ne_2"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Hybrid_Index","Year","Distance","Classification","Locality"))
IND13 <- as.data.frame(c(master_meth[1],"D_Ne_Y14",master_meth[53], master_meth[54],"0.6457692","2","546","1","Ne_2"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Hybrid_Index","Year","Distance","Classification","Locality"))
IND14 <- as.data.frame(c(master_meth[1],"D_Ne_Y31",master_meth[57], master_meth[58],"0.7511844","2","546","1","Ne_2"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Hybrid_Index","Year","Distance","Classification","Locality"))
IND15 <- as.data.frame(c(master_meth[1],"D_Ne_Y36",master_meth[61], master_meth[62],"0.6105241","2","546","1","Ne_2"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Hybrid_Index","Year","Distance","Classification","Locality"))
IND16 <- as.data.frame(c(master_meth[1],"D_Ne_Y40",master_meth[65], master_meth[66],"0.6080164","2","546","1","Ne_2"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Hybrid_Index","Year","Distance","Classification","Locality"))
IND17 <- as.data.frame(c(master_meth[1],"D_Ne_Y42",master_meth[69], master_meth[70],"0.5530574","2","546","1","Ne_2"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Hybrid_Index","Year","Distance","Classification","Locality"))
IND18 <- as.data.frame(c(master_meth[1],"D_Rb_Y07",master_meth[73], master_meth[74],"0.6271691","2","513","1","Rb_2"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Hybrid_Index","Year","Distance","Classification","Locality"))
IND19 <- as.data.frame(c(master_meth[1],"D_Rb_Y08",master_meth[77], master_meth[78],"0.5266439","2","513","1","Rb_1"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Hybrid_Index","Year","Distance","Classification","Locality"))
IND20 <- as.data.frame(c(master_meth[1],"D_Rb_Y13",master_meth[81], master_meth[82],"0.7101439","2","513","1","Rb_1"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Hybrid_Index","Year","Distance","Classification","Locality"))
IND21 <- as.data.frame(c(master_meth[1],"D_Rb_Y15",master_meth[85], master_meth[86],"0.4707096","2","513","1","Rb_1"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Hybrid_Index","Year","Distance","Classification","Locality"))
IND22 <- as.data.frame(c(master_meth[1],"S_Up_H77",master_meth[89], master_meth[90],"1","1","1476","2","Up_1"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Hybrid_Index","Year","Distance","Classification","Locality"))
IND23 <- as.data.frame(c(master_meth[1],"S_Up_H80",master_meth[93], master_meth[94],"1","1","1476","2","Up_1"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Hybrid_Index","Year","Distance","Classification","Locality"))
IND24 <- as.data.frame(c(master_meth[1],"S_Up_H81",master_meth[97], master_meth[98],"1","1","1476","2","Up_1"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Hybrid_Index","Year","Distance","Classification","Locality"))

head(IND4)
head(IND24)

#Bind them all together and order it by window
datmeth <- rbind(IND1,IND2,IND3,IND4,IND5,IND6,IND7,IND8,IND9,IND10,IND11,IND12,IND13,IND14,IND15,IND16,IND17,IND18,IND19,IND20,IND21,IND22,IND23,IND24)
datmeth <- datmeth[order(datmeth$CGI),]

head(datmeth)

#sanity check 
nrow(datmeth)/24

#plot correlations between variables
dats <- read.table("metadata_hz.txt",header=TRUE)
text <- cor(dats$Hybrid_Index,dats$Year,method='pearson')
ggplot(dats,aes(x=Year,y=Hybrid_Index,col=Hybrid_Index))+geom_point()+theme_classic()+
  scale_color_gradient(low='black',high='grey90')+
  geom_smooth(method='lm')+ggtitle('Distance to Hybrid Index')+
  geom_text(aes(-Inf,Inf),label=paste("Pearson cor = ",round(text,2)),vjust=2,hjust=-.15,cex=4,col='black')


#save some files we might need later 
write.csv(datmeth,file="HZ_GLMM_Input_5x.csv",row.names=FALSE)
write.csv(master_meth,file="HZ_FULL_5x.csv",row.names=FALSE)

#now we can start here 
datmeth <- read.csv("HZ_GLMM_Input_5x.csv",header=TRUE)
head(datmeth)
nrow(datmeth)/24
#set these as factor variables 
datmeth$Hybrid_Index <- as.numeric(datmeth$Hybrid_Index)
datmeth$Year <- as.factor(datmeth$Year)
datmeth$Locality <- as.factor(datmeth$Locality)
datmeth$ID <- as.factor(datmeth$ID)
datmeth$Classification <- as.factor(datmeth$Classification)

#and scale the continuous variables
datmeth <- datmeth %>% mutate_at(c('Hybrid_Index', 'Distance'), ~(scale(.) %>% as.vector))

##### MODEL #####
DMA_GLMM_pValues <- NULL
DMA_GLMM_Estimates <- NULL
DMA_GLMM_RandomEffect_Var <- NULL
DMA_GLMM_SEs <- NULL
DMA_GLMM_DispersionStat <- NULL
DMA_GLMM_Singular <- NULL

for (i in seq(1, nrow(datmeth), by = 24)) {
  tryCatch({
    cat('\n Running model for CpG site: ',((i-1)/24))
    # model each site
    CpGsite <- datmeth[i:(i + 23), ]
    # run the model examining relationship between hybrid index and locality 
    glmm_CpG <- glmer(cbind(MethCounts, UnMethCounts) ~
                        Hybrid_Index + (1|Year) ,data = CpGsite, family = "binomial",
                      control = glmerControl(optimizer = "bobyqa", boundary.tol = 0.01,
                                             optCtrl = list(maxfun = 2e+08)))
    DMA_GLMM_Singularf <- isSingular(glmm_CpG)
    sum_mod <- summary(glmm_CpG)
    
    ## For random effect 
    GLMM_RandomEffect_Var <- as.data.frame(VarCorr(glmm_CpG))$sdcor
    
    # calculate dispersion statistic (following AF Zuur, JM Hilbe
    # & EN Ieno. A beginner’s guide to GLM and GLMM with R - a
    # frequentist and Bayesian perspective for
    # ecologists(Highland Statistics Ltd., 2013).)
    E1 <- residuals(glmm_CpG)
    # number of parameters: fixed effects + 1 random effect
    p1 <- length(fixef(glmm_CpG)) + 1
    GLMM_DispersionStat <- sum(E1^2)/(nrow(CpGsite) - p1)
    # format estimates to add estimates for each loop to
    # estimates from previous loops
    
    GLMM_pValues_f <- data.frame(site = as.character(CpGsite[1,1]),Hybrid_Index = sum_mod$coefficients[2,4])
    GLMM_SEs_f <- data.frame(site = as.character(CpGsite[1,1]),Hybrid_Index = sum_mod$coefficients[2,2])
    GLMM_Estimates_f <- data.frame(site = as.character(CpGsite[1,1]),Hybrid_Index = sum_mod$coefficients[2,1])
    GLMM_DispersionStat_f <- data.frame(site = as.character(CpGsite[1,1]), DispStat = GLMM_DispersionStat)
    GLMM_RandomEffect_Var_f <- data.frame(site = as.character(CpGsite[1,1]), ID = GLMM_RandomEffect_Var[1])
    # add estimates for each loop to estimates from previous
    # loops
    DMA_GLMM_pValues <- rbind(DMA_GLMM_pValues, GLMM_pValues_f)
    DMA_GLMM_Estimates <- rbind(DMA_GLMM_Estimates, GLMM_Estimates_f)
    DMA_GLMM_SEs <- rbind(DMA_GLMM_SEs, GLMM_SEs_f)
    DMA_GLMM_RandomEffect_Var <- rbind(DMA_GLMM_RandomEffect_Var,GLMM_RandomEffect_Var_f)
    DMA_GLMM_DispersionStat <- rbind(DMA_GLMM_DispersionStat,GLMM_DispersionStat_f)
    DMA_GLMM_Singular <- rbind(DMA_GLMM_Singular,DMA_GLMM_Singularf)
    
    #Save errors to this output
  }, error=function(e){cat(unique(CpGsite[,1]), "\n")
  }
  )
}

# 3. Save loop results
save(DMA_GLMM_pValues, DMA_GLMM_Estimates, DMA_GLMM_SEs, 
     DMA_GLMM_DispersionStat, DMA_GLMM_RandomEffect_Var,DMA_GLMM_Singular, file = "HZ_DMA_CrowResults_5x.RData")

#should be the samem otherwise we have models that didn't converge 
nrow(DMA_GLMM_pValues)
nrow(datmeth)/24

```

### Volcanos

```##
library(methylKit)
library(lme4)
library(Matrix)
library(emmeans)
library(multcomp)
library(tidyr)
library(ggplot2)
library(psych)
library(gridExtra)
library(HDInterval)
library(dplyr)
library(data.table)
library(stringr)
library(LICORS)
library(RColorBrewer)
library(viridis)

setwd("E:/Research/scratch/crow_hybrid_paper/CGI/HZ/dplyr/")

load("HZ_DMA_CrowResults_5x.RData")

#check out some of the data
head(DMA_GLMM_Estimates)
dim(DMA_GLMM_Estimates)
head(DMA_GLMM_pValues)
head(DMA_GLMM_SEs)


#let's check out the distributio of 5mC contrasts
multi.hist(DMA_GLMM_Estimates[,2:2],nrow=NULL,ncol=NULL,breaks=21)

#plot dispersion boxplot, we will exclude any sites outside the 99% interval
disp_limits <- hdi(DMA_GLMM_DispersionStat$DispStat,credMass = 0.99)
disp_limits

#lower       upper 
#0.002117176 0.459554568

jpeg("plots/HZ_dispersion_plot.jpg",units="in",res=300, height=6, width=5)
theme_set(theme_classic(base_size = 16))
dplot <- ggplot(DMA_GLMM_DispersionStat,aes(y=DispStat,x=0))+geom_violin()+theme_classic(base_size = 16)+
  geom_jitter(width=0.1)+
  geom_text(aes(y = (disp_limits[2]+0.05),x=0.2, label = "Upper Limit"))+
  geom_hline(yintercept = disp_limits[2],linetype="dashed",col="red") +
  geom_text(aes(y = (disp_limits[1]+0.05),x=0.2, label = "Lower Limit"))+
  geom_hline(yintercept = disp_limits[1],linetype="dashed",col="red") +
  xlab("")+ylab("Dispersion Statistic")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
grid.arrange(dplot)
dev.off()

#subset bad sites
bad_sites <- subset(DMA_GLMM_DispersionStat,DispStat < disp_limits[1] | DispStat > disp_limits[2])
nrow(bad_sites)

# calculate 10% differences, now we care about direction so we re-do this, since before it was absolute value 
meth <- read.csv("HZ_FULL_5x.csv",header=TRUE)
dim(meth)
head(meth)

names(DMA_GLMM_pValues) <- c('site','pHI')
names(DMA_GLMM_Estimates) <- c('site','eHI')
names(DMA_GLMM_SEs) <- c('site','sHI')
DMA_RF_S <- cbind(DMA_GLMM_RandomEffect_Var,DMA_GLMM_Singular)
names(DMA_RF_S) <- c('site','RandomEffect','IsSingular')
master_meth1 <- merge(meth, DMA_GLMM_pValues, by='site')
master_meth2 <- merge(master_meth1, DMA_GLMM_Estimates, by='site')
master_meth3 <- merge(master_meth2, DMA_GLMM_SEs, by='site')
master_meth4 <- merge(master_meth3, DMA_GLMM_DispersionStat, by='site')
master_meth <- merge(master_meth4, DMA_RF_S, by='site')
head(master_meth)
nrow(master_meth)

#### Volcano with bonferonni
pthresh <- (0.05/nrow(master_meth))
master_meth <- master_meth %>% mutate(OD = ifelse(DispStat < disp_limits[1] | DispStat > disp_limits[2], "OVERDISPERSED","OKAY"),
                                      SING = ifelse(IsSingular == FALSE, "OKAY","SINGULAR"),
                                      DMR = ifelse((pHI < pthresh & abs(species_H.C) > 0.05) , "SPECIES","UNKNOWN"))

master_meth <- subset(master_meth, OD == 'OKAY')
master_meth %>% dplyr::count(DMR)

### With FDR, first select sites that are significant for each variable, but non-significant for the other variables 
pval <- 0.05
master_meth$fHI <- p.adjust(master_meth$pHI,method='BH',n=nrow(master_meth))

master_meth <- master_meth %>% mutate(OD = ifelse(DispStat < disp_limits[1] | DispStat > disp_limits[2], "OVERDISPERSED","OKAY"),
                                      SING = ifelse(IsSingular == FALSE, "OKAY","SINGULAR"),
                                      DMR = ifelse((fHI < pval & abs(species_H.C) > 0.05) , "SPECIES","UNKNOWN"))

#remove overdispersed sites 
master_meth %>% dplyr::count(OD)
master_meth %>% dplyr::count(SING)
master_meth <- subset(master_meth, OD == 'OKAY')
master_meth %>% dplyr::count(DMR)


#save for manhattan
write.csv(master_meth,file='HZ_Manhattan_5x_bonferonni.csv')

#hunt for colors
display.brewer.all(colorblindFriendly=TRUE) 
cols <- brewer.pal(5,'Paired')
#cols <- viridis(7)
shapes <- c(16,1)

#plot SPECIES
p1 <- ggplot(master_meth, aes(x=species_H.C, y=threshold(-log10(pHI),max=quantile(-log10(pHI),0.99))))+
  geom_point(aes(colour = DMR,shape=DMR), size=4,show.legend = FALSE) +
  scale_colour_manual(values = c("UNKNOWN"= "black","SPECIES"=cols[3]))+
  theme_classic(base_size = 16)+xlab("5mC % Difference")+
  ggtitle("Species")+
  scale_shape_manual(values=shapes)+
  geom_hline(yintercept=-log10(0.05),lwd=1,lty=2,col="maroon")+
  geom_vline(xintercept=c(-0.05,0.05),lwd=1.5,col="black")+
  ylab('-log10(P)')
p1

#save it 
jpeg("plots/HZ_DMA_Volcano_5x_bonferonni.jpg",units="in",res=300, height=5, width=7)
theme_set(theme_classic(base_size = 16))
grid.arrange(p1,nrow=1)
dev.off()

#save file
write.csv(master_meth,file='HZ_DMA_DMRS_5x_bonferonni.csv',row.names=FALSE)

```

## WGBS 5x

```
library(methylKit)
library(lme4)
library(Matrix)
library(emmeans)
library(multcomp)
library(matrixStats)
library(ggplot2)
library(gridExtra)
library(HDInterval)
library(dplyr)

setwd("E:/Research/scratch/crow_hybrid_paper/CGI/WGBS/dplyr/")

#import all the data and average by CGI
tab <- read.table('D_Ko_C31_BL_ADL_M.CGI.bed',header=FALSE);  tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);C31_BL_M <- tab2[,c(7,1,2,3,8,9,10)]; names(C31_BL_M) <- c('site','chr','start','end','numCs1','numTs1','Sites1');C31_BL_M$C31_BL_M <- (C31_BL_M$numCs1/(C31_BL_M$numCs1+C31_BL_M$numTs1))
tab <- read.table('D_Ko_C45_BL_ADL_F.CGI.bed',header=FALSE);  tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);C45_BL_F <- tab2[,c(7,1,2,3,8,9,10)]; names(C45_BL_F) <- c('site','chr','start','end','numCs2','numTs2','Sites2');C45_BL_F$C45_BL_F <- (C45_BL_F$numCs2/(C45_BL_F$numCs2+C45_BL_F$numTs2))
tab <- read.table('S_Up_H59_BL_ADL_M.CGI.bed',header=FALSE);  tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);H59_BL_M <- tab2[,c(7,1,2,3,8,9,10)]; names(H59_BL_M) <- c('site','chr','start','end','numCs3','numTs3','Sites3');H59_BL_M$H59_BL_M <- (H59_BL_M$numCs3/(H59_BL_M$numCs3+H59_BL_M$numTs3))
tab <- read.table('S_Up_H59_M_ADL_M.CGI.bed',header=FALSE);  tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);H59_M_M <- tab2[,c(7,1,2,3,8,9,10)]; names(H59_M_M) <- c('site','chr','start','end','numCs4','numTs4','Sites4');H59_M_M$H59_M_M <- (H59_M_M$numCs4/(H59_M_M$numCs4+H59_M_M$numTs4))
tab <- read.table('S_Up_H60_BL_ADL_F.CGI.bed',header=FALSE);  tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);H60_BL_F <- tab2[,c(7,1,2,3,8,9,10)]; names(H60_BL_F) <- c('site','chr','start','end','numCs5','numTs5','Sites5');H60_BL_F$H60_BL_F <- (H60_BL_F$numCs5/(H60_BL_F$numCs5+H60_BL_F$numTs5))
tab <- read.table('S_Up_H60_M_ADL_F.CGI.bed',header=FALSE);  tab <- subset(tab, (V5 + V6) > 5); tab1 <- tab %>% group_by(V7) %>% mutate(numCs = sum(V5), numTs = sum(V6), Sites = n()); tab2 <- as.data.frame(tab1[!duplicated(tab1[ , c('V7')]),]);H60_M_F <- tab2[,c(7,1,2,3,8,9,10)]; names(H60_M_F) <- c('site','chr','start','end','numCs6','numTs6','Sites6');H60_M_F$H60_M_F <- (H60_M_F$numCs6/(H60_M_F$numCs6+H60_M_F$numTs6))

head(C31_BL_M)

master <- Reduce(function(x,y) merge(x = x, y = y, by=c('site','chr','start','end')),
                 list(C31_BL_M,C45_BL_F,H59_BL_M,H59_M_M,H60_BL_F,H60_M_F))
head(master)

#calculate coverage and CpGs
cov <- master %>% mutate(minsite1 = ifelse(Sites1 < 3, NA, Sites1),
                         minsite2 = ifelse(Sites2 < 3, NA, Sites2),
                         minsite3 = ifelse(Sites3 < 3, NA, Sites3),
                         minsite4 = ifelse(Sites4 < 3, NA, Sites4),
                         minsite5 = ifelse(Sites5 < 3, NA, Sites5),
                         minsite6 = ifelse(Sites6 < 3, NA, Sites6))



#count NAs
colSums(is.na(cov))

#filter according to missingness
sitedat <- cov[rowSums(is.na(cov[grepl('^minsite', names(cov))])) == 0, ]
nrow(sitedat)

master <- sitedat

#define groups
master$MALE <- rowMeans(master[grepl('_M$', names(master))],na.rm=TRUE) #M=1
master$FEMALE <- rowMeans(master[grepl('_F$', names(master))],na.rm=TRUE) #F=0
master$BLOOD <- rowMeans(master[grepl('_BL_', names(master))],na.rm=TRUE) #BL=0
master$SPLEEN <- rowMeans(master[grepl('_M_', names(master))],na.rm=TRUE) #M=1
master$HOODED <- rowMeans(master[grepl('^H', names(master))],na.rm=TRUE) #H=0
master$CARRION <- rowMeans(master[grepl('^C', names(master))],na.rm=TRUE) #C=1
compare <- master[,c(1,(ncol(master)-4-1):ncol(master))]
head(compare)
#now we must do pairwise comparisons between our 3 treatments for % 5mC
compare[8] <- compare[2]-compare[3]
compare[9] <- compare[4]-compare[5]
compare[10] <- compare[6]-compare[7]
head(compare)
win_diff <- compare[,c(8,9,10)]

#merge the total percent diff with the original frame
master_meth <- cbind(master,win_diff)
colnames(master_meth)[colnames(master_meth) == 'MALE.1'] <- 'sex_M.F'
colnames(master_meth)[colnames(master_meth) == 'BLOOD.1'] <- 'tissue_BL.M'
colnames(master_meth)[colnames(master_meth) == 'HOODED.1'] <- 'species_H.C'
nrow(master_meth)
names(master_meth)

#now subset for GLMM # C=1 H=0 / M=1 F=0 / BL=0 M=1
IND1 <- as.data.frame(c(master_meth[1],"C31",master_meth[5], master_meth[6],"1","1","0"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Species","Sex","Tissue"))
IND2 <- as.data.frame(c(master_meth[1],"C45",master_meth[9], master_meth[10],"1","0","0"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Species","Sex","Tissue"))
IND3 <- as.data.frame(c(master_meth[1],"H59",master_meth[13], master_meth[14],"0","1","0"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Species","Sex","Tissue"))
IND4 <- as.data.frame(c(master_meth[1],"H59",master_meth[17], master_meth[18],"0","1","1"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Species","Sex","Tissue"))
IND5 <- as.data.frame(c(master_meth[1],"H60",master_meth[21], master_meth[22],"0","0","0"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Species","Sex","Tissue"))
IND6 <- as.data.frame(c(master_meth[1],"H60",master_meth[25], master_meth[26],"0","0","1"),col.names=c("CGI","ID","MethCounts","UnMethCounts","Species","Sex","Tissue"))

head(IND4)
head(IND1)

#Bind them all together and order it by window
datmeth <- rbind(IND1,IND2,IND3,IND4,IND5,IND6)
datmeth <- datmeth[order(datmeth$CGI),]

head(datmeth)

#sanity check 
nrow(datmeth)/6

#save some files we might need later 
write.csv(datmeth,file="WGBS_GLMM_Input_5x.csv",row.names=FALSE)
write.csv(master_meth,file="WGBS_FULL_5x.csv",row.names=FALSE)

#now we can start here 
datmeth <- read.csv("WGBS_GLMM_Input_5x.csv",header=TRUE)
head(datmeth)
#set these as factor variables 
datmeth$Species <- as.factor(datmeth$Species)
datmeth$Tissue <- as.factor(datmeth$Tissue)
datmeth$Sex <- as.factor(datmeth$Sex)
datmeth$ID <- as.factor(datmeth$ID)

##### MODEL #####
DMA_GLMM_pValues <- NULL
DMA_GLMM_Estimates <- NULL
DMA_GLMM_SEs <- NULL
DMA_GLMM_RandomEffect_Var <- NULL
DMA_GLMM_DispersionStat <- NULL
DMA_GLMM_Singular <- NULL

for (i in seq(1, nrow(datmeth), by = 6)) {
  tryCatch({
    cat('\n Running model for CpG site: ',((i-1)/6))
    # model each site
    CpGsite <- datmeth[i:(i + 5), ]
    # run the model and get the contrasts between sex, stage, species
    glmm_CpG <- glmer(cbind(MethCounts, UnMethCounts) ~
                        Sex + Tissue + Species + (1|ID) ,data = CpGsite, family = "binomial",
                      control = glmerControl(optimizer = "bobyqa", boundary.tol = 0.01,
                                             optCtrl = list(maxfun = 2e+09)))
    DMA_GLMM_Singularf <- isSingular(glmm_CpG)
    sum_mod <- summary(glmm_CpG)
    
    ## 1 random effects ( (1 | ID) )
    GLMM_RandomEffect_Var <- as.data.frame(VarCorr(glmm_CpG))$sdcor
    
    # calculate dispersion statistic (following AF Zuur, JM Hilbe
    # & EN Ieno. A beginner’s guide to GLM and GLMM with R - a
    # frequentist and Bayesian perspective for
    # ecologists(Highland Statistics Ltd., 2013).)
    E1 <- residuals(glmm_CpG)
    # number of parameters: fixed effects + 1 random effect
    p1 <- length(fixef(glmm_CpG)) + 1
    GLMM_DispersionStat <- sum(E1^2)/(nrow(CpGsite) - p1)
    # format estimates to add estimates for each loop to
    # estimates from previous loops
    
    # TO GET THE COLUMN HEADERS WITH THE CONTRASTS, check the excel file. First I look at the glmm_CpG_contrast object, which lists the models in order of 1-8, for this comparison. From here, make the 28 pairwise comparison charts, and substitute the 0 and 1 with your variables
    GLMM_Estimates_f <- data.frame(site = as.character(CpGsite[1,1]),Sex1 = sum_mod$coefficients[2,1], Tissue1 = sum_mod$coefficients[3,1], Species1 = sum_mod$coefficients[4,1])
    GLMM_SEs_f <- data.frame(site = as.character(CpGsite[1,1]),Sex1 = sum_mod$coefficients[2,2], Tissue1 = sum_mod$coefficients[3,2], Species1 = sum_mod$coefficients[4,2])
    GLMM_pValues_f <- data.frame(site = as.character(CpGsite[1,1]),Sex1 = sum_mod$coefficients[2,4], Tissue1 = sum_mod$coefficients[3,4], Species1 = sum_mod$coefficients[4,4])
    GLMM_RandomEffect_Var_f <- data.frame(site = as.character(CpGsite[1,
                                                                      1]), ID = GLMM_RandomEffect_Var[1])
    GLMM_DispersionStat_f <- data.frame(site = as.character(CpGsite[1,
                                                                    1]), DispStat = GLMM_DispersionStat)
    # add estimates for each loop to estimates from previous
    # loops
    DMA_GLMM_pValues <- rbind(DMA_GLMM_pValues, GLMM_pValues_f)
    DMA_GLMM_Estimates <- rbind(DMA_GLMM_Estimates, GLMM_Estimates_f)
    DMA_GLMM_SEs <- rbind(DMA_GLMM_SEs, GLMM_SEs_f)
    DMA_GLMM_RandomEffect_Var <- rbind(DMA_GLMM_RandomEffect_Var,GLMM_RandomEffect_Var_f)
    DMA_GLMM_DispersionStat <- rbind(DMA_GLMM_DispersionStat,GLMM_DispersionStat_f)
    DMA_GLMM_Singular <- rbind(DMA_GLMM_Singular,DMA_GLMM_Singularf)
    
    #Save errors to this output
  }, error=function(e){cat(unique(CpGsite[,1]), "\n")
  }
  )
}

# 3. Save loop results
save(DMA_GLMM_pValues, DMA_GLMM_Estimates, DMA_GLMM_SEs, DMA_GLMM_RandomEffect_Var,
     DMA_GLMM_DispersionStat, DMA_GLMM_Singular, file = "WGBS_DMA_CrowResults_5x.RData")

hist(DMA_GLMM_pValues$Species1)
table(DMA_GLMM_Singular)['FALSE']

#should be the samem otherwise we have models that didn't converge 
nrow(DMA_GLMM_pValues)
nrow(datmeth)/6 

```

### Volcanos

```
library(methylKit)
library(lme4)
library(Matrix)
library(emmeans)
library(multcomp)
library(tidyr)
library(ggplot2)
library(psych)
library(gridExtra)
library(HDInterval)
library(dplyr)
library(data.table)
library(stringr)
library(matrixStats)
library(LICORS)
library(viridis)
library(RColorBrewer)

setwd("E:/Research/scratch/crow_hybrid_paper/CGI/WGBS/dplyr")

load("WGBS_DMA_CrowResults_5x.RData")

#check out some of the data
head(DMA_GLMM_Estimates)
dim(DMA_GLMM_Estimates)
head(DMA_GLMM_pValues)
head(DMA_GLMM_SEs)


#let's check out the distributio of 5mC contrasts
multi.hist(DMA_GLMM_Estimates[,2:4],nrow=1,ncol=NULL,breaks=21)

#plot dispersion boxplot, we will exclude any sites outside the 99% interval
disp_limits <- hdi(DMA_GLMM_DispersionStat$DispStat,credMass = 0.99)
disp_limits

#lower        upper 
#4.228215e-02 3.381379e+03 

jpeg("plots/WGBS_dispersion_plot.jpg",units="in",res=300, height=6, width=5)
theme_set(theme_classic(base_size = 16))
dplot <- ggplot(DMA_GLMM_DispersionStat,aes(y=DispStat,x=0))+geom_violin()+theme_classic(base_size = 16)+
  geom_jitter(width=0.1)+
  geom_text(aes(y = (disp_limits[2]+0.05),x=0.2, label = "Upper Limit"))+
  geom_hline(yintercept = disp_limits[2],linetype="dashed",col="red") +
  geom_text(aes(y = (disp_limits[1]+0.05),x=0.2, label = "Lower Limit"))+
  geom_hline(yintercept = disp_limits[1],linetype="dashed",col="red") +
  xlab("")+ylab("Dispersion Statistic")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
grid.arrange(dplot)
dev.off()

#subset bad sites
bad_sites <- subset(DMA_GLMM_DispersionStat,DispStat < disp_limits[1] | DispStat > disp_limits[2])
nrow(bad_sites)

# calculate 10% differences, now we care about direction so we re-do this, since before it was absolute value 
meth <- read.csv("WGBS_FULL_5x.csv",header=TRUE)
head(meth)

names(DMA_GLMM_pValues) <- c('site','pSex','pTissue','pSpecies')
names(DMA_GLMM_Estimates) <- c('site','eSex','eTissue','eSpecies')
names(DMA_GLMM_SEs) <- c('site','sSex','sTissue','sSpecies')
DMA_RF_S <- cbind(DMA_GLMM_RandomEffect_Var,DMA_GLMM_Singular)
names(DMA_RF_S) <- c('site','RandomEffect','IsSingular')
master_meth1 <- merge(meth, DMA_GLMM_pValues, by='site')
master_meth2 <- merge(master_meth1, DMA_GLMM_Estimates, by='site')
master_meth3 <- merge(master_meth2, DMA_GLMM_SEs, by='site')
master_meth4 <- merge(master_meth3, DMA_RF_S, by='site')
master_meth <- merge(master_meth4, DMA_GLMM_DispersionStat, by='site')
head(master_meth)
nrow(master_meth)

#### Volcano with bonferonni
pthresh <- (0.05/nrow(master_meth))
master_meth <- master_meth %>% mutate(OD = ifelse(DispStat < disp_limits[1] | DispStat > disp_limits[2], "OVERDISPERSED","OKAY"),
                                      SING = ifelse(IsSingular == FALSE, "OKAY","SINGULAR"),
                                      DMR = ifelse(pSex < pthresh & abs(sex_M.F) > 0.05 & pSpecies > pthresh & pTissue > pthresh, "SEX", 
                                                   ifelse(pSpecies < pthresh & abs(species_H.C) > 0.05 & pSex > pthresh & pTissue > pthresh,"SPECIES",
                                                          ifelse(pTissue < pthresh & abs(tissue_BL.M) > 0.05 & pSex > pthresh & pSpecies > pthresh ,"TISSUE",
                                                                 ifelse(pSex < pthresh & pSpecies < pthresh & abs(sex_M.F) > 0.05 & abs(species_H.C) > 0.05 & pTissue > pthresh ,"SEX_SPECIES", 
                                                                        ifelse(pSex < pthresh & pTissue < pthresh & abs(sex_M.F) > 0.05 & abs(tissue_BL.M) > 0.05 &  pSex > pthresh ,"SEX_TISSUE",
                                                                               ifelse(pSpecies < pthresh & pTissue < pthresh & abs(tissue_BL.M) > 0.05 & abs(species_H.C) > 0.05 &  pSex > pthresh ,"SPECIES_TISSUE",
                                                                                      ifelse(pSpecies < pthresh & pTissue < pthresh &  pSex < pthresh &
                                                                                               abs(tissue_BL.M) > 0.05 & abs(species_H.C) > 0.05 & abs(sex_M.F) > 0.05,"ALL","UNKNOWN"))))))))


master_meth <- subset(master_meth, OD == 'OKAY')
master_meth %>% dplyr::count(DMR)

### Or with FDR
pval <- 0.05
master_meth$fSex <- p.adjust(master_meth$pSex,method='BH',n=nrow(master_meth))
master_meth$fTissue <- p.adjust(master_meth$pTissue,method='BH',n=nrow(master_meth))
master_meth$fSpecies <- p.adjust(master_meth$pSpecies,method='BH',n=nrow(master_meth))
master_meth <- master_meth %>% mutate(OD = ifelse(DispStat < disp_limits[1] | DispStat > disp_limits[2], "OVERDISPERSED","OKAY"),
                                      SING = ifelse(IsSingular == FALSE, "OKAY","SINGULAR"),
                                      DMR = ifelse(fSex < pval & abs(sex_M.F) > 0.05 & fSpecies > pval & fTissue > pval, "SEX", 
                                                   ifelse(fSpecies < pval & abs(species_H.C) > 0.05 & fSex > pval & fTissue > pval,"SPECIES",
                                                          ifelse(fTissue < pval & abs(tissue_BL.M) > 0.05 & fSex > pval & fSpecies > pval ,"TISSUE",
                                                                 ifelse(fSex < pval & fSpecies < pval & abs(sex_M.F) > 0.05 & abs(species_H.C) > 0.05 & fTissue > pval ,"SEX_SPECIES", 
                                                                        ifelse(fSex < pval & fTissue < pval & abs(sex_M.F) > 0.05 & abs(tissue_BL.M) > 0.05 &  fSex > pval ,"SEX_TISSUE",
                                                                               ifelse(fSpecies < pval & fTissue < pval & abs(tissue_BL.M) > 0.05 & abs(species_H.C) > 0.05 &  fSex > pval ,"SPECIES_TISSUE",
                                                                                      ifelse(fSpecies < pval & fTissue < pval &  fSex < pval &
                                                                                               abs(tissue_BL.M) > 0.05 & abs(species_H.C) > 0.05 & abs(sex_M.F) > 0.05,"ALL","UNKNOWN"))))))))

master_meth %>% dplyr::count(OD)
master_meth %>% dplyr::count(SING)
master_meth <- subset(master_meth, OD == 'OKAY')
master_meth %>% dplyr::count(DMR)

#save for manhattan
write.csv(master_meth,file='WGBS_Manhattan_5x_bonferonni.csv')

#hunt for colors
#display.brewer.all(colorblindFriendly=TRUE) 
cols <- brewer.pal(8,'Set2')
#cols <- inferno(5)
shapes <- c(15,16,17,18,16,16,1)

#plot sex
p1 <- ggplot(master_meth, aes(x=sex_M.F,y=threshold(-log10(pSex),max=quantile(-log10(pSex),0.99)))) +
  geom_point(aes(colour = DMR,pch= DMR),size=4,show.legend = FALSE) +
  scale_colour_manual(values = c("UNKNOWN"= "black","SEX"=cols[5],"SPECIES"=cols[3],"SEX_SPECIES"=cols[4], "TISSUE"=cols[2],"SEX_TISSUE"=cols[6],"SPECIES_TISSUE"=cols[7],"ALL"=cols[8]))+  theme_classic(base_size = 16)+xlab("5mC % Difference")+
  scale_shape_manual(values=shapes)+
  geom_hline(yintercept=-log10(pthresh),col="maroon",lwd=1.5,lty='dashed')+
  geom_vline(xintercept=c(-0.05,0.05),lwd=1.5,col="black")+
  ggtitle("Sex")+
  ylab('-log10(p)')
p1

#plot species
p2 <- ggplot(master_meth, aes(x=species_H.C,y=threshold(-log10(pSpecies),max=quantile(-log10(pSpecies),0.99)))) +
  geom_point(aes(colour = DMR,pch= DMR),size=4,show.legend = FALSE) +
  scale_colour_manual(values = c("UNKNOWN"= "black","SEX"=cols[5],"SPECIES"=cols[3],"SEX_SPECIES"=cols[4], "TISSUE"=cols[2],"SEX_TISSUE"=cols[6],"SPECIES_TISSUE"=cols[7],"ALL"=cols[8]))+  theme_classic(base_size = 16)+xlab("5mC % Difference")+
  theme_classic(base_size = 16)+xlab("5mC % Difference")+
  scale_shape_manual(values=shapes)+
  geom_hline(yintercept=-log10(pthresh),col="maroon",lwd=1.5,lty='dashed')+
  geom_vline(xintercept=c(-0.05,0.05),lwd=1.5,col="black")+
  ggtitle("Species")+
  ylab('-log10(p)')
p2

#plot tissue
p3 <- ggplot(master_meth, aes(x=tissue_BL.M,y=threshold(-log10(pTissue),max=quantile(-log10(pTissue),0.99)))) +
  geom_point(aes(colour = DMR,pch= DMR),size=4,show.legend = FALSE) +
  scale_colour_manual(values = c("UNKNOWN"= "black","SEX"=cols[5],"SPECIES"=cols[3],"SEX_SPECIES"=cols[4], "TISSUE"=cols[2],"SEX_TISSUE"=cols[6],"SPECIES_TISSUE"=cols[7],"ALL"=cols[8]))+  theme_classic(base_size = 16)+xlab("5mC % Difference")+
  theme_classic(base_size = 16)+xlab("5mC % Difference")+
  scale_shape_manual(values=shapes)+
  geom_hline(yintercept=-log10(pthresh),col="maroon",lwd=1.5,lty='dashed')+
  geom_vline(xintercept=c(-0.05,0.05),lwd=1.5,col="black")+
  ggtitle("Tissue")+
  ylab('-log10(p)')
p3


#save it 
jpeg("plots/WGBS_DMA_Volcano_5x_bonferonni.jpg",units="in",res=300, height=5, width=12)
theme_set(theme_classic(base_size = 16))
grid.arrange(p1,p2,p3,nrow=1)
dev.off()

#save file
write.csv(master_meth,file='WGBS_DMA_DMRS_5x_bonferonni.csv',row.names=FALSE)
```

## Ultimate Overlap

```
setwd("E:/Research/scratch/crow_hybrid_paper/CGI/all")

library(methylKit)
library(lme4)
library(Matrix)
library(emmeans)
library(multcomp)
library(matrixStats)
library(ggplot2)
library(gridExtra)
library(HDInterval)
library(dplyr)
library(stringr)
library(maditr)
library(viridis)

#import truth set DMRs
truth_dmr_CG <- read.csv('CG_DMA_DMRs_5x_FDR.csv',header=TRUE)
head(truth_dmr_CG)

truth_dmr_WGBS <- read.csv('WGBS_DMA_DMRs_5x_FDR.csv',header=TRUE)
head(truth_dmr_WGBS)

truth <- merge(truth_dmr_WGBS,truth_dmr_CG,by='site',all=TRUE)

hz <- read.csv('HZ_DMA_DMRs_5x_FDR.csv',header=TRUE)

dmr_min <- merge(hz,truth,by='site')
dmr_min$DMR.x <- str_replace_na(dmr_min$DMR.x, replacement = "MISSING")
dmr_min$DMR.y <- str_replace_na(dmr_min$DMR.y, replacement = "MISSING")

### try full dataset overlap overlap 
dmr_min %>% dplyr::count(DMR,DMR.x,DMR.y)
dmr_min <- dmr_min %>% mutate(Overlap = ifelse(DMR == "SPECIES" & (DMR.x == "MISSING" | DMR.x == "UNKNOWN") &
                              (DMR.y == "SEX_SPECIES" | DMR.y == "SPECIES" | DMR.y == "SPECIES_STAGE" | DMR.y == "ALL" ) , "SPECIES",
                              ifelse(DMR == "SPECIES" & (DMR.y == "MISSING" | DMR.y == "UNKNOWN") &
                                       (DMR.x == "SEX_SPECIES" | DMR.x == "SPECIES" | DMR.x == "SPECIES_TISSUE" | DMR.x == "ALL" ) , "SPECIES",
                                     ifelse(DMR == "UNKNOWN" & (DMR.y == "MISSING" | DMR.y == "UNKNOWN") &
                                              (DMR.x == "SEX_SPECIES" | DMR.x == "SEX" | DMR.x == "ALL" ) , "SEX",
                                            ifelse(DMR == "UNKNOWN" & (DMR.x == "MISSING" | DMR.x == "UNKNOWN") &
                                                     (DMR.y == "SEX_SPECIES" | DMR.y == "SEX" | DMR.y == "ALL" ) , "SEX",
                                               ifelse(DMR == "UNKNOWN" & 
                                                        (DMR.y == "SPECIES_STAGE" | DMR.y == "STAGE") & 
                                                        (DMR.x == "UNKNOWN") , "STAGE",
                                                      ifelse(DMR == "UNKNOWN" & 
                                                               (DMR.y == "UNKNOWN") & 
                                                               (DMR.x == "TISSUE" | DMR.x == "SPECIES_TISSUE") , "TISSUE",
                                                             ifelse(DMR == "UNKNOWN" & 
                                                                      (DMR.y == "UNKNOWN") & 
                                                                      (DMR.x == "UNKNOWN") , "UNKNOWN","UNCLASSIFIED"))))))))

dmr_min %>% dplyr::count(Overlap)

dmr_min <- dmr_min[!grepl('UNCLASSIFIED',dmr_min$Overlap),]
nrow(dmr_min)
dmr_min %>% dplyr::count(Overlap)
#bootstraps
nboots <- 100
HZ_species <- NULL
HZ_sex <- NULL
HZ_tissue <- NULL
HZ_stage <- NULL
HZ_other <- NULL
HZ_pval <- NULL
HZ_means <- NULL
HZ_noint <- NULL

counter = 0

for (j in seq(1, nboots)) {
  cat('\n Running iteration: ',counter,'\n')
  #obtain stratified sample
  counter=counter+1
  strat_sample <- dmr_min %>%
    group_by(Overlap) %>%
    sample_n(size=4,replace=FALSE)
  
  #re-form model 
  datmeth <- read.csv("HZ_GLMM_Input_5x.csv",header=TRUE)
  datmeth$site <- datmeth$CGI
  head(datmeth)
  nrow(datmeth)/24
  
  parents <- str_subset(unique(datmeth$ID),'D_Ko|S_Up')
  hybrids <- str_subset(unique(datmeth$ID),'D_Ko|S_Up',negate=TRUE)
  hyb3 <- sample(hybrids,3,replace=FALSE)
  full <- c(parents,hyb3)
  
  #only keep windows identified
  mod_input1 <- datmeth %>% filter_at(.vars = vars(site),
                                      .vars_predicate = any_vars(str_detect(. , paste0("^(", paste(strat_sample$site, collapse = "$|^"), ")"))))
  
  mod_input2 <- mod_input1[grepl(paste0(full,collapse="$|^"),mod_input1$ID),]
  #keep levels 
  hzmrg <- strat_sample[grepl('^site|Overlap',names(strat_sample))]
  mod_input <- merge(mod_input2,hzmrg)
  nrow(mod_input)/9
  datmeth <- mod_input
  
  #set these as factor variables 
  datmeth$Classification <- as.factor(datmeth$Classification)
  datmeth$Locality <- as.factor(datmeth$Locality)
  datmeth$Year <- as.factor(datmeth$Year)
  datmeth$ID <- as.factor(datmeth$ID)
  
  #and scale the continuous variables
  datmeth <- datmeth %>% mutate_at(c('Hybrid_Index'), ~(scale(.) %>% as.vector))
  
  plotDat <- NULL 
  tryCatch({
    for (i in seq(1, nrow(datmeth), by = 9)) {
      
      CpGsite <- datmeth[i:(i + 8), ]
      
      glmm_CpG <- glmer(cbind(MethCounts, UnMethCounts) ~
                          Classification + (1 | Year) ,data = CpGsite, family = "binomial",
                        control = glmerControl(optimizer = "bobyqa", boundary.tol = 0.01,
                                               optCtrl = list(maxfun = 2e+08)))
      dats <- emmip(glmm_CpG, ~ Classification , type = "response", CIs = TRUE,plotit=FALSE) 
      site <- CpGsite[1:3,1]
      level <- CpGsite[1:3,11]
      group <- dats[1]
      mean <- dats[2]
      se <- dats[3]
      loops <- cbind(site,level,group,mean,se)
      plotDat <- rbind(plotDat,loops)
      
    }
    
    plotDat
    plotDat_LongSE <- dcast(plotDat, site ~ Classification, value.var = 'SE')
    plotDat_LongMEAN <- dcast(plotDat, site ~ Classification, value.var = 'yvar')
    names(plotDat_LongMEAN) <- c('site','Corone','Hybrid','Cornix')
    plotDat_LongMEAN$ParentalMidPoint <- ((plotDat_LongMEAN$Cornix + plotDat_LongMEAN$Corone)/2)
    plotDat_LongMEAN$HybDistance <- (plotDat_LongMEAN$Hybrid - plotDat_LongMEAN$ParentalMidPoint)
    
    #add simple intermediacy factor 
    plotDat_LongMEAN <- plotDat_LongMEAN %>% mutate(Intermediate = ifelse(Hybrid > Corone & Hybrid < Cornix | 
                                                                            Hybrid < Corone & Hybrid > Cornix, "INTERMEDIATE", "OTHER"))
    hyb_level <- plotDat_LongMEAN[,c('site','Intermediate','HybDistance')]
    for_tab <- merge(hyb_level,plotDat,by='site')
    k_tab <- unique(for_tab[,c('site','Intermediate','level','HybDistance')])
    
    #calculate distance between hybrids and the parental midpoint 
    hyb_means <- as.data.frame(k_tab %>% group_by(level) %>% dplyr::summarize(Mean = mean(HybDistance)))
    hyb_means$ITERATION <- counter
    
    #test significant differences in counts, if there are no intermediate events, this will add a column of zeroes
    con1 <- table(k_tab$level,k_tab$Intermediate)
    con1
    
    tryCatch({
      res <- fisher.test(con1)
      con2 <- table(k_tab$level,k_tab$HybDistance)
      con2
      
      spec <-  data.frame(ITERATION = counter, LEVEL = 'SPECIES',INTERMEDIATE = con1[2,1], MIXED = con1[2,2])
      stag <-  data.frame(ITERATION = counter, LEVEL = 'STAGE',INTERMEDIATE = con1[3,1], MIXED = con1[3,2])
      tiss <-  data.frame(ITERATION = counter, LEVEL = 'TISSUE',INTERMEDIATE = con1[4,1], MIXED = con1[4,2])
      sexx <-  data.frame(ITERATION = counter, LEVEL = 'SEX',INTERMEDIATE = con1[1,1], MIXED = con1[1,2])
      other <-  data.frame(ITERATION = counter, LEVEL = 'UNKNOWN',INTERMEDIATE = con1[5,1], MIXED = con1[5,2])
      pvals <-  data.frame(pvalue = res$p.value)
      
      
      #bind loops
      HZ_species <- rbind(HZ_species, spec)
      HZ_stage <- rbind(HZ_stage, stag)
      HZ_sex <- rbind(HZ_sex, sexx)
      HZ_tissue <- rbind(HZ_tissue, tiss)
      HZ_other <- rbind(HZ_other, other)
      HZ_pval <- rbind(HZ_pval,pvals)
      HZ_means <- rbind(HZ_means,hyb_means)
      
    }, error=function(e){HZ_noint <- rbind(HZ_noint,counter);cat(counter,"BAD \n")})
    
  }, error=function(e){cat(counter,": FAILED \n")})
  
  
}

HZ_distance <- HZ_means
HZ_distance$level <- factor(HZ_distance$level, levels = c('UNKNOWN','SEX','STAGE','TISSUE', 'SPECIES'))
a <- ggplot(HZ_distance,aes(x=level,fill=level, y=Mean))+
  geom_boxplot(show.legend=FALSE)+ylab('Distance to Parentals')+xlab('')+
  scale_fill_manual(values=viridis(6))+
  stat_boxplot(geom ='errorbar', width = 0.25) +
  geom_hline(yintercept=0,lty=2,lwd=1,col="black")+
  theme_classic()
#coord_cartesian(ylim=c(-1,1))
a


bootdat <- rbind(HZ_species,HZ_stage,HZ_sex,HZ_other)
bootdat

meltboot <- melt(bootdat,id.vars = c('LEVEL','ITERATION'))

meltboot$LEVEL <- factor(meltboot$LEVEL, levels = c('UNKNOWN','SEX','STAGE', 'SPECIES'))
b <- ggplot(meltboot,aes(x=LEVEL,fill=variable,y=value))+
  geom_boxplot()+ylab('Number of Replicates')+xlab('')+
  scale_fill_manual(values=c('purple','gray'))+
  theme_classic()
b

c <- ggplot(HZ_pval, aes(x=pvalue))+
  geom_histogram(aes(y = ..density..),bins=5) +
  geom_density(alpha = 0.1, fill = "turquoise")+
  ylab("Density")+xlab("p-Value")+
  theme_classic()
c

#save it 
jpeg("plots/FULL-BootstrapDistance.jpg",units="in",res=300, height=5, width=7)
theme_set(theme_classic(base_size = 16))
grid.arrange(a,nrow=1)
dev.off()


#save it 
jpeg("plots/FULL-Bootstrap.jpg",units="in",res=300, height=5, width=7)
theme_set(theme_classic(base_size = 16))
grid.arrange(b,nrow=1)
dev.off()

#save it 
jpeg("plots/FULL-Bootstrap_pvals.jpg",units="in",res=300, height=5, width=4)
theme_set(theme_classic(base_size = 16))
grid.arrange(c,nrow=1)
dev.off()
```

## Hybrid Intermediacy

```
library(methylKit)
library(lme4)
library(Matrix)
library(emmeans)
library(multcomp)
library(tidyr)
library(ggplot2)
library(psych)
library(gridExtra)
library(HDInterval)
library(dplyr)
library(data.table)
library(stringr)
library(LICORS)
library(RColorBrewer)
library(viridis)

setwd("E:/Research/scratch/crow_hybrid_paper/CGI/all")

#import truth set DMRs
truth_dmr_CG <- read.csv('CG_DMA_DMRs_5x_FDR.csv',header=TRUE)
head(truth_dmr_CG)

truth_dmr_WGBS <- read.csv('WGBS_DMA_DMRs_5x_FDR.csv',header=TRUE)
head(truth_dmr_WGBS)

truth <- merge(truth_dmr_WGBS,truth_dmr_CG,by='site',all=TRUE)

hz <- read.csv('HZ_DMA_DMRs_5x_FDR.csv',header=TRUE)

dmr_min <- merge(hz,truth,by='site')
dmr_min$DMR.x <- str_replace_na(dmr_min$DMR.x, replacement = "MISSING")
dmr_min$DMR.y <- str_replace_na(dmr_min$DMR.y, replacement = "MISSING")

### try full dataset overlap overlap 
dmr_min %>% dplyr::count(DMR,DMR.x,DMR.y)
dmr_min <- dmr_min %>% mutate(Overlap = ifelse(DMR == "SPECIES" & (DMR.x == "MISSING" | DMR.x == "UNKNOWN") &
                                                 (DMR.y == "SEX_SPECIES" | DMR.y == "SPECIES" | DMR.y == "SPECIES_STAGE" | DMR.y == "ALL" ) , "SPECIES",
                                               ifelse(DMR == "SPECIES" & (DMR.y == "MISSING" | DMR.y == "UNKNOWN") &
                                                        (DMR.x == "SEX_SPECIES" | DMR.x == "SPECIES" | DMR.x == "SPECIES_TISSUE" | DMR.x == "ALL" ) , "SPECIES",
                                                      ifelse(DMR == "UNKNOWN" & (DMR.y == "MISSING" | DMR.y == "UNKNOWN") &
                                                               (DMR.x == "SEX_SPECIES" | DMR.x == "SEX" | DMR.x == "ALL" ) , "SEX",
                                                             ifelse(DMR == "UNKNOWN" & (DMR.x == "MISSING" | DMR.x == "UNKNOWN") &
                                                                      (DMR.y == "SEX_SPECIES" | DMR.y == "SEX" | DMR.y == "ALL" ) , "SEX",
                                                                    ifelse(DMR == "UNKNOWN" & 
                                                                             (DMR.y == "SPECIES_STAGE" | DMR.y == "STAGE") & 
                                                                             (DMR.x == "UNKNOWN") , "STAGE",
                                                                           ifelse(DMR == "UNKNOWN" & 
                                                                                    (DMR.y == "UNKNOWN") & 
                                                                                    (DMR.x == "TISSUE" | DMR.x == "SPECIES_TISSUE") , "TISSUE",
                                                                                  ifelse(DMR == "UNKNOWN" & 
                                                                                           (DMR.y == "UNKNOWN") & 
                                                                                           (DMR.x == "UNKNOWN") , "UNKNOWN","UNCLASSIFIED"))))))))

dmr_min %>% dplyr::count(Overlap)


#save for manhattan
metadat <- read.table('metadata_hz.txt',header=TRUE)

hz_dats <- dmr_min[grepl('^site$|^S_|^D_|^Overlap',names(dmr_min))]
hz_dats <- hz_dats[!grepl('UNCLASSIFIED',hz_dats$Overlap),]
hz_dats %>% dplyr::count(Overlap)

#randomly sample the number of windows according to the number of species windows
strat_sample <- hz_dats %>%
  group_by(Overlap) %>%
  sample_n(size=5,replace=FALSE)

panel_labs <- unique(strat_sample$site)
plot_labs <- rep(seq(1,5),5)
labs <- cbind(panel_labs,plot_labs)
colnames(labs)[colnames(labs) == c('panel_labs','plot_labs')] <- c('site','siteID')
plotdis <- merge(strat_sample,labs,by='site')

hz.m1 = melt(strat_sample, id.vars = c('site','Overlap'),
             measure.vars = names(hz[grepl('^S_|^D_',names(hz_dmrs))]))

colnames(metadat)[colnames(metadat) == 'Library'] <- 'variable'
hzp <- merge(hz.m1,metadat)
hzp2 <- merge(hzp,labs,by='site')



hzp2$Classification <- as.factor(hzp2$Classification)

hyb <- ggplot(hzp2,aes(x=Classification,y=value,fill=Classification))+
  geom_boxplot(varwidth=FALSE,alpha=0.75)+
  stat_boxplot(geom ='errorbar', width = 0.6) +
  theme_bw()+facet_grid(Overlap~siteID)+
  scale_fill_manual(values=c('black','purple','grey60'))+ylab('% 5mC')+
  coord_cartesian(ylim=c(0,1))+theme(axis.text.x=element_blank(),axis.title.x=element_blank())

hyb

#save it 
jpeg("plots/HZ-Intermediacy.jpg",units="in",res=300, height=8, width=10)
theme_set(theme_classic(base_size = 16))
grid.arrange(hyb,nrow=1)
dev.off()

```

## CG-HZ Overlap

```
setwd("E:/Research/scratch/crow_hybrid_paper/CGI/all")

library(methylKit)
library(lme4)
library(Matrix)
library(emmeans)
library(multcomp)
library(matrixStats)
library(ggplot2)
library(gridExtra)
library(HDInterval)
library(dplyr)
library(stringr)
library(maditr)
library(viridis)

#import truth set DMRs
truth_dmr_CG <- read.csv('CG_DMA_DMRs_5x_bonferonni.csv',header=TRUE)
head(truth_dmr_CG)

hz <- read.csv('HZ_DMA_DMRs_5x_bonferonni.csv',header=TRUE)


hz_cg <- merge(hz,truth_dmr_CG,by='site')

### try CG-HZ overlap 
hz_cg %>% dplyr::count(DMR.x,DMR.y)
dmr_min <- hz_cg
dmr_min <- dmr_min %>% mutate(Overlap = ifelse(DMR.x == "SPECIES" & 
                                                 (DMR.y == "SEX_SPECIES" | DMR.y == "SPECIES" | DMR.y == "SPECIES_STAGE" | DMR.y == "ALL"), "SPECIES",
                                               ifelse(DMR.x == "UNKNOWN" & 
                                                        (DMR.y == "SEX_SPECIES" | DMR.y == "SEX" | DMR.y == "ALL" ) , "SEX",
                                                      ifelse(DMR.x == "UNKNOWN" & 
                                                               (DMR.y == "SPECIES_STAGE" | DMR.y == "STAGE") , "STAGE",
                                                             ifelse(DMR.x == "UNKNOWN" & 
                                                                      (DMR.y == "UNKNOWN") , "UNKNOWN","UNCLASSIFIED")))))

dmr_min %>% dplyr::count(Overlap)

#bootstraps
nboots <- 100
HZ_species <- NULL
HZ_sex <- NULL
HZ_stage <- NULL
HZ_unclassified <- NULL
HZ_other <- NULL
HZ_pval <- NULL
HZ_means <- NULL
HZ_noint <- NULL

counter = 0

for (j in seq(1, nboots)) {
  cat('\n Running iteration: ',counter,'\n')
  #obtain stratified sample
  counter=counter+1
  strat_sample <- dmr_min %>%
    group_by(Overlap) %>%
    sample_n(size=2,replace=FALSE)
  
  #re-form model 
  datmeth <- read.csv("HZ_GLMM_Input_5x.csv",header=TRUE)
  datmeth$site <- datmeth$CGI
  head(datmeth)
  nrow(datmeth)/24
  
  parents <- str_subset(unique(datmeth$ID),'D_Ko|S_Up')
  hybrids <- str_subset(unique(datmeth$ID),'D_Ko|S_Up',negate=TRUE)
  hyb3 <- sample(hybrids,3,replace=FALSE)
  full <- c(parents,hyb3)
  
  #only keep windows identified
  mod_input1 <- datmeth %>% filter_at(.vars = vars(site),
                                      .vars_predicate = any_vars(str_detect(. , paste0("^(", paste(strat_sample$site, collapse = "$|^"), ")"))))
  
  mod_input2 <- mod_input1[grepl(paste0(full,collapse="$|^"),mod_input1$ID),]
  #keep levels 
  hzmrg <- strat_sample[grepl('^site|Overlap',names(strat_sample))]
  mod_input <- merge(mod_input2,hzmrg)
  nrow(mod_input)/9
  datmeth <- mod_input
  
  #set these as factor variables 
  datmeth$Classification <- as.factor(datmeth$Classification)
  datmeth$Locality <- as.factor(datmeth$Locality)
  datmeth$Year <- as.factor(datmeth$Year)
  datmeth$ID <- as.factor(datmeth$ID)
  
  #and scale the continuous variables
  datmeth <- datmeth %>% mutate_at(c('Hybrid_Index'), ~(scale(.) %>% as.vector))
  
  plotDat <- NULL 
  tryCatch({
    for (i in seq(1, nrow(datmeth), by = 9)) {
    
    CpGsite <- datmeth[i:(i + 8), ]
    
    glmm_CpG <- glmer(cbind(MethCounts, UnMethCounts) ~
                        Classification + (1 | Year) ,data = CpGsite, family = "binomial",
                      control = glmerControl(optimizer = "bobyqa", boundary.tol = 0.01,
                                             optCtrl = list(maxfun = 2e+08)))
    dats <- emmip(glmm_CpG, ~ Classification , type = "response", CIs = TRUE,plotit=FALSE) 
    site <- CpGsite[1:3,1]
    level <- CpGsite[1:3,11]
    group <- dats[1]
    mean <- dats[2]
    se <- dats[3]
    loops <- cbind(site,level,group,mean,se)
    plotDat <- rbind(plotDat,loops)
    
    }
    
  plotDat
  plotDat_LongSE <- dcast(plotDat, site ~ Classification, value.var = 'SE')
  plotDat_LongMEAN <- dcast(plotDat, site ~ Classification, value.var = 'yvar')
  names(plotDat_LongMEAN) <- c('site','Corone','Hybrid','Cornix')
  plotDat_LongMEAN$ParentalMidPoint <- ((plotDat_LongMEAN$Cornix + plotDat_LongMEAN$Corone)/2)
  plotDat_LongMEAN$HybDistance <- (plotDat_LongMEAN$Hybrid - plotDat_LongMEAN$ParentalMidPoint)
  
  #add simple intermediacy factor 
  plotDat_LongMEAN <- plotDat_LongMEAN %>% mutate(Intermediate = ifelse(Hybrid > Corone & Hybrid < Cornix | 
                                                                          Hybrid < Corone & Hybrid > Cornix, "INTERMEDIATE", "OTHER"))
  hyb_level <- plotDat_LongMEAN[,c('site','Intermediate','HybDistance')]
  for_tab <- merge(hyb_level,plotDat,by='site')
  k_tab <- unique(for_tab[,c('site','Intermediate','level','HybDistance')])
  
  #calculate distance between hybrids and the parental midpoint 
  hyb_means <- as.data.frame(k_tab %>% group_by(level) %>% dplyr::summarize(Mean = mean(HybDistance)))
  hyb_means$ITERATION <- counter
  
  #test significant differences in counts, if there are no intermediate events, this will add a column of zeroes
  con1 <- table(k_tab$level,k_tab$Intermediate)
  con1
  
  tryCatch({
    res <- fisher.test(con1)
    con2 <- table(k_tab$level,k_tab$HybDistance)
    con2
    
    spec <-  data.frame(ITERATION = counter, LEVEL = 'SPECIES',INTERMEDIATE = con1[2,1], MIXED = con1[2,2])
    stag <-  data.frame(ITERATION = counter, LEVEL = 'STAGE',INTERMEDIATE = con1[3,1], MIXED = con1[3,2])
    fals <-  data.frame(ITERATION = counter, LEVEL = 'UNCLASSIFIED',INTERMEDIATE = con1[4,1], MIXED = con1[4,2])
    sexx <-  data.frame(ITERATION = counter, LEVEL = 'SEX',INTERMEDIATE = con1[1,1], MIXED = con1[1,2])
    other <-  data.frame(ITERATION = counter, LEVEL = 'UNKNOWN',INTERMEDIATE = con1[5,1], MIXED = con1[5,2])
    pvals <-  data.frame(pvalue = res$p.value)
    
    
    #bind loops
    HZ_species <- rbind(HZ_species, spec)
    HZ_stage <- rbind(HZ_stage, stag)
    HZ_sex <- rbind(HZ_sex, sexx)
    HZ_unclassified <- rbind(HZ_unclassified, fals)
    HZ_other <- rbind(HZ_other, other)
    HZ_pval <- rbind(HZ_pval,pvals)
    HZ_means <- rbind(HZ_means,hyb_means)
    
  }, error=function(e){HZ_noint <- rbind(HZ_noint,counter);cat(counter,"BAD \n")})
  
    }, error=function(e){cat(counter,": FAILED \n")})
  
  
}

HZ_distance <- HZ_means
HZ_distance <- HZ_distance[!grepl('UNCLASSIFIED',HZ_distance$level),]
HZ_distance$level <- factor(HZ_distance$level, levels = c('UNKNOWN','SEX','STAGE', 'SPECIES'))
a <- ggplot(HZ_distance,aes(x=level,fill=level, y=Mean))+
  geom_boxplot(show.legend=FALSE)+ylab('Distance to Parentals')+xlab('')+
  scale_fill_manual(values=viridis(6))+
  stat_boxplot(geom ='errorbar', width = 0.25) +
  geom_hline(yintercept=0,lty=2,lwd=1,col="black")+
  theme_classic()
#coord_cartesian(ylim=c(-1,1))
a


  bootdat <- rbind(HZ_species,HZ_stage,HZ_sex,HZ_other)
bootdat

meltboot <- melt(bootdat,id.vars = c('LEVEL','ITERATION'))

meltboot$LEVEL <- factor(meltboot$LEVEL, levels = c('UNKNOWN','SEX','STAGE', 'SPECIES'))
b <- ggplot(meltboot,aes(x=LEVEL,fill=variable,y=value))+
  geom_boxplot()+ylab('Number of Replicates')+xlab('')+
  scale_fill_manual(values=c('purple','gray'))+
  theme_classic()
b

c <- ggplot(HZ_pval, aes(x=pvalue))+
  geom_histogram(aes(y = ..density..),bins=5) +
  geom_density(alpha = 0.1, fill = "turquoise")+
  ylab("Density")+xlab("p-Value")+
  theme_classic()
c

#save it 
jpeg("plots/HZ-CG-BootstrapDistance.jpg",units="in",res=300, height=5, width=7)
theme_set(theme_classic(base_size = 16))
grid.arrange(a,nrow=1)
dev.off()


#save it 
jpeg("plots/HZ-CG-Bootstrap.jpg",units="in",res=300, height=5, width=7)
theme_set(theme_classic(base_size = 16))
grid.arrange(b,nrow=1)
dev.off()

#save it 
jpeg("plots/HZ-CG-Bootstrap_pvals.jpg",units="in",res=300, height=5, width=4)
theme_set(theme_classic(base_size = 16))
grid.arrange(c,nrow=1)
dev.off()
```

## WGBS-HZ Overlap

```
setwd("E:/Research/scratch/crow_hybrid_paper/CGI/all")

library(methylKit)
library(lme4)
library(Matrix)
library(emmeans)
library(multcomp)
library(matrixStats)
library(ggplot2)
library(gridExtra)
library(HDInterval)
library(dplyr)
library(stringr)
library(maditr)
library(viridis)

#import truth set DMRs
truth_dmr_WGBS <- read.csv('WGBS_DMA_DMRs_5x_bonferonni.csv',header=TRUE)
head(truth_dmr_WGBS)

hz <- read.csv('HZ_DMA_DMRs_5x_bonferonni.csv',header=TRUE)


hz_wgbs <- merge(hz,truth_dmr_WGBS,by='site')

### try WGBS-HZ overlap 
dmr_min <- hz_wgbs
hz_wgbs %>% dplyr::count(DMR.x,DMR.y)
dmr_min <- dmr_min %>% mutate(Overlap = ifelse(DMR.x == "SPECIES" & 
                                                 (DMR.y == "SEX_SPECIES" | DMR.y == "SPECIES" | DMR.y == "SPECIES_TISSUE" | DMR.y == "ALL"), "SPECIES",
                                               ifelse(DMR.x == "UNKNOWN" & 
                                                        (DMR.y == "SEX_SPECIES" | DMR.y == "SEX" | DMR.y == "ALL" ) , "SEX",
                                                      ifelse(DMR.x == "UNKNOWN" & 
                                                               (DMR.y == "SPECIES_TISSUE" | DMR.y == "TISSUE") , "TISSUE",
                                                             ifelse(DMR.x == "UNKNOWN" & 
                                                                      (DMR.y == "UNKNOWN") , "UNKNOWN","UNCLASSIFIED")))))

dmr_min %>% dplyr::count(Overlap)

#bootstraps
nboots <- 100
HZ_species <- NULL
HZ_sex <- NULL
HZ_tissue <- NULL
HZ_unclassified <- NULL
HZ_other <- NULL
HZ_pval <- NULL
HZ_means <- NULL
HZ_noint <- NULL

counter = 0

for (j in seq(1, nboots)) {
  cat('\n Running iteration: ',counter,'\n')
  #obtain stratified sample
  counter=counter+1
  strat_sample <- dmr_min %>%
    group_by(Overlap) %>%
    sample_n(size=2,replace=FALSE)
  
  #re-form model 
  datmeth <- read.csv("HZ_GLMM_Input_5x.csv",header=TRUE)
  datmeth$site <- datmeth$CGI
  head(datmeth)
  nrow(datmeth)/24
  
  parents <- str_subset(unique(datmeth$ID),'D_Ko|S_Up')
  hybrids <- str_subset(unique(datmeth$ID),'D_Ko|S_Up',negate=TRUE)
  hyb3 <- sample(hybrids,3,replace=FALSE)
  full <- c(parents,hyb3)
  
  #only keep windows identified
  mod_input1 <- datmeth %>% filter_at(.vars = vars(site),
                                      .vars_predicate = any_vars(str_detect(. , paste0("^(", paste(strat_sample$site, collapse = "$|^"), ")"))))
  
  mod_input2 <- mod_input1[grepl(paste0(full,collapse="$|^"),mod_input1$ID),]
  #keep levels 
  hzmrg <- strat_sample[grepl('^site|Overlap',names(strat_sample))]
  mod_input <- merge(mod_input2,hzmrg)
  nrow(mod_input)/9
  datmeth <- mod_input
  
  #set these as factor variables 
  datmeth$Classification <- as.factor(datmeth$Classification)
  datmeth$Locality <- as.factor(datmeth$Locality)
  datmeth$Year <- as.factor(datmeth$Year)
  datmeth$ID <- as.factor(datmeth$ID)
  
  #and scale the continuous variables
  datmeth <- datmeth %>% mutate_at(c('Hybrid_Index'), ~(scale(.) %>% as.vector))
  
  plotDat <- NULL 
  tryCatch({
    for (i in seq(1, nrow(datmeth), by = 9)) {
      CpGsite <- datmeth[i:(i + 8), ]
      
      glmm_CpG <- glmer(cbind(MethCounts, UnMethCounts) ~
                          Classification + (1 | Year) ,data = CpGsite, family = "binomial",
                        control = glmerControl(optimizer = "bobyqa", boundary.tol = 0.01,
                                               optCtrl = list(maxfun = 2e+08)))
      dats <- emmip(glmm_CpG, ~ Classification , type = "response", CIs = TRUE,plotit=FALSE) 
      site <- CpGsite[1:3,1]
      level <- CpGsite[1:3,11]
      group <- dats[1]
      mean <- dats[2]
      se <- dats[3]
      loops <- cbind(site,level,group,mean,se)
      plotDat <- rbind(plotDat,loops)
    }
  
  plotDat
  plotDat_LongSE <- dcast(plotDat, site ~ Classification, value.var = 'SE')
  plotDat_LongMEAN <- dcast(plotDat, site ~ Classification, value.var = 'yvar')
  names(plotDat_LongMEAN) <- c('site','Corone','Hybrid','Cornix')
  plotDat_LongMEAN$ParentalMidPoint <- ((plotDat_LongMEAN$Cornix + plotDat_LongMEAN$Corone)/2)
  plotDat_LongMEAN$HybDistance <- (plotDat_LongMEAN$Hybrid - plotDat_LongMEAN$ParentalMidPoint)
  
  #add simple intermediacy factor 
  plotDat_LongMEAN <- plotDat_LongMEAN %>% mutate(Intermediate = ifelse(Hybrid > Corone & Hybrid < Cornix | 
                                                                          Hybrid < Corone & Hybrid > Cornix, "INTERMEDIATE", "OTHER"))
  hyb_level <- plotDat_LongMEAN[,c('site','Intermediate','HybDistance')]
  for_tab <- merge(hyb_level,plotDat,by='site')
  k_tab <- unique(for_tab[,c('site','Intermediate','level','HybDistance')])
  
  #calculate distance between hybrids and the parental midpoint 
  hyb_means <- as.data.frame(k_tab %>% group_by(level) %>% dplyr::summarize(Mean = mean(HybDistance)))
  hyb_means$ITERATION <- counter
  
  #test significant differences in counts, if there are no intermediate events, this will add a column of zeroes
  con1 <- table(k_tab$level,k_tab$Intermediate)
  con1
  
  tryCatch({
    res <- fisher.test(con1)
    con2 <- table(k_tab$level,k_tab$HybDistance)
    con2
    
    spec <-  data.frame(ITERATION = counter, LEVEL = 'SPECIES',INTERMEDIATE = con1[2,1], MIXED = con1[2,2])
    tiss <-  data.frame(ITERATION = counter, LEVEL = 'TISSUE',INTERMEDIATE = con1[3,1], MIXED = con1[3,2])
    fals <-  data.frame(ITERATION = counter, LEVEL = 'UNCLASSIFIED',INTERMEDIATE = con1[4,1], MIXED = con1[4,2])
    sexx <-  data.frame(ITERATION = counter, LEVEL = 'SEX',INTERMEDIATE = con1[1,1], MIXED = con1[1,2])
    other <-  data.frame(ITERATION = counter, LEVEL = 'UNKNOWN',INTERMEDIATE = con1[5,1], MIXED = con1[5,2])
    pvals <-  data.frame(pvalue = res$p.value)
    
    
    #bind loops
    HZ_species <- rbind(HZ_species, spec)
    HZ_tissue <- rbind(HZ_tissue, tiss)
    HZ_sex <- rbind(HZ_sex, sexx)
    HZ_unclassified <- rbind(HZ_unclassified, fals)
    HZ_other <- rbind(HZ_other, other)
    HZ_pval <- rbind(HZ_pval,pvals)
    HZ_means <- rbind(HZ_means,hyb_means)
    
  }, error=function(e){HZ_noint <- rbind(HZ_noint,counter);cat(counter,"BAD \n")})
  
    }, error=function(e){cat(counter,": FAILED \n")})
  
  
}

HZ_distance <- HZ_means
HZ_distance <- HZ_distance[!grepl('UNCLASSIFIED',HZ_distance$level),]
HZ_distance$level <- factor(HZ_distance$level, levels = c('UNKNOWN','SEX','TISSUE', 'SPECIES'))
a <- ggplot(HZ_distance,aes(x=level,fill=level, y=Mean))+
  geom_boxplot(show.legend=FALSE)+ylab('Distance to Parentals')+xlab('')+
  scale_fill_manual(values=viridis(6))+
  stat_boxplot(geom ='errorbar', width = 0.25) +
  geom_hline(yintercept=0,lty=2,lwd=1,col="black")+
  theme_classic()
#coord_cartesian(ylim=c(-1,1))
a


bootdat <- rbind(HZ_species,HZ_tissue,HZ_sex,HZ_other)
bootdat

meltboot <- melt(bootdat,id.vars = c('LEVEL','ITERATION'))

meltboot$LEVEL <- factor(meltboot$LEVEL, levels = c('UNKNOWN','SEX','TISSUE', 'SPECIES'))
b <- ggplot(meltboot,aes(x=LEVEL,fill=variable,y=value))+
  geom_boxplot()+ylab('Number of Replicates')+xlab('')+
  scale_fill_manual(values=c('purple','gray'))+
  theme_classic()
b

c <- ggplot(HZ_pval, aes(x=pvalue))+
  geom_histogram(aes(y = ..density..),bins=5) +
  geom_density(alpha = 0.1, fill = "turquoise")+
  ylab("Density")+xlab("p-Value")+
  theme_classic()
c

#save it 
jpeg("plots/HZ-CG-BootstrapDistance.jpg",units="in",res=300, height=5, width=7)
theme_set(theme_classic(base_size = 16))
grid.arrange(a,nrow=1)
dev.off()


#save it 
jpeg("plots/HZ-CG-Bootstrap.jpg",units="in",res=300, height=5, width=7)
theme_set(theme_classic(base_size = 16))
grid.arrange(b,nrow=1)
dev.off()

#save it 
jpeg("plots/HZ-CG-Bootstrap_pvals.jpg",units="in",res=300, height=5, width=4)
theme_set(theme_classic(base_size = 16))
grid.arrange(c,nrow=1)
dev.off()
```

# 
