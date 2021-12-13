### 'Broodmix' neonicotinoid resilience analysis
# For review of Bartlett et al.
# Analysis by Lewis J Bartlett
# lewis.bartlett@uga.edu

####

# Read in data

DD <- 'D:/Google Drive/PostDoc2/BroodMix/Empirical Paper/Datasets/Cleaned/'

DFs <- list.files(path = DD,
                  pattern = NULL, all.files = FALSE,
                  full.names = FALSE, recursive = FALSE,
                  ignore.case = FALSE, include.dirs = FALSE)

for(DF in DFs){
  
  tmp <- read.csv(file = paste0(DD,DF),
                  header = T,
                  stringsAsFactors = F)
  
  assign(x = strsplit(DF, split = '.csv')[[1]][1],
         value = tmp)
  
  rm(tmp)
  
}

rm(list = c('DD','DF','DFs'))

SimData <- read.csv(file = 'D:/Google Drive/PostDoc2/BroodMix/Empirical Paper/Datasets/PolyMaxPredictions/SimData.csv',
                     header = T)

#SimData <- read.csv(file = 'C:/Users/ljb87745/Google Drive/PostDoc2/BroodMix/Empirical Paper/Datasets/PolyMaxPredictions/SimData.csv',
#                    header = T)

SimData$Emphasis <- SimData$Emphasis/max(SimData$Emphasis)

# Data loaded and ready for arrangement for analysis

#
## Assemble necessary analysis data 'on the fly' per task/response variable
# 

# load necessary packages

library(afex)

# see https://m-clark.github.io/mixed-models-with-R/extensions.html as a good example of how nested + crossed effects work

# Bind together reference sheets

AlabamaTreatRef$Site <- 'Alabama'
DelawareTreatRef$Site <- 'Delaware'
GeorgiaTreatRef$Site <- 'Georgia'

TreatRefAll <- rbind(AlabamaTreatRef, DelawareTreatRef, GeorgiaTreatRef)

# Let's begin.

## Aggression

AlabamaAggression$Site <- 'Alabama'
DelawareAggression$Site <- 'Delaware'
GeorgiaAggression$Site <- 'Georgia'

AggDF <- rbind(AlabamaAggression, DelawareAggression, GeorgiaAggression)

AggDF$Neonic <- NA
AggDF$Mix <- NA
AggDF$UID <- NA


for(N in 1:NROW(AggDF)){
  
  AggDF$Neonic[N] <- TreatRefAll$Neonic[which(TreatRefAll$Colony == AggDF$Colony[N] & TreatRefAll$Site == AggDF$Site[N])]
  
  AggDF$Mix[N] <- TreatRefAll$Mixed[which(TreatRefAll$Colony == AggDF$Colony[N] & TreatRefAll$Site == AggDF$Site[N])]

  AggDF$UID <- paste0(AggDF$Site, AggDF$Colony)
  AggDF$SD <- paste0(AggDF$Site, AggDF$Day)
  
}


AggMod <- mixed(Recruitment ~ BeesT0 + Neonic*Mix + (1|Site) + (1|Site:UID) + (1|Site:SD),
                 method = 'KR',
                 data = na.exclude(AggDF)
)
nice(AggMod)
summary(AggMod)

# No strong evidence of effects


## Brood Survival

AlabamaBroodSurvivorshipL2$Survived <- AlabamaBroodSurvivorshipL2$Total - AlabamaBroodSurvivorshipL2$Died
DelawareBroodSurvivorshipL2$Survived <- DelawareBroodSurvivorshipL2$Total - DelawareBroodSurvivorshipL2$Died
GeorgiaBroodSurvivorshipL2$Died <- GeorgiaBroodSurvivorshipL2$Total - GeorgiaBroodSurvivorshipL2$Survived

AlabamaBroodSurvivorshipL2$Site <- 'Alabama'
DelawareBroodSurvivorshipL2$Site <- 'Delaware'
GeorgiaBroodSurvivorshipL2$Site <- 'Georgia'

BSDF <- rbind(AlabamaBroodSurvivorshipL2, DelawareBroodSurvivorshipL2, GeorgiaBroodSurvivorshipL2)

BSDF$DS <- cbind(BSDF$Died, BSDF$Survived)

BSDF <- na.exclude(BSDF)

BSDF$Neonic <- NA
BSDF$Mix <- NA
BSDF$UID <- NA

for(N in 1:NROW(BSDF)){
  
  BSDF$Neonic[N] <- TreatRefAll$Neonic[which(TreatRefAll$Colony == BSDF$Colony[N] & TreatRefAll$Site == BSDF$Site[N])]
  
  BSDF$Mix[N] <- TreatRefAll$Mixed[which(TreatRefAll$Colony == BSDF$Colony[N] & TreatRefAll$Site == BSDF$Site[N])]
  
  BSDF$UID <- paste0(BSDF$Site, BSDF$Colony)
  BSDF$SD <- paste0(BSDF$Site, BSDF$Day)
}


BSMod <- mixed(DS ~ Neonic*Mix + (1|Site) + (1|Site:UID) + (1|Site:SD),
                family = 'binomial',
                method = 'LRT',
                data = na.exclude(BSDF)
)
nice(BSMod)
summary(BSMod)

## Comb Construction

DelawareCombCon$Site <- 'Delaware'
GeorgiaCombCon$Site <- 'Georgia'

CCDF <- rbind(DelawareCombCon, GeorgiaCombCon)

CCDF <- na.exclude(CCDF)

CCDF$Neonic <- NA
CCDF$Mix <- NA
CCDF$UID <- NA

for(N in 1:NROW(CCDF)){
  
  CCDF$Neonic[N] <- TreatRefAll$Neonic[which(TreatRefAll$Colony == CCDF$Colony[N] & TreatRefAll$Site == CCDF$Site[N])]
  
  CCDF$Mix[N] <- TreatRefAll$Mixed[which(TreatRefAll$Colony == CCDF$Colony[N] & TreatRefAll$Site == CCDF$Site[N])]
  
  CCDF$UID <- paste0(CCDF$Site, CCDF$Colony)
  CCDF$SD <- paste0(CCDF$Site, CCDF$Day)
}

CCMod <- mixed(CombArea ~ Neonic*Mix + (1|Site) + (1|Site:UID) + (1|Site:SD),
                family = 'poisson',
                method = 'LRT',
                data = na.exclude(CCDF)
)
nice(CCMod)
summary(CCMod)

## Mite Drop

AlabamaMiteDrop$Site <- 'Alabama'
DelawareMiteDrop$Site <- 'Delaware'
GeorgiaMiteDrop$Site <- 'Georgia'

MDDF <- rbind(AlabamaMiteDrop, DelawareMiteDrop, GeorgiaMiteDrop)

MDDF <- MDDF[which(MDDF$Day > 0),]

MDDF$Neonic <- NA
MDDF$Mix <- NA
MDDF$UID <- NA


for(N in 1:NROW(MDDF)){
  
  MDDF$Neonic[N] <- TreatRefAll$Neonic[which(TreatRefAll$Colony == MDDF$Colony[N] & TreatRefAll$Site == MDDF$Site[N])]
  
  MDDF$Mix[N] <- TreatRefAll$Mixed[which(TreatRefAll$Colony == MDDF$Colony[N] & TreatRefAll$Site == MDDF$Site[N])]
  
  MDDF$UID <- paste0(MDDF$Site, MDDF$Colony)
  MDDF$SD <- paste0(MDDF$Site, MDDF$Day)
  
}

MDDF$MDR <- MDDF$MiteDrop/MDDF$Time

MDDF$MDR2 <- ceiling((MDDF$MiteDrop/MDDF$Time)*24)

MDMod <- mixed(MDR2 ~ Neonic*Mix + (1|Site) + (1|Site:UID) + (1|Site:SD),
                family = 'poisson',
                method = 'LRT',
                data = na.exclude(MDDF)
)
nice(MDMod)
summary(MDMod)


## Mite Wash

AlabamaMiteWash$Site <- 'Alabama'
DelawareMiteWash$Site <- 'Delaware'

MWDF <- rbind(AlabamaMiteWash, DelawareMiteWash)

MWDF <- MWDF[which(MWDF$Day > 0),]

MWDF$Neonic <- NA
MWDF$Mix <- NA
MWDF$UID <- NA


for(N in 1:NROW(MWDF)){
  
  MWDF$Neonic[N] <- TreatRefAll$Neonic[which(TreatRefAll$Colony == MWDF$Colony[N] & TreatRefAll$Site == MWDF$Site[N])]
  
  MWDF$Mix[N] <- TreatRefAll$Mixed[which(TreatRefAll$Colony == MWDF$Colony[N] & TreatRefAll$Site == MWDF$Site[N])]
  
  if(!is.na(MWDF$Mix[N])){
    if(MWDF$Mix[N] == 1){
      
      if(MWDF$Day[N] > 60){
        MWDF$MixEmph[N] <- 0
      }else{
        
        MWDF$MixEmph[N] <- SimData$Emphasis[which(SimData$Colony == MWDF$Colony[N] & SimData$Day == MWDF$Day[N] & SimData$State == MWDF$Site[N] & SimData$Task == 'Guard')]
      }
    }
  }
  
  MWDF$UID <- paste0(MWDF$Site, MWDF$Colony)
  MWDF$SD <- paste0(MWDF$Site, MWDF$Day)
  
}

MWDF$MPB <- MWDF$MiteWash/MWDF$Bees

MWDF$PMI <- ceiling((MWDF$MiteWash/MWDF$Bees)*100)

MWMod <- mixed(PMI ~ Neonic*Mix + (1|Site) + (1|Site:UID) + (1|Site:SD),
                family = 'poisson',
                method = 'LRT',
                data = na.exclude(MWDF)
)
nice(MWMod)
summary(MWMod)


## Pollen

DelawarePollen$Site <- 'Delaware'
GeorgiaPollen$Site <- 'Georgia'

GeorgiaPollen$Multiplier <- 1

PollDF <- rbind(DelawarePollen, GeorgiaPollen)

PollDF$Neonic <- NA
PollDF$Mix <- NA
PollDF$UID <- NA

for(N in 1:NROW(PollDF)){
  
  PollDF$Neonic[N] <- TreatRefAll$Neonic[which(TreatRefAll$Colony == PollDF$Colony[N] & TreatRefAll$Site == PollDF$Site[N])]
  
  PollDF$Mix[N] <- TreatRefAll$Mixed[which(TreatRefAll$Colony == PollDF$Colony[N] & TreatRefAll$Site == PollDF$Site[N])]
  
  PollDF$UID <- paste0(PollDF$Site, PollDF$Colony)
  PollDF$SD <- paste0(PollDF$Site, PollDF$Day)
  
}

PollDF$PPD <- PollDF$PollenWeightG/PollDF$Multiplier

PollMod <- mixed(PPD ~ Neonic*Mix + (1|Site) + (1|Site:UID) + (1|Site:SD),
                  method = 'KR',
                  data = na.exclude(PollDF)
)
nice(PollMod)
summary(PollMod)

## Queen Survival

DelawareQueenSurv$Site <- 'Delaware'

QSDF <- rbind(DelawareQueenSurv)

QSDF$Neonic <- NA
QSDF$Mix <- NA
QSDF$UID <- NA

for(N in 1:NROW(QSDF)){
  
  QSDF$Neonic[N] <- TreatRefAll$Neonic[which(TreatRefAll$Colony == QSDF$Colony[N] & TreatRefAll$Site == QSDF$Site[N])]
  
  QSDF$Mix[N] <- TreatRefAll$Mixed[which(TreatRefAll$Colony == QSDF$Colony[N] & TreatRefAll$Site == QSDF$Site[N])]
  
  QSDF$UID <- paste0(QSDF$Site, QSDF$Colony)
  
}

QSMod <- glm(QueenLoss ~ Neonic*Mix,
              family = 'binomial',
              data = QSDF)


library('car')

Anova(QSMod, type = 3, test.statistic = 'Wald')

## Colony Sizes [brood area]

AlabamaSize$Site <- 'Alabama'
DelawareSize$Site <- 'Delaware'
GeorgiaSize$Site <- 'Georgia'

AlabamaSize <- AlabamaSize[,c('Colony','Day','Brood','Bees','Assessor','Site')]

BrADF <- rbind(AlabamaSize, DelawareSize,GeorgiaSize)

BrADF <- BrADF[which(BrADF$Day > 0),]

BrADF$Neonic <- NA
BrADF$Mix <- NA
BrADF$UID <- NA


for(N in 1:NROW(BrADF)){
  
  BrADF$Neonic[N] <- TreatRefAll$Neonic[which(TreatRefAll$Colony == BrADF$Colony[N] & TreatRefAll$Site == BrADF$Site[N])]
  
  BrADF$Mix[N] <- TreatRefAll$Mixed[which(TreatRefAll$Colony == BrADF$Colony[N] & TreatRefAll$Site == BrADF$Site[N])]

  BrADF$UID <- paste0(BrADF$Site, BrADF$Colony)
  BrADF$SD <- paste0(BrADF$Site, BrADF$Day)
  
}

BrAMod <- mixed(Brood ~ Neonic*Mix + (1|Site) + (1|Site:UID) + (1|Site:SD) + (1|Site:Assessor),
                 method = 'KR',
                 data = na.exclude(BrADF)
)
nice(BrAMod)
summary(BrAMod)


#### Let's plot the 'main finding' graphs

CompRes <- function(CUID, Blocking, RV){
  
  TempRVV <- vector(mode='numeric', length = NROW(unique(Blocking$SD[which(Blocking$UID == CUID)])))
  
  for(J in 1:NROW(unique(Blocking$SD[which(Blocking$UID == CUID)]))){
    
    TempBlock <- unique(Blocking$SD[which(Blocking$UID == CUID)])[J]
    
    TempRVV[J] <- (RV[which(Blocking$UID == CUID & Blocking$SD == TempBlock)] - mean(RV[which(Blocking$UID != CUID & Blocking$SD == TempBlock)]))/mean(RV[which(Blocking$SD == TempBlock)])
    
  }
  
  CR <- mean(TempRVV)
  
  return(CR)
  
}

# MDR graph

MDDF2 <- na.exclude(MDDF)

MDG <- data.frame('Colony' = unique(MDDF2$UID))
MDG$Neonic <- NA
MDG$Mix <- NA
MDG$CompositeResidualMDR2 <- NA
MDG$TCode <- NA

for(N in 1:NROW(MDG)){
  
  MDG$Neonic[N] <- unique(MDDF2$Neonic[which(MDDF2$UID == MDG$Colony[N])])
  MDG$Mix[N] <- unique(MDDF2$Mix[which(MDDF2$UID == MDG$Colony[N])])
  
  MDG$CompositeResidualMDR2[N] <- CompRes(CUID = MDG$Colony[N], Blocking = MDDF2[,c('UID', 'SD')], RV = MDDF2$MDR)
  
  MDG$TCode[N] <- paste0('N',MDG$Neonic[N],'M',MDG$Mix[N])
  
}

boxplot(MDG$CompositeResidualMDR2 ~ MDG$TCode)

# MW graph

MWDF2 <- na.exclude(MWDF)

MWG <- data.frame('Colony' = unique(MWDF2$UID))
MWG$Neonic <- NA
MWG$Mix <- NA
MWG$CompositeResidualMWR2 <- NA
MWG$TCode <- NA

for(N in 1:NROW(MWG)){
  
  MWG$Neonic[N] <- unique(MWDF2$Neonic[which(MWDF2$UID == MWG$Colony[N])])
  MWG$Mix[N] <- unique(MWDF2$Mix[which(MWDF2$UID == MWG$Colony[N])])
  
  MWG$CompositeResidualMWR2[N] <- CompRes(CUID = MWG$Colony[N], Blocking = MWDF2[,c('UID', 'SD')], RV = MWDF2$MPB)
  
  MWG$TCode[N] <- paste0('N',MWG$Neonic[N],'M',MWG$Mix[N])
  
}

boxplot(MWG$CompositeResidualMWR2 ~ MWG$TCode)

# BS graph

#BSDF$PropDead <- BSDF$Died/BSDF$Total

BSDF$PropSurv <- BSDF$Survived/BSDF$Total

BSDF2 <- na.exclude(BSDF)

BSG <- data.frame('Colony' = unique(BSDF2$UID))
BSG$Neonic <- NA
BSG$Mix <- NA
BSG$CompositeResidualBSR2 <- NA
BSG$TCode <- NA

for(N in 1:NROW(BSG)){
  
  BSG$Neonic[N] <- unique(BSDF2$Neonic[which(BSDF2$UID == BSG$Colony[N])])
  BSG$Mix[N] <- unique(BSDF2$Mix[which(BSDF2$UID == BSG$Colony[N])])
  
  #BSG$CompositeResidualBSR2[N] <- CompRes(CUID = BSG$Colony[N], Blocking = BSDF2[,c('UID', 'SD')], RV = BSDF2$PropDead)
  
  BSG$CompositeResidualBSR2[N] <- CompRes(CUID = BSG$Colony[N], Blocking = BSDF2[,c('UID', 'SD')], RV = BSDF2$PropSurv)
    
  BSG$TCode[N] <- paste0('N',BSG$Neonic[N],'M',BSG$Mix[N])
  
}

boxplot(BSG$CompositeResidualBSR2 ~ BSG$TCode)

#### Okay now make *nice* graphs

# Quick transparency function for plotting
Transpa <- function(color, percent) {
  
  rgb.val <- col2rgb(color)
  
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100-percent)*255/100)
  
  return(t.col)
  
}

# Set plotting colours
ColRef <- data.frame(Treatment = c('N0M0','N0M1','N1M0','N1M1'), Col =  c('pink3','blue3','red3','purple3'))

# Makes the plots we would like
par(mar=c(5,8,2,2))

## Mite Drop

MDG$TCode <- factor(MDG$TCode, levels = c('N0M0','N1M0','N0M1','N1M1'))

boxplot(MDG$CompositeResidualMDR2 ~ MDG$TCode, 
        main = NA, ylab = 'Difference from Typical Colony - Mite Drop', xlab = 'Treatment', 
        border = 'transparent', 
        cex.axis = 1.3, cex.lab = 1.5, outline = TRUE, lwd = 1.2,
        boxlty = 1, whisklty = 0, staplelty = 1, boxwex = 0.00, 
        col = sapply(X = as.character(ColRef$Col), FUN = Transpa, percent = 100))

abline(h = 0 , lty = 3)

stripchart(MDG$CompositeResidualMDR2 ~ MDG$TCode,
           col = sapply(X = as.character(ColRef$Col[match(ColRef$Treatment, levels(MDG$TCode))]), FUN = Transpa, percent = 40),
           vertical = T, add = T, pch = 4, cex = 1.75, 
           method = 'jitter', lwd = 2.5)

# (this is not real analysis it is for plotting purposes only)
library(emmeans)

PlotMod <- glm(MDG$CompositeResidualMDR2 ~ MDG$TCode,
               family = 'gaussian')

PlotCIs <- as.data.frame(emmeans(PlotMod, specs =c('TCode')))

for(L in 1:NROW(PlotCIs)){
  
  segments(x0 = L, x1 = L, y0 = PlotCIs$asymp.LCL[L], y1 = PlotCIs$asymp.UCL[L],
           col = as.character(ColRef$Col[which(as.character(ColRef$Treatment) == as.character(PlotCIs$TCode[L]))]),
           lwd = 4)
  
}



## Mite Wash

MWG$TCode <- factor(MWG$TCode, levels = c('N0M0','N1M0','N0M1','N1M1'))

boxplot(MWG$CompositeResidualMWR2 ~ MWG$TCode, 
        main = NA, ylab = 'Difference from Typical Colony - Mite Wash', xlab = 'Treatment', 
        border = 'transparent', 
        cex.axis = 1.3, cex.lab = 1.5, outline = TRUE, lwd = 1.2,
        boxlty = 1, whisklty = 0, staplelty = 1, boxwex = 0.00, 
        col = sapply(X = as.character(ColRef$Col), FUN = Transpa, percent = 100))

abline(h = 0 , lty = 3)

stripchart(MWG$CompositeResidualMWR2 ~ MWG$TCode,
           col = sapply(X = as.character(ColRef$Col[match(ColRef$Treatment, levels(MWG$TCode))]), FUN = Transpa, percent = 40),
           vertical = T, add = T, pch = 4, cex = 1.75, 
           method = 'jitter', lwd = 2.5)

# (this is not real analysis it is for plotting purposes only)
library(emmeans)

PlotMod <- glm(MWG$CompositeResidualMWR2 ~ MWG$TCode,
               family = 'gaussian')

PlotCIs <- as.data.frame(emmeans(PlotMod, specs =c('TCode')))

for(L in 1:NROW(PlotCIs)){
  
  segments(x0 = L, x1 = L, y0 = PlotCIs$asymp.LCL[L], y1 = PlotCIs$asymp.UCL[L],
           col = as.character(ColRef$Col[which(as.character(ColRef$Treatment) == as.character(PlotCIs$TCode[L]))]),
           lwd = 4)
  
}


## Brood Death

boxplot(BSG$CompositeResidualBSR2 ~ BSG$TCode, 
        main = NA, ylab = 'Difference from Typical Colony - Brood Survival', xlab = 'Treatment', 
        border = 'transparent', 
        cex.axis = 1.3, cex.lab = 1.5, outline = TRUE, lwd = 1.2,
        boxlty = 1, whisklty = 0, staplelty = 1, boxwex = 0.00, 
        col = sapply(X = as.character(ColRef$Col), FUN = Transpa, percent = 100))

abline(h = 0 , lty = 3)

stripchart(BSG$CompositeResidualBSR2 ~ BSG$TCode,
           col = sapply(X = as.character(ColRef$Col), FUN = Transpa, percent = 40),
           vertical = T, add = T, pch = 4, cex = 1.75, 
           method = 'jitter', lwd = 2.5)

# (this is not real analysis it is for plotting purposes only)
library(emmeans)

PlotMod <- glm(BSG$CompositeResidualBSR2 ~ BSG$TCode,
               family = 'gaussian')

PlotCIs <- as.data.frame(emmeans(PlotMod, specs =c('TCode')))

for(L in 1:NROW(PlotCIs)){
  
  segments(x0 = L, x1 = L, y0 = PlotCIs$asymp.LCL[L], y1 = PlotCIs$asymp.UCL[L],
           col = as.character(ColRef$Col[which(as.character(ColRef$Treatment) == as.character(PlotCIs$TCode[L]))]),
           lwd = 4)
  
}







