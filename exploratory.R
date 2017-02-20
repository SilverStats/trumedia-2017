########################################################
# Code for the TruMedia MLB Hackathon                  #
# 20 February 2017                                     #
# Steven Silverman                                     #
#                                                      #
# See the associated .Rmd file for all graphics code   #
# Contact: Steven@SilverStats.com                      #
########################################################

setwd("~/Analytics articles/TruMedia/trumedia-2017")

library(plyr)
library(dplyr)
library(scoring)
library(mgcv)

# first attempt, using the World Series data to play around
sample <- read.csv("2016-WS.csv")

# 0/1 vector of whether the batter swung at the pitch
sample$swing <- sample$pitchResult %in% c("SS", "F", "FT", "FB", "MB",
                                          "IP", "CI")

sample$brier <- (sample$swing - sample$probCalledStrike)^2


grouped <- (group_by(sample, batter) %>% summarize(count = n(), avgBrier = mean(brier))
                                     %>% arrange(avgBrier))

scores <- brierscore(swing ~ probCalledStrike, data=sample)


# actual code
data.2014 <- read.csv("2014.csv")
data.2015 <- read.csv("2015.csv")
train <- rbind(data.2014, data.2015)

data.2016 <- read.csv("2016.csv")

# 0/1 vector of whether the batter swung at the pitch
train$swing <- train$pitchResult %in% c("SS", "F", "FT", "FB", "MB",
                                        "IP", "CI")
train$brier <- (train$swing - train$probCalledStrike)^2


pitch_cutoff <- 1500 # on the order of 400 plate appearances, which will be good enough

# create some grouped data frames to get a feel for the Brier scores
groupedBatter <- (group_by(train, batter) %>% summarize(count = n(), avgBrier = mean(brier, na.rm=T))
                    %>% filter(count >= pitch_cutoff) %>% arrange(avgBrier))
groupedPitcher <- (group_by(train, pitcher) %>% summarize(count = n(), avgBrier = mean(brier, na.rm=T))
                   %>% filter(count >= pitch_cutoff) %>% arrange(desc(avgBrier)))




# processing code for getting a few aggregate measures
# some of which I used in last year's hackathon as well

# pitch groupings--used Brooks Baseball/PitchInfo definitions
fastballTypes <- c("FA", "FT", "FF", "FC", "SI")
breakingTypes <- c("CU", "SL", "KC", "GY")
offspeedTypes <- c("CH", "FS", "SC", "EP", "FO")
otherTypes <- c("KN", "PO", "IN", "AB", "AS", "UN")

veloMean <- function(velo, type, group) {  # again, using Brooks definitions
    if (group == 0) {
        bool <- type %in% fastballTypes
    }
    if (group == 1) {
        bool <- type %in% breakingTypes
    }
    if (group == 2) {
        bool <- type %in% offspeedTypes
    }
    if (group == 3) {
        return(NA)
    }
    velos <- velo[bool]
    return(mean(velos, na.rm = TRUE))
}

herf <- function(pitchType) { # Herfindahl index for pitch arsenal
    types <- table(pitchType)
    return(sum((types/sum(types))^2, na.rm = TRUE))
}

# movement: adds columns with x and z movement (in feet) to a data frame
addMovement <- function(dataframe) {
    
    # use kinematic equations
    
    # calculate time ball takes to travel from 50ft to home plate
    t <- ((-dataframe$vy0 - sqrt(dataframe$vy0^2+2*dataframe$ay*(-50)))
          / dataframe$ay)
    
    # calculate actual vs. predicted movement; find difference
    xd <- dataframe$vx0*t
    x.move <- (dataframe$px-dataframe$x0) - xd
    zd <- dataframe$vz0*t - 16*t^2 # subtract out gravity
    z.move <- (dataframe$pz - dataframe$z0) - zd
    
    dataframe$mx <- x.move
    dataframe$mz <- z.move
    
    return(dataframe)
    
}

spin <- function(spinRate, type, group) { # group = 0 for FB, 1 for breaking, 2 for offspeed
    if (group == 0) {
        bool <- type %in% fastballTypes
    }
    if (group == 1) {
        bool <- type %in% breakingTypes
    }
    if (group == 2) {
        bool <- type %in% offspeedTypes
    }
    if (group == 3) { # other
        return(NA)
    }
    return(mean(spinRate[bool], na.rm = TRUE))
}

move <- function(movement.pitch, type, group) { # group = 0 for FB, 1 for breaking, 2 for offspeed
    if (group == 0) {
        bool <- type %in% fastballTypes
    }
    if (group == 1) {
        bool <- type %in% breakingTypes
    }
    if (group == 2) {
        bool <- type %in% offspeedTypes
    }
    if (group == 3) {
        return(NA)
    }
    return(mean(abs(movement.pitch[bool]), na.rm = TRUE))
}

train <- addMovement(train) # add the movement columns


totalsPitcher <- (group_by(train, pitcher) %>% 
                  summarize(count = n(),
                            fbVelo = veloMean(releaseVelocity, pitchType, 0),
                            brVelo = veloMean(releaseVelocity, pitchType, 1),
                            offVelo = veloMean(releaseVelocity, pitchType, 2),
                            rx = mean(x0, na.rm = TRUE), # release x (left-right)
                            rz = mean(z0, na.rm = TRUE), # release z (up-down)
                            sdX = var(x0, na.rm = TRUE), # release standard deviation
                            sdZ = var(z0, na.rm = TRUE),
                            herf = herf(pitchType),
                            fbSpin = spin(spinRate, pitchType, 0),
                            breakSpin = spin(spinRate, pitchType, 1),
                            offSpin = spin(spinRate, pitchType, 2),
                            fbx = move(mx, pitchType, 0), # horizontal movement
                            fbz = move(mz, pitchType, 0), # vertical movement
                            brx = move(mx, pitchType, 1),
                            brz = move(mz, pitchType, 1),
                            offx = move(mx, pitchType, 2),
                            offz = move(mz, pitchType, 2),
                            avgBrier = mean(brier, na.rm=T)) %>% 
                      filter(count >= pitch_cutoff) %>% arrange(avgBrier))


# fit lots of GAMs, trying a bunch of variable combinations; pick the best one with CV
# I initially wrote code like this for a class; this is a modified version

formula.response <- "avgBrier ~"

always.gam <- c("s(herf)","s(fbVelo)")
sometimes.gam <- c("s(fbx)","s(fbz)","s(brx)","s(brz)","s(fbSpin)","s(breakSpin)",
                   "s(fbVelo,offVelo)", "s(fbx,brx,offx)", "s(fbz,brz,offz)")

formula.always.gam <- paste(always.gam, collapse = " + ")

n <- length(sometimes.gam)
l <- rep(list(0:1), n)

combinations <- expand.grid(l)

models.gam <- list()
formulas.gam <- as.formula(paste(formula.response, formula.always.gam))

models.gam[[1]] <- gam(formulas.gam, data=totalsPitcher)

for (index in 2:nrow(combinations)) { # ignore first row and do manually
    print(index)
    formula.temporary <- paste(sometimes.gam[as.logical(combinations[index,])],
                               collapse = " + ")
    formula.final <- as.formula(paste(formula.response, formula.always.gam, " + ",
                                      formula.temporary))
    formulas.gam <- c(formulas.gam, formula.final)
    models.gam[[index]] <- gam(formula.final, data=totalsPitcher)
}

CV.gam <- sapply(models.gam, function(model) { model$gcv.ubre })
best.gam <- models.gam[[which.min(CV.gam)]]
best.formula.gam <- formulas.gam[[which.min(CV.gam)]]

# manually type out the formula so I don't have to examine the object every time
#best.formula.gam2 <- "avgBrier ~ s(herf) + s(fbVelo) + + s(fbx) + s(fbz) + 
#                                 s(brx) + s(fbSpin) + s(breakSpin)  + s(fbVelo, offVelo) + 
#                                 s(fbx, brx, offx) + s(fbz, brz, offz)"
#best.formula.gam2 <- as.formula(best.formula.gam2)

plot(best.gam, scale=0, se=2, shade=TRUE, select=1)
plot(best.gam, scale=0, se=2, shade=TRUE, select=2)
plot(best.gam, scale=0, se=2, shade=TRUE, select=3)
plot(best.gam, scale=0, se=2, shade=TRUE, select=4)
plot(best.gam, scale=0, se=2, shade=TRUE, select=5)
plot(best.gam, scale=0, se=2, shade=TRUE, select=6)
plot(best.gam, scale=0, se=2, shade=TRUE, select=7)

plot(best.gam, scale=0, se=2, shade=TRUE, select=8, cex.lab=2, cex.main=2,
     cex.axis=2.5)

# Predictions on test data set (2016; training was 2014-15)

pitch_cutoff_test <- 750 # lower cutoff for one year

test <- addMovement(data.2016)

test$swing <- test$pitchResult %in% c("SS", "F", "FT", "FB", "MB",
                                        "IP", "CI")
test$brier <- (test$swing - test$probCalledStrike)^2

testPitcher <- (group_by(test, pitcher) %>% 
                      summarize(count = n(),
                                fbVelo = veloMean(releaseVelocity, pitchType, 0),
                                brVelo = veloMean(releaseVelocity, pitchType, 1),
                                offVelo = veloMean(releaseVelocity, pitchType, 2),
                                rx = mean(x0, na.rm = TRUE), # release x (left-right)
                                rz = mean(z0, na.rm = TRUE), # release z (up-down)
                                sdX = var(x0, na.rm = TRUE), # release standard deviation
                                sdZ = var(z0, na.rm = TRUE),
                                herf = herf(pitchType),
                                fbSpin = spin(spinRate, pitchType, 0),
                                breakSpin = spin(spinRate, pitchType, 1),
                                offSpin = spin(spinRate, pitchType, 2),
                                fbx = move(mx, pitchType, 0),
                                fbz = move(mz, pitchType, 0),
                                brx = move(mx, pitchType, 1),
                                brz = move(mz, pitchType, 1),
                                offx = move(mx, pitchType, 2),
                                offz = move(mz, pitchType, 2),
                                avgBrier = mean(brier, na.rm=T)) %>% 
                      filter(count >= pitch_cutoff_test) %>% arrange(avgBrier))

predicted <- predict(best.gam, newdata=testPitcher)
mae <- mean(abs(predicted-testPitcher$avgBrier), na.rm=T)


# Save off relevant objects to avoid re-reading CSV files and redoing calculations

save(best.gam.formula, "gamFormula.RData")
save(best.gam, "gam.RData")
save(totalsPitcher, "totalsPitcher.RData")
save(groupedPitcher, "groupedPitcher.RData")
save(groupedBatter, "groupedBatter.RData")
save(testPitcher, "testPitcher.RData")