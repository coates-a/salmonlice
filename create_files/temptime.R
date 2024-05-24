##### Assign seasonal temperatures on farms ####
##### Calculate daily growth rate on farms ####

library(ggplot2)
library(dplyr)
library(readr)


farms <- read_csv("src/aza_temp/farms.week1.csv") # Farm data
seasons <- read_csv("src/aza_temp/seasons.csv") # Weeks and corresponding 'seasons'
lice <- read_csv("dseed/data/lice2012_2020.csv") # Barentswatch data with temperatures

lice <- lice %>% # Only look at farms that appear in model
  filter(farmID %in% farms$farmID)

lice <- subset(lice, week<53) # Remove few days in week 53

lice <- merge(lice, seasons, by="week") # Assign each week a season

mtemp <- aggregate(temp ~ farmID + season, data=lice, FUN=mean) # Average temp per season
# There are some farms/seasons with NO TEMP DATA that don't appear here

# Create list of ALL farms and season combinations
season.farm <- data.frame(farmID = c(rep(farms$farmID, each=10)), # 10 'seasons' in year
                             season = c(rep(c("w1","w2","w3",
                                               "sp",
                                               "s1","s2","s3","s4","s5",
                                              "w4"), by=i)),
                          lat = c(rep(farms$lat, each=10)))
mtemp.all <- merge(season.farm, mtemp, by=c("farmID","season"), all=T) # Add mean temps per season

## Create new data frame to fill missing temperatures
mtemp2 <- mtemp.all

## Function that calculates the temperature at the nearest farm (by latitude) to a farm with no temp data

fun_nearfarm <- function(mtemp2, ff) { # Where mtemp2 is dataset, and ff is a row of that dataset
  near.farm.season <- mtemp2[!is.na(mtemp2$temp),] # Remove NA temps (don't want to look at farm (including ff farm) with no temp data)
  near.farm.season <- subset(near.farm.season, season == mtemp2[ff,"season"]) # Only looking at temps for the season that matches ff
  near.farm.season$diff <- abs(mtemp2[ff,"lat"] - # Calculate absolute difference between latitude of ff farm
                                 near.farm.season$lat) # and other farm latitudes
  near.farm <- near.farm.season[which.min(near.farm.season$diff),] # Find farm with smallest difference in latitude
  return(near.farm$temp) # Return the temp (during ff season) of that closest farm
  }

# Then fill all NAs in farm X temp X season dataset
for(ff in 1:nrow(mtemp2)){
  ifelse((is.na(mtemp2[ff,"temp"])), # If temp is NA
         
         mtemp2[ff,"temp"] <- fun_nearfarm(mtemp2, ff), # Run above function, input the temp at the nearest farm
         
         mtemp2[ff,"temp"] <- mtemp2[ff,"temp"]) # If not NA, temp is the same
}
 
# Turn into matrix?
mtemp2 <- subset(mtemp2, select = -lat) # Drop latitude column
temp.m <- reshape(mtemp2, idvar = "farmID", timevar = "season", direction = "wide")

write.csv(temp.m, "src/aza_temp/temp.t.csv")

# Need to rename columns

##################################################################

temp.m <- as.matrix(read_csv("src/aza_temp/temp.t.csv"))

#### Convert into development rate ####

# Just temp data (farmId as rownames)
rownames(temp.m)<-temp.m[,1]
temp.m <- temp.m[,-1]

# From Hamre et al. 2019
# Constants (from Hamre, average of male and female)
b <- 0.000581
c <- 0.0094805
d <- 0.0047395
# Calculate daily development rate at each farm & season according to temp 
s.m <- (b*temp.m^2)+(c*temp.m)+d # Eq from Hamre

# Convert from daily transition rate to weekly transition rate
# If daily RETENTION = (1-s), then weekly retention = (1-s)^7
# Therefore weekly TRANSITION = 1-(1-s)^7
s.w <- (1-(1 - s.m)^7)

write.csv(s.w, "src/aza_temp/temp.sw.csv")



###########################################################

### Convert into daily larval production rate

# From Johnsen et al. 2020
# Eggs per string (temp dep)
b0 <- 5.6
b1 <- -0.43
b2 <- -0.78
fT <- exp(b0 + b1 * (temp.m/10) + b2 * (log(temp.m/10))^2)

# Days between hatches
a0 <- 4.85e-4
a1 <- 8.667e-3
a2 <- 3.75e-3
DE = 0.25*(5/(a0 * temp.m^2 + a1 * temp.m + a2))

# Daily number of larvae per Adult
f <- fT/DE
# Weekly number of larvae
f.w <- f * 7

write.csv(f.w, "src/aza_temp/temp.fw.csv")

##### PLOTS #####
temp.t <- temp.t[,-1]
temp <- melt(temp.t)
colnames(temp) <- c("season","temp")

temp <- merge(temp, seasons, by="season", all=T)

temp.plot <- ggplot()+
  geom_point(data=lice, aes(x=week, y=temp),
             size=0.1, col="grey")+
  geom_point(data=temp, aes(x=week, y=temp),
             size=0.9, col="red")+
  scale_x_continuous(name="Week")+
  scale_y_continuous(name=expression(Temperature~(degree~C)))+
  theme(axis.text.x = element_text(size=14), #text(angle = 45, size = 28, vjust=1, hjust=1),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x = element_text(size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        plot.title = element_text(size=18, margin = margin(t = 0, r = 0, b = 5, l = 0)),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.key= element_rect(fill = NA),
        legend.text = element_text(size = 11),
        axis.line = element_line(colour = "gray50", size = 0.6, linetype = "solid"),
        panel.border = element_rect(colour = "black", fill=NA,size=1.5))

ggsave(file = sprintf("outputs/aza_temp/temp.supp.png"), 
       temp.plot, # Save as png
       width = 10, height = 5, dpi = 700, units = "in", device='png')


##### DEV & REPROD @ CONSTANT TEMP (6 DEGREES)

temp.m <- as.matrix(read_csv("src/aza_temp/temp.t.csv"))

#### Convert into development rate ####

# Just temp data (farmId as rownames)
rownames(temp.m)<-temp.m[,1]
temp.m <- temp.m[,-1]

# Constant temp of 6 degrees throughout
temp.c <- temp.m
temp.c[]<- 6

# From Hamre et al. 2019
# Constants (from Hamre, average of male and female)
b <- 0.000581
c <- 0.0094805
d <- 0.0047395
# Calculate daily development rate at each farm & season according to temp 
s.c <- (b*temp.c^2)+(c*temp.c)+d # Eq from Hamre

# Convert from daily transition rate to weekly transition rate
# If daily RETENTION = (1-s), then weekly retention = (1-s)^7
# Therefore weekly TRANSITION = 1-(1-s)^7
s.c <- (1-(1 - s.c)^7)

write.csv(s.c, "src/aza_temp/temp6.sw.csv")



###########################################################

### Convert into daily larval production rate

# From Johnsen et al. 2020
# Eggs per string (temp dep)
b0 <- 5.6
b1 <- -0.43
b2 <- -0.78
fT.c <- exp(b0 + b1 * (temp.c/10) + b2 * (log(temp.c/10))^2)

# Days between hatches
a0 <- 4.85e-4
a1 <- 8.667e-3
a2 <- 3.75e-3
DE.c = 0.25*(5/(a0 * temp.c^2 + a1 * temp.c + a2))

# Daily number of larvae per Adult
f.c <- fT.c/DE.c
# Weekly number of larvae
f.c <- f.c * 7

write.csv(f.c, "src/aza_temp/temp6.fw.csv")
