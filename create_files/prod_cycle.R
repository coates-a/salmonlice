#### CREATE MATRIX ASSIGNING EACH FARM AS 'ACTIVE' OR 'INACTIVE' (FALLOW) OVER A 2-YEAR PRODUCTION CYCLE ####

library(readr)
library(dplyr)
library(reshape)

lice <- read_csv("src/2loci/create_files/lice_2012_2022.csv") # Complete farm infection data for 2012-2022, from Barents Watch 
farms <- read_csv("src/2loci/farms.week1.csv") # Farm parameter data set used in simulation 

lice <- lice %>%
  filter(farmID %in% farms$farmID) # Only look at farms used in the simulation 

lice$yr <- lice$year + (lice$week - 1)/52 # Turn Week into continuous variable (proportion of year)

lice.prod <- lice[c("farmID","yr","no.fish")] # Subset farm ID, week, and whether farm has 'no fish' (estimated from Barents Watch)
lice.prod <- arrange(lice.prod, yr) # Reorder by yr and ID, to match simulation layout
lice.prod <- arrange(lice.prod, farmID)
lice.prod <- as.data.frame(lice.prod) # Save as dataset

### There are some cases where farms have just 1 or a few weeks of 'no fish' in the middle of production (or occassionally vice versa)
# Most likely a mistake with Barents Watch (e.g. farms just don't report for a few weeks)
# (although it could be due to farms moving fish around, etc)

#################################################################
#### IGNORE ISOLATED PERIODS (only 1-3 weeks) OF AN ACTIVITY ####

# If activity is different to previous week AND to 3 weeks later
for(nn in 2:nrow(lice.prod)){
  
  ifelse(lice.prod[nn,3] != lice.prod[nn-1,3] & # If activity in a week if different to week before 
      lice.prod[nn,3] != lice.prod[nn + 3,3] & # AND 3 WEEKS AFTER
           lice.prod[nn,1] == lice.prod[nn-1,1], # And it's the same farm as previous row
         
         lice.prod[nn,3] <- lice.prod[nn-1,3], # Change activity to be the same as the previous week
         
         lice.prod[nn,3] <- lice.prod[nn,3]) # Otherwise leave the same
}
# If activity is different to previous week AND to 3 weeks later
for(nn in 2:nrow(lice.prod)){
  
  ifelse(lice.prod[nn,3] != lice.prod[nn-1,3] & # If activity in a week if different to week before 
           lice.prod[nn,3] != lice.prod[nn + 2,3] & # AND 2 WEEKS AFTER
           lice.prod[nn,1] == lice.prod[nn-1,1], # And it's the same farm as previous row
         
         lice.prod[nn,3] <- lice.prod[nn-1,3], # Change activity to be the same as the previous week
         
         lice.prod[nn,3] <- lice.prod[nn,3]) # Otherwise leave the same
}
# If activity is different to previous week AND to following week
for(nn in 2:nrow(lice.prod)){
  
  ifelse(lice.prod[nn,3] != lice.prod[nn-1,3] & # If activity in a week if different to week before 
           lice.prod[nn,3] != lice.prod[nn + 1,3] & # AND 1 WEEKS AFTER
           lice.prod[nn,1] == lice.prod[nn-1,1], # And it's the same farm as previous row
         
         lice.prod[nn,3] <- lice.prod[nn-1,3], # Change activity to be the same as the previous week
         
         lice.prod[nn,3] <- lice.prod[nn,3]) # Otherwise leave the same
}


###############################################################
#### Group consecutive weeks of active or inactive together ####

lice.prod$grp <- 1

for(nn in 2:nrow(lice.prod)){
  
  ifelse(lice.prod[nn,3] == lice.prod[nn-1,3] & # If no.fish is the same as the previous week
     lice.prod[nn,1] == lice.prod[nn-1,1], # And it's the same farm
     
     lice.prod[nn,4] <- lice.prod[nn-1,4], # Make the group for that week the same as for the previous week
     
     lice.prod[nn,4] <- lice.prod[nn-1,4]+1) # Assign it to the next group, if a new activity or a new farm
}

lice.prod$wk <- 1
lice.ag <- aggregate(wk ~ farmID + no.fish + grp, data=lice.prod, FUN=sum) # Add up the number of weeks in each group
lice.ag$activity <- ifelse(lice.ag$no.fish=="Ja", # Assign activity names as 'active' or 'nofish'
                           "nofish",
                           "active")



########################################################################
#### DETERMINE WHAT TIME OF YEAR THE PRODUCTION CYCLE STARTS AND ENDS ####

lice.min <- aggregate(yr ~ farmID + no.fish + grp, data=lice.prod, FUN=min) # Subset the earliest instance of each group
colnames(lice.min)[4] <- "start" 
lice.ag <- merge(lice.ag, lice.min[c(1,3,4)], by=c("farmID","grp")) # Add start time for each group
# When does cycle end? 
lice.max <- aggregate(yr ~ farmID + no.fish + grp, data=lice.prod, FUN=max) # Subset the latest instance of each group
colnames(lice.max)[4] <- "end" 
lice.ag <- merge(lice.ag, lice.max[c(1,3,4)], by=c("farmID","grp")) # Add end time for each group


#####################################################################
#### REMOVE THE FIRST AND LAST GROUP FOR EACH FARM ####
# On Jan 1st 2012, farms are already in the middle of a cycle (starting pre-2012) OR have not have been established yet
# On Dec 31st 2022, farms are mid-cycle, OR have been permanently closed

min.start <- data.frame(farmID = unique(lice.ag$farmID), # List of the earliest start dates for each farm
                  min.start = tapply(lice.ag$start, lice.ag$farmID, min)) # Most begin at 2012, but some farms are NAs until later
lice.ag <- merge(lice.ag, min.start, by="farmID") # Add earliest start date to data set
lice.ag <- subset(lice.ag, lice.ag$start != lice.ag$min.start) # Remove groups where the start of a cycle IS THE SAME AS the earliest (first) start date

max.end <- data.frame(farmID = unique(lice.ag$farmID), # List of the latest end dates for each farm
                  max.end = tapply(lice.ag$end, lice.ag$farmID, max)) # Most end at 2022.9, but some farms are NAs by then
lice.ag <- merge(lice.ag, max.end, by="farmID") 
lice.ag <- subset(lice.ag, lice.ag$end != lice.ag$max.end) # Remove groups where the end of a cycle IS THE SAME AS the latest (last) end date



#############################################################################
#### CALCULATE THE MEAN NUMBER OF WEEKS ACTIVE AND INACTIVE FOR EACH FARM ####

lice.ag.m <- aggregate(wk ~ farmID + activity, data=lice.ag, FUN=mean)


############################################################
#### CALCULATE TIME OF THE VERY FIRST STOCKING OF FISH #####

stock1 <- subset(lice.ag, activity=="active") # Subset when farms are active
stock1 <- data.frame(farmID = unique(stock1$farmID), 
                         grp = tapply(stock1$grp, stock1$farmID, min)) # Subset by the earliest active group for each farm
stock1 <- merge(stock1, lice.ag[c("grp","farmID","start","end","min.start")])

stock1$start <- stock1$start - stock1$min.start # At what time of the year is each farm first stocked
stock1$end <- stock1$end - stock1$min.start # At what time does the stocking period end 

# AT THE BEGINNING OF THE SIMULATED PRODUCTION CYCLE, ALL FARMS ARE FIRST STOCKED WITHIN 2 YEARS OF EACH OTHER
stock1$start1 <- ifelse(stock1$start<2, # If the start year is <2 (i.e. first stocked 2012-2013)
                        stock1$start, # Keep values same.
                        # If the year is >2 (first stocked after 2013) AND
                        ifelse(abs(floor(stock1$start)/2) %% 1 == 0, # If the year (rounded down) is even (no decimal point),
                               abs(stock1$start) %% 1, # Then just keep decimal point (i.e. '0.x')
                               abs(stock1$start) %% 1 + 1)) # If odd year, make '1.x'

lice.ag.m <- merge(lice.ag.m, stock1[c("farmID","start1")], by="farmID")
lice.ag.m$start.ts <- lice.ag.m$start1 * 52 # convert decimal years to weekly time-step
lice.ag.m$wk <- round(lice.ag.m$wk) # round the lengths of cycles to the nearest whole week

#####################################
#### MAKE SAME STRUCTURE AS N  #####

# Reshape data set to match the N matrix in the Simulation
lice.ag.m <- reshape(lice.ag.m, v.names="wk", # Columns = time-steps
              timevar="activity",idvar="farmID", direction="wide") # Rows = farms

farms.ag <- merge(farms["farmID"], lice.ag.m, by="farmID",all=T) # Add missing farm sites not included in Barents Watch data

farms.ag <- as.matrix(farms.ag) # Convert to matrix
farms.ag[,"start.ts"] <- farms.ag[,"start.ts"] + 1 # add 1 so ts=0 is 1st column of matrix


#########################################
#### CREATE PRODUCTION CYCLE MATRIX #####

prodcy <- matrix(0, nrow=nrow(farms.ag), # 1 row for each farm
                 ncol=52*150) # Make matrix longer much longer than needed -> will subset timeframe in 'parameters' script
row.names(prodcy) <- c(farms.ag[,"farmID"])


farms.ag[,"start.ts"] <- ifelse(is.na(farms.ag[,"start.ts"]), # Can't have NA in start.ts 
             1, # Just make 1 for NAs
             farms.ag[,"start.ts"])

write.csv(farms.ag, 'src/2loci/create_files/prodcy.summary.csv')

#########################################################
#### FOR EACH FARM & TIME-STEP, ASSIGN FARM ACTIVITY ####

# '1' for active, '0' for no fish

for(rr in 1:nrow(prodcy)){
  
  ifelse(is.na(farms.ag[rr,"wk.nofish"]), # If no farm data (NA) --> assume broodstock/juveniles --> constant fish
         cycs <- rep(1, # Just make a string of 1s
                     times=104), 
  
  cycs <- rep(c # Cycling for 1 farm
              (rep(1, times=c(farms.ag[rr,"wk.active"])), # 1 for weeks active
        rep(0, times=c(farms.ag[rr,"wk.nofish"]))), # 0 for no fish
      times=20)) # repeat production cycle a number of times
  ### Repeat more times than needed 
  
  #
  prodcy[rr, farms.ag[rr,"start.ts"]: # Insert vector of cycles starting from initial stocking time
           c(farms.ag[rr,"start.ts"] + length(cycs) - 1)] <- cycs # For length of cycles
  
}
  
#### SAVE ####

write.csv(prodcy, file="src/2loci/prod.cycle.csv")




  



###################################
#### Farms not in Barentswatch ####
length(unique(lice$farmID)) # I'm using 10 farms that don't appear in barentswatch dataset...

farms.unique <- farms %>%
  filter(!farmID %in% lice$farmID) # Filter the 10 farms I'm using NOT found in Barentswatch

View(farms.unique)

ggplot()+
  geom_point(data=farms, aes(y=Dwcon, x=biomass), col="grey")+
  geom_point(data=farms.unique, aes(y=Dwcon, x=biomass), col="red")
# These 10 farms all have pretty low biomass (325t). 
# Pretty standard looking spread for connectivity. (None of these farms have the very high conn values)

ggplot()+
  geom_point(data=farms, aes(x=lon, y=lat), col="grey")+
  geom_point(data=farms.unique, aes(x=lon, y=lat), col="red")
# 3 farms in Zone2, 4 farms in Zone6

#################################################################

### Plot all farms - active vs inactive over time ####
active <- subset(lice, no.fish=="Nei")

active.plot <- ggplot(active, aes(x=yr, y= as.factor(farmID)))+
  geom_point(size=1.5, col="red")+
  geom_vline(xintercept=c(2012:2022), lty=2)

ggsave(file = "src/2loci/farm_prod/active.plot.pdf", 
       active.plot, # Save as png
       width = 30, height = 40, dpi = 1000, units = "in", device='pdf')

################
