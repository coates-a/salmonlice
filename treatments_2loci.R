##### LOG EACH TIME & PLACE DELOUSING TREATMENT IS APPLIED ####
treat.log <- matrix(0, nrow=i, ncol=ts) # Build matrix to record treatments 
rownames(treat.log) <- farms$farmID # Row for each farm
colnames(treat.log) <- c(0:(ts-1)) # Column for each time-step

####################################
#### Run treatment log function ####

treatments <- treat.log_fun(N, l, V, X, Y, Zdisc, MAX, treat.log, prod.cycle)

####################################

write.csv(treatments, file = sprintf("outputs/treatments.all_%s_%s.csv", scenarios, sim), row.names = F) # Save output

#### Summarise treatment frequency through TIME ####

xtreats <- ifelse(treatments == "X", 
                  1,
                  0)
treat.x <- data.frame(treat.x = colSums(xtreats), # Count number of X treatments (all farms) per time-step
                      t = colnames(treatments))
treat.x$yr <- rep(c(0:(years-1)), each = 52) # Add year
tts <- c(seq(0,51,1)/52) # Fraction of the year represented by time-step
treat.x$tt <- rep(c(tts), times=years)
treat.x$tt <- treat.x$tt + treat.x$yr # Add continuous measure of time
treat.x <- treat.x[,-c(2,3)] # Drop t and yr columns (too cluttered when merging dataframes)
# Count number of discrete Y treatments
ytreats <- ifelse(treatments == "Y", 
                  1,
                  0)
treat.y <- data.frame(treat.y = colSums(ytreats), # Count number of Y treatments (all farms) per time-step
                      t = colnames(treatments))
treat.y$yr <- rep(c(0:(years-1)), each = 52) # Add year
tts <- c(seq(0,51,1)/52) # Fraction of the year represented by time-step
treat.y$tt <- rep(c(tts), times=years)
treat.y$tt <- treat.y$tt + treat.y$yr # Add continuous measure of time
treat.y <- treat.y[,-c(2,3)] # Drop t and yr columns (too cluttered when merging dataframes)

treat.t <- merge(treat.x, treat.y, by="tt") # Merge 2 treatment datasets

# Count number of discrete Z treatments
ztreats <- ifelse(treatments == "Z", 
                  1,
                  0)
treat.z <- data.frame(treat.z = colSums(ztreats), # Count number of Y treatments (all farms) per time-step
                      t = colnames(treatments))
treat.z$yr <- rep(c(0:(years-1)), each = 52) # Add year
tts <- c(seq(0,51,1)/52) # Fraction of the year represented by time-step
treat.z$tt <- rep(c(tts), times=years)
treat.z$tt <- treat.z$tt + treat.z$yr # Add continuous measure of time
treat.z <- treat.z[,-c(2,3)] # Drop t and yr columns (too cluttered when merging dataframes)

treat.t <- merge(treat.t, treat.z, by="tt") # Merge 2 treatment datasets

# Treats per year BY ZONE

treats.zone <- matrix(nrow=nrow(ztreats), # Create matrix to fill
                 ncol=years+1)
colnames(treats.zone) <- c(c(1:years),"zone") # Column for each year, plus one with Zone ID
treats.zone[,years+1] <- c(farms$zone)

for(nn in 1:years){
  
  treats.zone[,nn] <- c(rowSums(ztreats[,c((0+(52*(nn-1))): # For each year (block of 52 weeks), sum number of treatments per farm
                                        (52+(52*(nn-1))))]))
  
}

treats.zone <- aggregate(. ~ zone, data=as.data.frame(treats.zone), FUN=sum) # Sum all farms in each zone

treats.zone <- melt(treats.zone, # Melt data into single column
               id.vars = "zone",
               variable.name = "year",
               value.name = "freq")
colnames(treats.zone) <- c('zone','year',sim)

# write.csv(treats.zone, file = sprintf("outputs/treat.zone/treat.z_%s_%s.csv", scenarios, sim), row.names = F)

##############

# Treats PER FARM PER YEAR

treats.farm <- matrix(nrow=nrow(ztreats), # Create matrix to fill
                      ncol=years+1)
colnames(treats.farm) <- c(c(1:years),"farm") # Column for each year, plus one with farm ID
treats.farm[,years+1] <- c(farms$farm)

for(nn in 1:years){
  
  treats.farm[,nn] <- c(rowSums(ztreats[,c((0+(52*(nn-1))): # For each year (block of 52 weeks), sum number of treatments per farm
                                             (52+(52*(nn-1))))]))
  
}

treats.farm <- as.data.frame(treats.farm)
treats.farm <- melt(treats.farm, # Melt data into single column
                    id.vars = "farm",
                    variable.name = "year",
                    value.name = "freq")
colnames(treats.farm) <- c('farm','year',sim)

write.csv(treats.farm, file = sprintf("outputs/treat.f.y/treat.f.y_%s_%s.csv", scenarios, sim), row.names = F)



################
# Count harvests
harvs <- ifelse(treatments == "H", 
                1,
                0)
treat.h <- data.frame(treat.h = colSums(harvs), # Count number of X treatments (all farms) per time-step
                      t = colnames(treatments))
treat.h$yr <- rep(c(0:(years-1)), each = 52) # Add year
tts <- c(seq(0,51,1)/52) # Fraction of the year represented by time-step
treat.h$tt <- rep(c(tts), times=years)
treat.h$tt <- treat.h$tt + treat.h$yr # Add continuous measure of time
treat.h <- treat.h[,-c(2,3)] # 

treat.t <- merge(treat.t, treat.h, by="tt") # Merge with harvest

treat.t <- melt(treat.t, # Melt data into single column
                id.vars = "tt")
colnames(treat.t)[c(2,3)] <- c("treat","freq")

tr.t <- treat.t # Create new dataframe to save (mapPlot.R will still use treat.t)
colnames(tr.t)[3] <- sim # Assign name of simulation
write.csv(tr.t, file = sprintf("outputs/treat.t/treat.t_%s_%s.csv", scenarios, sim), row.names = F) # Save to treat.t folder to compare simulations


#### Number of treatments frequency by FARM (over whole simulation) ####

ntreats <- ifelse(treatments == "X" | treatments == "Y" | treatments == "Z", 
                1,
                0)
treat.f <- data.frame(treats = rowSums(ntreats), # Count treatments per farm
                      farmID = rownames(treatments))
colnames(treat.f)[1] <- sim
# write.csv(treat.f, file = sprintf("outputs/treat.f/treat.f_%s_%s.csv", scenarios, sim), row.names = F) # Save to treat.f folder (to compare sims)

#### Number of harvests by FARM (over whole simulation) ####

nharv <- ifelse(treatments == "H", 
                  1,
                  0)
harv.f <- data.frame(harv = rowSums(nharv),
                      farmID = rownames(treatments))
colnames(harv.f)[1] <- sim
# write.csv(harv.f, file = sprintf("outputs/harv.f/harv.f_%s_%s.csv", scenarios, sim), row.names = F) # Save to harv.f folder (to compare sims)

# Treats PER FARM PER YEAR

harv.farm <- matrix(nrow=nrow(nharv), # Create matrix to fill
                      ncol=years+1)
colnames(harv.farm) <- c(c(1:years),"farm") # Column for each year, plus one with farm ID
harv.farm[,years+1] <- c(farms$farm)

for(nn in 1:years){
  
  harv.farm[,nn] <- c(rowSums(nharv[,c((0+(52*(nn-1))): # For each year (block of 52 weeks), sum number of treatments per farm
                                             (52+(52*(nn-1))))]))
  
}

harv.farm <- as.data.frame(harv.farm)
harv.farm <- melt(harv.farm, # Melt data into single column
                    id.vars = "farm",
                    variable.name = "year",
                    value.name = "harv")
colnames(harv.farm) <- c('farm','year',sim)

write.csv(harv.farm, file = sprintf("outputs/treat.f.y/harv.f.y_%s_%s.csv", scenarios, sim), row.names = F)


