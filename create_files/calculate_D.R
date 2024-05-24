
#### Calculate Dispersal Matrix (D) ####

library(readr)
install.packages("reshape2")
library(reshape2)

farms <- read_csv("src/2loci/farms.week1.csv") 

##### FOR FUTURE: #####
## Write function that goes through all the datasets (dispersal at different dates) and creates D matrix for each


dispersal <- read.csv('dseed/data/Dispersal sims/AprJul2014.csv', header = TRUE) # Different dispersal data for different times! Try a few...
dispersal <- dispersal[!(dispersal$source==10173),] # Remove the connections from the node that is not a sea farm. Node name = 10173
dispersal <- dispersal[!(dispersal$target==10173),]

total_part <- 122 * 24 # Total number of particles released per farm (hourly over 120 days in WINTER, 122 in SPRING))
dispersal$prob <- ((dispersal$weight/total_part)) # proportion of lice from farm, to farm 

## Need to include ALL FARMS (from farms data), 
# including those that don't apear in this D (since they might appear in the other D)

dispersal <- rbind(dispersal, # Add new data.frame (with missing farms) to dispersal
                   data.frame(source = setdiff(farms$farmID, dispersal$source), # Find all IDs in farms that DON'T appear in dispersal
                target = setdiff(farms$farmID, dispersal$source), # Build into data frame that matches dispersal layout
                weight=0,
                prob=0))
rownames(dispersal)<-c(1:nrow(dispersal)) # fix row numbers

# Need dcast to include all combinations of source & target (even if they don't appear in dispersal)
# Include ALL farm IDs as levelled factors
dispersal$source <-factor(dispersal$source)
dispersal$target <-factor(dispersal$target)
IDs <- unique(c(levels(dispersal$source),levels(dispersal$target)))
dispersal$source <-factor(dispersal$source, levels=IDs)
dispersal$target <-factor(dispersal$target, levels=IDs)

# Turn into matrix (or matrix-style data frame)

D <- dcast(dispersal, source~target, value.var="prob", fill=0, drop=FALSE)
D <- D[,-1] # Remove header column
D <- as.matrix(D) # Save as matrix

Ds14 <- t(D) # Transpose D so cols = source, rows = target

# Mean D values (each cell of matrix) for all years
Ds <- (Ds14 + Ds13 + Ds12 + Ds11 + Ds10 +Ds09)/6

# To save
write.csv(Ds, file = 'dseed/data/Ds.csv', row.names = F)



############################################################################################

##### QUANTIFY CONNECTIVITY OF FARMS ######

Dw <- read_csv("src/2loci/Dw.csv")
Dw <- as.matrix(Dw)
Ds <- read_csv("src/2loci/Ds.csv")
Ds <- as.matrix(Ds)

# Indicator of INCOMING connection strength:
# Sum of all D values (excluding self-recruitment) for lice going TO a farm
DwX <- Dw
diag(DwX) <- 0 # exclude self-recruitment
Dw.in <- rowSums(DwX) # Sum all the dispersal probabilities TO farm

DsX <- Ds
diag(DsX) <- 0 # exclude self-recruitment
Ds.in <- rowSums(DsX)

# Add to farms dataset
farms$Dw.in <- Dw.in
farms$Ds.in <- Ds.in

# The NUMBER of other farms sending ANY lice to a farm
Dw.f.in <- rowSums(DwX!=0) # Count number of farms with D>0 (excluding self-recruiment)
Ds.f.in <- rowSums(DsX!=0)

farms$Dw.f.in <- Dw.f.in
farms$Ds.f.in <- Ds.f.in


# Indicator of OUTGOING connection strength:
# Sum of all D values (excluding self-recruitment) for lice coming FROM a farm
Dw.out <- colSums(DwX)
Ds.out <- colSums(DsX)

# Add to farms dataset
farms$Dw.out <- Dw.out
farms$Ds.out <- Ds.out

# The NUMBER of other farms that each farm sends any lice to
Dw.f.out <- colSums(DwX!=0) # Count number of farms with D>0 (excluding self-recruiment)
Ds.f.out <- colSums(DsX!=0)

farms$Dw.f.out <- Dw.f.out
farms$Ds.f.out <- Ds.f.out

# Calculate the strength of SELF-RECRUITMENT
Dw.self <- diag(Dw) # Count number of farms with D>0 (excluding self-recruiment)
Ds.self <- diag(Ds)

farms$Dw.self <- Dw.self
farms$Ds.self <- Ds.self

write.csv(farms, "src/2loci/farms.week1.csv")


##################################################

### Can also explore other measures of connectivity...?

# What about 2nd order connections?
# Of the other farms that farm_i sends lice to, how many connections do THOSE farms have?
# (excluding/including those farms that farm_i already disperses to)

# Can make a function for this, but for now might not need...

# E.g. FLUX: connectivity strength scaled by the size of farms
# (assuming that bigger farms produce more lice, therefore have a stronger infection pressure on other farms)
fish <- matrix(c(farms$fs), nrow=537, ncol=1) # vector of farm sizes
fluxw <- matrix(0, nrow=537, ncol=1) # empty vector for flux in winter...
fluxs <- matrix(0, nrow=537, ncol=1) # .. and summer

fluxw[,1] <- DwX %*% fish
fluxs[,1] <- DsX %*% fish

farms$fluxw <- fluxw
farms$fluxs <- fluxs


#### OUTGOING FLUX # 
# Estimated number of lice going from each farm (based on farm size), to each farm
fs.m <- matrix(rep(farms$fs, times=537), nrow=537, ncol=537, byrow=T)

Dw.outflux <- DwX * fs.m
Ds.outflux <- DsX * fs.m

# Add to farms dataset
farms$Dw.outflux <- colSums(Dw.outflux)
farms$Ds.outflux <- colSums(Ds.outflux)


write.csv(farms, "src/2loci/farms.week1.csv")

#######
# PLOTS

library(mapdata)
norwaymap <- map_data("worldHires", "Norway", xlim = c(4,11), ylim = c(58,66.5)) # Save map of Norway (within lat & long of sim)

# Dw vs Ds
ggplot(farms, 
       aes(x=Dw.out, y=Ds.out))+
  geom_point()+
  geom_abline(slope=1)+
  stat_smooth(alpha=0)
# Farms with strong Dw tend to also have strong Ds
# (but overall Dw stronger than Ds)
# Same pattern for INCOMING D

#D.in vs D.out
ggplot(farms, 
       aes(x=Dw.in, y=Dw.out))+
  geom_point()+
  geom_abline(slope=1)+
  stat_smooth(alpha=0)
