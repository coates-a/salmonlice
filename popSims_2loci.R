
# Call functions
source("src/2loci/popFun_2loci.R") # Population function 
source("src/2loci/MFun_2loci.R") # Function for building transition matrix M for each season 
source("src/2loci/treatmentsFun_2loci.R") # Function for recording treatments 


###### Create empty transition matrix G #######
# Build skeleton of stage transition matrix, S, for 1 genotype, 1 farm
S <- matrix(c(rep(0,5), 1, rep(0, 3), 1, 1, rep(0, 3), 1, 1), # Add 1 to cells to be filled by values for development rate (s)
             nrow=4, ncol=4, byrow=T)

G <- kronecker(diag(1,g), S) # Expand S into G (stage transition matrix for ALL genotypes on 1 farm)

G[c(seq(1,z,4)), c(seq(4,z,4))] <- 1 # Add 1 to cells to be filled by values for fecundity (f) (corresponding to Larvae rows, Adult columns)

#### Create vector of season names for every time-step of simulation ####
V <- rep(rep(c("w1","w2","w3","sp", # 10x 'seasons'in 1 year 
               "s1","s2","s3","s4","s5","w4"), # Temperature, larval dispersal, and lice limits are unique for each season
             times=c(5,5,5,6, # Each 'season' = 5 weeks (5 time-steps) except for spring (6 weeks with lower lice limits),
                     5,5,6,5,5,5)), # and s3 (for total 52 weeks in year)
         times = years) # Repeat seasonal cycle over years
V <- matrix(V,nrow=1)

#### Create vector of lice limits for each farm ####
l <- matrix(0, nrow=i, ncol=3)
# Multiple number of fish on farm (fs) by adult abundance (adults/fish) limit (k_) to calculate maximum NUMBER of adults allowed 
# Delousing treatment/forced harvest applied when adults exceed limit
l[,1] <- farms$fs * kws # Winter/summer limit that triggers treatment 
l[,2] <- farms$fs * ksp # Spring limit that triggers treatment 
l[,3] <- farms$fs * kmax # Limit that triggers forced harvest
colnames(l) <- c("winter/summer", "spring", "max")

#### Population matrix (N) ####
N <- matrix(0, nrow=size, ncol=ts) # To be filled by number of lice
rownames(N) <- rep(farms$farmID, each=z) # A row for every life stage, genotype and farm
colnames(N) <- c(0:(ts-1)) # A column for every time-step

#### Calculate number of lice at start of simulation (N_t=0) ####
# For each farm, multiply lice abundance (lice/fish) by number of fish (fs)

# New script: adjust starting lice numbers
ad.start <- matrix(0, nrow = i, ncol = 4)
pa.start <- matrix(0, nrow = i, ncol = 4)
ch.start <- matrix(0, nrow = i, ncol = 4)

# Starting lice abundance based on Barents Watch data
ad.start[, 1] <- farms$female * farms$fs # Assume number of adults = number adult females
pa.start[, 1] <- farms$motile * farms$fs # Assume all motiles = pre-adult
ch.start[, 1] <- farms$sessile * farms$fs

# Use lice numbers calculated after the 2-year spin-up period (sim = mech only)
ad.start[, 2] <- farms$ad.y2 * farms$fs # Assume number of adults = number adult females
pa.start[, 2] <- farms$pa.y2 * farms$fs # Assume all motiles = pre-adult
ch.start[, 2] <- farms$ch.y2 * farms$fs

# Constant lice abundance #1
ad.start[, 3] <- 0.2 * farms$fs # Assume number of adults = number adult females
pa.start[, 3] <- 0.2 * farms$fs # Assume all motiles = pre-adult
ch.start[, 3] <- 0.2 * farms$fs

# Constant lice abundance #2
ad.start[, 4] <- 0.1 * farms$fs # Assume number of adults = number adult females
pa.start[, 4] <- 0.1 * farms$fs # Assume all motiles = pre-adult
ch.start[, 4] <- 0.1 * farms$fs

adult.x <- ad.start[, lice.start] # column of starting lice numbetrs given by lice.start in paras.csv
preadult.x <- pa.start[, lice.start]
chalimus.x <- ch.start[, lice.start]

# Prev script
# farms$adult.x <- farms$female * farms$fs # Assume number of adults = number adult females
# farms$preadult.x <- farms$motile * farms$fs # Assume all motiles = pre-adult
# farms$chalimus.x <- farms$sessile * farms$fs


# Create a marker (r) for the 1st row in N for each farm (corresponding to RR Larvae)
r <- seq(from=1, to=size, by=z)
# Fill N_t=0 with starting lice numbers
for(rr in seq(3, z, 4)){ # Adults in every 4th row from row 4 (3+1)
   N[r+rr, 1] <- adult.x #farms$adult.x # Adults in every 4th row, from rows 4 - 12
}
for(rr in seq(2, z, 4)){
   N[r+rr, 1] <- preadult.x #farms$preadult.x # Pre-adults in every 4th row, from rows 3 - 11
}
for(rr in seq(1, z, 4)){
   N[r+rr, 1] <- chalimus.x #farms$chalimus.x # Chalimus in every 4th row, from rows 2 - 10
}


# Create a vector of all the starting genotype proportions that matches each stage and farm
start.prop <- c(ifelse(rep(farms$lat >= lwr.latR & farms$lat <= upr.latR & 
                              farms$lat >= lwr.latT & farms$lat <= upr.latT, # if farm within both R AND T range
                           each=z), # (Need to repeat for each geno and stage)
                       # Include starting props of RSUU, SSTU and RSTU
                       c(rep(c(0, 0, 0, 0,# No RRTT, RSTT, SSTT, RRTU
                               prop.RSTU, 
                               prop.SSTU, 
                               0, # No RRUU
                               prop.RSUU, 
                               1-prop.RSTU-prop.SSTU-prop.RSUU), # Remainder SSUU
                             each=j)), # repeat for each stage, so matches structure of N
                       
                       ifelse(rep(farms$lat >= lwr.latR & farms$lat <= upr.latR &
                                     (farms$lat < lwr.latT | farms$lat > upr.latT), # If farm ONLY within R range,
                                  each=z),
                              
                              rep(c(c(rep(0, times=7), 
                                      prop.RSUU, # Starting props of RSUU
                                      1-prop.RSUU)), each=j),
                              
                              ifelse(rep(farms$lat >= lwr.latT & farms$lat <= upr.latT &
                                            (farms$lat < lwr.latR | farms$lat > upr.latR), # If ONLY within T range 
                                         each=z),
                                     
                                     rep(c(c(rep(0, times=5), 
                                             prop.SSTU, # Starting props of SSTU
                                             0, 0, 
                                             1-prop.SSTU)), each=j),
                                     # If farm outside starting range of R AND T
                                     rep(c(c(rep(0, times=8), 
                                             1)), # 100% SSUU
                                         each=j))))) 


N[,1] <- N[,1]*start.prop # Multiply total lice numbers by starting genotype proportions to get starting lice numbers


PC <- prod.cycle[,1] # Vector of whether farms are active (1) or inactive (0) at t=1
PC <- matrix(rep(PC, each=z), # Repeat for each genotype & stage on farm
             nrow=size,
             ncol=1)
PC[c(seq(1, size, by=4)),] <- 1 # Put 1 in all Larva rows (larvae in environment are unaffected by farm status)

N[,1] <- N[,1] * PC # Multiply farms by 0 if inactive; by 1 in active


# Include larvae at start of simulation
N.ad <- c(N[seq(4, size, 4) ,1]) # Vector of N adult lice (each genotype)
f.w1 <- c(f.temp[,"w1"]) # Vector of reproductive output in first week of year
f.w1 <- rep(f.w1, each=g) # Repeat fecundity for each genotype
N.la <- N.ad * f.w1 # For now, just assume the genotype proportions of larvae are ths SAME as in adults
N[seq(1, size, 4) ,1] <- N.la # Add larvae based on output of adults already on farms (without gene mixing)



#####################################
########### RUN FUNCTION ############

popdyn <- pop_2loci(N, l, V, X, Y, Zdisc, MAX, prod.cycle)

####################################

#### Edit output ####

colnames(popdyn) <-  c("farm","t","N") # Assign column names
popdyn$abund <- popdyn$N/(rep(farms$fs, each=z)) # Calculate abundance (lice per fish) by dividing number of lice by number of fish on farm (fs)
popdyn$stage <- rep(c("la","ch","pa","ad"), times = g*i*ts) # Label life stages 
popdyn$prodcy <- rep(c(prod.cycle), each=z) # Log whether farm active (1) or fallow (0) for each farm and t 
                                            # (goes ONE COLUMN AT A TIME, so should match up with popdyn...)

popdyn$genoR <- rep(c("RR","RS","SS"), # Genotype at R/S locus
                    each=j, # R/S genotype repeats for each life stage
                    times = 3*i*ts) # and 3 times on each farm

popdyn$genoT <- rep(c(rep(c("TT","TU","UU"), # Genotype at T/U locus
                          each=j*3)), # each T/U genotype repeats for each life stage and 3 times in each farm
                    times = i*ts) 

popdyn$geno <- do.call(paste0, popdyn[c("genoR", "genoT")]) # Full genotype (both loci)



# Turn t into years
popdyn$yr <- rep(c(0:(years-1)), each=i*z*52) # Year of simulation (52 time-steps per year)
tts <- c(seq(0,51,1)/52) # Fraction of the year represented by each week
popdyn$tt <- rep(c(tts), each=i*z, times=years) # Repeat fraction for each year
popdyn$tt <- popdyn$tt + popdyn$yr # Add fraction of year to year of simulation for a continuous measure of time (e.g. tt = 2.5 means two and a half years)

N.stage <- aggregate(N ~ farm + stage + t, data=popdyn, FUN=sum) # Calculate the total number of lice in each stage (across all genotypes)
data.table::setnames(N.stage,'N','N.stage') # Set colname of number of lice from N to N.stage
popdyn <- merge(x=popdyn, y=N.stage, by=c("stage","farm","t"), all.y=F) # Merge with dataset

popdyn$gen.prop <- popdyn$N / popdyn$N.stage # Calculate the proportion of each genotype (per stage and farm)

abund.stage <- aggregate(abund ~ farm + stage + t, data=popdyn, FUN=sum) # Calculate the abudnance of lice in each stage (across all genotypes)
data.table::setnames(abund.stage,'abund','abund.stage')
popdyn <- merge(x=popdyn, y=abund.stage, by=c("stage","farm","t"), all.y=F) # Merge with dataset

##################### HERE ###############################

### CALCULATE GENE FREQUENCIES
# pR (frequency R)
# Total number each of RR, RS & SS lice
genoR.stage <- aggregate(N ~ stage + farm + t + genoR, data=popdyn, FUN = sum) # Sum each of the RR, RS and SS genotypes (across T/U genotypes)
data.table::setnames(genoR.stage,'N','genoR.stage') # Set colname of number of lice from N to genoR.stage

# Calculate total number of R alleles for each genotype
genoR.stage$nR <- genoR.stage$genoR.stage * ifelse(genoR.stage$genoR == "RR", # For each louse,
                                                   2, # 2x R alleles in RR
                                                   ifelse(genoR.stage$genoR == "RS", 
                                                          1, # 1x in RS
                                                          0)) # 0 in SS
nR <- aggregate(nR ~ stage + farm + t, data=genoR.stage, FUN = sum) # Now add up to total number of R alleles at each t, farm and stage
popdyn <- merge(x=popdyn, y=nR, by=c("stage","farm","t"), all.y=F) # Merge with dataset (nR will repeat for all genotypes in a stage and farm)
# Gene frequency of R allele (pR) for each stage & farm
popdyn$pR <- popdyn$nR / (2*popdyn$N.stage) # Number of R alleles / total number of alleles (2 per lice)


# pT (frequency T)
# Total number each of RR, RS & SS lice
genoT.stage <- aggregate(N ~ stage + farm + t + genoT, data=popdyn, FUN = sum) # Sum each of the RR, RS and SS genotypes (across T/U genotypes)
data.table::setnames(genoT.stage,'N','genoT.stage') # Set colname of number of lice from N to genoR.stage

# Calculate total number of R alleles for each genotype
genoT.stage$nT <- genoT.stage$genoT.stage * ifelse(genoT.stage$genoT == "TT", # For each louse,
                                                   2, # 2x T alleles in TT
                                                   ifelse(genoT.stage$genoT == "TU", 
                                                          1, # 1x in TU
                                                          0)) # 0 in UU
nT <- aggregate(nT ~ stage + farm + t, data=genoT.stage, FUN = sum) # Now add up to total number of R alleles at each t, farm and stage
popdyn <- merge(x=popdyn, y=nT, by=c("stage","farm","t"), all.y=F) # Merge with dataset (nR will repeat for all genotypes in a stage and farm)
# Gene frequency of R allele (pR) for each stage & farm
popdyn$pT <- popdyn$nT / (2*popdyn$N.stage) # Number of R alleles / total number of alleles (2 per lice)



#################################
# Merge with farm info (from farms.csv)
farms$farm <- farms$farmID # Ensure same column names
farms$H <- H.vector # Distribution of continuous strategies across farms (in this simulation)
##  ** For now, just using H AFTER new distribution is added (ie ignoring H.start)
popdyn <- merge(x=popdyn, y=farms[,c("lat", "lon","farm","H","zone")], by="farm", all.y=TRUE) # Match farms to their latitude & longitude, production zone, and whether or not using continuous strategy

#### Save output ####
write.csv(popdyn, file = sprintf('outputs/popdyn_2loci_%s_%s.csv', scenarios, sim), row.names = F)

#### Calculate the number of treatments applied through simulation ####
source("src/2loci/treatments_2loci.R") 

#### Create Summary data set for farms ####
farm.sum <- merge(x=farms, y=treat.f, by="farmID", all=T) # Add number of treatments per farm (calculated from 'treatments_3.R')
# Calculate the average gene frequencies (pR and pT) in adults on each farm in the final year of simulation
p.farms <- subset(popdyn, yr==max(yr) & # Subset final year
                     stage=="ad" & geno=="RRTT") # Only need to look at 1 genotype (since p is repeated)

# Average R gene frequency (pR) in adults across whole final year of sim, for each farm
ifelse(length(na.omit(p.farms$pR)) > 0, # p = NA (can't calculate mean) when there are no adults 
       pR.end <- aggregate(pR ~ farm, FUN=mean, # Calculate mean p over final year for each farm
                          data = p.farms), 
       pR.end <- data.frame(farm = farms$farm, # If all values of p are NA (i.e. NO ADULt LICE), create separate dataframe
                           pR = rep(NaN, times=i)))
colnames(pR.end)[1] <- "farmID"
farm.sum <- merge(x=farm.sum, y=pR.end, by="farmID", all=T) # Add mean pR in final year to farm summary

# Average T gene frequency (pT) in adults across whole final year of sim, for each farm
ifelse(length(na.omit(p.farms$pT)) > 0, # p = NA (can't calculate mean) when there are no adults 
       pT.end <- aggregate(pT ~ farm, FUN=mean, # Calculate mean p over final year for each farm
                           data = p.farms), 
       pT.end <- data.frame(farm = farms$farm, # If all values of p are NA (i.e. NO ADULt LICE), create separate dataframe
                            pT = rep(NaN, times=i)))
colnames(pT.end)[1] <- "farmID"
farm.sum <- merge(x=farm.sum, y=pT.end, by="farmID", all=T) # Add mean pR in final year to farm summary


#### Save farm summary ####
write.csv(farm.sum, file = sprintf("outputs/farm_summary_%s_%s.csv", scenarios, sim), row.names = T, col.names = T )

### CALCULATE N and GENE FREQUENCY (p) IN WHOLE METAPOPULATION ####

# Look just at adults for now
ads <- subset(popdyn, stage=="ad") 

gf <- aggregate(N ~ geno + genoR + genoT + tt, data=ads, FUN = sum) # Sum all adults of each genotype across all farms
Ntot <- aggregate(N ~ tt, data=ads, FUN=sum) # Total number of adults in sim
data.table::setnames(Ntot,'N','Ntot')
gf <- merge(x=gf, y=Ntot, by="tt", all.y=F)
gf$geno.prop <- gf$N / gf$Ntot # Proportion of genotypes in whole metapopulation
# Number of R alleles in whole adult metapopulation
gf$nR <- gf$N * ifelse(gf$genoR == "RR", # For each louse,
                       2, # 2x R alleles in RR
                       ifelse(gf$genoR == "RS", 
                              1, # 1x in RS
                              0)) # 0 in SS
nRtot <- aggregate(nR ~ tt, data=gf, FUN=sum)
data.table::setnames(nRtot,'nR','nRtot')
gf <- merge(x=gf, y=nRtot, by="tt")
gf$pR <- gf$nRtot / (2* gf$Ntot) # p = number R alleles / 2*number adults
# Number of T alleles in whole adult metapopulation
gf$nT <- gf$N * ifelse(gf$genoT == "TT", # For each louse,
                       2, # 2x T alleles in TT
                       ifelse(gf$genoT == "TU", 
                              1, # 1x in TU
                              0)) # 0 in UU
nTtot <- aggregate(nT ~ tt, data=gf, FUN=sum)
data.table::setnames(nTtot,'nT','nTtot')
gf <- merge(x=gf, y=nTtot, by="tt")
gf$pT <- gf$nTtot / (2* gf$Ntot) # p = number T alleles / 2*number adults


gf.p <- subset(gf, geno=="SSUU") # Just need to plot for 1 genotype, since p values are repeated for each genotype
gf.p <- gf.p[,c("tt","pR","pT")]
gf.p <- melt(gf.p, # Melt data into single column
             id.vars = "tt")
colnames(gf.p)[c(2,3)]<-c('locus', 'p')

colnames(gf.p)[3] <- sim
write.csv(gf.p, file = sprintf('outputs/p/p_%s_%s.csv', scenarios, sim), row.names = F) # Save back into src with new data

##### Calculate MEAN ABUNDACE across all farms #####

m.ab <- subset(popdyn, geno=="SSUU" & stage=="ad" &
                 prodcy==1) # Only average ACTIVE farms
m.ab <- aggregate(abund.stage ~ tt, data=m.ab, FUN = mean)
# Save mean abundance (across farms) in a separate folder - will combine all files in folder to compare abundance across sims
colnames(m.ab)[2] <- sim
write.csv(m.ab, file = sprintf('outputs/abund/mab_%s_%s.csv', scenarios, sim), row.names = F) # Save back into src with new data


#### Create plots ####
# source("src/2loci/mapPlot_2loci.R") 


