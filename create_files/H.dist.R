#### ASSIGN THE DISTRIBUTION OF CONTINUOUS/DISCRETE STRATEGIES ACROSS ALL FARMS ####
#### Create data set containing different distributions, to be used for different simulations
#### Here, distribution of Cont strategies is kept CONSTANT THROUGH TIME 

### Can have up to 3 types of strategies in the simulation: x, y, z 
### x and y can be either continuous or discrete (type is assigned in 'paras.csv')
### z is always discrete (only to be used if x and y are both continuous)
### Individual farms can use 3 strategies alone or in combination

library(readr)
library(tidyverse)
# Build upon existing H:
H <- read_csv("src/2loci/H.csv")
i <- length(H$farmID) # Number of farms

# Or, starting new:
farms <- read_csv("src/2loci/farms.week1.csv") # Data set of farm information
i <- length(farms$farmID) # Number of farms

H <- farms[,c("farmID","farm",
              "lon","lat","zone", # latitude, longitude, production zone
              "biomass", # farm size (biomass produced)
              "Dw.in","Ds.in", # strength of incoming lice connections (winter & summer)
              "Dw.out","Ds.out", # strength of outgoing lice connections
              "Dw.self","Ds.self", # strength of self-recruitment
              "Dw.f.in","Ds.f.in", # number of incoming lice connections
              "Dw.f.out","Ds.f.out", # number of incoming lice connections
              "avtemp")] # Average temperature

### Each new column of H = Distribution of strategies FOR ONE SCENARIO 
### (each row of H corresponds to 1 farm)
#### Edit the following script as necessary to produce the required distributions ####

#### ALL FARMS USING THE SAME STRATEGIES ####
H$all_x <- "x" 
H$all_y <- "y" 
H$all_xy <- "xy"
H$all_xz <- "xz"
H$all_yz <- "yz"
H$all_xyz <- "xyz"


########## 2 STRATEGIES #######

### COMBOS TO START WITH:
# - all XYZ
# - XZ vs YZ
# - Z vs XYZ

#### RANDOM DISTRIBUTION OF STRATEGIES AT SET PROPORTIONS ####

#### Create Function to randomly distribute 2 strategies ####
prop_2strats <- function(strat1, prop1, strat2){
  
  sample(rep(c(strat1, # 'sample' randomly distributes the two strategies
               strat2), 
             times=c(ceiling(i*prop1), # Proportion of strat1 = prop1
                     floor(i*(1-prop1))))) # ceiling/floor to round up and down (of i is an odd number)
}

# XZ vs YZ
H$r10xz_r90yz <- prop_2strats(strat1="xz", prop1=0.1, strat2="yz")
H$r25xz_r75yz <- prop_2strats(strat1="xz", prop1=0.25, strat2="yz")
H$r50xz_r50yz <- prop_2strats(strat1="xz", prop1=0.5, strat2="yz")
H$r75xz_r25yz <- prop_2strats(strat1="xz", prop1=0.75, strat2="yz")
H$r90xz_r10yz <- prop_2strats(strat1="xz", prop1=0.9, strat2="yz")

# Z vs XYZ
H$r10z_r90xyz <- prop_2strats(strat1="z", prop1=0.1, strat2="xyz")
H$r25z_r75xyz <- prop_2strats(strat1="z", prop1=0.25, strat2="xyz")
H$r50z_r50xyz <- prop_2strats(strat1="z", prop1=0.5, strat2="xyz")
H$r75z_r25xyz <- prop_2strats(strat1="z", prop1=0.75, strat2="xyz")
H$r90z_r10xyz <- prop_2strats(strat1="z", prop1=0.9, strat2="xyz")

# Z vs XZ
H$r50z_r50xz <- prop_2strats(strat1="z", prop1=0.5, strat2="xz")

# XZ vs XYZ
H$r50xz_r50xyz <- prop_2strats(strat1="xz", prop1=0.5, strat2="xyz")

## Run different random distributions ##
# 'sample' creates different distribution each time (no need to set seed) 
H$r50xz_r50yz_1 <- prop_2strats(strat1="xz", prop1=0.5, strat2="yz")
H$r75xz_r25yz_1 <- prop_2strats(strat1="xz", prop1=0.75, strat2="yz")
H$r90xz_r10yz_1 <- prop_2strats(strat1="xz", prop1=0.9, strat2="yz")
H$r50xz_r50yz_2 <- prop_2strats(strat1="xz", prop1=0.5, strat2="yz")
H$r75xz_r25yz_2 <- prop_2strats(strat1="xz", prop1=0.75, strat2="yz")
H$r90xz_r10yz_2 <- prop_2strats(strat1="xz", prop1=0.9, strat2="yz")
H$r50xz_r50yz_3 <- prop_2strats(strat1="xz", prop1=0.5, strat2="yz")
H$r75xz_r25yz_3 <- prop_2strats(strat1="xz", prop1=0.75, strat2="yz")
H$r90xz_r10yz_3 <- prop_2strats(strat1="xz", prop1=0.9, strat2="yz")

H$r10z_r90xyz_1 <- prop_2strats(strat1="z", prop1=0.1, strat2="xyz")
H$r25z_r75xyz_1 <- prop_2strats(strat1="z", prop1=0.25, strat2="xyz")
H$r50z_r50xyz_1 <- prop_2strats(strat1="z", prop1=0.5, strat2="xyz")
H$r75z_r25xyz_1 <- prop_2strats(strat1="z", prop1=0.75, strat2="xyz")
H$r90z_r10xyz_1 <- prop_2strats(strat1="z", prop1=0.9, strat2="xyz")
H$r10z_r90xyz_2 <- prop_2strats(strat1="z", prop1=0.1, strat2="xyz")
H$r25z_r75xyz_2 <- prop_2strats(strat1="z", prop1=0.25, strat2="xyz")
H$r50z_r50xyz_2 <- prop_2strats(strat1="z", prop1=0.5, strat2="xyz")
H$r75z_r25xyz_2 <- prop_2strats(strat1="z", prop1=0.75, strat2="xyz")
H$r90z_r10xyz_2 <- prop_2strats(strat1="z", prop1=0.9, strat2="xyz")
H$r10z_r90xyz_3 <- prop_2strats(strat1="z", prop1=0.1, strat2="xyz")
H$r25z_r75xyz_3 <- prop_2strats(strat1="z", prop1=0.25, strat2="xyz")
H$r50z_r50xyz_3 <- prop_2strats(strat1="z", prop1=0.5, strat2="xyz")
H$r75z_r25xyz_3 <- prop_2strats(strat1="z", prop1=0.75, strat2="xyz")
H$r90z_r10xyz_3 <- prop_2strats(strat1="z", prop1=0.9, strat2="xyz")


#########################################################
#### DISTRIBUTED STRATEGIES ACCORDING TO FARM VALUES ####

# Create Function for 2 strategies where farms are ordered by values (e.g. size, connectivity)

ordered_2strats <- function(lowest, # Strategy used on farms with the lowest values (eg. smallest farms)
                            proplow, # proportion of farms using 'lowest' strategy
                            highest){ # strategy used on farms with highest values (eg. biggest farms)
  
  (rep(c(lowest,
         highest), 
       times=c(ceiling(i*proplow), 
               floor(i*(1-proplow))))) # ceiling/floor to round up and down (of i is an odd number)
}



##########################################################
#### DISTRIBUTED ACCORDING TO SIZE OF FARM (BIOMASS) #####

H <- arrange(H, biomass) # Sort by farm size (ascending)

# xz VS yz (1 continuous strategy added to 1 discrete strategy at some sites)
H$small10_xz_big90_yz <- ordered_2strats("xz", 0.1, "yz") # Smallest 10% using xz, largest 90% using yz
H$small25_xz_big75_yz <- ordered_2strats("xz", 0.25, "yz")
H$small50_xz_big50_yz <- ordered_2strats("xz", 0.5, "yz") # Smallest 50% using xz, largest 50% using yz
H$small75_xz_big25_yz <- ordered_2strats("xz", 0.75, "yz")
H$small90_xz_big10_yz <- ordered_2strats("xz", 0.9, "yz")
# Switch xz and yz (smallest vs biggest farms) -> for when xz and yz have different effects
H$small10_yz_big90_xz <- ordered_2strats("yz", 0.1, "xz") # Smallest 10% using yz, largest 90% using xz
H$small25_yz_big75_xz <- ordered_2strats("yz", 0.25, "xz")
H$small50_yz_big50_xz <- ordered_2strats("yz", 0.5, "xz") # Smallest 50% using yz, largest 50% using xz
H$small75_yz_big25_xz <- ordered_2strats("yz", 0.75, "xz")
H$small90_yz_big10_xz <- ordered_2strats("yz", 0.9, "xz")

# xyz VS z (1 continuous strategy added to 1 discrete strategy at some sites)
H$small10_xyz_big90_z <- ordered_2strats("xyz", 0.1, "z") # Smallest 10% using xyz, largest 90% using z
H$small25_xyz_big75_z <- ordered_2strats("xyz", 0.25, "z")
H$small50_xyz_big50_z <- ordered_2strats("xyz", 0.5, "z") # Smallest 50% using xyz, largest 50% using z
H$small75_xyz_big25_z <- ordered_2strats("xyz", 0.75, "z")
H$small90_xyz_big10_z <- ordered_2strats("xyz", 0.9, "z")
# Switch xyz and z (smallest vs biggest farms) -> for when xyz and z have different effects
H$small10_z_big90_xyz <- ordered_2strats("z", 0.1, "xyz") # Smallest 10% using z, largest 90% using xyz
H$small25_z_big75_xyz <- ordered_2strats("z", 0.25, "xyz")
H$small50_z_big50_xyz <- ordered_2strats("z", 0.5, "xyz") # Smallest 50% using z, largest 50% using xyz
H$small75_z_big25_xyz <- ordered_2strats("z", 0.75, "xyz")
H$small90_z_big10_xyz <- ordered_2strats("z", 0.9, "xyz")



###########################################################
#### DISTRIBUTED ACCORDING TO FARM LARVAL CONNECTIVITY ####

##### Sort by strength of INCOMING WINTER connections #####
H <- arrange(H, Dw.in) 

# xz VS yz (1 continuous strategy added to 1 discrete strategy at some sites)
H$leastDwin10_xz_mostDwin90_yz <- ordered_2strats("xz", 0.1, "yz") 
H$leastDwin25_xz_mostDwin75_yz <- ordered_2strats("xz", 0.25, "yz")
H$leastDwin50_xz_mostDwin50_yz <- ordered_2strats("xz", 0.5, "yz") 
H$leastDwin75_xz_mostDwin25_yz <- ordered_2strats("xz", 0.75, "yz")
H$leastDwin90_xz_mostDwin10_yz <- ordered_2strats("xz", 0.9, "yz")
# Switch xz and yz
H$leastDwin10_yz_mostDwin90_xz <- ordered_2strats("yz", 0.1, "xz")  
H$leastDwin25_yz_mostDwin75_xz <- ordered_2strats("yz", 0.25, "xz")
H$leastDwin50_yz_mostDwin50_xz <- ordered_2strats("yz", 0.5, "xz") 
H$leastDwin75_yz_mostDwin25_xz <- ordered_2strats("yz", 0.75, "xz")
H$leastDwin90_yz_mostDwin10_xz <- ordered_2strats("yz", 0.9, "xz")

# xyz VS z (2 discrete strategies)
H$leastDwin10_xyz_mostDwin90_z <- ordered_2strats("xyz", 0.1, "z") # least connected 10% using xyz, most connected 90% using z
H$leastDwin25_xyz_mostDwin75_z <- ordered_2strats("xyz", 0.25, "z")
H$leastDwin50_xyz_mostDwin50_z <- ordered_2strats("xyz", 0.5, "z") # least connected 50% using xyz, most connected 50% using z
H$leastDwin75_xyz_mostDwin25_z <- ordered_2strats("xyz", 0.75, "z")
H$leastDwin90_xyz_mostDwin10_z <- ordered_2strats("xyz", 0.9, "z")
# Switch xyz and z -> for when xyz and z have different effects
H$leastDwin10_z_mostDwin90_xyz <- ordered_2strats("z", 0.1, "xyz") 
H$leastDwin25_z_mostDwin75_xyz <- ordered_2strats("z", 0.25, "xyz")
H$leastDwin50_z_mostDwin50_xyz <- ordered_2strats("z", 0.5, "xyz") 
H$leastDwin75_z_mostDwin25_xyz <- ordered_2strats("z", 0.75, "xyz")
H$leastDwin90_z_mostDwin10_xyz <- ordered_2strats("z", 0.9, "xyz")


##### Sort by strength of INCOMING SUMMER connections #####
H <- arrange(H, Ds.in) 

# xz VS yz (1 continuous strategy added to 1 discrete strategy at some sites)
H$leastDsin10_xz_mostDsin90_yz <- ordered_2strats("xz", 0.1, "yz") 
H$leastDsin25_xz_mostDsin75_yz <- ordered_2strats("xz", 0.25, "yz")
H$leastDsin50_xz_mostDsin50_yz <- ordered_2strats("xz", 0.5, "yz") 
H$leastDsin75_xz_mostDsin25_yz <- ordered_2strats("xz", 0.75, "yz")
H$leastDsin90_xz_mostDsin10_yz <- ordered_2strats("xz", 0.9, "yz")
# Switch xz and yz
H$leastDsin10_yz_mostDsin90_xz <- ordered_2strats("yz", 0.1, "xz")  
H$leastDsin25_yz_mostDsin75_xz <- ordered_2strats("yz", 0.25, "xz")
H$leastDsin50_yz_mostDsin50_xz <- ordered_2strats("yz", 0.5, "xz") 
H$leastDsin75_yz_mostDsin25_xz <- ordered_2strats("yz", 0.75, "xz")
H$leastDsin90_yz_mostDsin10_xz <- ordered_2strats("yz", 0.9, "xz")

# xyz VS z (2 discrete strategies)
H$leastDsin10_xyz_mostDsin90_z <- ordered_2strats("xyz", 0.1, "z") # least connected 10% using xyz, most connected 90% using z
H$leastDsin25_xyz_mostDsin75_z <- ordered_2strats("xyz", 0.25, "z")
H$leastDsin50_xyz_mostDsin50_z <- ordered_2strats("xyz", 0.5, "z") # least connected 50% using xyz, most connected 50% using z
H$leastDsin75_xyz_mostDsin25_z <- ordered_2strats("xyz", 0.75, "z")
H$leastDsin90_xyz_mostDsin10_z <- ordered_2strats("xyz", 0.9, "z")
# Switch xyz and z -> for when xyz and z have different effects
H$leastDsin10_z_mostDsin90_xyz <- ordered_2strats("z", 0.1, "xyz") 
H$leastDsin25_z_mostDsin75_xyz <- ordered_2strats("z", 0.25, "xyz")
H$leastDsin50_z_mostDsin50_xyz <- ordered_2strats("z", 0.5, "xyz") 
H$leastDsin75_z_mostDsin25_xyz <- ordered_2strats("z", 0.75, "xyz")
H$leastDsin90_z_mostDsin10_xyz <- ordered_2strats("z", 0.9, "xyz")

##### Sort by strength of OUTGOING WINTER connections #####
H <- arrange(H, Dw.out) 

# xz VS yz (1 continuous strategy added to 1 discrete strategy at some sites)
H$leastDwout10_xz_mostDwout90_yz <- ordered_2strats("xz", 0.1, "yz") 
H$leastDwout25_xz_mostDwout75_yz <- ordered_2strats("xz", 0.25, "yz")
H$leastDwout50_xz_mostDwout50_yz <- ordered_2strats("xz", 0.5, "yz") 
H$leastDwout75_xz_mostDwout25_yz <- ordered_2strats("xz", 0.75, "yz")
H$leastDwout90_xz_mostDwout10_yz <- ordered_2strats("xz", 0.9, "yz")
# Switch xz and yz
H$leastDwout10_yz_mostDwout90_xz <- ordered_2strats("yz", 0.1, "xz")  
H$leastDwout25_yz_mostDwout75_xz <- ordered_2strats("yz", 0.25, "xz")
H$leastDwout50_yz_mostDwout50_xz <- ordered_2strats("yz", 0.5, "xz") 
H$leastDwout75_yz_mostDwout25_xz <- ordered_2strats("yz", 0.75, "xz")
H$leastDwout90_yz_mostDwout10_xz <- ordered_2strats("yz", 0.9, "xz")

# xyz VS z (2 discrete strategies)
H$leastDwout10_xyz_mostDwout90_z <- ordered_2strats("xyz", 0.1, "z") # least connected 10% usoutg xyz, most connected 90% usoutg z
H$leastDwout25_xyz_mostDwout75_z <- ordered_2strats("xyz", 0.25, "z")
H$leastDwout50_xyz_mostDwout50_z <- ordered_2strats("xyz", 0.5, "z") # least connected 50% usoutg xyz, most connected 50% usoutg z
H$leastDwout75_xyz_mostDwout25_z <- ordered_2strats("xyz", 0.75, "z")
H$leastDwout90_xyz_mostDwout10_z <- ordered_2strats("xyz", 0.9, "z")
# Switch xyz and z -> for when xyz and z have different effects
H$leastDwout10_z_mostDwout90_xyz <- ordered_2strats("z", 0.1, "xyz") 
H$leastDwout25_z_mostDwout75_xyz <- ordered_2strats("z", 0.25, "xyz")
H$leastDwout50_z_mostDwout50_xyz <- ordered_2strats("z", 0.5, "xyz") 
H$leastDwout75_z_mostDwout25_xyz <- ordered_2strats("z", 0.75, "xyz")
H$leastDwout90_z_mostDwout10_xyz <- ordered_2strats("z", 0.9, "xyz")

# Tweak proportion refugia:
H$leastDwout85_xyz_mostDwout15_z <- ordered_2strats("xyz", 0.85, "z")
H$leastDwout95_xyz_mostDwout5_z <- ordered_2strats("xyz", 0.95, "z")

##### Sort by strength of OUTGOING SUMMER connections #####
H <- arrange(H, Ds.out) 

# xz VS yz (1 continuous strategy added to 1 discrete strategy at some sites)
H$leastDsout10_xz_mostDsout90_yz <- ordered_2strats("xz", 0.1, "yz") 
H$leastDsout25_xz_mostDsout75_yz <- ordered_2strats("xz", 0.25, "yz")
H$leastDsout50_xz_mostDsout50_yz <- ordered_2strats("xz", 0.5, "yz") 
H$leastDsout75_xz_mostDsout25_yz <- ordered_2strats("xz", 0.75, "yz")
H$leastDsout90_xz_mostDsout10_yz <- ordered_2strats("xz", 0.9, "yz")
# Switch xz and yz
H$leastDsout10_yz_mostDsout90_xz <- ordered_2strats("yz", 0.1, "xz")  
H$leastDsout25_yz_mostDsout75_xz <- ordered_2strats("yz", 0.25, "xz")
H$leastDsout50_yz_mostDsout50_xz <- ordered_2strats("yz", 0.5, "xz") 
H$leastDsout75_yz_mostDsout25_xz <- ordered_2strats("yz", 0.75, "xz")
H$leastDsout90_yz_mostDsout10_xz <- ordered_2strats("yz", 0.9, "xz")

# xyz VS z (2 discrete strategies)
H$leastDsout10_xyz_mostDsout90_z <- ordered_2strats("xyz", 0.1, "z") # least connected 10% usoutg xyz, most connected 90% usoutg z
H$leastDsout25_xyz_mostDsout75_z <- ordered_2strats("xyz", 0.25, "z")
H$leastDsout50_xyz_mostDsout50_z <- ordered_2strats("xyz", 0.5, "z") # least connected 50% usoutg xyz, most connected 50% usoutg z
H$leastDsout75_xyz_mostDsout25_z <- ordered_2strats("xyz", 0.75, "z")
H$leastDsout90_xyz_mostDsout10_z <- ordered_2strats("xyz", 0.9, "z")
# Switch xyz and z -> for when xyz and z have different effects
H$leastDsout10_z_mostDsout90_xyz <- ordered_2strats("z", 0.1, "xyz") 
H$leastDsout25_z_mostDsout75_xyz <- ordered_2strats("z", 0.25, "xyz")
H$leastDsout50_z_mostDsout50_xyz <- ordered_2strats("z", 0.5, "xyz") 
H$leastDsout75_z_mostDsout25_xyz <- ordered_2strats("z", 0.75, "xyz")
H$leastDsout90_z_mostDsout10_xyz <- ordered_2strats("z", 0.9, "xyz")



##### Sort by strength of NUMBER of OUTGOING SUMMER connections #####
H <- arrange(H, Ds.f.out) 

# xz VS yz (1 continuous strategy added to 1 discrete strategy at some sites)
H$leastDsf.out10_xz_mostDsf.out90_yz <- ordered_2strats("xz", 0.1, "yz") 
H$leastDsf.out25_xz_mostDsf.out75_yz <- ordered_2strats("xz", 0.25, "yz")
H$leastDsf.out50_xz_mostDsf.out50_yz <- ordered_2strats("xz", 0.5, "yz") 
H$leastDsf.out75_xz_mostDsf.out25_yz <- ordered_2strats("xz", 0.75, "yz")
H$leastDsf.out90_xz_mostDsf.out10_yz <- ordered_2strats("xz", 0.9, "yz")
# Switch xz and yz
H$leastDsf.out10_yz_mostDsf.out90_xz <- ordered_2strats("yz", 0.1, "xz")  
H$leastDsf.out25_yz_mostDsf.out75_xz <- ordered_2strats("yz", 0.25, "xz")
H$leastDsf.out50_yz_mostDsf.out50_xz <- ordered_2strats("yz", 0.5, "xz") 
H$leastDsf.out75_yz_mostDsf.out25_xz <- ordered_2strats("yz", 0.75, "xz")
H$leastDsf.out90_yz_mostDsf.out10_xz <- ordered_2strats("yz", 0.9, "xz")

# xyz VS z (2 discrete strategies)
H$leastDsf.out10_xyz_mostDsf.out90_z <- ordered_2strats("xyz", 0.1, "z") # least connected 10% usf.outg xyz, most connected 90% usf.outg z
H$leastDsf.out25_xyz_mostDsf.out75_z <- ordered_2strats("xyz", 0.25, "z")
H$leastDsf.out50_xyz_mostDsf.out50_z <- ordered_2strats("xyz", 0.5, "z") # least connected 50% usf.outg xyz, most connected 50% usf.outg z
H$leastDsf.out75_xyz_mostDsf.out25_z <- ordered_2strats("xyz", 0.75, "z")
H$leastDsf.out90_xyz_mostDsf.out10_z <- ordered_2strats("xyz", 0.9, "z")
# Switch xyz and z -> for when xyz and z have different effects
H$leastDsf.out10_z_mostDsf.out90_xyz <- ordered_2strats("z", 0.1, "xyz") 
H$leastDsf.out25_z_mostDsf.out75_xyz <- ordered_2strats("z", 0.25, "xyz")
H$leastDsf.out50_z_mostDsf.out50_xyz <- ordered_2strats("z", 0.5, "xyz") 
H$leastDsf.out75_z_mostDsf.out25_xyz <- ordered_2strats("z", 0.75, "xyz")
H$leastDsf.out90_z_mostDsf.out10_xyz <- ordered_2strats("z", 0.9, "xyz")


##### Sort by strength of NUMBER of OUTGOING WINTER connections #####
H <- arrange(H, Dw.f.out) 

# xz VS yz (1 continuous strategy added to 1 discrete strategy at some sites)
H$leastDwf.out10_xz_mostDwf.out90_yz <- ordered_2strats("xz", 0.1, "yz") 
H$leastDwf.out25_xz_mostDwf.out75_yz <- ordered_2strats("xz", 0.25, "yz")
H$leastDwf.out50_xz_mostDwf.out50_yz <- ordered_2strats("xz", 0.5, "yz") 
H$leastDwf.out75_xz_mostDwf.out25_yz <- ordered_2strats("xz", 0.75, "yz")
H$leastDwf.out90_xz_mostDwf.out10_yz <- ordered_2strats("xz", 0.9, "yz")
# Switch xz and yz
H$leastDwf.out10_yz_mostDwf.out90_xz <- ordered_2strats("yz", 0.1, "xz")  
H$leastDwf.out25_yz_mostDwf.out75_xz <- ordered_2strats("yz", 0.25, "xz")
H$leastDwf.out50_yz_mostDwf.out50_xz <- ordered_2strats("yz", 0.5, "xz") 
H$leastDwf.out75_yz_mostDwf.out25_xz <- ordered_2strats("yz", 0.75, "xz")
H$leastDwf.out90_yz_mostDwf.out10_xz <- ordered_2strats("yz", 0.9, "xz")

# xyz VS z (2 discrete strategies)
H$leastDwf.out10_xyz_mostDwf.out90_z <- ordered_2strats("xyz", 0.1, "z") # least connected 10% usf.outg xyz, most connected 90% usf.outg z
H$leastDwf.out25_xyz_mostDwf.out75_z <- ordered_2strats("xyz", 0.25, "z")
H$leastDwf.out50_xyz_mostDwf.out50_z <- ordered_2strats("xyz", 0.5, "z") # least connected 50% usf.outg xyz, most connected 50% usf.outg z
H$leastDwf.out75_xyz_mostDwf.out25_z <- ordered_2strats("xyz", 0.75, "z")
H$leastDwf.out90_xyz_mostDwf.out10_z <- ordered_2strats("xyz", 0.9, "z")
# Switch xyz and z -> for when xyz and z have different effects
H$leastDwf.out10_z_mostDwf.out90_xyz <- ordered_2strats("z", 0.1, "xyz") 
H$leastDwf.out25_z_mostDwf.out75_xyz <- ordered_2strats("z", 0.25, "xyz")
H$leastDwf.out50_z_mostDwf.out50_xyz <- ordered_2strats("z", 0.5, "xyz") 
H$leastDwf.out75_z_mostDwf.out25_xyz <- ordered_2strats("z", 0.75, "xyz")
H$leastDwf.out90_z_mostDwf.out10_xyz <- ordered_2strats("z", 0.9, "xyz")


##### Sort by strength of NUMBER of OUTGOING SUMMER connections #####
H <- arrange(H, Ds.f.in) 

# xz VS yz (1 continuous strategy added to 1 discrete strategy at some sites)
H$leastDsf.in10_xz_mostDsf.in90_yz <- ordered_2strats("xz", 0.1, "yz") 
H$leastDsf.in25_xz_mostDsf.in75_yz <- ordered_2strats("xz", 0.25, "yz")
H$leastDsf.in50_xz_mostDsf.in50_yz <- ordered_2strats("xz", 0.5, "yz") 
H$leastDsf.in75_xz_mostDsf.in25_yz <- ordered_2strats("xz", 0.75, "yz")
H$leastDsf.in90_xz_mostDsf.in10_yz <- ordered_2strats("xz", 0.9, "yz")
# Switch xz and yz
H$leastDsf.in10_yz_mostDsf.in90_xz <- ordered_2strats("yz", 0.1, "xz")  
H$leastDsf.in25_yz_mostDsf.in75_xz <- ordered_2strats("yz", 0.25, "xz")
H$leastDsf.in50_yz_mostDsf.in50_xz <- ordered_2strats("yz", 0.5, "xz") 
H$leastDsf.in75_yz_mostDsf.in25_xz <- ordered_2strats("yz", 0.75, "xz")
H$leastDsf.in90_yz_mostDsf.in10_xz <- ordered_2strats("yz", 0.9, "xz")

# xyz VS z (2 discrete strategies)
H$leastDsf.in10_xyz_mostDsf.in90_z <- ordered_2strats("xyz", 0.1, "z") # least connected 10% usf.ing xyz, most connected 90% usf.ing z
H$leastDsf.in25_xyz_mostDsf.in75_z <- ordered_2strats("xyz", 0.25, "z")
H$leastDsf.in50_xyz_mostDsf.in50_z <- ordered_2strats("xyz", 0.5, "z") # least connected 50% usf.ing xyz, most connected 50% usf.ing z
H$leastDsf.in75_xyz_mostDsf.in25_z <- ordered_2strats("xyz", 0.75, "z")
H$leastDsf.in90_xyz_mostDsf.in10_z <- ordered_2strats("xyz", 0.9, "z")
# Switch xyz and z -> for when xyz and z have different effects
H$leastDsf.in10_z_mostDsf.in90_xyz <- ordered_2strats("z", 0.1, "xyz") 
H$leastDsf.in25_z_mostDsf.in75_xyz <- ordered_2strats("z", 0.25, "xyz")
H$leastDsf.in50_z_mostDsf.in50_xyz <- ordered_2strats("z", 0.5, "xyz") 
H$leastDsf.in75_z_mostDsf.in25_xyz <- ordered_2strats("z", 0.75, "xyz")
H$leastDsf.in90_z_mostDsf.in10_xyz <- ordered_2strats("z", 0.9, "xyz")


##### Sort by strength of NUMBER of OUTGOING WINTER connections #####
H <- arrange(H, Dw.f.in) 

# xz VS yz (1 continuous strategy added to 1 discrete strategy at some sites)
H$leastDwf.in10_xz_mostDwf.in90_yz <- ordered_2strats("xz", 0.1, "yz") 
H$leastDwf.in25_xz_mostDwf.in75_yz <- ordered_2strats("xz", 0.25, "yz")
H$leastDwf.in50_xz_mostDwf.in50_yz <- ordered_2strats("xz", 0.5, "yz") 
H$leastDwf.in75_xz_mostDwf.in25_yz <- ordered_2strats("xz", 0.75, "yz")
H$leastDwf.in90_xz_mostDwf.in10_yz <- ordered_2strats("xz", 0.9, "yz")
# Switch xz and yz
H$leastDwf.in10_yz_mostDwf.in90_xz <- ordered_2strats("yz", 0.1, "xz")  
H$leastDwf.in25_yz_mostDwf.in75_xz <- ordered_2strats("yz", 0.25, "xz")
H$leastDwf.in50_yz_mostDwf.in50_xz <- ordered_2strats("yz", 0.5, "xz") 
H$leastDwf.in75_yz_mostDwf.in25_xz <- ordered_2strats("yz", 0.75, "xz")
H$leastDwf.in90_yz_mostDwf.in10_xz <- ordered_2strats("yz", 0.9, "xz")

# xyz VS z (2 discrete strategies)
H$leastDwf.in10_xyz_mostDwf.in90_z <- ordered_2strats("xyz", 0.1, "z") # least connected 10% usf.ing xyz, most connected 90% usf.ing z
H$leastDwf.in25_xyz_mostDwf.in75_z <- ordered_2strats("xyz", 0.25, "z")
H$leastDwf.in50_xyz_mostDwf.in50_z <- ordered_2strats("xyz", 0.5, "z") # least connected 50% usf.ing xyz, most connected 50% usf.ing z
H$leastDwf.in75_xyz_mostDwf.in25_z <- ordered_2strats("xyz", 0.75, "z")
H$leastDwf.in90_xyz_mostDwf.in10_z <- ordered_2strats("xyz", 0.9, "z")
# Switch xyz and z -> for when xyz and z have different effects
H$leastDwf.in10_z_mostDwf.in90_xyz <- ordered_2strats("z", 0.1, "xyz") 
H$leastDwf.in25_z_mostDwf.in75_xyz <- ordered_2strats("z", 0.25, "xyz")
H$leastDwf.in50_z_mostDwf.in50_xyz <- ordered_2strats("z", 0.5, "xyz") 
H$leastDwf.in75_z_mostDwf.in25_xyz <- ordered_2strats("z", 0.75, "xyz")
H$leastDwf.in90_z_mostDwf.in10_xyz <- ordered_2strats("z", 0.9, "xyz")


##### Sort by OUTGOING WINTER FLUX (Connectivity scaled by source farm size) #####
H <- arrange(H, Dw.outflux) 

# xz VS yz (1 continuous strategy added to 1 discrete strategy at some sites)

H$leastDw.outflux90_xyz_most10_z <- ordered_2strats("xyz", 0.9, "z")



#############################################################
#### DISTRIBUTED ACCORDING TO AVERAGE YEARLY TEMPERATURE ####

H <- arrange(H, avtemp) 

# xz VS yz (2 discrete strategies)
H$cold10_xz_hot90_yz <- ordered_2strats("xz", 0.1, "yz") # least connected 10% using xz, most connected 90% using yz
H$cold25_xz_hot75_yz <- ordered_2strats("xz", 0.25, "yz")
H$cold50_xz_hot50_yz <- ordered_2strats("xz", 0.5, "yz") # least connected 50% using xz, most connected 50% using yz
H$cold75_xz_hot25_yz <- ordered_2strats("xz", 0.75, "yz")
H$cold90_xz_hot10_yz <- ordered_2strats("xz", 0.9, "yz")
# Switch xz and yz -> for when xz and yz have different effects
H$cold10_yz_hot90_xz <- ordered_2strats("yz", 0.1, "xz") 
H$cold25_yz_hot75_xz <- ordered_2strats("yz", 0.25, "xz")
H$cold50_yz_hot50_xz <- ordered_2strats("yz", 0.5, "xz") 
H$cold75_yz_hot25_xz <- ordered_2strats("yz", 0.75, "xz")
H$cold90_yz_hot10_xz <- ordered_2strats("yz", 0.9, "xz")

# xyz VS z (1 continuous strategy added to 1 discrete strategy at some sites)
H$cold10_xyz_hot90_z <- ordered_2strats("xyz", 0.1, "z") 
H$cold25_xyz_hot75_z <- ordered_2strats("xyz", 0.25, "z")
H$cold50_xyz_hot50_z <- ordered_2strats("xyz", 0.5, "z") 
H$cold75_xyz_hot25_z <- ordered_2strats("xyz", 0.75, "z")
H$cold90_xyz_hot10_z <- ordered_2strats("xyz", 0.9, "z")
# Switch xyz and z
H$cold10_z_hot90_xyz <- ordered_2strats("z", 0.1, "xyz")  
H$cold25_z_hot75_xyz <- ordered_2strats("z", 0.25, "xyz")
H$cold50_z_hot50_xyz <- ordered_2strats("z", 0.5, "xyz") 
H$cold75_z_hot25_xyz <- ordered_2strats("z", 0.75, "xyz")
H$cold90_z_hot10_xyz <- ordered_2strats("z", 0.9, "xyz")

#####################################################################

#### Alternate by latitude ####
H$altlat_xz_yz <- ifelse(floor(H$lat) %% 2 == 0, # if latitude (rounded down) is even
                       "xz",
                       "yz")
# Swap xz & yz (lowest latitude is EVEN)
H$altlat_yz_xz <- ifelse(floor(H$lat) %% 2 == 0, # if latitude (rounded down) is even
                       "yz",
                       "xz")
## xyz vs z
H$altlat_xyz_z <- ifelse(floor(H$lat) %% 2 == 0, # if latitude (rounded down) is even
                       "xyz",
                       "z")
H$altlat_z_xyz <- ifelse(floor(H$lat) %% 2 == 0, # if latitude (rounded down) is even
                       "z",
                       "xyz")


#### Alternate byz PRODUCTION ZONE ####
H$altzone_xz_yz <- ifelse(floor(H$zone) %% 2 == 0, # if zone is even
                       "yz",
                       "xz") # Zone 1 = xz
# Swap xz & yz (Zone 1 = yz)
H$altzone_yz_xz <- ifelse(floor(H$zone) %% 2 == 0,  
                       "xz",
                       "yz")
## xyz vs z
H$altzone_xyz_z <- ifelse(floor(H$zone) %% 2 == 0, 
                        "z",
                        "xyz")
H$altzone_z_xyz <- ifelse(floor(H$zone) %% 2 == 0, 
                        "xyz",
                        "z")



### IMPORTANT! ###
H <- arrange(H, farmID) # Need to rearrange order of data back to order of farmID, so matches with original farms.csv data set!
write.csv(H, file='src/2loci/H.csv')

##########################################

######### 3 STRATEGIES ######

#### Create Function to randomly distribute 3 strategies ####
prop_3strats <- function(strat1, prop1, 
                         strat2, prop2,
                         strat3){
  sample(rep(c(strat1, # sample randomly distributes the three strategies
               strat2,
               strat3), 
             times=c(ceiling(i*prop1),
                     ceiling(i*prop2),
                     ceiling(i*(1-prop1-prop2)))) # round up
  )[1:i] # Then ensure vector is #i values long 
}

# Eg
# Equal proportions X, Y, XY
H$r33xy_r33x_r33y <- prop_3strats("xy", prop1=0.333,
                                  "x", prop2=0.333, 
                                  "y")
# Equal proportions XZ, YZ, XYZ
H$r33xyz_r33xz_r33yz <- prop_3strats("xyz", prop1=0.333,
                                     "xz", prop2=0.333, 
                                     "yz")
##
# Equal proportions XZ, YZ, Z
H$r33xz_r33yz_r33z <- prop_3strats("xz", prop1=0.333,
                                     "yz", prop2=0.333, 
                                     "z")
# Equal proportions Z, XZ, XYZ
H$r33z_r33xz_r33xyz <- prop_3strats("z", prop1=0.333,
                                     "xz", prop2=0.333, 
                                     "xyz")

# Unequal x, y, xy
H$r20xy_r40x_r40y <- prop_3strats("xy", prop1=0.2,
                                  "x", prop2=0.4, 
                                  "y")







############################################################################

#### Distribute 2 strategy types within 1 Production Zone (rest of area just use 1)

H <- arrange(H, farmID) # Make sure H is back in original order

pz.2strat <- function(zone1, # Zone receiving new strategy distribution
                  regstrat, # Regular strategy used throughout area
                  newstrat, # New strategy added to Zone
                    propnewstrat){ # Proportion of farms in Zone using newstrat

pz <- subset(H, zone==zone1)["farmID"] # Subset just farm IDs in that Zone
n_pz <- nrow(pz) # Total number of farms in zone

pz$dist <- sample(rep(c(regstrat, # Randomly distribute 2 strategies through zone
             newstrat), 
           times=c(floor(n_pz*(1-propnewstrat)), # according to proportion of new strat
                   ceiling(n_pz*propnewstrat)))) # Round up in number of newstrat (for zones with very few farms)

pz <- merge(H["farmID"], pz, by="farmID", all=T) # Add full list of farms (across all zones)

pz$dist[is.na(pz$dist)] <- regstrat # Replace NAs with regstrat

pz <- arrange(pz, farmID) # Make sure back in original order (by farm ID)

return(c(pz$dist)) # Just take distribution column out of new data frame
}

# E.g. Focus on Zone 3
# Most farms using XZ, some Zone 3 using YZ
H$xz_zone3.r25yz <- pz.2strat(zone1 = 3,
                             regstrat = "xz", # Most farms using XZ
                             newstrat = "yz", propnewstrat = 0.25) # Except 25% of farms in Zone 3 using YZ
H$xz_zone3.r50yz <- pz.2strat(zone1 = 3,
                              regstrat = "xz", # Most farms using XZ
                              newstrat = "yz", propnewstrat = 0.5) # Except 50% of farms in Zone 3 using YZ
H$xz_zone3.r75yz <- pz.2strat(zone1 = 3,
                              regstrat = "xz", # Most farms using XZ
                              newstrat = "yz", propnewstrat = 0.75) # Except 75% of farms in Zone 3 using YZ
H$xz_zone3.r100yz <- pz.2strat(zone1 = 3,
                              regstrat = "xz", # Most farms using XZ
                              newstrat = "yz", propnewstrat = 1) # Except 100% of farms in Zone 3 using YZ

# Most farms using Z, some Zone 3 using XYZ
H$z_zone3.r50xyz <- pz.2strat(zone1 = 3,
                              regstrat = "z", # Most farms using Z
                              newstrat = "xyz", propnewstrat = 0.5) # Except 50% of farms in Zone 3 using XYZ
H$z_zone3.r75xyz <- pz.2strat(zone1 = 3,
                              regstrat = "z", # Most farms using Z
                              newstrat = "xyz", propnewstrat = 0.75)
# VS. Most farms using XYZ, some Zone 3 using just Z
H$xyz_zone3.r50z <- pz.2strat(zone1 = 3,
                              regstrat = "xyz", # Most farms using Z
                              newstrat = "z", propnewstrat = 0.5) # Except 50% of farms in Zone 3 using XYZ
H$xyz_zone3.r75z <- pz.2strat(zone1 = 3,
                              regstrat = "xyz", # Most farms using Z
                              newstrat = "z", propnewstrat = 0.75) # Except 75% of farms in Zone 3 using XYZ                             


### Repeat for others...
H$x_zone1_r25xy <- pz.2strat(zone1 = 1, # 25% of farms in Zone 1 using xy, rest using x
                              regstrat = "x",
                              newstrat = "xy",
                              propnewstrat = 0.25)
H$x_zone1_r50xy <- pz.2strat(zone1 = 1,
                             regstrat = "x",
                             newstrat = "xy",
                             propnewstrat = 0.5)
H$x_zone1_r75xy <- pz.2strat(zone1 = 1,
                             regstrat = "x",
                             newstrat = "xy",
                             propnewstrat = 0.75)
# Z 2
H$x_zone2_r25xy <- pz.2strat(zone1 = 2, # 25% of farms in Zone 2 using xy, rest using x
                             regstrat = "x",
                             newstrat = "xy",
                             propnewstrat = 0.25)
H$x_zone2_r50xy <- pz.2strat(zone1 = 2,
                             regstrat = "x",
                             newstrat = "xy",
                             propnewstrat = 0.5)
H$x_zone2_r75xy <- pz.2strat(zone1 = 2,
                             regstrat = "x",
                             newstrat = "xy",
                             propnewstrat = 0.75)
# Z 3
H$x_zone3_r25xy <- pz.2strat(zone1 = 3, # 25% of farms in Zone 2 using xy, rest using x
                             regstrat = "x",
                             newstrat = "xy",
                             propnewstrat = 0.25)
H$x_zone3_r50xy <- pz.2strat(zone1 = 3,
                             regstrat = "x",
                             newstrat = "xy",
                             propnewstrat = 0.5)
H$x_zone3_r75xy <- pz.2strat(zone1 = 3,
                             regstrat = "x",
                             newstrat = "xy",
                             propnewstrat = 0.75)
# Z 4
H$x_zone4_r25xy <- pz.2strat(zone1 = 4, # 25% of farms in Zone 2 using xy, rest using x
                             regstrat = "x",
                             newstrat = "xy",
                             propnewstrat = 0.25)
H$x_zone4_r50xy <- pz.2strat(zone1 = 4,
                             regstrat = "x",
                             newstrat = "xy",
                             propnewstrat = 0.5)
H$x_zone4_r75xy <- pz.2strat(zone1 = 4,
                             regstrat = "x",
                             newstrat = "xy",
                             propnewstrat = 0.75)
# Z 5
H$x_zone5_r25xy <- pz.2strat(zone1 = 5, # 25% of farms in Zone 2 using xy, rest using x
                             regstrat = "x",
                             newstrat = "xy",
                             propnewstrat = 0.25)
H$x_zone5_r50xy <- pz.2strat(zone1 = 5,
                             regstrat = "x",
                             newstrat = "xy",
                             propnewstrat = 0.5)
H$x_zone5_r75xy <- pz.2strat(zone1 = 5,
                             regstrat = "x",
                             newstrat = "xy",
                             propnewstrat = 0.75)

# Z 6
H$x_zone6_r25xy <- pz.2strat(zone1 = 6, # 25% of farms in Zone 2 using xy, rest using x
                             regstrat = "x",
                             newstrat = "xy",
                             propnewstrat = 0.25)
H$x_zone6_r50xy <- pz.2strat(zone1 = 6,
                             regstrat = "x",
                             newstrat = "xy",
                             propnewstrat = 0.5)
H$x_zone6_r75xy <- pz.2strat(zone1 = 6,
                             regstrat = "x",
                             newstrat = "xy",
                             propnewstrat = 0.75)

# Z 7
H$x_zone7_r25xy <- pz.2strat(zone1 = 7, # 25% of farms in Zone 2 using xy, rest using x
                             regstrat = "x",
                             newstrat = "xy",
                             propnewstrat = 0.25)
H$x_zone7_r50xy <- pz.2strat(zone1 = 7,
                             regstrat = "x",
                             newstrat = "xy",
                             propnewstrat = 0.5)
H$x_zone7_r75xy <- pz.2strat(zone1 = 7,
                             regstrat = "x",
                             newstrat = "xy",
                             propnewstrat = 0.75)

# Z 8
H$x_zone8_r25xy <- pz.2strat(zone1 = 8, # 25% of farms in Zone 2 using xy, rest using x
                             regstrat = "x",
                             newstrat = "xy",
                             propnewstrat = 0.25)
H$x_zone8_r50xy <- pz.2strat(zone1 = 8,
                             regstrat = "x",
                             newstrat = "xy",
                             propnewstrat = 0.5)
H$x_zone8_r75xy <- pz.2strat(zone1 = 8,
                             regstrat = "x",
                             newstrat = "xy",
                             propnewstrat = 0.75)

# Z 9
H$x_zone9_r25xy <- pz.2strat(zone1 = 9, # 25% of farms in Zone 2 using xy, rest using x
                             regstrat = "x",
                             newstrat = "xy",
                             propnewstrat = 0.25)
H$x_zone9_r50xy <- pz.2strat(zone1 = 9,
                             regstrat = "x",
                             newstrat = "xy",
                             propnewstrat = 0.5)
H$x_zone9_r75xy <- pz.2strat(zone1 = 9,
                             regstrat = "x",
                             newstrat = "xy",
                             propnewstrat = 0.75)

# Z 10
H$x_zone10_r25xy <- pz.2strat(zone1 = 10, # 25% of farms in Zone 2 using xy, rest using x
                             regstrat = "x",
                             newstrat = "xy",
                             propnewstrat = 0.25)
H$x_zone10_r50xy <- pz.2strat(zone1 = 10,
                             regstrat = "x",
                             newstrat = "xy",
                             propnewstrat = 0.5)
H$x_zone10_r75xy <- pz.2strat(zone1 = 10,
                             regstrat = "x",
                             newstrat = "xy",
                             propnewstrat = 0.75)

# Z 11
H$x_zone11_r25xy <- pz.2strat(zone1 = 11, # 25% of farms in Zone 2 using xy, rest using x
                              regstrat = "x",
                              newstrat = "xy",
                              propnewstrat = 0.25)
H$x_zone11_r50xy <- pz.2strat(zone1 = 11,
                              regstrat = "x",
                              newstrat = "xy",
                              propnewstrat = 0.5)
H$x_zone11_r75xy <- pz.2strat(zone1 = 11,
                              regstrat = "x",
                              newstrat = "xy",
                              propnewstrat = 0.75)

# Z 12
H$x_zone12_r25xy <- pz.2strat(zone1 = 12, # 25% of farms in Zone 2 using xy, rest using x
                              regstrat = "x",
                              newstrat = "xy",
                              propnewstrat = 0.25)
H$x_zone12_r50xy <- pz.2strat(zone1 = 12,
                              regstrat = "x",
                              newstrat = "xy",
                              propnewstrat = 0.5)
H$x_zone12_r75xy <- pz.2strat(zone1 = 12,
                              regstrat = "x",
                              newstrat = "xy",
                              propnewstrat = 0.75)


# Z 13
H$x_zone13_r25xy <- pz.2strat(zone1 = 13, # 25% of farms in Zone 2 using xy, rest using x
                              regstrat = "x",
                              newstrat = "xy",
                              propnewstrat = 0.25)
H$x_zone13_r50xy <- pz.2strat(zone1 = 13,
                              regstrat = "x",
                              newstrat = "xy",
                              propnewstrat = 0.5)
H$x_zone13_r75xy <- pz.2strat(zone1 = 13,
                              regstrat = "x",
                              newstrat = "xy",
                              propnewstrat = 0.75)



##############################################

### IMPORTANT! ###
H <- arrange(H, farmID) # Need to rearrange order of data back to order of farmID, so matches with original farms.csv data set!
write.csv(H, file='src/2loci/H.csv')

#################################################################


# Create Function for 3 strategies, with a random mix of 2 strategies in subset of farms
# Eg x & y are randomly distributed among the 25% most connected farms, the remaining 75% use xy
ordered_3strats <- function(mix1, # 1st strategy used in subset (e.g. smallest farms)
                            prop1, # Total proportion of farms using mix1
                            mix2, # 2nd strategy used in subset
                            prop2,# Total proportion of farms using mix2
                            rest1){ # strategy used by rest of farms
  
  # Create vector with
  c(
    # 1. randomised mix of first 2 strategies:
    c(sample(rep(c(mix1, # randomly distribute mix1 and mix2
                   mix2), 
                 times=c(ceiling(i*prop1), # according to their respective proportions
                         ceiling(i*prop2))))),
    # 2. remainder as third strategies
    c(rep(rest1, times=ceiling(i*(1-prop1-prop2))))
  )[1:i] # ensure vector is same length as i
  
  # 'mix' will be assigned to first rows of H, so need to ensure H is ordered correctly   
}

# Smallest farms using strategy 1 or strategy 2, rest using strategy
H$small_25xzORyz_z <- ordered_3strats("xz", 0.125, "yz", 0.125, "z") # 25% smallest using EITHER xz or xy, rest z
# BIGGEST farms using strategy 1 or 2 -> reorder!
H <- arrange(H, desc(biomass))
# Smallest farms using strategy 1 or strategy 2, rest using strategy
H$big_25xzORyz_z <- ordered_3strats("xz", 0.125, "yz", 0.125, "z") # 25% BIGGEST using EITHER xz or xy, rest z


H$leastD_25xzORyz_z <- ordered_3strats("xz", 0.125, "yz", 0.125, "z") # 25% smallest using EITHER xz or xy, rest z
# BIGGEST farms using strategy 1 or 2 -> reorder!
H <- arrange(H, desc(Dwcon))
# Smallest farms using strategy 1 or strategy 2, rest using strategy
H$mostD_25xzORyz_z <- ordered_3strats("xz", 0.125, "yz", 0.125, "z") # 25% BIGGEST using EITHER xz or xy, rest z



