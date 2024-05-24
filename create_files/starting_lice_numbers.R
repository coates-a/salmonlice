### Estimate starting lice infestations ####
install.packages("dplyr")
library(dplyr)
library(readr)

### FROM 2012 - 2022 DATA:

lice <- read_csv("~/Evolution model/dseed/src/2loci/create_files/lice_2012_2022.csv") 

# Subset only farms used in simulation
farms <- read.csv('src/2loci/farms.week1.csv', header = TRUE) # Add on to existing farm dataset
lice <- lice %>%
  filter(farmID %in% farms$farmID)

#### Explore infestation rates through time ####
lice$wk <- (lice$week-1)/52
lice$tt <- lice$year+lice$wk
lice.yrs <- (aggregate(cbind(female,motile,sessile) ~ tt, mean, data=lice)) 
ggplot(lice.yrs, aes(x=tt, y=female))+
  geom_line()+
  scale_x_continuous(breaks=seq(2012,2022,by=1))

# Much higher lice levels & more variation 2012-2014
# Use data from 2015 onwards (more relevant to today)

lice <- subset(lice, year>=2015)
# Since simulations begin at Week 1:
lice.wk1 <- subset(lice, week==1)
# Mean Week 1 infestations over 2015-2022 for each farm
lice.wk1 <- (aggregate(cbind(female,motile,sessile)~farmID, mean, data=lice.wk1)) 

farms <- merge(farms, lice.wk1, by="farmID", all=TRUE) # Merge datasets
farms[is.na(farms)] <- 0  # Replace NAs with 0 (assume no infestation on those farms)

write.csv(farms, file = 'src/2loci/farms.week1.csv', row.names = F) # Save as CSV
