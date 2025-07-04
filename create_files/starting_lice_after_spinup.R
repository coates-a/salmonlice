# Use the number of lice in the 'control' sim (all farms using mech only) *AFTER THE 2-YEAR SPIN-UP*
# as the starting numbers for other sims

# Load'allz' sim:
noge <- read_csv("outputs/popdyn_2loci_new_s_noGE.csv")

noge.y2 <- subset(noge, tt==2 & # year 2 (end up spin up)
                    geno=='SSUU') # can isolate 1 genotype since looking at total abundance per stage (abund.stage)
noge.y2 <- noge.y2[,c("farm","stage","abund.stage")] # 

y2 <- dcast(noge.y2, farm ~ stage, value.var = "abund.stage" )
colnames(y2) <- c("farm","ad.y2","ch.y2","la.y2","pa.y2")

farms <- merge(farms, y2, by="farm", all=T) # Add to farms dataset

write.csv(farms, 'src/farms.csv')

