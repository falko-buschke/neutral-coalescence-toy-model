setwd("SET WORKING DIRECTORY")   # Set working directory    

# Open bird community dataset
# The data set are available from: https://github.com/falko-buschke/Madagascar-neutral-model
birds <-  read.table("Mada_endemics.txt",header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
birds <- birds[,-1] #remove column with species names
richness <- colSums(birds)


# Open grid information dataset
mad.grid <-   read.table("Mada_grid.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
coord <- cbind(mad.grid$Long,mad.grid$Lat)
energy <- mad.grid$NDVI


# create a distance matrix (measures distance in kilometers)
distance.1 <- as.matrix(dist(coord, method = "euclidean", diag = TRUE, upper = TRUE))
distance <- distance.1/50  # convert to "grid units"

para.a <- 0.3      # parameter for short-distance dispersal
para.b <- 0.2      # parameter for long-distance dispersal
K.par <- 350       # parameter for habitat caacity

# calculate the dispersal kernel
dispersal <- para.a^(distance) + (para.b^2/(distance^2 + para.b^2))
diag(dispersal) <- 1



conversion <- length(energy)* (energy/sum(energy))     # this number is the correction factor equivalent to Muneerpeerakul et al. (s2008)

S <- dim(birds)[1]      # number of species in metacommunity
max.rich <- S           # maxmum number of species in the coalescence model
M <- dim(birds)[2]      # number of sites
K <- ceiling(K.par * conversion)     # number of individuals per sample (average K is 700)


species <-  (seq(1:S))    # just a vector of species
samples <- seq(1:M)       # a vector of sites names
ccapacity <- K    

# add disrectional dispersal by making dispersal more likely from a quadrat with high habitat capacity
# to a quadrat with low habitat capacity
imm.null <- t( apply(dispersal,1,function(dispersal) dispersal*ccapacity))
denom <- colSums(imm.null)
# this is the dispersal kernel to be used in the simulations
disp <- t( apply(imm.null,1,function(imm.null) imm.null/denom))

# Creat a vector of the site identification for each bird unit.
sites <- rep(1:M,K)

J <- sum(K) #  total number of bird units

sp.name <- 1:J  # Start by assigning each bird unit a unique ID

# This is a probability vector that determines whether a bird unit has already coalesed. 
ind.prob <- rep(1,J)

# This is simply to set the counter to display the number of time-steps (starts at 1)
x <- 1


# Start the loop and continue until the number of species in the simulation matches the observed species richness
while(length(unique(sp.name))> max.rich) {

# Choose a bird unit at random
ind.id <- sample(ind,1,prob=ind.prob)
# Identify the site location of selected bird unit
site.id <- sites[ind.id]  
# Idntify the source of selected bird unit based on the dispersal kernel
source.id <- sample(1:M,1,prob=disp[site.id,])

# Identify a subset of bird units that occur in the source quadrat
migr.id <- ind[which(sites==source.id)]

# Select a subset of bird units that DO NOT already share the same species ID
imm.id <- migr.id[which(sp.name[migr.id]!=sp.name[ind.id])]

# Select an parent/immigrant from the source quadrat
sist.sp <- sample(imm.id,1)

# Assign the same species ID as the parent bird unit to the randomly selected focal bird unit 
# AND all other bird units that already share that species ID
sp.name[which(sp.name==ind.id)] <- sp.name[sist.sp]

# Since the focal bird unit already has an assigned species ID, it is excluded from the model from now on.
ind.prob[ind.id] <- 0

# Report the progress of the simulation 
x <- x+1
if (x/1000 == floor(x/1000))
{print(x)
flush.console()
} 
}

# Re-format data into a site by species matrix
Patches <- as.matrix(table(sp.name,sites))
