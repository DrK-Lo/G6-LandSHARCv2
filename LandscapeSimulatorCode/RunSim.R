setwd("LandscapeSimulatorCode/")

# LandSHARC simulates one locus with two alleles

outDir = "output" 
LocusType = "NEUTRAL"                       
InitialPopMat = read.table("Demographies/northern_refuge") 
pMat = 0.2 # initial allele frequency on landscape, or a matrix of allele frequencies
Kmat = 1000 # carrying capacity of each deme, or a matrix of carrying capacity of each cell
u=1e-06 # mutation rate
r_growth = 0.5 # population growth rate
DispersalKernel = as.matrix(read.table("MigMat.txt"))
Gen = 1000 # number of generations
GenFreqMatOut = seq(100,1000,by=100) #a vector of generations when to output the allele frequency landscape 
GenPopMatOut = c(10,20, 50, 75, seq(100,1000,by=100)) #a vector of generations when to output the allele frequency landscape 
MoviePfreq= TRUE # output a movie of allele frequencies?
MoviePGenFreq = 100 # how often to output frame to the movie, every X generations
AlleeDensity = 2 # below this density, deme goes extinct
MoviePopSize = FALSE # output a movie of population size?
MoviePopGenFreq = 100 # # how often to output frame to the movie, every X generations
ENVI_matrix = NA  
s_high = NA
s_low=-s_high

LandSHARC_IBD(outDir, LocusType, InitialPopMat, pMat, # a value or matrix
                      Kmat, u=1e-06, r_growth, 
                      DispersalKernel, Gen, GenFreqMatOut=GEN, GenPopMatOut=GEN,
                      MoviePfreq=FALSE, MoviePGenFreq = 100, 
                      AlleeDensity = 2, 
                      MoviePopSize = FALSE, MoviePopGenFreq = 100,
                      ENVI_matrix = NULL,  
                      s_high = NULL, s_low=-s_high)
