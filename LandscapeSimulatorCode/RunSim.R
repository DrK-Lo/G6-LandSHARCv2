
# LandSHARC simulates one locus with two alleles
setwd("/Users/katie/Desktop/CurrResearch/G6-LandSHARCv2/LandscapeSimulatorCode")
source("0Code_ALL_Select_Cv2.R")

## Example 1: northern ref
#InitialPopMat = as.matrix(read.table("Demographies/northern_refuge"))
#Kmat = 50 # carrying capacity of each deme, or a matrix of carrying capacity of each cell

## Example 2: hole of unsuitable habitat in middle of landscape
InitialPopMat <- matrix(100, 500, 500)
Kmat <- matrix(10, dim(InitialPopMat)[1], dim(InitialPopMat)[2])
dim(Kmat)
Kmat[300:400, 300:400] <- 0

outDir = "output" #will use or create this directory relative to wd
LocusType = "NEUTRAL"                       
pMat = 0.5 # initial allele frequency on landscape, or a matrix of allele frequencies
u=1e-06 # mutation rate
r_growth = 0.5 # population growth rate
DispersalKernel = as.matrix(read.table("MigMat.txt"))
Gen = 1000 # number of generations
GenPFreqMatOut = Gen #c(10,20, 50, 75, seq(100,1000,by=100)) #a vector of generations when to output the allele frequency landscape 
GenPopMatOut = Gen #c(10,20, 50, 75, seq(100,1000,by=100)) #a vector of generations when to output the allele frequency landscape 
AlleeDensity = 2 # below this density, deme goes extinct
MoviePopSize = TRUE # output a movie of population size?
MoviePfreq= TRUE # output a movie of allele frequencies?
MovieGenFreq = 50 # # how often to output frame to the movie, every X generations
ENVI_matrix = NA  
s_high = NA
s_low=-s_high

setwd("/Users/katie/Desktop/CurrResearch/G6-LandSHARCv2/LandscapeSimulatorCode")
LandSHARC_IBD(outDir=outDir, LocusType=LocusType, 
              InitialPopMat=InitialPopMat, pMat=pMat, # a value or matrix
                      Kmat, u=1e-06, r_growth=r_growth, 
                      DispersalKernel=DispersalKernel, Gen=Gen, 
                    GenPFreqMatOut=GenPFreqMatOut, 
                    GenPopMatOut=GenPopMatOut,
                      AlleeDensity = 2, 
                     MoviePfreq=TRUE,  MoviePopSize = TRUE, MovieGenFreq = MovieGenFreq,
                      ENVI_matrix = NULL,  
                      s_high = NULL, s_low=-s_high)
