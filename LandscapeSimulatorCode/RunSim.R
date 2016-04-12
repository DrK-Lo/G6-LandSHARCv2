setwd("/Users/katie/Desktop/CurrResearch/LandSHARCv2/LandscapeSimulatorCode/")

outDir = "output" 
LocusType = "NEUTRAL"                       
InitialPopMat = read.table("") 
pMat = 0.2 
Kmat = 1000
u=1e-06, 
r_growth = 1 
DispersalKernel = as.matrix(read.table(""))
Gen = 100
GenFreqMatOut = GEN 
GenPopMatOut = GEN
MoviePfreq=FALSE 
MoviePGenFreq = 100
AlleeDensity = 2
MoviePopSize = FALSE 
MoviePopGenFreq = 100
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
