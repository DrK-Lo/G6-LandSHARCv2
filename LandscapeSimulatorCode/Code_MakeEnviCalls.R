#Katie Lotterhos
#July 2012
#Make several environments for a landscape

CodeDir <- "~/Dropbox/MetapopWestgrid/MetapopCodeWestgrid/"
setwd(CodeDir)
source("Code_AllSourceFiles.R")

headDir <- paste(CodeDir, "Environments", sep="")#"~/Dropbox/MetapopCode/Environments"

X_Dist <- 270
Y_Dist <- 270
SD <- 1

BETA <- matrix( c(1, 0, 0,
				1, 0.01, 0,
				1, 0, 0.01), ncol=3, byrow=TRUE)
#spatial trend in environmental data 
											#ie (1,0.1,0) has a trend in x-direction
											#strength of trend will be affected by RANGE

RANGE <- c(10, 30, 60, 120,  200)			#the distance at which covariance 
											#in the environment goes to 0


#BETA <- matrix(c(1, 0.1,0, 1, 0, 0.1), ncol=3, byrow=TRUE)
#RANGE <- 60

SILL <- 1									#scale of variance in the environment
REPS <- 1:10

num_b <- nrow(BETA)


#b <- 4
#r <- 11
	for (b in 1:num_b){
		for (r in 1: length(RANGE)){
			for (i in REPS){
			MakeEnvi(X_Dist, Y_Dist, SD, BETA[b,], RANGE[r], SILL, i, headDir)
		}
	}
}