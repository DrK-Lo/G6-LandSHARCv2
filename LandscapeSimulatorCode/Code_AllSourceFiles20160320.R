	
if (!("sp" %in% installed.packages())){install.packages("sp", dependencies=TRUE)}
if (!("spacetime" %in% installed.packages())){install.packages("spacetime", dependencies=TRUE)}
if (!("zoo" %in% installed.packages())){install.packages("zoo", dependencies=TRUE)}
if (!("xts" %in% installed.packages())){install.packages("xts", dependencies=TRUE)}
if (!("gstat" %in% installed.packages())){install.packages("gstat", dependencies=TRUE)}
if (!("fields" %in% installed.packages())){install.packages("fields", dependencies=TRUE)}
if (!("animation" %in% installed.packages())){install.packages("animation", dependencies=TRUE)}
if (!("spam" %in% installed.packages())){install.packages("spam", dependencies=TRUE)}
if (!("gdata" %in% installed.packages())){install.packages("gdata", dependencies=TRUE)}
if (!("gtools" %in% installed.packages())){install.packages("gtools", dependencies=TRUE)}
if (!("gmodels" %in% installed.packages())){install.packages("gmodels", dependencies=TRUE)}
if (!("gplots" %in% installed.packages())){install.packages("gplots", dependencies=TRUE)}
	
			library(sp)
				library(spacetime)
				library(zoo)
				library(xts)
				library(gstat)
				library(fields)
			library(spam)
			library(gdata)
				library(gtools) 
				library(gmodels)
				library(gplots)   
	

###################### DRIFT IN C ####################	
if(file.exists("src/drift.o")){system(paste("rm src/drift.o"))}
if(file.exists("src/drift.so")){system(paste("rm src/drift.so"))}	
system("R CMD SHLIB src/drift.c")

	dyn.load("src/drift.so")		
	drift_C<- function(Ne, P_freqMAT){
		out<- .C("drift", as.integer(Ne), as.double(P_freqMAT), as.integer(length(P_freqMAT)))
		return(out[[2]])
	}


###################### MUTATION ####################	
if(file.exists("src/mutation.o")){system(paste("rm src/mutation.o"))}
if(file.exists("src/mutation.so")){system(paste("rm src/mutation.so"))}	
system("R CMD SHLIB src/mutation.c")
	dyn.load("src/mutation.so")

	mutation_C <- function(u, P_freqMAT){
			l <- length(P_freqMAT)
			out<- .C("mutation", as.double(u), as.double(P_freqMAT), as.integer(l))
		return(out[[2]])
	}


###################### SELECTION ####################
if(file.exists("src/selection.o")){system(paste("rm src/selection.o"))}
if(file.exists("src/selection.so")){system(paste("rm src/selection.so"))}
system("R CMD SHLIB src/selection.c")
dyn.load("src/selection.so")

selection_C <- function(p_vect, s_vect){
		l <- length(p_vect)
		out <- .C("selection", as.double(p_vect), as.double(s_vect), as.integer(l))
		return(out[[1]])
}


###################### MIGRATION IN C ####################
if(file.exists("src/migration.o")){system(paste("rm src/migration.o"))}
if(file.exists("src/migration.so")){system(paste("rm src/migration.so"))}

	system("R CMD SHLIB src/migration.c")
	dyn.load("src/migration.so")

	migration_C <- function(P, xdim, ydim, MigMat){
		 MigMatDim <- dim(MigMat)[1]
		 Extend <- floor(MigMatDim/2)
		 pnew <- rep(0, (xdim+Extend*2)*(ydim+Extend*2))
		 MigVect <- as.vector(MigMat)
	
		 out <- .C("migration", as.double(P), as.integer(xdim), as.integer(ydim),as.double(pnew), as.double(MigVect), as.integer(MigMatDim), as.integer(Extend)) 
		
		 pnew <- matrix(out[[4]], (xdim+Extend*2), (ydim+Extend*2))
		
		 pnew2 <- as.vector(pnew[(Extend+1):(xdim+Extend), (Extend+1):(ydim+Extend)])
			
		# pnew3 <- pnew2/P_scaleVECT	
		 return(pnew2)
		 }


############### OTHER FUNCTIONS ##############################
ConvertEnviToIslands <- function(s_VECT, IslandSize, X_Demes, Y_Demes){
	x = seq(IslandSize,X_Demes, by=IslandSize)
	y = seq(IslandSize,Y_Demes, by=IslandSize)
	
	out <- matrix(NA, length(x), length(y))
	s_MAT <- matrix(s_VECT, X_Demes, Y_Demes)
	
	for (i in 1:length(x)){
		for (j in 1:length(y)){
			out[i,j] = mean(s_MAT[(x[i]-4):x[i], (y[j]-4): y[j]])
			
		}
	}
	return(as.vector(out))
}


ConvertVectToLS <- function(X_Demes, Y_Demes, P_Vect, IslandSize){
	k <- 0
	P_freqMAT <- matrix(NA, X_Demes,Y_Demes)
	X_Demes_Is <- X_Demes/IslandSize
	Y_Demes_Is <- Y_Demes/IslandSize
	for (i in 1:X_Demes_Is){ 
		for (j in 1:Y_Demes_Is){
			k <- k + 1
			P_freqMAT[((i-1)*IslandSize+1):(i*IslandSize), ((j-1)*IslandSize+1):(j*IslandSize)] <-  P_Vect[k]
			
		}
	}
	return(P_freqMAT)
}				

Get_selVECT <- function(s_low,s_high, Envi_VECT){
	if (s_high==0 & s_low==0){s_VECT=0}else{
		if(s_high==-s_low){
			s_slope <- (s_high - s_low)/(4*sd(Envi_VECT))
			#slope is determined by 4 SD (ie ~ 95% CI)
			#intercept is 0 when s_low=s_high
			s_VECT <- s_slope*Envi_VECT
			}else{
				s_slope <- (s_high - s_low)/(4*sd(Envi_VECT))
				s_int <- s_high-2*sd(Envi_VECT)*s_slope 
				s_VECT <- s_slope*Envi_VECT + s_int
			} #end nested if
	}		#end if
	return(s_VECT)
} #end function

getRunName <- function(){
	paste(round(unclass(Sys.time())), round(runif(1)*100000), sep="_")
}

plotAlleleFreqImage <- function(X_Locs, Y_Locs, P_plotMAT, gen){
    P_plotMAT=matrix(as.vector(P_plotMAT), nrow=length(X_Locs), ncol=length(Y_Locs))    
	image.plot(X_Locs, Y_Locs, P_plotMAT, main=
               paste("Allele freqs:", gen, "gens"), 
              breaks=seq(0,1,0.01), col=two.colors(n=100, start="darkblue", end="darkred", middle=grey(0.8), alpha=1.0),cex.main=0.8)	
}


GetColNamesEnvi <- function(){return(c("runname", "X_Dist", "Y_Dist", "X_Dist_Adj", "Y_Dist_Adj", "X_Demes","Y_Demes", "Gamma", "SD", "CellSize", "BETA_1", "BETA_X", "BETA_Y", "RANGE", "SILL", "REPID", "NumLociBase"))}

GetColNamesLandscapeMeta <- function(){return(c("MODEL_TYPE", "LOCUS_TYPE", "runname", "X_Demes", "Y_Demes", "Scale", "u", "GEN", "LandscapeMatGen", "FSTgen", "MOVIE", "MOVIEgen", "m", "IslandSize", "Island_N", "DispersalKernel_FilePath", "K_All", "GEN_Refugia", "r_growth", "Refugia_StartMat_FilePath", "AlleeDensity", "MOVIE_Density", "Refugia_Selection", "ENVI_ID", "ENVI_DIR", "s_low", "s_high",  "p_start", "p_LS_end", "He_LS_start", "He_LS_end", "FST_All", "Landscape_Corr", "TotalT"))}


GetLandscape <- function(DemogSubFolder, gen_foc, runname ){
		LandscapePath <- paste(DemogSubFolder,"/",runname,"Movie/",runname,"Matrix_AlleleFreqs_Gen",gen_foc, sep="")
		#NOTE: path may be different for a different type of demography (my bad)
		LandscapeMat <- as.matrix(read.table(LandscapePath, header=FALSE)) 
	return(LandscapeMat)
}

GetEnviMeta <- function(ENVI_DIR){
	EnviFile <- read.table(paste(ENVI_DIR, "All_Environments_Metadata.txt", sep=""))
	colnames(EnviFile)<- GetColNamesEnvi()
    BETA <- paste(EnviFile$BETA1, EnviFile$BETA_X, EnviFile$BETA_Y, sep="_")
    Envis <- cbind(EnviFile, BETA)
	return(Envis)
}

GetEnviInfo <- function(Index, Envi_Meta){
	return(c(as.character(Envi_Meta$runname[Index]), Envi_Meta$BETA_X[Index], Envi_Meta$BETA_Y[Index], Envi_Meta$RANGE[Index], Envi_Meta$REPID[Index]))
}

GetEnvi <-function(runname, headDir){
	file=paste(headDir, runname,"EnviMat.txt",sep="")
	Envi<-data.matrix(read.table(file,header=FALSE))
	colnames(Envi)=NULL
	return(t(Envi))
}

SetOrCreateDir <- function(Dir){
	if (file.exists(Dir)){
		setwd(file.path(Dir))
	} else {
		dir.create((Dir))
		setwd(file.path(Dir))
	}
}

#############################################################
######Tic Toc functions for timing runs
#######################################################
#############################################################
tic <- function(gcFirst = TRUE, type=c("elapsed", "user.self", "sys.self")){
   type <- match.arg(type)
   assign(".type", type, envir=baseenv())
   if(gcFirst) gc(FALSE)
   tic <- proc.time()[type]         
   assign(".tic", tic, envir=baseenv())
   invisible(tic)
}

toc <- function(){
   type <- get(".type", envir=baseenv())
   toc <- proc.time()[type]
   tic <- get(".tic", envir=baseenv())
   print(toc - tic)
   invisible(toc)
}

