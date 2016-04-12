## To Do:
## Input matrices for environment, starting population size, carrying capacity,
## unsuitable habitat, starting allele frequency chosen from beta
source("Code_AllSourceFiles20160320.R")
LandSHARC_IBD <- function(outDir, LocusType, 
                      InitialPopMat, 
                      pMat, # a value or matrix
                      Kmat, # a value or matrix
                      u=1e-06, r_growth, 
                      DispersalKernel, 
                      # need to check if p_start can be vector
                      Gen, GenFreqMatOut=GEN, GenPopMatOut=GEN,
                      MoviePfreq=FALSE, MoviePGenFreq = 100, 
                      AlleeDensity = 2, 
                      MoviePopSize = FALSE, MoviePopGenFreq = 100,
                      ENVI_matrix = NULL,  
                      s_high = NULL, s_low=-s_high, runname=NA){

#############################################################		
########### ARE INPUT FILES THERE ?????? #########################
##################################################################
 	if (is.null(outDir)==TRUE){print ("Missing the head directory for output: outDir")}
 	if (is.null(LocusType)==TRUE){print("Please enter locus type (NEUTRAL, SELECTION)")	}
 	if (is.null(InitialPopMat)==TRUE){print("Missing Initial population matrix: InitialPopMat")}			
 	if (is.null(u)==TRUE){print ("Missing mutation rate: u")}
 	if (is.null(pMat)==TRUE){print ("Missing starting allele freq: pMat")} 	
 	if (is.null(Gen)==TRUE){print ("Missing Gen")}	
# 	if (LocusType == "SELECTION") {
# 			if (is.null(ENVI_ID)==TRUE){print ("Missing environment file: ENVI_ID")}
# 			if (is.null(ENVI_DIR)==TRUE){print ("Missing directory for the ENVI file: ENVI_DIR")}
# 			if (is.null(s_low)==TRUE){print ("Missing s_low")}
# 			if (is.null(s_high)==TRUE){print ("Missing s_high")}
# 		}

	#########################################
	#Get name for this run
	#########################################
	if(is.na(runname)){runname <- getRunName()}
  if (file.exists(outDir)){
    setwd(outDir); print(c("Moving to existing output directory", getwd()))
  } else {
    dir.create(outDir); ; print(c("Creating output directory", getwd()))
  }

	 MovieDirname <- paste(outDir,"/", runname,"Movie",sep="")
    dir.create(MovieDirname)
		
	#########################################
	####If LocusType=="SELECTION", Input Environment
	#########################################	
if (LocusType=="SELECTION"){
		ENVI_ID <- as.character(ENVI_ID)
		Envi <- GetEnvi(ENVI_ID, ENVI_DIR)
		EnviFile <- read.table(paste(ENVI_DIR, "All_Environments_Metadata.txt", sep=""))
		colnames(EnviFile) <- GetColNamesEnvi()
		Index<- which(EnviFile$runname==ENVI_ID)
		ENVI_Info <- c(Index, EnviFile$BETA_X[Index],  EnviFile$BETA_Y[Index], EnviFile$RANGE[Index])
	
		Envi_VECT = as.vector(Envi)
	
		if(dim(Envi)[1]!=X_Demes){print("Environment Matrix not same size as Landscape (x)")}
		if(dim(Envi)[2]!=Y_Demes){print("Environment Matrix not same size as Landscape (y)")}
		
		#Input relationship between selection and environment
		#######################################
		s_VECT <- Get_selVECT(s_low, s_high, Envi_VECT)	
		min_s <- min(s_VECT)
		max_s <- max(s_VECT)
		
		if(length(s_VECT)!=(X_Demes*Y_Demes)){print("Environment not same size as landscape"); break}
		
		#Write selection matrix to file
		#########################################	
		setwd(MovieDirname)
		write(as.vector(s_VECT), file = paste(runname, "_SelectionMatrix", sep=""), ncol=Y_Demes)
		
		pdf(file=paste(runname,"_SelectionLandscape_Envi",ENVI_ID,".pdf",sep=""), width=8, height=3, bg="white")
		par(mfrow=c(1,2), mar=c(2,2,2,6), cex.axis=0.5)
		image.plot(X_Locs,Y_Locs, Envi, main=paste(ENVI_ID), cex.main=0.7)	
		image.plot(X_Locs,Y_Locs, matrix(s_VECT, ncol=Y_Demes), main=paste( "s_min=",round(min_s,5), " ,s_max=", round(max_s,5),  sep=""), cex.main=0.7, col=two.colors(n=200, start="blue", end="red", middle="white", alpha=1.0))	
		dev.off()

 	}# end if selection
		
	#########################################
	#Print Run Metadata
	#########################################	
  X_Demes <- nrow(InitialPopMat)
  Y_Demes <- ncol(InitialPopMat)
 		Metadatafile <- paste(outDir,"/",runname,"Metadata",sep="")
 		write(runname, file= Metadatafile) 
 		write(paste("Date=",date()), file= Metadatafile, append=TRUE)
 		write(paste("LocusType=",LocusType), file= Metadatafile, append=TRUE)
 		write(paste("X_Demes=",X_Demes), file= Metadatafile, append=TRUE)
 		write(paste("Y_Demes=",Y_Demes), file= Metadatafile, append=TRUE)
    write(paste("InitialPopMat=", paste(InitialPopMat, collapse=",")), file= Metadatafile, append=TRUE)
 		write(paste("u=",u), file= Metadatafile, append=TRUE)			
 		write(paste("pMat=",paste(pMat, collapse=",")), file= Metadatafile, append=TRUE)
 		write(paste("Gen=",gen), file= Metadatafile, append=TRUE)
 		write(paste("LandscapeMatGen=",LandscapeMatGen), file= Metadatafile, append=TRUE)
 		write(paste("DispersalKernel=", paste(DispersalKernel, collapse="")), file= Metadatafile, append=TRUE)
 		write(paste("Kmat=",paste(Kmat, collapse=",")), file= Metadatafile, append=TRUE)
 		write(paste("r_growth=",r_growth), file= Metadatafile, append=TRUE)
   	write(paste("AlleeDensity=",AlleeDensity), file= Metadatafile, append=TRUE)

  	write(paste("MoviePfreq=",MoviePfreq), file= Metadatafile, append=TRUE)
    write(paste("MoviePGenFreq=",MoviePGenFreq), file=Metadatafile, append=TRUE)
    write(paste("MoviePopSize=",MoviePopSize), file=Metadatafile, append=TRUE)
    write(paste("MoviePopGenFreq=",MoviePopGenFreq), file=Metadatafile, append=TRUE)
# 		write(paste("ENVI_ID=",ENVI_ID), file= Metadatafile, append=TRUE)
# 		write(paste("ENVI_DIR=",ENVI_DIR), file= Metadatafile, append=TRUE)
# 		write(paste("s_low=",s_low), file= Metadatafile, append=TRUE)
# 		write(paste("s_high=",s_high), file= Metadatafile, append=TRUE)
# 		
#############################################################			
	#########################################
	#Setup starting population and allele frequencies
	#########################################			
#############################################################	
	
	### ISLAND MODEL
	##############################
# 	if (MODEL_TYPE=="ISLAND"){
# 		TotIslandSize <- IslandSize*IslandSize
# 		NumIslands <- (X_Demes*Y_Demes)/(IslandSize*IslandSize)	
# 		N_vect <- rep(Island_N, NumIslands)
# 		vectorsize <- NumIslands	
# 		OccupiedVect <- rep(1, vectorsize)
# 		P_freqMAT <- rep(p_start, vectorsize)
# 		DensityPerKM <- round(sum(N_vect)/(X_Dist_Adj*Y_Dist_Adj), 3)
# 		
# 		#Need to convert selection-vector to same size as island vector
# 		if (LocusType=="SELECTION"){
# 		s_VECT <- ConvertEnviToIslands(s_VECT, IslandSize, X_Demes, Y_Demes)
# 		Envi_VECT <- ConvertEnviToIslands(Envi_VECT, IslandSize, X_Demes, Y_Demes)
# 		}
# 	}

		vectorsize<- X_Demes*Y_Demes
		MigMat <- as.matrix(read.table(DispersalKernel_FilePath))
		MigVect <- as.vector(MigMat)
		K_Vect <- rep(K_All, vectorsize)

			OccupiedVect <- InitialPopVect
			OccupiedVect[which(InitialPopVect>0)]=1
			N_vect <- InitialPopVect
			P_scaleVECT <- migration_C(OccupiedVect, X_Demes, Y_Demes, MigMat)
        # this is doing some thing weird
			P_freqMAT <- rep(p_start, vectorsize)*OccupiedVect
  	PopVect2<-InitialPopVect
		AtEQ <- FALSE
		EqGen<-NA
		
#############################################################		
	#########################################
	# Loop over generations and let evolution happen
	##########################################
	
		for (gen in 1:GEN){
			if (sum(is.na(PopVect2))>0){print("Error: NAs in Abundance vector"); break;}

			if (AtEQ==FALSE){
				Last_PopVect <- PopVect2
				Last_Occupied <- Last_PopVect
				Last_Occupied[Last_PopVect>0]=1
			
				Fecundity_Vect <- PopVect2 + r_growth*PopVect2*(K_All-PopVect2)/K_All #logistic
								
				PopVect <-migration_C(Fecundity_Vect, X_Demes, Y_Demes, MigMat)  #this is the next generation
				PopVect[PopVect<AlleeDensity]=0	
				if(sum(PopVect)==0){print("Sorry you went extinct!"); break;}
				
				OccupiedVect <- PopVect
				OccupiedVect[PopVect>0]=1	 #These cells become occupied in this time step
					
				if(gen>100 & sum(Last_PopVect==round(PopVect))==(X_Demes*Y_Demes) & AtEQ == FALSE){AtEQ=TRUE; EqGen=gen; print(c("Eq has been reached at gen =", gen))} #if at equilibrium, tell it
			}#end if AtEQ 
			
			if (AtEQ==TRUE & gen==(EqGen+1)){
			Fecundity_Vect <- PopVect + r_growth*PopVect*(K_Vect-PopVect)/K_Vect
			PopVect <- migration_C(Fecundity_Vect, X_Demes, Y_Demes, MigMat)
			}#end if
		
			PopVect2 <- round(PopVect)			# round pop Vect

		############################################################
		#Step 1: Mutation occurs on reproducing demes
		############################################################
		mutVECT <- mutation_C(u,P_freqMAT)*Last_Occupied

		############################################################
		#Step 2: Dispersal/migration into new demes
		############################################################
	# In this case, we calculate the number of offpsring with allele 1 that migrate and where they go (mutVECT*Fecundity_Vect)
	# And then we scale by the new population size in newly occupied demes (2nd line of code)
	#don't use the rounded PopVect2 here, as it will not account correctly for allele type
		
		migVECT<-migration_C(mutVECT*Fecundity_Vect, X_Demes, Y_Demes, MigMat)*OccupiedVect
		migVECT[OccupiedVect==0] = 0	# #this step removes alleles from demes with no individuals.  multiplying doesn't work.
		migVECT[OccupiedVect==1] = migVECT[OccupiedVect==1]/PopVect[OccupiedVect==1]
			#this step scales allele freqs by deme size
				# the scaling must be done with the "fractional" population matrix, if it is rounded off you can get allele freqs > 1
		
		############################################################
		#Step 3: Selection on offspring
		############################################################
		if(LocusType=="SELECTION"){
				selVECT <- selection_C( migVECT, s_VECT)*OccupiedVect
			}
		
		if(LocusType=="NEUTRAL"){
			selVECT <- migVECT
		}
			
		############################################################
		#Step 4: Drift (ie population regulation)
		############################################################
		 driftVECT <- drift_C(PopVect2, selVECT)*OccupiedVect
		
		############################################################
		#Step 5: Reassign allele freqs to next generation
		############################################################	
		P_freqMAT <- driftVECT	

		############################################################
		#Write desired outputs for Refugia model
		############################################################

			#Write pdf Image of allele frequencies to a file
			##################################	
			if (gen%%MOVIEgen==0 & MOVIE==TRUE){
				 setwd(MovieDirname)
				
				 P_plot = migVECT
		          P_plot[OccupiedVect==0]=NA
		          P_plotMAT <- matrix(P_plot, X_Demes, Y_Demes)
				pdf(file=paste(runname,"ImageAfterMig_",gen, ".pdf",sep=""), width=4, height=3, bg="white")
	             par(mar=c(3,3,1,1))
	             plotAlleleFreqImage(X_Locs, Y_Locs, P_plotMAT, gen, FST_All)	#need to edit this function
	             dev.off()			
			}	
			
			#Write pdfs Image of density to a file for refugia model
			##################################	
			if (MOVIE_Density==TRUE & gen%%MOVIEgen==0){
				Density <- PopVect2
				Density[OccupiedVect==0]=NA
				Density_plot <- matrix(Density, X_Demes, Y_Demes)
			pdf(file=paste(runname,"ImageDensity_",gen, ".pdf",sep=""), width=4, height=3 ,bg="white");  par(mar=c(3,3,1,1))
			image.plot(Density_plot, x=X_Locs, y=Y_Locs,  main=paste(runname, "Density at gen", gen),cex.main=0.8, xlab="", ylab="")
			dev.off()	
			}
			
			#Write matrices of landscape to a file
			##################################		
			if (gen%%LandscapeMatGen==0){
				setwd(MovieDirname) 
	            write(as.vector(migVECT), file=paste(runname,"Matrix_AlleleFreqs_Gen", gen,sep=""), ncolumns=Y_Demes)
			}
				
		}#end output if Island or IBD

} # end loop through gen

return()
  # return p freq, pop size, 


# ############################################################				
# ######################## Write Ouputs ####################################
# ############################################################			
# 	#################################
# 	#Write time elapsed to metadata
# 	#################################	
# 		setwd(MovieDirname)
# 			write(" ", file=Metadatafile, append=TRUE)
# 			write(paste("Number generations elapsed =", gen), file=Metadatafile, append=TRUE)
# 
# 			p_LS_end <- mean(migVECT)
# 			He_LS_start	<- 	p_start*(1-p_start)	
# 			He_LS_end <- mean(migVECT*(1-migVECT))
# 
# 			write(paste("p_LS_end=",p_LS_end), file= Metadatafile, append=TRUE)
# 			write(paste("He_LS_start=",p_LS_end), file= Metadatafile, append=TRUE)
# 			write(paste("He_LS_end=",p_LS_end), file= Metadatafile, append=TRUE)
# 			endT <- toc()
# 			TotalT <- round(endT-startT)
# 			write(paste("Time Elapsed (sec)",TotalT), file=Metadatafile, append=TRUE)
# 
# 	print(c(runname, p_start))
# 
# 	#################################
# 	#Write main output to a file in main directory
# 	#################################	
# 		setwd(headDir)	
# 		allVars <- c(MODEL_TYPE, LocusType, runname, X_Demes, Y_Demes, Scale, u, GEN, LandscapeMatGen, FSTgen, MOVIE, MOVIEgen, m, IslandSize, Island_N, DispersalKernel_FilePath, K_All, GEN_Refugia, r_growth, Refugia_StartMat_FilePath, AlleeDensity, MOVIE_Density, Refugia_Selection, ENVI_ID, ENVI_DIR, s_low, s_high,  p_start, p_LS_end, He_LS_start, He_LS_end, FST_All, Landscape_Corr, TotalT, "\n", sep=" ")
# 		write(allVars, FileName_LandscapeMetadata , append=TRUE, ncolumns=length(allVars))

}#end function
