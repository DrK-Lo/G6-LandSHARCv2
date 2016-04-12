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
                      Gen, GenPFreqMatOut=Gen, GenPopMatOut=Gen,
                      MoviePfreq=FALSE, MoviePopSize = FALSE, MovieGenFreq = 100, 
                      AlleeDensity = 2, 
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
  print(runname)
  if (file.exists(outDir)){
    setwd(outDir); print(c("Moving to existing output directory", getwd()))
  } else {
    dir.create(outDir); setwd(outDir); print(c("Creating output directory", getwd()))
  }

	 MovieDirname <- paste(runname,"Movie",sep="")
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
		#setwd(MovieDirname)
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
 		Metadatafile <- paste(runname,"Metadata",sep="")
 		write(runname, file= Metadatafile) 
 		write(paste("Date=",date()), file= Metadatafile, append=TRUE)
 		write(paste("LocusType=",LocusType), file= Metadatafile, append=TRUE)
 		write(paste("X_Demes=",X_Demes), file= Metadatafile, append=TRUE)
 		write(paste("Y_Demes=",Y_Demes), file= Metadatafile, append=TRUE)
    write(paste("InitialPopMat=", paste(InitialPopMat, collapse=",")), file= Metadatafile, append=TRUE)
 		write(paste("u=",u), file= Metadatafile, append=TRUE)			
 		write(paste("pMat=",paste(pMat, collapse=",")), file= Metadatafile, append=TRUE)
 		write(paste("Gen=",Gen), file= Metadatafile, append=TRUE)
 		write(paste("DispersalKernel=", paste(DispersalKernel, collapse="")), file= Metadatafile, append=TRUE)
 		write(paste("Kmat=",paste(Kmat, collapse=",")), file= Metadatafile, append=TRUE)
 		write(paste("r_growth=",r_growth), file= Metadatafile, append=TRUE)
   	write(paste("AlleeDensity=",AlleeDensity), file= Metadatafile, append=TRUE)

  	write(paste("MoviePfreq=",MoviePfreq), file= Metadatafile, append=TRUE)
    write(paste("MovieGenFreq=",MovieGenFreq), file=Metadatafile, append=TRUE)
    write(paste("MoviePopSize=",MoviePopSize), file=Metadatafile, append=TRUE)
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
		vectorsize<- X_Demes*Y_Demes
    MigMat <- as.matrix(DispersalKernel)
    if(length(Kmat)==1){K_Vect <- rep(Kmat, vectorsize)}else{
      K_Vect <- as.vector(Kmat) # careful here, could be subject to recycling
    }
    InitialPopVect <- as.vector(InitialPopMat)
		OccupiedVect <- InitialPopVect
		OccupiedVect[which(InitialPopVect>0)]=1
    N_vect <- InitialPopVect
		  # P_scaleVECT <- migration_C(OccupiedVect, X_Demes, Y_Demes, MigVect)
        # this is doing some thing weird
    PopVect2<-InitialPopVect
		AtEQ <- FALSE
		EqGen<-NA
    if(length(pMat)==1){P_freqMAT <- rep(pMat, vectorsize)*OccupiedVect}else{
      P_freqMAT <- pMat
    }
		
#############################################################		
	#########################################
	# Loop over generations and let evolution happen
	##########################################
	
  for (gen in 1:Gen){
			if (sum(is.na(PopVect2))>0){print("Error: NAs in Abundance vector"); break;}

			if (AtEQ==FALSE){
				Last_PopVect <- PopVect2
				Last_Occupied <- Last_PopVect
				Last_Occupied[Last_PopVect>0]=1
			
        growth <- (r_growth*PopVect2*(K_Vect-PopVect2)/K_Vect)
        growth[growth<0 & K_Vect==0] <- 0
          # here the infinite values indicate 0 carrying capacity
				  # any adults in those demes have 0 reproductive value
        Fecundity_Vect <- PopVect2 +  growth#logistic
				Fecundity_Vect[K_Vect==0] <- 0
        
				PopVect <-migration_C(Fecundity_Vect, X_Demes, Y_Demes, MigMat)  #this is the next generation
				PopVect[PopVect<AlleeDensity | K_Vect==0]=0	
          # any juveniles that migrate to patches with 0 carrying capacity bite the dust
				if(sum(PopVect)==0){print("Sorry you went extinct!"); break;}
				
				OccupiedVect <- PopVect
				OccupiedVect[PopVect>0]=1	 #These cells become occupied in this time step
					
				if(gen>100 & sum(Last_PopVect==round(PopVect))==(X_Demes*Y_Demes) & AtEQ == FALSE){AtEQ=TRUE; EqGen=gen; print(c("Demographic Eq has been reached at gen =", gen))} #if at equilibrium, tell it
			}#end if AtEQ 
		
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
			if (gen%%MovieGenFreq==0 & MoviePfreq==TRUE){
				      P_plot = migVECT
		          P_plot[OccupiedVect==0]=NA
		          P_plotMAT <- matrix(P_plot, X_Demes, Y_Demes)
				png(file=paste(MovieDirname,"/",runname,"ImageAfterMig_",gen, ".png",sep=""), 
            width=4, height=3, bg="white", units="in",res=200)
	             par(mar=c(3,3,1,1))
	             plotAlleleFreqImage(1:X_Demes, 1:Y_Demes, P_plotMAT, gen)	#need to edit this function
	             dev.off()			
			}	
			
			#Write pdfs Image of density to a file for refugia model
			##################################	
			if (MoviePopSize ==TRUE & gen%%MovieGenFreq==0){
				Density <- PopVect2
				Density[OccupiedVect==0]=NA
				Density_plot <- matrix(Density, X_Demes, Y_Demes)
			png(file=paste(MovieDirname,"/",runname,"ImageDensity_",gen, ".png",sep=""), 
          width=4, height=3, units="in",res=200 ,bg="white");  
      par(mar=c(3,3,1,1))
			image.plot(Density_plot, x=1:X_Demes, y=1:Y_Demes,  main=paste(runname, "Density at gen", gen),cex.main=0.8, xlab="", ylab="")
			dev.off()	
			}
			
			#Write matrices of landscape to a file
			##################################		
    		if (gen %in% GenPFreqMatOut){
           write(as.vector(migVECT), file=paste(runname,"Matrix_AlleleFreqs_Gen_", gen,sep=""), ncolumns=Y_Demes)
    		}
        if (gen %in% GenPopMatOut){
          write(as.vector(PopVect2), file=paste(runname,"Matrix_PopSize_Gen_", gen,sep=""), ncolumns=Y_Demes)
        }
  #if(gen%%10==0){print(gen)}
  #print(gen)
  } # end loop through gen
  setwd("..")
}#end function
