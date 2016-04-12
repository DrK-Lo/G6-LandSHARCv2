setwd("/Users/katie/Desktop/CurrResearch/LandSHARCv2/LandscapeSimulatorCode/")
X_Demes <- 10
Y_Demes <- 10
envimat <- matrix(1:10,X_Demes,Y_Demes)
Kmat <- matrix(1000, X_Demes, Y_Demes)
Kmat[5,5] <- NA
Kmat
#source("Code_AllSourceFiles20160320.R")
p_start=0.4
u=1e-06
LOCUS_TYPE = "NEUTRAL" #"NEUTRAL", "SELECTION"

MODEL_TYPE= "Refugia" # "ISLAND", "IBD", "Refugia"
    # if ISLAND, specify 
        # m
    # if IBD, specify
        # MigMat
    # if Refugia, specify
      Refugia_Selection = FALSE
      GEN <- 100
      MigMat <- as.matrix(read.table("MigMat.txt", colClasses = "numeric"))
      InitialPopMat <- matrix(0, X_demes, Y_demes)
      InitialPopMat[1,1] <- InitialPopMat[1,3] <- InitialPopMat[3,1] <-1000
      class(MigMat)
      InitialPopMat
      r_growth = 0.5
      AlleeDensity <- 2


s_high=0.1
s_low=-0.1
Envi_VECT = as.vector(envimat)
  	s_VECT <- Get_selVECT(s_low, s_high, Envi_VECT)	
		min_s <- min(s_VECT)
		max_s <- max(s_VECT)
vectorsize <- length(Envi_VECT)	
MigVect <- as.vector(MigMat)
InitialPopVect <- as.vector(InitialPopMat)
N_vect <- InitialPopVect
OccupiedVect <- rep(0, length(InitialPopVect))
OccupiedVect[which(InitialPopVect>0)]=1
OccupiedVectNA <- OccupiedVect
OccupiedVectNA[OccupiedVect==0]<-NA
KVect <- as.vector(Kmat)


  ### REFUGIA MODEL SPECIFICS
	####################################			
	if (MODEL_TYPE=="Refugia"){ ?????
		P_freqMAT <- rep(p_start, vectorsize)*OccupiedVect
    P_scaleVECT <- migration_C(OccupiedVect, X_Demes, Y_Demes, MigMat)
	}

##### LOOP THROUGH GENERATIONS

if (MODEL_TYPE=="Refugia"){
		PopVect2<-InitialPopVect
		AtEQ <- FALSE
		EqGen<-NA
	
		for (gen in 1:GEN){
			if (sum(is.na(PopVect2))>0){print("Error: NAs in Abundance vector"); break;}

			if (AtEQ==FALSE){
				Last_PopVect <- PopVect2
				Last_Occupied <- Last_PopVect
				Last_Occupied[Last_PopVect>0]=1
				
				
				Fecundity_Vect <- PopVect2 + r_growth*PopVect2*(KVect-PopVect2)/KVect #logistic
				Fecundity_Vect[is.na(Fecundity_Vect)] <- 0
				PopVect <- migration_C(Fecundity_Vect, X_Demes, Y_Demes, MigMat)  #this is the next generation
				PopVect[PopVect<AlleeDensity]=0	
        PopVect[is.na(KVect)]=0 #NAs in KVect are for unsuitable habitat
				if(sum(PopVect)==0){print("Sorry you went extinct!"); break;}
				
        
				OccupiedVect <- rep(0, length(PopVect))
				OccupiedVect[PopVect>0]=1	 #These cells become occupied in this time step

				if(gen>100 & sum(Last_PopVect==round(PopVect))==(X_Demes*Y_Demes) & AtEQ == FALSE){
          AtEQ=TRUE; 
          EqGen=gen; 
          print(c("Eq has been reached at gen =", gen))
          } #if at equilibrium, tell it
			}#end if AtEQ 
			
			if (AtEQ==TRUE & gen==(EqGen+1)){
			Fecundity_Vect <- PopVect + r_growth*PopVect*(KVect-PopVect)/KVect
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
		sourceN <- mutVECT*Fecundity_Vect
		migVECT<-migration_C(sourceN, X_Demes, Y_Demes, MigMat)*OccupiedVect
		migVECT[OccupiedVect<=0] = 0	# #this step removes alleles from demes with no individuals.  multiplying doesn't work.
		migVECT[OccupiedVect==1] = migVECT[OccupiedVect==1]/PopVect[OccupiedVect==1]
			#this step scales allele freqs by deme size
				# the scaling must be done with the "fractional" population matrix, if it is rounded off you can get allele freqs > 1
		
		
		############################################################
		#Step 3: Selection on offspring
		############################################################
		if(LOCUS_TYPE=="SELECTION"){
				selVECT <- selection_C( migVECT, s_VECT)*OccupiedVect
			}
		
		if(LOCUS_TYPE=="NEUTRAL"){
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

}

