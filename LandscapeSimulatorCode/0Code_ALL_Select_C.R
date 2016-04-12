
ALLSelect <- function(CodeDir, headDir, subfolder, FileName_ResultsFST, FSTgen=100, FileName_LandscapeMetadata, MODEL_TYPE, LOCUS_TYPE, X_Demes, Y_Demes, Scale=1, u=1e-06, p_start = runif(1), GEN, LandscapeMatGen=GEN, MOVIE=FALSE , MOVIEgen = GEN/2, DispersalKernel_FilePath = NA, K_All = NA	, GEN_Refugia = NA, r_growth = NA	, Refugia_StartMat_FilePath = NA, AlleeDensity = 2, MOVIE_Density = FALSE	, Refugia_Selection = NA , m = NA	, IslandSize = NA, Island_N = NA, ENVI_ID =  NA, ENVI_DIR = NA,  s_high = NA, s_low=-s_high){
	
	startT <- tic()
	
	
#############################################################		
	########### ARE INPUT FILES THERE ?????? #########################
	##################################################################
	if (is.null(CodeDir)==TRUE){print ("Missing the directory for code: CodeDir")}
	if (is.null(headDir)==TRUE){print ("Missing the head directory for output: headDir")}
	if (is.null(subfolder)==TRUE){print ("Missing the subfolder within the head directory for output: subfolder")}
	if (is.null(FileName_ResultsFST)==TRUE){print ("Missing the FST results filename")}
	if (is.null(FileName_LandscapeMetadata)==TRUE){print ("Missing the LandscapeMetadata filename")}
	
	
	if (is.null(MODEL_TYPE)==TRUE){print("Please enter  model type (IBD, REFUGIA, ISLAND)")	}
	if (is.null(LOCUS_TYPE)==TRUE){print("Please enter locus type (NEUTRAL, SELECTION)")	}
	if (is.null(X_Demes)==TRUE){print("Missing Landscape Size (X)")}	
	if (is.null(Y_Demes)==TRUE){print("Missing Landscape Size (Y)")}
		
	if (is.null(u)==TRUE){print ("Missing mutation rate: u")}
	if (is.null(p_start)==TRUE){print ("Missing starting allele freq: p_start")}
	
	if (is.null(GEN)==TRUE){print ("Missing GEN")}
	
	
	if (LOCUS_TYPE == "SELECTION") {
			if (is.null(ENVI_ID)==TRUE){print ("Missing environment file: ENVI_ID")}
			if (is.null(ENVI_DIR)==TRUE){print ("Missing directory for the ENVI file: ENVI_DIR")}
			if (is.null(s_low)==TRUE){print ("Missing s_low")}
			if (is.null(s_high)==TRUE){print ("Missing s_high")}
		}



#############################################################	
	#########################################
	#Get name for this run
	#########################################
#############################################################		
	runname <- getRunName() 
	SetOrCreateDir(headDir)
	 if (file.exists(paste(headDir, subfolder, sep=""))){}else{dir.create(paste(headDir, subfolder, sep=""))}
	   MovieDirname <- paste(headDir,subfolder,"/", runname,"Movie",sep="")
        dir.create(MovieDirname)
	
#############################################################			
	#########################################
	#Determine scaling for population based on Scale
	#########################################	
#############################################################		
		X_Dist_Adj <- X_Demes*Scale
		Y_Dist_Adj <- Y_Demes*Scale
		X_Locs <- seq(0,X_Dist_Adj, by=Scale)
		Y_Locs <- seq(0,Y_Dist_Adj, by=Scale)
		
#############################################################		
	#########################################
	####If LOCUS_TYPE=="SELECTION", Input Environment
	#########################################	
#############################################################	
if (LOCUS_TYPE=="SELECTION"){
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

 	}
#############################################################			
	#########################################
	#Print Run Metadata
	#########################################	
#############################################################	
		setwd(MovieDirname)
		Metadatafile <- paste(runname,"Metadata",sep="")
		write(runname, file= Metadatafile) 
		write(paste("Date=",date()), file= Metadatafile, append=TRUE)
		write(paste("CodeDir=",CodeDir), file= Metadatafile, append=TRUE)
		write(paste("headDir=",headDir), file= Metadatafile, append=TRUE)
		write(paste("subfolder=",subfolder), file= Metadatafile, append=TRUE) 
		write(paste("FileName_ResultsFST=",FileName_ResultsFST), file= Metadatafile, append=TRUE)			
		write(paste("FSTgen=",FSTgen), file= Metadatafile, append=TRUE)
		write(paste("FileName_LandscapeMetadata=",FileName_LandscapeMetadata), file= Metadatafile, append=TRUE) 
		write(paste("MODEL_TYPE=",MODEL_TYPE), file= Metadatafile, append=TRUE)
		write(paste("LOCUS_TYPE=",LOCUS_TYPE), file= Metadatafile, append=TRUE)
		write(paste("X_Demes=",X_Demes), file= Metadatafile, append=TRUE)
		write(paste("Y_Demes=",Y_Demes), file= Metadatafile, append=TRUE)
		write(paste("Scale=",Scale), file= Metadatafile, append=TRUE)
		write(paste("u=",u), file= Metadatafile, append=TRUE)			
		write(paste("p_start=",p_start), file= Metadatafile, append=TRUE)
		write(paste("GEN=",GEN), file= Metadatafile, append=TRUE)
		write(paste("LandscapeMatGen=",LandscapeMatGen), file= Metadatafile, append=TRUE)
		write(paste("MOVIE=",MOVIE), file= Metadatafile, append=TRUE)
		write(paste("MOVIEgen=",MOVIEgen), file= Metadatafile, append=TRUE)
		write(paste("DispersalKernel_FilePath=",DispersalKernel_FilePath), file= Metadatafile, append=TRUE)
		write(paste("K_All=",K_All), file= Metadatafile, append=TRUE)
		write(paste("GEN_Refugia=",GEN_Refugia), file= Metadatafile, append=TRUE)
		write(paste("r_growth=",r_growth), file= Metadatafile, append=TRUE)
		write(paste("Refugia_StartMat_FilePath=",Refugia_StartMat_FilePath), file= Metadatafile, append=TRUE)
		write(paste("AlleeDensity=",AlleeDensity), file= Metadatafile, append=TRUE)
		write(paste("MOVIE_Density=",MOVIE_Density), file= Metadatafile, append=TRUE)
		write(paste("Refugia_Selection=",Refugia_Selection), file= Metadatafile, append=TRUE)
		write(paste("m=",m), file= Metadatafile, append=TRUE)
		write(paste("IslandSize=",IslandSize), file= Metadatafile, append=TRUE)
		write(paste("Island_N=",Island_N), file= Metadatafile, append=TRUE)
		write(paste("ENVI_ID=",ENVI_ID), file= Metadatafile, append=TRUE)
		write(paste("ENVI_DIR=",ENVI_DIR), file= Metadatafile, append=TRUE)
		write(paste("s_low=",s_low), file= Metadatafile, append=TRUE)
		write(paste("s_high=",s_high), file= Metadatafile, append=TRUE)
		
#############################################################			
	#########################################
	#Setup starting population and allele frequencies
	#########################################			
#############################################################	
	
	### ISLAND MODEL
	##############################
	if (MODEL_TYPE=="ISLAND"){
		TotIslandSize <- IslandSize*IslandSize
		NumIslands <- (X_Demes*Y_Demes)/(IslandSize*IslandSize)	
		N_vect <- rep(Island_N, NumIslands)
		vectorsize <- NumIslands	
		OccupiedVect <- rep(1, vectorsize)
		P_freqMAT <- rep(p_start, vectorsize)
		DensityPerKM <- round(sum(N_vect)/(X_Dist_Adj*Y_Dist_Adj), 3)
		
		#Need to convert selection-vector to same size as island vector
		if (LOCUS_TYPE=="SELECTION"){
		s_VECT <- ConvertEnviToIslands(s_VECT, IslandSize, X_Demes, Y_Demes)
		Envi_VECT <- ConvertEnviToIslands(Envi_VECT, IslandSize, X_Demes, Y_Demes)
		}
	}
		
	### IBD MODEL OR REFUGIA MODEL
	####################################		
	if (MODEL_TYPE=="IBD" | MODEL_TYPE=="Refugia"){
		vectorsize<- X_Demes*Y_Demes
		MigMat <- as.matrix(read.table(DispersalKernel_FilePath))
		MigVect <- as.vector(MigMat)
		DensityPerKM <- round((K_All*X_Demes*Y_Demes)/(X_Dist_Adj*Y_Dist_Adj), 3)
		K_Vect <- rep(K_All, vectorsize)
	}

	### IBD MODEL SPECIFICS
	####################################			
	if (MODEL_TYPE=="IBD"){
		InitialPopVect <- K_Vect
	}

	### REFUGIA MODEL SPECIFICS
	####################################			
	if (MODEL_TYPE=="Refugia"){
		#rename variables (naming system set by IBD and Island)
			Post_Refugia_Gen <- GEN
			GEN <- GEN_Refugia
		
		#Setup Intial population
			InitialPopVect <- as.vector(as.matrix(scan(Refugia_StartMat_FilePath, quiet=TRUE)))	

	}
	
		### IBD MODEL OR REFUGIA MODEL
	####################################		
	if (MODEL_TYPE=="IBD" | MODEL_TYPE=="Refugia"){
			OccupiedVect <- InitialPopVect
			OccupiedVect[which(InitialPopVect>0)]=1
			N_vect <- InitialPopVect
			P_scaleVECT <- migration_C(OccupiedVect, X_Demes, Y_Demes, MigMat)

			P_freqMAT <- rep(p_start, vectorsize)*OccupiedVect
		
	}
		
#############################################################		
	#########################################
	# Loop over generations and let evolution happen
	# For the refugia model, this is the first set of generations
	# before population expansion is allowed to occur
	##########################################
#############################################################################
	for (gen in 1:GEN){
		
		############################################################
		#Step 1: Mutation
		############################################################
		mutVECT <- mutation_C(u,P_freqMAT)*OccupiedVect
		
		############################################################
		#Step 2: Dispersal/migration
		############################################################
		if(MODEL_TYPE=="ISLAND"){migVECT<- (1-m)*mutVECT + m*mean(mutVECT)}
		if(MODEL_TYPE=="IBD" | MODEL_TYPE=="Refugia" ){migVECT<-migration_C(mutVECT, X_Demes, Y_Demes, MigMat)*OccupiedVect
				migVECT[OccupiedVect==1] = migVECT[OccupiedVect==1]/P_scaleVECT[OccupiedVect==1]
			}
		
		############################################################
		#Step 3: Selection
		############################################################
		if(LOCUS_TYPE=="SELECTION"){
			if (length(s_VECT)!=length(migVECT)){
				print("Vector lengths not equal"); break}
			
			if(MODEL_TYPE=="ISLAND" | MODEL_TYPE=="IBD"){
				selVECT <- selection_C(migVECT, s_VECT)*OccupiedVect
			}
			#note that just applying OneLocusSelection to two vecotrs will not give correct answer
			
			if(MODEL_TYPE=="Refugia" & Refugia_Selection==TRUE){
				selVECT <- selection_C(migVECT, s_VECT)*OccupiedVect
			}
			if(MODEL_TYPE=="Refugia" & Refugia_Selection==FALSE){
				selVECT <- migVECT
			}
		}	
		
		
		if(LOCUS_TYPE=="NEUTRAL"){
			selVECT <- migVECT
		}
			
		############################################################
		#Step 4: Drift
		############################################################
		if (length(which(is.na(selVECT)==TRUE))>1){print(c("NAs in selVECT gen", gen)); return(list(mut=mutVECT, mig=migVECT, sel=selVECT, drift=driftVECT)); break}
		
		 driftVECT <- drift_C(N_vect, selVECT)*OccupiedVect
		
				if (length(which(is.na(driftVECT)==TRUE))>1){print("NAs in driftVECT"); return(list(mut=mutVECT, mig=migVECT, sel=selVECT, drift=driftVECT)); break}
		############################################################
		#Write desired outputs
		############################################################
		if (MODEL_TYPE=="Refugia" & gen%%FSTgen==0){print(c("Refugia generation:",gen))}
			
		if (MODEL_TYPE=="ISLAND" | MODEL_TYPE == "IBD"){	
		
			#Write FST of landscape to a file
			##################################			
			if (gen%%FSTgen==0){
				print(gen)
				#print(c(max(selVECT), min(selVECT)))
				FST_All<- WC_FST_LargeSample(migVECT)
				setwd(headDir)	
				
				
				if (LOCUS_TYPE=="SELECTION"){
					test <- cor.test(Envi_VECT, migVECT)
					Landscape_Corr <- test$est
				}else{Landscape_Corr <-NA}
				
				
				fst_out = c(runname, gen, FST_All, Landscape_Corr)
				write(fst_out, FileName_ResultsFST , append=TRUE, ncolumns=length(fst_out))
				setwd(MovieDirname)
				write(fst_out, Metadatafile , append=TRUE, ncolumns=length(fst_out))
				}
				
			#Write pdf Image of allele frequencies to a file
			##################################	
			if (gen%%MOVIEgen==0 & MOVIE==TRUE){
				 setwd(MovieDirname)
				if(MODEL_TYPE=="ISLAND"){P_plotMAT <- ConvertVectToLS(X_Demes, Y_Demes, migVECT, IslandSize)
					}else{
						P_plotMAT = migVECT
		          		P_plotMAT[OccupiedVect=0]=NA
		          }
				pdf(file=paste(runname,"ImageAfterMig_",gen, ".pdf",sep=""), width=4, height=3, bg="white")
	             par(mar=c(3,3,1,1))
	             plotAlleleFreqImage(X_Locs, Y_Locs, P_plotMAT, gen, FST_All)	#need to edit this function
	             dev.off()			
			}	
			
			
			#Write matrices of allele frequencies to a file
			##################################		
			if (gen%%LandscapeMatGen==0){
				setwd(MovieDirname) 
				if(MODEL_TYPE=="ISLAND"){P_mat <- ConvertVectToLS(X_Demes, Y_Demes, migVECT, IslandSize)
					}else{
						P_mat <- migVECT
						}
	            write(as.vector(P_mat), file=paste(runname,"Matrix_AlleleFreqs_Gen", gen,sep=""), ncolumns=Y_Demes)
			}
				
		}#end output if Island or IBD
	
		
		############################################################
		#Step 4: Reassign allele freqs to next generation
		############################################################	
		P_freqMAT <- driftVECT	
			
}#end loop through generations		
	#if this was the refugia model, this was the time in the refugia
############################################################				
################### Post-refugia loops #########################################
############################################################		


if (MODEL_TYPE=="Refugia"){
		PopVect2<-InitialPopVect
		AtEQ <- FALSE
		EqGen<-NA
	
		for (gen in 1:Post_Refugia_Gen){
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

			#Write FST of landscape to a file
			##################################			
			if (gen%%FSTgen==0){
				print(gen)
				FST_All<- WC_FST_LargeSample(migVECT[OccupiedVect==1])
				setwd(headDir)	
				
				
				if (LOCUS_TYPE=="SELECTION"){
					test <- cor.test(Envi_VECT[OccupiedVect==1], migVECT[OccupiedVect==1])
					Landscape_Corr <- test$est
				}else{Landscape_Corr <-NA}
				
				
				fst_out = c(runname, gen, FST_All, Landscape_Corr)
				write(fst_out, FileName_ResultsFST , append=TRUE, ncolumns=length(fst_out))
				
				}
				
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


############################################################				
######################## Write Ouputs ####################################
############################################################			
	#################################
	#Write time elapsed to metadata
	#################################	
		setwd(MovieDirname)
			write(" ", file=Metadatafile, append=TRUE)
			write(paste("Number generations elapsed =", gen), file=Metadatafile, append=TRUE)

			p_LS_end <- mean(migVECT)
			He_LS_start	<- 	p_start*(1-p_start)	
			He_LS_end <- mean(migVECT*(1-migVECT))

			write(paste("p_LS_end=",p_LS_end), file= Metadatafile, append=TRUE)
			write(paste("He_LS_start=",p_LS_end), file= Metadatafile, append=TRUE)
			write(paste("He_LS_end=",p_LS_end), file= Metadatafile, append=TRUE)
			endT <- toc()
			TotalT <- round(endT-startT)
			write(paste("Time Elapsed (sec)",TotalT), file=Metadatafile, append=TRUE)

	print(c(runname, p_start))

	#################################
	#Write main output to a file in main directory
	#################################	
		setwd(headDir)	
		allVars <- c(MODEL_TYPE, LOCUS_TYPE, runname, X_Demes, Y_Demes, Scale, u, GEN, LandscapeMatGen, FSTgen, MOVIE, MOVIEgen, m, IslandSize, Island_N, DispersalKernel_FilePath, K_All, GEN_Refugia, r_growth, Refugia_StartMat_FilePath, AlleeDensity, MOVIE_Density, Refugia_Selection, ENVI_ID, ENVI_DIR, s_low, s_high,  p_start, p_LS_end, He_LS_start, He_LS_end, FST_All, Landscape_Corr, TotalT, "\n", sep=" ")
		write(allVars, FileName_LandscapeMetadata , append=TRUE, ncolumns=length(allVars))

}#end function
