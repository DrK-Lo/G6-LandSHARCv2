MakeRefugia <- function(N_start, Side, Loc, X_Locs, Y_Locs, K_cell, CellSize){
	
	#Refugia are spread out equally over a 20km area, so would no work in corners of cell	
	#	N_start<-10000 #starting pop size for first refugia
	#	Side <- "bottom", "left", "top" , "right"
	# 	Loc = 181 on that Side is in km
	#   K_cell is carrying capacity for that cell
	
	
	#Insert general starting location in grid units (in km) for each refugia
	
	Loc_start <- which(abs(X_Locs-Loc)==min(abs(X_Locs-Loc)))
	#gives starting cell
	#good for square matrix
	
		X_cells <- 1:length(X_Locs)
		Y_cells <- 1:length(Y_Locs)
	
	
	#Figure out how many cells will be filled in each refugia
		numcells <- round(N_start/K_cell)
		remainder <- N_start%%K_cell
		N_eachCell <- rep(K_cell, numcells)	
		if (remainder>=1){numcells_tot=numcells+1;N_eachCell <-c(N_eachCell, remainder)}else{numcells_tot=numcells}	
	
	
	# If the starting pop size fits in one cell	
		if (numcells==0){x_start_locs_i <- start_x
			y_start_locs_i <- start_y
			break;
		}
	
	
	# Stack cells on top of each other over a 20km distance
		numCellsPerRow <- round(20/CellSize)
		numRows <- floor(numcells_tot/numCellsPerRow)
		xtraCells <- numcells_tot - numRows*numCellsPerRow
		
	#Cell Locations for "rows" (taken from perspective of side)
	Central_Locs <- rep( (Loc_start - round(numCellsPerRow/2)) : (Loc_start - 1- round(numCellsPerRow/2) + numCellsPerRow) ,numRows )
		
	Row_Locs <- rep(1:numRows, each=numCellsPerRow)
	
	
	if(xtraCells>0){
	Xtra_Locs <- (Loc_start - round(xtraCells/2)) : (Loc_start - round(xtraCells/2) + xtraCells -1)
	Xtra_RowLoc <- 	rep(numRows+1, xtraCells)
		Central_Locs <- c(Central_Locs, Xtra_Locs)
		Row_Locs <- c(Row_Locs, Xtra_RowLoc)
		}
		
	
	if (Side=="bottom"){start_x <- Central_Locs; start_y=Row_Locs}
	if (Side=="left"){start_x <- Row_Locs; start_y=Central_Locs}
	if (Side=="top"){start_x <- Central_Locs; start_y= max(Y_cells)+1-Row_Locs}
	if (Side=="right"){start_x <- max(X_cells)+1-Row_Locs; start_y=Central_Locs}
		
		
	AllLocs <- cbind(start_x,start_y,N_eachCell)
	
	return(AllLocs)
	
}