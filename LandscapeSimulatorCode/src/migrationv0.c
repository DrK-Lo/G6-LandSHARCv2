#include <R.h> 
#include <Rdefines.h>
#include <Rmath.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//system("R CMD SHLIB src/migration.c")
//dyn.load("src/migration.so")

void migration(double *P_freqMAT, int *xdim, int *ydim, double *pnew, double *MigVect, int *MigVectDim, int *xtra)
{
	int k;			//index to loop through each deme
	int row=0; int i; //index to loop through rows
	int col=0; int j; //index to loop through columns

	int length=*xdim * *ydim;

	//double MigVect[]={0,0,0.0003234533,0.0007520515,0.0009963064,0.0007520515,0.0003234533,0,0,0,0.0005676783,0.002316482,0.005385981,0.007135266,0.005385981,0.002316482,0.0005676783,0,0.0003234533,0.002316482,0.009452692,0.02197817,0.02911634,0.02197817,0.009452692,0.002316482,0.0003234533,0.0007520515,0.005385981,0.02197817,0.05110077,0.06769752,0.05110077,0.02197817,0.005385981,0.0007520515,0.0009963064,0.007135266,0.02911634,0.06769752,0.08968464,0.06769752,0.02911634,0.007135266,0.0009963064,0.0007520515,0.005385981,0.02197817,0.05110077,0.06769752,0.05110077,0.02197817,0.005385981,0.0007520515,0.0003234533,0.002316482,0.009452692,0.02197817,0.02911634,0.02197817,0.009452692,0.002316482,0.0003234533,0,0.0005676783,0.002316482,0.005385981,0.007135266,0.005385981,0.002316482,0.0005676783,0,0,0,0.0003234533,0.0007520515,0.0009963064,0.0007520515,0.0003234533,0,0};
	//int MigVectDim=9;
	//int xtra=4;
	int xdim_ext = *xdim + *xtra*2;
	int ydim_ext = *ydim + *xtra*2;
	int index=0;
	int migIndex=0;

	double focalfreq;
	int focalCol; int focalRow; int focalInd;
	int localRow; int localCol;


	//loop through each deme in the P_freqMAT (which is input as a vector)
	for (k=0; k<length; k++){
		focalfreq = P_freqMAT[k];
		//loop through columns in MigVect (if the dimension of MigVect changes, so does the -4 here)
		//start at upper left of migration matrix, move through columns until end of row, etc.
		//in this case, the center of the migration matrix (MigVect[40]) is the focal cell

		localRow = floor(k / *xdim);
		localCol = k - localRow * *xdim;
		//focalInd = *ydim * localRow + localCol;


		focalRow= localRow + *xtra;
		focalCol= localCol + *xtra ;
		focalInd= xdim_ext * focalRow + focalCol;
	//	Rprintf("%d %d %d %d %d %d %d \n", localRow, localCol, focalInd, k, focalRow, focalCol, focalInd);
		migIndex=0;

	//	pnew[focalInd]= focalInd;

		//loop through all demes that will exchange migrants with focal deme in X direction
		for (i = -*xtra; i< (*xtra+1); i++){
			//loop through all demes that will exchange migrants with focal deme in Y direction
			for (j = -*xtra; j< (*xtra+1); j++) {
				col = focalCol + j;
				row = focalRow + i;

				index = row * xdim_ext + col;
		//	Rprintf("%d %d %d %d %d %d \n",focalRow, focalCol, row, col, index, migIndex);

				pnew[index] += focalfreq * MigVect[migIndex];

			migIndex +=1;
			}	//end loop through cols
		} //end loop through rows



	}// end loop through all

	//pnew is subset when returned to R
}
