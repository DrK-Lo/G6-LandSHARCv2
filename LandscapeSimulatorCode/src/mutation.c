#include <R.h> 
#include <Rdefines.h>
#include <Rmath.h> 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void mutation(double *u, double *P_freqMAT, int *length)
{
	int i;
	for (i=0; i<*length; i++) {
		P_freqMAT[i] = (1-*u) * P_freqMAT[i] + *u * (1-P_freqMAT[i]);
	}	
}