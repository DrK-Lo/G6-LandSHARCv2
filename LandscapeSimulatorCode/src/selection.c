#include <R.h> 
#include <Rdefines.h>
#include <Rmath.h> 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void selection( double *P_freq,double *sel_VECT, int *length)
{
	int i;
	for (i=0; i<*length; i++) {
		P_freq[i] = P_freq[i] * (1 +sel_VECT[i]) / (P_freq[i] * (1 + sel_VECT[i]) + (1 - P_freq[i]));
	}	
}
