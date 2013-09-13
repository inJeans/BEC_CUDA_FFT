/*Functions for basic matrix operations, integrating to 2D and normalisation.*/

//Included libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
//#include <cufftw.h>

//Designated header file with input parameters.
#include "BEC_3D.h"



//Computes the weighted sum of two matrices.
void weight_sum(double a, fftw_complex *q, double b, fftw_complex *r)
{
	int i;
	
	for(i = 0; i < Zpoints * Ypoints * Xpoints; i++) {
		*(q + i) = *(q + i) * a + *(r + i) * b;
	}
}


//Creates a memory copy of the given matrix.
fftw_complex *copy(fftw_complex *q)
{
	int i;
	fftw_complex *p;
	
	p = fftw_malloc(Xpoints * Ypoints * Zpoints * sizeof(fftw_complex));
	
	for(i = 0; i < Zpoints * Ypoints * Xpoints; i++) {
		*(p + i) = *(q + i);
	}
	return p;
}


//Computes the modulus squared of a 3D complex wavefunction.
double *modulus3(fftw_complex *q)
{
	int i;
	double *p;
	
	p = malloc(Xpoints * Ypoints * Zpoints * sizeof(double));
	
	for(i = 0; i < Xpoints * Ypoints * Zpoints; i++) {
		*(p + i) = pow(cabs(*(q + i)), 2);
	}
	return (double *) p;
}


//Integrates a 3D wavefunction through the Z direction onto the XY plane.
double *integrateZ(double *q)
{
	int i, j, k;
	double *p, sum;
	
	sum = 0;
	
	p = malloc(Xpoints * Ypoints * sizeof(double));
	
	for(i = 0; i < Ypoints; i++) {	
		for(j = 0; j < Xpoints; j++) {
			for(k = 0; k < Zpoints; k++) {
				sum += *(q + k * Ypoints * Xpoints + i * Xpoints + j);
			}
			*(p + i * Xpoints + j) = sum * DELTAz;
			sum = 0;	
		}
	}
	
	return p;
}


//integrates a 3D wavefunction through the X direction onto the YZ plane.
double *integrateX(double *q)
{
	int i, j, k;
	double *p, sum;
	
	sum = 0;
	
	p = malloc(Ypoints * Zpoints * sizeof(double));
	
	for(i = 0; i < Zpoints; i++) {	
		for(j = 0; j < Ypoints; j++) {
			for(k = 0; k < Xpoints; k++) {
				sum += *(q + i * Ypoints * Xpoints + j * Xpoints + k);
			}
			*(p + i * Ypoints + j) = sum * DELTAx;
			sum = 0;
		}
	}
	
	return p;
}


//Computes the phase at a set Z level (k increment) across the XY plane.
double *phase(fftw_complex *q)
{
	int i, j, k;
	double *p;
	
	k = Zpoints/2;
	
	p = malloc(Xpoints * Ypoints * sizeof(double));
	
	for(i = 0; i < Ypoints; i++) {	
		for(j = 0; j < Xpoints; j++) {
			*(p + i * Xpoints + j) = carg(*(q + k * Ypoints * Xpoints + i * Xpoints + j));
		}
	}
	
	return p;
}


//Normalises a complex 3D wavefunction.
void norm(fftw_complex *q)
{
	int i;
	fftw_complex sum;
	sum = 0;

	for(i = 0; i < Xpoints * Ypoints * Zpoints; i++) {	
		sum += pow(cabs(*(q + i)), 2);
	}
	
	sum = sum * DELTAx * DELTAy * DELTAz;
	
	for(i = 0; i < Xpoints * Ypoints * Zpoints; i++) {	
		*(q + i) = *(q + i) / sqrt(sum);
	}
}

