/* Creates an array containing a normalised initial wavefunction from the Thomas-Fermi approximation.*/

//Included libraries.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
//#include <cufftw.h>

//Designated header file with input parameters.
#include "BEC_3D.h"

//Takes the input of an empty wavefunction, two empty potentials and the x, y, z grid vectors.
void initial_cond(fftw_complex *w, double *v1, double *v2, double x[], double y[], double z[])
{
	int i, j, k;
	double utf, V, X, Y, Z, Xt, Yt, Zt;
	fftw_complex zero;

	zero = 0.0 + _Complex_I*0.0;


	//Filling the x, y, z grid vectors.
	x[0] = xSTART;
	y[0] = ySTART;
	z[0] = zSTART;
	
	for (i = 0; i < Xpoints; i++) {
		x[i] = x[0] + i * DELTAx;
	}
	
	for (i = 0; i < Ypoints; i++) {
		y[i] = y[0] + i * DELTAy;
	}
	
	for (i = 0; i < Zpoints; i++) {
		z[i] = z[0] + i * DELTAz;
	}
	
	//Calculating the Thomas-Fermi chemical energy.
	utf = pow(15 * C * L / (64 * PI), 0.4);
	
	for(i = 0; i < Zpoints; i++) {
		for(j = 0; j < Ypoints; j++) {
			for(k = 0; k < Xpoints; k++) {
				
				X = x[k];
				Y = y[j];
				Z = z[i];
				
				//Filling the potential matrix for the level potential.
				*(v2 + i * Ypoints * Xpoints + j * Xpoints + k) = 0.25 * (X * X + Y * Y + L * L * Z * Z);

				Xt = X;
				Yt = Y * cos(theta) + Z * sin(theta);
				Zt = -Y * sin(theta) + Z * cos(theta);
				
				V = 0.25 * (Xt * Xt + Yt * Yt + L * L * Zt * Zt);
				
				//Filling the potential matrix for the tilted potential.
				*(v1 + i * Ypoints * Xpoints + j * Xpoints + k) = V;
				
				//Filling the initial condition with the Thomas-Fermi approximation.
				if ((utf - V) < 0)
				  //*(w + i * Ypoints * Xpoints + j * Xpoints + k) = (fftw_complex) 0.0;
				  *(w + i * Ypoints * Xpoints + j * Xpoints + k) = 0.0;
				else					
				  *(w + i * Ypoints * Xpoints + j * Xpoints + k) = (fftw_complex) pow((utf - V) / C,0.5);
				  //*(w + i * Ypoints * Xpoints + j * Xpoints + k) = make_cuDoubleComplex( pow((utf - V) / C,0.5), 0.0 );
				  
			}
		}
	}
}
