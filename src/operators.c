/* File contains key operators including the linear (diffusion) term and non-linear (inc. potential)*/

//Included libraries.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
//#include <cufftw.h>

//Designated header file with input parameters.
#include "BEC_3D.h"


/*Calculates the diffusion coefficients for progressing in Imaginary Time.
Computed once at the start of the algorithm.
Input is a pointer to a single empty matrix.*/
void DCi(fftw_complex *p)
{
	int l, m, n;
	double k_lmn2;

	//Filling a grid with diffusion coefficients. Each section of the rectangular prism is considered seperately.	
	for(n = 0; n < Zpoints; n++) {
			
		for(m = 0; m < Ypoints; m++) {
				
			for(l = 0; l < Xpoints; l++) {
									
				if(n <= Zpoints / 2 && m <= Ypoints / 2 && l <= Xpoints / 2)
					k_lmn2 = pow(2 * PI * l / (DELTAx * Xpoints), 2) + pow(2 * PI * m / (DELTAy * Ypoints), 2) + pow(2 * PI * n / (DELTAz * Zpoints), 2);
					
				else
				if(n <= Zpoints / 2 && m <= Ypoints / 2 && l > Xpoints / 2)
					k_lmn2 = pow( - 2 * PI * (Xpoints - l) / (DELTAx * Xpoints), 2) + pow(2 * PI * m / (DELTAy * Ypoints), 2) + pow(2 * PI * n / (DELTAz * Zpoints), 2);
					
				else
				if(n <= Zpoints / 2 && m > Ypoints / 2 && l <= Xpoints / 2)
					k_lmn2 = pow(2 * PI * l / (DELTAx * Xpoints), 2) + pow( - 2 * PI * (Ypoints - m) / (DELTAy * Ypoints), 2) + pow(2 * PI * n / (DELTAz * Zpoints), 2);
					
				else
				if(n <= Zpoints / 2 && m > Ypoints / 2 && l > Xpoints / 2)
					k_lmn2 = pow( - 2 * PI * (Xpoints - l) / (DELTAx * Xpoints), 2) + pow( - 2 * PI * (Ypoints - m) / (DELTAy * Ypoints), 2) + pow(2 * PI * n / (DELTAz * Zpoints), 2);
						
						
						
						
				else
				if(n > Zpoints / 2 && m <= Ypoints / 2 && l <= Xpoints / 2)
					k_lmn2 = pow(2 * PI * l / (DELTAx * Xpoints), 2) + pow(2 * PI * m / (DELTAy * Ypoints), 2) + pow( - 2 * PI * (Zpoints - n) / (DELTAz * Zpoints), 2);
					
				else
				if(n > Zpoints / 2 && m <= Ypoints / 2 && l > Xpoints / 2)
					k_lmn2 = pow( - 2 * PI * (Xpoints - l) / (DELTAx * Xpoints), 2) + pow(2 * PI * m / (DELTAy * Ypoints), 2) + pow( - 2 * PI * (Zpoints - n) / (DELTAz * Zpoints), 2);
					
				else
				if(n > Zpoints / 2 && m > Ypoints / 2 && l <= Xpoints / 2)
					k_lmn2 = pow(2 * PI * l / (DELTAx * Xpoints), 2) + pow( - 2 * PI * (Ypoints - m) / (DELTAy * Ypoints), 2) + pow( - 2 * PI * (Zpoints - n) / (DELTAz * Zpoints), 2);
					
				else
				if(n > Zpoints / 2 && m > Ypoints / 2 && l > Xpoints / 2)
					k_lmn2 = pow( - 2 * PI * (Xpoints - l) / (DELTAx * Xpoints), 2) + pow( - 2 * PI * (Ypoints - m) / (DELTAy * Ypoints), 2) + pow( - 2 * PI * (Zpoints - n) / (DELTAz * Zpoints), 2);	
				
			*(p + n * Ypoints * Xpoints + m * Xpoints + l) = cexp(- (DELTAt / 2) * k_lmn2) / (Xpoints * Ypoints * Zpoints);
			}
		}			
	}		
}


/*Calculates the diffusion coefficients for progressing in real time.
Computed once after progression in Imaginary Time is complete.
Input is a pointer to a single empty matrix.*/
void DCr(fftw_complex *p)
{
	int l, m, n;
	double k_lmn2;
		
	//Filling the grid with diffusion coefficients. Each section of the rectangular prism is considered seperately.	
	for(n = 0; n < Zpoints; n++) {
			
		for(m = 0; m < Ypoints; m++) {
				
			for(l = 0; l < Xpoints; l++) {
									
				if(n <= Zpoints / 2 && m <= Ypoints / 2 && l <= Xpoints / 2)
					k_lmn2 = pow(2 * PI * l / (DELTAx * Xpoints), 2) + pow(2 * PI * m / (DELTAy * Ypoints), 2) + pow(2 * PI * n / (DELTAz * Zpoints), 2);
					
				else
				if(n <= Zpoints / 2 && m <= Ypoints / 2 && l > Xpoints / 2)
					k_lmn2 = pow( - 2 * PI * (Xpoints - l) / (DELTAx * Xpoints), 2) + pow(2 * PI * m / (DELTAy * Ypoints), 2) + pow(2 * PI * n / (DELTAz * Zpoints), 2);
					
				else
				if(n <= Zpoints / 2 && m > Ypoints / 2 && l <= Xpoints / 2)
					k_lmn2 = pow(2 * PI * l / (DELTAx * Xpoints), 2) + pow( - 2 * PI * (Ypoints - m) / (DELTAy * Ypoints), 2) + pow(2 * PI * n / (DELTAz * Zpoints), 2);
					
				else
				if(n <= Zpoints / 2 && m > Ypoints / 2 && l > Xpoints / 2)
					k_lmn2 = pow( - 2 * PI * (Xpoints - l) / (DELTAx * Xpoints), 2) + pow( - 2 * PI * (Ypoints - m) / (DELTAy * Ypoints), 2) + pow(2 * PI * n / (DELTAz * Zpoints), 2);
						
						
						
						
				else
				if(n > Zpoints / 2 && m <= Ypoints / 2 && l <= Xpoints / 2)
					k_lmn2 = pow(2 * PI * l / (DELTAx * Xpoints), 2) + pow(2 * PI * m / (DELTAy * Ypoints), 2) + pow( - 2 * PI * (Zpoints - n) / (DELTAz * Zpoints), 2);
					
				else
				if(n > Zpoints / 2 && m <= Ypoints / 2 && l > Xpoints / 2)
					k_lmn2 = pow( - 2 * PI * (Xpoints - l) / (DELTAx * Xpoints), 2) + pow(2 * PI * m / (DELTAy * Ypoints), 2) + pow( - 2 * PI * (Zpoints - n) / (DELTAz * Zpoints), 2);
					
				else
				if(n > Zpoints / 2 && m > Ypoints / 2 && l <= Xpoints / 2)
					k_lmn2 = pow(2 * PI * l / (DELTAx * Xpoints), 2) + pow( - 2 * PI * (Ypoints - m) / (DELTAy * Ypoints), 2) + pow( - 2 * PI * (Zpoints - n) / (DELTAz * Zpoints), 2);
					
				else
				if(n > Zpoints / 2 && m > Ypoints / 2 && l > Xpoints / 2)
					k_lmn2 = pow( - 2 * PI * (Xpoints - l) / (DELTAx * Xpoints), 2) + pow( - 2 * PI * (Ypoints - m) / (DELTAy * Ypoints), 2) + pow( - 2 * PI * (Zpoints - n) / (DELTAz * Zpoints), 2);	
				
			*(p + n * Ypoints * Xpoints + m * Xpoints + l) = cexp(- I * (DELTAt / 2) * k_lmn2) / (Xpoints * Ypoints * Zpoints);
			}
		}			
	}		
}

/*Non linear operator, containing the potential for Imaginary Time Progression.
Input is pointers to the wavefunction, potential and the x,y,z grid vectors.*/
void Ni(fftw_complex *w, double *v, double x[], double y[], double z[])
{
	int i, j, k;
	
	//Iterating through the grid.
	for(i = 0; i < Zpoints; i++) {
		for(j = 0; j < Ypoints; j++) {
			for(k = 0; k < Xpoints; k++) {
				
				*(w + i * Ypoints * Xpoints + j * Xpoints + k) = -DELTAt * (*(v + i * Ypoints * Xpoints + j * Xpoints + k) + C * pow(cabs(*(w + i * Ypoints * Xpoints + j * Xpoints + k)), 2)) * 
				(*(w + i * Ypoints * Xpoints + j * Xpoints + k));
			}
		}
	}
}




/*Non linear operator, containing the potential for Real Time Progression.
Input is pointers to the wavefunction, tilted and un-tilted potential, the step, time and x,y,z grid vectors.*/
void Nr(fftw_complex *w, double *v1, double *v2, int step, double s, double x[], double y[], double z[])
{
	int i, j, k;
	double T, tt, Xt, Yt, Zt, Tt, ct, st, L2;	
	
	//Iterating in Real Time prior to the tilt.
	if (step >= Sground + Svortex && step < Sground + Svortex + Sevolve1) {
		for(i = 0; i < Zpoints; i++) {
			for(j = 0; j < Ypoints; j++) {
				for(k = 0; k < Xpoints; k++) {

					*(w + i * Ypoints * Xpoints + j * Xpoints + k) = -DELTAt * I * (*(v1 + i * Ypoints * Xpoints + j * Xpoints + k) + C * pow(cabs(*(w + i * Ypoints * Xpoints + j * Xpoints + k)), 2)) * 
					(*(w + i * Ypoints * Xpoints + j * Xpoints + k));
				}
			}
		}
	}
	
	//Iterating in Real Time during the tilt.
	else
	if (step >= Sground + Svortex + Sevolve1 && step < Sground + Svortex + Sevolve1 + Stilt) {

		//Pre-computing time dependend parameters for each step.
		tt = 0;
		
		if(Ttilt == 0)
			Tt = 1;
		else
			Tt = Ttilt;
		
		T = (s - (Tground + Tvortex + Tevolve1)) / (Tt);

		tt = theta * (sin(2 * PI * T) / (2 * PI) - T + 1);
		
		st = sin(tt);
		
		ct = cos(tt);
		
		L2 = L * L;
		
		for(i = 0; i < Zpoints; i++) {
			for(j = 0; j < Ypoints; j++) {
				for(k = 0; k < Xpoints; k++) {
					
					//Implementing the tilt through rotation matrices, using pre-computed cosines.
					Xt = x[k];
					Yt = y[j] * ct + z[i] * st;
					Zt = -y[j] * st + z[i] * ct;
					
					*(w + i * Ypoints * Xpoints + j * Xpoints + k) = -DELTAt * I * (0.25 * (Xt * Xt + Yt * Yt + L2 * Zt * Zt) + C * pow(cabs(*(w + i * Ypoints * Xpoints + j * Xpoints + k)), 2)) * 
					(*(w + i * Ypoints * Xpoints + j * Xpoints + k));
				}
			}
		}
	}
	
	//Iterating in real time after completion of the tilt.
	else {
		for(i = 0; i < Zpoints; i++) {
			for(j = 0; j < Ypoints; j++) {
				for(k = 0; k < Xpoints; k++) {
					
					*(w + i * Ypoints * Xpoints + j * Xpoints + k) = -DELTAt * I * (*(v2 + i * Ypoints * Xpoints + j * Xpoints + k) + C * pow(cabs(*(w + i * Ypoints * Xpoints + j * Xpoints + k)), 2)) * 
					(*(w + i * Ypoints * Xpoints + j * Xpoints + k));
				}					
			}
		}
	}
}




/*Function for imprinting a vortex.
Input is pointers to the wavefunction the x,y,a grid vectors and the chosen x0,y0 position of the vortex to be imprinted.*/
void vortex(fftw_complex *q, double x[], double x0, double y[], double y0, double z[], int sign)
{
	int i, j, k;
	double X, Y;
	
	for(i = 0; i < Ypoints; i++) {
		for(j = 0; j < Xpoints; j++) {
			for(k = 0; k < Zpoints; k++) {
				
				//Incoroporating the initial tilt.
				X = x[j];
				Y = y[i] * cos(theta) + z[k] * sin(theta);
				
				//Applying the phase shift to imprint the vortex.
				*(q + k * Ypoints * Xpoints + i * Xpoints + j) = 
				(*(q + k * Ypoints * Xpoints + i * Xpoints + j)) * cexp(I * carg((X - x0) + sign * I * (Y - y0 * cos(theta))));
			}
		}
	}	
}


