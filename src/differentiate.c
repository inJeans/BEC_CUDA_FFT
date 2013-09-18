/* File contains key operators including the linear (diffusion) term and non-linear (inc. potential)*/

//Included libraries.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
//#include <fftw3.h>
#include <cufftw.h>

//Designated header file with input parameters.
#include "BEC_3D.h"



/*Linear (diffusion) operator. Utilised throughout the algorithm.
Function input is the pointers to the wavefunction and coefficient matrix.
This is where I assume you would implement the CUDA FFT algorithm*/
void D(fftw_complex *in, fftw_complex *c)
{
	fftw_plan plan, plan2;
	int l, m, n;
	
	//Setting up the backward and forward transforms.
	plan = fftw_plan_dft_3d(Zpoints, Ypoints, Xpoints, in, in, FFTW_BACKWARD, FFTW_MEASURE);
	plan2 = fftw_plan_dft_3d(Zpoints, Ypoints, Xpoints, in, in, FFTW_FORWARD, FFTW_MEASURE);
		
	//Executing the backward transform.
	fftw_execute(plan);
	
	//Multiplying the transformed wavefunction by the previously calculated coefficients (depends on if we are in ITP or RTP).		
	for(n = 0; n < Zpoints; n++) {
			
		for(m = 0; m < Ypoints; m++) {
				
			for(l = 0; l < Xpoints; l++) {
					
			  //*(in + n * Ypoints * Xpoints + m * Xpoints + l) = *(in + n * Ypoints * Xpoints + m * Xpoints + l) * (*(c + n * Ypoints * Xpoints + m * Xpoints + l));
			  *(in + n * Ypoints * Xpoints + m * Xpoints + l)[0] = *(in + n * Ypoints * Xpoints + m * Xpoints + l)[0] * (*(c + n * Ypoints * Xpoints + m * Xpoints + l)[0]) - *(in + n * Ypoints * Xpoints + m * Xpoints + l)[1] * (*(c + n * Ypoints * Xpoints + m * Xpoints + l)[1]);
			  *(in + n * Ypoints * Xpoints + m * Xpoints + l)[1] = *(in + n * Ypoints * Xpoints + m * Xpoints + l)[0] * (*(c + n * Ypoints * Xpoints + m * Xpoints + l)[1]) + *(in + n * Ypoints * Xpoints + m * Xpoints + l)[1] * (*(c + n * Ypoints * Xpoints + m * Xpoints + l)[0]);
			}
		}			
	}		
	
	//Executing the forward transform.
	fftw_execute(plan2);

	fftw_destroy_plan(plan);
	fftw_destroy_plan(plan2);
}
