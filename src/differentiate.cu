/* File contains key operators including the linear (diffusion) term and non-linear (inc. potential)*/

//Included libraries.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <complex.h>
//#include <fftw3.h>
#include <cufftw.h>

//Designated header file with input parameters.
extern "C" {
#include "BEC_3D.h"
}

__global__ void multcoeff( fftw_complex *d_in , fftw_complex *d_c) {
  
    int idx = threadIdx.x + blockIdx.x*gridDim.x;

    if ( idx > Zpoints ) {
	return;
    }

    for(int m = 0; m < Ypoints; m++) {	
			for(int l = 0; l < Xpoints; l++) {

			  *(d_in + idx * Ypoints * Xpoints + m * Xpoints + l)[0] = *(d_in + idx * Ypoints * Xpoints + m * Xpoints + l)[0] * (*(d_c + idx * Ypoints * Xpoints + m * Xpoints + l)[0]) - *(d_in + idx * Ypoints * Xpoints + m * Xpoints + l)[1] * (*(d_c + idx * Ypoints * Xpoints + m * Xpoints + l)[1]);
			  *(d_in + idx * Ypoints * Xpoints + m * Xpoints + l)[1] = *(d_in + idx * Ypoints * Xpoints + m * Xpoints + l)[0] * (*(d_c + idx * Ypoints * Xpoints + m * Xpoints + l)[1]) + *(d_in + idx * Ypoints * Xpoints + m * Xpoints + l)[1] * (*(d_c + idx * Ypoints * Xpoints + m * Xpoints + l)[0]);
			}
		}
      
    return;
}


/*Linear (diffusion) operator. Utilised throughout the algorithm.
Function input is the pointers to the wavefunction and coefficient matrix.
This is where I assume you would implement the CUDA FFT algorithm*/
extern "C" void D(fftw_complex *in, fftw_complex *c)
{
	cufftHandle plan;
	
	cudaMemcpy( d_in, in, sizeof(fftw_complex) * Xpoints * Ypoints * Zpoints, cudaMemcpyHostToDevice );
	cudaMemcpy( d_c , c , sizeof(fftw_complex) * Xpoints * Ypoints * Zpoints, cudaMemcpyHostToDevice );

	//Setting up the backward and forward transforms.
	//plan = fftw_plan_dft_3d(Zpoints, Ypoints, Xpoints, in, in, FFTW_BACKWARD, FFTW_MEASURE);
	//plan2 = fftw_plan_dft_3d(Zpoints, Ypoints, Xpoints, in, in, FFTW_FORWARD, FFTW_MEASURE);	
	
	cufftPlan3d( &plan, Xpoints, Ypoints, Zpoints, CUFFT_Z2Z );

	//Executing the backward transform.
	//fftw_execute(plan);
	cufftExecZ2Z( plan, (cufftDoubleComplex*) d_in, (cufftDoubleComplex*) d_in, CUFFT_INVERSE );

	//Multiplying the transformed wavefunction by the previously calculated coefficients (depends on if we are in ITP or RTP).		
        multcoeff <<<(Zpoints + 15)/16, 16>>> ( d_in, d_c );		
	
	//Executing the forward transform.
	//fftw_execute(plan2);
	cufftExecZ2Z( plan, (cufftDoubleComplex*) d_in, (cufftDoubleComplex*) d_in, CUFFT_FORWARD );

	//fftw_destroy_plan(plan);
	//fftw_destroy_plan(plan2);
	cufftDestroy( plan );

	cudaMemcpy( in, d_in, sizeof(fftw_complex) * Xpoints * Ypoints * Zpoints, cudaMemcpyDeviceToHost );
	cudaMemcpy( c , d_c , sizeof(fftw_complex) * Xpoints * Ypoints * Zpoints, cudaMemcpyDeviceToHost );
}
