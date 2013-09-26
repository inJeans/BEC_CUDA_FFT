/* Main body of BEC evolution algorithm. */

//Included libraries.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>									//This library is declared before fftw3.h
#include <fftw3.h>										//The FFTW library
//#include <cufftw.h>

//Designated header file with input parameters.
#include "BEC_3D.h"



int main()
{
        cudaMalloc( (void**) &d_in, sizeof(fftw_complex) * Xpoints * Ypoints * Zpoints );
	cudaMalloc( (void**) &d_c,  sizeof(fftw_complex) * Xpoints * Ypoints * Zpoints );
/*-----------------------------Declaring variables & setting the initial conditions-----------------------------*/
	
	int i, j, k;
	double *mod, *out, *out2, *ph, xgrid[Xpoints] = {}, ygrid[Ypoints] = {}, zgrid[Zpoints] = {}; 
	double *tgrid, *V1, *V2;
	double t; //Can we make t just a real value?
	fftw_complex *psi, *psiK, *psiI, *dconsti, *dconstr;	

	//Allocating memory for the wavefunction, potential matrices and time step grid.
	psi = fftw_malloc(sizeof(fftw_complex) * Xpoints * Ypoints * Zpoints);
	tgrid = malloc(sizeof(double) * (Stotal));
	V1 = malloc(sizeof(double) * Xpoints * Ypoints * Zpoints);
	V2 = malloc(sizeof(double) * Xpoints * Ypoints * Zpoints);
	
	//Running the initial condition code to compute the initial condition, level and tilted potential and fill the x, y, z grid vectors.
	initial_cond(psi, V1, V2, xgrid, ygrid, zgrid);
	
	//Filling the time step grid.
	for (j = 0; j < Stotal; j++)	{
		*(tgrid + j) = T_zero + j * DELTAt;
	}
	

	/*Printing the plotting data and initial condition to file.
	 There are 3 files. One prints the wavefunction modulus integrated through Z.
	 One prints the wavefunction modulus integrated through Y.
	 One prints the phase map of the XY plane through the Z-origin.*/
	FILE *outputXY;
	FILE *outputYZ;
	FILE *phaseout;
	outputXY = fopen("outputXY.txt", "w");
	outputYZ = fopen("outputYZ.txt", "w");
	phaseout = fopen("phaseXY.txt", "w");
	
	{
		//Printing integer constants to each file.
		fprintf(outputXY, "%d %d %d %d %d %d %f\n%d\n%d\n%d\n", Sprint, Sground/printer, Svortex/printer, Sevolve1/printer, Stilt/printer, Sevolve2/printer, printer*DELTAt, Xpoints, Ypoints, Zpoints);
		fprintf(outputYZ, "%d %d %d %d %d %d %f\n%d\n%d\n%d\n", Sprint, Sground/printer, Svortex/printer, Sevolve1/printer, Stilt/printer, Sevolve2/printer, printer*DELTAt, Xpoints, Ypoints, Zpoints);
		fprintf(phaseout, "%d %d %d %d %d %d %f\n%d\n%d\n%d\n", Sprint, Sground/printer, Svortex/printer, Sevolve1/printer, Stilt/printer, Sevolve2/printer, printer*DELTAt, Xpoints, Ypoints, Zpoints);
	
		//Printing the xgrid to each file.
		for (i = 0; i < Xpoints; i++)
		{
			fprintf(outputXY, "%lf ", xgrid[i]);
		}	
		fprintf(outputXY, "\n");
	
		for (i = 0; i < Xpoints; i++)
		{
			fprintf(outputYZ, "%lf ", xgrid[i]);
		}	
		fprintf(outputYZ, "\n");
	
		for (i = 0; i < Xpoints; i++)
		{
			fprintf(phaseout, "%lf ", xgrid[i]);
		}	
		fprintf(phaseout, "\n");
	
	
		//Printing the ygrid to each file.
		for (j = 0; j < Ypoints; j++)
		{
			fprintf(outputXY, "%lf ", ygrid[j]);
		}	
		fprintf(outputXY, "\n");
	
		for (j = 0; j < Ypoints; j++)
		{
			fprintf(outputYZ, "%lf ", ygrid[j]);
		}	
		fprintf(outputYZ, "\n");
	
		for (j = 0; j < Ypoints; j++)
		{
			fprintf(phaseout, "%lf ", ygrid[j]);
		}	
		fprintf(phaseout, "\n");
	
	
		//Printing the zgrid to each file.
		for (k = 0; k < Zpoints; k++)
		{
			fprintf(outputXY, "%lf ", zgrid[k]);
		}	
		fprintf(outputXY, "\n");	
	
		for (k = 0; k < Zpoints; k++)
		{
		fprintf(outputYZ, "%lf ", zgrid[k]);
		}	
		fprintf(outputYZ, "\n");
	
		for (k = 0; k < Zpoints; k++)
		{
			fprintf(phaseout, "%lf ", zgrid[k]);
		}	
		fprintf(phaseout, "\n");

		mod = modulus3(psi);
		out = integrateZ(mod);
	
	
		//Printing the initial condition to each file.
		for (i = 0; i < Ypoints; i++)
		{
			for (j = 0; j < Xpoints; j++)
			{
				fprintf(outputXY, "%lf ", *(out + i * Xpoints + j)); //Write the data to a temporary file
			}
			fprintf(outputXY, "\n");
		}	
		fprintf(outputXY, "\n");
	
		free(out);
	

		out2 = integrateX(mod);
	
		for (i = 0; i < Zpoints; i++)
		{
			for (j = 0; j < Ypoints; j++)
			{
				fprintf(outputYZ, "%lf ", *(out2 + i * Ypoints + j)); //Write the data to a temporary file
			}
			fprintf(outputYZ, "\n");
		}	
		fprintf(outputYZ, "\n");
	
		free(out2);
		free(mod);
		
		ph = phase(psi);
	
		for (i = 0; i < Ypoints; i++)
		{
			for (j = 0; j < Xpoints; j++)
			{
				fprintf(phaseout, "%lf ", *(ph + i * Xpoints + j)); //Write the data to a temporary file
			}
			fprintf(phaseout, "\n");
		}	
		fprintf(phaseout, "\n");
	
		free(ph);
	}

/*----------------------------Running the efficient Runge Kutta 4th order method in ITP-----------------------------------*/

	//Filling the diffusion coefficients for Imaginary Time Progression.
	dconsti = fftw_malloc(sizeof(fftw_complex) * Xpoints * Ypoints * Zpoints);
	DCi(dconsti);

	//Progressing in Imaginary Time
	for (k = 0; k < Sground + Svortex; k++) 
	{
		
		//Imprinting vortices after finding the ground state.
		if (k == Sground + 1)
		{
			vortex(psi, xgrid, 0, ygrid, 0, zgrid, 1);
			vortex(psi, xgrid, -6, ygrid, 1.5, zgrid, -1);
			vortex(psi, xgrid, -2, ygrid, 4, zgrid, 1);
			vortex(psi, xgrid, 3, ygrid, -6, zgrid, -1);
			vortex(psi, xgrid, -2, ygrid, -4, zgrid, 1);
			
			vortex(psi, xgrid, 3, ygrid, 5.666, zgrid, -1);
			vortex(psi, xgrid, 3.5, ygrid, 0, zgrid, 1);
			vortex(psi, xgrid, 7, ygrid, 1, zgrid, -1);
			vortex(psi, xgrid, -3.5, ygrid, -3, zgrid, 1);
			vortex(psi, xgrid, 4.5, ygrid, -3, zgrid, 1);
		}
		
		psiK = copy(psi);
			
		D(psi, dconsti);
		
		psiI = copy(psi);
		
		Ni(psiK, V1, xgrid, ygrid, zgrid);
		
		D(psiK, dconsti);
	
		weight_sum(1, psi, 0.1666666667, psiK);
		
		weight_sum(0.5, psiK, 1, psiI);
		
		Ni(psiK, V1, xgrid, ygrid, zgrid);
		
		weight_sum(1, psi, 0.3333333333, psiK);
		
		weight_sum(0.5, psiK, 1, psiI);
		
		Ni(psiK, V1, xgrid, ygrid, zgrid);
		
		weight_sum(1, psi, 0.3333333333, psiK);
		
		weight_sum(1, psiK, 1, psiI);
		fftw_free(psiI);
		
		D(psiK, dconsti);
		
		D(psi, dconsti);
		
		Ni(psiK, V1, xgrid, ygrid, zgrid);
		
		weight_sum(1, psi, 0.1666666667, psiK);
		fftw_free(psiK);
		
		//Normalisation
		norm(psi);
		
		//The algorithm now iterates and the data printed at set increments.
		if (k % printer == 0)
		{	
			//Taking the modulus and integrating through the Z axis.	
			mod = modulus3(psi);
			out = integrateZ(mod);
		
			//Printing the data to a text file.
			for (i = 0; i < Ypoints; i++)
			{
				for (j = 0; j < Xpoints; j++)
				{
					fprintf(outputXY, "%lf ", *(out + i * Xpoints + j));
				}
				fprintf(outputXY, "\n");
			}
		
			fprintf(outputXY, "\n");		
			free(out);		
			
			//Integrating through the X-axis.
			out2 = integrateX(mod);
		
			//Printing the data to a text file.
			for (i = 0; i < Zpoints; i++)
			{
				for (j = 0; j < Ypoints; j++)
				{
					fprintf(outputYZ, "%lf ", *(out2 + i * Ypoints + j));
				}
				fprintf(outputYZ, "\n");
			}
		
			fprintf(outputYZ, "\n");
			
			free(out2);
			free(mod);
			
			//Computing the phase.
			ph = phase(psi);
			
			//Printing the data to a text file.
			for (i = 0; i < Ypoints; i++)
			{
				for (j = 0; j < Xpoints; j++)
				{
					fprintf(phaseout, "%lf ", *(ph + i * Xpoints + j)); //Write the data to a temporary file
				}
				fprintf(phaseout, "\n");
			}	
			fprintf(phaseout, "\n");
	
			free(ph);		
		}	

		//Printing the step number.
		printf("%d\n", k);
	}
	
/*----------------------------Running the efficient Runge Kutta 4th order method in RTP-----------------------------------*/	
	
	//Freeing the imaginary diffusion coefficients and filling the coefficients for Real Time Progression.
	free(dconsti);
	dconstr = fftw_malloc(sizeof(fftw_complex) * Xpoints * Ypoints * Zpoints);
	DCr(dconstr);	
	
	
	//Progressing in Real Time
	for (k = Sground + Svortex; k < Stotal; k++)
	{

		t = tgrid[k];
		
		psiK = copy(psi);

		D(psi, dconstr);
		
		psiI = copy(psi);
		
		Nr(psiK, V1, V2, k, t, xgrid, ygrid, zgrid);
		
		D(psiK, dconstr);
		
		weight_sum(1, psi, 0.1666666667, psiK);
		
		weight_sum(0.5, psiK, 1, psiI);
		
		t = t + DELTAt * 0.5;
		
		Nr(psiK, V1, V2, k, t, xgrid, ygrid, zgrid);
		
		weight_sum(1, psi, 0.3333333333, psiK);
		
		weight_sum(0.5, psiK, 1, psiI);
		
		Nr(psiK, V1, V2, k, t, xgrid, ygrid, zgrid);
		
		weight_sum(1, psi, 0.3333333333, psiK);
		
		weight_sum(1, psiK, 1, psiI);
		fftw_free(psiI);
		
		D(psiK, dconstr);
			
		D(psi, dconstr);
		
		t = t + DELTAt * 0.5;
		
		Nr(psiK, V1, V2, k, t, xgrid, ygrid, zgrid);
		
		weight_sum(1, psi, 0.1666666667, psiK);
		fftw_free(psiK);
		
		//Normalisation is not required in RTP.
		
		//The algorithm now iterates and prints for set steps. Same as above.
		if (k % printer == 0)
		{
			mod = modulus3(psi);
			out = integrateZ(mod);

			for (i = 0; i < Ypoints; i++)
			{
				for (j = 0; j < Xpoints; j++)
				{
					fprintf(outputXY, "%lf ", *(out + i * Xpoints + j));
				}
				fprintf(outputXY, "\n");
			}
		
			fprintf(outputXY, "\n");		
			free(out);		
		
		
			out2 = integrateX(mod);
					
			for (i = 0; i < Zpoints; i++)
			{
				for (j = 0; j < Ypoints; j++)
				{
					fprintf(outputYZ, "%lf ", *(out2 + i * Ypoints + j));
				}
				fprintf(outputYZ, "\n");
			}
		
			fprintf(outputYZ, "\n");
			
			free(out2);
			free(mod);
				
			ph = phase(psi);
	
			for (i = 0; i < Ypoints; i++)
			{
				for (j = 0; j < Xpoints; j++)
				{
					fprintf(phaseout, "%lf ", *(ph + i * Xpoints + j)); 
				}
				fprintf(phaseout, "\n");
			}	
			fprintf(phaseout, "\n");
	
			free(ph);		
		}	
		
			//Printing the step.
			printf("%d\n", k);		
	}	
	

	//Closing the output files.
	fclose(outputXY);
	fclose(outputYZ);
	fclose(phaseout);

	cudaFree( d_in );
	cudaFree( d_c  );

	return 0;	
}


