/*Header file were parameters, grid size and number of time steps can be set for the experiment */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <complex.h>									//This library is declared before fftw3.h
//#include <fftw3.h>										//The FFTW library
//#include <cufftw.h>
#include <cuda.h>
#include <cuda_runtime.h>

#define PI 3.141592654									//Defining pi

#define Xpoints 32										//Number of x increments
#define xSTART - 5 * PI									//Start point for the x grid
#define xEND 5 * PI										//End point for the x grid

#define Ypoints 32										//Number of y increments
#define ySTART - 5 * PI									//Start point for the y grid
#define yEND 5 * PI										//End point for the y grid

#define Zpoints 32										//Number of z increments
#define zSTART - 2 * PI									//Start point for the z grid
#define zEND 2 * PI										//End point for the z grid

//Defining the grid increment sizes from the spacial arrangement above.
#define DELTAx (double) ((xEND - xSTART) / (Xpoints - 1))
#define DELTAy (double) ((yEND - ySTART) / (Ypoints - 1))
#define DELTAz (double) ((zEND - zSTART) / (Zpoints - 1))

#define Tground 0.2										//Time to evolve to ground state in ITP.
#define Tvortex 0										//Time for vortex to develop in ITP.
#define Tevolve1 0.2									//Time to evolve in real time.
#define Ttilt 1											//Time for tilting
#define Tevolve2 0.5									//Time to evolve after tilt
#define T_zero 0										//Start time
#define DELTAt 0.0005									//Time increment

#define C 4660 											//Non-linear constant C
#define L 11.25											//Ratio of harmonic frequencies Lambda
#define theta 20 * PI / 180								//Tilt angle theta

#define printer 10										//Number of steps between prints

//Defining the number of steps for each part of the experiment from the times above.
#define Sground (int) (Tground / DELTAt)
#define Svortex (int) (Tvortex / DELTAt)
#define Sevolve1 (int) (Tevolve1 / DELTAt)
#define Stilt (int) (Ttilt / DELTAt)
#define Sevolve2 (int) (Tevolve2 / DELTAt)
#define Stotal (Sground + Svortex + Sevolve1 + Stilt + Sevolve2)
#define Sprint ((Sground + Svortex + Sevolve1 + Stilt + Sevolve2) / printer)

// Device pointers
fftw_complex *d_in, *d_c;

//Operators. See operators.c for function descriptions.
void initial_cond(fftw_complex *w, double *v1, double *v2, double x[], double y[], double z[]);

void DCi(fftw_complex *p);

void DCr(fftw_complex *p);

void D(fftw_complex *in, fftw_complex *c);

void Ni(fftw_complex *w, double *v, double x[], double y[], double z[]);

void Nr(fftw_complex *w, double *v1, double *v2, int step, double s, double x[], double y[], double z[]);

void vortex(fftw_complex *q, double x[], double x0, double y[], double y0, double z[], int sign);

//Matrix functions. See matrix_functions.c for descriptions.
void weight_sum(double a, fftw_complex *q, double b, fftw_complex *r);

fftw_complex *copy(fftw_complex *q);

double *modulus3(fftw_complex *q);

double *integrateZ(double *q);

double *integrateX(double *q);

void norm(fftw_complex *q);

double *phase(fftw_complex *q);
