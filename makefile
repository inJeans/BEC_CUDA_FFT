CUDA_INC = -I/usr/local/cuda-5.5/include
CUDA_LIB = -L/usr/local/cuda-5.5/lib
FFTW_INC = -I/usr/local/include
FFTW_LIB = -L/usr/local/lib

all: BEC_groundstate

BEC_groundstate: BEC_evolve.o initial_cond.o operators.o matrix_functions.o differentiate.o
	gcc -std=c99 -o $@ $(FFTW_INC) $(CUDA_INC) $(FFTW_LIB) $(CUDA_LIB) BEC_evolve.o initial_cond.o operators.o matrix_functions.o differentiate.o -lfftw3 -lcufft -lcufftw -lm

BEC_evolve.o: BEC_evolve.c 
	gcc -std=c99 -c $(FFTW_INC) BEC_evolve.c

initial_cond.o: initial_cond.c
	gcc -std=c99 -c $(FFTW_INC) initial_cond.c

operators.o: operators.c
	gcc -std=c99 -c $(FFTW_INC) operators.c

matrix_functions.o: matrix_functions.c
	gcc -std=c99 -c $(FFTW_INC) matrix_functions.c

differentiate.o: differentiate.c
	gcc -std=c99 -c $(CUDA_INC) differentiate.c

clean:
	rm -rf *.o BEC_groundstate

#BEC_groudstate: BEC_evolve.c initial_cond.c operators.c matrix_functions.c
#	gcc -std=c99 -o BEC_evolve -I/usr/local/cuda-5.5/include -L/usr/local/cuda-5.5/lib BEC_evolve.c initial_cond.c operators.c matrix_functions.c -lcufft -lcufftw -lm
