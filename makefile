UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin) #If building on an OSX system
    CUDA_INC = -I/Developer/NVIDIA/CUDA-5.5/include
    CUDA_LIB = -L/Developer/NVIDIA/CUDA-5.5/lib
else                     #If building on a Linux system
    CUDA_INC = -I/usr/local/cuda-5.5/include
    UNAME_P := $(shell uname -p)
    ifeq ($(UNAME_P),x86_64)  #If 64 bit
        CUDA_LIB = -L/usr/local/cuda-5.5/lib64
    else
        CUDA_LIB = -L/usr/local/cuda-5.5/lib
    endif
endif

BUILDDIR = bin/
OBJDIR  = obj/

FFTW_INC = -I/usr/local/include
FFTW_LIB = -L/usr/local/lib

all: BEC_groundstate

BEC_groundstate: BEC_evolve.o initial_cond.o operators.o matrix_functions.o differentiate.o
	gcc -std=c99 -o $(addprefix $(BUILDDIR), $@) $(FFTW_INC) $(CUDA_INC) $(FFTW_LIB) $(CUDA_LIB) $(addprefix $(BUILDDIR)$(OBJDIR), $?) -lcufft -lcufftw -lm

BEC_evolve.o: BEC_evolve.c 
	gcc -std=c99 -c $(FFTW_INC) BEC_evolve.c -o $(addprefix $(BUILDDIR)$(OBJDIR), $@)

initial_cond.o: initial_cond.c
	gcc -std=c99 -c $(FFTW_INC) initial_cond.c -o $(addprefix $(BUILDDIR)$(OBJDIR), $@)

operators.o: operators.c
	gcc -std=c99 -c $(FFTW_INC) operators.c -o $(addprefix $(BUILDDIR)$(OBJDIR), $@)

matrix_functions.o: matrix_functions.c
	gcc -std=c99 -c $(FFTW_INC) matrix_functions.c -o $(addprefix $(BUILDDIR)$(OBJDIR), $@)

differentiate.o: differentiate.c
	gcc -std=c99 -c $(CUDA_INC) differentiate.c -o $(addprefix $(BUILDDIR)$(OBJDIR), $@)

clean:
	rm -rf $(addprefix $(BUILDDIR)$(OBJDIR), *.o) $(addprefix $(BUILDDIR), BEC_groundstate)

#BEC_groudstate: BEC_evolve.c initial_cond.c operators.c matrix_functions.c
#	gcc -std=c99 -o BEC_evolve -I/usr/local/cuda-5.5/include -L/usr/local/cuda-5.5/lib BEC_evolve.c initial_cond.c operators.c matrix_functions.c -lcufft -lcufftw -lm
