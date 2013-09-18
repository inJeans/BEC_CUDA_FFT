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
OBJDIR   = obj/
SRCDIR   = src/

FFTW_INC = -I/usr/local/include
FFTW_LIB = -L/usr/local/lib

INCLUDE = -I include/

all: BEC_groundstate

BEC_groundstate: BEC_evolve.o initial_cond.o operators.o matrix_functions.o differentiate.o
	gcc -std=c99 -o $(addprefix $(BUILDDIR), $@) $(INCLUDE) $(FFTW_INC) $(FFTW_LIB) $(CUDA_LIB) $(addprefix $(BUILDDIR)$(OBJDIR), $?) -lcufft -lcufftw -lm
#	gcc -std=c99 -o $(addprefix $(BUILDDIR), $@) $(INCLUDE) $(FFTW_INC) $(FFTW_LIB) $(CUDA_LIB) $(addprefix $(BUILDDIR)$(OBJDIR), $?) -lfftw3 -lm

BEC_evolve.o: $(addprefix $(SRCDIR), BEC_evolve.c ) 
	gcc -std=c99 -c $(INCLUDE) $(FFTW_INC) $(addprefix $(SRCDIR), BEC_evolve.c ) -o $(addprefix $(BUILDDIR)$(OBJDIR), $@)

initial_cond.o: $(addprefix $(SRCDIR), initial_cond.c )
	gcc -std=c99 -c $(INCLUDE) $(FFTW_INC) $(addprefix $(SRCDIR), initial_cond.c ) -o $(addprefix $(BUILDDIR)$(OBJDIR), $@)

operators.o: $(addprefix $(SRCDIR), operators.c )
	gcc -std=c99 -c $(INCLUDE) $(FFTW_INC) $(addprefix $(SRCDIR), operators.c ) -o $(addprefix $(BUILDDIR)$(OBJDIR), $@)

matrix_functions.o: $(addprefix $(SRCDIR), matrix_functions.c )
	gcc -std=c99 -c $(INCLUDE) $(FFTW_INC) $(addprefix $(SRCDIR), matrix_functions.c ) -o $(addprefix $(BUILDDIR)$(OBJDIR), $@)

differentiate.o: $(addprefix $(SRCDIR), differentiate.c )
	gcc -std=c99 -c $(INCLUDE) $(FFTW_INC) $(addprefix $(SRCDIR), differentiate.c ) -o $(addprefix $(BUILDDIR)$(OBJDIR), $@)

clean:
	rm -rf $(addprefix $(BUILDDIR)$(OBJDIR), *.o) $(addprefix $(BUILDDIR), BEC_groundstate)

#BEC_groudstate: BEC_evolve.c initial_cond.c operators.c matrix_functions.c
#	gcc -std=c99 -o BEC_evolve -I/usr/local/cuda-5.5/include -L/usr/local/cuda-5.5/lib BEC_evolve.c initial_cond.c operators.c matrix_functions.c -lcufft -lcufftw -lm
