# Triples makefile

PROGRAM=sd

CC = cc
CPPFLAGS = -I../ddi-onesided/src -I../ddi-onesided/include
CFLAGS = "-O2"

NVCC = nvcc
NVFLAGS = -arch sm_20 -Xptxas=-v
CUDA_LIBS = -lcublas -lcudart

FCFLAGS = $(FCFLAGS_CCE) -DGPU

OBJ = streaming_dgemm_ddi.o

all: $(PROGRAM).x 

$(PROGRAM).x : $(OBJ) Makefile
	$(CC) -o $(PROGRAM).x $(OBJ) -L../ddi-onesided -lddi $(CUDA_LIBS)

$(PROGRAM).x-pat : triples.x
	@make clean
	@make $(PROGRAM).x
	@pat_build -f -u $(PROGRAM).x $(PROGRAM).x-pat

%.o : %.c Makefile
	$(CC) $(CPPFLAGS) $(CCFLAGS) -c $*.c

%.o : %.cu Makefile
	$(NVCC) $(CPPFLAGS) $(NVFLAGS) -c $*.cu

remove_objects:
	@echo Removing Object Files ...
	@rm -f *.o
	@rm -f *.mod
	@rm -f *.lst
	@rm -f *.opt
	@rm -f *.cg
	@rm -f *.ptx
	@rm -f *.cub

clean: remove_objects
	@echo Removing Binaries ...
	@rm -f *.x
