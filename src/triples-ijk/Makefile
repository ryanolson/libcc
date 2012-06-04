# Triples makefile

PROGRAM=triples
PROGRAM_CUDA=triples-cuda
PROGRAM_ACC=triples-acc

FC = ftn
FCFLAGS_PGI = -i8 -fastsse -Munroll -tp shanghai-64 -Mfpmisalign
FCFLAGS_CCE = -sdefault64 -target=linux -rm -rd -O2,ipa5 -Onoomp
FCFLAGS_CCE = -sdefault64 -rm -rd -O2,ipa5 -Onoomp

NVCC = nvcc
NVINCLUDE = -I../ddi-onesided/include
NVFLAGS = -arch sm_20 -Xptxas=-v -O2
CUDA_LIBS = -lcublas -lcudart

FCFLAGS = $(FCFLAGS_CCE) -DGPU

MOD = common_cc.o
OBJ = ddcc_aux.o t1wt3.o t3squa.o triples.o ddcc_t_ijj_iij_big.o
OBJ90 = ijj_tuple.o iij_tuple.o 
OBJCPU = ijk_tuple.o ddcc_t_ijk_big.o
OBJACC = ijk_tuple_acc.o ddcc_t_ijk_big.o
OBJCUDA = cuda_ijk_kernels.o \
	      cuda_t1wt3_kernels.o \
	      cuda_ijk_driver.o \
		  cuda_triples_driver.o \
	      $(OBJCPU)

all: $(PROGRAM).x 

$(PROGRAM).x : $(MOD) $(OBJ) $(OBJ90) $(OBJCPU) Makefile
	$(FC) $(FCFLAGS) -Onoacc -o $(PROGRAM).x $(MOD) $(OBJ) $(OBJ90) $(OBJCPU) -L../ddi-onesided -lddi

$(PROGRAM_ACC).x : $(MOD) $(OBJ) $(OBJ90) $(OBJACC) Makefile
	$(FC) $(FCFLAGS) -o $(PROGRAM_ACC).x $(MOD) $(OBJ) $(OBJ90) $(OBJACC) -L../ddi-onesided -lddi

$(PROGRAM_CUDA).x : $(MOD) $(OBJ) $(OBJ90) $(OBJCUDA) Makefile
	$(FC) $(FCFLAGS) -Onoacc -o $(PROGRAM_CUDA).x $(MOD) $(OBJ) $(OBJ90) $(OBJCUDA) -L../ddi-onesided -lddi $(CUDA_LIBS)

$(PROGRAM).x-pat : triples.x
	@make clean
	@make $(PROGRAM).x
	@pat_build -f -u $(PROGRAM).x $(PROGRAM).x-pat

%.o : %.f common_cc.f90 Makefile
	$(FC) $(FCFLAGS) -c $*.f

%.o : %.cu cuda_kernels.h
	$(NVCC) $(NVINCLUDE) $(NVFLAGS) -c $*.cu

# triples.o : triples.f90

%.o : %.f90 common_cc.f90 Makefile
	$(FC) $(FCFLAGS) -c $*.f90

%.o : %.F90 common_cc.f90 Makefile
	$(FC) $(FCFLAGS) -c $*.F90

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
