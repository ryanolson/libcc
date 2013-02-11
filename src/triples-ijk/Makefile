# Triples makefile

PROGRAM=triples
PROGRAM_CUDA=triples-cuda
PROGRAM_ACC=triples-acc

FC = ftn
FCFLAGS_PGI = -i8 -fastsse -Munroll -tp shanghai-64 -Mfpmisalign
#FCFLAGS_CCE = -sdefault64 -target=linux -rm -Oipa5 -Onoomp
FCFLAGS_CCE = -sdefault64 -rmd -Onoomp

NVCC = nvcc
NVFLAGS = -arch sm_20 -Xptxas=-v
CUDA_LIBS = -lcublas -lcudart

# use for CUDA version
#FCFLAGS = $(FCFLAGS_CCE) -DGPU -DUSE_CUDA -F

# use for OpenACC version
FCFLAGS = $(FCFLAGS_CCE) -DGPU  -DUSE_OPEN_ACC -F -hacc_model=fast_addr 

# use for CPU version
#FCFLAGS = $(FCFLAGS_CCE) -F

MOD = common_cc.o
#OBJ = ddcc_aux.o t3squa.o triples.o ddcc_t_ijj_iij_big.o
OBJ = ddcc_aux.o ddcc_t_ijj_iij_big.o common_functions.o program.o 
#OBJ90 = ijj_tuple.o iij_tuple.o
OBJ90 = 
OBJCPU = ijk_tuple.o ddcc_t_ijk_big.o t1wt3.o triples_cpu.o t3squa.o ijj_tuple.o iij_tuple.o
OBJACC = ijk_tuple_acc.o ddcc_t_ijk_big.o t1wt3_acc.o triples.o iij_tuple_acc.o ijj_tuple_acc.o t3squa_acc.o dgemm_async.o
#OBJACC = ijk_tuple_acc.o ddcc_t_ijk_big.o t1wt3_acc.o triples.o
OBJCUDA = ijk_tuple_cuda_kernels.o triples.o \
	      t1wt3_ijk_cuda_kernels.o \
	      ddcc_t_ijk_big_cuda.o ijk_gpu_driver.o \
	      ijk_tuple.o ddcc_t_ijk_big.o t1wt3.o t3squa.o ijj_tuple.o iij_tuple.o

all: $(PROGRAM_ACC).x 

$(PROGRAM).x : $(MOD) $(OBJ) $(OBJ90) $(OBJCPU) Makefile
	$(FC) $(FCFLAGS) -Onoacc -o $(PROGRAM).x $(MOD) $(OBJ) $(OBJ90) $(OBJCPU) -L../ddi-onesided -lddi

$(PROGRAM_ACC).x : $(MOD) $(OBJ) $(OBJ90) $(OBJACC) Makefile
	$(FC) $(FCFLAGS) -o $(PROGRAM_ACC).x $(MOD) $(OBJ) $(OBJ90) $(OBJACC) -L../ddi-onesided -lddi $(usemylib)

$(PROGRAM_CUDA).x : $(MOD) $(OBJ) $(OBJ90) $(OBJCUDA) Makefile
	$(FC) $(FCFLAGS) -Onoacc -o $(PROGRAM_CUDA).x $(MOD) $(OBJ) $(OBJ90) $(OBJCUDA) -L../ddi-onesided -lddi $(CUDA_LIBS)

$(PROGRAM).x-pat : triples.x
	@make clean
	@make $(PROGRAM).x
	@pat_build -f -u $(PROGRAM).x $(PROGRAM).x-pat

%.o : %.f common_cc.f90 Makefile
	$(FC) $(FCFLAGS) -c $*.f

%.o : %.cu cuda_kernels.h
	$(NVCC) $(NVFLAGS) -c $*.cu

# triples.o : triples.f90

%.o : %.f90 common_cc.f90 Makefile
	$(FC) $(FCFLAGS) -c $*.f90

%.o : %.F90 common_cc.f90 Makefile
	$(FC) $(FCFLAGS) -c $*.F90

%.o : %.cpp Makefile
	CC -hacc_model=fast_addr -c $*.cpp

remove_objects:
	@echo Removing Object Files ...
	@rm -f *.o
	@rm -f *.mod
	@rm -f *.lst
	@rm -f *.opt
	@rm -f *.cg
	@rm -f *.ptx
	@rm -f *.cub
	@rm -f *.cu

clean: remove_objects
	@echo Removing Binaries ...
	@rm -f *.x
