# Triples makefile

FC = ftn
FCFLAGS_PGI = -i8 -fastsse -Munroll -tp shanghai-64 -Mfpmisalign
FCFLAGS_CCE = -sdefault64 -target=linux -rm -rd -O2,ipa5

FCFLAGS = $(FCFLAGS_CCE)

MOD = common_cc.o
OBJ = triples.o ddcc_aux.o ddcc_t_ijj_iij_big.o ddcc_t_ijk_big.o tranmd.o trant3.o t3squa.o
OBJ90 = ijk_tuple.o ijj_tuple.o iij_tuple.o # t3squa-cache.o

triples.x : $(MOD) $(OBJ) $(OBJ90) Makefile
	$(FC) $(FCFLAGS) -o triples.x $(MOD) $(OBJ) $(OBJ90) -L./ -lddi

%.o : %.f triples.f90 Makefile
	$(FC) $(FCFLAGS) -c $*.f

# triples.o : triples.f90

%.o : %.f90 triples.f90 Makefile
	$(FC) $(FCFLAGS) -c $*.f90

clean:
	@echo Cleaning ...
	@rm -f *.o
	@rm -f *.mod
	@rm -f *.lst
	@rm -f *.x
