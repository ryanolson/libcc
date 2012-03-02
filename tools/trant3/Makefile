# Triples makefile

PROGRAM=trant3

FC = ftn
FCFLAGS_PGI = -i8 -fastsse -Munroll -tp shanghai-64 -Mfpmisalign
FCFLAGS_CCE = -sdefault64 -target=linux -rm -rx -O2,ipa5

FCFLAGS = $(FCFLAGS_CCE)

OBJ = trant3.o

$(PROGRAM).x : $(OBJ) 
	$(FC) $(FCFLAGS) -o $(PROGRAM).x $(OBJ) 


%.o : %.f $(MOD) Makefile
	$(FC) $(FCFLAGS) -c $*.f

%.o : %.f90 $(MOD) Makefile
	$(FC) $(FCFLAGS) -c $*.f90

clean:
	@echo Cleaning ...
	@rm -f *.o
	@rm -f *.mod
	@rm -f *.lst
	@rm -f *.x
