
OBJS = transforms.o ccaux_subset.o

t.x : unit_tests.F90 $(OBJS)
	ftn -o t.x unit_tests.F90 $(OBJS)

%.o : %.F90
	ftn -c $*.F90

%.o : %.f
	ftn -c $*.f

clean :
	@echo "Deleting *.o *.x"
	@rm -f *.o
	@rm -f *.x
