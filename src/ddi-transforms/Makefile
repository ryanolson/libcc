
OBJS = transforms.o

t.x : test.F90 $(OBJS)
	ftn -o t.x test.F90 $(OBJS)

%.o : %.F90
	ftn -c $*.F90
