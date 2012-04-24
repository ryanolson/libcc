  GA_CONFIG := /lus/scratch/ryan/local/bin/ga-config
ifneq ($(GA_CONFIG),)
CPPFLAGS := $(shell $(GA_CONFIG) --cppflags)
  FFLAGS := $(shell $(GA_CONFIG) --fflags)
  CFLAGS := $(shell $(GA_CONFIG) --cflags)
    LIBS := $(shell $(GA_CONFIG) --libs)
   FLIBS := $(shell $(GA_CONFIG) --flibs)
 LDFLAGS := $(shell $(GA_CONFIG) --ldflags)
      CC := $(shell $(GA_CONFIG) --cc)
     F77 := $(shell $(GA_CONFIG) --f77)
endif

DEPS = ddi.h ddi_armci.h
OBJS = ddi_armci.o ddi_armci_init.o ddi_armci_get.o ddi_armci_put.o ddi_armci_acc.o ddi_armci_getacc.o \
       ddi.o ddi_comm_bcast.o ddi_comm_create.o ddi_comm_group.o ddi_comm_gsum.o ddi_comm_nnode.o      \
	   ddi_comm_nproc.o ddi_comm_recv.o ddi_comm_send.o ddi_comm_sync.o ddi_get.o ddi_put.o ddi_acc.o  \
	   ddi_getacc.o ddi_distrib.o ddi_subpatch.o ddi_index.o ddi_create.o ddi_destroy.o ddi_finalize.o \
	   ddi_memory.o ddi_isend.o ddi_irecv.o ddi_wait.o ddi_signals.o ddi_id.o ddi_timer.o ddi_gdlb.o   \
	   ddi_dlb.o ddi_dlb_proc.o ddi_init.o ddi_runtime.o ddi_util.o mmath.o ddi_arr.o ddi_arr_acc.o    \
	   ddi_arr_add.o ddi_arr_dot.o ddi_arr_scalar.o ddi_arr_select.o ddi_comm.o ddi_smp.o std_system.o \
	   ddi_send.o debug.o ddi_fortran.o

CPPFLAGS += -DLINUX64 -DDDI_MPI -DDDI_ARMCI -DUSE_MPI_BARRIER -DMAX_NODES=4096 -DMAX_SMP_PROCS=1 -I./
CFLAGS = " "

libddi.a: tags
	ar cr libddi.a $(OBJS)


tags: $(OBJS) $(DEPS)
	ctags *


%.o : %.c $(DEPS)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $*.c


ddi_fortran.o : ddi_fortran.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -DINT_SIZE=long -D_UNDERSCORES=1  -c $*.c

clean:
	@echo "Deleting *.o *.a ..."
	@rm -rf *.o *.a
