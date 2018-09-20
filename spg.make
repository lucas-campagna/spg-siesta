OBJS += spglib_f08.o\
        siesta2bt.o

COMP_LIBS+= spg-siesta/libsymspg.a

LDFLAGS += -fopenmp -lm


siesta2bt.o: siesta2bt.f90 spglib_f08.o


spglig_f08.o: spglib_f08.f90


	



