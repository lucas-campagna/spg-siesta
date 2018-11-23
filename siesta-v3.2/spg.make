OBJS += spglib_f08.o\
        siesta2bt.o

COMP_LIBS+= spg-siesta/spglib/_build/libsymspg.a

LDFLAGS += -fopenmp -lm


siesta2bt.o: siesta2bt.f90 spglib_f08.o fdf/fdf_mod.o precision.o siesta_geom.o sys.o atomlist.o m_spin.o sparse_matrices.o siesta_options.o m_energies.o bands.o reclat.o


spglig_f08.o: spglib_f08.f90

