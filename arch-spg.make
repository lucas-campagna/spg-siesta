.PHONY: clean-spg donwload-spg install-spg

OBJS += spglib_f08.o\
        siesta2bt.o

LIBS += -fopenmp -lm

SPGPATH=$(VPATH)/spg

install-spg: download-spg
	export LIBS='-lgomp'
	export CFLAGS='-fopenmp'
        ( cd $(SPGPATH)/spglib ; mkdir _build ; cd _build ; cmake .. ; make )

lib-spg: lib libsymspg.a

siesta2bt.o: spglib_f08.o

libsymspg.a: $(SPGPATH)/spglib $(SPGPATH)/siesta2bt.f90
        @ ( cd $(SPGPATH) ; cp siesta2bt.f90 .. ; spglib/example/spglib_f08.f90 .. ) 
        cp $</_build/libsymspg.a .

$(SPGPATH)/siesta.f90: $(SPGPATH)/siesta2bt_spg_findsym.f90
	@ ( cd $(SPGPATH) ; cp siesta2bt_spg_findsym.f90 siesta2bt.f90 )


donwload-spg: 
        @ ( cd $(SPGPATH) ; test -d spglib || git clone https://github.com/atztogo/spglib.git )

clean-spg:
        rm -fr $(SPGPATH)/spglib

$(SPGPATH)/spglib: install-spg 

$(SPGPATH)/siesta2bt_%.f90:
	git clone git@github.com:lucas-campagna/spg-siesta.git
		
