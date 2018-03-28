.PHONY: donwload-spg install-spg clean-download-spg clean-install-spg clean-spg

export LIBS='-lgomp'
export CFLAGS='-fopenmp'

OBJS += spglib_f08.o\
        siesta2bt.o

LIBS += -fopenmp -lm

VPATH ?= $(shell pwd)
SPGPATH=$(VPATH)/spg

install-spg: $(SPGPATH)/spglib
	@ echo Target: $@
	@ ( cd $(SPGPATH)/spglib ; mkdir -p _build ; cd _build ; cmake .. ; make )

siesta2bt.o: spglib_f08.o
	@ echo Target: $@

libsymspg.a: $(SPGPATH)/spglib/_build/libsymspg.a $(SPGPATH)/spglib/example/spglib_f08.f90 $(SPGPATH)/siesta2bt.f90 $(SPGPATH)/siesta_analysis_spg.F
	@ echo Target: $@
	@ ( cd $(SPGPATH) ; cp siesta2bt.f90 .. ; \
	cp spglib/example/spglib_f08.f90 .. ; \
	cp ../siesta_analysis.F siesta_analysis_original.F ; \
	cp siesta_analysis_spg.F ../siesta_analysis.F )
	cp $(SPGPATH)/spglib/_build/libsymspg.a .

$(SPGPATH)/siesta2bt.f90: $(SPGPATH)/siesta2bt_spg_findsym.f90
	@ echo Target: $@
	cp $< $@


download-spg: 
	@ echo Target: $@
	@ ( cd $(SPGPATH) ; test -d spglib || git clone https://github.com/atztogo/spglib.git )

$(SPGPATH)/spglib: download-spg 
	@ echo Target: $@

$(SPGPATH)/spglib/_build/libsymspg.a $(SPGPATH)/spglib/example/spglib_f08.f90: install-spg
	@ echo Target: $@

$(SPGPATH)/siesta2bt_*.f90 $(SPGPATH)/siesta_analysis_spg.F:
	@ echo Target: $@
	@ ( cd $(SPGPATH) ; git clone git@github.com:lucas-campagna/spg-siesta.git )
	mv spg-siesta/* .
	rm -d spg-siesta


clean-spg: clean-install-spg
	rm -f $(VPATH)/siesta2bt.f90 libsymspg.a  $(VPATH)/spglib_f08.f90
	cd $(SPGPATH) ; cp siesta_analysis_original.F ../siesta_analysis.F

clean-install-spg:
	rm -fr $(SPGPATH)/spglib/_build

clean-download-spg:
	rm -fr $(SPGPATH)/spglib

