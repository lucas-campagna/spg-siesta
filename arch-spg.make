.PHONY: donwload-spg install-spg clean-download-spg clean-install-spg clean-spg

OBJS += spglib_f08.o\
        siesta2bt.o

LIBS += -fopenmp -lm

#VPATH ?= $(shell pwd)
SPGPATH=$(shell find $(VPATH) -name arch-spg.make -printf "%h\n" | tail -n 1)


siesta2bt.o: spglib_f08.o

install-spg: $(SPGPATH)/spglib
	@ echo Target: $@
	@ ( cd $(SPGPATH); \
	 echo "export LIBS='-lgomp' ; export CFLAGS='-fopenmp'; cd spglib; mkdir -p _build ; cd _build ; cmake .. ; make" > script.sh; \
	 chmod +x $(SPGPATH)/script.sh; \
	 $(SPGPATH)/script.sh; \
	 rm $(SPGPATH)/script.sh )


libsymspg.a: $(SPGPATH)/spglib/_build/libsymspg.a $(SPGPATH)/spglib/example/spglib_f08.f90 $(SPGPATH)/siesta2bt.f90 $(SPGPATH)/siesta_analysis_spg.F
	@ echo Target: $@
	@ ( cd $(SPGPATH) ; cp siesta2bt.f90 .. ; \
	cp spglib/example/spglib_f08.f90 .. ; \
	cp ../siesta_analysis.F siesta_analysis_original.F ; \
	cp siesta_analysis_spg.F ../siesta_analysis.F )
	@ cp $(SPGPATH)/spglib/_build/libsymspg.a .

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
	cd $(SPGPATH) ; cp siesta_analysis_original.F ../siesta_analysis.F
	rm -f $(VPATH)/siesta2bt.f90 libsymspg.a  $(VPATH)/spglib_f08.f90 $(VPATH)/libsymspg.a $(SPGPATH)/siesta_analysis_original.F 

clean-install-spg:
	rm -fr $(SPGPATH)/spglib/_build

clean-download-spg:
	rm -fr $(SPGPATH)/spglib

