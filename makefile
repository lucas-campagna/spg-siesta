FC=gfortran

SRC=$(wildcard *.f90)
OBJ_SRC=$(patsubst %.f90,%.o,$(SRC))
#MOD_SRC=$(patsubst %.f90,%.mod,$(SRC))

OBJ=$(OBJ_SRC) 
MOD= $(wildcard ../*.mod)

LIBS=

EXE=spg.x

.PHONY: debug clean

all: $(OBJ)
	$(FC) -o $(EXE) $^ $(LIBS)

$(OBJ_SRC): %.o: %.f90
	$(FC) -c $< $(MOD) -L../

debug:
	@ echo "SRC = $(SRC)"
	@ echo "OBJ = $(OBJ)"
	@ echo "OBJ_SRC = $(OBJ_SRC)"
	@ echo "MOD_SRC = $(MOD_SRC)"
	@ echo "MOD = $(MOD)"

clean:
	rm -f $(EXE) $(OBJ_SRC) $(MOD_SRC)

