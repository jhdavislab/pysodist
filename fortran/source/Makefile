#SHELL = /bin/sh
.SUFFIXES: 
.SUFFIXES: .f .o

PGM = isodist

LOCDIR = ./
SRC = $(wildcard $(LOCDIR)*.f)
OBJ = $(addsuffix .o, $(basename $(notdir $(SRC))))
#JMT- FC = g77
FC = gfortran

#JMT- FFLAGS = -O3
FFLAGS = -O3 -ffixed-line-length-none

$(PGM): $(OBJ) 
	$(FC) $(FFLAGS) $(OBJ) -o $(PGM)
		
$(OBJ): $(PGM).h

.f.o: $(FC) $(FFLAGS) %.f

install:
	cp isodist ../.
	
clean:   
	rm isodist
	rm *.o


